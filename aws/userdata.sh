#!/bin/bash
set -euxo pipefail
exec > >(tee -a /var/log/imp-bootstrap.log) 2>&1

REGION="${REGION:-us-east-1}"
S3_BUCKET="${S3_BUCKET}"
SCAN_KIND="${SCAN_KIND:-all}"
GIT_REF="${GIT_REF:-master}"
N_THREADS="${N_THREADS:-$(nproc)}"

if [ -n "${AWS_ACCESS_KEY_ID:-}" ]; then
  mkdir -p /home/ubuntu/.aws
  cat > /home/ubuntu/.aws/credentials <<EOF
[default]
aws_access_key_id=${AWS_ACCESS_KEY_ID}
aws_secret_access_key=${AWS_SECRET_ACCESS_KEY}
aws_session_token=${AWS_SESSION_TOKEN}
EOF
  cat > /home/ubuntu/.aws/config <<EOF
[default]
region=${REGION}
EOF
  chown -R ubuntu:ubuntu /home/ubuntu/.aws
  chmod 600 /home/ubuntu/.aws/credentials
fi

apt-get update -y
apt-get install -y git ffmpeg awscli wget

cd /home/ubuntu
REPO_URL="${REPO_URL:-https://github.com/okezue/ImpAgingSim.git}"
REPO_DIR="${REPO_DIR:-ImpAgingSim}"
sudo -u ubuntu git clone "${REPO_URL}" "${REPO_DIR}" || (cd "${REPO_DIR}" && sudo -u ubuntu git fetch --all && sudo -u ubuntu git reset --hard origin/${GIT_REF})
cd "${REPO_DIR}"
sudo -u ubuntu git checkout "${GIT_REF}"

if [ ! -d /home/ubuntu/miniconda3 ]; then
  echo "=== Installing miniconda ==="
  sudo -u ubuntu wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
  sudo -u ubuntu bash /tmp/miniconda.sh -b -p /home/ubuntu/miniconda3
fi
CONDA=/home/ubuntu/miniconda3/bin/conda
PY=/home/ubuntu/miniconda3/envs/imp/bin/python

if [ ! -d /home/ubuntu/miniconda3/envs/imp ]; then
  echo "=== Creating imp env with openmm + CUDA 12.9 (matches openmm 8.5 conda-forge requirement; <13.0 driver max) ==="
  sudo -u ubuntu ${CONDA} create -n imp --override-channels -c conda-forge -y \
    python=3.11 'openmm>=8.5' 'cuda-version=12.9' numpy scipy pandas matplotlib
fi

echo "=== GPU + CUDA platform check ==="
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv 2>&1 | head -3 || true
sudo -u ubuntu ${PY} -c "import openmm as mm; print('platforms:', [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])" 2>&1

PLATFORM_FLAG="--platform CUDA"
if ! sudo -u ubuntu ${PY} -c "import openmm as mm;sys=mm.System();sys.addParticle(1.0);i=mm.LangevinMiddleIntegrator(1.0,1.0,0.005);ctx=mm.Context(sys,i,mm.Platform.getPlatformByName('CUDA'));print('CUDA OK',ctx.getPlatform().getPropertyValue(ctx,'DeviceName'))" 2>&1 | tee -a /var/log/imp-bootstrap.log | grep -q "CUDA OK"; then
  echo "WARNING: CUDA not available, falling back to CPU"
  PLATFORM_FLAG="--platform CPU"
fi
echo "PLATFORM_FLAG=${PLATFORM_FLAG}"

mkdir -p /home/ubuntu/${REPO_DIR}/output/melt
chown -R ubuntu:ubuntu /home/ubuntu/${REPO_DIR}

run_kappa(){
  sudo -u ubuntu ${PY} -m melt.kappa_scan \
    --scan_id kscan_aws \
    --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
    --seeds 1 2 3 4 5 \
    --n_chains 144 --chain_length 40 --box_size 22.0 \
    --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
    --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
    --grid_size 56 ${PLATFORM_FLAG}
}
run_temperature(){
  sudo -u ubuntu ${PY} -m melt.temperature_scan \
    --scan_id tscan_aws \
    --T_quenches 0.25 0.35 0.5 0.65 0.8 1.0 1.3 1.7 \
    --sequences random block correlated \
    --seeds 1 2 3 4 5 \
    --n_chains 144 --chain_length 40 --box_size 22.0 \
    --lj_eps_AB 0.1 \
    --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
    --grid_size 56 ${PLATFORM_FLAG}
}
run_big(){
  for seq in correlated random block; do
    sudo -u ubuntu ${PY} -m melt.big_run \
      --sequence ${seq} \
      --n_chains 400 --chain_length 50 --box_size 32.0 \
      --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
      --n_steps 500000 --equilibration 60000 --snapshot_interval 5000 \
      --grid_size 80 \
      --run_id big_${seq} ${PLATFORM_FLAG}
  done
  for d in output/melt/big/big_*; do
    sudo -u ubuntu ${PY} -m melt.viz "$d" || true
  done
}
run_smoke(){
  sudo -u ubuntu ${PY} -m melt.run \
    --sequence correlated --n_chains 32 --chain_length 20 --box_size 12.0 \
    --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
    --equilibration 2000 --n_steps 8000 --snapshot_interval 500 \
    --grid_size 24 --compute_density --save_trajectory \
    --out output/melt --run_id aws_smoke ${PLATFORM_FLAG}
  sudo -u ubuntu ${PY} -m melt.viz output/melt/aws_smoke || true
}
run_rerun(){
  sudo -u ubuntu env PY="${PY}" PLATFORM_FLAG="${PLATFORM_FLAG}" \
       S3_BUCKET="${S3_BUCKET}" REGION="${REGION}" \
       bash aws/rerun_corrected.sh
}
run_seed_extension(){
  sudo -u ubuntu env PY="${PY}" PLATFORM_FLAG="${PLATFORM_FLAG}" \
       S3_BUCKET="${S3_BUCKET}" REGION="${REGION}" \
       bash aws/seed_extension.sh
}

case "${SCAN_KIND}" in
  smoke) run_smoke ;;
  kappa) run_kappa ;;
  temperature) run_temperature ;;
  big) run_big ;;
  rerun) run_rerun ;;
  seed_extension) run_seed_extension ;;
  all) run_kappa; run_temperature; run_big ;;
  *) echo "unknown SCAN_KIND=${SCAN_KIND}"; exit 2 ;;
esac

cd /home/ubuntu/${REPO_DIR}
sudo -u ubuntu aws s3 sync output/ "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_$(hostname -s)/" --region "${REGION}"

shutdown -h +5
