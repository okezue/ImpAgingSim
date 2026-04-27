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
apt-get install -y python3-pip git ffmpeg awscli

cd /home/ubuntu
sudo -u ubuntu git clone https://github.com/okezue/ImpAgingSim.git || (cd ImpAgingSim && sudo -u ubuntu git pull)
cd ImpAgingSim
sudo -u ubuntu git checkout "${GIT_REF}"
sudo -u ubuntu pip3 install --user -r requirements.txt

echo "=== GPU + CUDA platform check ==="
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv 2>&1 | head -3 || true
sudo -u ubuntu python3 -c "import openmm as mm; print('platforms:', [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])" 2>&1

PLATFORM_FLAG="--platform CUDA"
if ! sudo -u ubuntu python3 -c "import openmm as mm;p=mm.Platform.getPlatformByName('CUDA');print('CUDA OK')" 2>&1 | tee -a /var/log/imp-bootstrap.log | grep -q "CUDA OK"; then
  echo "WARNING: CUDA not available, falling back to CPU"
  PLATFORM_FLAG="--platform CPU"
fi
echo "PLATFORM_FLAG=${PLATFORM_FLAG}"

mkdir -p /home/ubuntu/ImpAgingSim/output/melt
chown -R ubuntu:ubuntu /home/ubuntu/ImpAgingSim

run_kappa(){
  sudo -u ubuntu python3 -m melt.kappa_scan \
    --scan_id kscan_aws \
    --kappas 0.0 0.2 0.4 0.6 0.8 1.0 \
    --seeds 1 2 3 4 \
    --n_chains 96 --chain_length 30 --box_size 18.0 \
    --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
    --n_steps 150000 --equilibration 20000 --snapshot_interval 1500 \
    --grid_size 48 ${PLATFORM_FLAG}
}
run_temperature(){
  sudo -u ubuntu python3 -m melt.temperature_scan \
    --scan_id tscan_aws \
    --T_quenches 0.3 0.5 0.7 1.0 1.5 2.0 \
    --sequences random block correlated \
    --seeds 1 2 3 \
    --n_chains 96 --chain_length 30 --box_size 18.0 \
    --lj_eps_AB 0.1 \
    --n_steps 150000 --equilibration 20000 --snapshot_interval 1500 \
    --grid_size 48 ${PLATFORM_FLAG}
}
run_big(){
  sudo -u ubuntu python3 -m melt.big_run \
    --sequence correlated \
    --n_chains 256 --chain_length 30 --box_size 24.0 \
    --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
    --n_steps 400000 --equilibration 50000 --snapshot_interval 4000 \
    --grid_size 64 ${PLATFORM_FLAG}
  for d in output/melt/big/big_*; do
    sudo -u ubuntu python3 -m melt.viz "$d" || true
  done
}
run_smoke(){
  sudo -u ubuntu python3 -m melt.run \
    --sequence correlated --n_chains 16 --chain_length 12 --box_size 8.0 \
    --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
    --equilibration 1000 --n_steps 3000 --snapshot_interval 200 \
    --grid_size 16 --compute_density --save_trajectory \
    --out output/melt --run_id aws_smoke ${PLATFORM_FLAG}
  sudo -u ubuntu python3 -m melt.viz output/melt/aws_smoke || true
}

case "${SCAN_KIND}" in
  smoke) run_smoke ;;
  kappa) run_kappa ;;
  temperature) run_temperature ;;
  big) run_big ;;
  all) run_kappa; run_temperature; run_big ;;
  *) echo "unknown SCAN_KIND=${SCAN_KIND}"; exit 2 ;;
esac

cd /home/ubuntu/ImpAgingSim
sudo -u ubuntu aws s3 sync output/ "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_$(hostname -s)/" --region "${REGION}"

shutdown -h +5
