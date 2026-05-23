#!/bin/bash
set -euo pipefail

REGION="${REGION:-us-east-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-g5.2xlarge}"
AMI_ID="${AMI_ID:-ami-0ff290337e78c83bf}"
KEY_NAME="${KEY_NAME:-okezue}"
KEY_FILE="${KEY_FILE:-$HOME/Downloads/okezue.pem}"
S3_BUCKET="${S3_BUCKET:?must export S3_BUCKET=...}"
SCAN_KIND="${SCAN_KIND:-all}"
GIT_REF="${GIT_REF:-master}"
SECURITY_GROUP="${SECURITY_GROUP:-default}"
SUBNET_ID="${SUBNET_ID:-}"
TAG_NAME="${TAG_NAME:-imp-aging-$(date +%Y%m%d-%H%M%S)}"
STS_DURATION="${STS_DURATION:-43200}"
SKIP_CONFIRM="${SKIP_CONFIRM:-no}"

if [ ! -f "${KEY_FILE}" ]; then
  echo "ERROR: KEY_FILE not found: ${KEY_FILE}" >&2
  exit 1
fi
chmod 600 "${KEY_FILE}"

echo "=== Launch summary ==="
echo "  region:        ${REGION}"
echo "  instance type: ${INSTANCE_TYPE}"
echo "  AMI:           ${AMI_ID}"
echo "  key:           ${KEY_NAME} (${KEY_FILE})"
echo "  s3 bucket:     ${S3_BUCKET}"
echo "  scan kind:     ${SCAN_KIND}"
echo "  git ref:       ${GIT_REF}"
echo "  security grp:  ${SECURITY_GROUP}"
echo "  subnet:        ${SUBNET_ID:-<default>}"
echo "  tag:           ${TAG_NAME}"
echo "  STS duration:  ${STS_DURATION}s"
echo "======================"
if [ "${SKIP_CONFIRM}" != "yes" ]; then
  read -p "Proceed and create instance? (yes/no) " confirm
  if [ "${confirm}" != "yes" ]; then
    echo "Aborted."
    exit 0
  fi
fi

echo "Acquiring STS session token (${STS_DURATION}s)..."
STS_JSON=$(aws sts get-session-token --duration-seconds "${STS_DURATION}" --output json)
STS_KEY=$(echo "${STS_JSON}" | python3 -c "import sys,json;print(json.load(sys.stdin)['Credentials']['AccessKeyId'])")
STS_SECRET=$(echo "${STS_JSON}" | python3 -c "import sys,json;print(json.load(sys.stdin)['Credentials']['SecretAccessKey'])")
STS_TOKEN=$(echo "${STS_JSON}" | python3 -c "import sys,json;print(json.load(sys.stdin)['Credentials']['SessionToken'])")
echo "  got temp creds expiring: $(echo ${STS_JSON} | python3 -c "import sys,json;print(json.load(sys.stdin)['Credentials']['Expiration'])")"

USERDATA_FILE="$(dirname "$0")/userdata.sh"
RENDERED=$(mktemp)
trap "rm -f ${RENDERED}" EXIT
{
  echo "#!/bin/bash"
  echo "export REGION='${REGION}'"
  echo "export S3_BUCKET='${S3_BUCKET}'"
  echo "export SCAN_KIND='${SCAN_KIND}'"
  echo "export GIT_REF='${GIT_REF}'"
  echo "export AWS_ACCESS_KEY_ID='${STS_KEY}'"
  echo "export AWS_SECRET_ACCESS_KEY='${STS_SECRET}'"
  echo "export AWS_SESSION_TOKEN='${STS_TOKEN}'"
  tail -n +2 "${USERDATA_FILE}"
} > "${RENDERED}"

EXTRA_ARGS=()
if [ -n "${SUBNET_ID}" ]; then
  EXTRA_ARGS+=(--subnet-id "${SUBNET_ID}")
fi
EXTRA_ARG_STR=""
if [ ${#EXTRA_ARGS[@]} -gt 0 ]; then
  EXTRA_ARG_STR="${EXTRA_ARGS[@]}"
fi

ROOT_VOLUME_GB="${ROOT_VOLUME_GB:-100}"
INSTANCE_ID=$(aws ec2 run-instances \
  --region "${REGION}" \
  --image-id "${AMI_ID}" \
  --instance-type "${INSTANCE_TYPE}" \
  --key-name "${KEY_NAME}" \
  --security-groups "${SECURITY_GROUP}" \
  --user-data "file://${RENDERED}" \
  --instance-initiated-shutdown-behavior terminate \
  --block-device-mappings "[{\"DeviceName\":\"/dev/sda1\",\"Ebs\":{\"VolumeSize\":${ROOT_VOLUME_GB},\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=${TAG_NAME}},{Key=Project,Value=ImpAgingSim}]" \
  ${EXTRA_ARG_STR} \
  --query 'Instances[0].InstanceId' --output text)

echo "Launched: ${INSTANCE_ID}"
echo "Polling for public IP..."
while true; do
  IP=$(aws ec2 describe-instances --region "${REGION}" --instance-ids "${INSTANCE_ID}" \
       --query 'Reservations[0].Instances[0].PublicIpAddress' --output text 2>/dev/null || echo "None")
  if [ "${IP}" != "None" ] && [ -n "${IP}" ]; then break; fi
  sleep 5
done
echo "Public IP: ${IP}"
echo
echo "Connect: ssh -i ${KEY_FILE} ubuntu@${IP}"
echo "Tail logs: ssh -i ${KEY_FILE} ubuntu@${IP} 'sudo tail -f /var/log/imp-bootstrap.log'"
echo "Results will sync to: s3://${S3_BUCKET}/"
echo "Instance auto-terminates 5 min after job completes."
