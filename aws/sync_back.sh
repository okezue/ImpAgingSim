#!/bin/bash
set -euo pipefail
S3_BUCKET="${S3_BUCKET:?must export S3_BUCKET=...}"
REGION="${REGION:-us-east-1}"
LOCAL_DIR="${LOCAL_DIR:-output/aws}"
mkdir -p "${LOCAL_DIR}"
echo "Syncing s3://${S3_BUCKET}/ -> ${LOCAL_DIR}/ ..."
aws s3 sync "s3://${S3_BUCKET}/" "${LOCAL_DIR}/" --region "${REGION}"
echo "Done."
du -sh "${LOCAL_DIR}"/* 2>/dev/null || true
