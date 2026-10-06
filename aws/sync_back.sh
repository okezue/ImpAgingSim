#!/bin/bash
set -euo pipefail
S3_BUCKET="${S3_BUCKET:?must export S3_BUCKET=...}"
REGION="${REGION:-us-east-1}"
LOCAL_DIR="${LOCAL_DIR:-output/aws}"
S3_PREFIX="${S3_PREFIX:-}"
SOURCE="s3://${S3_BUCKET}/${S3_PREFIX#/}"
mkdir -p "${LOCAL_DIR}"
echo "Syncing ${SOURCE} -> ${LOCAL_DIR}/ ..."
aws s3 sync "${SOURCE}" "${LOCAL_DIR}/" --region "${REGION}"
echo "Done."
du -sh "${LOCAL_DIR}"/* 2>/dev/null || true
