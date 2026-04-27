# AWS deployment for ImpAgingSim melt experiments

End-to-end pipeline: launch a GPU/CPU EC2 instance, sync code from GitHub, run the full set of multi-chain copolymer-melt scans (kappa scan + temperature scan + big production run), upload all results to S3, terminate.

## Prerequisites

1. AWS CLI configured: `aws configure` (credentials for IAM user with EC2 + S3 permissions)
2. SSH key pair `okezue` registered in your AWS account, private key at `~/Downloads/okezue.pem` (chmod 600)
3. An S3 bucket you own to receive results
4. A security group allowing SSH (port 22) inbound from your IP

## Files

| file | purpose |
|---|---|
| `launch.sh` | Provisions an EC2 instance and kicks off the experiment campaign. Asks for confirmation before billing. |
| `userdata.sh` | EC2 boot script: installs deps, clones repo, runs scans, syncs results to S3, shuts down. |
| `sync_back.sh` | Pulls completed results from S3 to your laptop. |
| `Dockerfile` | Optional containerized form (CUDA 12.4) for ECS / Batch / Sagemaker. |

## Usage

### Launch a full campaign

```bash
export S3_BUCKET=okezue-imp-aging-results          # your S3 bucket
export REGION=us-east-1                            # your region
export INSTANCE_TYPE=g5.2xlarge                    # GPU: $1.21/hr; or c6i.4xlarge for CPU at $0.68/hr
export AMI_ID=ami-0e2c8caa4b6378d8c                # Ubuntu 22.04 LTS in us-east-1
export KEY_NAME=okezue                             # AWS-side key name
export KEY_FILE=~/Downloads/okezue.pem
export SCAN_KIND=all                               # all | kappa | temperature | big
export GIT_REF=master                              # branch / tag / commit to run
export SECURITY_GROUP=default
# optional:
# export SUBNET_ID=subnet-xxxxxxxx
# export INSTANCE_PROFILE=ImpAgingSimRunner       # IAM role granting S3 write

bash aws/launch.sh
```

The script prints a summary, then prompts `Proceed and create instance? (yes/no)` before doing anything that costs money.

### What runs on the instance

`SCAN_KIND=all` runs three back-to-back campaigns:

1. **Kappa scan** — 6 kappa values × 4 seeds = 24 runs at 96 chains × 30 beads, 150k BD steps each. ~6-10 GPU-hours.
2. **Temperature scan** — 6 T values × 3 sequences × 3 seeds = 54 runs same size. ~12-20 GPU-hours.
3. **Big run** — single 256 chains × 30 beads, 400k steps, full position+grid recording. ~4-8 GPU-hours.

Total: ~25-40 GPU-hours on a g5.2xlarge → ~$30-50 in compute.

You can override via `SCAN_KIND=kappa`, `temperature`, or `big`.

### Pull results back

```bash
export S3_BUCKET=okezue-imp-aging-results
bash aws/sync_back.sh
```

Results land in `output/aws/<date>_<hostname>/melt/` mirroring the on-instance layout.

### Generate animations from a completed run

```bash
python3 -m melt.viz output/aws/.../big_correlated_*/   # produces density_slice.gif, polymer_3d.gif, Sk_evolution.gif, final_3d.png
python3 -m melt.analyze output/aws/.../scans/kscan_aws/   # cross-run plots + summary.csv
```

### Cost & safety

- Instance is launched with `--instance-initiated-shutdown-behavior terminate`. The userdata calls `shutdown -h +5` after the S3 sync, so the instance self-terminates ~5 min after the job completes. **No long-running zombie costs.**
- If you cancel the SSH session, the simulation keeps running (it's invoked from userdata, not from your shell).
- To force-kill: `aws ec2 terminate-instances --instance-ids i-XXXX --region $REGION`
- Check status: `aws ec2 describe-instances --filters "Name=tag:Project,Values=ImpAgingSim" --region $REGION --query 'Reservations[].Instances[].[InstanceId,State.Name,PublicIpAddress]' --output table`

### Resuming a failed run

The userdata is idempotent for git-clone (does `git pull` if repo exists). To resume after a partial failure, SSH in, `cd ImpAgingSim`, and re-run whichever scan command from `aws/userdata.sh` you need — no AMI rebuild needed.
