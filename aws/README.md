# Run the paper's static simulations on AWS

These scripts run the parameter sweeps and fixed-density validation used in
the manuscript, *Independent control of sequence-correlation amplitude tunes
post-quench composition fluctuations in A/B heteropolymer melts*.
The run-to-figure mapping and archived inputs are in
[RESULTS_MANIFEST.md](RESULTS_MANIFEST.md).

## Run in an existing environment

From the repository root, with the dependencies in `requirements.txt` installed:

```bash
PY=python PLATFORM_FLAG="--platform CUDA" bash aws/rerun_corrected.sh
PY=python PLATFORM_FLAG="--platform CUDA" bash aws/fixed_density.sh
```

The first command runs 744 simulations: the four amplitude controls, the
persistence/amplitude grid, the off-stoichiometric diagnostic and the
sequence-class/cross-attraction scan. The second runs the 30 matched
fixed-density size comparisons. Use `--platform CPU` when CUDA is unavailable.
The scripts write to `output/melt/scans_corrected/` and `output/melt/fd_v2/`,
respectively. Existing completed runs are reused; fixed-density runs additionally
verify their recorded hashes. Set `S3_BUCKET` to upload results after a campaign.

## Launch an EC2 instance

Configure the AWS CLI, an SSH key pair, an S3 bucket, and a security group that
permits SSH from your address. Choose an Ubuntu GPU AMI valid in your region
with drivers compatible with the CUDA environment created by `userdata.sh`.
The launcher requires these configuration values rather than assuming an old AMI
or account-specific key. Set `GIT_REF` to the reviewed commit or tag to run.

```bash
export REGION=us-east-1
export INSTANCE_TYPE=g5.2xlarge
export AMI_ID=your-regional-ubuntu-gpu-ami
export KEY_NAME=your-ec2-key-name
export KEY_FILE=/path/to/your-key.pem
export SECURITY_GROUP=your-security-group
export S3_BUCKET=your-paper-results-bucket
export GIT_REF=your-reviewed-commit-or-tag
export SCAN_KIND=rerun
bash aws/launch.sh
```

`launch.sh` displays the configuration and asks before creating the instance.
`userdata.sh` installs the simulation environment, checks CUDA availability,
runs the selected workload, uploads results, and schedules shutdown after a
successful upload. EC2 is configured to terminate on shutdown. A failure before
that final step can leave the instance running; inspect its status and log.

| `SCAN_KIND` | Workload |
|---|---|
| `smoke` | Small setup/trajectory check |
| `kappa` | Baseline static amplitude scan |
| `rerun` | The 744-run static parameter suite |
| `fixed_density` | The 30-run fixed-density size comparison |
| `all` | Static suite followed by fixed-density comparison |

Results and logs use date/hostname prefixes in your bucket. The launcher prints
the instance ID, SSH command and bootstrap-log command. To terminate an instance:

```bash
aws ec2 terminate-instances --instance-ids YOUR_INSTANCE_ID --region "$REGION"
```

## Retrieve and analyse results

Select a prefix containing the paper campaign you need:

```bash
S3_BUCKET=your-paper-results-bucket S3_PREFIX=your-campaign-prefix \
  LOCAL_DIR=output/aws/paper bash aws/sync_back.sh

python -m melt.analyze output/melt/scans_corrected/fig1_baseline
python -m melt.fixed_density_size_scan --analyze \
  --out output/melt/fd_v2 --campaign-id fixed_density_pi099_v2
```

The checked-in completed fixed-density data are under
`output/melt/fixed_density_size/fixed_density_pi099_v2/`; preserve them when
running a new campaign. Rebuild the submission figures from the retained source
tables with `python scripts/paper/figures.py --out output/paper_figures`.
The RPA correction instructions are under [`docs/corrections/scientific_reports_rpa/`](../docs/corrections/scientific_reports_rpa/).
