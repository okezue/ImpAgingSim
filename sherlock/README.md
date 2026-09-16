# Sherlock deployment for the eps_AB sweep

Runs `python -m melt.epsab_scan` on Stanford's Sherlock cluster (SLURM) as a GPU
job array, aggregates it as a CPU job, and pulls the results back to the laptop.
Code and venv live in `$HOME`; all outputs live in `$SCRATCH/impagingsim/`.

Sherlock reference: <https://www.sherlock.stanford.edu/docs/>

## Files

| file | where it runs | purpose |
|---|---|---|
| `setup_env.sh` | `sh_dev` node, once | python module + venv at `~/venvs/imp`, installs `requirements.txt`, `openmm[cuda12]`, `pytest` |
| `env.sh` | sourced by every job/shell | loads module, activates venv, detects the OpenMM platform, sets `CAMPAIGN_ROOT`, `cd`s into the repo |
| `smoke.sh` | `sh_dev -g 1` | unit tests + a two-run campaign (4000 steps) |
| `grid.sh` | sourced | turns `CAMPAIGN_ID`, `EPS_ABS`, `SEEDS`, `SIZES`, `EXTRA_ARGS` into driver flags, shared by the three scripts below |
| `epsab_manifest.sh` | login node | writes `manifest.json` once, before the array |
| `epsab_sweep.sbatch` | `gpu` partition, job array | one shard of the sweep per array task |
| `epsab_analyze.sbatch` | `normal` partition | `--analyze` |
| `sync_back.sh` | laptop | rsync from the data transfer nodes into `output/sherlock/<campaign>/` |

## Environment variables

| variable | default | used by |
|---|---|---|
| `CAMPAIGN_ID` | `epsab_kappa05` | manifest, sweep, analyze |
| `EPS_ABS`, `SEEDS`, `SIZES` | unset = driver defaults (37 points 1.0 -> 0.1, seeds 1 2 3 4, size 144) | manifest, sweep, analyze |
| `EXTRA_ARGS` | empty; appended verbatim, e.g. `"--n-steps 500000"` or `"--recover-interrupted"` | manifest, sweep, analyze |
| `CAMPAIGN_PLATFORM` | `CUDA`; `--platform` value used where there is no GPU to probe (manifest, analyze). Not part of the design hash, so it never has to match the GPU tasks | manifest, analyze |
| `OPENMM_PLATFORM` | probed by `env.sh` (CUDA, else OpenCL, else CPU); preset to skip the probe | sweep, smoke |
| `CAMPAIGN_ROOT` | `$SCRATCH/impagingsim` | everything on Sherlock |
| `REPO_DIR` | `$HOME/ImpAgingSim` | everything on Sherlock |
| `VENV` | `$HOME/venvs/imp` | `setup_env.sh`, `env.sh` |
| `PY_MODULE` | `python/3.12.1`, else best of `ml spider python` in 3.10 - 3.12 | `setup_env.sh`, `env.sh` |
| `ALLOW_INCOMPLETE` | `0`; `1` adds `--allow-incomplete-analysis` | analyze |
| `SUNETID`, `CAMPAIGN`, `FULL` | see step 7 | `sync_back.sh` (laptop) |

Design flags (`EPS_ABS`, `SEEDS`, `SIZES`, anything in `EXTRA_ARGS` such as
`--n-steps`) must be identical for manifest, sweep and analyze: designs are
hashed and the driver rejects a mismatch. A changed grid needs a new
`CAMPAIGN_ID`.

## 0. Before you start

- A Sherlock account (PI-sponsored) and your SUNet ID.
- The `sherlock/` kit committed and pushed. The driver refuses to write a
  manifest unless `git status` is clean, so nothing may be edited or copied
  into the clone on Sherlock; change code locally, push, `git pull` there.
- Nothing under `output/sherlock/` is gitignored yet; add `output/sherlock/`
  to `.gitignore` before syncing results to the laptop if you also run
  campaigns locally.

## 1. One-time local SSH config (laptop)

Sherlock uses SUNet password + Duo; SSH public keys are not accepted. Multiplexing
keeps one authenticated connection open for 4 h so `ssh`, `scp`, `rsync` do not
re-prompt. Add to `~/.ssh/config`:

```
Host login.sherlock.stanford.edu dtn.sherlock.stanford.edu
    User <sunetid>
    ControlMaster auto
    ControlPath ~/.ssh/%C
    ControlPersist 4h
```

`%C` (not `%r@%h:%p`) is required on macOS to avoid `unix_listener: ... too long`.

Host key fingerprints to accept on first connection:

```
RSA   SHA256:T1q1Tbq8k5XBD5PIxvlCfTxNMi1ORWwKNRPeZPXUfJA
ECDSA SHA256:eB0bODKdaCWtPgv0pYozsdC5ckfcBFVOxeMwrNKdkmg
```

Optional, skips the password (Duo still prompts): Kerberos. Put
<https://web.stanford.edu/dept/its/support/kerberos/dist/krb5.conf> at
`/etc/krb5.conf`, run `kinit <sunetid>@stanford.edu`, and add to the host block:

```
    GSSAPIAuthentication yes
    GSSAPIDelegateCredentials yes
```

Check the multiplexed master: `ssh -O check login.sherlock.stanford.edu`.

## 2. First login and checks

```bash
ssh login.sherlock.stanford.edu
sh_part                            # partitions: gpu, normal, dev, ...
node_feat -p gpu | grep GPU_        # GPU features usable with -C, e.g. GPU_MEM:16GB, GPU_BRD:TESLA
ml spider python                    # available python modules (expect python/3.12.1)
echo $SCRATCH $HOME                 # /scratch/users/<sunetid>  /home/users/<sunetid>
git clone https://github.com/okezue/ImpAgingSim.git ~/ImpAgingSim
cd ~/ImpAgingSim && git status      # must be clean
```

Login nodes are for editing and submitting only. Everything below that computes
runs on `sh_dev` or through `sbatch`.

## 3. One-time environment setup (compute node)

```bash
sh_dev                              # 1 core, 4 GB, 1 h; enough for pip
cd ~/ImpAgingSim
bash sherlock/setup_env.sh
exit
```

The script refuses to run on a login node (set `FORCE_LOGIN_NODE=1` to override).
It is idempotent: rerun it to upgrade packages. It prints
`python -m openmm.testInstallation` and the platform list; on a node without a
GPU, CUDA/OpenCL failures in that output are expected. The venv is ~1-2 GB in
`$HOME`; the pip cache is redirected to `$SCRATCH/.cache/pip`.

The module actually used is recorded in `~/venvs/imp/sherlock_python_module` so
`env.sh` always loads the interpreter the venv was built with.

## 4. Smoke test (MIG GPU slice)

```bash
sh_dev -g 1                         # lightweight MIG slice, usually instant
cd ~/ImpAgingSim
bash sherlock/smoke.sh
exit
```

Expected: tests pass, `env: ... platform=CUDA`, two runs complete under
`$SCRATCH/impagingsim/smoke/smoke/runs/`, `analysis/` is listed at the end.
Total a few minutes, mostly OpenMM kernel compilation on first use.

For a full GPU instead of a slice: `sh_dev -c 4 -m 8GB -g 1 -p gpu`.
For the CPU path: plain `sh_dev` (expect `platform=CPU`).

## 5. Stage 1 campaign: default 37 x 4 grid (148 runs)

### 5a. Manifest (login node, once per campaign)

```bash
cd ~/ImpAgingSim && git pull && git status
bash sherlock/epsab_manifest.sh
```

Writes `$SCRATCH/impagingsim/epsab_kappa05/manifest.json`, pins the git
commit, prints the run count. Rerunning is a no-op if the design is unchanged.

### 5b. Submit the array

```bash
sbatch sherlock/epsab_sweep.sbatch                 # 8 tasks, -p gpu -G 1 -c 4 --mem=16G -t 02:00:00
```

Task `K` of `S` runs plan indices `i` with `i % S == K` (`--shard-index`,
`--shard-count` from `$SLURM_ARRAY_TASK_ID` / `$SLURM_ARRAY_TASK_COUNT`), so the
array size is free to change:

```bash
sbatch --array=0-15 sherlock/epsab_sweep.sbatch    # 16 shards, ~9 runs each
sbatch --array=0-15%8 sherlock/epsab_sweep.sbatch  # at most 8 running at once
```

Grid overrides travel as environment variables; they must match the manifest:

```bash
sbatch --export=ALL,CAMPAIGN_ID=epsab_kappa05,SEEDS="1 2 3 4" sherlock/epsab_sweep.sbatch
# equivalent (sbatch exports the caller's environment by default):
SEEDS="1 2 3 4" sbatch sherlock/epsab_sweep.sbatch
```

Submit from the repo root: `slurm-<jobid>_<task>.out` logs are written there and
`*.out` is gitignored, so they do not dirty the tree.

### 5c. Monitor

```bash
squeue -u $USER                                          # queue state; ST=PD pending, R running
squeue -u $USER --start                                  # estimated start times
tail -f slurm-<jobid>_0.out                              # one shard's log
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS     # after completion, per task
find $SCRATCH/impagingsim/epsab_kappa05/runs -name completion.json | wc -l   # completed runs / 148
ls $SCRATCH/impagingsim/epsab_kappa05/.locks             # runs currently in flight
du -sh $SCRATCH/impagingsim/epsab_kappa05
scancel <jobid>                                          # whole array; scancel <jobid>_3 for one task
```

The driver's own report (light, fine on a login node; same grid variables):

```bash
source sherlock/env.sh
python -m melt.epsab_scan --status --out $CAMPAIGN_ROOT --campaign-id epsab_kappa05   # n_completed, missing_indices
```

### 5d. Recover and resubmit (idempotent)

A task killed at the time limit or by a node failure leaves one run's lock file
and staging directory behind. Once `squeue -u $USER` shows no `epsab` tasks:

```bash
sbatch --export=ALL,EXTRA_ARGS="--recover-interrupted" sherlock/epsab_sweep.sbatch
```

Completed runs are checksum-verified and skipped, so resubmitting the array any
number of times only does the missing work. Without `--recover-interrupted` the
driver stops at a stale lock with `lock exists for ...`.

## 6. Analysis job

```bash
sbatch sherlock/epsab_analyze.sbatch                                 # -p normal -c 4 --mem=16G -t 01:00:00
sbatch --dependency=afterok:<array jobid> sherlock/epsab_analyze.sbatch   # queue it right after the sweep
sbatch --export=ALL,ALLOW_INCOMPLETE=1 sherlock/epsab_analyze.sbatch      # aggregate whatever is complete
```

Outputs go to `$SCRATCH/impagingsim/epsab_kappa05/analysis/`; the log
`slurm-<jobid>.out` lists them. The job passes the same grid variables and
`--platform $CAMPAIGN_PLATFORM` (no dynamics run; the flag only has to match the
manifest).

## 7. Pull results back (laptop)

```bash
SUNETID=<sunetid> bash sherlock/sync_back.sh            # manifest, analysis/, runs/*/{completion,meta}.json, snapshots.csv
SUNETID=<sunetid> FULL=1 bash sherlock/sync_back.sh     # everything incl. mode_amplitudes.npz (~3-4 GB)
SUNETID=<sunetid> CAMPAIGN=epsab_stage2 bash sherlock/sync_back.sh
```

Transfers go through `dtn.sherlock.stanford.edu` (rsync/scp only, no shell);
`$SCRATCH` resolves to `/scratch/users/<sunetid>`. Results land in
`output/sherlock/<campaign>/`. `.staging/`, `.locks/`, `.interrupted/` are never
pulled. `$SCRATCH` is purged after inactivity: pull `FULL=1` once the campaign
is done, or copy it to `$GROUP_HOME`/`$OAK` on Sherlock.

## 8. Stage 2: refinement

Stage 1 locates the transition; Stage 2 refines it with a denser `eps_AB`
grid, more seeds and a second system size, under a new campaign id. Example for
a transition between 0.5 and 0.3:

```bash
# login node
export CAMPAIGN_ID=epsab_stage2
export EPS_ABS="0.50 0.49 0.48 0.47 0.46 0.45 0.44 0.43 0.42 0.41 0.40 0.39 0.38 0.37 0.36 0.35 0.34 0.33 0.32 0.31 0.30"
export SEEDS="1 2 3 4 5 6 7 8"
export SIZES="144 288"
bash sherlock/epsab_manifest.sh                          # 21 x 8 x 2 = 336 runs
sbatch --array=0-23%8 sherlock/epsab_sweep.sbatch        # ~14 runs per task; exported env carries the grid
sbatch --dependency=afterok:<jobid> sherlock/epsab_analyze.sbatch
```

Or all on the `sbatch` line with `--export=ALL,CAMPAIGN_ID=epsab_stage2,EPS_ABS="...",SEEDS="...",SIZES="144 288"`.
Same variables again for `sync_back.sh` (`CAMPAIGN=epsab_stage2`). Do not reuse
`epsab_kappa05`: the manifest hash differs and the driver will refuse.

## 9. Expected wall time and storage

Per 144-chain run: ~30 s of GPU time plus ~1 min of Python for the Fourier-mode
recording (measured on an AWS A10G); 288 chains roughly 2-3x. Sherlock GPU
types vary, keep a 2x margin. Output per run ~15-25 MB (`mode_amplitudes.npz`).

| campaign | runs | GPU-hours | wall time, 8 tasks | storage in `$SCRATCH` |
|---|---|---|---|---|
| smoke (2 runs, 4000 steps) | 2 | < 0.1 | 2-5 min interactive | < 0.1 GB |
| Stage 1, 37 x 4 x {144} | 148 | 2-4 | 20-45 min per task + queue | 3-4 GB |
| Stage 2 example, 21 x 8 x {144, 288} | 336 | 12-18 | 24 tasks x ~45-60 min | ~10 GB |
| venv | | | | 1-2 GB in `$HOME` (15 GB quota) |

Queue wait on `gpu` is the unknown; `squeue -u $USER --start` shows estimates.

## 10. Troubleshooting

| symptom | fix |
|---|---|
| `env: ... platform=OpenCL` or `CPU` on a GPU node | CUDA platform failed. `nvidia-smi`, then `python -m openmm.testInstallation` in `sh_dev -g 1` for the error. Rerun `setup_env.sh` if `openmm-cuda-12` is missing. The sweep falls back to OpenCL (then CPU) automatically and the log shows a `WARNING`; OpenCL results are valid, CPU is ~10-50x slower, so `scancel` and fix the install instead. |
| `git tree is dirty` | `git status`. Untracked files count. Delete them or `git pull`; never edit on Sherlock. A manual `pytest` leaves `.pytest_cache/`: `rm -rf .pytest_cache` (`smoke.sh` avoids it with `-p no:cacheprovider`). |
| `current commit ... differs from campaign commit` | Code changed after the manifest was written. `git checkout <creation commit from manifest.json>` to finish the campaign, or start a new `CAMPAIGN_ID`. |
| `lock exists for <run>` / `staging directory exists` | A task died mid-run. Confirm nothing is running (`squeue -u $USER`), then resubmit with `EXTRA_ARGS="--recover-interrupted"` (step 5d). |
| `no manifest at ...` in the job log | Run `bash sherlock/epsab_manifest.sh` with the same `CAMPAIGN_ID`/grid, then resubmit. |
| `campaign design differs from existing manifest` | Grid or `EXTRA_ARGS` differ from the manifest. Use the same variables, or a new `CAMPAIGN_ID`. |
| tasks pending for hours | Request less (`-t 01:00:00`, fewer concurrent tasks with `--array=0-7%4`) or target a less busy GPU type with `-C` from `node_feat -p gpu \| grep GPU_`, e.g. `sbatch -C GPU_MEM:16GB ...`. |
| `sbatch: error: ... GPU` / job rejected | GPU jobs need both `-p gpu` and `-G 1`; the script sets them, do not override one without the other. |
| task ends `DUE TO TIME LIMIT` | Fine: resubmit (step 5d). Or `sbatch -t 04:00:00 ...`, max 2 days on `gpu`. |
| `Out Of Memory` in `sacct` | `sbatch --mem=32G sherlock/epsab_sweep.sbatch`; 288-chain runs need more than 144. |
| `python: command not found` / wrong version in a job | `env.sh` could not load the module. `ml spider python`, then `PY_MODULE=python/<version> bash sherlock/setup_env.sh` on `sh_dev` and export `PY_MODULE` for jobs. |
| Duo prompt on every command | Multiplexing not active: check `~/.ssh/config` block in step 1, `ssh -O check login.sherlock.stanford.edu`. |
| `unix_listener: ... too long` | Use `ControlPath ~/.ssh/%C`. |
| `$HOME` full | `du -sh ~/venvs ~/.cache`; the venv belongs in `$HOME`, everything else in `$SCRATCH`. |
| smoke fails with a commit mismatch | `rm -rf $SCRATCH/impagingsim/smoke` (or `SMOKE_ID=smoke2 bash sherlock/smoke.sh`). |
