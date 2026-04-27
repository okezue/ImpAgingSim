from __future__ import annotations
import argparse,os,time
from types import SimpleNamespace
from .run import execute

def make_cfg(scan_root,scan_id,seq,seed,base):
    c=SimpleNamespace(**vars(base))
    c.sequence=seq
    c.seed=seed
    c.out=os.path.join(scan_root,scan_id)
    c.run_id=f"{seq}_s{seed}"
    return c

def parse_args():
    p=argparse.ArgumentParser(description="Sequence-correlation scan for copolymer melt quench")
    p.add_argument("--out",type=str,default="output/melt/scans")
    p.add_argument("--scan_id",type=str,default=None)
    p.add_argument("--sequences",type=str,nargs="+",
                   default=["random","block","alternating","correlated"])
    p.add_argument("--seeds",type=int,nargs="+",default=[1,2,3])
    p.add_argument("--n_chains",type=int,default=24)
    p.add_argument("--chain_length",type=int,default=16)
    p.add_argument("--box_size",type=float,default=10.0)
    p.add_argument("--f_A",type=float,default=0.5)
    p.add_argument("--block_length",type=int,default=4)
    p.add_argument("--kappa",type=float,default=0.7)
    p.add_argument("--pi",type=float,default=0.9)
    p.add_argument("--bond_k",type=float,default=200.0)
    p.add_argument("--bond_r0",type=float,default=1.0)
    p.add_argument("--lj_eps_AA",type=float,default=1.0)
    p.add_argument("--lj_eps_BB",type=float,default=1.0)
    p.add_argument("--lj_eps_AB",type=float,default=0.1)
    p.add_argument("--lj_sigma",type=float,default=1.0)
    p.add_argument("--lj_cutoff",type=float,default=2.5)
    p.add_argument("--T_equilibrate",type=float,default=5.0)
    p.add_argument("--T_quench",type=float,default=0.7)
    p.add_argument("--temperature",type=float,default=1.0)
    p.add_argument("--friction",type=float,default=1.0)
    p.add_argument("--dt",type=float,default=0.005)
    p.add_argument("--n_steps",type=int,default=20000)
    p.add_argument("--equilibration",type=int,default=5000)
    p.add_argument("--snapshot_interval",type=int,default=400)
    p.add_argument("--platform",type=str,default="Reference")
    p.add_argument("--grid_size",type=int,default=24)
    return p.parse_args()

def main():
    a=parse_args()
    a.compute_density=True
    scan_id=a.scan_id or f"scan_{int(time.time())}"
    print(f"Scan id: {scan_id}")
    print(f"Sequences: {a.sequences}")
    print(f"Seeds: {a.seeds}")
    print(f"Quench: T={a.T_equilibrate} -> T={a.T_quench}, eps_AB={a.lj_eps_AB}")
    print(f"System: {a.n_chains} chains x {a.chain_length} beads = {a.n_chains*a.chain_length} beads, box={a.box_size}")
    print(f"Steps: {a.equilibration} eq + {a.n_steps} prod, snapshot every {a.snapshot_interval}\n")
    n_total=len(a.sequences)*len(a.seeds);k=0
    t0=time.time()
    for seq in a.sequences:
        for seed in a.seeds:
            k+=1
            cfg=make_cfg(a.out,scan_id,seq,seed,a)
            t1=time.time()
            rd=execute(cfg)
            print(f"  [{k}/{n_total}] {seq} seed={seed} -> {rd}  ({time.time()-t1:.1f}s)")
    print(f"\nScan complete in {time.time()-t0:.1f}s -> {os.path.join(a.out,scan_id)}")

if __name__=="__main__":
    main()
