from __future__ import annotations
import argparse,os,time
from types import SimpleNamespace
from .run import execute

def parse_args():
    p=argparse.ArgumentParser(description="Kappa scan for correlated copolymer melt")
    p.add_argument("--out",type=str,default="output/melt/scans")
    p.add_argument("--scan_id",type=str,default=None)
    p.add_argument("--kappas",type=float,nargs="+",default=[0.0,0.2,0.4,0.6,0.8,1.0])
    p.add_argument("--seeds",type=int,nargs="+",default=[1,2,3])
    p.add_argument("--n_chains",type=int,default=48)
    p.add_argument("--chain_length",type=int,default=24)
    p.add_argument("--box_size",type=float,default=14.0)
    p.add_argument("--f_A",type=float,default=0.5)
    p.add_argument("--block_length",type=int,default=4)
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
    p.add_argument("--n_steps",type=int,default=80000)
    p.add_argument("--equilibration",type=int,default=10000)
    p.add_argument("--snapshot_interval",type=int,default=800)
    p.add_argument("--platform",type=str,default=None)
    p.add_argument("--grid_size",type=int,default=48)
    p.add_argument("--save_trajectory",action="store_true")
    p.add_argument("--save_density_grids",action="store_true")
    return p.parse_args()

def main():
    a=parse_args()
    a.compute_density=True
    scan_id=a.scan_id or f"kscan_{int(time.time())}"
    print(f"Kappa scan: id={scan_id}, kappas={a.kappas}, seeds={a.seeds}")
    print(f"System: {a.n_chains}x{a.chain_length}={a.n_chains*a.chain_length} beads, box={a.box_size}")
    print(f"Quench: T={a.T_equilibrate}->{a.T_quench}, eps_AB={a.lj_eps_AB}, n_steps={a.n_steps}\n")
    n_total=len(a.kappas)*len(a.seeds);k=0
    t0=time.time()
    for kappa in a.kappas:
        for seed in a.seeds:
            k+=1
            cfg=SimpleNamespace(**vars(a))
            cfg.sequence="correlated"
            cfg.kappa=float(kappa)
            cfg.seed=int(seed)
            cfg.out=os.path.join(a.out,scan_id)
            cfg.run_id=f"kappa{kappa:.2f}_s{seed}"
            cached=os.path.join(cfg.out,cfg.run_id,"structure_factor.npz")
            if os.path.exists(cached):
                print(f"  [{k}/{n_total}] kappa={kappa} seed={seed} (cached, skipped)");continue
            t1=time.time()
            rd=execute(cfg)
            print(f"  [{k}/{n_total}] kappa={kappa} seed={seed} -> {rd}  ({time.time()-t1:.1f}s)")
    print(f"\nKappa scan complete in {time.time()-t0:.1f}s -> {os.path.join(a.out,scan_id)}")

if __name__=="__main__":
    main()
