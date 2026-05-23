from __future__ import annotations
import argparse,os,time
from types import SimpleNamespace
from .run import execute

def make_cfg(scan_root,scan_id,base,sequence,seed,run_id,**overrides):
    c=SimpleNamespace(**vars(base))
    c.sequence=sequence
    c.seed=int(seed)
    c.out=os.path.join(scan_root,scan_id)
    c.run_id=run_id
    for k,v in overrides.items():setattr(c,k,v)
    return c

def parse_args():
    p=argparse.ArgumentParser(description="Sequence-correlation scans for copolymer melt quench")
    p.add_argument("--kind",type=str,default="sequence",
                   choices=["sequence","pi_kappa","fA_kappa","epsAB"],
                   help="sequence: vary --sequences at fixed (kappa,pi). pi_kappa: grid over --pis x --kappas. fA_kappa: grid over --f_As x --kappas. epsAB: vary --eps_ABs at fixed (kappa,pi).")
    p.add_argument("--out",type=str,default="output/melt/scans")
    p.add_argument("--scan_id",type=str,default=None)
    p.add_argument("--sequences",type=str,nargs="+",
                   default=["random","block","alternating","correlated"])
    p.add_argument("--kappas",type=float,nargs="+",default=None,help="for kind=pi_kappa/fA_kappa")
    p.add_argument("--pis",type=float,nargs="+",default=None,help="for kind=pi_kappa")
    p.add_argument("--f_As",type=float,nargs="+",default=None,help="for kind=fA_kappa")
    p.add_argument("--eps_ABs",type=float,nargs="+",default=None,help="for kind=epsAB")
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
    p.add_argument("--lj_eps_core",type=float,default=1.0)
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

def iter_grid(a):
    if a.kind=="sequence":
        for seq in a.sequences:
            for seed in a.seeds:
                yield make_cfg(a.out,a.scan_id,a,seq,seed,f"{seq}_s{seed}")
    elif a.kind=="pi_kappa":
        assert a.kappas and a.pis,"--kappas and --pis required"
        for k in a.kappas:
            for pi in a.pis:
                for seed in a.seeds:
                    yield make_cfg(a.out,a.scan_id,a,"correlated",seed,
                                   f"kappa{k:.2f}_pi{pi:.3f}_s{seed}",kappa=float(k),pi=float(pi))
    elif a.kind=="fA_kappa":
        assert a.kappas and a.f_As,"--kappas and --f_As required"
        for k in a.kappas:
            for fA in a.f_As:
                for seed in a.seeds:
                    yield make_cfg(a.out,a.scan_id,a,"correlated",seed,
                                   f"kappa{k:.2f}_fA{fA:.2f}_s{seed}",kappa=float(k),f_A=float(fA))
    elif a.kind=="epsAB":
        assert a.eps_ABs,"--eps_ABs required"
        for e in a.eps_ABs:
            for seq in a.sequences:
                for seed in a.seeds:
                    yield make_cfg(a.out,a.scan_id,a,seq,seed,
                                   f"{seq}_eps{e:.5f}_s{seed}",lj_eps_AB=float(e))
    else:
        raise ValueError(f"unknown kind {a.kind}")

def main():
    a=parse_args()
    a.compute_density=True
    scan_id=a.scan_id or f"scan_{a.kind}_{int(time.time())}"
    a.scan_id=scan_id
    print(f"Scan kind={a.kind} id={scan_id}")
    cfgs=list(iter_grid(a))
    print(f"{len(cfgs)} runs total")
    t0=time.time()
    for i,cfg in enumerate(cfgs,1):
        t1=time.time()
        out_check=os.path.join(cfg.out,cfg.run_id,"structure_factor.npz")
        if os.path.exists(out_check):
            print(f"  [{i}/{len(cfgs)}] {cfg.run_id} (cached, skipped)");continue
        rd=execute(cfg)
        print(f"  [{i}/{len(cfgs)}] {cfg.run_id} -> {rd} ({time.time()-t1:.1f}s)")
    print(f"\nScan complete in {time.time()-t0:.1f}s -> {os.path.join(a.out,scan_id)}")

if __name__=="__main__":
    main()
