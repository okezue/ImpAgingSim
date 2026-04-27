from __future__ import annotations
import argparse,os,time
from types import SimpleNamespace
from .run import execute

def parse_args():
    p=argparse.ArgumentParser(description="Single big production run with full recording")
    p.add_argument("--out",type=str,default="output/melt/big")
    p.add_argument("--run_id",type=str,default=None)
    p.add_argument("--seed",type=int,default=42)
    p.add_argument("--sequence",type=str,default="correlated",
                   choices=["random","block","alternating","correlated"])
    p.add_argument("--n_chains",type=int,default=200)
    p.add_argument("--chain_length",type=int,default=30)
    p.add_argument("--box_size",type=float,default=22.0)
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
    p.add_argument("--n_steps",type=int,default=200000)
    p.add_argument("--equilibration",type=int,default=20000)
    p.add_argument("--snapshot_interval",type=int,default=2000)
    p.add_argument("--platform",type=str,default=None)
    p.add_argument("--grid_size",type=int,default=64)
    p.add_argument("--no_save_trajectory",action="store_true")
    p.add_argument("--no_save_grids",action="store_true")
    return p.parse_args()

def main():
    a=parse_args()
    a.compute_density=True
    a.save_trajectory=not a.no_save_trajectory
    a.save_density_grids=not a.no_save_grids
    a.run_id=a.run_id or f"big_{a.sequence}_{int(time.time())}"
    print(f"Big run: id={a.run_id}")
    print(f"System: {a.n_chains}x{a.chain_length}={a.n_chains*a.chain_length} beads, box={a.box_size}")
    print(f"Quench: T={a.T_equilibrate}->{a.T_quench}, eps_AB={a.lj_eps_AB}")
    print(f"Steps: {a.equilibration} eq + {a.n_steps} prod, snapshot every {a.snapshot_interval}")
    print(f"Save: traj={a.save_trajectory}, grids={a.save_density_grids}\n")
    t0=time.time()
    rd=execute(a)
    print(f"\nBig run complete in {time.time()-t0:.1f}s -> {rd}")

if __name__=="__main__":
    main()
