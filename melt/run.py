from __future__ import annotations
import argparse,os,time
import numpy as np
from .model import MeltParams,RunParams
from .sequences import generate_per_chain
from .box import init_chains_in_box
from .integrator import (build_openmm_system,make_langevin_integrator,make_context,run_simulation,
                         kinetic_temperature_kelvin,kinetic_tstar,minimize_energy,
                         tstar_to_kelvin,KB_KJMOLK)
from .io import make_run_dir,write_meta,open_csv,append_row,compute_mean_Rg,save_structure_factors,save_trajectory,save_density_grids
from .density import density_field_A,density_field_B
from .observables import compute_all_observables

def execute(cfg):
    eps_ref=float(cfg.lj_eps_AA)
    T_eq_star=getattr(cfg,"T_equilibrate",None)
    if T_eq_star is None:T_eq_star=cfg.temperature
    T_q_star=getattr(cfg,"T_quench",None)
    if T_q_star is None:T_q_star=cfg.temperature
    T_eq_K=tstar_to_kelvin(T_eq_star,eps_ref)
    T_q_K=tstar_to_kelvin(T_q_star,eps_ref)
    mp=MeltParams(n_chains=cfg.n_chains,chain_length=cfg.chain_length,box_size=cfg.box_size,
                  bond_k=cfg.bond_k,bond_r0=cfg.bond_r0,
                  lj_eps_AA=cfg.lj_eps_AA,lj_eps_BB=cfg.lj_eps_BB,lj_eps_AB=cfg.lj_eps_AB,
                  lj_eps_core=getattr(cfg,"lj_eps_core",1.0),
                  lj_sigma=cfg.lj_sigma,lj_cutoff=cfg.lj_cutoff,
                  temperature=T_q_K,friction=cfg.friction,dt=cfg.dt)
    rp=RunParams(n_steps=cfg.n_steps,equilibration_steps=cfg.equilibration,
                 snapshot_interval=cfg.snapshot_interval,seed=cfg.seed)
    run_id=cfg.run_id or f"run_{int(time.time())}_{cfg.sequence}_s{cfg.seed}"
    rd=make_run_dir(cfg.out,run_id)
    if os.path.exists(os.path.join(rd,"structure_factor.npz")) and getattr(cfg,"skip_if_cached",True):
        print(f"  cached: {rd} (skipping)")
        return rd
    rng=np.random.default_rng(cfg.seed)
    types=generate_per_chain(cfg.sequence,cfg.n_chains,cfg.chain_length,cfg.f_A,
                             cfg.block_length,cfg.kappa,cfg.pi,rng)
    N=cfg.n_chains*cfg.chain_length
    pos=init_chains_in_box(cfg.n_chains,cfg.chain_length,cfg.box_size,cfg.bond_r0,rng)
    write_meta(rd,mp,rp,cfg.sequence,types,extra={"T_equilibrate_kelvin":float(T_eq_K),
                                                  "T_quench_kelvin":float(T_q_K),
                                                  "T_equilibrate_star":float(T_eq_star),
                                                  "T_quench_star":float(T_q_star),
                                                  "eps_ref_kjmol":eps_ref,
                                                  "k_B_kjmolK":KB_KJMOLK})
    sys=build_openmm_system(types,mp)
    integ=make_langevin_integrator(T_eq_K,cfg.friction,cfg.dt,seed=cfg.seed)
    ctx=make_context(sys,integ,pos,cfg.box_size,platform_name=cfg.platform)
    minimize_energy(ctx)
    if cfg.equilibration>0:
        integ.step(cfg.equilibration)
    if T_q_K!=T_eq_K:
        from openmm import unit
        integ.setTemperature(T_q_K*unit.kelvin)
    f,w=open_csv(rd)
    sf_steps=[];sf_k=None;sf_AA=[];sf_BB=[];sf_AB=[]
    traj_steps=[];traj_pos=[]
    grid_steps=[];grid_phiA=[];grid_phiB=[]
    save_traj=getattr(cfg,"save_trajectory",False)
    save_grids=getattr(cfg,"save_density_grids",False)
    try:
        def cb(step,p_arr,v_arr,pe,ke):
            nonlocal sf_k
            T_inst_K=kinetic_temperature_kelvin(ke,N)
            T_inst_star=kinetic_tstar(ke,N,eps_ref)
            rg=compute_mean_Rg(p_arr,cfg.n_chains,cfg.chain_length)
            kstar=xi=Speak=float("nan")
            if cfg.compute_density:
                phi_A=density_field_A(p_arr,types,cfg.box_size,cfg.grid_size)
                phi_B=density_field_B(p_arr,types,cfg.box_size,cfg.grid_size)
                obs=compute_all_observables(phi_A,phi_B,cfg.box_size)
                kstar=obs["k_star"];xi=obs["xi_AA"];Speak=obs["S_AA_peak"]
                if sf_k is None:
                    sf_k=obs["k"]
                sf_steps.append(int(step))
                sf_AA.append(obs["S_AA"]);sf_BB.append(obs["S_BB"]);sf_AB.append(obs["S_AB"])
                if save_grids:
                    grid_steps.append(int(step));grid_phiA.append(phi_A);grid_phiB.append(phi_B)
            if save_traj:
                traj_steps.append(int(step));traj_pos.append(p_arr.copy())
            append_row(w,run_id,cfg.sequence,mp,cfg.kappa,cfg.pi,cfg.f_A,step,pe,ke,rg,T_inst_star,
                       k_star=kstar,xi_AA=xi,S_AA_peak=Speak)
            f.flush()
        run_simulation(ctx,cfg.n_steps,cfg.snapshot_interval,callback=cb)
    finally:
        f.close()
    if cfg.compute_density and sf_k is not None:
        save_structure_factors(rd,sf_steps,sf_k,sf_AA,sf_BB,sf_AB)
    if save_traj and traj_pos:
        save_trajectory(rd,traj_steps,traj_pos,types,cfg.box_size)
    if save_grids and grid_phiA:
        save_density_grids(rd,grid_steps,grid_phiA,grid_phiB)
    return rd

def parse_args():
    p=argparse.ArgumentParser(description="Multi-chain A/B copolymer melt BD (reduced units)")
    p.add_argument("--out",type=str,default="output/melt")
    p.add_argument("--run_id",type=str,default=None)
    p.add_argument("--seed",type=int,default=12345)
    p.add_argument("--n_chains",type=int,default=10)
    p.add_argument("--chain_length",type=int,default=20)
    p.add_argument("--box_size",type=float,default=10.0)
    p.add_argument("--sequence",type=str,default="random",
                   choices=["random","block","alternating","correlated"])
    p.add_argument("--f_A",type=float,default=0.5)
    p.add_argument("--block_length",type=int,default=4)
    p.add_argument("--kappa",type=float,default=0.7)
    p.add_argument("--pi",type=float,default=0.9)
    p.add_argument("--bond_k",type=float,default=100.0)
    p.add_argument("--bond_r0",type=float,default=1.0)
    p.add_argument("--lj_eps_AA",type=float,default=1.0)
    p.add_argument("--lj_eps_BB",type=float,default=1.0)
    p.add_argument("--lj_eps_AB",type=float,default=0.5)
    p.add_argument("--lj_eps_core",type=float,default=1.0,
                   help="WCA core depth, shared across all pairs (excluded volume)")
    p.add_argument("--lj_sigma",type=float,default=1.0)
    p.add_argument("--lj_cutoff",type=float,default=2.5)
    p.add_argument("--temperature",type=float,default=1.0,help="reduced T* (kB T / eps_AA)")
    p.add_argument("--T_equilibrate",type=float,default=None,help="reduced T*")
    p.add_argument("--T_quench",type=float,default=None,help="reduced T*")
    p.add_argument("--friction",type=float,default=1.0)
    p.add_argument("--dt",type=float,default=0.005)
    p.add_argument("--n_steps",type=int,default=10000)
    p.add_argument("--equilibration",type=int,default=1000)
    p.add_argument("--snapshot_interval",type=int,default=100)
    p.add_argument("--platform",type=str,default=None)
    p.add_argument("--compute_density",action="store_true")
    p.add_argument("--grid_size",type=int,default=32)
    p.add_argument("--save_trajectory",action="store_true")
    p.add_argument("--save_density_grids",action="store_true")
    return p.parse_args()

def main():
    a=parse_args()
    rd=execute(a)
    print(f"Run complete: {rd}")

if __name__=="__main__":
    main()
