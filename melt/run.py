from __future__ import annotations
import argparse,os,time
import numpy as np
from .model import MeltParams,RunParams
from .sequences import generate_per_chain,generate_per_chain_exact_total
from .box import init_chains_in_box,relax_overlaps
from .integrator import (build_openmm_system,make_langevin_integrator,make_context,run_simulation,
                         kinetic_temperature_kelvin,kinetic_tstar,minimize_energy,
                         tstar_to_kelvin,KB_KJMOLK)
from .io import make_run_dir,write_meta,open_csv,append_row,compute_mean_Rg,save_structure_factors,save_trajectory,save_density_grids
from .density import density_field_A,density_field_B
from .observables import compute_all_observables
from .direct_structure import DirectStructureFactorRecorder,early_and_final_snapshot_steps
from .modes import ModeAmplitudeRecorder
from .marks import MarkDynamics,MarkRecorder,local_B_counts,step_marks
from .clusters import condensation_summary,CONDENSATION_FIELDS
from .integrator import update_types_in_context

def make_rng_streams(seed,split=False):
    """Return sequence and placement RNGs, splitting only for matched exact-balance studies."""
    if split:
        root=np.random.SeedSequence(int(seed))
        sequence_seed,placement_seed=root.spawn(2)
        return (
            np.random.default_rng(sequence_seed),
            np.random.default_rng(placement_seed),
            {
                "method":"numpy_SeedSequence_spawn",
                "root_seed":int(seed),
                "sequence_spawn_key":[int(x) for x in sequence_seed.spawn_key],
                "placement_spawn_key":[int(x) for x in placement_seed.spawn_key],
                "openmm_integrator_seed":int(seed),
            },
        )
    shared=np.random.default_rng(seed)
    return shared,shared,{"method":"legacy_shared_numpy_stream","root_seed":int(seed)}

def execute(cfg):
    # Reduced temperature unit: the shared WCA core energy, which is 1 for every campaign to
    # date; eps_AA is only a fallback for configurations that do not set a core energy.
    eps_ref=float(getattr(cfg,"lj_eps_core",None) or cfg.lj_eps_AA)
    if eps_ref<=0.0:
        raise ValueError("the reduced temperature reference energy must be positive")
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
    compute_density=bool(getattr(cfg,"compute_density",False))
    compute_direct=bool(getattr(cfg,"compute_direct_structure_factor",False))
    record_modes=bool(getattr(cfg,"record_modes",False))
    mode_interval=int(getattr(cfg,"mode_interval",cfg.snapshot_interval) or cfg.snapshot_interval)
    mode_q_max=float(getattr(cfg,"mode_q_max",1.5))
    if record_modes and (mode_interval<1 or cfg.n_steps%mode_interval!=0):
        raise ValueError("mode_interval must be positive and divide n_steps")
    marks_dynamic=bool(getattr(cfg,"marks_dynamic",False))
    mark_params=None
    if marks_dynamic:
        mark_params=MarkDynamics(k_off=float(cfg.mark_k_off),k_on=float(cfg.mark_k_on),
                                 k_fb=float(getattr(cfg,"mark_k_fb",0.0)),
                                 r_c=float(getattr(cfg,"mark_r_c",1.5)),
                                 n_half=float(getattr(cfg,"mark_n_half",6.0)),
                                 hill=float(getattr(cfg,"mark_hill",2.0)),
                                 interval_steps=int(getattr(cfg,"mark_interval",200)))
        if cfg.n_steps%mark_params.interval_steps!=0:
            raise ValueError("mark_interval must divide n_steps")
    condensation_interval=int(getattr(cfg,"condensation_interval",0) or 0)
    if condensation_interval<0 or (condensation_interval and cfg.n_steps%condensation_interval!=0):
        raise ValueError("condensation_interval must be nonnegative and divide n_steps")
    if condensation_interval and condensation_interval%int(cfg.snapshot_interval)!=0:
        raise ValueError("condensation_interval must be a multiple of snapshot_interval")
    required_cache=[]
    if compute_density:
        required_cache.append(os.path.join(rd,"structure_factor.npz"))
    if compute_direct:
        required_cache.append(os.path.join(rd,"direct_structure_factor.npz"))
    if record_modes:
        required_cache.append(os.path.join(rd,"mode_amplitudes.npz"))
    if marks_dynamic:
        required_cache.append(os.path.join(rd,"marks.npz"))
    if condensation_interval:
        required_cache.append(os.path.join(rd,"condensation.csv"))
    if required_cache and all(os.path.exists(path) for path in required_cache) and getattr(cfg,"skip_if_cached",True):
        print(f"  cached: {rd} (skipping)")
        return rd
    exact_global_composition=bool(getattr(cfg,"exact_global_composition",False))
    sequence_rng,placement_rng,random_streams=make_rng_streams(
        cfg.seed,split=exact_global_composition
    )
    if exact_global_composition:
        balance_max_attempts=int(getattr(cfg,"sequence_balance_max_attempts",100000))
        types,balance_attempts=generate_per_chain_exact_total(
            cfg.sequence,cfg.n_chains,cfg.chain_length,cfg.f_A,
            cfg.block_length,cfg.kappa,cfg.pi,sequence_rng,balance_max_attempts
        )
        sequence_ensemble={
            "method":"independent_per_chain_rejection_conditioned_on_exact_global_count",
            "n_chains_drawn_independently_per_attempt":int(cfg.n_chains),
            "acceptance_attempts":int(balance_attempts),
            "max_attempts":balance_max_attempts,
            "target_A_count":int(cfg.n_chains*cfg.chain_length*cfg.f_A),
            "realized_A_count":int(np.sum(types,dtype=np.int64)),
            "exact_global_f_A":float(np.mean(types)),
            "constructed_complement_pairs":False,
        }
    else:
        types=generate_per_chain(cfg.sequence,cfg.n_chains,cfg.chain_length,cfg.f_A,
                                 cfg.block_length,cfg.kappa,cfg.pi,sequence_rng)
        sequence_ensemble={
            "method":"independent_per_chain",
            "n_independent_base_chains":int(cfg.n_chains),
            "exact_global_f_A":None,
        }
    N=cfg.n_chains*cfg.chain_length
    pos=init_chains_in_box(
        cfg.n_chains,cfg.chain_length,cfg.box_size,cfg.bond_r0,placement_rng
    )
    pos,overlap_info=relax_overlaps(pos,cfg.box_size,cfg.chain_length)
    direct_recorder=None
    direct_extra={"enabled":False}
    if compute_direct:
        direct_early_steps,direct_final_steps=early_and_final_snapshot_steps(
            cfg.n_steps,
            cfg.snapshot_interval,
            getattr(cfg,"direct_early_frames",0),
            getattr(cfg,"direct_final_frames",5),
        )
        direct_steps=np.concatenate((direct_early_steps,direct_final_steps))
        direct_windows=np.asarray(
            ["early"]*len(direct_early_steps)+["final"]*len(direct_final_steps)
        )
        direct_recorder=DirectStructureFactorRecorder(
            cfg.box_size,
            getattr(cfg,"direct_q_max",1.5),
            direct_steps,
            chunk_size=getattr(cfg,"direct_chunk_size",128),
            window_labels=direct_windows,
        )
        direct_extra={
            "enabled":True,
            "q_max":float(getattr(cfg,"direct_q_max",1.5)),
            "final_frames":int(getattr(cfg,"direct_final_frames",5)),
            "early_frames":int(getattr(cfg,"direct_early_frames",0)),
            "early_steps":[int(x) for x in direct_early_steps],
            "final_steps":[int(x) for x in direct_final_steps],
            "n_modes":int(len(direct_recorder.modes.q)),
            "n_shells":int(len(direct_recorder.modes.shell_q)),
            "q_units":"nm^-1",
            "normalization":"S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total",
            "composition_channel":"S_psi_psi^(N)=S_CC=S_AA+S_BB-2*S_AB",
            "primary_channel":"S_psi_psi^(N)/2=S_CC/2",
        }
    mode_recorder=None
    mode_extra={"enabled":False}
    if record_modes:
        mode_recorder=ModeAmplitudeRecorder(cfg.box_size,mode_q_max,
                                            dtype=getattr(cfg,"mode_dtype","complex64"),
                                            allow_type_changes=marks_dynamic)
        mode_extra={
            "enabled":True,
            "interval_steps":mode_interval,
            "n_frames":int(cfg.n_steps//mode_interval),
            "q_max":mode_q_max,
            "n_modes":int(len(mode_recorder.modes.q)),
            "n_shells":int(len(mode_recorder.modes.shell_q)),
            "dtype":str(mode_recorder.dtype),
            "q_units":"nm^-1",
            "amplitude_definition":"rho_alpha(q)=sum_{j in alpha} exp(i q.r_j)",
            "normalization":"S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total",
        }
    write_meta(rd,mp,rp,cfg.sequence,types,extra={"T_equilibrate_kelvin":float(T_eq_K),
                                                  "T_quench_kelvin":float(T_q_K),
                                                  "T_equilibrate_star":float(T_eq_star),
                                                  "T_quench_star":float(T_q_star),
                                                  "eps_ref_kjmol":eps_ref,
                                                  "k_B_kjmolK":KB_KJMOLK,
                                                  "sequence_ensemble":sequence_ensemble,
                                                  "sequence_parameters":{
                                                      "f_A":float(cfg.f_A),
                                                      "kappa":float(cfg.kappa),
                                                      "pi":float(cfg.pi),
                                                      "block_length":int(cfg.block_length),
                                                  },
                                                  "density_analysis":{
                                                      "enabled":compute_density,
                                                      "grid_size":int(cfg.grid_size) if compute_density else None,
                                                  },
                                                  "requested_openmm_platform":cfg.platform,
                                                  "random_streams":random_streams,
                                                  "initial_overlap_relaxation":overlap_info,
                                                  "direct_structure_factor":direct_extra,
                                                  "mode_amplitudes":mode_extra,
                                                  "marks":({"dynamic":True,**mark_params.as_dict(),
                                                            "mark_rng_seed":int(cfg.seed)+7919}
                                                           if marks_dynamic else {"dynamic":False}),
                                                  "condensation":{"enabled":bool(condensation_interval),
                                                                  "interval_steps":condensation_interval,
                                                                  "r_c":float(getattr(cfg,"condensation_r_c",1.5)),
                                                                  "n_dense":int(getattr(cfg,"condensation_n_dense",6))}})
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
    state={"types":np.asarray(types).astype(np.int8).copy()}
    mark_recorder=MarkRecorder(mark_params,N) if marks_dynamic else None
    mark_rng=np.random.default_rng(int(cfg.seed)+7919) if marks_dynamic else None
    cond_rows=[]
    cond_r_c=float(getattr(cfg,"condensation_r_c",1.5));cond_n_dense=int(getattr(cfg,"condensation_n_dense",6))
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
                phi_A=density_field_A(p_arr,state["types"],cfg.box_size,cfg.grid_size)
                phi_B=density_field_B(p_arr,state["types"],cfg.box_size,cfg.grid_size)
                obs=compute_all_observables(phi_A,phi_B,cfg.box_size)
                kstar=obs["k_star"];xi=obs["xi_AA"];Speak=obs["S_AA_peak"]
                if sf_k is None:
                    sf_k=obs["k"]
                sf_steps.append(int(step))
                sf_AA.append(obs["S_AA"]);sf_BB.append(obs["S_BB"]);sf_AB.append(obs["S_AB"])
                if save_grids:
                    grid_steps.append(int(step));grid_phiA.append(phi_A);grid_phiB.append(phi_B)
            if direct_recorder is not None:
                direct_recorder.observe(step,p_arr,state["types"])
            if condensation_interval and step%condensation_interval==0:
                row=condensation_summary(p_arr,state["types"],cfg.box_size,r_c=cond_r_c,n_dense=cond_n_dense)
                cond_rows.append({"step":int(step),**row})
            if save_traj:
                traj_steps.append(int(step));traj_pos.append(p_arr.copy())
            append_row(w,run_id,cfg.sequence,mp,cfg.kappa,cfg.pi,cfg.f_A,step,pe,ke,rg,T_inst_star,
                       k_star=kstar,xi_AA=xi,S_AA_peak=Speak)
            f.flush()
        mode_cb=None
        if mode_recorder is not None:
            def mode_cb(step,p_arr):
                mode_recorder.observe(step,p_arr,state["types"])
        mark_cb=None
        if marks_dynamic:
            dt_marks=float(cfg.dt)*mark_params.interval_steps
            def mark_cb(step,p_arr):
                n_B=local_B_counts(p_arr,state["types"],cfg.box_size,mark_params.r_c)
                new_types,n_on,n_off=step_marks(state["types"],n_B,mark_params,dt_marks,mark_rng)
                mark_recorder.observe(step,new_types,n_on,n_off,n_B)
                if n_on or n_off:
                    state["types"]=new_types
                    update_types_in_context(sys,ctx,new_types)
        run_simulation(ctx,cfg.n_steps,cfg.snapshot_interval,callback=cb,
                       mode_interval=mode_interval if mode_recorder is not None else None,
                       mode_callback=mode_cb,
                       mark_interval=mark_params.interval_steps if marks_dynamic else None,
                       mark_callback=mark_cb)
    finally:
        f.close()
    if cfg.compute_density and sf_k is not None:
        save_structure_factors(rd,sf_steps,sf_k,sf_AA,sf_BB,sf_AB)
    if direct_recorder is not None:
        direct_recorder.save(os.path.join(rd,"direct_structure_factor.npz"))
    if mode_recorder is not None:
        mode_recorder.save(os.path.join(rd,"mode_amplitudes.npz"),
                           mode_interval_steps=mode_interval,
                           production_steps=int(cfg.n_steps),
                           dt_ps=float(cfg.dt),
                           kappa=float(cfg.kappa),
                           pi=float(cfg.pi),
                           f_A=float(cfg.f_A),
                           lj_eps_AB=float(cfg.lj_eps_AB),
                           T_quench_star=float(T_q_star),
                           seed=int(cfg.seed))
    if mark_recorder is not None:
        mark_recorder.save(os.path.join(rd,"marks.npz"),production_steps=int(cfg.n_steps),
                           dt_ps=float(cfg.dt),initial_types=np.asarray(types).astype(np.int8),
                           seed=int(cfg.seed),T_quench_star=float(T_q_star))
    if cond_rows:
        import csv as _csv
        with open(os.path.join(rd,"condensation.csv"),"w",newline="") as ch:
            cw=_csv.DictWriter(ch,fieldnames=list(CONDENSATION_FIELDS),lineterminator="\n")
            cw.writeheader()
            for row in cond_rows:
                cw.writerow({k:row.get(k,"") for k in CONDENSATION_FIELDS})
    if save_traj and traj_pos:
        save_trajectory(rd,traj_steps,traj_pos,state["types"],cfg.box_size)
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
    p.add_argument("--exact_global_composition",action="store_true",
                   help="rejection-sample independent chain sets until the global A count is exact")
    p.add_argument("--sequence_balance_max_attempts",type=int,default=100000,
                   help="maximum full-set draws for --exact_global_composition")
    p.add_argument("--compute_direct_structure_factor",action="store_true",
                   help="compute particle-level partial spectra at reciprocal-box modes")
    p.add_argument("--direct_q_max",type=float,default=1.5,
                   help="largest |q| included in direct particle-level spectra")
    p.add_argument("--direct_final_frames",type=int,default=5,
                   help="number of final production snapshots used for direct spectra")
    p.add_argument("--direct_early_frames",type=int,default=0,
                   help="number of preceding snapshots stored as a convergence block")
    p.add_argument("--direct_chunk_size",type=int,default=128,
                   help="number of q vectors evaluated per memory-bounded chunk")
    p.add_argument("--record_modes",action="store_true",
                   help="store exact rho_A(q,t), rho_B(q,t) at every box mode for the dynamic structure factor")
    p.add_argument("--mode_interval",type=int,default=None,
                   help="steps between recorded mode frames (default: snapshot_interval)")
    p.add_argument("--mode_q_max",type=float,default=1.5,
                   help="largest |q| recorded in mode_amplitudes.npz")
    p.add_argument("--mode_dtype",type=str,default="complex64",choices=["complex64","complex128"])
    p.add_argument("--marks_dynamic",action="store_true",
                   help="let marks turn over (k_off) and be written with reader-writer feedback during the run")
    p.add_argument("--mark_k_off",type=float,default=1e-3,help="B->A turnover rate per tau")
    p.add_argument("--mark_k_on",type=float,default=1e-3,help="basal A->B writing rate per tau")
    p.add_argument("--mark_k_fb",type=float,default=0.0,help="feedback gain: extra A->B rate at saturating local B count")
    p.add_argument("--mark_r_c",type=float,default=1.5,help="neighbor radius for the local B count")
    p.add_argument("--mark_n_half",type=float,default=6.0,help="Hill half-saturation B-neighbor count")
    p.add_argument("--mark_hill",type=float,default=2.0)
    p.add_argument("--mark_interval",type=int,default=200,help="steps between mark updates")
    p.add_argument("--condensation_interval",type=int,default=0,
                   help="steps between condensation summaries (0 disables; multiple of snapshot_interval)")
    p.add_argument("--condensation_r_c",type=float,default=1.5)
    p.add_argument("--condensation_n_dense",type=int,default=6)
    return p.parse_args()

def main():
    a=parse_args()
    rd=execute(a)
    print(f"Run complete: {rd}")

if __name__=="__main__":
    main()
