import sys, time
import numpy as np
sys.path.insert(0,"/Users/okezuebell/Documents/GitHub/ImpAgingSim")
from melt.model import MeltParams
from melt.sequences import generate_correlated
from melt.box import init_chains_in_box
from melt.integrator import (build_openmm_system, make_langevin_integrator,
                             make_context, minimize_energy, get_state_arrays)
from melt.density import density_field_A, density_field_B

n_chains=144; N=40; L=17.0
dt=0.005; friction=1.0
T_eq=5.0; T_q=0.7
n_eq=30000; n_prod=250000
snap_eq=600; snap_prod=2500
grid_size=56
pi=0.99
seed=1

mp=MeltParams(n_chains=n_chains, chain_length=N, box_size=L,
              bond_k=100.0, bond_r0=1.0,
              lj_eps_AA=1.0, lj_eps_BB=1.0, lj_eps_AB=0.1,
              lj_sigma=1.0, lj_cutoff=2.5,
              temperature=T_q, friction=friction, dt=dt)

platform=None
try:
    import openmm
    from openmm import Platform
    for nm in ("OpenCL","CUDA","CPU"):
        try:
            Platform.getPlatformByName(nm); platform=nm; break
        except Exception: continue
except Exception: platform="CPU"
print(f"platform: {platform}", flush=True)

out={}
for kap in (0.0, 0.5, 1.0):
    print(f"=== kappa={kap} ===", flush=True)
    rng=np.random.default_rng(seed)
    types=generate_correlated(n_chains*N, kap, pi, rng)
    pos=init_chains_in_box(n_chains, N, L, 1.0, rng)
    sys_=build_openmm_system(types, mp)
    integ=make_langevin_integrator(T_eq, friction, dt, seed=seed)
    ctx=make_context(sys_, integ, pos, L, platform_name=platform)
    minimize_energy(ctx)
    steps_log=[]; phiA_log=[]; phiB_log=[]; phase_log=[]
    p_arr,_,_,_=get_state_arrays(ctx)
    phiA_log.append(density_field_A(p_arr, types, L, grid_size))
    phiB_log.append(density_field_B(p_arr, types, L, grid_size))
    steps_log.append(0); phase_log.append("min")
    t0=time.time()
    for s in range(n_eq//snap_eq):
        integ.step(snap_eq)
        p_arr,_,_,_=get_state_arrays(ctx)
        phiA_log.append(density_field_A(p_arr, types, L, grid_size))
        phiB_log.append(density_field_B(p_arr, types, L, grid_size))
        steps_log.append((s+1)*snap_eq); phase_log.append("eq")
    print(f"  eq done in {time.time()-t0:.1f}s", flush=True)
    integ.setTemperature(T_q)
    t1=time.time()
    for s in range(n_prod//snap_prod):
        integ.step(snap_prod)
        p_arr,_,_,_=get_state_arrays(ctx)
        phiA_log.append(density_field_A(p_arr, types, L, grid_size))
        phiB_log.append(density_field_B(p_arr, types, L, grid_size))
        steps_log.append(n_eq+(s+1)*snap_prod); phase_log.append("prod")
    print(f"  prod done in {time.time()-t1:.1f}s", flush=True)
    out[kap]={"steps":np.asarray(steps_log),"phase":np.asarray(phase_log),
              "phi_A":np.asarray(phiA_log,dtype=np.float32),
              "phi_B":np.asarray(phiB_log,dtype=np.float32)}

np.savez_compressed("/tmp/prequench_grids_v2.npz",
    k0_steps=out[0.0]["steps"], k0_phase=out[0.0]["phase"],
    k0_phi_A=out[0.0]["phi_A"], k0_phi_B=out[0.0]["phi_B"],
    k05_steps=out[0.5]["steps"], k05_phase=out[0.5]["phase"],
    k05_phi_A=out[0.5]["phi_A"], k05_phi_B=out[0.5]["phi_B"],
    k1_steps=out[1.0]["steps"], k1_phase=out[1.0]["phase"],
    k1_phi_A=out[1.0]["phi_A"], k1_phi_B=out[1.0]["phi_B"],
    n_eq=n_eq, T_eq=T_eq, T_q=T_q, L=L, grid_size=grid_size)
print("saved /tmp/prequench_grids_v2.npz")
