from __future__ import annotations
import numpy as np
try:
    import openmm as mm
    from openmm import unit
    HAS_OPENMM=True
except ImportError:
    HAS_OPENMM=False
from .box import bonded_pairs

def _require_openmm():
    if not HAS_OPENMM:
        raise ImportError("openmm not installed. Run: pip install openmm>=8.0")

def build_openmm_system(types,mp):
    _require_openmm()
    N=len(types)
    sys=mm.System()
    for _ in range(N):
        sys.addParticle(1.0)
    L=mp.box_size
    sys.setDefaultPeriodicBoxVectors(mm.Vec3(L,0,0),mm.Vec3(0,L,0),mm.Vec3(0,0,L))
    bf=mm.HarmonicBondForce()
    for i,j in bonded_pairs(mp.n_chains,mp.chain_length):
        bf.addBond(int(i),int(j),mp.bond_r0,mp.bond_k)
    sys.addForce(bf)
    expr=("4*eps*((sig/r)^12-(sig/r)^6);"
          "eps=epsAA*tA1*tA2+epsBB*(1-tA1)*(1-tA2)+epsAB*(tA1*(1-tA2)+(1-tA1)*tA2);"
          f"epsAA={mp.lj_eps_AA};epsBB={mp.lj_eps_BB};epsAB={mp.lj_eps_AB};"
          f"sig={mp.lj_sigma}")
    nbf=mm.CustomNonbondedForce(expr)
    nbf.addPerParticleParameter("tA")
    nbf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nbf.setCutoffDistance(mp.lj_cutoff*mp.lj_sigma)
    types_py=[1.0 if int(t)==1 else 0.0 for t in types]
    for k in range(N):
        nbf.addParticle([types_py[k]])
    bonds_list=bonded_pairs(mp.n_chains,mp.chain_length)
    nbf.createExclusionsFromBonds([(int(i),int(j)) for i,j in bonds_list],1)
    sys.addForce(nbf)
    return sys

def make_langevin_integrator(temperature,friction,dt,seed=None):
    _require_openmm()
    integ=mm.LangevinMiddleIntegrator(temperature,friction,dt)
    if seed is not None:
        integ.setRandomNumberSeed(int(seed))
    return integ

def minimize_energy(ctx,tolerance=10.0,max_iterations=0):
    _require_openmm()
    mm.LocalEnergyMinimizer.minimize(ctx,tolerance,max_iterations)

def make_context(system,integrator,positions,box_size,platform_name=None):
    _require_openmm()
    if platform_name is None:
        ctx=mm.Context(system,integrator)
    else:
        plat=mm.Platform.getPlatformByName(platform_name)
        ctx=mm.Context(system,integrator,plat)
    ctx.setPeriodicBoxVectors(mm.Vec3(box_size,0,0),mm.Vec3(0,box_size,0),mm.Vec3(0,0,box_size))
    ctx.setPositions(positions.tolist())
    return ctx

def get_state_arrays(ctx):
    _require_openmm()
    st=ctx.getState(getPositions=True,getEnergy=True,getVelocities=True)
    pos=np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    vel=np.asarray(st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond))
    pe=float(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    ke=float(st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
    return pos,vel,pe,ke

def kinetic_temperature(ke,n_particles,k_B=0.0083144626):
    return 2.0*ke/(3.0*n_particles*k_B)

def run_simulation(ctx,n_steps,snapshot_interval,callback=None):
    _require_openmm()
    integ=ctx.getIntegrator()
    n_snap=n_steps//snapshot_interval
    for s in range(n_snap):
        integ.step(snapshot_interval)
        if callback is not None:
            step=(s+1)*snapshot_interval
            pos,vel,pe,ke=get_state_arrays(ctx)
            callback(step,pos,vel,pe,ke)
