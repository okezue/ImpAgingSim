from __future__ import annotations
import numpy as np
try:
    import openmm as mm
    from openmm import unit
    HAS_OPENMM=True
except ImportError:
    HAS_OPENMM=False
from .box import bonded_pairs

KB_KJMOLK=0.00831446261815324
TWO_SIXTH=2.0**(1.0/6.0)

def _require_openmm():
    if not HAS_OPENMM:
        raise ImportError("openmm not installed. Run: pip install openmm>=8.0")

def tstar_to_kelvin(tstar,eps_kjmol=1.0):
    return float(tstar)*float(eps_kjmol)/KB_KJMOLK

def kelvin_to_tstar(tk,eps_kjmol=1.0):
    return float(tk)*KB_KJMOLK/float(eps_kjmol)

def build_openmm_system(types,mp):
    """Build melt system with shared WCA repulsive core + separately tunable attractive tail.

    Repulsive core (all pairs): 4*eps_core*[(sig/r)^12-(sig/r)^6]+eps_core for r<2^(1/6)*sig, else 0.
    Attractive tail (type-dependent): 4*eps_att*[(sig/r)^12-(sig/r)^6] for r>=2^(1/6)*sig, with eps_att
    set by the type pair (eps_AA, eps_BB, eps_AB). Tail truncated and shifted to zero at lj_cutoff.

    This separates excluded volume (universal) from incompatibility (per-pair). The previous
    single-LJ form scaled both core and tail with eps_AB, making the eps_AB=10^-4 scan equivalent
    to removing A-B excluded volume.
    """
    _require_openmm()
    N=len(types)
    sys=mm.System()
    for _ in range(N):
        sys.addParticle(1.0*unit.amu)
    L=float(mp.box_size)
    sys.setDefaultPeriodicBoxVectors(mm.Vec3(L,0,0)*unit.nanometer,
                                     mm.Vec3(0,L,0)*unit.nanometer,
                                     mm.Vec3(0,0,L)*unit.nanometer)
    bf=mm.HarmonicBondForce()
    bonds=bonded_pairs(mp.n_chains,mp.chain_length)
    for i,j in bonds:
        bf.addBond(int(i),int(j),float(mp.bond_r0)*unit.nanometer,
                   float(mp.bond_k)*unit.kilojoule_per_mole/unit.nanometer**2)
    sys.addForce(bf)

    eps_core=float(getattr(mp,"lj_eps_core",1.0))
    sig=float(mp.lj_sigma)
    rc_wca=TWO_SIXTH*sig
    core_expr=(f"step({rc_wca}-r)*(4*{eps_core}*((sig/r)^12-(sig/r)^6)+{eps_core});"
               f"sig={sig}")
    core=mm.CustomNonbondedForce(core_expr)
    core.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    core.setCutoffDistance(rc_wca*unit.nanometer)
    for _ in range(N):
        core.addParticle([])
    core.createExclusionsFromBonds([(int(i),int(j)) for i,j in bonds],1)
    sys.addForce(core)

    rc=float(mp.lj_cutoff)*sig
    shift_AA=4.0*float(mp.lj_eps_AA)*((sig/rc)**12-(sig/rc)**6)
    shift_BB=4.0*float(mp.lj_eps_BB)*((sig/rc)**12-(sig/rc)**6)
    shift_AB=4.0*float(mp.lj_eps_AB)*((sig/rc)**12-(sig/rc)**6)
    att_expr=(f"step(r-{rc_wca})*(4*eps*((sig/r)^12-(sig/r)^6)-shift);"
              "eps=epsAA*tA1*tA2+epsBB*(1-tA1)*(1-tA2)+epsAB*(tA1*(1-tA2)+(1-tA1)*tA2);"
              "shift=shiftAA*tA1*tA2+shiftBB*(1-tA1)*(1-tA2)+shiftAB*(tA1*(1-tA2)+(1-tA1)*tA2);"
              f"epsAA={mp.lj_eps_AA};epsBB={mp.lj_eps_BB};epsAB={mp.lj_eps_AB};"
              f"shiftAA={shift_AA};shiftBB={shift_BB};shiftAB={shift_AB};"
              f"sig={sig}")
    att=mm.CustomNonbondedForce(att_expr)
    att.addPerParticleParameter("tA")
    att.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    att.setCutoffDistance(rc*unit.nanometer)
    types_py=[1.0 if int(t)==1 else 0.0 for t in types]
    for k in range(N):
        att.addParticle([types_py[k]])
    att.createExclusionsFromBonds([(int(i),int(j)) for i,j in bonds],1)
    sys.addForce(att)
    return sys

def make_langevin_integrator(temperature_kelvin,friction,dt,seed=None):
    """temperature_kelvin must be in Kelvin. Use tstar_to_kelvin() to convert from reduced T*."""
    _require_openmm()
    integ=mm.LangevinMiddleIntegrator(float(temperature_kelvin)*unit.kelvin,
                                      float(friction)/unit.picosecond,
                                      float(dt)*unit.picosecond)
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
    L=float(box_size)
    ctx.setPeriodicBoxVectors(mm.Vec3(L,0,0)*unit.nanometer,
                              mm.Vec3(0,L,0)*unit.nanometer,
                              mm.Vec3(0,0,L)*unit.nanometer)
    ctx.setPositions((positions*unit.nanometer).tolist() if isinstance(positions,np.ndarray) else positions.tolist())
    return ctx

def get_state_arrays(ctx):
    _require_openmm()
    st=ctx.getState(getPositions=True,getEnergy=True,getVelocities=True)
    pos=np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    vel=np.asarray(st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond))
    pe=float(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    ke=float(st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
    return pos,vel,pe,ke

def kinetic_temperature_kelvin(ke_kjmol,n_particles):
    return 2.0*ke_kjmol/(3.0*n_particles*KB_KJMOLK)

def kinetic_tstar(ke_kjmol,n_particles,eps_kjmol=1.0):
    return kelvin_to_tstar(kinetic_temperature_kelvin(ke_kjmol,n_particles),eps_kjmol)

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
