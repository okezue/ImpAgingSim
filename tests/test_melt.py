from __future__ import annotations
import math
import numpy as np
import pytest
from melt.model import MeltParams,RunParams,lj_pair_energy
from melt.sequences import (generate_random,generate_block,generate_alternating,
                            generate_correlated,generate_per_chain,valid_pi_range,autocorrelation)
from melt.box import init_chains_in_box,apply_pbc,bonded_pairs,chain_index
from melt.density import grid_density_cic,density_field_A,density_field_B,density_contrast
from melt.observables import (structure_factor_3d,cross_structure_factor_3d,k_grid,
                              spherical_average,find_peak,domain_length,fit_ornstein_zernike,
                              compute_all_observables)
from melt.integrator import tstar_to_kelvin,kelvin_to_tstar,KB_KJMOLK
try:
    import openmm
    HAS_OPENMM=True
except ImportError:
    HAS_OPENMM=False
needs_openmm=pytest.mark.skipif(not HAS_OPENMM,reason="openmm not installed")

class TestSequences:
    def test_random_composition(self):
        rng=np.random.default_rng(0)
        s=generate_random(20000,0.3,rng)
        assert abs(s.mean()-0.3)<0.01
        assert set(np.unique(s)).issubset({0,1})
    def test_alternating(self):
        s=generate_alternating(10)
        assert list(s)==[1,0,1,0,1,0,1,0,1,0]
    def test_block_period(self):
        s=generate_block(20,4,0.5)
        assert s[0]==1 and s[3]==1 and s[4]==0 and s[7]==0
        assert s[8]==1
    def test_correlated_lambda_2pi_minus_1_at_symmetric(self):
        rng=np.random.default_rng(3)
        pi=0.85
        s=generate_correlated(40000,kappa=1.0,pi=pi,f_A=0.5,rng=rng)
        ac=autocorrelation(s,kmax=3)
        assert abs(ac[1]-(2*pi-1))<0.02
    def test_correlated_lambda_invariant_in_fA(self):
        rng=np.random.default_rng(7)
        pi=0.9
        for fA in (0.3,0.5,0.7):
            s=generate_correlated(40000,kappa=1.0,pi=pi,f_A=fA,rng=rng)
            ac=autocorrelation(s,kmax=2)
            sm=ac[0]
            lag1_norm=(ac[1]-(2*fA-1)**2)/(sm-(2*fA-1)**2) if abs(sm-(2*fA-1)**2)>1e-6 else 0
            assert abs(lag1_norm-(2*pi-1))<0.03,f"fA={fA} got {lag1_norm}"
    def test_correlated_composition_matches_fA(self):
        rng=np.random.default_rng(11)
        for fA in (0.2,0.5,0.7):
            s=generate_correlated(20000,kappa=1.0,pi=0.9,f_A=fA,rng=rng)
            assert abs(s.mean()-fA)<0.02
    def test_correlated_invalid_pi_raises(self):
        rng=np.random.default_rng(0)
        with pytest.raises(ValueError):
            generate_correlated(100,kappa=1.0,pi=0.1,f_A=0.2,rng=rng)
    def test_valid_pi_range(self):
        pmin,pmax=valid_pi_range(0.5)
        assert pmin==0.0 and pmax==1.0
        pmin,_=valid_pi_range(0.2)
        assert abs(pmin-0.375)<1e-9
    def test_per_chain_no_cross_chain_correlation(self):
        rng=np.random.default_rng(13)
        n_chains=200;L=50
        s=generate_per_chain("correlated",n_chains,L,f_A=0.5,block_length=4,
                             kappa=1.0,pi=0.95,rng=rng)
        boundary_corr=float(np.mean([(2*int(s[(c+1)*L])-1)*(2*int(s[(c+1)*L-1])-1)
                                     for c in range(n_chains-1)]))
        assert abs(boundary_corr)<0.1,f"cross-chain correlation = {boundary_corr}"
    def test_per_chain_within_chain_correlation_preserved(self):
        rng=np.random.default_rng(17)
        pi=0.9
        s=generate_per_chain("correlated",200,100,f_A=0.5,block_length=4,
                             kappa=1.0,pi=pi,rng=rng)
        ac=autocorrelation(s,kmax=3,per_chain_length=100)
        assert abs(ac[1]-(2*pi-1))<0.03

class TestTemperatureUnits:
    def test_tstar_to_kelvin_at_unit_eps(self):
        assert abs(tstar_to_kelvin(0.7,1.0)-0.7/KB_KJMOLK)<1e-6
    def test_tstar_to_kelvin_round_trip(self):
        for tstar in (0.1,0.7,1.0,5.0):
            for eps in (0.5,1.0,2.0):
                Tk=tstar_to_kelvin(tstar,eps)
                assert abs(kelvin_to_tstar(Tk,eps)-tstar)<1e-9
    def test_tstar_0_7_gives_84_19_K(self):
        assert abs(tstar_to_kelvin(0.7,1.0)-84.187)<0.05

class TestBox:
    def test_init_shape(self):
        rng=np.random.default_rng(0)
        pos=init_chains_in_box(5,10,10.0,1.0,rng)
        assert pos.shape==(50,3)
    def test_bonded_pairs_count(self):
        assert len(bonded_pairs(5,10))==5*9
    def test_bonded_pairs_no_inter_chain(self):
        for i,j in bonded_pairs(3,10):
            assert chain_index(i,10)==chain_index(j,10)
    def test_pbc_wraps(self):
        L=10.0
        pos=np.array([[11.5,-1.0,5.0]])
        wp=apply_pbc(pos,L)
        assert abs(wp[0,0]-1.5)<1e-9
        assert abs(wp[0,1]-9.0)<1e-9

class TestEnergyHelper:
    def test_lj_min_at_2to1over6_sigma(self):
        rmin=2**(1/6)
        assert abs(lj_pair_energy(rmin,1.0,1.0)-(-1.0))<1e-9
    def test_lj_zero_at_sigma(self):
        assert abs(lj_pair_energy(1.0,1.0,1.0))<1e-9

@needs_openmm
class TestForceFieldSeparation:
    """Verify the WCA core is shared across AA, BB, AB pairs and the attractive tail is
    separately tunable. The previous single-LJ form failed this test because lowering
    lj_eps_AB also weakened the cross-pair excluded volume."""
    def _pair_energy(self,types,r,eps_AB,eps_core=1.0):
        from melt.integrator import (build_openmm_system,make_langevin_integrator,
                                     make_context,get_state_arrays)
        mp=MeltParams(n_chains=2,chain_length=1,box_size=20.0,bond_k=0.0,bond_r0=1.0,
                      lj_eps_AA=1.0,lj_eps_BB=1.0,lj_eps_AB=eps_AB,lj_eps_core=eps_core,
                      lj_sigma=1.0,lj_cutoff=2.5)
        pos=np.array([[5.0,5.0,5.0],[5.0+r,5.0,5.0]])
        sys=build_openmm_system(np.asarray(types,dtype=np.int8),mp)
        integ=make_langevin_integrator(1.0,1.0,0.001,seed=0)
        ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
        _,_,pe,_=get_state_arrays(ctx)
        return float(pe)
    def test_wca_core_shared_across_pair_types(self):
        r=0.95
        eAA=self._pair_energy([1,1],r,eps_AB=0.1)
        eAB=self._pair_energy([1,0],r,eps_AB=0.1)
        eBB=self._pair_energy([0,0],r,eps_AB=0.1)
        assert abs(eAA-eAB)<1e-3 and abs(eAA-eBB)<1e-3
    def test_wca_core_unaffected_by_eps_AB(self):
        r=0.95
        e_eAB_large=self._pair_energy([1,0],r,eps_AB=0.5)
        e_eAB_tiny=self._pair_energy([1,0],r,eps_AB=1e-4)
        assert abs(e_eAB_large-e_eAB_tiny)<1e-3
    def test_attractive_tail_scales_with_eps_AB(self):
        r=1.5
        e1=self._pair_energy([1,0],r,eps_AB=1.0)
        e2=self._pair_energy([1,0],r,eps_AB=0.1)
        assert e2/e1<0.5

@needs_openmm
class TestIntegrator:
    def test_build_system(self):
        from melt.integrator import build_openmm_system
        mp=MeltParams(n_chains=2,chain_length=5,box_size=10.0)
        types=np.array([1,0,1,0,1,1,0,1,0,1],dtype=np.int8)
        s=build_openmm_system(types,mp)
        assert s.getNumParticles()==10
        assert s.getNumForces()==3
    def test_kinetic_temperature_matches_reduced_target(self):
        from melt.integrator import (build_openmm_system,make_langevin_integrator,make_context,
                                     minimize_energy,kinetic_tstar,get_state_arrays,tstar_to_kelvin)
        T_star=1.0
        T_K=tstar_to_kelvin(T_star,1.0)
        mp=MeltParams(n_chains=5,chain_length=10,box_size=20.0,temperature=T_K,
                      friction=1.0,dt=0.005)
        rng=np.random.default_rng(0)
        types=generate_random(50,0.5,rng)
        pos=init_chains_in_box(5,10,20.0,1.0,rng)
        sys=build_openmm_system(types,mp)
        integ=make_langevin_integrator(T_K,1.0,0.005,seed=0)
        ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
        minimize_energy(ctx)
        integ.step(3000)
        Ts=[]
        for _ in range(40):
            integ.step(200)
            _,_,_,ke=get_state_arrays(ctx)
            Ts.append(kinetic_tstar(ke,50,1.0))
        assert abs(float(np.mean(Ts))-T_star)<0.25

class TestDensity:
    def test_cic_mass_conservation(self):
        rng=np.random.default_rng(0)
        N=500;L=10.0
        pos=rng.uniform(0,L,size=(N,3))
        g=grid_density_cic(pos,np.ones(N),L,16)
        assert abs(g.sum()-N)<1e-9
    def test_AB_partition(self):
        rng=np.random.default_rng(0)
        N=200;L=10.0
        pos=rng.uniform(0,L,size=(N,3))
        types=(rng.random(N)<0.5).astype(np.int8)
        gA=density_field_A(pos,types,L,16)
        gB=density_field_B(pos,types,L,16)
        assert abs((gA+gB).sum()-N)<1e-9

class TestStructureFactor:
    def test_sinusoidal_pattern_peaks_at_known_k(self):
        L=10.0;G=32;n=4
        x=np.arange(G)*L/G
        X,Y,Z=np.meshgrid(x,x,x,indexing="ij")
        kw=2.0*np.pi*n/L
        phi=np.cos(kw*X)
        S3=structure_factor_3d(phi,L)
        k,Sk=spherical_average(S3,L,n_bins=G)
        kstar,_=find_peak(k,Sk,k_min=0.5*kw)
        assert abs(kstar-kw)/kw<0.1
    def test_domain_length_inverts_k_star(self):
        kstar=2.0*np.pi/3.5
        assert abs(domain_length(kstar)-3.5)<1e-9
    def test_oz_fit_recovers_xi(self):
        xi_true=2.0;S0=10.0
        k=np.linspace(0.05,1.5,30)
        S=S0/(1.0+(k*xi_true)**2)
        xi_fit,S0_fit=fit_ornstein_zernike(k,S,k_max=1.5)
        assert abs(xi_fit-xi_true)<0.05

class TestChi4Normalization:
    """chi_4 = N * Var_runs(Q) requires N = system size, not n_runs. Regression test
    for the patch that fixed a ~1920x underestimate on the corrected aging campaign."""
    def test_chi4_uses_n_particles(self):
        from melt.dynamics import chi4_from_overlap_realizations
        rng=np.random.default_rng(0)
        n_runs=3;n_lags=20;N=5760
        Q=rng.uniform(0.2,0.8,size=(n_runs,n_lags))
        chi4=chi4_from_overlap_realizations(Q,n_particles=N)
        expected=N*np.var(Q,axis=0,ddof=1)
        np.testing.assert_allclose(chi4,expected,rtol=1e-12)
    def test_chi4_rejects_nonpositive_N(self):
        from melt.dynamics import chi4_from_overlap_realizations
        Q=np.array([[0.5,0.4],[0.3,0.2]])
        with pytest.raises(ValueError):
            chi4_from_overlap_realizations(Q,n_particles=0)
    def test_chi4_single_run_returns_nan(self):
        from melt.dynamics import chi4_from_overlap_realizations
        Q=np.array([[0.5,0.4,0.3]])
        out=chi4_from_overlap_realizations(Q,n_particles=100)
        assert np.all(np.isnan(out))

@needs_openmm
class TestReproducibility:
    def test_same_seed_same_trajectory(self):
        from melt.integrator import (build_openmm_system,make_langevin_integrator,
                                     make_context,get_state_arrays,tstar_to_kelvin)
        def run_one(seed):
            mp=MeltParams(n_chains=3,chain_length=8,box_size=10.0,dt=0.005,
                          temperature=tstar_to_kelvin(1.0,1.0))
            rng=np.random.default_rng(seed)
            types=generate_random(24,0.5,rng)
            pos=init_chains_in_box(3,8,10.0,1.0,rng)
            sys=build_openmm_system(types,mp)
            integ=make_langevin_integrator(mp.temperature,1.0,0.005,seed=seed)
            ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
            integ.step(500)
            p,_,_,_=get_state_arrays(ctx)
            return p
        a=run_one(7);b=run_one(7)
        np.testing.assert_allclose(a,b,atol=1e-5)
