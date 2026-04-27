from __future__ import annotations
import math
import numpy as np
import pytest
from melt.model import MeltParams,RunParams,lj_pair_energy
from melt.sequences import generate_random,generate_block,generate_alternating,generate_correlated,autocorrelation
from melt.box import init_chains_in_box,apply_pbc,bonded_pairs,chain_index
from melt.density import grid_density_cic,density_field_A,density_field_B,density_contrast
from melt.observables import structure_factor_3d,cross_structure_factor_3d,k_grid,spherical_average,find_peak,domain_length,fit_ornstein_zernike,compute_all_observables
try:
    import openmm
    HAS_OPENMM=True
except ImportError:
    HAS_OPENMM=False
needs_openmm=pytest.mark.skipif(not HAS_OPENMM,reason="openmm not installed")

class TestSequences:
    def test_random_composition(self):
        rng=np.random.default_rng(0)
        N=20000
        s=generate_random(N,0.3,rng)
        assert abs(s.mean()-0.3)<0.01
        assert set(np.unique(s)).issubset({0,1})
    def test_alternating(self):
        s=generate_alternating(10)
        assert list(s)==[1,0,1,0,1,0,1,0,1,0]
    def test_block_period(self):
        s=generate_block(20,4,0.5)
        assert s[0]==1 and s[3]==1 and s[4]==0 and s[7]==0
        assert s[8]==1
    def test_correlated_extreme_kappa1_pi1(self):
        rng=np.random.default_rng(1)
        s=generate_correlated(50,kappa=1.0,pi=1.0,rng=rng)
        assert len(set(s))==1
    def test_correlated_pi05_uncorrelated(self):
        rng=np.random.default_rng(2)
        N=20000
        s=generate_correlated(N,kappa=1.0,pi=0.5,rng=rng)
        ac=autocorrelation(s,kmax=5)
        assert abs(ac[1])<0.05
    def test_correlated_persistence_decay(self):
        rng=np.random.default_rng(3)
        N=20000
        pi=0.85
        s=generate_correlated(N,kappa=1.0,pi=pi,rng=rng)
        ac=autocorrelation(s,kmax=5)
        expect=(2*pi-1)**1
        assert abs(ac[1]-expect)<0.05
    def test_reproducibility(self):
        rng1=np.random.default_rng(42)
        rng2=np.random.default_rng(42)
        a=generate_random(100,0.5,rng1)
        b=generate_random(100,0.5,rng2)
        np.testing.assert_array_equal(a,b)

class TestBox:
    def test_init_shape(self):
        rng=np.random.default_rng(0)
        pos=init_chains_in_box(5,10,L:=10.0,1.0,rng)
        assert pos.shape==(50,3)
    def test_bonded_pairs_count(self):
        bp=bonded_pairs(5,10)
        assert len(bp)==5*9
    def test_bonded_pairs_no_inter_chain(self):
        bp=bonded_pairs(3,10)
        for i,j in bp:
            assert chain_index(i,10)==chain_index(j,10)
    def test_bonded_pairs_adjacent(self):
        bp=bonded_pairs(2,5)
        for i,j in bp:
            assert abs(j-i)==1
    def test_pbc_wraps(self):
        L=10.0
        pos=np.array([[11.5,-1.0,5.0],[20.0,9.5,3.0]])
        wp=apply_pbc(pos,L)
        assert 0<=wp[0,0]<L and 0<=wp[0,1]<L
        assert abs(wp[0,0]-1.5)<1e-9
        assert abs(wp[0,1]-9.0)<1e-9

class TestEnergyHelper:
    def test_lj_min_at_2to1over6_sigma(self):
        sig=1.0;eps=1.0
        rmin=2**(1/6)*sig
        e=lj_pair_energy(rmin,eps,sig)
        assert abs(e-(-eps))<1e-9
    def test_lj_zero_at_sigma(self):
        e=lj_pair_energy(1.0,1.0,1.0)
        assert abs(e)<1e-9

@needs_openmm
class TestIntegrator:
    def test_build_system(self):
        from melt.integrator import build_openmm_system
        mp=MeltParams(n_chains=2,chain_length=5,box_size=10.0)
        types=np.array([1,0,1,0,1,1,0,1,0,1],dtype=np.int8)
        s=build_openmm_system(types,mp)
        assert s.getNumParticles()==10
        assert s.getNumForces()==2
    def test_equilibrium_temperature(self):
        from melt.integrator import build_openmm_system,make_langevin_integrator,make_context,minimize_energy,kinetic_temperature,get_state_arrays
        mp=MeltParams(n_chains=5,chain_length=10,box_size=20.0,temperature=1.0,friction=1.0,dt=0.005)
        rng=np.random.default_rng(0)
        types=generate_random(50,0.5,rng)
        pos=init_chains_in_box(5,10,20.0,1.0,rng)
        sys=build_openmm_system(types,mp)
        integ=make_langevin_integrator(mp.temperature,mp.friction,mp.dt,seed=0)
        ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
        minimize_energy(ctx)
        integ.step(3000)
        Ts=[]
        for _ in range(40):
            integ.step(200)
            _,_,_,ke=get_state_arrays(ctx)
            Ts.append(kinetic_temperature(ke,50))
        Tm=float(np.mean(Ts))
        assert abs(Tm-1.0)<0.25

@needs_openmm
class TestEnergyConsistency:
    def test_lj_AA_value_matches_helper(self):
        from melt.integrator import build_openmm_system,make_langevin_integrator,make_context,get_state_arrays
        mp=MeltParams(n_chains=2,chain_length=1,box_size=20.0,bond_k=0.0,bond_r0=1.0,
                      lj_eps_AA=1.0,lj_eps_BB=1.0,lj_eps_AB=0.5,lj_sigma=1.0,lj_cutoff=4.0)
        types=np.array([1,1],dtype=np.int8)
        r=1.5
        pos=np.array([[5.0,5.0,5.0],[5.0+r,5.0,5.0]])
        sys=build_openmm_system(types,mp)
        integ=make_langevin_integrator(1.0,1.0,0.001,seed=0)
        ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
        _,_,pe,_=get_state_arrays(ctx)
        expected=lj_pair_energy(r,1.0,1.0)
        assert abs(pe-expected)<1e-4
    def test_lj_AB_value(self):
        from melt.integrator import build_openmm_system,make_langevin_integrator,make_context,get_state_arrays
        mp=MeltParams(n_chains=2,chain_length=1,box_size=20.0,
                      lj_eps_AA=1.0,lj_eps_BB=1.0,lj_eps_AB=0.5,lj_sigma=1.0,lj_cutoff=4.0)
        types=np.array([1,0],dtype=np.int8)
        r=1.5
        pos=np.array([[5.0,5.0,5.0],[5.0+r,5.0,5.0]])
        sys=build_openmm_system(types,mp)
        integ=make_langevin_integrator(1.0,1.0,0.001,seed=0)
        ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
        _,_,pe,_=get_state_arrays(ctx)
        expected=lj_pair_energy(r,0.5,1.0)
        assert abs(pe-expected)<1e-4

class TestDensity:
    def test_cic_mass_conservation(self):
        rng=np.random.default_rng(0)
        N=500;L=10.0
        pos=rng.uniform(0,L,size=(N,3))
        w=np.ones(N)
        g=grid_density_cic(pos,w,L,16)
        assert abs(g.sum()-N)<1e-9
    def test_uniform_density_constant(self):
        L=8.0;G=8
        pos=np.array([[(i+0.5)*L/G,(j+0.5)*L/G,(k+0.5)*L/G]
                      for i in range(G) for j in range(G) for k in range(G)])
        w=np.ones(len(pos))
        g=grid_density_cic(pos,w,L,G)
        assert np.allclose(g,1.0,atol=1e-9)
    def test_AB_partition(self):
        rng=np.random.default_rng(0)
        N=200;L=10.0
        pos=rng.uniform(0,L,size=(N,3))
        types=(rng.random(N)<0.5).astype(np.int8)
        gA=density_field_A(pos,types,L,16)
        gB=density_field_B(pos,types,L,16)
        assert abs(gA.sum()-types.sum())<1e-9
        assert abs(gB.sum()-(N-types.sum()))<1e-9
        assert abs((gA+gB).sum()-N)<1e-9

class TestStructureFactor:
    def test_uniform_density_is_flat(self):
        rng=np.random.default_rng(0)
        L=10.0;G=16
        phi=rng.normal(loc=1.0,scale=0.01,size=(G,G,G))
        S3=structure_factor_3d(phi,L)
        k,Sk=spherical_average(S3,L)
        finite=np.isfinite(Sk)&(k>0)
        assert finite.sum()>3
        Sm=Sk[finite].mean()
        assert Sm<1e-2
    def test_sinusoidal_pattern_peaks_at_known_k(self):
        L=10.0;G=32
        n=4
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
        xi=domain_length(kstar)
        assert abs(xi-3.5)<1e-9
    def test_oz_fit_recovers_xi(self):
        xi_true=2.0;S0=10.0
        k=np.linspace(0.05,1.5,30)
        S=S0/(1.0+(k*xi_true)**2)
        xi_fit,S0_fit=fit_ornstein_zernike(k,S,k_max=1.5)
        assert abs(xi_fit-xi_true)<0.05
        assert abs(S0_fit-S0)<0.5
    def test_compute_all_observables_shape(self):
        rng=np.random.default_rng(1)
        N=300;L=10.0;G=16
        pos=rng.uniform(0,L,size=(N,3))
        types=(rng.random(N)<0.5).astype(np.int8)
        phi_A=density_field_A(pos,types,L,G)
        phi_B=density_field_B(pos,types,L,G)
        out=compute_all_observables(phi_A,phi_B,L)
        assert "k" in out and "S_AA" in out and "xi_AA" in out
        assert len(out["k"])==len(out["S_AA"])

@needs_openmm
class TestReproducibility:
    def test_same_seed_same_trajectory(self):
        from melt.integrator import build_openmm_system,make_langevin_integrator,make_context,get_state_arrays
        def run_one(seed):
            mp=MeltParams(n_chains=3,chain_length=8,box_size=10.0,dt=0.005)
            rng=np.random.default_rng(seed)
            types=generate_random(24,0.5,rng)
            pos=init_chains_in_box(3,8,10.0,1.0,rng)
            sys=build_openmm_system(types,mp)
            integ=make_langevin_integrator(1.0,1.0,0.005,seed=seed)
            ctx=make_context(sys,integ,pos,mp.box_size,platform_name="Reference")
            integ.step(500)
            p,_,_,_=get_state_arrays(ctx)
            return p
        a=run_one(7)
        b=run_one(7)
        np.testing.assert_allclose(a,b,atol=1e-5)
