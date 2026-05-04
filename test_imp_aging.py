import math
import numpy as np
import pytest
from imp_aging import (
    ModelParams,RunParams,
    generate_eta_iid,generate_eta_correlated,generate_sigma_markov,
    make_contact_pair_indices,dist_sq_matrix,
    delta_energy_move,metropolis_sweep,init_random_walk,
    run_aging_trajectory,compute_energy,compute_Rg,compute_ncontacts,
    contact_vector_from_dist_sq,compute_observables_for_trajectory,
    log_lags,build_snapshot_times,lj_rmin,contact_cutoff,
)

def full_energy(pos,params,eta,sqrt_eps):
    N=pos.shape[0]
    E=0.0
    for i in range(N):
        for j in range(i+1,N):
            d2=float(np.sum((pos[i]-pos[j])**2))
            d2=max(d2,1e-12)
            inv2=1.0/d2;inv6=inv2**3;inv12=inv6**2
            E+=params.R*inv12-params.A*inv6+sqrt_eps*eta[i,j]*inv6
            if j==i+1:
                E+=params.h*d2
    return E

class TestEnergyFunction:
    def test_delta_energy_consistent_with_full(self):
        N=5
        p=ModelParams(N=N,A=3.8,R=2.0,h=1.0)
        rng=np.random.default_rng(42)
        eta=generate_eta_iid(N,rng)
        se=math.sqrt(6.0)
        pos=init_random_walk(N,1.0,rng)
        E0=full_energy(pos,p,eta,se)
        i=2
        new_p=pos[i]+rng.uniform(-0.25,0.25,size=3)
        dE=delta_energy_move(pos,i,new_p,p,eta,se)
        pos2=pos.copy()
        pos2[i]=new_p
        E1=full_energy(pos2,p,eta,se)
        assert abs(dE-(E1-E0))<1e-8

    def test_delta_energy_zero_for_no_move(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(7)
        eta=generate_eta_iid(N,rng)
        pos=init_random_walk(N,1.0,rng)
        dE=delta_energy_move(pos,1,pos[1].copy(),p,eta,math.sqrt(6.0))
        assert abs(dE)<1e-12

    def test_two_beads_known_distance(self):
        N=2
        p=ModelParams(N=N,A=3.8,R=2.0,h=1.0)
        eta=np.zeros((2,2))
        pos=np.array([[0.0,0.0,0.0],[1.0,0.0,0.0]])
        pos-=pos.mean(axis=0)
        E=full_energy(pos,p,eta,0.0)
        d2=1.0
        exp=p.R*1.0-p.A*1.0+p.h*d2
        assert abs(E-exp)<1e-12

    def test_delta_energy_boundary_beads(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(99)
        eta=generate_eta_iid(N,rng)
        se=math.sqrt(6.0)
        pos=init_random_walk(N,1.0,rng)
        for i in [0,N-1]:
            E0=full_energy(pos,p,eta,se)
            np_i=pos[i]+rng.uniform(-0.1,0.1,size=3)
            dE=delta_energy_move(pos,i,np_i,p,eta,se)
            pos2=pos.copy();pos2[i]=np_i
            E1=full_energy(pos2,p,eta,se)
            assert abs(dE-(E1-E0))<1e-8

    def test_compute_energy_matches_full(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(55)
        eta=generate_eta_iid(N,rng)
        se=math.sqrt(6.0)
        pos=init_random_walk(N,1.0,rng)
        E1=full_energy(pos,p,eta,se)
        E2=compute_energy(pos,p,eta,se)
        assert abs(E1-E2)<1e-6

class TestDisorderIID:
    def test_mean_var(self):
        N=200
        rng=np.random.default_rng(1)
        eta=generate_eta_iid(N,rng)
        iu=np.triu_indices(N,1)
        v=eta[iu]
        assert abs(np.mean(v))<0.1
        assert abs(np.var(v)-1.0)<0.15

    def test_symmetry(self):
        N=20
        rng=np.random.default_rng(2)
        eta=generate_eta_iid(N,rng)
        assert np.allclose(eta,eta.T)

    def test_diagonal_zero(self):
        N=20
        rng=np.random.default_rng(3)
        eta=generate_eta_iid(N,rng)
        assert np.all(np.diag(eta)==0.0)

class TestDisorderCorrelated:
    def test_mean_var(self):
        N=200
        rng=np.random.default_rng(10)
        eta,sig=generate_eta_correlated(N,pi=0.9,kappa=0.7,rng=rng)
        iu=np.triu_indices(N,1)
        v=eta[iu]
        assert abs(np.mean(v))<0.15
        assert abs(np.var(v)-1.0)<0.2

    def test_symmetry(self):
        N=20
        rng=np.random.default_rng(11)
        eta,_=generate_eta_correlated(N,pi=0.9,kappa=0.7,rng=rng)
        assert np.allclose(eta,eta.T)

    def test_diagonal_zero(self):
        N=20
        rng=np.random.default_rng(12)
        eta,_=generate_eta_correlated(N,pi=0.9,kappa=0.7,rng=rng)
        assert np.all(np.diag(eta)==0.0)

    def test_sigma_persistence(self):
        N=10000
        rng=np.random.default_rng(13)
        pi=0.9
        sig=generate_sigma_markov(N,pi,rng)
        same=(sig[:-1]==sig[1:]).mean()
        assert abs(same-pi)<0.02

    def test_sigma_values(self):
        N=100
        rng=np.random.default_rng(14)
        sig=generate_sigma_markov(N,0.5,rng)
        assert set(np.unique(sig)).issubset({-1,1})

class TestObservables:
    def test_contact_vector_identifies_contacts(self):
        pos=np.array([[0,0,0],[0.5,0,0],[10,0,0]],dtype=float)
        dsq=dist_sq_matrix(pos)
        pi_arr=np.array([0,0,1],dtype=np.int32)
        pj_arr=np.array([1,2,2],dtype=np.int32)
        cv=contact_vector_from_dist_sq(dsq,pi_arr,pj_arr,rc2=1.0)
        assert cv[0]==True
        assert cv[1]==False
        assert cv[2]==False

    def test_Q_identical_configs(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(20)
        pos=init_random_walk(N,lj_rmin(p),rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos.copy(),1:pos.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        Q=res["Q"]
        dsq=dist_sq_matrix(pos)
        cv=contact_vector_from_dist_sq(dsq,pi,pj,rc2)
        exp_q=float(np.mean(cv&cv))
        assert abs(Q[0,0]-exp_q)<1e-12

    def test_Q_different_configs(self):
        N=8
        p=ModelParams(N=N)
        rng=np.random.default_rng(21)
        pos1=init_random_walk(N,lj_rmin(p),rng)
        pos2=init_random_walk(N,lj_rmin(p),rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos1.copy(),1:pos2.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        Q=res["Q"]
        dsq1=dist_sq_matrix(pos1);dsq2=dist_sq_matrix(pos2)
        c1=contact_vector_from_dist_sq(dsq1,pi,pj,rc2)
        c2=contact_vector_from_dist_sq(dsq2,pi,pj,rc2)
        exp_q=float(np.mean(c1&c2))
        assert abs(Q[0,0]-exp_q)<1e-12

    def test_chi4_nonneg(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(22)
        eta=generate_eta_iid(N,rng)
        pi_a,pj_a=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        n_pairs=len(pi_a)
        Q_list=[]
        for t in range(4):
            rng2=np.random.default_rng(100+t)
            pos=init_random_walk(N,1.0,rng2)
            snaps={0:pos.copy(),1:pos+rng2.normal(0,0.1,pos.shape)}
            res=compute_observables_for_trajectory(snaps,[0],[1],pi_a,pj_a,rc2)
            Q_list.append(res["Q"])
        Qa=np.stack(Q_list)
        Qm=Qa.mean(axis=0)
        Q2m=(Qa**2).mean(axis=0)
        chi4=n_pairs*(Q2m-Qm**2)
        assert np.all(chi4>=-1e-12)

    def test_D4_zero_identical(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(23)
        pos=init_random_walk(N,1.0,rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos.copy(),1:pos.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        assert abs(res["D4"][0,0])<1e-12

class TestContactPairIndices:
    def test_excludes_nearest_neighbors(self):
        N=8
        pi,pj=make_contact_pair_indices(N)
        diffs=pj-pi
        assert np.all(diffs>1)

    def test_correct_count(self):
        N=8
        pi,pj=make_contact_pair_indices(N)
        exp=N*(N-1)//2-(N-1)
        assert len(pi)==exp

    def test_all_i_less_j(self):
        N=6
        pi,pj=make_contact_pair_indices(N)
        assert np.all(pi<pj)

class TestReproducibility:
    def test_same_seed_same_result(self):
        N=5
        p=ModelParams(N=N)
        eta=generate_eta_iid(N,np.random.default_rng(50))
        se=math.sqrt(6.0)
        rng1=np.random.default_rng(77)
        pos1=init_random_walk(N,1.0,rng1)
        metropolis_sweep(pos1,1.0,0.25,p,eta,se,rng1)
        rng2=np.random.default_rng(77)
        pos2=init_random_walk(N,1.0,rng2)
        metropolis_sweep(pos2,1.0,0.25,p,eta,se,rng2)
        assert np.allclose(pos1,pos2)

    def test_different_seeds_differ(self):
        N=5
        p=ModelParams(N=N)
        eta=generate_eta_iid(N,np.random.default_rng(50))
        se=math.sqrt(6.0)
        rng1=np.random.default_rng(77)
        pos1=init_random_walk(N,1.0,rng1)
        metropolis_sweep(pos1,1.0,0.25,p,eta,se,rng1)
        rng2=np.random.default_rng(88)
        pos2=init_random_walk(N,1.0,rng2)
        metropolis_sweep(pos2,1.0,0.25,p,eta,se,rng2)
        assert not np.allclose(pos1,pos2)

class TestGeometryHelpers:
    def test_dist_sq_matrix(self):
        pos=np.array([[0,0,0],[3,4,0]],dtype=float)
        ds=dist_sq_matrix(pos)
        assert abs(ds[0,1]-25.0)<1e-12
        assert abs(ds[1,0]-25.0)<1e-12
        assert abs(ds[0,0])<1e-12

    def test_dist_sq_matrix_symmetric(self):
        rng=np.random.default_rng(30)
        pos=rng.normal(size=(8,3))
        ds=dist_sq_matrix(pos)
        assert np.allclose(ds,ds.T)

    def test_log_lags_unique_sorted(self):
        lags=log_lags(1000,20)
        assert len(lags)==len(set(lags))
        assert np.all(np.diff(lags)>0)

    def test_log_lags_includes_zero(self):
        lags=log_lags(1000,20)
        assert lags[0]==0
    def test_log_lags_bounds(self):
        lags=log_lags(500,15)
        assert lags[0]==0
        assert lags[-1]<=500

    def test_log_lags_integer(self):
        lags=log_lags(100,10)
        assert lags.dtype in [np.int32,np.int64,int]

    def test_build_snapshot_times_includes_all(self):
        tw=[0,10]
        lags=[1,5]
        st=build_snapshot_times(tw,lags)
        for tw_v in tw:
            assert tw_v in st
            for l in lags:
                assert tw_v+l in st

class TestMetropolisSweep:
    def test_acceptance_rate_range(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(40)
        eta=generate_eta_iid(N,rng)
        pos=init_random_walk(N,1.0,rng)
        ar=metropolis_sweep(pos,1.0,0.25,p,eta,math.sqrt(6.0),rng)
        assert 0.0<=ar<=1.0

    def test_high_temp_high_acceptance(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(41)
        eta=generate_eta_iid(N,rng)
        pos=init_random_walk(N,lj_rmin(p),rng)
        for _ in range(20):
            metropolis_sweep(pos,0.01,0.01,p,eta,0.0,rng)
        rates=[]
        for _ in range(10):
            ar=metropolis_sweep(pos,0.01,0.01,p,eta,0.0,rng)
            rates.append(ar)
        assert np.mean(rates)>0.5

    def test_com_near_zero(self):
        N=8
        p=ModelParams(N=N)
        rng=np.random.default_rng(42)
        eta=generate_eta_iid(N,rng)
        pos=init_random_walk(N,1.0,rng)
        for _ in range(50):
            metropolis_sweep(pos,5.0,0.25,p,eta,math.sqrt(6.0),rng)
        com=pos.mean(axis=0)
        assert np.linalg.norm(com)<1e-10

class TestFullTrajectory:
    def test_returns_all_snapshots(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(60)
        eta=generate_eta_iid(N,rng)
        st=[0,1,3,5]
        result=run_aging_trajectory(
            params=p,eta=eta,epsilon=6.0,
            beta0=0.05,beta=1.0,pre_sweeps=5,
            meas_sweeps=5,step_size=0.25,
            snapshot_times=st,rng=rng,bond_length=1.0,
        )
        snaps=result["snapshots"]
        for t in st:
            assert t in snaps

    def test_snapshot_shape(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(61)
        eta=generate_eta_iid(N,rng)
        st=[0,2]
        result=run_aging_trajectory(
            params=p,eta=eta,epsilon=6.0,
            beta0=0.05,beta=1.0,pre_sweeps=3,
            meas_sweeps=2,step_size=0.25,
            snapshot_times=st,rng=rng,bond_length=1.0,
        )
        snaps=result["snapshots"]
        for t in st:
            assert snaps[t].shape==(N,3)

    def test_snapshots_are_copies(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(62)
        eta=generate_eta_iid(N,rng)
        st=[0,1,2]
        result=run_aging_trajectory(
            params=p,eta=eta,epsilon=6.0,
            beta0=0.05,beta=1.0,pre_sweeps=3,
            meas_sweeps=2,step_size=0.25,
            snapshot_times=st,rng=rng,bond_length=1.0,
        )
        snaps=result["snapshots"]
        assert not np.shares_memory(snaps[0],snaps[1])

    def test_returns_energy_and_rg(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(63)
        eta=generate_eta_iid(N,rng)
        st=[0,1]
        result=run_aging_trajectory(
            params=p,eta=eta,epsilon=6.0,
            beta0=0.05,beta=1.0,pre_sweeps=3,
            meas_sweeps=1,step_size=0.25,
            snapshot_times=st,rng=rng,bond_length=1.0,
        )
        assert "energy" in result
        assert "Rg" in result
        for t in st:
            assert t in result["energy"]
            assert t in result["Rg"]
            assert isinstance(result["energy"][t],float)
            assert result["Rg"][t]>=0
class TestNewObservables:
    def test_msd_zero_for_identical(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(70)
        pos=init_random_walk(N,1.0,rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos.copy(),1:pos.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        assert abs(res["MSD"][0,0])<1e-12
    def test_msd_positive_for_different(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(71)
        pos1=init_random_walk(N,1.0,rng)
        pos2=pos1+rng.normal(0,0.5,pos1.shape)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos1.copy(),1:pos2.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        assert res["MSD"][0,0]>0
    def test_alpha2_zero_for_identical(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(72)
        pos=init_random_walk(N,1.0,rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        snaps={0:pos.copy(),1:pos.copy()}
        res=compute_observables_for_trajectory(snaps,[0],[1],pi,pj,rc2)
        assert abs(res["alpha2"][0,0])<1e-12
    def test_Rg_positive(self):
        N=8
        rng=np.random.default_rng(73)
        pos=init_random_walk(N,1.0,rng)
        rg=compute_Rg(pos)
        assert rg>0
    def test_Rg_collinear(self):
        pos=np.array([[float(i),0,0] for i in range(5)])
        pos-=pos.mean(axis=0)
        rg=compute_Rg(pos)
        exp=math.sqrt(np.mean(np.sum(pos**2,axis=1)))
        assert abs(rg-exp)<1e-10
    def test_ncontacts_nonneg(self):
        N=8
        p=ModelParams(N=N)
        rng=np.random.default_rng(74)
        pos=init_random_walk(N,lj_rmin(p),rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        nc=compute_ncontacts(pos,pi,pj,rc2)
        assert nc>=0
        assert isinstance(nc,int)
    def test_trajectory_returns_ncontacts(self):
        N=5
        p=ModelParams(N=N)
        rng=np.random.default_rng(75)
        eta=generate_eta_iid(N,rng)
        pi,pj=make_contact_pair_indices(N)
        rc2=contact_cutoff(p)**2
        st=[0,1]
        result=run_aging_trajectory(
            params=p,eta=eta,epsilon=6.0,
            beta0=0.05,beta=1.0,pre_sweeps=3,
            meas_sweeps=1,step_size=0.25,
            snapshot_times=st,rng=rng,bond_length=1.0,
            pair_i=pi,pair_j=pj,rc2=rc2,
        )
        assert "n_contacts" in result
        for t in st:
            assert t in result["n_contacts"]
class TestKWWFitting:
    def test_kww_pure_exponential(self):
        from extended_analysis import fit_kww
        lag=np.logspace(0,4,50)
        qn=np.exp(-lag/500.0)
        tau,beta,r2=fit_kww(lag,qn)
        assert abs(beta-1.0)<0.1
        assert abs(tau-500.0)<50
        assert r2>0.99
    def test_kww_stretched(self):
        from extended_analysis import fit_kww
        lag=np.logspace(0,4,50)
        qn=np.exp(-(lag/300.0)**0.5)
        tau,beta,r2=fit_kww(lag,qn)
        assert abs(beta-0.5)<0.1
        assert r2>0.95
    def test_kww_nan_for_insufficient(self):
        from extended_analysis import fit_kww
        lag=np.array([1.0,2.0])
        qn=np.array([0.9,0.8])
        tau,beta,r2=fit_kww(lag,qn)
        assert math.isnan(tau)
