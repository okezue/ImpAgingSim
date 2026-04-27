from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class MeltParams:
    n_chains:int=10
    chain_length:int=20
    box_size:float=10.0
    bond_k:float=100.0
    bond_r0:float=1.0
    lj_eps_AA:float=1.0
    lj_eps_BB:float=1.0
    lj_eps_AB:float=0.5
    lj_sigma:float=1.0
    lj_cutoff:float=2.5
    temperature:float=1.0
    friction:float=1.0
    dt:float=0.005

@dataclass(frozen=True)
class RunParams:
    n_steps:int=10000
    equilibration_steps:int=1000
    snapshot_interval:int=100
    seed:int=12345

def lj_pair_energy(r,eps,sigma):
    sr=sigma/r
    sr6=sr**6
    return 4.0*eps*(sr6*sr6-sr6)
