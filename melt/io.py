from __future__ import annotations
import csv,json,os
from dataclasses import asdict
import numpy as np

CSV_FIELDS=["run_id","sequence_type","n_chains","chain_length","f_A","kappa","pi",
            "step","total_energy","kinetic_energy","mean_Rg","temperature_inst",
            "k_star","xi_AA","S_AA_peak"]

def make_run_dir(out_root,run_id):
    d=os.path.join(out_root,run_id)
    os.makedirs(d,exist_ok=True)
    return d

def write_meta(run_dir,mp,rp,sequence_type,sequence,extra=None):
    meta={"melt_params":asdict(mp),"run_params":asdict(rp),
          "sequence_type":sequence_type,"sequence":sequence.tolist(),
          "n_particles":int(mp.n_chains*mp.chain_length)}
    try:
        import openmm
        meta["openmm_version"]=openmm.version.short_version
    except Exception:
        meta["openmm_version"]=None
    if extra:
        meta.update(extra)
    with open(os.path.join(run_dir,"meta.json"),"w") as f:
        json.dump(meta,f,indent=2,sort_keys=True)

def open_csv(run_dir):
    p=os.path.join(run_dir,"snapshots.csv")
    new=not os.path.exists(p) or os.path.getsize(p)==0
    f=open(p,"a",newline="")
    w=csv.writer(f)
    if new:
        w.writerow(CSV_FIELDS)
    return f,w

def append_row(writer,run_id,sequence_type,mp,kappa,pi,f_A,step,pe,ke,rg,T_inst,
               k_star=float("nan"),xi_AA=float("nan"),S_AA_peak=float("nan")):
    writer.writerow([run_id,sequence_type,int(mp.n_chains),int(mp.chain_length),
                     float(f_A),float(kappa),float(pi),
                     int(step),float(pe),float(ke),float(rg),float(T_inst),
                     float(k_star),float(xi_AA),float(S_AA_peak)])

def save_structure_factors(run_dir,steps,ks,S_AA,S_BB,S_AB):
    p=os.path.join(run_dir,"structure_factor.npz")
    np.savez_compressed(p,steps=np.asarray(steps),k=np.asarray(ks),
                        S_AA=np.asarray(S_AA),S_BB=np.asarray(S_BB),S_AB=np.asarray(S_AB))
    return p

def save_trajectory(run_dir,steps,positions,types,box_size):
    p=os.path.join(run_dir,"trajectory.npz")
    np.savez_compressed(p,steps=np.asarray(steps),
                        positions=np.asarray(positions,dtype=np.float32),
                        types=np.asarray(types,dtype=np.int8),
                        box_size=float(box_size))
    return p

def save_density_grids(run_dir,steps,phi_A_grids,phi_B_grids):
    p=os.path.join(run_dir,"density_grids.npz")
    np.savez_compressed(p,steps=np.asarray(steps),
                        phi_A=np.asarray(phi_A_grids,dtype=np.float32),
                        phi_B=np.asarray(phi_B_grids,dtype=np.float32))
    return p

def compute_mean_Rg(pos,n_chains,chain_length):
    rgs=[]
    for c in range(n_chains):
        s=c*chain_length;e=s+chain_length
        cm=pos[s:e].mean(axis=0)
        rgs.append(float(np.sqrt(np.mean(np.sum((pos[s:e]-cm)**2,axis=1)))))
    return float(np.mean(rgs))
