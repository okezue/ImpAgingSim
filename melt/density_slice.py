"""High-resolution real-space density slices from saved trajectories.

The corrected aging campaign saved trajectory.npz files at every fig5_aging condition.
Density slices written during the run used grid_size=56 (matching the structure-factor
grid), which is fine for FFT-based S(k) extraction but undersampled for publication
triptychs.

This module re-grids a saved trajectory at an arbitrary grid_size (default 128) and
returns the symmetric A-excess field (phi_A - phi_B) / (phi_A + phi_B) through a z-slab,
lightly Gaussian-smoothed. The output is a 2D float array suitable for imshow with
cmap='RdBu_r' and symmetric vmin/vmax.

Usage:
    from melt.density_slice import slab_from_trajectory
    slab, L = slab_from_trajectory("output/.../trajectory.npz",
                                   n_grid=128, slab_width_sigma=3.0,
                                   smooth_sigma_voxels=2.0, last_n_snapshots=8)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from .density import density_field_A,density_field_B

def slab_from_grids(phi_A,phi_B,L,slab_width_sigma=3.0,smooth_sigma_voxels=2.0):
    """Project density grids onto a z-slab and return (phi_A - phi_B) / (phi_A + phi_B).

    phi_A, phi_B : (G, G, G) density grids (after CIC deposition).
    L            : box edge length in sigma.
    """
    from scipy.ndimage import gaussian_filter
    G=phi_A.shape[0]
    h=L/G
    width_vox=max(2,int(round(slab_width_sigma/h)))
    z0=G//2-width_vox//2
    A=phi_A[:,:,z0:z0+width_vox].mean(2)
    B=phi_B[:,:,z0:z0+width_vox].mean(2)
    A=gaussian_filter(A,smooth_sigma_voxels)
    B=gaussian_filter(B,smooth_sigma_voxels)
    tot=A+B
    return np.where(tot>1e-6,(A-B)/np.maximum(tot,1e-9),0.0)

def slab_from_trajectory(traj_npz,n_grid=128,slab_width_sigma=3.0,smooth_sigma_voxels=2.0,
                         last_n_snapshots=8):
    """Re-grid positions from a saved trajectory.npz at chosen n_grid, average the last
    last_n_snapshots frames, project on a z-slab, return symmetric A-excess slice and L.
    """
    d=np.load(str(traj_npz),allow_pickle=False)
    pos=d["positions"];types=d["types"];L=float(d["box_size"])
    n=pos.shape[0]
    k=min(int(last_n_snapshots),n)
    A_sum=np.zeros((n_grid,n_grid,n_grid),dtype=np.float64)
    B_sum=np.zeros_like(A_sum)
    for t in range(n-k,n):
        A_sum+=density_field_A(pos[t],types,L,n_grid)
        B_sum+=density_field_B(pos[t],types,L,n_grid)
    A_sum/=k;B_sum/=k
    slab=slab_from_grids(A_sum,B_sum,L,slab_width_sigma,smooth_sigma_voxels)
    return slab,L

def triptych_from_runs(run_paths,n_grid=128,slab_width_sigma=3.0,smooth_sigma_voxels=2.0,
                       last_n_snapshots=8):
    """Compute slabs for a list of run directories (typically kappa=0, 0.5, 1.0).
    Each entry in run_paths is a Path or string pointing to a directory containing
    trajectory.npz. Returns (slabs, L, vmax_99) where slabs is a list and vmax_99
    is the 99th percentile of |slab| across all entries for symmetric color scaling.
    """
    slabs=[];L=None
    for rp in run_paths:
        rp=Path(rp)
        tjp=rp/"trajectory.npz"
        if not tjp.exists():raise FileNotFoundError(tjp)
        s,L=slab_from_trajectory(tjp,n_grid=n_grid,slab_width_sigma=slab_width_sigma,
                                 smooth_sigma_voxels=smooth_sigma_voxels,
                                 last_n_snapshots=last_n_snapshots)
        slabs.append(s)
    vmax_99=float(max(np.percentile(np.abs(s),99) for s in slabs)) if slabs else 1.0
    return slabs,L,vmax_99

def _cli():
    import argparse
    p=argparse.ArgumentParser(description="Render high-resolution density slice from a saved trajectory.npz")
    p.add_argument("trajectory",help="path to trajectory.npz")
    p.add_argument("--n_grid",type=int,default=128)
    p.add_argument("--slab_width_sigma",type=float,default=3.0)
    p.add_argument("--smooth_sigma_voxels",type=float,default=2.0)
    p.add_argument("--last_n_snapshots",type=int,default=8)
    p.add_argument("--out",type=str,default=None,help="optional .npz path to save the slab")
    a=p.parse_args()
    slab,L=slab_from_trajectory(a.trajectory,n_grid=a.n_grid,slab_width_sigma=a.slab_width_sigma,
                                smooth_sigma_voxels=a.smooth_sigma_voxels,
                                last_n_snapshots=a.last_n_snapshots)
    print(f"slab shape={slab.shape}, L={L:.3f} sigma, range=[{slab.min():.3f}, {slab.max():.3f}]")
    if a.out:
        np.savez_compressed(a.out,slab=slab,L=L,n_grid=a.n_grid,
                            slab_width_sigma=a.slab_width_sigma,
                            smooth_sigma_voxels=a.smooth_sigma_voxels)
        print(f"saved -> {a.out}")

if __name__=="__main__":
    _cli()
