import numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter
from matplotlib import animation

d=np.load("/tmp/prequench_grids_v2.npz",allow_pickle=True)
L=float(d["L"]); n_eq=int(d["n_eq"]); T_eq=float(d["T_eq"]); T_q=float(d["T_q"])

def grids(k):
    return d[f"{k}_phi_A"],d[f"{k}_phi_B"],d[f"{k}_steps"],d[f"{k}_phase"]

A0,B0,steps,phase=grids("k0"); A05,B05,_,_=grids("k05"); A1,B1,_,_=grids("k1")
T,G,_,_=A0.shape
w=max(2,int(G*0.18)); z0=G//2-w//2

def slab(phiA,phiB,t):
    A=phiA[t,:,:,z0:z0+w].mean(2); B=phiB[t,:,:,z0:z0+w].mean(2)
    A=gaussian_filter(A,1.6); B=gaussian_filter(B,1.6)
    tot=A+B
    return np.where(tot>1e-6,(A-B)/np.maximum(tot,1e-9),0.0)

vmax=max(np.percentile(np.abs(slab(A1,B1,T-1)),98),0.6)

plt.rcParams.update({"font.size":12,"savefig.facecolor":"white","figure.facecolor":"white"})
fig,axes=plt.subplots(1,3,figsize=(13,5.6))
ims=[]
for ax,(phiA,phiB),lab in zip(axes,[(A0,B0),(A05,B05),(A1,B1)],
                              [r"$\kappa=0$",r"$\kappa=0.5$",r"$\kappa=1$"]):
    im=ax.imshow(slab(phiA,phiB,0).T,origin="lower",extent=[0,L,0,L],
        cmap="RdBu_r",vmin=-vmax,vmax=vmax,aspect="equal",interpolation="bilinear")
    ax.set_title(lab,fontsize=14,pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    ims.append(im)
sup=fig.suptitle("",fontsize=12,y=0.985)
cax=fig.add_axes([0.22,0.13,0.56,0.030])
cb=fig.colorbar(ims[0],cax=cax,orientation="horizontal")
cb.set_label(r"local A excess $(\phi_A-\phi_B)/(\phi_A+\phi_B)$",fontsize=10)
cb.ax.tick_params(labelsize=9)
fig.subplots_adjust(left=0.03,right=0.97,top=0.84,bottom=0.26,wspace=0.06)

def upd(i):
    ims[0].set_data(slab(A0,B0,i).T)
    ims[1].set_data(slab(A05,B05,i).T)
    ims[2].set_data(slab(A1,B1,i).T)
    ph=phase[i]
    if ph=="min":
        tag=f"step 0, post-minimization (T=0)"
    elif ph=="eq":
        tag=f"step {steps[i]:,}, equilibration at T={T_eq:.1f}"
    else:
        tag=f"step {steps[i]:,}, quenched at T={T_q:.1f}"
    sup.set_text(tag)
    return ims+[sup]

ani=animation.FuncAnimation(fig,upd,frames=range(T),interval=80,blit=False)
out=Path.home()/"Downloads/melt_quench_evolution.gif"
ani.save(out,writer=animation.PillowWriter(fps=12),dpi=100)
print(f"saved {out}, frames={T}")
