from __future__ import annotations
import argparse,os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter,FFMpegWriter
from .density import density_field_A,density_field_B

def density_slice_movie(traj_npz,outpath,grid_size=48,axis=2,fps=10,dpi=100):
    d=np.load(traj_npz)
    pos=d["positions"];types=d["types"];L=float(d["box_size"]);steps=d["steps"]
    n_frames=pos.shape[0]
    fig,ax=plt.subplots(1,2,figsize=(11,5))
    phi_A0=density_field_A(pos[0],types,L,grid_size)
    phi_B0=density_field_B(pos[0],types,L,grid_size)
    sl_A=phi_A0.mean(axis=axis);sl_B=phi_B0.mean(axis=axis)
    vmax=max(sl_A.max(),sl_B.max())*1.5
    im0=ax[0].imshow(sl_A.T,origin="lower",cmap="Reds",vmin=0,vmax=vmax,
                     extent=[0,L,0,L])
    im1=ax[1].imshow(sl_B.T,origin="lower",cmap="Blues",vmin=0,vmax=vmax,
                     extent=[0,L,0,L])
    ax[0].set_title("φ_A density (z-projection)");ax[1].set_title("φ_B density (z-projection)")
    for a in ax:a.set_xlabel("x");a.set_ylabel("y")
    suptitle=fig.suptitle(f"step {steps[0]}")
    plt.tight_layout()
    def update(i):
        phi_A=density_field_A(pos[i],types,L,grid_size).mean(axis=axis)
        phi_B=density_field_B(pos[i],types,L,grid_size).mean(axis=axis)
        im0.set_data(phi_A.T);im1.set_data(phi_B.T)
        suptitle.set_text(f"step {steps[i]}")
        return im0,im1,suptitle
    anim=FuncAnimation(fig,update,frames=n_frames,interval=1000//fps,blit=False)
    if outpath.endswith(".mp4"):
        try:
            anim.save(outpath,writer=FFMpegWriter(fps=fps),dpi=dpi)
        except Exception:
            outpath=outpath.replace(".mp4",".gif")
            anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    else:
        anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    plt.close(fig)
    return outpath

def Sk_evolution_movie(sf_npz,outpath,fps=10,dpi=100):
    d=np.load(sf_npz)
    steps=d["steps"];k=d["k"];SAA=d["S_AA"];SBB=d["S_BB"];SAB=d["S_AB"]
    fig,ax=plt.subplots(figsize=(7,5))
    lA,=ax.plot(k,SAA[0],"r-",label="S_AA(k)",lw=2)
    lB,=ax.plot(k,SBB[0],"b-",label="S_BB(k)",lw=2)
    lAB,=ax.plot(k,SAB[0],"g-",label="S_AB(k)",lw=2)
    ax.set_xscale("log");ax.set_yscale("symlog",linthresh=1e-4)
    ymax=max(SAA.max(),SBB.max())*1.5
    ax.set_ylim(-ymax*0.1,ymax)
    ax.set_xlabel("k");ax.set_ylabel("S(k)");ax.legend(frameon=False)
    title=ax.set_title(f"step {steps[0]}")
    plt.tight_layout()
    def update(i):
        lA.set_ydata(SAA[i]);lB.set_ydata(SBB[i]);lAB.set_ydata(SAB[i])
        title.set_text(f"step {steps[i]}")
        return lA,lB,lAB,title
    anim=FuncAnimation(fig,update,frames=len(steps),interval=1000//fps,blit=False)
    if outpath.endswith(".mp4"):
        try:
            anim.save(outpath,writer=FFMpegWriter(fps=fps),dpi=dpi)
        except Exception:
            outpath=outpath.replace(".mp4",".gif")
            anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    else:
        anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    plt.close(fig)
    return outpath

def polymer_3d_snapshot(traj_npz,outpath,frame=-1,dpi=140):
    d=np.load(traj_npz)
    pos=d["positions"][frame];types=d["types"];L=float(d["box_size"]);step=int(d["steps"][frame])
    fig=plt.figure(figsize=(7,7))
    ax=fig.add_subplot(111,projection="3d")
    A=types==1
    ax.scatter(pos[A,0],pos[A,1],pos[A,2],c="red",s=8,label="A",alpha=0.7)
    ax.scatter(pos[~A,0],pos[~A,1],pos[~A,2],c="blue",s=8,label="B",alpha=0.7)
    ax.set_xlim(0,L);ax.set_ylim(0,L);ax.set_zlim(0,L)
    ax.set_xlabel("x");ax.set_ylabel("y");ax.set_zlabel("z")
    ax.set_title(f"step {step}");ax.legend(loc="upper right")
    plt.tight_layout();plt.savefig(outpath,dpi=dpi);plt.close()
    return outpath

def polymer_3d_movie(traj_npz,outpath,fps=10,dpi=80,stride=1):
    d=np.load(traj_npz)
    pos=d["positions"];types=d["types"];L=float(d["box_size"]);steps=d["steps"]
    A=types==1
    fig=plt.figure(figsize=(7,7))
    ax=fig.add_subplot(111,projection="3d")
    pA0=pos[0][A];pB0=pos[0][~A]
    sA=ax.scatter(pA0[:,0],pA0[:,1],pA0[:,2],c="red",s=10,alpha=0.7,label="A")
    sB=ax.scatter(pB0[:,0],pB0[:,1],pB0[:,2],c="blue",s=10,alpha=0.7,label="B")
    ax.set_xlim(0,L);ax.set_ylim(0,L);ax.set_zlim(0,L)
    ax.legend(loc="upper right")
    title=ax.set_title(f"step {steps[0]}")
    plt.tight_layout()
    frames=list(range(0,len(steps),stride))
    def update(i):
        pi=pos[i]
        sA._offsets3d=(pi[A,0],pi[A,1],pi[A,2])
        sB._offsets3d=(pi[~A,0],pi[~A,1],pi[~A,2])
        title.set_text(f"step {steps[i]}")
        return sA,sB,title
    anim=FuncAnimation(fig,update,frames=frames,interval=1000//fps,blit=False)
    if outpath.endswith(".mp4"):
        try:
            anim.save(outpath,writer=FFMpegWriter(fps=fps),dpi=dpi)
        except Exception:
            outpath=outpath.replace(".mp4",".gif")
            anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    else:
        anim.save(outpath,writer=PillowWriter(fps=fps),dpi=dpi)
    plt.close(fig)
    return outpath

def main():
    p=argparse.ArgumentParser()
    p.add_argument("run_dir",type=str)
    p.add_argument("--outdir",type=str,default=None)
    p.add_argument("--fps",type=int,default=10)
    p.add_argument("--grid_size",type=int,default=48)
    args=p.parse_args()
    od=args.outdir or args.run_dir
    traj_p=os.path.join(args.run_dir,"trajectory.npz")
    sf_p=os.path.join(args.run_dir,"structure_factor.npz")
    written=[]
    if os.path.exists(traj_p):
        out=density_slice_movie(traj_p,os.path.join(od,"density_slice.gif"),
                                grid_size=args.grid_size,fps=args.fps)
        written.append(out)
        out=polymer_3d_snapshot(traj_p,os.path.join(od,"final_3d.png"))
        written.append(out)
        out=polymer_3d_movie(traj_p,os.path.join(od,"polymer_3d.gif"),fps=args.fps)
        written.append(out)
    if os.path.exists(sf_p):
        out=Sk_evolution_movie(sf_p,os.path.join(od,"Sk_evolution.gif"),fps=args.fps)
        written.append(out)
    for w in written:
        print(f"  wrote {w}")

if __name__=="__main__":
    main()
