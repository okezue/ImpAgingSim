#!/usr/bin/env python3
from __future__ import annotations
import os,math,argparse,warnings
from typing import Dict,List,Tuple,Optional
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
matplotlib.rcParams.update({
    'font.size':11,'axes.labelsize':13,'axes.titlesize':13,
    'legend.fontsize':9,'figure.dpi':200,'savefig.dpi':200,
    'figure.figsize':(7,5),'lines.linewidth':1.5,
    'lines.markersize':5,'axes.grid':True,'grid.alpha':0.3,
})
ENS_STYLE={"iid":("tab:blue","o","-"),"correlated":("tab:red","s","--")}
EPS_COLORS={0.0:"tab:green",3.0:"tab:orange",6.0:"tab:purple"}
def kww(t,tau,beta):
    return np.exp(-(t/tau)**beta)
def fit_kww(lag,qn,p0=None):
    m=(lag>0)&np.isfinite(qn)&(qn>0)&(qn<=1.0)
    x,y=lag[m].astype(float),qn[m].astype(float)
    if len(x)<4:
        return np.nan,np.nan,np.nan
    if p0 is None:
        p0=[float(x[len(x)//2]),0.5]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt,_=curve_fit(kww,x,y,p0=p0,bounds=([1e-3,0.01],[1e8,1.5]),maxfev=5000)
        yhat=kww(x,*popt)
        ss_res=np.sum((y-yhat)**2)
        ss_tot=np.sum((y-y.mean())**2)
        r2=1-ss_res/ss_tot if ss_tot>0 else np.nan
        return float(popt[0]),float(popt[1]),float(r2)
    except:
        return np.nan,np.nan,np.nan
def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--csv",required=True)
    p.add_argument("--outdir",default="analysis/extended")
    p.add_argument("--n_bootstrap",type=int,default=1000)
    return p.parse_args()
def load_and_aggregate(csv_path):
    df=pd.read_csv(csv_path)
    gcols=["ensemble","epsilon","tw","lag"]
    agg_spec={"Q":("Q_mean","mean"),"chi4":("chi4","mean"),"D4":("D4_mean","mean")}
    if "MSD_mean" in df.columns:
        agg_spec["MSD"]=("MSD_mean","mean")
    if "alpha2_mean" in df.columns:
        agg_spec["alpha2"]=("alpha2_mean","mean")
    if "energy_mean" in df.columns:
        agg_spec["energy"]=("energy_mean","mean")
    if "Rg_mean" in df.columns:
        agg_spec["Rg"]=("Rg_mean","mean")
    if "n_contacts_mean" in df.columns:
        agg_spec["n_contacts"]=("n_contacts_mean","mean")
    agg=df.groupby(gcols).agg(**agg_spec).reset_index()
    key=["ensemble","epsilon","tw"]
    if 0 in agg["lag"].values:
        q0=agg[agg["lag"]==0][key+["Q"]].rename(columns={"Q":"Q0"})
    else:
        ref=agg.loc[agg.groupby(key)["lag"].idxmin()][key+["Q"]].rename(columns={"Q":"Q0"})
        q0=ref
    agg=agg.merge(q0,on=key,how="left")
    agg["Q_norm"]=np.where(agg["Q0"]>0,agg["Q"]/agg["Q0"],np.nan)
    return df,agg
def analysis_kww(agg,outdir):
    rows=[]
    for (ens,eps),g in agg.groupby(["ensemble","epsilon"]):
        for tw,s in g.groupby("tw"):
            s=s.sort_values("lag")
            lag=s["lag"].to_numpy(float)
            qn=s["Q_norm"].to_numpy(float)
            tau_k,beta_k,r2_k=fit_kww(lag,qn)
            rows.append({"ensemble":ens,"epsilon":eps,"tw":tw,
                         "tau_kww":tau_k,"beta_kww":beta_k,"r2_kww":r2_k})
    kdf=pd.DataFrame(rows)
    kdf.to_csv(os.path.join(outdir,"kww_fits.csv"),index=False)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    epsilons=sorted(kdf["epsilon"].unique())
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            sub=kdf[(kdf["ensemble"]==ens)&(kdf["epsilon"]==eps)&(kdf["tw"]>0)]
            sub=sub.sort_values("tw")
            if len(sub)==0: continue
            ax.plot(sub["tw"],sub["beta_kww"],color=c,marker=m,linestyle=ls,label=ens)
        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)")
        ax.set_ylabel(r"$\beta_{KWW}$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
        ax.set_ylim(0,1.2)
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_kww_beta_vs_tw.png"))
    plt.savefig(os.path.join(outdir,"Fig_kww_beta_vs_tw.pdf"))
    plt.close()
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            sub=kdf[(kdf["ensemble"]==ens)&(kdf["epsilon"]==eps)&(kdf["tw"]>0)]
            sub=sub[np.isfinite(sub["tau_kww"])].sort_values("tw")
            if len(sub)==0: continue
            ax.plot(sub["tw"],sub["tau_kww"],color=c,marker=m,linestyle=ls,label=ens)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)")
        ax.set_ylabel(r"$\tau_{KWW}$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_kww_tau_vs_tw.png"))
    plt.savefig(os.path.join(outdir,"Fig_kww_tau_vs_tw.pdf"))
    plt.close()
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for i,eps in enumerate(epsilons):
        for ens,(c,m,ls) in ENS_STYLE.items():
            for tw in sorted(agg["tw"].unique()):
                s=agg[(agg["ensemble"]==ens)&(agg["epsilon"]==eps)&(agg["tw"]==tw)].sort_values("lag")
                lag=s["lag"].to_numpy(float)
                qn=s["Q_norm"].to_numpy(float)
                mp=lag>0
                if not np.any(mp): continue
                row_idx=0 if ens=="iid" else 1
                ax=axes[row_idx,i]
                alpha=0.3+0.7*(list(sorted(agg["tw"].unique())).index(tw)/max(1,len(agg["tw"].unique())-1))
                ax.plot(lag[mp],qn[mp],color=c,alpha=alpha,marker=m,markersize=3,linestyle="-",linewidth=0.8)
                kr=kdf[(kdf["ensemble"]==ens)&(kdf["epsilon"]==eps)&(kdf["tw"]==tw)]
                if len(kr)>0 and np.isfinite(kr.iloc[0]["tau_kww"]):
                    tau_k=kr.iloc[0]["tau_kww"]
                    beta_k=kr.iloc[0]["beta_kww"]
                    lf=np.logspace(np.log10(max(1,lag[mp].min())),np.log10(lag[mp].max()),100)
                    ax.plot(lf,kww(lf,tau_k,beta_k),color="black",alpha=alpha*0.7,linewidth=0.6,linestyle="--")
        for row in range(2):
            axes[row,i].set_xscale("log")
            axes[row,i].set_xlabel("lag")
            axes[row,i].set_ylabel(r"$\tilde{Q}$")
            axes[row,i].set_title(rf"{'iid' if row==0 else 'corr'}, $\epsilon={eps:.0f}$")
            axes[row,i].set_ylim(0,1.1)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_kww_fit_curves.png"))
    plt.close()
    return kdf
def analysis_msd(agg,outdir):
    if "MSD" not in agg.columns:
        print("  MSD not in data, skipping")
        return None
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    epsilons=sorted(agg["epsilon"].unique())
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            for tw in sorted(agg["tw"].unique()):
                s=agg[(agg["ensemble"]==ens)&(agg["epsilon"]==eps)&(agg["tw"]==tw)].sort_values("lag")
                lag=s["lag"].to_numpy(float)
                msd=s["MSD"].to_numpy(float)
                mp=lag>0
                if not np.any(mp): continue
                alpha=0.3+0.7*(list(sorted(agg["tw"].unique())).index(tw)/max(1,len(agg["tw"].unique())-1))
                ax.plot(lag[mp],msd[mp],color=c,alpha=alpha,marker=m,markersize=2,linestyle=ls,linewidth=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("lag"); ax.set_ylabel("MSD")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
        x_ref=np.logspace(0,4,50)
        ax.plot(x_ref,0.01*x_ref,"k:",alpha=0.3,label=r"$\sim t$")
        ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_MSD_vs_lag.png"))
    plt.savefig(os.path.join(outdir,"Fig_MSD_vs_lag.pdf"))
    plt.close()
    rows=[]
    for (ens,eps),g in agg.groupby(["ensemble","epsilon"]):
        for tw,s in g.groupby("tw"):
            s=s.sort_values("lag")
            lag=s["lag"].to_numpy(float)
            msd=s["MSD"].to_numpy(float)
            mp=(lag>0)&(msd>0)
            if np.sum(mp)<3: continue
            sl,ic,rv,_,_=linregress(np.log10(lag[mp]),np.log10(msd[mp]))
            rows.append({"ensemble":ens,"epsilon":eps,"tw":tw,"msd_exponent":sl,"msd_r2":rv**2})
    msd_df=pd.DataFrame(rows)
    msd_df.to_csv(os.path.join(outdir,"msd_exponents.csv"),index=False)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            sub=msd_df[(msd_df["ensemble"]==ens)&(msd_df["epsilon"]==eps)&(msd_df["tw"]>0)]
            sub=sub.sort_values("tw")
            if len(sub)==0: continue
            ax.plot(sub["tw"],sub["msd_exponent"],color=c,marker=m,linestyle=ls,label=ens)
        ax.axhline(y=1.0,color="gray",linestyle="--",alpha=0.5,label="diffusive")
        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$"); ax.set_ylabel(r"MSD exponent $\alpha$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_MSD_exponent_vs_tw.png"))
    plt.savefig(os.path.join(outdir,"Fig_MSD_exponent_vs_tw.pdf"))
    plt.close()
    return msd_df
def analysis_alpha2(agg,outdir):
    if "alpha2" not in agg.columns:
        print("  alpha2 not in data, skipping")
        return None
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    epsilons=sorted(agg["epsilon"].unique())
    rows=[]
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            for tw in sorted(agg["tw"].unique()):
                s=agg[(agg["ensemble"]==ens)&(agg["epsilon"]==eps)&(agg["tw"]==tw)].sort_values("lag")
                lag=s["lag"].to_numpy(float)
                a2=s["alpha2"].to_numpy(float)
                mp=lag>0
                if not np.any(mp): continue
                alpha=0.3+0.7*(list(sorted(agg["tw"].unique())).index(tw)/max(1,len(agg["tw"].unique())-1))
                ax.plot(lag[mp],a2[mp],color=c,alpha=alpha,marker=m,markersize=2,linestyle=ls,linewidth=0.8)
                pk_idx=np.argmax(a2[mp])
                rows.append({"ensemble":ens,"epsilon":eps,"tw":tw,
                             "alpha2_peak":float(a2[mp][pk_idx]),"alpha2_tpeak":float(lag[mp][pk_idx])})
        ax.set_xscale("log")
        ax.set_xlabel("lag"); ax.set_ylabel(r"$\alpha_2(t_w,t)$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_alpha2_vs_lag.png"))
    plt.savefig(os.path.join(outdir,"Fig_alpha2_vs_lag.pdf"))
    plt.close()
    a2df=pd.DataFrame(rows)
    a2df.to_csv(os.path.join(outdir,"alpha2_peaks.csv"),index=False)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            sub=a2df[(a2df["ensemble"]==ens)&(a2df["epsilon"]==eps)&(a2df["tw"]>0)]
            sub=sub.sort_values("tw")
            if len(sub)==0: continue
            ax.plot(sub["tw"],sub["alpha2_peak"],color=c,marker=m,linestyle=ls,label=ens)
        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$"); ax.set_ylabel(r"$\alpha_2^*$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_alpha2_peak_vs_tw.png"))
    plt.savefig(os.path.join(outdir,"Fig_alpha2_peak_vs_tw.pdf"))
    plt.close()
    return a2df
def analysis_energy_rg(agg,outdir):
    has_e="energy" in agg.columns
    has_rg="Rg" in agg.columns
    has_nc="n_contacts" in agg.columns
    if not (has_e or has_rg or has_nc):
        print("  No energy/Rg/n_contacts in data, skipping")
        return None
    epsilons=sorted(agg["epsilon"].unique())
    ncols=sum([has_e,has_rg,has_nc])
    fig,axes=plt.subplots(1,ncols,figsize=(5*ncols,4.5))
    if ncols==1: axes=[axes]
    col_idx=0
    snap=agg[agg["lag"]==agg["lag"].min()].copy()
    if has_e:
        ax=axes[col_idx]; col_idx+=1
        for ens,(c,m,ls) in ENS_STYLE.items():
            for eps in epsilons:
                sub=snap[(snap["ensemble"]==ens)&(snap["epsilon"]==eps)].sort_values("tw")
                if len(sub)==0: continue
                ax.plot(sub["tw"],sub["energy"],color=EPS_COLORS.get(eps,c),
                        marker=m,linestyle=ls,label=f"{ens} $\\epsilon$={eps:.0f}")
        ax.set_xscale("log"); ax.set_xlabel(r"$t_w$"); ax.set_ylabel("Energy")
        ax.set_title("Energy at snapshot"); ax.legend(fontsize=7)
    if has_rg:
        ax=axes[col_idx]; col_idx+=1
        for ens,(c,m,ls) in ENS_STYLE.items():
            for eps in epsilons:
                sub=snap[(snap["ensemble"]==ens)&(snap["epsilon"]==eps)].sort_values("tw")
                if len(sub)==0: continue
                ax.plot(sub["tw"],sub["Rg"],color=EPS_COLORS.get(eps,c),
                        marker=m,linestyle=ls,label=f"{ens} $\\epsilon$={eps:.0f}")
        ax.set_xscale("log"); ax.set_xlabel(r"$t_w$"); ax.set_ylabel(r"$R_g$")
        ax.set_title("Radius of gyration"); ax.legend(fontsize=7)
    if has_nc:
        ax=axes[col_idx]; col_idx+=1
        for ens,(c,m,ls) in ENS_STYLE.items():
            for eps in epsilons:
                sub=snap[(snap["ensemble"]==ens)&(snap["epsilon"]==eps)].sort_values("tw")
                if len(sub)==0: continue
                ax.plot(sub["tw"],sub["n_contacts"],color=EPS_COLORS.get(eps,c),
                        marker=m,linestyle=ls,label=f"{ens} $\\epsilon$={eps:.0f}")
        ax.set_xscale("log"); ax.set_xlabel(r"$t_w$"); ax.set_ylabel("Contact count")
        ax.set_title("Number of contacts"); ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_energy_Rg_contacts_vs_tw.png"))
    plt.savefig(os.path.join(outdir,"Fig_energy_Rg_contacts_vs_tw.pdf"))
    plt.close()
    return snap
def analysis_rel_chi4_diff(agg,outdir):
    epsilons=sorted(agg["epsilon"].unique())
    rows=[]
    for (ens,eps),g in agg.groupby(["ensemble","epsilon"]):
        for tw,s in g.groupby("tw"):
            s=s.sort_values("lag")
            lag=s["lag"].to_numpy(float)
            chi4=s["chi4"].to_numpy(float)
            mp=lag>0
            if not np.any(mp): continue
            pk=np.argmax(chi4[mp])
            rows.append({"ensemble":ens,"epsilon":eps,"tw":tw,
                         "chi4_star":float(chi4[mp][pk]),"t_star":float(lag[mp][pk])})
    ts=pd.DataFrame(rows)
    sep_rows=[]
    for eps in epsilons:
        if eps==0: continue
        iid=ts[(ts["ensemble"]=="iid")&(ts["epsilon"]==eps)].set_index("tw")
        corr=ts[(ts["ensemble"]=="correlated")&(ts["epsilon"]==eps)].set_index("tw")
        common=iid.index.intersection(corr.index)
        for tw in common:
            c4i=iid.loc[tw,"chi4_star"]
            c4c=corr.loc[tw,"chi4_star"]
            if c4c>0:
                sep_rows.append({"epsilon":eps,"tw":tw,"rel_diff_pct":(c4i-c4c)/c4c*100})
    sep=pd.DataFrame(sep_rows)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for eps,c in [(3.0,"tab:orange"),(6.0,"tab:purple")]:
        sub=sep[sep["epsilon"]==eps].sort_values("tw")
        if len(sub)==0: continue
        axes[0].plot(sub["tw"],sub["rel_diff_pct"],marker="o",color=c,label=rf"$\epsilon={eps:.0f}$")
    axes[0].axhline(0,color="gray",linestyle="--",alpha=0.5)
    axes[0].set_xscale("log"); axes[0].set_xlabel(r"$t_w$")
    axes[0].set_ylabel(r"$\Delta\chi_4^*/\chi_{4,corr}^*$ [%]")
    axes[0].set_title("Relative heterogeneity enhancement (iid vs corr)")
    axes[0].legend()
    for i,eps in enumerate(epsilons):
        ax=axes[1]
        for ens,(c,m,ls) in ENS_STYLE.items():
            sub=ts[(ts["ensemble"]==ens)&(ts["epsilon"]==eps)&(ts["tw"]>0)].sort_values("tw")
            if len(sub)==0: continue
            ax.plot(sub["tw"],sub["chi4_star"],color=EPS_COLORS.get(eps,c),
                    marker=m,linestyle=ls,label=f"{ens} $\\epsilon$={eps:.0f}")
    axes[1].set_xscale("log"); axes[1].set_xlabel(r"$t_w$")
    axes[1].set_ylabel(r"$\chi_4^*$"); axes[1].set_title(r"$\chi_4^*$ vs age")
    axes[1].legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_chi4_star_and_separation.png"))
    plt.savefig(os.path.join(outdir,"Fig_chi4_star_and_separation.pdf"))
    plt.close()
    return ts,sep
def analysis_d4_scaling(agg,outdir):
    epsilons=sorted(agg["epsilon"].unique())
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for i,eps in enumerate(epsilons):
        ax=axes[i]
        for ens,(c,m,ls) in ENS_STYLE.items():
            for tw in sorted(agg["tw"].unique()):
                s=agg[(agg["ensemble"]==ens)&(agg["epsilon"]==eps)&(agg["tw"]==tw)].sort_values("lag")
                lag=s["lag"].to_numpy(float)
                d4=s["D4"].to_numpy(float)
                mp=(lag>0)&(d4>0)
                if not np.any(mp): continue
                alpha=0.3+0.7*(list(sorted(agg["tw"].unique())).index(tw)/max(1,len(agg["tw"].unique())-1))
                ax.plot(lag[mp],d4[mp],color=c,alpha=alpha,marker=m,markersize=2,linestyle=ls,linewidth=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("lag"); ax.set_ylabel(r"$D_4(t_w,t)$")
        ax.set_title(rf"$\epsilon={eps:.0f}$")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir,"Fig_D4_vs_lag.png"))
    plt.savefig(os.path.join(outdir,"Fig_D4_vs_lag.pdf"))
    plt.close()
def print_findings(kdf,msd_df,a2df,ts,sep):
    print("\n"+"="*70)
    print("EXTENDED ANALYSIS FINDINGS")
    print("="*70)
    print("\n--- KWW Stretched Exponential Fits ---")
    if kdf is not None:
        for eps in sorted(kdf["epsilon"].unique()):
            for ens in ["iid","correlated"]:
                sub=kdf[(kdf["ensemble"]==ens)&(kdf["epsilon"]==eps)&(kdf["tw"]>0)]
                sub=sub[np.isfinite(sub["beta_kww"])]
                if len(sub)==0: continue
                b_mean=sub["beta_kww"].mean()
                b_late=sub[sub["tw"]==sub["tw"].max()]["beta_kww"]
                b_late_v=float(b_late.iloc[0]) if len(b_late)>0 else np.nan
                print(f"  eps={eps:.0f} {ens}: <beta>={b_mean:.3f}, beta(tw_max)={b_late_v:.3f}")
    print("\n--- MSD Subdiffusion Exponents ---")
    if msd_df is not None:
        for eps in sorted(msd_df["epsilon"].unique()):
            for ens in ["iid","correlated"]:
                sub=msd_df[(msd_df["ensemble"]==ens)&(msd_df["epsilon"]==eps)&(msd_df["tw"]>0)]
                if len(sub)==0: continue
                a_late=sub[sub["tw"]==sub["tw"].max()]["msd_exponent"]
                if len(a_late)>0:
                    print(f"  eps={eps:.0f} {ens}: alpha(tw_max)={float(a_late.iloc[0]):.3f}")
    print("\n--- Non-Gaussianity Peak ---")
    if a2df is not None:
        for eps in sorted(a2df["epsilon"].unique()):
            for ens in ["iid","correlated"]:
                sub=a2df[(a2df["ensemble"]==ens)&(a2df["epsilon"]==eps)&(a2df["tw"]>0)]
                if len(sub)==0: continue
                pk_late=sub[sub["tw"]==sub["tw"].max()]["alpha2_peak"]
                if len(pk_late)>0:
                    print(f"  eps={eps:.0f} {ens}: alpha2_peak(tw_max)={float(pk_late.iloc[0]):.3f}")
    print("\n--- Ensemble Separation (chi4*) ---")
    if sep is not None and len(sep)>0:
        for eps in sorted(sep["epsilon"].unique()):
            sub=sep[sep["epsilon"]==eps]
            late=sub[sub["tw"]==sub["tw"].max()]
            if len(late)>0:
                print(f"  eps={eps:.0f}: late-age (iid-corr)/corr = {float(late.iloc[0]['rel_diff_pct']):+.1f}%")
    print("\n"+"="*70)
def main():
    args=parse_args()
    os.makedirs(args.outdir,exist_ok=True)
    print(f"Loading {args.csv}")
    df,agg=load_and_aggregate(args.csv)
    print(f"  {len(df)} raw rows, {len(agg)} aggregated rows")
    print(f"  Columns: {list(df.columns)}")
    print("\n1. KWW stretched exponential fitting...")
    kdf=analysis_kww(agg,args.outdir)
    print("   Saved kww_fits.csv + 3 figures")
    print("\n2. MSD analysis...")
    msd_df=analysis_msd(agg,args.outdir)
    print("\n3. Non-Gaussianity parameter...")
    a2df=analysis_alpha2(agg,args.outdir)
    print("\n4. Energy / Rg / contact count...")
    analysis_energy_rg(agg,args.outdir)
    print("\n5. D4 structural displacement...")
    analysis_d4_scaling(agg,args.outdir)
    print("\n6. chi4* separation analysis...")
    ts,sep=analysis_rel_chi4_diff(agg,args.outdir)
    print_findings(kdf,msd_df,a2df,ts,sep)
    print(f"\nAll outputs saved to {args.outdir}")
if __name__=="__main__":
    main()
