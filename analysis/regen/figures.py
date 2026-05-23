from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator,MaxNLocator
from .observables import (rpa_predictor_finite,rpa_predictor_asymptotic,fit_loglog,
                          aging_slope_bootstrap,collapse_predictor_comparison)

KAPPA_COLORS={0.0:"#1f77b4",0.15:"#3a8fbe",0.3:"#69a8c5",0.45:"#9ec0c9",
              0.5:"#cf9a78",0.6:"#e6864c",0.75:"#ec6e2a",0.8:"#e85a18",
              0.9:"#dc3a14",1.0:"#b2182b"}
SEQ_COLORS={"random":"#1f77b4","correlated":"#d62728","block":"#2ca02c","alternating":"#9467bd"}
COND_COLORS={"baseline":"#386cb0","dense":"#f0027f","soft":"#4daf4a","short_chain":"#984ea3"}
COND_MARKERS={"baseline":"o","dense":"s","soft":"D","short_chain":"^"}

def apply_style():
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":8,
        "axes.labelsize":8,"axes.titlesize":9,"axes.titleweight":"bold",
        "xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,
        "axes.linewidth":0.8,"xtick.major.width":0.8,"ytick.major.width":0.8,
        "xtick.major.size":3,"ytick.major.size":3,
        "axes.grid":True,"grid.alpha":0.25,"grid.linewidth":0.4,
        "savefig.dpi":400,"savefig.bbox":"tight","savefig.pad_inches":0.05,
        "figure.facecolor":"white","savefig.facecolor":"white",
    })

def kappa_color(k):
    try:return KAPPA_COLORS[float(round(k,2))]
    except Exception:return "#888888"

def panel_label(ax,letter,x=-0.18,y=1.05):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top",ha="left")

def fig1_kappa_robustness(df_by_condition:dict[str,pd.DataFrame],out:Path,N=40,b=1.0):
    """Fig 1. The kappa effect is not an artifact of one parameter setting.
    Panels: (A) C vs kappa across 4 model variants with RPA prediction overlaid for baseline.
    (B) raw S_AA(k*) showing why we use C=S_AA-S_AB. (C) Domain length xi vs kappa.
    (D) Percent energy drop relative to kappa=0 within each condition."""
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.6))
    A,B,C,D=axes[0,0],axes[0,1],axes[1,0],axes[1,1]
    for cond,d in df_by_condition.items():
        g=d.groupby("kappa")
        mu=g["C"].mean();se=g["C"].sem();ks=mu.index.to_numpy()
        c=COND_COLORS.get(cond,"#666"); mk=COND_MARKERS.get(cond,"o")
        A.errorbar(ks,mu.values,yerr=se.values,fmt=mk+"-",color=c,markersize=4,
                   linewidth=1.2,capsize=2.5,label=f"{cond} (n={int(d.groupby('kappa').size().min())})")
        mu_s=g["S_AA_peak"].mean();se_s=g["S_AA_peak"].sem()
        B.errorbar(ks,mu_s.values,yerr=se_s.values,fmt=mk+"-",color=c,markersize=4,linewidth=1.2,capsize=2.5)
        mu_xi=g["xi"].mean();se_xi=g["xi"].sem()
        C.errorbar(ks,mu_xi.values,yerr=se_xi.values,fmt=mk+"-",color=c,markersize=4,linewidth=1.2,capsize=2.5)
        mu_E=g["E"].mean()
        E0=mu_E[0.0] if 0.0 in mu_E.index else mu_E.iloc[0]
        dE=100.0*(mu_E.values-E0)/abs(E0) if E0!=0 else mu_E.values*0
        D.plot(ks,dE,mk+"-",color=c,markersize=4,linewidth=1.2)
    base=df_by_condition.get("baseline")
    if base is not None and len(base):
        pi_base=float(base["pi"].dropna().iloc[0]) if base["pi"].notna().any() else 0.99
        N_base=int(base["N"].dropna().iloc[0]) if base["N"].notna().any() else N
        ks_th=np.linspace(0,1,40)
        pred=np.array([rpa_predictor_finite(k,pi_base,N=N_base,k=0.0,b=b) for k in ks_th])
        c0=float(base.groupby("kappa")["C"].mean().get(0.0,np.nan))
        if np.isfinite(c0):
            scale=(base.groupby("kappa")["C"].mean().max()-c0)/(pred.max()-1.0) if pred.max()>1 else 1
            A.plot(ks_th,c0+scale*(pred-1.0),"--",color="black",linewidth=1.0,alpha=0.7,
                   label=fr"RPA $1+2\kappa^2\lambda/(1-\lambda)$, $\pi={pi_base:.3g}$")
    A.set_xlabel(r"sequence-correlation amplitude $\kappa$");A.set_ylabel(r"contrast $C=\max_k[S_{AA}-S_{AB}]$")
    A.set_title("A. Contrast vs $\\kappa$ across model variants",loc="left")
    A.legend(loc="upper left",frameon=True,framealpha=0.9)
    panel_label(A,"A")
    B.set_xlabel(r"$\kappa$");B.set_ylabel(r"raw $S_{AA}(k_\star)$")
    B.set_title("B. Raw $S_{AA}$ tracks $\\kappa$ only in dense melt",loc="left");panel_label(B,"B")
    C.set_xlabel(r"$\kappa$");C.set_ylabel(r"domain length $\xi=2\pi/k_\star\ [\sigma]$")
    C.set_title("C. Domain length vs $\\kappa$",loc="left");panel_label(C,"C")
    D.set_xlabel(r"$\kappa$");D.set_ylabel(r"$\Delta E$ vs $\kappa{=}0$ (%)")
    D.set_title("D. Energy drop confirms ordering",loc="left");panel_label(D,"D")
    D.axhline(0,color="grey",linewidth=0.4,linestyle=":")
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig2_blockiness_corrected(df_pi_kappa:pd.DataFrame,df_xi_at_pi:pd.DataFrame,out:Path,N=40,b=1.0):
    """Fig 2. Persistence and amplitude combine multiplicatively (corrected RPA).
    Panels: (A) heatmap C across (pi, kappa). (B) Excess contrast vs corrected RPA
    predictor with old-vs-new R^2 comparison. (C) Domain length plateau at pi=0.99."""
    apply_style()
    fig=plt.figure(figsize=(8.5,3.2))
    gs=fig.add_gridspec(1,3,width_ratios=[1.05,1.6,0.95],wspace=0.45)
    A=fig.add_subplot(gs[0,0]);B=fig.add_subplot(gs[0,1]);C=fig.add_subplot(gs[0,2])
    pivot=df_pi_kappa.groupby(["pi","kappa"])["C"].mean().unstack("kappa")
    pivot=pivot[pivot.index<1.0].sort_index()
    im=A.imshow(pivot.values,origin="lower",aspect="auto",
                extent=[pivot.columns.min(),pivot.columns.max(),pivot.index.min(),pivot.index.max()],
                cmap="viridis")
    A.set_xlabel(r"$\kappa$");A.set_ylabel(r"$\pi$");A.set_title("A. C across $(\\pi,\\kappa)$",loc="left")
    cb=fig.colorbar(im,ax=A,shrink=0.85,pad=0.02);cb.set_label("C",fontsize=7);cb.ax.tick_params(labelsize=6)
    panel_label(A,"A")
    cmp=collapse_predictor_comparison(df_pi_kappa,N=N,b=b)
    if len(cmp):
        x_corr=cmp["pred_corrected_minus_1"].values;y=cmp["excess"].values
        x_asy=cmp["pred_asymptotic"].values
        s_c,b_c,r2_c,lo_c,hi_c=fit_loglog(x_corr,y)
        s_a,b_a,r2_a,lo_a,hi_a=fit_loglog(x_asy,y)
        cvals=cmp["kappa"].values
        sc=B.scatter(x_corr,y,c=cvals,cmap="plasma",s=22,edgecolor="black",linewidth=0.3,label="corrected predictor")
        if np.isfinite(s_c):
            xs=np.logspace(np.log10(max(np.nanmin(x_corr),1e-3)),np.log10(np.nanmax(x_corr)),50)
            B.plot(xs,10**(s_c*np.log10(xs)+b_c),"-",color="black",linewidth=1.0,
                   label=fr"corrected: slope={s_c:.2f} [{lo_c:.2f},{hi_c:.2f}], $R^2$={r2_c:.2f}")
        if np.isfinite(s_a):
            xs=np.logspace(np.log10(max(np.nanmin(x_asy),1e-3)),np.log10(np.nanmax(x_asy)),50)
            B.plot(xs,10**(s_a*np.log10(xs)+b_a),"--",color="grey",linewidth=1.0,
                   label=fr"old asymptotic: slope={s_a:.2f} [{lo_a:.2f},{hi_a:.2f}], $R^2$={r2_a:.2f}")
        B.set_xscale("log");B.set_yscale("log")
        cb2=fig.colorbar(sc,ax=B,shrink=0.85,pad=0.02);cb2.set_label(r"$\kappa$",fontsize=7);cb2.ax.tick_params(labelsize=6)
    B.set_xlabel("RPA predictor (corrected: $2\\kappa^2\\Sigma(1-\\ell/N)\\lambda^\\ell$; dashed: $\\kappa^2/[2(1-\\pi)]$)")
    B.set_ylabel(r"excess contrast $C-C_{\kappa=0}$")
    B.set_title("B. Corrected predictor collapses the data",loc="left")
    B.legend(loc="upper left",frameon=True,framealpha=0.9,fontsize=6)
    panel_label(B,"B")
    if len(df_xi_at_pi):
        g=df_xi_at_pi.groupby("kappa")["xi"]
        mu=g.mean();se=g.sem()
        C.errorbar(mu.index,mu.values,yerr=se.values,fmt="o-",color="#1b9e77",markersize=4,capsize=2.5)
    C.set_xlabel(r"$\kappa$");C.set_ylabel(r"$\xi\ [\sigma]$")
    C.set_title(r"C. Domain length saturates at $\pi=0.99$",loc="left");panel_label(C,"C")
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig3_composition(df_fA_kappa:pd.DataFrame,df_at_kappa07:pd.DataFrame,df_high_pi:pd.DataFrame,out:Path):
    """Fig 3. Composition asymmetry under unified lambda generator.
    The discontinuity at f_A=0.5 from the previous code path is gone: the (f_A, kappa)
    heatmap should now be smooth across f_A=0.5."""
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.6))
    A,B,C,D=axes[0,0],axes[0,1],axes[1,0],axes[1,1]
    pivot=df_fA_kappa.groupby(["f_A","kappa"])["C"].mean().unstack("kappa")
    im=A.imshow(pivot.values,origin="lower",aspect="auto",
                extent=[pivot.columns.min(),pivot.columns.max(),pivot.index.min(),pivot.index.max()],
                cmap="viridis")
    A.set_xlabel(r"$\kappa$");A.set_ylabel(r"$f_A$");A.set_title("A. C across $(f_A,\\kappa)$ at $\\pi=0.95$",loc="left")
    cb=fig.colorbar(im,ax=A,shrink=0.85,pad=0.02);cb.set_label("C",fontsize=7);cb.ax.tick_params(labelsize=6)
    panel_label(A,"A")
    if len(df_at_kappa07):
        g=df_at_kappa07.groupby("f_A")["C"]
        mu=g.mean();se=g.sem();fA=mu.index.to_numpy()
        B.errorbar(fA,mu.values,yerr=se.values,fmt="o-",color="#d95f02",markersize=4,linewidth=1.2,capsize=2.5,
                   label=fr"$\kappa=0.7$, $\pi=0.95$ (n={int(df_at_kappa07.groupby('f_A').size().min())})")
        fA_th=np.linspace(0.1,0.9,80)
        base=4.0*fA_th*(1.0-fA_th)
        scale=float(np.nanmax(mu.values))/float(np.nanmax(base)) if np.nanmax(base)>0 else 1
        B.plot(fA_th,base*scale,"--",color="grey",linewidth=1.0,label=r"$4f_A(1-f_A)$ baseline (scaled)")
    B.set_xlabel(r"$f_A$");B.set_ylabel("contrast C")
    B.set_title("B. C vs $f_A$ at fixed $\\kappa=0.7$",loc="left")
    B.legend(loc="upper center",frameon=True,framealpha=0.9);panel_label(B,"B")
    pivot=df_fA_kappa.groupby(["f_A","kappa"])["C"].mean().unstack("kappa")
    if 0.0 in pivot.columns and 1.0 in pivot.columns:
        ratio=pivot[1.0]/pivot[0.0]
        C.plot(ratio.index,ratio.values,"o-",color="#7570b3",markersize=5,linewidth=1.4)
    C.set_xlabel(r"$f_A$");C.set_ylabel(r"$C_{\kappa=1}/C_{\kappa=0}$")
    C.set_title("C. $\\kappa$ amplification peaks at minority A",loc="left");panel_label(C,"C")
    if len(df_high_pi):
        g=df_high_pi.groupby("f_A")["C"]
        mu=g.mean();se=g.sem()
        D.errorbar(mu.index,mu.values,yerr=se.values,fmt="s-",color="#1b9e77",markersize=4,linewidth=1.2,capsize=2.5,
                   label=fr"$\pi=0.995$, $\epsilon_{{AB}}=0.05$")
    D.set_xlabel(r"$f_A$");D.set_ylabel("contrast C")
    D.set_title("D. High-persistence reproduction",loc="left")
    D.legend(loc="best",frameon=True,framealpha=0.9);panel_label(D,"D")
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig4_interactions_density(df_eps_moderate:pd.DataFrame,df_eps_high_pi:pd.DataFrame,
                              df_density:pd.DataFrame,density_slice_path:Path|None,out:Path):
    """Fig 4. Cross-attraction scan with corrected WCA core (probes incompatibility, not
    cross-pair excluded volume). At eps_AB=10^-4 cross pairs now still have WCA repulsion;
    only the attractive tail is suppressed."""
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.6))
    A,B,C,D=axes[0,0],axes[0,1],axes[1,0],axes[1,1]
    if len(df_eps_moderate):
        for seq,grp in df_eps_moderate.groupby("sequence"):
            g=grp.groupby("eps_AB")["C"]
            mu=g.mean();se=g.sem()
            A.errorbar(mu.index,mu.values,yerr=se.values,fmt="o-",
                       color=SEQ_COLORS.get(seq,"#444"),markersize=4,linewidth=1.2,capsize=2.5,label=seq)
        A.set_xscale("log")
    A.set_xlabel(r"$\epsilon_{AB}/\epsilon_{AA}$");A.set_ylabel("contrast C")
    A.set_title(r"A. $\epsilon_{AB}$ scan at $(\kappa,\pi)=(0.7,0.9)$",loc="left")
    A.legend(loc="best",frameon=True,framealpha=0.9);panel_label(A,"A")
    if len(df_eps_high_pi):
        g=df_eps_high_pi.groupby("eps_AB")["C"]
        mu=g.mean();se=g.sem()
        B.errorbar(mu.index,mu.values,yerr=se.values,fmt="s-",color="#d62728",markersize=4,linewidth=1.2,capsize=2.5,
                   label=fr"$\kappa=1$, $\pi=0.99$, WCA core fixed")
        B.set_xscale("log")
    B.set_xlabel(r"$\epsilon_{AB}/\epsilon_{AA}$");B.set_ylabel("contrast C")
    B.set_title("B. High-persistence: cross attraction is amplifier",loc="left")
    B.legend(loc="best",frameon=True,framealpha=0.9);panel_label(B,"B")
    if len(df_density):
        g1=df_density[df_density["kappa"]==1.0].groupby("rho")["C"]
        g0=df_density[df_density["kappa"]==0.0].groupby("rho")["C"]
        mu1=g1.mean();se1=g1.sem();mu0=g0.mean();se0=g0.sem()
        C.errorbar(mu1.index,mu1.values,yerr=se1.values,fmt="o-",color="#d62728",markersize=4,linewidth=1.2,capsize=2.5,label=r"$\kappa=1$")
        C.errorbar(mu0.index,mu0.values,yerr=se0.values,fmt="s-",color="#1f77b4",markersize=4,linewidth=1.2,capsize=2.5,label=r"$\kappa=0$")
        s,b,r2,lo,hi=fit_loglog(mu1.index.to_numpy(),mu1.values)
        if np.isfinite(s):
            xs=np.logspace(np.log10(mu1.index.min()),np.log10(mu1.index.max()),60)
            C.plot(xs,10**(s*np.log10(xs)+b),"--",color="black",linewidth=1.0,
                   label=fr"$\kappa=1$ fit: slope={s:.2f} [{lo:.2f},{hi:.2f}], $R^2$={r2:.2f}")
        C.set_xscale("log");C.set_yscale("log")
    C.set_xlabel(r"density $\rho\ [\sigma^{-3}]$");C.set_ylabel("contrast C")
    C.set_title("C. Density scaling",loc="left");C.legend(loc="best",frameon=True,framealpha=0.9);panel_label(C,"C")
    if density_slice_path is not None and Path(density_slice_path).exists():
        img=plt.imread(density_slice_path)
        D.imshow(img);D.axis("off")
    else:
        D.text(0.5,0.5,"density slice $\\kappa\\in\\{0,0.5,1\\}$\n(rendered when fig5_aging trajectories sync)",
               ha="center",va="center",transform=D.transAxes,fontsize=8,color="grey")
        D.axis("off")
    D.set_title(r"D. Real-space $(\phi_A-\phi_B)/(\phi_A+\phi_B)$ at $\rho=1.17$",loc="left");panel_label(D,"D")
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig5_architecture(df_block:pd.DataFrame,df_diblock_N:pd.DataFrame,df_ratio_N:pd.DataFrame,out:Path):
    """Fig 5. Deterministic block architecture amplifies the kappa effect.
    Panels: (A) C vs deterministic block length, (B) di-block C ~ N^alpha, (C) amplification ratio."""
    apply_style()
    fig,axes=plt.subplots(1,3,figsize=(9.0,3.2))
    A,B,C=axes
    if len(df_block):
        g=df_block.groupby("block_length")["C"] if "block_length" in df_block.columns else df_block.groupby("N")["C"]
        mu=g.mean();se=g.sem()
        A.errorbar(mu.index,mu.values,yerr=se.values,fmt="o-",color="#2ca02c",markersize=4,linewidth=1.2,capsize=2.5)
    A.set_xlabel(r"block length");A.set_ylabel("contrast C")
    A.set_title(r"A. Deterministic block scan at $N=40$",loc="left");panel_label(A,"A")
    if len(df_diblock_N):
        g=df_diblock_N.groupby("N")["C"]
        mu=g.mean();se=g.sem()
        B.errorbar(mu.index,mu.values,yerr=se.values,fmt="s",color="#d95f02",markersize=5,linewidth=0,capsize=2.5)
        s,b,r2,lo,hi=fit_loglog(mu.index.to_numpy(),mu.values)
        if np.isfinite(s):
            xs=np.logspace(np.log10(mu.index.min()),np.log10(mu.index.max()),60)
            B.plot(xs,10**(s*np.log10(xs)+b),"-",color="black",linewidth=1.0,
                   label=fr"slope={s:.2f} [{lo:.2f},{hi:.2f}], $R^2$={r2:.2f}")
        B.set_xscale("log");B.set_yscale("log")
        B.legend(loc="best",frameon=True,framealpha=0.9)
    B.set_xlabel(r"chain length $N$");B.set_ylabel("di-block contrast C")
    B.set_title("B. Di-block scaling",loc="left");panel_label(B,"B")
    if len(df_ratio_N):
        for pi,grp in df_ratio_N.groupby("pi"):
            g0=grp[grp["kappa"]==0.0].groupby("N")["C"].mean()
            g1=grp[grp["kappa"]==1.0].groupby("N")["C"].mean()
            common=g0.index.intersection(g1.index)
            ratio=g1.loc[common]/g0.loc[common]
            C.plot(common,ratio.values,"o-",markersize=4,linewidth=1.2,label=fr"$\pi={pi}$")
    C.set_xlabel(r"chain length $N$");C.set_ylabel(r"$C_{\kappa=1}/C_{\kappa=0}$")
    C.set_title(r"C. $\kappa$-amplification peaks at intermediate $N$",loc="left")
    C.legend(loc="best",frameon=True,framealpha=0.9);panel_label(C,"C")
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig6_aging_dynamics(dyn_by_kappa:dict,out:Path):
    """Fig 6. Two-time dynamical observables for the glass-aging claim.
    dyn_by_kappa: {kappa: {tw: {"fs": array_seeds_x_lags, "msd":..., "Q":..., "alpha2":...,
                                "lags_steps": array_lags}}}
    Panels: (A) F_s(k*, t; t_w) decay shape for kappa values at fixed t_w. (B) tau_alpha(t_w)
    vs t_w on log-log -- if it grows, that is aging. (C) chi_4(t; t_w) peak vs t_w. (D) MSD
    vs t showing ballistic-caged-diffusive regimes."""
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.6))
    A,B,C,D=axes[0,0],axes[0,1],axes[1,0],axes[1,1]
    kappas=sorted(dyn_by_kappa.keys())
    tw_ref=None
    for k in kappas:
        tws=sorted(dyn_by_kappa[k].keys())
        if not tws:continue
        tw_ref=tws[len(tws)//2]
        obs=dyn_by_kappa[k][tw_ref]
        fs=np.nanmean(obs["fs"],axis=0) if obs["fs"].ndim==2 else obs["fs"]
        lags=obs["lags_steps"]
        A.plot(lags,fs,color=kappa_color(k),linewidth=1.4,label=fr"$\kappa={k}$")
    A.axhline(1/np.e,color="grey",linestyle=":",linewidth=0.8,label=r"$1/e$")
    if tw_ref is not None:
        A.set_title(r"A. $F_s(k_\star,t;t_w\approx10^{%d})$"%int(round(np.log10(tw_ref))),loc="left")
    A.set_xscale("log");A.set_xlabel("BD lag steps");A.set_ylabel(r"$F_s$")
    A.set_ylim(-0.05,1.05);A.legend(loc="upper right",frameon=True,framealpha=0.9);panel_label(A,"A")
    for k in kappas:
        tws=sorted(dyn_by_kappa[k].keys())
        ta=[];tw_arr=[]
        for tw in tws:
            o=dyn_by_kappa[k][tw]
            if "tau_alpha" in o and np.isfinite(o["tau_alpha"]):
                ta.append(o["tau_alpha"]);tw_arr.append(tw)
        if len(tw_arr)>1:
            B.plot(tw_arr,ta,"o-",color=kappa_color(k),markersize=5,linewidth=1.4,label=fr"$\kappa={k}$")
            s,b,r2,lo,hi=fit_loglog(np.asarray(tw_arr),np.asarray(ta))
            if np.isfinite(s):
                xs=np.logspace(np.log10(min(tw_arr)),np.log10(max(tw_arr)),20)
                B.plot(xs,10**(s*np.log10(xs)+b),"--",color=kappa_color(k),alpha=0.5,linewidth=0.9)
    B.set_xscale("log");B.set_yscale("log");B.set_xlabel(r"waiting time $t_w$ (BD steps)")
    B.set_ylabel(r"$\tau_\alpha(t_w)$ (BD steps)")
    B.set_title(r"B. $\tau_\alpha$ growth = aging; flat = equilibrated",loc="left")
    B.legend(loc="best",frameon=True,framealpha=0.9);panel_label(B,"B")
    for k in kappas:
        tws=sorted(dyn_by_kappa[k].keys())
        peaks=[];tw_arr=[]
        for tw in tws:
            o=dyn_by_kappa[k][tw]
            ch=o.get("chi4")
            if ch is None or not np.isfinite(np.nanmax(ch)):continue
            peaks.append(float(np.nanmax(ch)));tw_arr.append(tw)
        if len(tw_arr)>1:
            C.plot(tw_arr,peaks,"o-",color=kappa_color(k),markersize=5,linewidth=1.4,label=fr"$\kappa={k}$")
    C.set_xscale("log");C.set_xlabel(r"$t_w$ (BD steps)");C.set_ylabel(r"$\chi_4^{\rm peak}(t_w)$")
    C.set_title("C. Cooperative-dynamics scale vs $t_w$",loc="left")
    C.legend(loc="best",frameon=True,framealpha=0.9);panel_label(C,"C")
    for k in kappas:
        tws=sorted(dyn_by_kappa[k].keys())
        if not tws:continue
        tw=tws[len(tws)//2]
        obs=dyn_by_kappa[k][tw]
        msd=np.nanmean(obs["msd"],axis=0) if obs["msd"].ndim==2 else obs["msd"]
        lags=obs["lags_steps"]
        m=lags>0
        D.plot(lags[m],msd[m],color=kappa_color(k),linewidth=1.4,label=fr"$\kappa={k}$")
    D.set_xscale("log");D.set_yscale("log");D.set_xlabel("BD lag steps")
    D.set_ylabel(r"MSD $\langle\Delta r^2\rangle\ [\sigma^2]$")
    D.set_title("D. MSD: ballistic / caged / diffusive",loc="left")
    D.legend(loc="best",frameon=True,framealpha=0.9);panel_label(D,"D")
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)

def fig7_atlas(df_pi_kappa:pd.DataFrame,df_fA_kappa:pd.DataFrame,df_eps_kappa:pd.DataFrame,out:Path):
    """Fig 7. Sequence-control atlas. Three (kappa x X) heat maps."""
    apply_style()
    fig,axes=plt.subplots(1,3,figsize=(9.5,3.2))
    for ax,df,xcol,xl,title,letter in zip(
            axes,
            [df_pi_kappa,df_fA_kappa,df_eps_kappa],
            ["pi","f_A","eps_AB"],
            [r"$\pi$",r"$f_A$",r"$\epsilon_{AB}/\epsilon_{AA}$"],
            ["A. $(\\kappa,\\pi)$ at $f_A=0.5$","B. $(\\kappa,f_A)$ at $\\pi=0.95$","C. $(\\kappa,\\epsilon_{AB})$ at $\\pi=0.99$"],
            ["A","B","C"]):
        if not len(df):continue
        d=df[df[xcol]<1.0] if xcol=="pi" else df
        pivot=d.groupby([xcol,"kappa"])["C"].mean().unstack("kappa")
        if xcol=="eps_AB":
            pivot=pivot.sort_index()
            im=ax.pcolormesh(pivot.columns.values,pivot.index.values,pivot.values,cmap="viridis",shading="auto")
            ax.set_yscale("log")
        else:
            im=ax.imshow(pivot.values,origin="lower",aspect="auto",
                         extent=[pivot.columns.min(),pivot.columns.max(),pivot.index.min(),pivot.index.max()],
                         cmap="viridis")
        cb=fig.colorbar(im,ax=ax,shrink=0.85,pad=0.02);cb.set_label("C",fontsize=7);cb.ax.tick_params(labelsize=6)
        ax.set_xlabel(r"$\kappa$");ax.set_ylabel(xl);ax.set_title(title,loc="left")
        panel_label(ax,letter)
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out);plt.close(fig)
