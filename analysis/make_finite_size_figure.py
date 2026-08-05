"""Finite-size figure: peak composition structure factor and its wavevector
versus system size at fixed bead density, for kappa = 0 and kappa = 1.

Reads the fixed-density campaign analysis and writes FigS9_finite_size.{pdf,png}
in the manuscript figures directory, matching the house style of make_figures.py.
"""
from __future__ import annotations
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parents[1]
FIG=HERE/'figures'
ROOT=Path(os.environ.get('IMP_V3_ROOT',HERE.parent))
CAMP=ROOT/'output/melt/fd_v2/fixed_density_pi099_v2'
ANA=CAMP/'analysis'
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    'font.size':9.5,'axes.labelsize':9.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,
    'legend.fontsize':9.5,'legend.markerscale':1.05,
    'figure.dpi':170,'savefig.dpi':520,
    'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans',
    'axes.linewidth':0.9,'xtick.major.width':0.8,'ytick.major.width':0.8,
    'pdf.fonttype':42,'ps.fonttype':42,
})

COL={'k0':'#3B6EA8','k1':'#7B3294'}

def panel(ax,label,size=12):
    ax.annotate(label,xy=(0,1),xycoords='axes fraction',xytext=(-12,7),
                textcoords='offset points',ha='right',va='bottom',
                fontsize=size,fontweight='bold',clip_on=False,annotation_clip=False)

def healthy_runs():
    """Runs whose final state is physically sane.

    Guards against runs that integrate unstably from an overlapping initial
    placement: those diverge to total energy ~1e27 and instantaneous T ~1e15 K
    against healthy values of ~-3e4 and 0.70. Six of the thirty runs failed
    this way before melt.box.relax_overlaps was introduced; all thirty pass
    with it. Kept as a standing check so a divergence is excluded rather than
    silently averaged in.
    """
    ok,bad=[],[]
    for d in sorted((CAMP/'runs').glob('*/')):
        s=pd.read_csv(d/'snapshots.csv')
        T=float(s.temperature_inst.iloc[-1]); E=float(s.total_energy.iloc[-1])
        Rg=float(s.mean_Rg.iloc[-1])
        (ok if (0.5<T<1.0 and E<0 and Rg<10) else bad).append(d.name)
    return set(ok),sorted(bad)

def main():
    ok,bad=healthy_runs()
    spec=pd.read_csv(ANA/'shell_spectra_by_seed.csv')
    spec=spec[spec.run_id.isin(ok)]
    rows=[]
    for (M,kap),g in spec.groupby(['n_chains','kappa']):
        m=g.groupby('q').S_psi_psi_over_2.agg(['mean','sem','count']).reset_index()
        m=m[m['count']==m['count'].max()]
        r=m.loc[m['mean'].idxmax()]
        rows.append(dict(n_chains=M,kappa=kap,n_seeds=int(r['count']),
                         k_star=r['q'],peak=r['mean'],peak_sem=r['sem'],
                         box_size=22.0*(M/144)**(1/3)))
    cond=pd.DataFrame(rows)
    seeds=spec.loc[spec.groupby(['run_id']).S_psi_psi_over_2.idxmax(),
                   ['run_id','n_chains','kappa','q','S_psi_psi_over_2']]

    fig,(axA,axB)=plt.subplots(1,2,figsize=(6.8,2.9),constrained_layout=True)
    for kap,col in [(0.0,COL['k0']),(1.0,COL['k1'])]:
        c=cond[np.isclose(cond.kappa,kap)].sort_values('n_chains')
        s=seeds[np.isclose(seeds.kappa,kap)]
        lab=fr'$\kappa={kap:g}$'
        # per-seed peaks, jittered, to show the real scatter behind the SEM
        jit=0.97 if np.isclose(kap,0.0) else 1.03
        axA.plot(s.n_chains*jit,s.S_psi_psi_over_2,'o',
                 ms=3.0,color=col,alpha=0.30,mew=0,zorder=1)
        axA.errorbar(c.n_chains,c.peak,yerr=c.peak_sem,
                     marker='o',ms=6,lw=1.8,capsize=3,color=col,label=lab,zorder=3)
        axB.plot(c.n_chains,c.k_star,marker='o',ms=6,lw=1.8,color=col,zorder=3)

    axA.set_yscale('log')
    axA.set_ylim(3,4e3)
    axA.set_ylabel(r'$S_{\psi\psi}(k_\star)/2$')
    axB.set_ylabel(r'$k_\star$ $(\sigma^{-1})$')
    axB.set_ylim(0,1.35)
    for ax in (axA,axB):
        ax.set_xscale('log')
        ax.set_xticks([144,288,576])
        ax.set_xticklabels(['144','288','576'])
        ax.minorticks_off()
        ax.set_xlabel('chains $M$ (fixed density)')
        ax.grid(alpha=0.18)
    # direct labels so identity is never colour-alone
    for kap,col,dy in [(1.0,COL['k1'],12),(0.0,COL['k0'],12)]:
        c=cond[np.isclose(cond.kappa,kap)].sort_values('n_chains').iloc[-1]
        axA.annotate(fr'$\kappa={kap:g}$',xy=(c.n_chains,c.peak),
                     xytext=(-6,dy),textcoords='offset points',
                     ha='right',color=col,fontsize=9.0)
    # seed count per point, since it is not five everywhere
    for _,r in cond.iterrows():
        axA.annotate(f'n={r.n_seeds:.0f}',xy=(r.n_chains,r.peak),xytext=(7,-9),
                     textcoords='offset points',fontsize=6.8,color='0.45')
    handles,labels=axA.get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,
               bbox_to_anchor=(0.50,1.00),fontsize=9.5)
    for ax,lab in zip((axA,axB),'AB'):
        panel(ax,lab)
    for ext in ('pdf','png'):
        fig.savefig(FIG/f'FigS9_finite_size.{ext}',bbox_inches='tight',
                    pad_inches=0.025,transparent=True)
    plt.close(fig)

    if bad:
        print(f'excluded {len(bad)} diverged runs:')
        for b in bad: print('  ',b)
    print('conditions (healthy seeds only)')
    for _,r in cond.sort_values(['kappa','n_chains']).iterrows():
        print(f"  M={int(r.n_chains):4d} kappa={r.kappa:g} L={r.box_size:6.2f} "
              f"n={r.n_seeds:.0f}  S={r.peak:8.2f} +/- {r.peak_sem:6.2f}"
              f"  k*={r.k_star:.4f}")
    cond.to_csv(ANA/'condition_summary_healthy.csv',index=False)
    print(f"wrote {FIG/'FigS9_finite_size.pdf'}")

if __name__=='__main__':
    main()
