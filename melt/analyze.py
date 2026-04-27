from __future__ import annotations
import argparse,glob,json,os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEQ_COLORS={"random":"#1f77b4","block":"#d62728","alternating":"#2ca02c","correlated":"#9467bd"}
SEQ_LABEL={"random":"random","block":"block","alternating":"alternating","correlated":"correlated"}

def load_scan(scan_dir):
    rows=[]
    sf_data={}
    for run_dir in sorted(glob.glob(os.path.join(scan_dir,"*"))):
        if not os.path.isdir(run_dir):continue
        csv_p=os.path.join(run_dir,"snapshots.csv")
        meta_p=os.path.join(run_dir,"meta.json")
        sf_p=os.path.join(run_dir,"structure_factor.npz")
        if not (os.path.exists(csv_p) and os.path.exists(meta_p)):continue
        with open(meta_p) as f:meta=json.load(f)
        df=pd.read_csv(csv_p)
        df["seed"]=meta["run_params"]["seed"]
        df["sequence"]=meta["sequence_type"]
        rows.append(df)
        if os.path.exists(sf_p):
            sf=np.load(sf_p)
            sf_data[(meta["sequence_type"],meta["run_params"]["seed"])]={
                "steps":sf["steps"],"k":sf["k"],"S_AA":sf["S_AA"],
                "S_BB":sf["S_BB"],"S_AB":sf["S_AB"]}
    if not rows:
        raise RuntimeError(f"No runs found in {scan_dir}")
    return pd.concat(rows,ignore_index=True),sf_data

def aggregate_by_seq(df,col):
    g=df.groupby(["sequence","step"])[col]
    mean=g.mean().reset_index().rename(columns={col:f"{col}_mean"})
    std=g.std().reset_index().rename(columns={col:f"{col}_std"})
    n=g.count().reset_index().rename(columns={col:f"{col}_n"})
    out=mean.merge(std,on=["sequence","step"]).merge(n,on=["sequence","step"])
    out[f"{col}_sem"]=out[f"{col}_std"]/np.sqrt(np.maximum(out[f"{col}_n"],1))
    return out

def plot_xi_t(df,outpath):
    agg=aggregate_by_seq(df,"S_AA_peak")
    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    for seq,sub in agg.groupby("sequence"):
        c=SEQ_COLORS.get(seq,"k")
        sub=sub.sort_values("step")
        ax[0].plot(sub["step"],sub["S_AA_peak_mean"],color=c,label=SEQ_LABEL.get(seq,seq),lw=2)
        ax[0].fill_between(sub["step"],
                           sub["S_AA_peak_mean"]-sub["S_AA_peak_sem"],
                           sub["S_AA_peak_mean"]+sub["S_AA_peak_sem"],color=c,alpha=0.2)
    ax[0].set_xlabel("steps after quench");ax[0].set_ylabel(r"$S_{AA}(k_\star)$")
    ax[0].set_title("peak segregation amplitude vs time")
    ax[0].legend(frameon=False)
    agg2=aggregate_by_seq(df,"xi_AA")
    for seq,sub in agg2.groupby("sequence"):
        c=SEQ_COLORS.get(seq,"k")
        sub=sub.sort_values("step")
        ax[1].plot(sub["step"],sub["xi_AA_mean"],color=c,label=SEQ_LABEL.get(seq,seq),lw=2)
        ax[1].fill_between(sub["step"],
                           sub["xi_AA_mean"]-sub["xi_AA_sem"],
                           sub["xi_AA_mean"]+sub["xi_AA_sem"],color=c,alpha=0.2)
    ax[1].set_xlabel("steps after quench");ax[1].set_ylabel(r"$\xi_{AA}=2\pi/k_\star$")
    ax[1].set_title("domain length vs time")
    ax[1].legend(frameon=False)
    plt.tight_layout();plt.savefig(outpath,dpi=140);plt.close()

def plot_final_Sk(sf_data,outpath,n_avg=3):
    seqs=sorted(set(s for s,_ in sf_data.keys()))
    fig,ax=plt.subplots(figsize=(7,5))
    for seq in seqs:
        keys=[(s,sd) for (s,sd) in sf_data.keys() if s==seq]
        if not keys:continue
        last_avgs=[]
        k_ref=None
        for s,sd in keys:
            d=sf_data[(s,sd)]
            if k_ref is None:k_ref=d["k"]
            last=d["S_AA"][-n_avg:].mean(axis=0)
            last_avgs.append(last)
        Sm=np.mean(np.stack(last_avgs,axis=0),axis=0)
        Ss=np.std(np.stack(last_avgs,axis=0),axis=0)/np.sqrt(len(last_avgs))
        c=SEQ_COLORS.get(seq,"k")
        ax.plot(k_ref,Sm,color=c,label=SEQ_LABEL.get(seq,seq),lw=2)
        ax.fill_between(k_ref,Sm-Ss,Sm+Ss,color=c,alpha=0.2)
    ax.set_xscale("log");ax.set_yscale("log")
    ax.set_xlabel("k");ax.set_ylabel(r"$S_{AA}(k)$")
    ax.set_title(f"final-state structure factor (avg of last {n_avg} snapshots)")
    ax.legend(frameon=False)
    plt.tight_layout();plt.savefig(outpath,dpi=140);plt.close()

def summary_table(df):
    g=df.groupby("sequence").agg(
        S_AA_peak_init=("S_AA_peak","first"),
        S_AA_peak_final=("S_AA_peak","last"),
        xi_init=("xi_AA","first"),
        xi_final=("xi_AA","last"),
        E_init=("total_energy","first"),
        E_final=("total_energy","last"),
        Rg_mean=("mean_Rg","mean"))
    g["S_AA_growth_x"]=g["S_AA_peak_final"]/g["S_AA_peak_init"]
    return g

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("scan_dir",type=str)
    ap.add_argument("--outdir",type=str,default=None)
    args=ap.parse_args()
    od=args.outdir or args.scan_dir
    df,sf=load_scan(args.scan_dir)
    plot_xi_t(df,os.path.join(od,"Fig_xi_vs_t.png"))
    if sf:
        plot_final_Sk(sf,os.path.join(od,"Fig_Sk_final.png"))
    tab=summary_table(df)
    tab_p=os.path.join(od,"summary.csv")
    tab.to_csv(tab_p)
    print(tab.to_string())
    print(f"\nWrote: Fig_xi_vs_t.png, Fig_Sk_final.png, summary.csv -> {od}")

if __name__=="__main__":
    main()
