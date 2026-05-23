from __future__ import annotations
import argparse,subprocess,sys,json
from pathlib import Path
import numpy as np
import pandas as pd
from .ingest import walk_chunk,chunk_dataframe
from .observables import (rpa_predictor_finite,aging_slope_bootstrap,dynamics_summary,
                          tau_alpha_per_tw,chi4_per_tw)
from .figures import (fig1_kappa_robustness,fig2_blockiness_corrected,fig3_composition,
                      fig4_interactions_density,fig5_architecture,fig6_aging_dynamics,fig7_atlas)

ROOT=Path(__file__).resolve().parents[2]
DEFAULT_S3=f"s3://okezue-imp-aging-results"
DEFAULT_LOCAL=ROOT/"output"/"aws_corrected"
FIGDIR=ROOT/"paper"/"figures_corrected"
CSVDIR=ROOT/"paper"/"source_data_corrected"

CHUNK_MAP={
    "fig1_baseline":"fig1_baseline",
    "fig1_dense":"fig1_dense",
    "fig1_soft":"fig1_soft",
    "fig1_short_chain":"fig1_short_chain",
    "fig2_pi_kappa":"fig2_pi_kappa",
    "fig3_fA_kappa":"fig3_fA_kappa",
    "fig4_epsAB":"fig4_epsAB",
    "fig5_aging":"fig5_aging",
}

def pull_s3(local:Path,bucket:str,prefix:str=""):
    local.mkdir(parents=True,exist_ok=True)
    src=f"s3://{bucket}/{prefix}".rstrip("/")+"/"
    print(f"syncing {src} -> {local}")
    subprocess.run(["aws","s3","sync",src,str(local),"--region","us-east-1"],check=True)

def load_chunk(local:Path,name:str):
    cd=local/name
    runs=walk_chunk(cd)
    if not runs:
        for sub in local.glob(f"*/scans_corrected/{name}"):
            runs=walk_chunk(sub);break
    if not runs:
        for sub in local.glob(f"*/{name}"):
            runs=walk_chunk(sub);break
    return chunk_dataframe(runs)

def save_csv(df:pd.DataFrame,name:str):
    CSVDIR.mkdir(parents=True,exist_ok=True)
    p=CSVDIR/f"{name}.csv";df.to_csv(p,index=False);print(f"  -> {p} ({len(df)} rows)")

def cmd_observables(args):
    local=Path(args.local)
    out={}
    for label,chunk in CHUNK_MAP.items():
        df=load_chunk(local,chunk)
        out[label]=df
        if len(df):save_csv(df,label)
        else:print(f"  warning: no runs found for {label}")
    return out

def cmd_dynamics(args):
    local=Path(args.local)
    aging_dir=local/"fig5_aging"
    if not aging_dir.exists():
        for c in local.glob(f"*/fig5_aging"):aging_dir=c;break
    runs=walk_chunk(aging_dir)
    if not runs:print("no aging runs found");return {}
    by_kappa={}
    Q_by_kappa_tw={}
    for r in runs:
        m=r["meta"];rd=r["dir"]
        name=rd.name
        try:
            parts=dict(p.split("=") if "=" in p else (p[:1],p[1:]) for p in name.split("_"))
            kappa=float(name.split("k")[1].split("_")[0]) if name.startswith("k") else float(m.get("kappa",np.nan))
            tw=int(m.get("run_params",{}).get("equilibration_steps",0))
        except Exception:
            kappa=np.nan;tw=0
        tjp=rd/"trajectory.npz"
        if not tjp.exists():continue
        L=float(m.get("melt_params",{}).get("box_size",0))
        sf=r["sf"];kvec=sf.get("k");SAA=sf.get("S_AA");SAB=sf.get("S_AB")
        if kvec is None or SAA is None:continue
        SAAm=SAA[-5:].mean(0);SABm=SAB[-5:].mean(0) if SAB is not None else np.zeros_like(SAAm)
        i=int(np.argmax((SAAm-SABm)[kvec>1.5*(2*np.pi/L)]))
        k_star=float(kvec[kvec>1.5*(2*np.pi/L)][i])
        snap=r["snap"]
        steps_per_lag=int(snap["step"].diff().dropna().iloc[0]) if snap is not None and len(snap)>1 else 1
        n_snap=len(sf["steps"])
        tw_indices=[0,n_snap//4,n_snap//2,3*n_snap//4][:max(1,n_snap-2)]
        dyn=dynamics_summary(tjp,L,k_star,tw_indices)
        for tw_idx,obs in dyn.items():
            obs["lags_steps"]=obs["lags"]*steps_per_lag
        ta=tau_alpha_per_tw(dyn,steps_per_lag)
        by_kappa.setdefault(kappa,{}).setdefault(tw,{})["msd"]=np.atleast_2d(dyn[tw_indices[0]]["msd"])
        by_kappa[kappa][tw]["fs"]=np.atleast_2d(dyn[tw_indices[0]]["fs"])
        by_kappa[kappa][tw]["alpha2"]=np.atleast_2d(dyn[tw_indices[0]]["alpha2"])
        by_kappa[kappa][tw]["Q"]=np.atleast_2d(dyn[tw_indices[0]]["Q_self"])
        by_kappa[kappa][tw]["lags_steps"]=dyn[tw_indices[0]]["lags_steps"]
        by_kappa[kappa][tw]["tau_alpha"]=ta[tw_indices[0]]
        Q_by_kappa_tw.setdefault(kappa,{}).setdefault(tw,[]).append(dyn[tw_indices[0]]["Q_self"])
    for k,d in by_kappa.items():
        for tw,o in d.items():
            qs=np.array(Q_by_kappa_tw[k][tw])
            o["chi4"]=qs.shape[0]*np.var(qs,axis=0,ddof=1) if qs.shape[0]>1 else np.full(qs.shape[-1],np.nan)
    np.savez_compressed(CSVDIR/"fig6_dynamics_cache.npz",
                        cache=json.dumps({str(k):{str(tw):{kk:vv.tolist() if isinstance(vv,np.ndarray) else float(vv)
                                                          for kk,vv in obs.items()} for tw,obs in d.items()}
                                          for k,d in by_kappa.items()}))
    return by_kappa

def cmd_figures(args,obs_dfs=None,dyn=None):
    if obs_dfs is None:
        obs_dfs={}
        for label in CHUNK_MAP:
            p=CSVDIR/f"{label}.csv"
            obs_dfs[label]=pd.read_csv(p) if p.exists() else pd.DataFrame()
    FIGDIR.mkdir(parents=True,exist_ok=True)
    print("rendering Fig 1...")
    df_by_cond={c:obs_dfs.get(f"fig1_{c}",pd.DataFrame()) for c in ["baseline","dense","soft","short_chain"]}
    fig1_kappa_robustness(df_by_cond,FIGDIR/"Fig1_kappa_robustness.pdf")
    print("rendering Fig 2...")
    pk=obs_dfs.get("fig2_pi_kappa",pd.DataFrame())
    df_xi=pk[(pk["pi"]>0.98)&(pk["pi"]<1.0)] if len(pk) else pk
    fig2_blockiness_corrected(pk,df_xi,FIGDIR/"Fig2_blockiness_corrected.pdf")
    print("rendering Fig 3...")
    fa=obs_dfs.get("fig3_fA_kappa",pd.DataFrame())
    df_at07=fa[(fa["kappa"]==0.7)] if len(fa) else fa
    df_hp=fa[(fa["pi"]>=0.99)] if len(fa) and "pi" in fa else fa
    fig3_composition(fa,df_at07,df_hp,FIGDIR/"Fig3_composition.pdf")
    print("rendering Fig 4...")
    eps=obs_dfs.get("fig4_epsAB",pd.DataFrame())
    df_em=eps[(eps["kappa"]==0.7)&(eps["pi"]<0.99)] if len(eps) else eps
    df_eh=eps[(eps["kappa"]==1.0)&(eps["pi"]>=0.99)] if len(eps) else eps
    df_den=eps
    slice_png=ROOT/"paper"/"figures"/"Fig4_density_slice.png"
    fig4_interactions_density(df_em,df_eh,df_den,slice_png,FIGDIR/"Fig4_interactions.pdf")
    print("rendering Fig 5...")
    fig5_architecture(pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),FIGDIR/"Fig5_architecture.pdf")
    print("rendering Fig 6 (dynamics)...")
    if dyn is None:
        cache=CSVDIR/"fig6_dynamics_cache.npz"
        if cache.exists():
            raw=json.loads(np.load(cache,allow_pickle=False)["cache"].item())
            dyn={float(k):{int(tw):{kk:(np.array(vv) if isinstance(vv,list) else vv) for kk,vv in obs.items()}
                          for tw,obs in d.items()} for k,d in raw.items()}
        else:dyn={}
    if dyn:fig6_aging_dynamics(dyn,FIGDIR/"Fig6_aging_dynamics.pdf")
    else:print("  no dynamics cache; run --dynamics first")
    print("rendering Fig 7 atlas...")
    fig7_atlas(pk,fa,eps,FIGDIR/"Fig7_atlas.pdf")
    print(f"figures -> {FIGDIR}")

def main():
    p=argparse.ArgumentParser(description="Regenerate corrected manuscript figures.")
    p.add_argument("--bucket",default="okezue-imp-aging-results")
    p.add_argument("--prefix",default="")
    p.add_argument("--local",default=str(DEFAULT_LOCAL))
    p.add_argument("--pull-s3",action="store_true")
    p.add_argument("--observables",action="store_true")
    p.add_argument("--dynamics",action="store_true")
    p.add_argument("--figures",action="store_true")
    p.add_argument("--all",action="store_true")
    a=p.parse_args()
    if a.all:a.pull_s3=a.observables=a.dynamics=a.figures=True
    if a.pull_s3:pull_s3(Path(a.local),a.bucket,a.prefix)
    obs=None;dyn=None
    if a.observables:obs=cmd_observables(a)
    if a.dynamics:dyn=cmd_dynamics(a)
    if a.figures:cmd_figures(a,obs,dyn)
    if not(a.pull_s3 or a.observables or a.dynamics or a.figures):
        p.print_help()

if __name__=="__main__":
    main()
