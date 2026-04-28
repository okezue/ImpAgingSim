from __future__ import annotations
import argparse,glob,json,os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_chunk_runs(scan_dir):
    rows=[]
    for d in sorted(glob.glob(scan_dir+"/*/")):
        sf_p=os.path.join(d,"structure_factor.npz")
        csv_p=os.path.join(d,"snapshots.csv")
        if not(os.path.exists(sf_p)and os.path.exists(csv_p)):continue
        rn=os.path.basename(d.rstrip('/'))
        try:
            kappa=float(rn.split('_T')[0].lstrip('k'))
            T_q=float(rn.split('_T')[1].split('_')[0])
            seed=int(rn.split('_s')[-1])
        except:continue
        sf=np.load(sf_p)
        k=sf['k'];SAA=sf['S_AA'][-5:].mean(axis=0);SAB=sf['S_AB'][-5:].mean(axis=0)
        contrast=SAA-SAB
        valid=(k>0.5)&np.isfinite(contrast)
        if not np.any(valid):continue
        idx=np.where(valid)[0]
        j_c=idx[np.argmax(contrast[idx])]
        j_s=idx[np.argmax(SAA[idx])]
        df=pd.read_csv(csv_p)
        rows.append({
            'kappa':kappa,'T_q':T_q,'seed':seed,
            'k_star_contrast':k[j_c],'contrast_peak':contrast[j_c],
            'SAA_peak':SAA[j_s],'SAA_at_kc':SAA[j_c],'SAB_at_kc':SAB[j_c],
            'E_final':df['total_energy'].iloc[-1] if len(df) else float('nan'),
            'Rg_final':df['mean_Rg'].iloc[-1] if len(df) else float('nan'),
        })
    return pd.DataFrame(rows)

def subblock_rg(traj_path,meta_path):
    td=np.load(traj_path)
    pos=td['positions'];types=td['types'];L=float(td['box_size'])
    with open(meta_path) as f:meta=json.load(f)
    n_chains=meta['melt_params']['n_chains']
    chain_length=meta['melt_params']['chain_length']
    seq_arr=np.array(meta['sequence'])
    last=pos[-1]
    Rg_full=[];Rg_A=[];Rg_B=[]
    for c in range(n_chains):
        s,e=c*chain_length,(c+1)*chain_length
        cp=last[s:e];ct=seq_arr[s:e]
        cm=cp.mean(axis=0)
        Rg_full.append(np.sqrt(np.mean(np.sum((cp-cm)**2,axis=1))))
        for mask,store in [(ct==1,Rg_A),(ct==0,Rg_B)]:
            if mask.sum()>1:
                cmm=cp[mask].mean(axis=0)
                store.append(np.sqrt(np.mean(np.sum((cp[mask]-cmm)**2,axis=1))))
    return float(np.mean(Rg_full)),float(np.mean(Rg_A)) if Rg_A else float('nan'),float(np.mean(Rg_B)) if Rg_B else float('nan')

def intra_globule_variance(traj_path,meta_path,grid_size=16,min_beads_per_cell=5):
    td=np.load(traj_path)
    pos=td['positions'];types=td['types'];L=float(td['box_size'])
    with open(meta_path) as f:meta=json.load(f)
    seq_arr=np.array(meta['sequence'])
    last=pos[-1]
    G=grid_size;cs=L/G
    bc=((last%L)//cs).astype(int)%G
    ci=bc[:,0]*G*G+bc[:,1]*G+bc[:,2]
    cA=np.zeros(G**3);cT=np.zeros(G**3)
    for i in range(len(ci)):
        cT[ci[i]]+=1
        if seq_arr[i]==1:cA[ci[i]]+=1
    occ=cT>=min_beads_per_cell
    if occ.sum()<3:
        return float('nan'),float('nan'),0
    Af=cA[occ]/cT[occ]
    var=float(np.var(Af))
    n_avg=float(cT[occ].mean())
    p=float(seq_arr.mean())
    rand_var=p*(1-p)/n_avg if n_avg>0 else float('nan')
    enh=var/rand_var if rand_var>0 else float('nan')
    return var,enh,int(occ.sum())

def plot_kt_heatmap(df,outpath,value_col='SAA_peak',cmap='viridis',title=None):
    pv=df.pivot_table(index='T_q',columns='kappa',values=value_col,aggfunc='mean')
    fig,ax=plt.subplots(figsize=(8,6))
    im=ax.imshow(pv.values,aspect='auto',origin='lower',cmap=cmap,
                 extent=[pv.columns.min(),pv.columns.max(),pv.index.min(),pv.index.max()])
    ax.set_xlabel('κ (correlation)');ax.set_ylabel('T_quench')
    ax.set_title(title or f'{value_col} across (κ, T)')
    plt.colorbar(im,ax=ax)
    plt.tight_layout();plt.savefig(outpath,dpi=140);plt.close()
    return outpath

def plot_kappa_curve(df,outpath,T_target=0.7,T_tol=0.05,value_col='SAA_peak',title=None):
    sel=df[(df['T_q']>T_target-T_tol)&(df['T_q']<T_target+T_tol)]
    agg=sel.groupby('kappa')[value_col].agg(['mean','std','count'])
    agg['sem']=agg['std']/np.sqrt(agg['count'])
    fig,ax=plt.subplots(figsize=(7,5))
    ax.errorbar(agg.index,agg['mean'],yerr=agg['sem'],marker='o',lw=2,capsize=3)
    ax.set_xlabel('κ');ax.set_ylabel(value_col)
    ax.set_title(title or f'{value_col} vs κ at T={T_target}')
    ax.grid(True,alpha=0.3)
    plt.tight_layout();plt.savefig(outpath,dpi=140);plt.close()
    return outpath

def main():
    p=argparse.ArgumentParser()
    p.add_argument("scan_dir",type=str)
    p.add_argument("--outdir",type=str,default=None)
    p.add_argument("--T_target",type=float,default=0.7)
    args=p.parse_args()
    od=args.outdir or os.path.dirname(args.scan_dir.rstrip('/'))
    os.makedirs(od,exist_ok=True)
    df=load_chunk_runs(args.scan_dir)
    if len(df)==0:
        print(f"No runs found in {args.scan_dir}");return
    df.to_csv(os.path.join(od,"runs.csv"),index=False)
    plot_kt_heatmap(df,os.path.join(od,"Fig_SAA_kT.png"),'SAA_peak','viridis','raw S_AA peak')
    plot_kt_heatmap(df,os.path.join(od,"Fig_contrast_kT.png"),'contrast_peak','magma','contrast S_AA-S_AB')
    plot_kt_heatmap(df,os.path.join(od,"Fig_E_kT.png"),'E_final','RdBu_r','final energy')
    plot_kappa_curve(df,os.path.join(od,"Fig_SAA_vs_kappa.png"),args.T_target,value_col='SAA_peak',title=f'raw S_AA peak vs κ at T={args.T_target}')
    plot_kappa_curve(df,os.path.join(od,"Fig_contrast_vs_kappa.png"),args.T_target,value_col='contrast_peak',title=f'contrast peak vs κ at T={args.T_target}')
    print(f"Wrote outputs to {od}")
    print(df.groupby('kappa')[['SAA_peak','contrast_peak','E_final']].mean().round(4).to_string())

if __name__=="__main__":
    main()
