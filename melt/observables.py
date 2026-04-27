from __future__ import annotations
import numpy as np

def fft_density(phi):
    return np.fft.fftn(phi)

def structure_factor_3d(phi,box_size,subtract_mean=True):
    G=phi.shape[0];V=float(box_size)**3
    p=phi-phi.mean() if subtract_mean else phi
    f=fft_density(p)
    return (np.abs(f)**2)*(V/(G**6))

def cross_structure_factor_3d(phi1,phi2,box_size,subtract_mean=True):
    G=phi1.shape[0];V=float(box_size)**3
    p1=phi1-phi1.mean() if subtract_mean else phi1
    p2=phi2-phi2.mean() if subtract_mean else phi2
    f1=fft_density(p1);f2=fft_density(p2)
    return np.real(f1*np.conj(f2))*(V/(G**6))

def k_grid(grid_size,box_size):
    G=int(grid_size);L=float(box_size)
    kx=2.0*np.pi*np.fft.fftfreq(G,d=L/G)
    KX,KY,KZ=np.meshgrid(kx,kx,kx,indexing="ij")
    return np.sqrt(KX*KX+KY*KY+KZ*KZ)

def spherical_average(S3d,box_size,n_bins=None,k_max=None):
    G=S3d.shape[0];L=float(box_size)
    K=k_grid(G,L)
    kmin=2.0*np.pi/L
    kmax_def=np.pi*G/L
    if k_max is None:
        k_max=kmax_def
    if n_bins is None:
        n_bins=G//2
    edges=np.linspace(0.0,k_max,n_bins+1)
    centers=0.5*(edges[1:]+edges[:-1])
    Kf=K.ravel();Sf=S3d.ravel()
    idx=np.digitize(Kf,edges)-1
    out=np.zeros(n_bins);cnt=np.zeros(n_bins)
    valid=(idx>=0)&(idx<n_bins)
    np.add.at(out,idx[valid],Sf[valid])
    np.add.at(cnt,idx[valid],1.0)
    cnt=np.where(cnt>0,cnt,1.0)
    return centers,out/cnt

def find_peak(k,S,k_min=None):
    mask=np.isfinite(S)
    if k_min is not None:
        mask=mask&(k>=k_min)
    if not np.any(mask):
        return float("nan"),float("nan")
    ki=np.where(mask)[0]
    j=int(ki[np.argmax(S[ki])])
    return float(k[j]),float(S[j])

def domain_length(k_star):
    if not np.isfinite(k_star) or k_star<=0:
        return float("nan")
    return float(2.0*np.pi/k_star)

def fit_ornstein_zernike(k,S,k_max=None):
    mask=np.isfinite(S)&(S>0)&(k>0)
    if k_max is not None:
        mask=mask&(k<=k_max)
    if mask.sum()<3:
        return float("nan"),float("nan")
    x=k[mask]**2
    y=1.0/S[mask]
    a,b=np.polyfit(x,y,1)
    if a<=0 or b<=0:
        return float("nan"),float("nan")
    S0=1.0/b
    xi=float(np.sqrt(a/b))
    return xi,float(S0)

def compute_all_observables(phi_A,phi_B,box_size,n_bins=None,k_min_factor=1.5):
    L=float(box_size);G=phi_A.shape[0]
    SA3=structure_factor_3d(phi_A,L)
    SB3=structure_factor_3d(phi_B,L)
    SAB3=cross_structure_factor_3d(phi_A,phi_B,L)
    k,SA=spherical_average(SA3,L,n_bins=n_bins)
    _,SB=spherical_average(SB3,L,n_bins=n_bins)
    _,SAB=spherical_average(SAB3,L,n_bins=n_bins)
    k_min=k_min_factor*(2.0*np.pi/L)
    kstar,Speak=find_peak(k,SA,k_min=k_min)
    xi=domain_length(kstar)
    return {"k":k,"S_AA":SA,"S_BB":SB,"S_AB":SAB,
            "k_star":kstar,"S_AA_peak":Speak,"xi_AA":xi}
