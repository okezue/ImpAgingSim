from __future__ import annotations
import numpy as np

def generate_random(N,f_A,rng):
    return (rng.random(N)<f_A).astype(np.int8)

def generate_block(N,block_length,f_A=0.5):
    s=np.empty(N,dtype=np.int8)
    nA=int(round(block_length*f_A*2))
    nA=max(1,nA)
    nB=max(1,2*block_length-nA)
    period=nA+nB
    for i in range(N):
        s[i]=1 if (i%period)<nA else 0
    return s

def generate_alternating(N):
    s=np.zeros(N,dtype=np.int8)
    s[::2]=1
    return s

def valid_pi_range(f_A):
    """Return (pi_min, pi_max=1) for which P(A->B)=2(1-pi)(1-f_A) and P(B->A)=2(1-pi)f_A
    both stay in [0,1]. pi_min = max(0, 1 - 1/(2 max(f_A,1-f_A)))."""
    fm=max(float(f_A),1.0-float(f_A))
    return max(0.0,1.0-1.0/(2.0*fm)),1.0

def generate_correlated(N,kappa,pi,f_A,rng):
    """Unified Markov correlated generator with eigenvalue lambda=2*pi-1 for all f_A.

    Transition probabilities: P(A->B)=2(1-pi)(1-f_A), P(B->A)=2(1-pi)f_A. This matches
    the manuscript formula (Eq. eq:markov) at all compositions; the symmetric f_A=0.5
    case reduces to P(stay)=pi. Detailed balance gives stationary P(A)=f_A. The second
    eigenvalue of the transition matrix is 2*pi-1 for every f_A.

    Replaces both the old symmetric generate_correlated (which was {-1,+1} only) and
    generate_correlated_biased (which used a=(1-pi) and gave lambda=pi for f_A!=0.5).

    Sequence is mixed with i.i.d. Bernoulli(f_A) draws: with probability kappa each bead
    inherits the Markov realization, with probability 1-kappa it is resampled i.i.d.

    Raises ValueError if pi is outside the valid range for this f_A; use valid_pi_range()
    to query the allowed window.
    """
    if not(0.0<=pi<=1.0):
        raise ValueError(f"pi in [0,1], got {pi}")
    if not(0.0<=kappa<=1.0):
        raise ValueError(f"kappa in [0,1], got {kappa}")
    if not(0.0<f_A<1.0):
        raise ValueError(f"f_A in (0,1), got {f_A}")
    p_AB=2.0*(1.0-pi)*(1.0-f_A)
    p_BA=2.0*(1.0-pi)*f_A
    if p_AB>1.0 or p_BA>1.0:
        pmin,_=valid_pi_range(f_A)
        raise ValueError(f"pi={pi} below valid_pi_min={pmin:.4f} for f_A={f_A}")
    z=np.empty(N,dtype=np.int8)
    z[0]=1 if rng.random()<f_A else 0
    for i in range(N-1):
        if z[i]==1:
            z[i+1]=0 if rng.random()<p_AB else 1
        else:
            z[i+1]=1 if rng.random()<p_BA else 0
    if kappa<1.0:
        flip=rng.random(N)>kappa
        rand=(rng.random(N)<f_A).astype(np.int8)
        z=np.where(flip,rand,z).astype(np.int8)
    return z

def generate_per_chain(kind,n_chains,chain_length,f_A,block_length,kappa,pi,rng):
    """Generate one sequence per chain independently, then concatenate.
    Eliminates the across-chain Markov bleed-through of the previous run.py code path."""
    out=np.empty(n_chains*chain_length,dtype=np.int8)
    for c in range(n_chains):
        if kind=="random":
            s=generate_random(chain_length,f_A,rng)
        elif kind=="block":
            s=generate_block(chain_length,block_length,f_A)
        elif kind=="alternating":
            s=generate_alternating(chain_length)
        elif kind=="correlated":
            s=generate_correlated(chain_length,kappa,pi,f_A,rng)
        else:
            raise ValueError(f"unknown sequence kind: {kind}")
        out[c*chain_length:(c+1)*chain_length]=s
    return out

def generate_per_chain_exact_total(kind,n_chains,chain_length,f_A,block_length,kappa,pi,rng,
                                   max_attempts=100000):
    """Sample independent per-chain sequences conditional on an exact global A count.

    Each attempt draws the full set with :func:`generate_per_chain`.  The first set whose
    total A count equals ``n_chains*chain_length*f_A`` is retained.  Thus the accepted law is
    the original independent-chain disorder ensemble conditioned only on canonical global
    composition; no chain is copied, complemented, or otherwise constructed from another.
    """
    M=int(n_chains);N=int(chain_length);attempt_limit=int(max_attempts)
    if M<1 or N<1:
        raise ValueError("n_chains and chain_length must be positive")
    target_float=M*N*float(f_A)
    target=int(round(target_float))
    if not np.isclose(target_float,target,rtol=0.0,atol=1e-12):
        raise ValueError("exact global composition requires an integer target A count")
    if attempt_limit<1:
        raise ValueError("max_attempts must be positive")
    for attempt in range(1,attempt_limit+1):
        sequence=generate_per_chain(kind,M,N,f_A,block_length,kappa,pi,rng)
        if int(np.sum(sequence,dtype=np.int64))==target:
            return sequence,attempt
    raise RuntimeError(
        f"failed to sample exact A count {target}/{M*N} after {attempt_limit} independent attempts"
    )

def autocorrelation(seq,kmax,per_chain_length=None):
    """Centred normalized autocorrelation. If per_chain_length is given, compute
    within-chain autocorrelation (averaged over chains) instead of treating the
    concatenated sequence as one long chain."""
    s=2*seq.astype(np.float64)-1.0
    if per_chain_length is None:
        N=len(s)
        out=np.empty(kmax+1)
        for k in range(kmax+1):
            out[k]=float(np.mean(s[:N-k]*s[k:])) if k<N else 0.0
        return out
    L=int(per_chain_length)
    n_chains=len(s)//L
    out=np.zeros(kmax+1)
    for c in range(n_chains):
        sc=s[c*L:(c+1)*L]
        for k in range(kmax+1):
            out[k]+=float(np.mean(sc[:L-k]*sc[k:])) if k<L else 0.0
    return out/max(n_chains,1)
