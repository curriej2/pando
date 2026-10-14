"""Probe for Betancourt §2.4 (notes §G.4.15): cost per independent draw, random-walk Metropolis vs a simple
static HMC, on standard Gaussians of dimension D (part 1) and on an unpreconditioned 2-D Gaussian with SD
ratio k (part 2). ESS of q_1 by batch means. Teaching probe, run 2026-10-13 from the login node (seconds);
printed numbers are the ones recorded in §G.4.15. Seeds fixed below.
"""
import numpy as np, time
def ess_bm(x, nb=50):                       # batch-means ESS
    n=len(x)//nb*nb; b=x[:n].reshape(nb,-1).mean(1); m=n//nb
    return x[:n].var()/ (m*b.var(ddof=1)) * n
def rwm(D, n, rng):
    s=2.38/np.sqrt(D); q=rng.standard_normal(D); lp=-0.5*q@q; out=np.empty(n); acc=0
    for t in range(n):
        qp=q+s*rng.standard_normal(D); lpp=-0.5*qp@qp
        if np.log(rng.random())<lpp-lp: q,lp=qp,lpp; acc+=1
        out[t]=q[0]
    return out, acc/n
def hmc(D, n, rng, c=1.0):
    eps=c*D**-0.25; q=rng.standard_normal(D); out=np.empty(n); acc=0; grads=0
    for t in range(n):
        p=rng.standard_normal(D); T=rng.uniform(0.8,1.6); L=max(1,int(np.ceil(T/eps)))
        q0,p0=q.copy(),p.copy(); H0=0.5*q@q+0.5*p@p
        p=p-0.5*eps*q
        for l in range(L):
            q=q+eps*p
            if l<L-1: p=p-eps*q
        p=p-0.5*eps*q; grads+=L
        if np.log(rng.random())<H0-(0.5*q@q+0.5*p@p): acc+=1
        else: q=q0
        out[t]=q[0]
    return out, acc/n, grads
rng=np.random.default_rng(5)
print("part 1: dimension D"); 
for D in (10,100,400,1000):
    t=time.time()
    x,a=rwm(D, 400000, rng); e=ess_bm(x[20000:])
    y,ah,g=hmc(D, 20000, rng); eh=ess_bm(y[1000:])
    print(f"D={D:5d}  RWM acc {a:.2f} ESS/eval {e/380000:.5f}  (evals per indep {380000/e:7.1f}; D/evals {D/(380000/e):.2f})"
          f" | HMC acc {ah:.2f} ESS/grad {eh/(g*19000/20000):.4f} (grads per indep {(g*19000/20000)/eh:6.1f})  {time.time()-t:.0f}s")

print("part 2: SD ratio k")
rng=np.random.default_rng(3)
for k in (1,5,20,50):
    sd=np.array([float(k),1.0])          # widest/narrowest SD ratio = k; no preconditioning
    # RWM: step tuned to the narrow direction (best single isotropic step found by small search)
    best=None
    for s in (0.5,1.0,1.5,2.0,3.0):
        q=np.zeros(2); lp=0; n=400000; out=np.empty(n)
        for t in range(n):
            qp=q+s*rng.standard_normal(2); lpp=-0.5*np.sum((qp/sd)**2)
            if np.log(rng.random())<lpp-lp: q,lp=qp,lpp
            out[t]=q[0]
        e=ess_bm(out[20000:]); c=380000/e
        best=c if best is None else min(best,c)
    # HMC, identity mass matrix, eps tuned to narrow direction, trajectory ~ wide scale
    eps=0.9; q=np.zeros(2); n=20000; out=np.empty(n); g=0
    for t in range(n):
        p=rng.standard_normal(2); L=max(1,int(np.ceil(rng.uniform(0.8,1.6)*k/eps)))
        q0=q.copy(); H0=0.5*np.sum((q/sd)**2)+0.5*p@p
        p=p-0.5*eps*q/sd**2
        for l in range(L):
            q=q+eps*p
            if l<L-1: p=p-eps*q/sd**2
        p=p-0.5*eps*q/sd**2; g+=L
        if np.log(rng.random())>=H0-(0.5*np.sum((q/sd)**2)+0.5*p@p): q=q0
        out[t]=q[0]
    eh=ess_bm(out[1000:]); ch=g*19/20/eh
    print(f"SD ratio {k:3d}: RWM evals/indep (wide coord) {best:8.1f}  -> /k^2 {best/k**2:.2f} | HMC grads/indep {ch:7.1f} -> /k {ch/k:.2f}")
