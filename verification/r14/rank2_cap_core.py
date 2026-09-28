"""Common Arb Zernike algebra for I2(m), m in {3,4,6}, step=2m.

This module only supplies algebra and the principal-tail inverse. A dimension
specific verifier must certify every finite and infinite inequality.
"""
import hashlib
import json
from pathlib import Path
from flint import arb, arb_mat, ctx

Z=arb(0)
M_ROOT=3
STEP=6
BOUNDARY_SCALE=arb(1)


def configure(data):
    global M_ROOT,STEP,BOUNDARY_SCALE
    m=int(data['m'])
    step=int(data['step'])
    need(m in (3,4,6) and step==2*m,'unsupported Weyl symmetry')
    scale=dyadic(data.get('boundary_scale','0x1.0000000000000p+0'))
    need(scale>0,'boundary field scale must be positive')
    M_ROOT=m
    STEP=step
    BOUNDARY_SCALE=scale

def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def acc(out, key, value):
    if not value.is_zero():
        out[key] = out.get(key, Z) + value


def add(*seq):
    out = {}
    for f in seq:
        for key, value in f.items():
            acc(out, key, value)
    return out


def scale(f, a):
    return {key:a*value for key,value in f.items() if not value.is_zero()}


def zw(f, bar=False):
    """Multiply a disk-polynomial expansion by w or conjugate(w)."""
    out = {}
    sign = -1 if bar else 1
    for (m,s),v in f.items():
        n = abs(m)
        d = n+2*s+1
        if sign*m >= 0:
            acc(out,(m+sign,s),v*(n+s+1)/d)
            if s:
                acc(out,(m+sign,s-1),v*s/d)
        else:
            acc(out,(m+sign,s),v*(n+s)/d)
            acc(out,(m+sign,s+1),v*(s+1)/d)
    return out


def shift(f,k,bar=False):
    for _ in range(k):
        f = zw(f,bar)
    return f


def hol(f,p,bar=False):
    """Multiply by sum p[a] w^a (or conjugate); p degrees are 0 mod STEP."""
    if not p:
        return {}
    deg = max(p)
    need(all(a>=0 and a%STEP==0 for a in p),'hol degrees not 0 mod STEP')
    out = scale(f,p[deg])
    for a in range(deg-STEP,-1,-STEP):
        out = shift(out,STEP,bar)
        if a in p:
            out = add(out,scale(f,p[a]))
    return out


def hol_psi(f,c):
    return zw(hol(f,{j-1:v for j,v in c.items()}))


def abs2(f,p):
    return hol(hol(f,p,True),p)


def sine(f):
    out = {}
    for (n,s),v in f.items():
        need(n>0 and n%STEP==M_ROOT,'sine-sector failure')
        acc(out,(n,s),-v/2)
        acc(out,(-n,s),v/2)
    return out


def unsine(f):
    out = {}
    for (m,s),v in f.items():
        need(m%STEP==M_ROOT,'unsine-sector failure')
        acc(out,(abs(m),s),-v if m>0 else v)
    return out


def K(f):
    out = {}
    for (n,s),v in f.items():
        need(n>0 and n%STEP==M_ROOT and s>=1,'K domain failure')
        d=n+2*s
        acc(out,(n,s-1),v/(4*d*(d+1)))
        acc(out,(n,s),-v/(2*d*(d+2)))
        acc(out,(n,s+1),v/(4*(d+1)*(d+2)))
    return out


def norm(f,rho):
    return sum((abs(v)*rho**abs(m) for (m,s),v in f.items()),Z)


def max_upper(seq):
    """Return a rigorous upper bound, never a midpoint comparison.

    Fail closed exactly like R^4's upper_max: a NaN or infinite interval (or
    an infinite upper endpoint) anywhere in the sequence aborts verification
    instead of being silently dropped by max().
    """
    bounds=[]
    for v in seq:
        need(v.is_finite(),'nonfinite interval in maximum')
        u=v.upper()
        need(u.is_finite(),'nonfinite upper endpoint in maximum')
        bounds.append(u)
    return max(bounds,default=Z)


def argmax_upper(seq):
    """Index of the largest upper endpoint (first one on ties), fail closed.

    Like max_upper, a NaN or infinite interval (or infinite upper endpoint)
    anywhere in the sequence aborts instead of being skipped by a comparison.
    """
    best=None
    for i,v in enumerate(seq):
        need(v.is_finite(),'nonfinite interval in maximum')
        u=v.upper()
        need(u.is_finite(),'nonfinite upper endpoint in maximum')
        if best is None or u>best_u:
            best,best_u=i,u
    need(best is not None,'empty maximum')
    return best


def hcoeff(c):
    powers={0:arb(1)}
    for _ in range(M_ROOT):
        next_power={}
        for a,aa in powers.items():
            for j,cj in c.items():
                next_power[a+j]=next_power.get(a+j,Z)+aa*cj
        powers=next_power
    return {(n,0):v/BOUNDARY_SCALE for n,v in powers.items() if not v.is_zero()}


def dyadic(h):
    a,b=float.fromhex(h).as_integer_ratio()
    v=arb(a)/b
    need(v.is_exact(),'input not exact dyadic')
    return v


class Inverse:
    """Inverse of the ball-principal infinite shape tail (m=3,4,6)."""
    def __init__(self,M,S,rho,b,Ahat):
        self.M,self.S,self.rho,self.b,self.Ahat=M,S,rho,b,Ahat
        self.ms=list(range(M_ROOT,M+1,STEP))
        self.js=list(range(1,1+STEP*len(self.ms),STEP))
        self.rows=[(n,s) for n in self.ms for s in range(1,S+1)]+[(n,0) for n in self.ms]
        self.pos={k:i for i,k in enumerate(self.rows)}
        self.N=len(self.rows)
        self.alpha=b**(M_ROOT+1)/BOUNDARY_SCALE
        self.cw=[rho**n if s else self.omega(n-M_ROOT+1) for n,s in self.rows]
        self.radial_cache={}

    def omega(self,j):
        return self.alpha*(j+M_ROOT)*self.rho**(j+M_ROOT-1)

    def radial_coeff(self,n):
        if n not in self.radial_cache:
            f=shift(shift({(n,0):arb(1)},M_ROOT),M_ROOT,True)
            self.radial_cache[n]={s:v for (k,s),v in f.items() if k==n}
            need(set(self.radial_cache[n])==set(range(M_ROOT+1)),
                 'radial-power support failure')
            need((arb(n+1)-arb(n+M_ROOT+1)*self.radial_cache[n][0]).contains(0),
                 'radial-power harmonic coefficient failure')
        return self.radial_cache[n]

    def split(self,f):
        vf=[arb(0) for _ in self.rows]
        dg={};harm={}
        for (n,s),v in f.items():
            if (n,s) in self.pos:
                vf[self.pos[n,s]]=v
            elif s==0:
                need(n>self.M and n%STEP==M_ROOT,'tail harmonic index invalid')
                harm[n]=v
            else:
                acc(dg,(n,s),v)
        eta={}
        run=arb(0)
        hi=max(harm,default=self.M)
        for n in range(hi,self.M,-STEP):
            # Before adding row n, run is eta_(n+m+1).
            if not run.is_zero():
                coeff=self.radial_coeff(n)
                for s in range(1,M_ROOT+1):
                    acc(dg,(n,s),self.alpha*(n+M_ROOT+1)*coeff[s]*run)
            run += harm.get(n,Z)/(self.alpha*(n+1))
            eta[n-M_ROOT+1]=run
            if n==self.M+STEP:
                vf[self.pos[self.M,0]] += self.alpha*(self.M+1)*run
                coeff=self.radial_coeff(self.M)
                for s in range(1,M_ROOT+1):
                    vf[self.pos[self.M,s]] += self.alpha*(self.M+M_ROOT+1)*coeff[s]*run
        return vf,eta,dg

    def tailnorm(self,eta,dg):
        return sum((self.omega(j)*abs(v) for j,v in eta.items()),Z)+norm(dg,self.rho)

    def apply_many(self,fs,subtract=None):
        pieces=[self.split(f) for f in fs]
        matrix=arb_mat([[item[0][i] for item in pieces] for i in range(self.N)])
        image=self.Ahat*matrix
        out=[]
        for q,(_,eta,dg) in enumerate(pieces):
            if subtract is not None:
                typ,idx=subtract[q]
                if typ=='finite':
                    image[idx,q]-=1
                elif typ=='shape':
                    eta[idx]=eta.get(idx,Z)-1
                elif typ=='g':
                    acc(dg,idx,-arb(1))
                else:
                    raise ValueError(typ)
            out.append(sum((abs(image[i,q])*self.cw[i] for i in range(self.N)),Z)+self.tailnorm(eta,dg))
        return out


def make_algebra(data):
    configure(data)
    M,S=data['M'],data['S']
    need(data.get('m')==M_ROOT and data.get('step')==STEP and M%STEP==M_ROOT and S>=3,
         'invalid m/step/M/S')
    ms=list(range(M_ROOT,M+1,STEP))
    js=list(range(1,1+STEP*len(ms),STEP))
    ng=len(ms)*S
    need(len(data['g'])==ng and len(data['c'])==len(js),'centre dimension mismatch')
    rows=[(n,s) for n in ms for s in range(1,S+1)]+[(n,0) for n in ms]
    g={k:dyadic(v) for k,v in zip(rows[:ng],data['g'])}
    c={j:dyadic(v) for j,v in zip(js,data['c'])}
    p={j-1:j*v for j,v in c.items()}
    V=add(K(g),hcoeff(c))
    F=add(g,unsine(abs2(sine(V),p)))
    return rows,g,c,p,V,F


def finite_columns(rows,g,c,p,V,S):
    ng=len(g);cols=[]
    for key in rows[:ng]:
        e={key:arb(1)}
        cols.append(add(e,unsine(abs2(sine(K(e)),p))))
    Aj=hol(sine(V),p,True)
    Bj=zw(abs2({(0,0):arb(1)},p))
    for _ in range(M_ROOT-1):
        Bj=hol_psi(Bj,c)
    for j in sorted(c):
        cols.append(unsine(add(scale(Aj,arb(2*j)),scale(Bj,-arb(M_ROOT)/BOUNDARY_SCALE))))
        Aj=shift(Aj,STEP);Bj=shift(Bj,STEP)
    return cols


def shape_column(j,c,p,V):
    need(j>=1 and j%STEP==1, 'invalid shape column')
    Aj=shift(hol(sine(V),p,True),j-1)
    Bj=zw(abs2({(0,0):arb(1)},p))
    for _ in range(M_ROOT-1):
        Bj=hol_psi(Bj,c)
    Bj=shift(Bj,j-1)
    return unsine(add(scale(Aj,arb(2*j)),scale(Bj,-arb(M_ROOT)/BOUNDARY_SCALE)))


def principal_shape(j,b):
    """Ball-principal shape derivative for j=1 mod STEP, j>=1+STEP."""
    need(j>=1+STEP and j%STEP==1, 'invalid principal shape index')
    n=j-M_ROOT-1
    q=shift(shift({(n,0):arb(1)},M_ROOT),M_ROOT,True)
    out={(j+M_ROOT-1,0):b**(M_ROOT+1)*(j+M_ROOT)/BOUNDARY_SCALE}
    for key,value in q.items():
        acc(out,key,-b**(M_ROOT+1)*j*value/BOUNDARY_SCALE)
    return out

