#!/usr/bin/env python3
"""Fail-closed Arb CAP verifier for the D6-symmetric Schiffer domain in R8.

The algebra and finite stages are diagnostics. The final stage checks every
frozen receipt, infinite-column estimate, nonlinear inequality, and geometry.
"""
import argparse
import hashlib
import json
from pathlib import Path
from flint import arb, arb_mat, ctx
import numpy as np

Z = arb(0)

M_ROOT = 3
STEP = 6
FROZEN_BUNDLE_SHA256 = 'c32ce586dafc332fdf8d05d74a90a98babe97a4275d86529946f482c06d5f21e'
FROZEN_FILES = [
    'center_r8_M45_S24.json', 'inverse_r8_M45_S24.npy',
    'finite_bound_r8.json', 'tail_inverse_bound_r8.json',
    'g_far_bound_r8.json', 'shape_near_r8_049_133.json',
    'shape_near_r8_139_217.json', 'far_shape_bound_r8.json',
    'finite_interval_r8.py', 'tail_inverse_r8.py', 'g_tail_r8.py',
    'shape_tail_r8.py', 'far_shape_r8.py',
    'check_identities_r8.py', 'identities_r8.json',
] + [f'g_near_r8_{a:03d}_{min(a+20,396):03d}.json'
     for a in range(0,396,20)]


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
    """Multiply by sum p[a] w^a (or conjugate); p degrees are 0 mod 6."""
    if not p:
        return {}
    deg = max(p)
    need(all(a>=0 and a%STEP==0 for a in p),'hol degrees not 0 mod 6')
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
    return {(n,0):v for n,v in powers.items() if not v.is_zero()}


def dyadic(h):
    a,b=float.fromhex(h).as_integer_ratio()
    v=arb(a)/b
    need(v.is_exact(),'input not exact dyadic')
    return v


class Inverse:
    """Inverse of the R8 ball-principal infinite shape tail."""
    def __init__(self,M,S,rho,b,Ahat):
        self.M,self.S,self.rho,self.b,self.Ahat=M,S,rho,b,Ahat
        self.ms=list(range(M_ROOT,M+1,STEP))
        self.js=list(range(1,1+STEP*len(self.ms),STEP))
        self.rows=[(n,s) for n in self.ms for s in range(1,S+1)]+[(n,0) for n in self.ms]
        self.pos={k:i for i,k in enumerate(self.rows)}
        self.N=len(self.rows)
        self.alpha=b**(M_ROOT+1)
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
        cols.append(unsine(add(scale(Aj,arb(2*j)),scale(Bj,-arb(M_ROOT)))))
        Aj=shift(Aj,STEP);Bj=shift(Bj,STEP)
    return cols


def shape_column(j,c,p,V):
    need(j>=1 and j%STEP==1, 'invalid shape column')
    Aj=shift(hol(sine(V),p,True),j-1)
    Bj=zw(abs2({(0,0):arb(1)},p))
    for _ in range(M_ROOT-1):
        Bj=hol_psi(Bj,c)
    Bj=shift(Bj,j-1)
    return unsine(add(scale(Aj,arb(2*j)),scale(Bj,-arb(M_ROOT))))


def principal_shape(j,b):
    """Ball-principal shape derivative for j=1 mod 6, j>=7."""
    need(j>=1+STEP and j%STEP==1, 'invalid principal shape index')
    n=j-M_ROOT-1
    q=shift(shift({(n,0):arb(1)},M_ROOT),M_ROOT,True)
    out={(j+M_ROOT-1,0):b**(M_ROOT+1)*(j+M_ROOT)}
    for key,value in q.items():
        acc(out,key,-b**(M_ROOT+1)*j*value)
    return out


def _read(name):
    return json.loads(Path(name).read_text())


def _interval(value):
    v=arb(value)
    need(v.is_finite() and v.upper().is_finite(),'nonfinite interval in receipt: '+str(value))
    return v


def frozen_digest():
    digest=hashlib.sha256()
    for name in FROZEN_FILES:
        data=Path(name).read_bytes()
        digest.update(name.encode('utf-8')+b'\0')
        digest.update(len(data).to_bytes(8,'big'))
        digest.update(data)
    return digest.hexdigest()


def final_gate(data,rows,g,c,p,V,F,rho,output):
    need(len(FROZEN_BUNDLE_SHA256)==64 and frozen_digest()==FROZEN_BUNDLE_SHA256,
         'stale or incomplete frozen R8 bundle')
    need(data['m']==3 and data['step']==6 and data['M']==45 and data['S']==24,
         'unsupported R8 centre')
    finite=_read('finite_bound_r8.json')
    need(finite['status']=='OPEN_COMPONENT' and finite['N']==len(rows)
         and finite['bits']>=128,'finite receipt mismatch')
    centre_hash=hashlib.sha256(Path('center_r8_M45_S24.json').read_bytes()).hexdigest()
    inverse_hash=hashlib.sha256(Path('inverse_r8_M45_S24.npy').read_bytes()).hexdigest()
    need(finite['centre_sha256']==centre_hash and finite['inverse_sha256']==inverse_hash,
         'finite receipt input mismatch')
    Anfin=_interval(finite['A_finite']).upper()
    inv_def=_interval(finite['inverse_defect']).upper()
    Y=_interval(finite['Y']).upper()
    Zfin=_interval(finite['Z_finite']).upper()
    need(Anfin<arb(964) and inv_def<arb(1),
         'finite approximate inverse not certified')
    tail=_read('tail_inverse_bound_r8.json')
    need(tail['status']=='OPEN_COMPONENT' and tail['M']==45 and tail['S']==24
         and tail['count']==100 and tail['far_start']==651,
         'principal tail inverse receipt mismatch')
    Atail=_interval(tail['A_tail_upper']).upper()
    need(Atail<arb(2303)/1000 and Anfin>Atail,
         'principal tail inverse or full-norm comparison failed')
    # The finite columns have been checked in one Arb pass. The g-tail file
    # names and [start,stop) fields must cover every one of the 396 near keys.
    far_g=_read('g_far_bound_r8.json')
    need(far_g['kind']=='far' and far_g['degree']==24
         and far_g['near_count']==396 and far_g['n_far']==75
         and far_g['s_far']==50,'g-far support mismatch')
    cur=0
    gnear=[]
    for a in range(0,396,20):
        item=_read(f'g_near_r8_{a:03d}_{min(a+20,396):03d}.json')
        need(item['kind']=='near' and item['start']==cur
             and item['stop']==min(a+20,396) and item['near_count']==396,
             'g-near column coverage gap')
        cur=item['stop']
        gnear.append(_interval(item['max']))
    need(cur==396,'g-near columns incomplete')
    Zgnear=max_upper(gnear)
    Zgfar=max_upper([_interval(far_g['far_n']),_interval(far_g['far_s'])])
    sh1=_read('shape_near_r8_049_133.json')
    sh2=_read('shape_near_r8_139_217.json')
    need(sh1['kind']=='near' and sh1['first']==49 and sh1['last']==133
         and sh1['count']==15 and sh2['kind']=='near'
         and sh2['first']==139 and sh2['last']==217 and sh2['count']==14,
         'shape-near column coverage gap')
    Zshnear=max_upper([_interval(sh1['max']),_interval(sh2['max'])])
    far_sh=_read('far_shape_bound_r8.json')
    need(far_sh['status']=='OPEN_COMPONENT' and far_sh['jstar']==223
         and far_sh['M']==45 and far_sh['min_freq']>45
         and _interval(far_sh['Q_identity_error'])<arb('1e-20'),
         'far-shape support or Q identity failed')
    Zshfar=_interval(far_sh['Z_shape_far']).upper()
    identities=_read('identities_r8.json')
    need(identities['status']=='IDENTITIES_PASS' and identities['bits']>=256
         and identities['K_cases']==12 and identities['radial_power_cases']==3
         and _interval(identities['ball_principal_max_error'])<arb('1e-60')
         and _interval(identities['Q_connection_max_error'])<arb('1e-60'),
         'identity receipt failed')
    Zbound=max_upper([Zfin,Zgnear,Zgfar,Zshnear,Zshfar])

    # H=g+|psi'|^2(Kg+Im psi^3) has degree five in (g,c).
    b=c[1]
    alpha=b**4
    P=sum((abs(x)*rho**a for a,x in p.items()),Z)
    C=sum((abs(x)*rho**j for j,x in c.items()),Z)
    Vn=norm(V,rho)
    theta0=1/(4*alpha*rho**2)
    theta1=1/(alpha*rho**3)
    kap=arb(1)/35
    Anq=arb(964)
    C2=Anq*(2*P*theta1*kap+3*P*P*C*theta0**2
             +6*P*C*C*theta1*theta0+theta1**2*Vn)
    C3=Anq*(theta1**2*kap+P*P*theta0**3
             +6*P*C*theta1*theta0**2+3*C*C*theta1**2*theta0)
    C4=Anq*(2*P*theta1*theta0**3+3*C*theta1**2*theta0**2)
    C5=Anq*theta1**2*theta0**3
    Yq=arb(153)/10000000
    Zq=arb(43)/100
    C2q=arb(1)/125
    C3q=arb(5)/10000000000
    C4q=arb(2)/1000000000000000000
    C5q=arb(3)/100000000000000000000000000
    need(Y<Yq and Zbound<Zq and C2<C2q and C3<C3q
         and C4<C4q and C5<C5q,'CAP majorant exceeded its rational threshold')
    radius=arb(1)/30000
    polynomial=Yq+(Zq-1)*radius+C2q*radius**2+C3q*radius**3+C4q*radius**4+C5q*radius**5
    derivative=Zq-1+2*C2q*radius+3*C3q*radius**2+4*C4q*radius**3+5*C5q*radius**4
    need(polynomial<0 and derivative<0,'radii polynomial failed')

    # Univalence on an enlarged disk and a global, unsampled convexity margin.
    ext=arb(21)/20
    need(ext<rho,'analytic extension radius exceeds coefficient radius')
    extension_margin=b-sum((j*abs(x)*ext**(j-1) for j,x in c.items() if j>=7),Z)-radius*theta1
    derivative_lower=b-sum((j*abs(x) for j,x in c.items() if j>=7),Z)-radius*theta1
    second_perturb=radius/(alpha*rho**2*arb.const_e()*rho.log())
    second_upper=sum((j*(j-1)*abs(x) for j,x in c.items() if j>=7),Z)+second_perturb
    curvature_lower=1-second_upper/derivative_lower
    c7_lower=abs(c[7])-radius/(alpha*10*rho**9)
    need(b>0 and extension_margin>0 and derivative_lower>0
         and curvature_lower>0 and c7_lower>0,
         'analytic non-disc convex geometry failed')

    out={'status':'PROVED','dimension':8,'bits':ctx.prec,'M':45,'S':24,
         'rho_weight':'11/10','Y':Y.str(30),'Z':Zbound.str(30),
         'Z_finite':Zfin.str(30),'Z_g_near':Zgnear.str(30),
         'Z_g_far':Zgfar.str(30),'Z_shape_near':Zshnear.str(30),
         'Z_shape_far':Zshfar.str(30),'A_finite':Anfin.str(30),
         'A_tail':Atail.str(30),'inverse_defect':inv_def.str(30),
         'C2':C2.str(30),'C3':C3.str(30),'C4':C4.str(30),'C5':C5.str(30),
         'radius':radius.str(30),'radii_polynomial':polynomial.str(30),
         'derivative_polynomial':derivative.str(30),
         'extension_margin':extension_margin.str(30),
         'curvature_lower':curvature_lower.str(30),
         'c7_lower':c7_lower.str(30),
         'frozen_bundle_sha256':FROZEN_BUNDLE_SHA256,
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(output).write_text(json.dumps(out,indent=2)+'\n')
    print('PROVED',flush=True)
    print('certificate:',output,flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--centre', default='center_r8.json')
    parser.add_argument('--stage', choices=['algebra', 'finite', 'final'], default='algebra')
    parser.add_argument('--bits', type=int, default=128)
    parser.add_argument('--rho', type=str, default='11/10')
    parser.add_argument('--inverse-file', help='freeze the finite inverse as binary64 .npy')
    parser.add_argument('--output',default='certificate_r8.json')
    args = parser.parse_args()
    ctx.prec = args.bits
    ctx.threads = 1
    data = json.loads(Path(args.centre).read_text())
    a, b = args.rho.split('/')
    rho = arb(int(a))/int(b)
    rows, g, c, p, V, F = make_algebra(data)
    print('algebra', len(rows), 'residual coefficients', len(F),
          'residual norm', norm(F, rho).str(15), flush=True)
    if args.stage == 'algebra':
        print('OPEN: finite inverse and infinite tail checks pending', flush=True)
        return
    if args.stage == 'final':
        need(args.rho=='11/10' and args.bits>=128,
             'unsupported weight or precision for final gate')
        final_gate(data,rows,g,c,p,V,F,rho,args.output)
        return
    cols = finite_columns(rows, g, c, p, V, data['S'])
    N = len(rows)
    Mf = arb_mat([[col.get(k, Z) for col in cols] for k in rows])
    Ainv = Mf.inv()
    if args.inverse_file:
        frozen = np.array([[float(Ainv[i,j].mid()) for j in range(N)] for i in range(N)])
        np.save(args.inverse_file, frozen)
        Ahat = arb_mat([[dyadic(float(frozen[i,j]).hex()) for j in range(N)] for i in range(N)])
    else:
        Ahat = arb_mat([[Ainv[i,j].mid() for j in range(N)] for i in range(N)])
    inv = Inverse(data['M'], data['S'], rho, c[1], Ahat)
    E = arb_mat(N,N)
    for i in range(N):
        E[i,i] = 1
    E = E-Ahat*Mf
    err = max_upper(sum((abs(E[i,j])*inv.cw[i] for i in range(N)),Z)/inv.cw[j]
                    for j in range(N))
    Y = inv.apply_many([F])[0].upper()
    print('finite inverse defect', err.str(15), 'Y', Y.str(15), flush=True)
    print('OPEN: g and shape tail bounds pending', flush=True)


if __name__ == '__main__':
    main()
