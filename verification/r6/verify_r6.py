#!/usr/bin/env python3
"""Fail-closed interval checker for the frozen R6 conformal certificate.

Partial algebra/finite stages print OPEN. The final stage checks the hashes,
coverage and Arb inequalities of every separately computed interval batch,
then the nonlinear radii polynomial and lifted geometry before PROVED.
Use reproduce_r6.sh to recompute the batch receipts from source.
"""
import argparse
import hashlib
import json
from pathlib import Path
from flint import arb, arb_mat, ctx

Z = arb(0)

# Frozen inputs/interval batches for this specific certificate. Each stage
# script can regenerate its corresponding receipt from the centre and inverse.
EXPECTED_SHA256 = {
    'center_r6.json': '83018f0975e0cd61f4168e2de8bb1946032265869b5bfa168773c5920198aaba',
    'inverse_r6_M58_S30.npy': '8bb1b34aa882e1dbbe6943cba4726ba4eb18aefa4150cb734665a85b9f1026a9',
    'finite_batches_r6.json': '138e166004e655a3e76f3ceea5e6d48fc99549734d6de28f573246176044c97c',
    'residual_bound_r6.json': 'b111bf5b26279810d582125f0d04fa07f07614d3e0bb60cf0119049319db556a',
    'tail_inverse_bound_r6.json': '393643e94df5e6d32189911ae61c10664e3658e6fff06d961931ef7d07f0433a',
    'g_tail_batches_r6.json': '194efaeb46b2762d86532cc9b91d6c6b652eef1ead0a879ceb082ed3cdd98375',
    'g_far_bound_r6.json': '8af49da5adc3aa3801c898949ea4ab8ff55c7ffa0ce94b51cd601b4ba300e474',
    'shape_tail_batches_r6.json': '7331746c80180533ce527565064aa39749ceca22172ae407f5f9d6a3531ce8ef',
    'far_shape_bound_r6.json': '6cf4938b1aaad6a1f40f7acc11ae3e602b615b6a9d469cb319deff6595153478',
    # Source of the interval stages that produce the receipts above
    # (reproduce_r6.sh regenerates the receipts from exactly these files).
    'finite_interval_r6.py': 'dec54274aa7b52940ef6cf9253102b5d12393629b2bd40434e2b21a526dc060c',
    'tail_inverse_r6.py': 'ea14f859d44a9f643263b73b10a4f95548c7b603cf6e7e3c3b6e109164c183ab',
    'g_tail_interval_r6.py': '6799e59083a354db0178e06dcbbc03f468f8f1edfdc12978fb85a3e2a0dca211',
    'shape_tail_interval_r6.py': '05566777c507e8c37ecc5e3fd20cbfda2010cb1b8a1d8825bede55a77ecbbb3c',
    'far_shape_r6.py': '99f1da9274008d3122bb5e7fffbdf949465f9d530639935c548674a005fe2804',
}


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
    """Multiply by sum p[a] w^a (or conjugate); p degrees are 0 mod 4."""
    if not p:
        return {}
    deg = max(p)
    need(all(a>=0 and a%4==0 for a in p),'hol degrees not 0 mod 4')
    out = scale(f,p[deg])
    for a in range(deg-4,-1,-4):
        out = shift(out,4,bar)
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
        need(n>0 and n%4==2,'sine-sector failure')
        acc(out,(n,s),-v/2)
        acc(out,(-n,s),v/2)
    return out


def unsine(f):
    out = {}
    for (m,s),v in f.items():
        need(m%4==2,'unsine-sector failure')
        acc(out,(abs(m),s),-v if m>0 else v)
    return out


def K(f):
    out = {}
    for (n,s),v in f.items():
        need(n>0 and n%4==2 and s>=1,'K domain failure')
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


def hcoeff(c):
    out = {}
    for j,cj in c.items():
        for k,ck in c.items():
            acc(out,(j+k,0),cj*ck/2)
    return out


def dyadic(h):
    a,b=float.fromhex(h).as_integer_ratio()
    v=arb(a)/b
    need(v.is_exact(),'input not exact dyadic')
    return v


class Inverse:
    """Exact inverse of the R6 ball-principal infinite shape tail."""
    def __init__(self,M,S,rho,b,Ahat):
        self.M,self.S,self.rho,self.b,self.Ahat=M,S,rho,b,Ahat
        self.ms=list(range(2,M+1,4))
        self.js=list(range(1,M,4))
        self.rows=[(n,s) for n in self.ms for s in range(1,S+1)]+[(n,0) for n in self.ms]
        self.pos={k:i for i,k in enumerate(self.rows)}
        self.N=len(self.rows)
        self.alpha=b**3/2
        self.cw=[rho**n if s else self.omega(n-1) for n,s in self.rows]

    def omega(self,j):
        return self.alpha*(j+2)*self.rho**(j+1)

    def split(self,f):
        vf=[arb(0) for _ in self.rows]
        dg={};harm={}
        for (n,s),v in f.items():
            if (n,s) in self.pos:
                vf[self.pos[n,s]]=v
            elif s==0:
                need(n>self.M and n%4==2,'tail harmonic index invalid')
                harm[n]=v
            else:
                acc(dg,(n,s),v)
        eta={}
        run=arb(0)
        hi=max(harm,default=self.M)
        for n in range(hi,self.M,-4):
            # Before adding row n, run is eta_(n+3).
            if not run.is_zero():
                acc(dg,(n,1),self.b**3*(n+3)/(n+4)*run)
                acc(dg,(n,2),self.b**3/(n+4)*run)
            run += harm.get(n,Z)/(self.alpha*(n+1))
            eta[n-1]=run
            if n==self.M+4:
                vf[self.pos[self.M,0]] += self.alpha*(self.M+1)*run
                vf[self.pos[self.M,1]] += self.b**3*(self.M+3)/(self.M+4)*run
                vf[self.pos[self.M,2]] += self.b**3/(self.M+4)*run
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
                else:
                    raise ValueError(typ)
            out.append(sum((abs(image[i,q])*self.cw[i] for i in range(self.N)),Z)+self.tailnorm(eta,dg))
        return out


def make_algebra(data):
    M,S=data['M'],data['S']
    need(M%4==2 and S>=2,'invalid M/S')
    ms=list(range(2,M+1,4));js=list(range(1,M,4));ng=len(ms)*S
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
    Bj=hol_psi(abs2({(1,0):arb(1)},p),c)
    for j in sorted(c):
        cols.append(scale(unsine(add(scale(Aj,arb(j)),scale(Bj,-arb(1)/2))),arb(2)))
        Aj=shift(Aj,4);Bj=shift(Bj,4)
    return cols


def _read(name):
    return json.loads(Path(name).read_text())


def _interval(value):
    v=arb(value)
    need(v.is_finite() and v.upper().is_finite(),'nonfinite interval in receipt: '+str(value))
    return v


def final_gate(data,rows,g,c,p,V,F,rho,output):
    """Combine every interval stage with analytic majorants and geometry gates."""
    for name,expected in EXPECTED_SHA256.items():
        actual=hashlib.sha256(Path(name).read_bytes()).hexdigest()
        need(actual==expected,'stale or changed certificate input: '+name)
    need(data['M']==58 and data['S']==30,'unsupported centre')
    from finite_interval_r6 import load_inverse
    Ahat=load_inverse('inverse_r6_M58_S30.npy',len(rows))
    inv=Inverse(data['M'],data['S'],rho,c[1],Ahat)
    Anfin=max_upper(sum((abs(Ahat[i,j])*inv.cw[i] for i in range(inv.N)),Z)/rho**rows[j][0]
                    for j in range(inv.N))
    Anq=arb(1402314)/1000
    need(Anfin<Anq,'finite approximate-inverse norm')

    finite=_read('finite_batches_r6.json')
    cur=0;finite_z=[];finite_inv=[]
    for item in finite:
        need(item['start']==cur and item['stop']>cur and item['N']==inv.N,'finite batch coverage')
        cur=item['stop'];finite_z.append(_interval(item['Zfin_max']))
        finite_inv.append(_interval(item['inverse_defect_max']))
    need(cur==inv.N,'finite columns incomplete')
    inv_def=max_upper(finite_inv)
    need(inv_def<arb(1),'finite inverse not certified')
    Zfin=max_upper(finite_z)

    residual=_read('residual_bound_r6.json')
    need(residual['kind']=='residual' and residual['N']==inv.N,'residual receipt mismatch')
    Y=_interval(residual['Y']).upper()

    tail=_read('tail_inverse_bound_r6.json')
    Atail=_interval(tail['A_tail_upper']).upper()
    Atailq=arb(3164016)/1000000
    need(Atail<Atailq,'principal tail inverse bound')

    gc=_read('g_tail_batches_r6.json')
    cur=0;gz=[]
    for item in gc:
        need(item['start']==cur and item['stop']>cur and item['near_count']==705,'g tail batch coverage')
        cur=item['stop'];gz.append(_interval(item['max']))
    need(cur==705,'g near tail incomplete')
    Zgnear=max_upper(gz)
    gfar=_read('g_far_bound_r6.json')
    need(gfar['degree']==24 and gfar['near_count']==705,'g far support mismatch')
    Zgfar=max_upper([_interval(gfar['far_n']),_interval(gfar['far_s'])])

    sh=_read('shape_tail_batches_r6.json')
    cur=61;sz=[]
    for item in sh:
        need(item['first']==cur and item['last']>=cur and item['count']==(item['last']-cur)//4+1,
             'shape tail batch coverage')
        need((item['last']-cur)%4==0,'shape batch step')
        cur=item['last']+4;sz.append(_interval(item['max']))
    need(cur==401,'shape near tail incomplete')
    Zshnear=max_upper(sz)
    shfar=_read('far_shape_bound_r6.json')
    need(shfar['jstar']==401 and shfar['M']==58 and shfar['min_freq']>58,
         'far shape support failure')
    Zshfar=_interval(shfar['Zshape_far']).upper()
    Zbound=max_upper([Zfin,Zgnear,Zgfar,Zshnear,Zshfar])

    # H=g+|psi'|^2(Kg+Im(psi^2)/2) is quartic in (g,c).
    b=c[1];alpha=b**3/2
    P=sum((abs(v)*rho**a for a,v in p.items()),Z)
    C=sum((abs(v)*rho**j for j,v in c.items()),Z)
    Vn=norm(V,rho)
    theta0=1/(3*alpha*rho)
    theta1=1/(alpha*rho*rho)
    kap=arb(1)/24
    C2=Anq*max_upper([2*P*theta1*kap,
                      2*P*theta1*C*theta0+P*P*theta0*theta0/2+theta1*theta1*Vn])
    C3=Anq*max_upper([theta1*theta1*kap,
                      P*theta1*theta0*theta0+theta1*theta1*C*theta0])
    C4=Anq*theta1*theta1*theta0*theta0/2
    Yq=arb(53)/1000000
    Zq=arb(47)/100
    C2q=arb(251)/1000
    C3q=arb(3)/10000000
    C4q=arb(3)/1000000000000000
    need(Y<Yq and Zbound<Zq and C2<C2q and C3<C3q and C4<C4q,
         'one of the interval majorants failed')
    radius=arb(1)/5000
    polynomial=Yq+(Zq-1)*radius+C2q*radius**2+C3q*radius**3+C4q*radius**4
    derivative=Zq-1+2*C2q*radius+3*C3q*radius**2+4*C4q*radius**3
    need(polynomial<0 and derivative<0,'radii polynomial failed')

    # Analytic extension, strict star shape, and a nonzero c5 for the true zero.
    c1lo=b-radius/inv.omega(1)
    ext=arb(21)/20
    dcenter=sum((j*abs(v)*ext**(j-1) for j,v in c.items() if j>=5),Z)
    mcenter=sum((abs(v) for j,v in c.items() if j>=5),Z)
    derivative_tail=dcenter+radius*theta1
    map_tail=mcenter+radius*theta0
    extension_margin=c1lo-derivative_tail
    star=(derivative_tail+map_tail)/(c1lo-map_tail)
    c5lo=abs(c[5])-radius/inv.omega(5)
    need(ext<rho and c1lo>0 and extension_margin>0 and star<1 and c5lo>0,
         'lifted-domain geometry failed')

    hashes=EXPECTED_SHA256
    result={'status':'PROVED','bits':ctx.prec,'M':58,'S':30,'rho_weight':'11/10',
            'Y':Y.str(30),'Z_finite':Zfin.str(30),'Z_g_near':Zgnear.str(30),
            'Z_g_far':Zgfar.str(30),'Z_shape_near':Zshnear.str(30),
            'Z_shape_far':Zshfar.str(30),'Z':Zbound.str(30),
            'inverse_defect':inv_def.str(30),'A_finite':Anfin.str(30),'A_tail':Atail.str(30),
            'C2':C2.str(30),'C3':C3.str(30),'C4':C4.str(30),
            'radius':radius.str(30),'radii_polynomial':polynomial.str(30),
            'derivative_polynomial':derivative.str(30),
            'extension_margin':extension_margin.str(30),'star_defect':star.str(30),
            'c5_lower':c5lo.str(30),'sha256':hashes,
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print('PROVED',flush=True)
    print('certificate:',output,flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--centre',default='center_r6.json')
    parser.add_argument('--stage',choices=['algebra','finite','final'],default='algebra')
    parser.add_argument('--bits',type=int,default=128)
    parser.add_argument('--rho',type=str,default='11/10')
    parser.add_argument('--output',default='certificate_r6.json')
    args=parser.parse_args()
    ctx.prec=args.bits;ctx.threads=1
    data=json.load(open(args.centre))
    a,b=args.rho.split('/');rho=arb(int(a))/int(b)
    rows,g,c,p,V,F=make_algebra(data)
    if args.stage=='final':
        need(args.rho=='11/10','unsupported analytic weight')
        final_gate(data,rows,g,c,p,V,F,rho,args.output)
        return
    print('algebra',len(rows),'residual coefficients',len(F),'residual norm',norm(F,rho).str(15),flush=True)
    if args.stage=='algebra':
        print('OPEN: finite inverse and infinite tail checks pending',flush=True)
        return
    cols=finite_columns(rows,g,c,p,V,data['S'])
    N=len(rows)
    Mf=arb_mat([[col.get(k,Z) for col in cols] for k in rows])
    Ainv=Mf.inv()
    Ahat=arb_mat([[Ainv[i,j].mid() for j in range(N)] for i in range(N)])
    inv=Inverse(data['M'],data['S'],rho,c[1],Ahat)
    E=arb_mat(N,N)
    for i in range(N):E[i,i]=1
    E=E-Ahat*Mf
    err=max_upper(sum((abs(E[i,j])*inv.cw[i] for i in range(N)),Z)/inv.cw[j] for j in range(N))
    Y=inv.apply_many([F])[0].upper()
    print('finite inverse defect',err.str(15),'Y',Y.str(15),flush=True)
    print('OPEN: g and shape tail bounds pending',flush=True)


if __name__=='__main__':
    main()
