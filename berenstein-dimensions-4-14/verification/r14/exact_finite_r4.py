"""Exact algebra plus Arb enclosure for finite R4 divided fixed-disc equations.

This evaluates finite centres and their finite columns without quadrature.
It is not an infinite-dimensional proof or a standalone PROVED verifier.
The disk-polynomial multiplication recurrence is given in Section 3.2.
"""
from proof_guard import require
import argparse
import json
import hashlib
from flint import arb, arb_mat, ctx

ZERO=arb(0)


def finite_upper(v):
    require((v.is_finite()), 'non-finite Arb enclosure')
    u=v.upper()
    require((u.is_finite()), 'non-finite Arb upper endpoint')
    return u


def maximum_upper(values):
    best=arb(0)
    for v in values:
        u=finite_upper(v)
        if u>best:best=u
    return best


def dyadic(s):
    num,den=float.fromhex(s).as_integer_ratio()
    return arb(num)/arb(den)


def acc(d,k,v):
    if not v.is_zero():d[k]=d.get(k,ZERO)+v


def add(*items):
    d={}
    for item in items:
        for k,v in item.items():acc(d,k,v)
    return d


def scale(d,a):return {k:a*v for k,v in d.items() if not v.is_zero()}


def z(f,conjugate=False):
    out={}; sign=-1 if conjugate else 1
    for (m,s),v in f.items():
        n=abs(m); den=n+2*s+1
        if sign*m>=0:
            acc(out,(m+sign,s),v*(n+s+1)/den)
            if s:acc(out,(m+sign,s-1),v*s/den)
        else:
            acc(out,(m+sign,s),v*(n+s)/den)
            acc(out,(m+sign,s+1),v*(s+1)/den)
    return out


def poly(f,p,conjugate=False):
    require((p and all(e>=0 and e%2==0 for e in p)), 'exact_finite_r4.py:63: proof gate failed')
    degree=max(p); out=scale(f,p[degree])
    for e in range(degree-2,-1,-2):
        out=z(z(out,conjugate),conjugate)
        if e in p:out=add(out,scale(f,p[e]))
    return out


def abs2(f,p):return poly(poly(f,p,True),p)


def sine(f):
    out={}
    for (m,s),v in f.items():
        require((m>0 and m%2==1), 'exact_finite_r4.py:77: proof gate failed')
        acc(out,(m,s),-v/2)
        acc(out,(-m,s),v/2)
    return out


def unsine(f):
    out={}
    for (m,s),v in f.items():
        require((m%2==1), 'exact_finite_r4.py:86: proof gate failed')
        acc(out,(abs(m),s),-v if m>0 else v)
    return out


def KD(f):
    out={}
    for (m,s),v in f.items():
        require((m>0 and m%2==1 and s>=0), 'exact_finite_r4.py:94: proof gate failed')
        if s==0:
            a=v/(4*(m+1)*(m+2))
            acc(out,(m,0),-a)
            acc(out,(m,1),a)
        else:
            d=m+2*s
            acc(out,(m,s-1),v/(4*d*(d+1)))
            acc(out,(m,s),-v/(2*d*(d+2)))
            acc(out,(m,s+1),v/(4*(d+1)*(d+2)))
    return out


def conv(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():acc(out,i+j,x*y)
    return out


def laur_p(c):return {j-1:j*v for j,v in c.items()}


def laur_char(j):return {k:arb(1) for k in range(-(j-1),j,2)}


def laur_h(c):
    return add(*(scale(laur_char(j),v) for j,v in c.items()))


def laur_q(g):
    return add(*(scale(laur_char(m),v/(2*(m+1))) for (m,s),v in g.items() if s==0))


def laur_a2(p):return conv(p,{-k:v for k,v in p.items()})


def boundary(q,h,a2):return add(conv(q,q),scale(conv(conv(h,h),a2),-1))


def to_char(l):
    require((all(k%2==0 for k in l)), 'exact_finite_r4.py:135: proof gate failed')
    kmax=max((abs(k) for k in l),default=0)
    return {k+1:l.get(k,ZERO)-l.get(k+2,ZERO) for k in range(0,kmax+1,2)}


def get_center(path):
    with open(path) as f:d=json.load(f)
    M,S=d['M'],d['S']; ms=list(range(1,M+1,2))
    keys=[(m,s) for m in ms for s in range(S+1)]
    require((len(keys)==len(d['g_hex']) and len(ms)==len(d['c_hex'])), 'exact_finite_r4.py:144: proof gate failed')
    g={k:dyadic(v) for k,v in zip(keys,d['g_hex'])}
    c={j:dyadic(v) for j,v in zip(ms,d['c_hex'])}
    return d,g,c,keys,ms


def evaluate(g,c):
    p=laur_p(c)
    F=add(g,unsine(abs2(sine(KD(g)),p)))
    q,h,a2=laur_q(g),laur_h(c),laur_a2(p)
    B=to_char(boundary(q,h,a2))
    return F,B


def columns(g,c,keys,ms,start=0,end=None):
    p=laur_p(c)
    q,h,a2=laur_q(g),laur_h(c),laur_a2(p)
    V=sine(KD(g))
    out=[]
    if end is None:end=len(keys)+len(ms)
    for idx in range(start,end):
        if idx<len(keys):
            k=keys[idx]; e={k:arb(1)}
            Fi=add(e,unsine(abs2(sine(KD(e)),p)))
            if k[1]==0:
                dq=scale(laur_char(k[0]),arb(1)/(2*(k[0]+1)))
                Bi=to_char(scale(conv(q,dq),2))
            else:Bi={}
        else:
            j=ms[idx-len(keys)]
            dp={j-1:arb(j)}
            mixed=add(poly(poly(V,p,True),dp),poly(poly(V,dp,True),p))
            Fi=unsine(mixed)
            dh=laur_char(j)
            da2=add(conv(dp,{-k:v for k,v in p.items()}),
                    conv(p,{-k:v for k,v in dp.items()}))
            Bi=to_char(add(scale(conv(conv(h,dh),a2),-2),
                           scale(conv(conv(h,h),da2),-1)))
        out.append((Fi,Bi))
    return out


def norm(f,b,rho):
    return (sum((rho**m*abs(v) for (m,s),v in f.items()),ZERO)
            +sum((rho**(j-1)*abs(v) for j,v in b.items()),ZERO))


def apply_tail_norm(fi,bi,keys,ms,b,rho,q0):
    ftail={k:v for k,v in fi.items() if k not in keys}
    btail={j:v for j,v in bi.items() if j not in ms}
    rhs=dict(btail)
    coupling=to_char(scale(conv(q0,laur_q(fi)),2))
    for j,v in coupling.items():
        if j not in ms:rhs[j]=rhs.get(j,ZERO)-v
    maxj=max(rhs,default=ms[-1])
    run=arb(0);ctail={}
    for j in range(maxj if maxj%2 else maxj-1,ms[-1],-2):
        run+=rhs.get(j,ZERO)/(b**3*(j+2))
        ctail[j]=-run
    tailnorm=(sum((rho**m*abs(v) for (m,s),v in ftail.items()),ZERO)
              +sum((b**3*(j+2)*rho**(j-1)*abs(v) for j,v in ctail.items()),ZERO))
    return tailnorm,ctail.get(ms[-1]+2,ZERO)


def main(args):
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,keys,ms=get_center(args.center)
    F,B=evaluate(g,c)
    rho=arb(11)/10
    print('dimension',len(keys)+len(ms),'interior support',len(F),'boundary support',len(B))
    print('full exact-expanded residual norm upper',finite_upper(norm(F,B,rho)))
    print('finite interior max upper',maximum_upper(abs(F.get(k,ZERO)) for k in keys))
    print('finite boundary max upper',maximum_upper(abs(B.get(j,ZERO)) for j in ms))
    if args.tail_column_batch:
        start,end=map(int,args.tail_column_batch.split(':'))
        require((0<=start<end<=len(keys)+len(ms) and args.output), 'exact_finite_r4.py:219: proof gate failed')
        cols=columns(g,c,keys,ms,start,end)
        weights=[rho**m for m,s in keys]+[c[1]**3*(j+2)*rho**(j-1) for j in ms]
        q0=laur_q(g)
        edge_norm=arb(54778)/1000 # certified by verify_finite_inverse_r4.py
        vals=[]
        for idx,(fi,bi) in enumerate(cols,start):
            tailnorm,first=apply_tail_norm(fi,bi,keys,ms,c[1],rho,q0)
            finite_cross=edge_norm*c[1]**3*(ms[-1]+2)*abs(first)
            vals.append(finite_upper((tailnorm+finite_cross)/weights[idx]))
        maxval=maximum_upper(vals);argmax=start+vals.index(maxval)
        top_m,top_e=maxval.man_exp()
        payload={'center':args.center,'center_sha256':hashlib.sha256(open(args.center,'rb').read()).hexdigest(),
                 'bits':args.bits,'start':start,'end':end,'dimension':len(weights),
                 'max_upper_dyadic':[int(top_m),int(top_e)],'argmax':argmax,
                 'column_upper_dyadic':[[int(q) for q in v.man_exp()] for v in vals]}
        with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
        print('tail finite-column batch',start,end,'max',maxval,'argmax',argmax)
        return
    if args.batch:
        start,end=map(int,args.batch.split(':'))
        require((0<=start<end<=len(keys)+len(ms) and args.output), 'exact_finite_r4.py:240: proof gate failed')
        cols=columns(g,c,keys,ms,start,end)
        entries=[]
        for j,(fi,bi) in enumerate(cols,start):
            for i,k in enumerate(keys):
                v=fi.get(k,ZERO)
                if not v.is_zero():
                    m,e=v.mid().man_exp();r,f=v.rad().man_exp()
                    entries.append([i,j,int(m),int(e),int(r),int(f)])
            for t,row_j in enumerate(ms):
                v=bi.get(row_j,ZERO)
                if not v.is_zero():
                    m,e=v.mid().man_exp();r,f=v.rad().man_exp()
                    entries.append([len(keys)+t,j,int(m),int(e),int(r),int(f)])
        payload={'center':args.center,'center_sha256':hashlib.sha256(open(args.center,'rb').read()).hexdigest(),
                 'bits':args.bits,'start':start,'end':end,'dimension':len(keys)+len(ms),'entries':entries}
        with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
        print('wrote batch',start,end,'entries',len(entries),args.output)
        return
    if args.columns:
        cols=columns(g,c,keys,ms)
        rows=keys+[(j,0) for j in ms]
        matrix=arb_mat([[col[0].get(k,ZERO) if i<len(keys) else col[1].get(k[0],ZERO)
                         for col in cols] for i,k in enumerate(rows)])
        print('finite Jacobian',matrix.nrows(),matrix.ncols())
        if args.interval_inverse:
            N=len(rows)
            enclosed=matrix.inv()
            A=arb_mat([[enclosed[i,j].mid() for j in range(N)] for i in range(N)])
            E=arb_mat([[arb(int(i==j)) for j in range(N)] for i in range(N)])-A*matrix
            rho=arb(11)/10; b=c[1]
            win=[rho**m for m,s in keys]+[b**3*(j+2)*rho**(j-1) for j in ms]
            wout=[rho**m for m,s in keys]+[rho**(j-1) for j in ms]
            normA=maximum_upper(sum((win[i]*abs(A[i,j]) for i in range(N)),ZERO)/wout[j]
                                for j in range(N))
            defect=maximum_upper(sum((win[i]*abs(E[i,j]) for i in range(N)),ZERO)/win[j]
                                 for j in range(N))
            if not (defect<arb(1)):
                raise RuntimeError(f'finite inverse defect >=1: {defect}')
            print('CERTIFIED FINITE INVERSE; norm A upper',normA,'defect upper',defect)
        if args.numeric_diagnostics:
            import numpy as np
            J=np.array([[float(matrix[i,j].mid()) for j in range(len(rows))]
                        for i in range(len(rows))])
            inverse=np.linalg.inv(J)
            b=float(c[1].mid()); rho_f=1.1
            wout=np.array([rho_f**m for m,s in keys]+[rho_f**(j-1) for j in ms])
            win=np.array([rho_f**m for m,s in keys]
                         +[b**3*(j+2)*rho_f**(j-1) for j in ms])
            normA=np.max(np.sum(np.abs(inverse)*win[:,None],axis=0)/wout)
            defect=np.max(np.sum(np.abs(np.eye(len(rows))-inverse@J)*wout[:,None],axis=0)/wout)
            print('finite weighted inverse 1-norm',normA,'midpoint inverse defect',defect)
            if args.numeric_diagnostics!='-':
                np.savez(args.numeric_diagnostics,J=J,inverse=inverse,win=win,wout=wout,
                         center=args.center)
        if args.output:
            with open(args.output,'w') as f:
                json.dump({'center':args.center,'bits':args.bits,'dimension':len(rows),
                           'residual_norm_upper':str(finite_upper(norm(F,B,rho)))},f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--center',default='finite_center_r4_M11_S10.json')
    p.add_argument('--bits',type=int,default=128)
    p.add_argument('--columns',action='store_true')
    p.add_argument('--numeric-diagnostics',default='')
    p.add_argument('--interval-inverse',action='store_true')
    p.add_argument('--batch',default='')
    p.add_argument('--tail-column-batch',default='')
    p.add_argument('--output',default='')
    main(p.parse_args())
