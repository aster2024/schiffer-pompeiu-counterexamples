"""Arb support/Wiener bound for A_f P_f J on all high-source columns.

Uses a certified per-output-column profile of the finite inverse. This is
one cross block only; it does not bound the full inverse or derivative defect.
"""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,q_of,to_char
from exact_finite_r4 import ZERO


def upper(x):
    require((x.is_finite()), 'bound_finite_cross_source_rank2.py:15: proof gate failed')
    y=x.upper()
    require((y.is_finite()), 'bound_finite_cross_source_rank2.py:17: proof gate failed')
    return y


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('profile')
    ap.add_argument('--bits',type=int,default=160)
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m,M,S=d['m'],d['M'],d['S'];rho=arb(102)/100;b=c[1]
    sha=hashlib.sha256(open(a.center,'rb').read()).hexdigest()
    with open(a.profile) as f:p=json.load(f)
    require((p['center_sha256']==sha and p['dimension']==len(keys)+len(js)), 'bound_finite_cross_source_rank2.py:32: proof gate failed')
    require(((p['rho_num'],p['rho_den'])==(102,100)), 'bound_finite_cross_source_rank2.py:33: proof gate failed')
    prof=[upper(arb(int(man))*arb(2)**int(exp))
          for man,exp in p['profile_upper_dyadic']]
    ng=len(keys)
    source={n:[prof[i*(S+1)+s] for s in range(S+1)]
            for i,n in enumerate(ms)}
    cumulative={}
    for n,row in source.items():
        vals=[arb(0)]*(S+2);best=arb(0)
        for s in range(S,-1,-1):
            if row[s]>best:best=row[s]
            vals[s]=best
        cumulative[n]=vals
    boundary={n:prof[ng+i] for i,n in enumerate(ms)}
    pp={j-1:j*v for j,v in c.items()}
    by_degree={}
    for ea,va in pp.items():
        for eb,vb in pp.items():
            degree=ea+eb
            by_degree[degree]=by_degree.get(degree,ZERO)+abs(va)*abs(vb)*rho**degree
    degrees=sorted(by_degree)
    dmax=degrees[-1]
    qchar=to_char(q_of(g,m),m)
    nmax=max(M+dmax,2*M+m)
    maxval=arb(0);best=None;bestparts=None;checked=0
    cache={}
    def rowmax(n,degree,smin):
        key=(n,degree,smin)
        if key in cache:return cache[key]
        val=arb(0)
        for outn in ms:
            if abs(outn-n)>degree:continue
            cand=cumulative[outn][min(smin,S+1)]
            if cand>val:val=cand
        cache[key]=val
        return val
    for n in range(m,nmax+1,2*m):
        for s in range(S+2+dmax):
            if n<=M and s<=S:continue
            if s==0:kappa=arb(1)/(2*(n+1)*(n+2))
            else:
                D=n+2*s;kappa=arb(1)/(D*(D+2))
            inside=arb(0)
            for degree in degrees:
                smin=max(0,s-1-degree)
                if smin>S:continue
                inside+=by_degree[degree]*rowmax(n,degree,smin)
            inside*=kappa
            bound=inside
            boundary_part=arb(0)
            if s==0:
                for qa,qv in qchar.items():
                    for outn in range(abs(qa-n)+m,qa+n-m+1,2*m):
                        if outn not in boundary:continue
                        boundary_part+=abs(qv)*boundary[outn]*rho**(outn-m-n)/(b*(n+1))
                bound+=boundary_part
            bound=upper(bound)
            if bound>maxval:
                maxval=bound;best=(n,s);bestparts=(upper(inside),upper(boundary_part))
            checked+=1
    require((maxval<arb(1)/5), 'bound_finite_cross_source_rank2.py:93: proof gate failed')
    print('FINITE_CROSS_SOURCE_BOUND_ONLY','m',m,'checked',checked,
          'p degree',dmax,'maximum upper',maxval,'argmax',best,
          'parts',bestparts,'cache entries',len(cache))


if __name__=='__main__':
    main()
