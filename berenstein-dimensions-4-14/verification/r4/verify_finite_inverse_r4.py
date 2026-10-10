#!/usr/bin/env python3
"""Assemble disjoint Arb Jacobian batches and certify the finite block inverse.

This is a finite-dimensional result only; infinite tails are not checked.
"""
import argparse
import hashlib
import json
from flint import arb,arb_mat,ctx
from exact_finite_r4 import get_center,evaluate,ZERO,laur_q,conv,to_char,scale


def main(args):
    ctx.prec=args.bits;ctx.threads=1
    center,g,c,keys,ms=get_center(args.center)
    n=len(keys)+len(ms)
    sha=hashlib.sha256(open(args.center,'rb').read()).hexdigest()
    cells=[[arb(0) for _ in range(n)] for _ in range(n)]
    spans=[];exp_cache={}
    def power(e):
        if e not in exp_cache:exp_cache[e]=arb(2)**e
        return exp_cache[e]
    for path in args.batches:
        with open(path) as f:d=json.load(f)
        assert d['center_sha256']==sha and d['bits']==args.bits and d['dimension']==n
        spans.append((d['start'],d['end']))
        for i,j,m,e,r,f in d['entries']:
            assert d['start']<=j<d['end'] and 0<=i<n
            assert cells[i][j].is_zero()
            cells[i][j]=arb(m)*power(e)+arb(0,arb(r)*power(f))
    assert sorted(spans)==[(k,min(k+args.batch_width,n)) for k in range(0,n,args.batch_width)]
    M=arb_mat(cells)
    enclosed=M.inv()
    A=arb_mat([[enclosed[i,j].mid() for j in range(n)] for i in range(n)])
    if args.check_arb_inverse:
        with open(args.check_arb_inverse) as f:saved=json.load(f)
        assert saved['center_sha256']==sha and saved['bits']==args.bits and saved['dimension']==n
        assert all(tuple(int(x) for x in saved['entries'][i][j])==
                   tuple(int(x) for x in A[i,j].mid().man_exp())
                   for i in range(n) for j in range(n))
    if args.save_arb_inverse:
        payload={'center':args.center,'center_sha256':sha,'bits':args.bits,'dimension':n,
                 'entries':[[[int(q) for q in A[i,j].mid().man_exp()] for j in range(n)]
                            for i in range(n)]}
        with open(args.save_arb_inverse,'w') as f:json.dump(payload,f,separators=(',',':'))
    if args.save_midpoint_inverse:
        import numpy as np
        np.save(args.save_midpoint_inverse,
                np.array([[float(A[i,j].mid()) for j in range(n)] for i in range(n)]))
    E=arb_mat([[arb(int(i==j)) for j in range(n)] for i in range(n)])-A*M
    rho=arb(11)/10;b=c[1]
    win=[rho**m for m,s in keys]+[b**3*(j+2)*rho**(j-1) for j in ms]
    wout=[rho**m for m,s in keys]+[rho**(j-1) for j in ms]
    normA=max((sum((win[i]*abs(A[i,j]) for i in range(n)),ZERO)/wout[j]).upper()
              for j in range(n))
    boundary_edge_column=sum((win[i]*abs(A[i,n-1]) for i in range(n)),ZERO).upper()
    if not (boundary_edge_column<arb(54778)/1000):
        raise RuntimeError('finite edge-column bound failed')
    defect=max((sum((win[i]*abs(E[i,j]) for i in range(n)),ZERO)/win[j]).upper()
               for j in range(n))
    F,B=evaluate(g,c)
    ftail={k:v for k,v in F.items() if k not in keys}
    btail={j:v for j,v in B.items() if j not in ms}
    rhs=dict(btail)
    coupling=to_char(scale(conv(laur_q(g),laur_q(F)),2))
    for j,v in coupling.items():
        if j not in ms:rhs[j]=rhs.get(j,ZERO)-v
    maxj=max(rhs,default=ms[-1])
    run=arb(0);ctail={}
    for j in range(maxj if maxj%2 else maxj-1,ms[-1],-2):
        run+=rhs.get(j,ZERO)/(b**3*(j+2))
        ctail[j]=-run
    finite_values=[F.get(k,ZERO) for k in keys]+[B.get(j,ZERO) for j in ms]
    finite_values[-1]-=b**3*(ms[-1]+2)*ctail.get(ms[-1]+2,ZERO)
    r_f=arb_mat([[v] for v in finite_values])
    correction=A*r_f
    Yf=sum((win[i]*abs(correction[i,0]) for i in range(n)),ZERO).upper()
    Ytail=(sum((rho**m*abs(v) for (m,s),v in ftail.items()),ZERO)
           +sum((b**3*(j+2)*rho**(j-1)*abs(v) for j,v in ctail.items()),ZERO)).upper()
    assert normA<arb(1665434)/1000
    assert (Yf+Ytail).upper()<arb(3038)/10**9
    if not (defect<arb(1)):
        raise RuntimeError('finite inverse failed')
    assert defect<arb(21)/10**35, 'finite inverse defect ceiling failed'
    print('CERTIFIED FINITE BLOCK ONLY')
    print('dimension',n,'bits',args.bits,'batches',spans)
    print('normA upper',normA)
    print('last boundary column weighted norm upper',boundary_edge_column)
    print('defect upper',defect)
    print('preconditioned finite residual upper',Yf)
    print('preconditioned tail residual upper',Ytail)
    print('preconditioned total residual upper',(Yf+Ytail).upper())


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--center',default='finite_center_r4_M31_S24.json')
    p.add_argument('--bits',type=int,default=128)
    p.add_argument('--batch-width',type=int,default=104)
    p.add_argument('--save-midpoint-inverse',default='')
    p.add_argument('--save-arb-inverse',default='')
    p.add_argument('--check-arb-inverse',default='')
    p.add_argument('batches',nargs='+')
    main(p.parse_args())
