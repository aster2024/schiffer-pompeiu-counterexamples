"""Fail-closed Arb check of the frozen finite Jacobian and dyadic inverse.

This checks the finite matrix estimates of Section 6.1.
Tail and nonlinear estimates are checked by the other programs of Section 6.
"""
from proof_guard import require
import argparse
import hashlib
import json
import time
import numpy as np
from flint import arb,arb_mat,ctx
from exact_finite_rank2 import center
from exact_finite_r4 import dyadic,ZERO


def finite_upper(v):
    if not v.is_finite():
        raise ValueError('non-finite Arb enclosure')
    u=v.upper()
    if not u.is_finite():
        raise ValueError('non-finite Arb upper endpoint')
    return u


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--center',required=True)
    ap.add_argument('--inverse',required=True)
    ap.add_argument('--batch-width',type=int,default=50)
    ap.add_argument('--bits',type=int,default=128)
    ap.add_argument('--rho-num',type=int,default=102)
    ap.add_argument('--rho-den',type=int,default=100)
    ap.add_argument('--save-profile',default='')
    ap.add_argument('batches',nargs='+')
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    t=time.monotonic()
    d,g,c,B,keys,ms,js=center(a.center)
    n=len(keys)+len(js)
    sha=hashlib.sha256(open(a.center,'rb').read()).hexdigest()
    cells=[[ZERO for _ in range(n)] for _ in range(n)]
    seen=set();spans=[];pow2={}
    def power(e):
        if e not in pow2:pow2[e]=arb(2)**e
        return pow2[e]
    for path in a.batches:
        with open(path) as f:receipt=json.load(f)
        require((receipt['center_sha256']==sha and receipt['bits']==a.bits), 'verify_finite_block_rank2.py:48: proof gate failed')
        require((receipt['dimension']==n), 'verify_finite_block_rank2.py:49: proof gate failed')
        lo,hi=receipt['start'],receipt['end'];spans.append((lo,hi))
        for row,col,m,e,r,f in receipt['entries']:
            require((0<=row<n and lo<=col<hi and (row,col) not in seen), 'verify_finite_block_rank2.py:52: proof gate failed')
            seen.add((row,col))
            cells[row][col]=arb(m)*power(e)+arb(0,arb(r)*power(f))
    require((sorted(spans)==[(k,min(k+a.batch_width,n)) for k in range(0,n,a.batch_width)]), 'verify_finite_block_rank2.py:55: proof gate failed')
    print('assembled',n,'entries',len(seen),'seconds',time.monotonic()-t,flush=True)
    J=arb_mat(cells)
    aa=np.load(a.inverse)
    require((aa.shape==(n,n) and np.isfinite(aa).all()), 'verify_finite_block_rank2.py:59: proof gate failed')
    A=arb_mat([[dyadic(float(aa[i,j]).hex()) for j in range(n)] for i in range(n)])
    print('converted inverse seconds',time.monotonic()-t,flush=True)
    E=arb_mat([[arb(int(i==j)) for j in range(n)] for i in range(n)])-A*J
    print('multiplied seconds',time.monotonic()-t,flush=True)
    rho=arb(a.rho_num)/a.rho_den;m=d['m']
    win=[rho**nn for nn,s in keys]+[(j+2*m)*rho**(j-1) for j in js]
    wout=[rho**nn for nn,s in keys]+[rho**(nn-m) for nn in ms]
    maxA=arb(0);maxE=arb(0);maxAcol=-1;maxEcol=-1;profile=[]
    for j in range(n):
        ca=finite_upper(sum((win[i]*abs(A[i,j]) for i in range(n)),ZERO)/wout[j])
        man,exp=ca.man_exp();profile.append([int(man),int(exp)])
        ce=finite_upper(sum((win[i]*abs(E[i,j]) for i in range(n)),ZERO)/win[j])
        if ca>maxA:maxA=ca;maxAcol=j
        if ce>maxE:maxE=ce;maxEcol=j
    if (d['m'],d['M'],d['S'])==(6,114,36):
        if not maxA<arb(82859) or not maxE<arb(3)/10**9:
            raise ValueError('R14 finite inverse ceiling fails')
    elif not maxE<arb(1):
        raise ValueError('finite inverse defect is not contractive')
    print('CERTIFIED FINITE BLOCK ONLY')
    print('dimension',n,'rho',rho,'normA upper',maxA,'argmax',maxAcol,
          'defect upper',maxE,'argmax',maxEcol,
          'seconds',time.monotonic()-t)
    if a.save_profile:
        data=dict(center_sha256=sha,
                  inverse_sha256=hashlib.sha256(open(a.inverse,'rb').read()).hexdigest(),
                  bits=a.bits,rho_num=a.rho_num,rho_den=a.rho_den,
                  dimension=n,profile_upper_dyadic=profile)
        with open(a.save_profile,'w') as f:json.dump(data,f,separators=(',',':'))
        print('saved profile',a.save_profile)


if __name__=='__main__':
    main()
