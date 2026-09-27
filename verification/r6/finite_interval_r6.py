"""Batched finite columns for the fail-closed R6 interval workbench.

Each command is bounded so the local DSW is not held by a long dense job.
Output is a batch of outward Arb upper bounds; a global certificate must take
the maximum over *all* batches and establish the infinite tail estimates.
"""
import argparse
import json
import numpy as np
from flint import arb, arb_mat, ctx
from verify_r6 import (Z, Inverse, abs2, add, hol, hol_psi, K, make_algebra,
                       max_upper, norm, scale, shift, sine, unsine)


def load_inverse(path, n):
    a=np.load(path)
    assert a.shape==(n,n) and np.isfinite(a).all()
    mat=[]
    for row in a:
        mat.append([arb(int(num))/int(den) for num,den in (float(x).as_integer_ratio() for x in row)])
    return arb_mat(mat)


def make_column(i,rows,g,c,p,V,ng):
    if i<ng:
        e={rows[i]:arb(1)}
        return add(e,unsine(abs2(sine(K(e)),p)))
    js=sorted(c)
    j=js[i-ng]
    A=shift(hol(sine(V),p,True),j-1)
    B=shift(hol_psi(abs2({(1,0):arb(1)},p),c),j-1)
    return scale(unsine(add(scale(A,arb(j)),scale(B,-arb(1)/2))),arb(2))


def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--start',type=int,default=0)
    pa.add_argument('--stop',type=int,default=0)
    pa.add_argument('--inverse',default='inverse_r6_M58_S30.npy')
    pa.add_argument('--centre',default='center_r6.json')
    pa.add_argument('--bits',type=int,default=128)
    pa.add_argument('--output')
    ar=pa.parse_args()
    ctx.prec=ar.bits;ctx.threads=1
    data=json.load(open(ar.centre));rows,g,c,p,V,F=make_algebra(data)
    N=len(rows);ng=len(g);rho=arb(11)/10
    Ahat=load_inverse(ar.inverse,N)
    inv=Inverse(data['M'],data['S'],rho,c[1],Ahat)
    if ar.stop==0:
        val=inv.apply_many([F])[0].upper()
        output={'kind':'residual','Y':val.str(30),'N':N,'F_norm':norm(F,rho).str(30)}
    else:
        assert 0<=ar.start<ar.stop<=N
        cols=[make_column(i,rows,g,c,p,V,ng) for i in range(ar.start,ar.stop)]
        raw=arb_mat([[col.get(key,Z) for col in cols] for key in rows])
        defect=Ahat*raw
        D=[]
        for q,i in enumerate(range(ar.start,ar.stop)):
            D.append(sum((abs(defect[k,q]-(1 if k==i else 0))*inv.cw[k] for k in range(N)),Z)/inv.cw[i])
        val=inv.apply_many(cols,[('finite',i) for i in range(ar.start,ar.stop)])
        output={'kind':'columns','start':ar.start,'stop':ar.stop,'inverse_defect_max':max_upper(D).str(30),
                'Zfin_max':max_upper(v/inv.cw[i] for i,v in zip(range(ar.start,ar.stop),val)).str(30),
                'N':N}
    print(json.dumps(output,indent=2),flush=True)
    if ar.output:
        with open(ar.output,'w') as f:json.dump(output,f,indent=2)


if __name__=='__main__':main()
