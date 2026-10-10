"""Arb bounds for A_f times every high-shape column that can reach finite rows.

The reported operator is A_f P_f D_c F at the frozen center. The analytic
support inequalities make all later columns exactly zero. This is one block
bound, not a complete inverse or nonlinear existence certificate.
"""
from proof_guard import require
import argparse
import hashlib
import json
import numpy as np
from flint import arb,arb_mat,ctx
from exact_finite_rank2 import center,columns,h_of,dh_of,a2_of,to_char
from exact_finite_r4 import ZERO,add,scale,conv,dyadic


def finite_upper(x):
    require((x.is_finite()), 'cert_finite_cross_shape_rank2.py:17: proof gate failed')
    y=x.upper()
    require((y.is_finite()), 'cert_finite_cross_shape_rank2.py:19: proof gate failed')
    return y


def boundary_shape(c,B,m,j,h,a2,p):
    dp={j-1:arb(j)}
    dh=dh_of(c,j,m,B)
    da2=add(conv(dp,{-k:v for k,v in p.items()}),
            conv(p,{-k:v for k,v in dp.items()}))
    raw=add(scale(conv(conv(h,dh),a2),-2),
            scale(conv(conv(h,h),da2),-1))
    return scale(to_char(raw,m),1/c[1])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('inverse')
    ap.add_argument('--start-index',type=int,required=True)
    ap.add_argument('--end-index',type=int,required=True)
    ap.add_argument('--bits',type=int,default=160)
    ap.add_argument('--output',required=True)
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m,M=d['m'],d['M'];rho=arb(102)/100
    p={j-1:j*v for j,v in c.items()};pmax=max(p)
    h=h_of(c,m,B);a2=a2_of(c)
    hmax=m*max(c)
    h2max=2*hmax-m
    first=js[-1]+2*m
    jzero=first
    while not (jzero+m-1-hmax-pmax>M and
               jzero-1-pmax-h2max>M and
               jzero-1-M-pmax>M):
        jzero+=2*m
    total=(jzero-first)//(2*m)
    require((0<=a.start_index<a.end_index<=total), 'cert_finite_cross_shape_rank2.py:56: proof gate failed')
    aa=np.load(a.inverse,mmap_mode='r')
    n=len(keys)+len(js)
    require((aa.shape==(n,n) and np.isfinite(aa).all()), 'cert_finite_cross_shape_rank2.py:59: proof gate failed')
    A=arb_mat([[dyadic(float(aa[i,k]).hex()) for k in range(n)] for i in range(n)])
    win=[rho**nn for nn,s in keys]+[(j+2*m)*rho**(j-1) for j in js]
    vals=[]
    for idx in range(a.start_index,a.end_index):
        j=first+2*m*idx
        if j-1-M-pmax>M:
            Fi={}
            Gi=boundary_shape(c,B,m,j,h,a2,p)
        else:
            Fi,Gi=columns(g,c,B,m,keys,js+[j],len(keys)+len(js),len(keys)+len(js)+1)[0]
        vec=arb_mat([[Fi.get(k,ZERO)] for k in keys]+
                    [[Gi.get(nn,ZERO)] for nn in ms])
        x=A*vec
        value=sum((win[i]*abs(x[i,0]) for i in range(n)),ZERO)/((j+2*m)*rho**(j-1))
        value=finite_upper(value)
        man,exp=value.man_exp()
        vals.append([j,int(man),int(exp)])
        print('finite cross shape j',j,'upper',value,flush=True)
    data=dict(center_sha256=hashlib.sha256(open(a.center,'rb').read()).hexdigest(),
              inverse_sha256=hashlib.sha256(open(a.inverse,'rb').read()).hexdigest(),
              bits=a.bits,m=m,M=M,jzero=jzero,first=first,total=total,
              start_index=a.start_index,end_index=a.end_index,
              column_upper_dyadic=vals)
    with open(a.output,'w') as f:json.dump(data,f,separators=(',',':'))
    print('wrote',a.output,'zero for all j>=',jzero)


if __name__=='__main__':
    main()
