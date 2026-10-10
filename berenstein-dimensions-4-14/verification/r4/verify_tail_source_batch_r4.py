#!/usr/bin/env python3
"""Arb enumeration of near-tail source columns from Section 5.2.

The remaining columns are bounded by verify_far_source_bound_r4.py."""
import argparse
import hashlib
import json
from flint import arb,arb_mat,ctx
from exact_finite_r4 import (get_center,ZERO,laur_p,laur_q,laur_char,conv,to_char,
                             scale,add,sine,unsine,KD,abs2)


def near_modes(M,S):
    return ([(m,s) for m in range(1,M+1,2) for s in range(S+1,S+32)]
            +[(m,s) for m in range(M+2,2*M,2) for s in range(0,S+32)])


def main(args):
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,keys,ms=get_center(args.center)
    M,S=d['M'],d['S'];modes=near_modes(M,S);start,end=map(int,args.batch.split(':'))
    assert 0<=start<end<=len(modes)
    with open(args.inverse) as f:inv=json.load(f)
    center_sha=hashlib.sha256(open(args.center,'rb').read()).hexdigest()
    assert inv['center_sha256']==center_sha and inv['bits']==args.bits
    N=len(keys)+len(ms);assert inv['dimension']==N
    powers={}
    def power(e):
        if e not in powers:powers[e]=arb(2)**e
        return powers[e]
    A=arb_mat([[arb(man)*power(exp) for man,exp in row] for row in inv['entries']])
    p=laur_p(c);q0=laur_q(g);b=c[1];rho=arb(11)/10
    win=[rho**m for m,s in keys]+[b**3*(j+2)*rho**(j-1) for j in ms]
    cols=[];tails=[]
    for m,s in modes[start:end]:
        e={(m,s):arb(1)}
        fi=add(e,unsine(abs2(sine(KD(e)),p)))
        if s==0:
            dq=scale(laur_char(m),arb(1)/(2*(m+1)))
            bi=to_char(scale(conv(q0,dq),2))
        else:bi={}
        ftail={k:v for k,v in fi.items() if k not in keys}
        btail={j:v for j,v in bi.items() if j not in ms}
        coupling=to_char(scale(conv(q0,laur_q(fi)),2))
        rhs=dict(btail)
        for j,v in coupling.items():
            if j not in ms:rhs[j]=rhs.get(j,ZERO)-v
        maxj=max(rhs,default=ms[-1]);run=arb(0);ctail={}
        for j in range(maxj if maxj%2 else maxj-1,ms[-1],-2):
            run+=rhs.get(j,ZERO)/(b**3*(j+2));ctail[j]=-run
        first=ctail.get(M+2,ZERO)
        ftail[(m,s)]=ftail.get((m,s),ZERO)-1
        tailnorm=(sum((rho**mm*abs(v) for (mm,ss),v in ftail.items()),ZERO)
                  +sum((b**3*(j+2)*rho**(j-1)*abs(v) for j,v in ctail.items()),ZERO))
        tails.append(tailnorm)
        finite=[fi.get(k,ZERO) for k in keys]+[bi.get(j,ZERO) for j in ms]
        finite[-1]-=b**3*(M+2)*first
        cols.append(finite)
    R=arb_mat([[col[i] for col in cols] for i in range(N)])
    X=A*R
    vals=[]
    for t,(m,s) in enumerate(modes[start:end]):
        fnorm=sum((win[i]*abs(X[i,t]) for i in range(N)),ZERO)
        vals.append(((fnorm+tails[t])/rho**m).upper())
    maxval=max(vals);argmax=start+vals.index(maxval)
    payload={'center':args.center,'center_sha256':center_sha,
             'inverse_sha256':hashlib.sha256(open(args.inverse,'rb').read()).hexdigest(),
             'bits':args.bits,'start':start,'end':end,'mode_count':len(modes),
             'max_upper_dyadic':[int(x) for x in maxval.man_exp()],
             'argmax':argmax,'argmax_mode':modes[argmax],
             'column_upper_dyadic':[[int(x) for x in v.man_exp()] for v in vals]}
    with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
    print('near source batch',start,end,'max',maxval,'mode',modes[argmax])


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--center',default='finite_center_r4_M31_S24.json')
    p.add_argument('--inverse',default='Ahat_arb_r4_M31_S24.json')
    p.add_argument('--bits',type=int,default=128)
    p.add_argument('--batch',required=True)
    p.add_argument('--output',required=True)
    main(p.parse_args())
