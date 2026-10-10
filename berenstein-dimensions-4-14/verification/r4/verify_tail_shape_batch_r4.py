#!/usr/bin/env python3
"""Arb enumeration of near-tail shape columns for the R4 divided system."""
import argparse
import hashlib
import json
from flint import arb,arb_mat,ctx
from exact_finite_r4 import (get_center,ZERO,laur_p,laur_q,laur_h,laur_a2,
                             laur_char,conv,to_char,scale,add,sine,unsine,
                             KD,poly)


def main(args):
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,keys,ms=get_center(args.center)
    M=d['M'];js=list(range(M+2,301,2));start,end=map(int,args.batch.split(':'))
    assert 0<=start<end<=len(js)
    with open(args.inverse) as f:inv=json.load(f)
    center_sha=hashlib.sha256(open(args.center,'rb').read()).hexdigest()
    assert inv['center_sha256']==center_sha and inv['bits']==args.bits
    N=len(keys)+len(ms);assert inv['dimension']==N
    powers={}
    def power(e):
        if e not in powers:powers[e]=arb(2)**e
        return powers[e]
    A=arb_mat([[arb(man)*power(exp) for man,exp in row] for row in inv['entries']])
    p=laur_p(c);q0=laur_q(g);h=laur_h(c);a2=laur_a2(p)
    V=sine(KD(g));b=c[1];rho=arb(11)/10
    win=[rho**m for m,s in keys]+[b**3*(j+2)*rho**(j-1) for j in ms]
    cols=[];tails=[]
    for j in js[start:end]:
        dp={j-1:arb(j)}
        fi=unsine(add(poly(poly(V,p,True),dp),poly(poly(V,dp,True),p)))
        dh=laur_char(j)
        da2=add(conv(dp,{-k:v for k,v in p.items()}),
                conv(p,{-k:v for k,v in dp.items()}))
        bi=to_char(add(scale(conv(conv(h,dh),a2),-2),
                       scale(conv(conv(h,h),da2),-1)))
        ftail={k:v for k,v in fi.items() if k not in keys}
        btail={k:v for k,v in bi.items() if k not in ms}
        coupling=to_char(scale(conv(q0,laur_q(fi)),2))
        rhs=dict(btail)
        for k,v in coupling.items():
            if k not in ms:rhs[k]=rhs.get(k,ZERO)-v
        maxj=max(rhs,default=ms[-1]);run=arb(0);ctail={}
        for k in range(maxj if maxj%2 else maxj-1,ms[-1],-2):
            run+=rhs.get(k,ZERO)/(b**3*(k+2));ctail[k]=-run
        first=ctail.get(M+2,ZERO)
        ctail[j]=ctail.get(j,ZERO)-1
        tailnorm=(sum((rho**m*abs(v) for (m,s),v in ftail.items()),ZERO)
                  +sum((b**3*(k+2)*rho**(k-1)*abs(v) for k,v in ctail.items()),ZERO))
        tails.append(tailnorm)
        finite=[fi.get(k,ZERO) for k in keys]+[bi.get(k,ZERO) for k in ms]
        finite[-1]-=b**3*(M+2)*first
        cols.append(finite)
    R=arb_mat([[col[i] for col in cols] for i in range(N)])
    X=A*R
    vals=[]
    for t,j in enumerate(js[start:end]):
        fnorm=sum((win[i]*abs(X[i,t]) for i in range(N)),ZERO)
        vals.append(((fnorm+tails[t])/(b**3*(j+2)*rho**(j-1))).upper())
    maxval=max(vals);argmax=start+vals.index(maxval)
    payload={'center':args.center,'center_sha256':center_sha,
             'inverse_sha256':hashlib.sha256(open(args.inverse,'rb').read()).hexdigest(),
             'bits':args.bits,'start':start,'end':end,'mode_count':len(js),
             'max_upper_dyadic':[int(x) for x in maxval.man_exp()],
             'argmax':argmax,'argmax_mode':js[argmax],
             'column_upper_dyadic':[[int(x) for x in v.man_exp()] for v in vals]}
    with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
    print('near shape batch',start,end,'max',maxval,'mode',js[argmax])


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--center',default='finite_center_r4_M31_S24.json')
    p.add_argument('--inverse',default='Ahat_arb_r4_M31_S24.json')
    p.add_argument('--bits',type=int,default=128)
    p.add_argument('--batch',required=True)
    p.add_argument('--output',required=True)
    main(p.parse_args())
