"""Arb coarse finite-input defect estimates of Section 6.4.

The resulting column bounds are combined with the last-column residual estimate."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,columns,q_of,to_char
from exact_finite_r4 import ZERO,conv,scale


def top(v):
    require((v.is_finite()), 'cert_finite_input_coarse_rank2.py:15: proof gate failed')
    x=v.upper();require((x.is_finite()), 'cert_finite_input_coarse_rank2.py:16: proof gate failed')
    return x


def pair(v):
    x=top(v);a,e=x.man_exp();return [int(a),int(e)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--start',type=int,required=True)
    ap.add_argument('--end',type=int,required=True)
    ap.add_argument('--output',required=True)
    a=ap.parse_args()
    ctx.prec=160;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m=d['m'];require(((m,d['M'],d['S'])==(6,114,36)), 'cert_finite_input_coarse_rank2.py:33: proof gate failed')
    require((0<=a.start<a.end<=len(keys)+len(js)), 'cert_finite_input_coarse_rank2.py:34: proof gate failed')
    rho=arb(102)/100;b=c[1]
    hb=b**m/B;alpha=hb*hb
    q0=q_of(g,m)
    k=arb(474)/1000
    E=arb(3)/10**9
    Cfactor=1/(1-k)
    Xs=arb(1)/5;Xc=arb(10);H=arb(250)
    records=[];needs=[];best=arb(0);besti=-1
    for i,(F,G) in enumerate(columns(g,c,B,m,keys,js,a.start,a.end),a.start):
        w=(rho**keys[i][0] if i<len(keys) else
           (js[i-len(keys)]+2*m)*rho**(js[i-len(keys)]-1))
        ft=sum((rho**n*abs(v) for (n,s),v in F.items() if (n,s) not in keys),ZERO)/w
        Q=to_char(scale(conv(q0,q_of(F,m)),2/b),m)
        rhs={n:G.get(n,ZERO)-Q.get(n,ZERO) for n in set(G)|set(Q)}
        ct={};run=arb(0)
        if rhs:
            for j in range(max(rhs)-m+1,js[-1],-2*m):
                run+=rhs.get(j+m-1,ZERO)/(alpha*(j+2*m))
                ct[j]=-run
        c0=sum(((j+2*m)*rho**(j-1)*abs(v) for j,v in ct.items()),ZERO)/w
        z=E+(1+Xs)*ft+((1+Xs)*H+(1+Xc))*Cfactor*c0
        z=top(z)
        if z>=arb(1):needs.append(i)
        if z>best:best=z;besti=i
        records.append(dict(i=i,ft=pair(ft),c0=pair(c0),coarseZ=pair(z)))
    data=dict(center_sha256=hashlib.sha256(open(a.center,'rb').read()).hexdigest(),
              start=a.start,end=a.end,dimension=len(keys)+len(js),
              records=records,needs_refinement=needs)
    with open(a.output,'w') as f:json.dump(data,f,separators=(',',':'))
    print('FINITE_INPUT_COARSE_ONLY','range',a.start,a.end,
          'max',best,'argmax',besti,'needs',needs,'wrote',a.output)


if __name__=='__main__':
    main()
