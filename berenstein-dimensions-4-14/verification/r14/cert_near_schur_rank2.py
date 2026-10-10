"""Arb column bounds for a finite, exactly covered high-shape Schur prefix.

This certifies only the specified columns, not a full Berenstein solution.
"""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,columns,q_of,to_char
from exact_finite_r4 import ZERO,conv,scale


def finite_upper(v):
    require((v.is_finite()), 'cert_near_schur_rank2.py:14: proof gate failed')
    y=v.upper()
    require((y.is_finite()), 'cert_near_schur_rank2.py:16: proof gate failed')
    return y


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--J',type=int,required=True)
    ap.add_argument('--start-index',type=int,required=True)
    ap.add_argument('--end-index',type=int,required=True)
    ap.add_argument('--rho-num',type=int,default=11)
    ap.add_argument('--rho-den',type=int,default=10)
    ap.add_argument('--bits',type=int,default=192)
    ap.add_argument('--output',required=True)
    a=ap.parse_args()
    ctx.prec=a.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m=d['m'];b=c[1]
    first=js[-1]+2*m
    require(((a.J-first)%(2*m)==0 and a.J>first), 'cert_near_schur_rank2.py:35: proof gate failed')
    total=(a.J-first)//(2*m)
    require((0<=a.start_index<a.end_index<=total), 'cert_near_schur_rank2.py:37: proof gate failed')
    rho=arb(a.rho_num)/a.rho_den
    hb=(arb(1)/2 if m==2 else arb(1))*b**m/B
    alpha=hb*hb
    q0=q_of(g,m)
    vals=[]
    for idx in range(a.start_index,a.end_index):
        j=first+2*m*idx
        F,G=columns(g,c,B,m,keys,js+[j],len(keys)+len(js),len(keys)+len(js)+1)[0]
        coupling=to_char(scale(conv(q0,q_of(F,m)),2/b),m)
        rhs={n:G.get(n,ZERO)-coupling.get(n,ZERO)
             for n in set(G)|set(coupling)}
        ct={};run=arb(0)
        if rhs:
            for jj in range(max(rhs)-m+1,js[-1],-2*m):
                run+=rhs.get(jj+m-1,ZERO)/(alpha*(jj+2*m))
                ct[jj]=-run
        norm=sum((((jj+2*m)*rho**(jj-1))*abs(ct.get(jj,ZERO)-(arb(1) if jj==j else ZERO))
                  for jj in set(ct)|{j}),ZERO)/((j+2*m)*rho**(j-1))
        val=finite_upper(norm)
        man,exp=val.man_exp()
        vals.append([j,int(man),int(exp)])
        print('near Schur j',j,'upper',val,flush=True)
    payload=dict(center=a.center,
                 center_sha256=hashlib.sha256(open(a.center,'rb').read()).hexdigest(),
                 J=a.J,rho_num=a.rho_num,rho_den=a.rho_den,bits=a.bits,
                 start_index=a.start_index,end_index=a.end_index,total=total,
                 columns=vals)
    with open(a.output,'w') as f:json.dump(payload,f,separators=(',',':'))
    print('wrote',a.output)


if __name__=='__main__':
    main()
