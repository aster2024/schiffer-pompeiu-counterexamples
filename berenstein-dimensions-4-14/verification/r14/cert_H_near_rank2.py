"""Arb source norm of near high-shape derivative columns."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center,columns
from exact_finite_r4 import ZERO


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--J',type=int,required=True)
    ap.add_argument('--start-index',type=int,required=True)
    ap.add_argument('--end-index',type=int,required=True)
    ap.add_argument('--output',required=True)
    a=ap.parse_args()
    ctx.prec=160;ctx.threads=1
    d,g,c,B,keys,ms,js=center(a.center)
    m=d['m'];rho=arb(102)/100;first=js[-1]+2*m
    require(((a.J-first)%(2*m)==0), 'cert_H_near_rank2.py:21: proof gate failed')
    total=(a.J-first)//(2*m)
    require((0<=a.start_index<a.end_index<=total), 'cert_H_near_rank2.py:23: proof gate failed')
    vals=[]
    for idx in range(a.start_index,a.end_index):
        j=first+2*m*idx
        F,G=columns(g,c,B,m,keys,js+[j],len(keys)+len(js),len(keys)+len(js)+1)[0]
        val=sum((rho**n*abs(v) for (n,s),v in F.items()),ZERO)/((j+2*m)*rho**(j-1))
        require((val.is_finite()), 'cert_H_near_rank2.py:29: proof gate failed')
        top=val.upper();require((top.is_finite() and top<arb(250)), 'cert_H_near_rank2.py:30: proof gate failed')
        man,exp=top.man_exp();vals.append([j,int(man),int(exp)])
        print('H norm j',j,'upper',top,flush=True)
    data=dict(center_sha256=hashlib.sha256(open(a.center,'rb').read()).hexdigest(),
              m=m,J=a.J,first=first,total=total,
              start_index=a.start_index,end_index=a.end_index,
              column_upper_dyadic=vals)
    with open(a.output,'w') as f:json.dump(data,f,separators=(',',':'))
    print('wrote',a.output)


if __name__=='__main__':
    main()
