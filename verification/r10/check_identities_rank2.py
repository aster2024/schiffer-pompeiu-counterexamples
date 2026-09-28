"""Exact K identities and high-precision principal/Q checks for R10/R14."""

import argparse
import json
from pathlib import Path
import sympy as sp
from flint import arb, ctx
import rank2_cap_core as v
from far_shape_rank2 import q_of_g, check_q


def check_k(m,step):
    r=sp.symbols('r',positive=True)
    count=0
    for n in (m,m+step,m+2*step):
        for s in range(1,5):
            S=lambda t:r**n*sp.jacobi(t,0,n,2*r*r-1)
            d=n+2*s
            k=S(s-1)/(4*d*(d+1))-S(s)/(2*d*(d+2))+S(s+1)/(4*(d+1)*(d+2))
            lap=sp.diff(k,r,2)+sp.diff(k,r)/r-n*n*k/r**2
            v.need(sp.cancel(lap-S(s))==0,'K Laplacian identity failed')
            v.need(sp.simplify(k.subs(r,1))==0,'K Dirichlet trace failed')
            v.need(sp.simplify(sp.diff(k,r).subs(r,1))==0,'K Neumann trace failed')
            count+=1
    return count


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    d=json.loads(Path(args.centre).read_text())
    m,step=d['m'],d['step']
    k_cases=check_k(m,step)
    ctx.prec=256
    ctx.threads=1
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    max_principal=arb(0)
    for j in (1+step,1+2*step,d['M']+m+1):
        cb={1:arb(1)}
        pb={0:arb(1)}
        Vb=v.hcoeff(cb)
        actual=v.shape_column(j,cb,pb,Vb)
        expected=v.principal_shape(j,arb(1))
        err=v.max_upper(abs(actual.get(key,v.Z)-expected.get(key,v.Z))
                        for key in actual.keys()|expected.keys())
        max_principal=v.max_upper([max_principal,err])
    v.need(max_principal<arb('1e-50'),'ball-principal identity failed')
    q_cases=0
    inv=v.Inverse(d['M'],d['S'],arb(11)/10,arb(1),None)
    for n in (m,m+step,d['M']):
        q=inv.radial_coeff(n)
        v.need(all(x>0 for x in q.values()),'radial-power coefficient not positive')
        v.need(abs(sum(q.values(),v.Z)-1).upper()<arb('1e-50'),'radial-power coefficients do not sum to one')
        q_cases+=1
    Q=q_of_g(g,d['M'],d['S'])
    q_error=check_q(g,Q)
    out={'status':'IDENTITIES_PASS','dimension':2*m+2,'bits':ctx.prec,
         'K_cases':k_cases,'radial_power_cases':q_cases,
         'ball_principal_max_error':max_principal.str(30),
         'Q_connection_max_error':q_error.str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
