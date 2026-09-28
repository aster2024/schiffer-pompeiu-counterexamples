"""Analytic far-shape bound for scaled R10/R14 fixed-disc CAPs."""

import argparse
import json
from pathlib import Path
from flint import arb, ctx
import rank2_cap_core as v


def q_of_g(g,M,S):
    Q={}
    for n in range(v.M_ROOT,M+1,v.STEP):
        C1=[{0:arb(1)}]
        C2=[{0:arb(1)}]
        for k in range(1,S):
            one={i:x*(k+n)/(k+n+1) for i,x in C1[-1].items()}
            one[k]=arb(2*k+n+1)/(k+n+1)
            two={i:x*(2*k+n+2)/(k+n+2) for i,x in one.items()}
            for i,x in C2[-1].items():
                v.acc(two,i,x*(k+n)/(k+n+2))
            C1.append(one)
            C2.append(two)
        for s in range(1,S+1):
            for t,x in C2[s-1].items():
                v.acc(Q,(n,t),g.get((n,s),v.Z)*x/(4*s*(s+1)))
    return Q


def check_q(g,Q):
    r2=v.shift(v.shift(v.sine(Q),1),1,True)
    r4=v.shift(v.shift(r2,1),1,True)
    recovered=v.unsine(v.add(v.sine(Q),v.scale(r2,-arb(2)),r4))
    target=v.K(g)
    diff=v.add(recovered,v.scale(target,-arb(1)))
    error=v.max_upper(abs(x) for x in diff.values())
    v.need(error<arb('1e-20'),'Q connection identity failed')
    return error


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--centre',required=True)
    p.add_argument('--tail-receipt',required=True)
    p.add_argument('--jstar',type=int,required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path(args.centre).read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    b=c[1]
    alpha=b**(v.M_ROOT+1)/v.BOUNDARY_SCALE
    v.need(args.jstar%v.STEP==1,'invalid far shape index')
    h=v.hcoeff(c)
    h0={(v.M_ROOT,0):b**v.M_ROOT/v.BOUNDARY_SCALE}
    Gamma1=v.add(v.hol(v.sine(h),pcoef,True),
                 v.scale(v.sine(h0),-b))
    a=v.abs2({(0,0):arb(1)},pcoef)
    Gamma0=a
    for _ in range(v.M_ROOT-1):
        Gamma0=v.hol_psi(Gamma0,c)
    Gamma0=v.add(v.scale(Gamma0,1/v.BOUNDARY_SCALE),
                 {(v.M_ROOT-1,0):-alpha})
    N1=v.norm(Gamma1,rho)
    N0=v.norm(Gamma0,rho)
    Q=q_of_g(g,M,S)
    q_error=check_q(g,Q)
    L=M+max(pcoef)
    T=v.shift(v.hol(v.sine(Q),pcoef,True),L)
    v.need(min(k for (k,t) in T)>=0,'Q shift support failure')
    min_freq=args.jstar-1-(v.M_ROOT*max(c)+max(pcoef))
    v.need(min_freq>M,'far shape output enters finite rows')
    nonW=(2*N1/rho**v.M_ROOT
          +v.M_ROOT*N0/(rho**(v.M_ROOT-1)*(args.jstar+v.M_ROOT)))/alpha
    wsum=sum((abs(x)*(2*t+1)*(2*t+3)*rho**(k-L-v.M_ROOT)
              for (k,t),x in T.items()),v.Z)
    Wbound=8*wsum/(alpha*(args.jstar-L)**2)
    tail=json.loads(Path(args.tail_receipt).read_text())
    Atail=arb(tail['A_tail_upper']).upper()
    total=Atail*(nonW+Wbound)
    out={'status':'OPEN_COMPONENT','dimension':2*v.M_ROOT+2,
         'jstar':args.jstar,'M':M,'L':L,'min_freq':min_freq,
         'Gamma1':N1.str(25),'Gamma0':N0.str(25),
         'Qnorm':v.norm(Q,rho).str(25),'Q_identity_error':q_error.str(25),
         'Wsum':wsum.str(25),'nonW_ratio':nonW.str(25),
         'W_ratio':Wbound.str(25),'Z_shape_far':total.upper().str(25)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
