"""Analytic far-shape bound for the R8 D6-symmetric CAP.

Uses K(g)=(1-r^2)^2 Q and a checked Zernike connection identity. Valid only
when every output frequency is beyond the finite cutoff.
"""

import argparse
import json
from flint import arb, ctx
import verify_r8 as v


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
    p.add_argument('--jstar',type=int,default=223)
    p.add_argument('--output',default='far_shape_bound_r8.json')
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.load(open('center_r8_M45_S24.json'))
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    b=c[1]
    alpha=b**(v.M_ROOT+1)
    v.need(args.jstar%v.STEP==1,'invalid far shape index')
    h=v.hcoeff(c)
    h0={(v.M_ROOT,0):b**v.M_ROOT}
    Gamma1=v.add(v.hol(v.sine(h),pcoef,True),
                 v.scale(v.sine(h0),-b))
    a=v.abs2({(0,0):arb(1)},pcoef)
    Gamma0=v.hol_psi(v.hol_psi(a,c),c)
    Gamma0=v.add(Gamma0,{(v.M_ROOT-1,0):-b**(v.M_ROOT+1)})
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
    Atail=arb('2.303')
    total=Atail*(nonW+Wbound)
    out={'status':'OPEN_COMPONENT','jstar':args.jstar,'M':M,'L':L,
         'min_freq':min_freq,'Gamma1':N1.str(25),'Gamma0':N0.str(25),
         'Qnorm':v.norm(Q,rho).str(25),'Q_identity_error':q_error.str(25),
         'Wsum':wsum.str(25),'nonW_ratio':nonW.str(25),
         'W_ratio':Wbound.str(25),'Z_shape_far':total.upper().str(25)}
    print(json.dumps(out,indent=2))
    with open(args.output,'w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')


if __name__=='__main__':
    main()
