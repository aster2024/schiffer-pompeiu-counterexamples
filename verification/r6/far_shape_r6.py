"""Conservative analytic far-shape bound using K(g)=(1-r^2)^2 Q.

This avoids the cancellation-heavy R4 far-window decomposition. It is valid
only once every output frequency of the tail columns lies above the finite
cutoff M; that support condition is checked explicitly.
"""
import argparse
import json
from flint import arb, ctx
from verify_r6 import Z, acc, add, hcoeff, K, hol, hol_psi, abs2, make_algebra, norm, scale, shift, sine


def q_of_g(g,M,S):
    Q={}
    for n in range(2,M+1,4):
        C1=[{0:arb(1)}];C2=[{0:arb(1)}]
        for k in range(1,S):
            one={i:v*(k+n)/(k+n+1) for i,v in C1[-1].items()}
            one[k]=arb(2*k+n+1)/(k+n+1)
            two={i:v*(2*k+n+2)/(k+n+2) for i,v in one.items()}
            for i,v in C2[-1].items():acc(two,i,v*(k+n)/(k+n+2))
            C1.append(one);C2.append(two)
        for s in range(1,S+1):
            for t,v in C2[s-1].items():
                acc(Q,(n,t),g.get((n,s),Z)*v/(4*s*(s+1)))
    return Q


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--jstar',type=int,default=233)
    ap.add_argument('--output',default='far_shape_bound_r6.json')
    args=ap.parse_args()
    ctx.prec=128;ctx.threads=1
    d=json.load(open('center_r6.json'))
    rows,g,c,p,V,F=make_algebra(d)
    M,S=d['M'],d['S'];rho=arb(11)/10;b=c[1]
    assert args.jstar%4==1
    h=hcoeff(c);h0={(2,0):b*b/2};alpha=b**3/2
    Gamma1=add(hol(sine(h),p,True),scale(sine(h0),-b))
    Gamma0=add(hol_psi(abs2({(0,0):arb(1)},p),c),{(1,0):-b**3})
    N1=norm(Gamma1,rho);N0=norm(Gamma0,rho)
    Q=q_of_g(g,M,S)
    # Q includes only sin(2 mod 4) frequencies; bar(p) can lower m by max(p).
    L=M+max(p)
    T=shift(hol(sine(Q),p,True),L)
    assert min(m for (m,t) in T)>=0
    # The non-W polynomial has harmonic frequencies at least j-1-minshift.
    # A conservative absolute bound: psi^2 terms max 2*max(c), p max max(p).
    min_freq=args.jstar-1-(2*max(c)+max(p))
    assert min_freq>M, (min_freq,M)
    nonW=(2*N1/rho**2+N0/(rho*(args.jstar+2)))/alpha
    # Bound 2 Re(j w^{j-1} bar(p) W), where W=(1-r^2)^2 Q.
    wsum=sum((abs(v)*(2*t+1)*(2*t+3)*rho**(m-L-2) for (m,t),v in T.items()),Z)
    Wbound=8*wsum/(alpha*(args.jstar-L)**2)
    Atail=arb('3.164016')
    total=Atail*(nonW+Wbound)
    out={'jstar':args.jstar,'M':M,'L':L,'min_freq':min_freq,
         'Gamma1':N1.str(25),'Gamma0':N0.str(25),'Qnorm':norm(Q,rho).str(25),
         'Wsum':wsum.str(25),'nonW_ratio':nonW.str(25),
         'W_ratio':Wbound.str(25),'Zshape_far':total.upper().str(25)}
    print(json.dumps(out,indent=2))
    json.dump(out,open(args.output,'w'),indent=2)


if __name__=='__main__':main()
