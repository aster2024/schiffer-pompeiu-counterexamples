"""Two-norm analytic far-shape bound with tiny low-frequency spill for R14.

Output above M uses the principal-tail inverse norm. The remaining exact
low-frequency coefficients use the full inverse norm. This permits a far
cutoff earlier than the crude full-support threshold.
"""

import argparse
import json
from pathlib import Path
from flint import arb, ctx
import rank2_cap_core as v
from far_shape_rank2 import q_of_g, check_q


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--jstar',type=int,required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args()
    ctx.prec=128
    ctx.threads=1
    d=json.loads(Path('center_r14_M114_S60.json').read_text())
    rows,g,c,pcoef,V,F=v.make_algebra(d)
    M=d['M']
    rho=arb(11)/10
    b=c[1]
    alpha=b**(v.M_ROOT+1)/v.BOUNDARY_SCALE
    v.need(args.jstar%v.STEP==1,'invalid far shape index')
    h=v.hcoeff(c)
    h0={(v.M_ROOT,0):b**v.M_ROOT/v.BOUNDARY_SCALE}
    G1=v.add(v.hol(v.sine(h),pcoef,True),
             v.scale(v.sine(h0),-b))
    a=v.abs2({(0,0):arb(1)},pcoef)
    G0=a
    for _ in range(v.M_ROOT-1):
        G0=v.hol_psi(G0,c)
    G0=v.add(v.scale(G0,1/v.BOUNDARY_SCALE),
             {(v.M_ROOT-1,0):-alpha})
    Q=q_of_g(g,M,d['S'])
    q_error=check_q(g,Q)
    L=M+max(pcoef)
    v.need(args.jstar>L,'far shape start must exceed Q shift')
    T=v.shift(v.hol(v.sine(Q),pcoef,True),L)
    v.need(min(k for (k,t) in T)>=0,'Q shift support failure')
    An=arb(json.loads(Path('finite_bound_r14_M114_S60.json').read_text())['A_finite']).upper()
    Atail=arb(json.loads(Path('tail_inverse_bound_r14_M114_S60.json').read_text())['A_tail_upper']).upper()
    def split_norm(f,offset):
        good=bad=arb(0)
        for (k,t),x in f.items():
            mass=abs(x)*rho**abs(k)
            if k+offset>M:good+=mass
            else:bad+=mass
        return good,bad
    N1g,N1b=split_norm(G1,args.jstar-1)
    N0g,N0b=split_norm(G0,args.jstar)
    wgood=wbad=arb(0)
    for (k,t),x in T.items():
        # Keep each exact output frequency in the two-boundary-factor estimate.
        # For j>=jstar, k+j-L grows, so its value at jstar bounds every
        # higher column. Replacing all denominators by (jstar-L)^2 is much
        # looser for the large k carrying most of the Q mass.
        mass=(abs(x)*(2*t+1)*(2*t+3)*rho**(k-L-v.M_ROOT)
              /(k+args.jstar-L)**2)
        if k+args.jstar-1-L>M:wgood+=mass
        else:wbad+=mass
    nonW=(2*(Atail*N1g+An*N1b)/rho**v.M_ROOT
          +v.M_ROOT*(Atail*N0g+An*N0b)
            /(rho**(v.M_ROOT-1)*(args.jstar+v.M_ROOT)))/alpha
    Wbound=8*(Atail*wgood+An*wbad)/alpha
    total=nonW+Wbound
    out={'status':'OPEN_COMPONENT','dimension':14,'jstar':args.jstar,
         'M':M,'L':L,'Q_identity_error':q_error.str(25),
         'Gamma1_good':N1g.str(25),'Gamma1_bad':N1b.str(25),
         'Gamma0_good':N0g.str(25),'Gamma0_bad':N0b.str(25),
         'Wsum_good':wgood.str(25),'Wsum_bad':wbad.str(25),
         'nonW_ratio':nonW.str(25),'W_ratio':Wbound.str(25),
         'Z_shape_far':total.upper().str(25)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
