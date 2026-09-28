"""Rigorous norm of the explicit R8 principal-tail inverse.

The finite list is checked column by column in Arb. The remaining tail uses
positivity and sum-one of the r^(n+6) Zernike coefficients.
"""

import json
import numpy as np
from flint import arb, arb_mat, ctx
import verify_r8 as v


def load_inverse(path):
    raw=np.load(path)
    return arb_mat([[v.dyadic(float(raw[i,j]).hex())
                     for j in range(raw.shape[1])] for i in range(raw.shape[0])])


def main():
    ctx.prec=128
    ctx.threads=1
    d=json.load(open('center_r8_M45_S24.json'))
    rows,g,c,p,V,F=v.make_algebra(d)
    M,S=d['M'],d['S']
    rho=arb(11)/10
    A=load_inverse('inverse_r8_M45_S24.npy')
    inv=v.Inverse(M,S,rho,c[1],A)
    q=inv.radial_coeff(M)
    vec=[v.Z for _ in rows]
    vec[inv.pos[M,0]]=inv.alpha*(M+1)
    for s in range(1,v.M_ROOT+1):
        vec[inv.pos[M,s]]=inv.alpha*(M+v.M_ROOT+1)*q[s]
    image=A*arb_mat([[x] for x in vec])
    H=sum((abs(image[i,0])*inv.cw[i] for i in range(inv.N)),v.Z)
    n0=M+v.STEP*101
    ns=list(range(M+v.STEP,n0,v.STEP))
    vals=inv.apply_many([{(n,0):arb(1)} for n in ns])
    ratios=[val/rho**n for n,val in zip(ns,vals)]
    tau=1/(1-rho**-v.STEP)
    # For an input at n>=n0: the shape response is at most tau; the total
    # nonharmonic correction at each lower k is alpha*m/(alpha*(n+1)), since
    # sum_{s>=1}(k+m+1)q_s(k)=m. The first-tail finite correction is H/alpha(n+1).
    far=tau+v.M_ROOT*tau*rho**-v.STEP/(n0+1)+H/(inv.alpha*(n0+1)*rho**n0)
    bound=v.max_upper([arb(1),far]+ratios)
    imax=v.argmax_upper(ratios)
    out={'status':'OPEN_COMPONENT','M':M,'S':S,'count':len(ns),
         'near_max':v.max_upper(ratios).str(30),'near_arg':ns[imax],
         'far_start':n0,'far_upper':far.upper().str(30),
         'A_tail_upper':bound.str(30),'H':H.str(30)}
    with open('tail_inverse_bound_r8.json','w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
