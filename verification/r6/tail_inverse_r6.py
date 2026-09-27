"""Rigorous column-norm bound for the explicit principal-tail inverse A."""
import argparse
import json
from flint import arb, arb_mat, ctx
from finite_interval_r6 import load_inverse
from verify_r6 import Inverse, make_algebra, max_upper, Z


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',default='tail_inverse_bound_r6.json')
    args=ap.parse_args()
    ctx.prec=128;ctx.threads=1
    d=json.load(open('center_r6.json'))
    rows,g,c,p,V,F=make_algebra(d)
    M,S=d['M'],d['S'];rho=arb(11)/10;b=c[1]
    Ahat=load_inverse('inverse_r6_M58_S30.npy',len(rows))
    inv=Inverse(M,S,rho,b,Ahat)
    vec=[Z for _ in rows]
    vec[inv.pos[M,0]]=inv.alpha*(M+1)
    vec[inv.pos[M,1]]=b**3*(M+3)/(M+4)
    vec[inv.pos[M,2]]=b**3/(M+4)
    v=Ahat*arb_mat([[x] for x in vec])
    H=sum((abs(v[i,0])*inv.cw[i] for i in range(inv.N)),Z)
    bounds=[]
    for n in range(M+4,M+404,4):
        eta=1/(inv.alpha*(n+1))
        shape=sum((inv.omega(k-1) for k in range(M+4,n+1,4)),Z)*eta
        gcoupling=sum((rho**k*b**3*((k+3)/(k+4)+1/(k+4)) for k in range(M+4,n,4)),Z)*eta
        fin=H*eta
        bounds.append((n,(shape+gcoupling+fin)/rho**n))
    # Beyond the checked interval: shape weights ratio <= tau, each lower g
    # coupling is <= (b^3/alpha)/(n+1)*rho^-4*tau = 2*tau*rho^-4/(n+1).
    # The finite correction decreases geometrically and by 1/(n+1).
    n0=M+404
    tau=1/(1-rho**-4)
    far=tau+2*tau*rho**-4/(n0+1)+H/(inv.alpha*(n0+1)*rho**n0)
    whole=max_upper([arb(1),far]+[v for _,v in bounds])
    result={'A_tail_upper':whole.str(30),'max_checked_n':max(bounds,key=lambda t:t[1].upper())[0],
            'max_checked_value':max_upper(v for _,v in bounds).str(30),
            'far_upper':far.upper().str(30),'H':H.str(30)}
    print(json.dumps(result,indent=2))
    json.dump(result,open(args.output,'w'),indent=2)


if __name__=='__main__':main()
