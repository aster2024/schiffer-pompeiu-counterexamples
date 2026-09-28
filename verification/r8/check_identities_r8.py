"""Independent exact and high-precision identity checks for R8 CAP formulas."""

import json
import sympy as sp
from flint import arb, ctx
import verify_r8 as v
from far_shape_r8 import q_of_g, check_q


def check_k():
    r=sp.symbols('r',positive=True)
    count=0
    for n in (3,9,15):
        for s in range(1,5):
            S=lambda t: r**n*sp.jacobi(t,0,n,2*r*r-1)
            d=n+2*s
            k=(S(s-1)/(4*d*(d+1))-S(s)/(2*d*(d+2))
               +S(s+1)/(4*(d+1)*(d+2)))
            lap=sp.diff(k,r,2)+sp.diff(k,r)/r-n*n*k/r**2
            v.need(sp.cancel(lap-S(s))==0,'K Laplacian identity failed')
            v.need(sp.simplify(k.subs(r,1))==0,'K Dirichlet trace failed')
            v.need(sp.simplify(sp.diff(k,r).subs(r,1))==0,'K Neumann trace failed')
            count+=1
    return count


def main():
    k_cases=check_k()
    ctx.prec=256
    ctx.threads=1
    max_principal=arb(0)
    for j in (7,13,49,121):
        c={1:arb(1)}
        p={0:arb(1)}
        V=v.hcoeff(c)
        actual=v.shape_column(j,c,p,V)
        expected=v.principal_shape(j,arb(1))
        err=v.max_upper(abs(actual.get(key,v.Z)-expected.get(key,v.Z))
                        for key in actual.keys()|expected.keys())
        max_principal=v.max_upper([max_principal,err])
    v.need(max_principal<arb('1e-60'),'ball-principal identity failed')
    q_cases=0
    for n in (3,9,45):
        q=v.Inverse(45,24,arb(11)/10,arb(1),None).radial_coeff(n)
        v.need(all(x>0 for x in q.values()),'radial-power coefficient not positive')
        v.need(abs(sum(q.values(),v.Z)-1).upper()<arb('1e-60'),'radial-power coefficients do not sum to one')
        q_cases+=1
    d=json.load(open('center_r8_M45_S24.json'))
    rows,g,c,p,V,F=v.make_algebra(d)
    Q=q_of_g(g,d['M'],d['S'])
    q_error=check_q(g,Q)
    out={'status':'IDENTITIES_PASS','bits':ctx.prec,'K_cases':k_cases,
         'radial_power_cases':q_cases,
         'ball_principal_max_error':max_principal.str(30),
         'Q_connection_max_error':q_error.str(30)}
    with open('identities_r8.json','w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
