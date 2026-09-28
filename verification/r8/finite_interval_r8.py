"""Arb finite Jacobian, frozen inverse defect, finite columns, residual."""

import hashlib
import json
from flint import arb, arb_mat, ctx
from tail_inverse_r8 import load_inverse
import verify_r8 as v


def sha(path):
    return hashlib.sha256(open(path,'rb').read()).hexdigest()


def main():
    ctx.prec=128
    ctx.threads=1
    centre='center_r8_M45_S24.json'
    inverse='inverse_r8_M45_S24.npy'
    d=json.load(open(centre))
    rows,g,c,p,V,F=v.make_algebra(d)
    N=len(rows)
    rho=arb(11)/10
    A=load_inverse(inverse)
    inv=v.Inverse(d['M'],d['S'],rho,c[1],A)
    cols=v.finite_columns(rows,g,c,p,V,d['S'])
    mat=arb_mat([[col.get(key,v.Z) for col in cols] for key in rows])
    E=arb_mat(N,N)
    for i in range(N):
        E[i,i]=1
    E=E-A*mat
    defect=v.max_upper(sum((abs(E[i,j])*inv.cw[i] for i in range(N)),v.Z)/inv.cw[j]
                       for j in range(N))
    Zfinite=v.max_upper(val/inv.cw[i] for i,val in enumerate(
        inv.apply_many(cols,subtract=[('finite',i) for i in range(N)])))
    Y=inv.apply_many([F])[0].upper()
    Anfin=v.max_upper(sum((abs(A[i,j])*inv.cw[i] for i in range(N)),v.Z)/rho**rows[j][0]
                      for j in range(N))
    out={'status':'OPEN_COMPONENT','N':N,'bits':ctx.prec,
         'centre_sha256':sha(centre),'inverse_sha256':sha(inverse),
         'inverse_defect':defect.str(30),'Z_finite':Zfinite.str(30),
         'Y':Y.str(30),'A_finite':Anfin.str(30)}
    with open('finite_bound_r8.json','w') as f:
        json.dump(out,f,indent=2)
        f.write('\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
