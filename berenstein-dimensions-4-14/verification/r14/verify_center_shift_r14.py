"""Arb transfer of the R14 derivative defect to the higher-precision center.

Requires the separately certified old-center Z<0.49 and A<100000 blocks.
"""
from proof_guard import require
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center
from exact_finite_r4 import ZERO


def upper(x):
    require((x.is_finite()), 'verify_center_shift_r14.py:13: proof gate failed')
    y=x.upper();require((y.is_finite()), 'verify_center_shift_r14.py:14: proof gate failed')
    return y


def main():
    ctx.prec=192;ctx.threads=1
    old='center_r14_M114_S36_arbrefined.json'
    new='center_r14_M186_S48_arbrefined2.json'
    d1,g1,c1,B1,*_=center(old)
    d2,g2,c2,B2,*_=center(new)
    require(((d1['m'],d2['m'])==(6,6) and B1==B2), 'verify_center_shift_r14.py:24: proof gate failed')
    rho=arb(102)/100;m=6
    dg=sum((rho**n*abs(g1.get((n,s),ZERO)-g2.get((n,s),ZERO))
            for n,s in set(g1)|set(g2)),ZERO)
    dc=sum(((j+2*m)*rho**(j-1)*abs(c1.get(j,ZERO)-c2.get(j,ZERO))
            for j in set(c1)|set(c2)),ZERO)
    distance=upper(dg+dc)
    require((distance<arb(7)/arb(10)**18), 'verify_center_shift_r14.py:31: proof gate failed')
    with open('r14_M114_nonlinear_bounds.json') as f:r=json.load(f)
    require((r['center_sha256']==hashlib.sha256(open(old,'rb').read()).hexdigest()), 'verify_center_shift_r14.py:33: proof gate failed')
    require(((r['rho_num'],r['rho_den'])==(102,100)), 'verify_center_shift_r14.py:34: proof gate failed')
    man,exp=r['derivative_sum_upper_dyadic']
    lip=arb(int(man))*arb(2)**int(exp)
    if not lip.is_finite() or not lip<arb(38):
        raise ValueError('old-center nonlinear derivative bound fails')
    Znew=arb(49)/100+arb(100000)*lip*distance
    require((upper(Znew)<arb(1)/2), 'verify_center_shift_r14.py:40: proof gate failed')
    print('CENTER_SHIFT_DERIVATIVE_ONLY','distance upper',distance,
          'old nonlinear derivative sum upper',lip,
          'new Z upper',upper(Znew))


if __name__=='__main__':
    main()
