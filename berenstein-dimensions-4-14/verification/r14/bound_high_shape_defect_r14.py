"""Rational all-column bound for R14 high-shape defect of 13-step tail inverse."""
from proof_guard import require
from fractions import Fraction as F


def main():
    k=F(474,1000)         # Schur K upper < 0.473562
    H=F(250)              # high-shape source derivative norm
    Xs=F(1,5)             # finite cross from high source
    Xc=F(10)              # finite cross from high shape
    steps=13
    z=((1+Xs)*H+(1+Xc))*k**steps
    require((z<F(1,50)), 'bound_high_shape_defect_r14.py:12: proof gate failed')
    print('HIGH_SHAPE_DEFECT_BOUND_ONLY','Neumann steps',steps,
          'whole column upper',z,'decimal',float(z))


if __name__=='__main__':
    main()
