"""Rational assembly of the all-column R14 high-source defect bound.

Requires independently checked component ceilings from the companion Arb
scripts. It does not cover finite input columns or high-shape columns.
"""
from fractions import Fraction as F


def main():
    delta=F(142,1000)       # verify_source_tail_rank2: < 0.141426
    qd=F(114,1000000)      # bound_QD_source_rank2: < 0.000113
    cinv=F(9)              # aggregate_schur_rank2: < 8.982
    h=F(250)               # aggregate_H_rank2: < 249.711
    xs=F(1,5)             # bound_finite_cross_source_rank2: < 0.172
    xc=F(10)              # aggregate_finite_cross_shape_rank2: < 9.184
    c=cinv*qd
    g=delta+h*c
    z=(1+xs)*g+(1+xc)*c
    if not z<F(49,100):
        raise ValueError('high-source defect exceeds old-center Z ceiling')
    print('HIGH_SOURCE_DEFECT_BOUND_ONLY','shape correction',c,
          'source correction',g,'whole column upper',z,'decimal',float(z))


if __name__=='__main__':
    main()
