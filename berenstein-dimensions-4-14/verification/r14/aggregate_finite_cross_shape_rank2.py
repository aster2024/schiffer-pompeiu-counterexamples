"""Read back exact coverage of high-shape to finite-block cross receipts."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('inverse')
    ap.add_argument('receipts',nargs='+')
    a=ap.parse_args()
    ctx.prec=192;ctx.threads=1
    raw=open(a.center,'rb').read();d=json.loads(raw)
    sha=hashlib.sha256(raw).hexdigest()
    invsha=hashlib.sha256(open(a.inverse,'rb').read()).hexdigest()
    m,M=d['m'],d['M']
    ms=list(range(m,M+1,2*m));js=list(range(1,1+2*m*len(ms),2*m))
    first=js[-1]+2*m
    jmax=js[-1];pmax=jmax-1;hmax=m*jmax;h2max=2*hmax-m
    jzero=first
    while not (jzero+m-1-hmax-pmax>M and
               jzero-1-pmax-h2max>M and
               jzero-1-M-pmax>M):
        jzero+=2*m
    total=(jzero-first)//(2*m)
    seen=set();best=arb(0);bestj=-1
    for path in a.receipts:
        with open(path) as f:r=json.load(f)
        require((r['center_sha256']==sha and r['inverse_sha256']==invsha), 'aggregate_finite_cross_shape_rank2.py:31: proof gate failed')
        require(((r['m'],r['M'],r['first'],r['jzero'],r['total'])==(m,M,first,jzero,total)), 'aggregate_finite_cross_shape_rank2.py:32: proof gate failed')
        require((r['bits']>=128 and len(r['column_upper_dyadic'])==r['end_index']-r['start_index']), 'aggregate_finite_cross_shape_rank2.py:33: proof gate failed')
        for idx,(j,man,exp) in enumerate(r['column_upper_dyadic'],r['start_index']):
            require((idx not in seen and j==first+2*m*idx), 'aggregate_finite_cross_shape_rank2.py:35: proof gate failed')
            seen.add(idx)
            val=arb(int(man))*arb(2)**int(exp)
            require((val.is_finite()), 'aggregate_finite_cross_shape_rank2.py:38: proof gate failed')
            if val>best:best=val;bestj=j
    require((seen==set(range(total))), 'aggregate_finite_cross_shape_rank2.py:40: proof gate failed')
    require((best<arb(10)), 'aggregate_finite_cross_shape_rank2.py:41: proof gate failed')
    print('FINITE_CROSS_SHAPE_BOUND_ONLY','m',m,'columns',total,
          'first',first,'zero_from',jzero,
          'maximum upper',best,'argmax',bestj)


if __name__=='__main__':
    main()
