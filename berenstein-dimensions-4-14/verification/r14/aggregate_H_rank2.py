"""Readback of near H columns plus the analytic far source bound."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx


def val(pair):
    x=arb(int(pair[0]))*arb(2)**int(pair[1])
    require((x.is_finite()), 'aggregate_H_rank2.py:10: proof gate failed')
    return x


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--far',required=True)
    ap.add_argument('near',nargs='+')
    a=ap.parse_args()
    ctx.prec=192;ctx.threads=1
    sha=hashlib.sha256(open(a.center,'rb').read()).hexdigest()
    with open(a.far) as f:far=json.load(f)
    require((far['center_sha256']==sha and 'Hbound_upper_dyadic' in far), 'aggregate_H_rank2.py:23: proof gate failed')
    J=far['J'];m=far['m']
    peak=val(far['Hbound_upper_dyadic']);where='far'
    seen=set();total=None;first=None
    for path in a.near:
        with open(path) as f:r=json.load(f)
        require((r['center_sha256']==sha and r['J']==J and r['m']==m), 'aggregate_H_rank2.py:29: proof gate failed')
        if total is None:total=r['total'];first=r['first']
        require((r['total']==total and r['first']==first), 'aggregate_H_rank2.py:31: proof gate failed')
        require((len(r['column_upper_dyadic'])==r['end_index']-r['start_index']), 'aggregate_H_rank2.py:32: proof gate failed')
        for idx,(j,man,exp) in enumerate(r['column_upper_dyadic'],r['start_index']):
            require((idx not in seen and j==first+2*m*idx), 'aggregate_H_rank2.py:34: proof gate failed')
            seen.add(idx)
            x=val([man,exp])
            if x>peak:peak=x;where=j
    require((seen==set(range(total)) and first+2*m*total==J), 'aggregate_H_rank2.py:38: proof gate failed')
    require((peak<arb(250)), 'aggregate_H_rank2.py:39: proof gate failed')
    print('HIGH_SHAPE_SOURCE_NORM_ONLY','near columns',total,
          'far start',J,'H upper',peak,'argmax',where)


if __name__=='__main__':
    main()
