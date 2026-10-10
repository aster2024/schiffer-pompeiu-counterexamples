"""Readback of all 380 finite input columns under the coarse tail bound."""
from proof_guard import require
import hashlib
import json
import glob
from flint import arb,ctx


def val(pair):
    x=arb(int(pair[0]))*arb(2)**int(pair[1])
    require((x.is_finite()), 'aggregate_finite_input_coarse_r14.py:10: proof gate failed')
    return x


def main():
    ctx.prec=192;ctx.threads=1
    name='center_r14_M114_S36_arbrefined.json'
    sha=hashlib.sha256(open(name,'rb').read()).hexdigest()
    files=glob.glob('r14_M114_finite_coarse_*.json')
    require((len(files)==8), 'aggregate_finite_input_coarse_r14.py:19: proof gate failed')
    seen=set();needs=[];peak=arb(0);peakidx=-1
    for path in files:
        with open(path) as f:r=json.load(f)
        require((r['center_sha256']==sha and r['dimension']==380), 'aggregate_finite_input_coarse_r14.py:23: proof gate failed')
        require((len(r['records'])==r['end']-r['start']), 'aggregate_finite_input_coarse_r14.py:24: proof gate failed')
        require((r['needs_refinement']==[q['i'] for q in r['records']
                                       if val(q['coarseZ'])>=arb(1)]), 'aggregate_finite_input_coarse_r14.py:25: proof gate failed')
        for q in r['records']:
            i=q['i'];require((r['start']<=i<r['end'] and i not in seen), 'aggregate_finite_input_coarse_r14.py:28: proof gate failed')
            seen.add(i)
            z=val(q['coarseZ'])
            if z>=arb(1):needs.append(i)
            elif z>peak:peak=z;peakidx=i
    require((seen==set(range(380)) and needs==[379]), 'aggregate_finite_input_coarse_r14.py:33: proof gate failed')
    require((peak<arb(1)/8), 'aggregate_finite_input_coarse_r14.py:34: proof gate failed')
    print('FINITE_INPUT_COARSE_COVERAGE_ONLY','columns',len(seen),
          'all-but-379 maximum',peak,'argmax',peakidx,
          'one remaining column',needs)


if __name__=='__main__':
    main()
