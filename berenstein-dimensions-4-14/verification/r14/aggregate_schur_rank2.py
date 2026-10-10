"""Fail-closed coverage/readback of near and far Schur certificates."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import arb,ctx
from exact_finite_rank2 import center


def dyad(pair):
    a,e=pair
    return arb(a)*arb(2)**e


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    ap.add_argument('--far',required=True)
    ap.add_argument('near',nargs='+')
    a=ap.parse_args()
    ctx.prec=256;ctx.threads=1
    sha=hashlib.sha256(open(a.center,'rb').read()).hexdigest()
    with open(a.center) as f:d=json.load(f)
    m=d['m'];ms=list(range(m,d['M']+1,2*m));js=list(range(1,1+2*m*len(ms),2*m))
    first=js[-1]+2*m
    with open(a.far) as f:far=json.load(f)
    require((far['center_sha256']==sha and far['m']==m), 'aggregate_schur_rank2.py:26: proof gate failed')
    require((far['rho_num']>far['rho_den']>0 and far['bits']>=128), 'aggregate_schur_rank2.py:27: proof gate failed')
    J=far['J'];require((J>first and (J-first)%(2*m)==0), 'aggregate_schur_rank2.py:28: proof gate failed')
    total=(J-first)//(2*m)
    require((far['kmin']+J>0), 'aggregate_schur_rank2.py:30: proof gate failed')
    far_upper=dyad(far['bound_upper_dyadic'])
    require((far_upper.is_finite() and far_upper<arb(1)), 'aggregate_schur_rank2.py:32: proof gate failed')
    spans=[];seen=[];peak=arb(0)
    for path in a.near:
        with open(path) as f:r=json.load(f)
        require((r['center_sha256']==sha and r['J']==J and r['total']==total), 'aggregate_schur_rank2.py:36: proof gate failed')
        require(((r['rho_num'],r['rho_den'])==(far['rho_num'],far['rho_den']) and r['bits']>=128), 'aggregate_schur_rank2.py:37: proof gate failed')
        spans.append((r['start_index'],r['end_index']))
        require((len(r['columns'])==r['end_index']-r['start_index']), 'aggregate_schur_rank2.py:39: proof gate failed')
        for idx,(j,man,exp) in enumerate(r['columns'],r['start_index']):
            require((j==first+2*m*idx), 'aggregate_schur_rank2.py:41: proof gate failed')
            val=dyad([man,exp])
            require((val.is_finite() and val<arb(1)), 'aggregate_schur_rank2.py:43: proof gate failed')
            if val>peak:peak=val
            seen.append(idx)
    require((sorted(seen)==list(range(total))), 'aggregate_schur_rank2.py:46: proof gate failed')
    require((sorted(spans)==sorted(set(spans))), 'aggregate_schur_rank2.py:47: proof gate failed')
    bound=peak if peak>far_upper else far_upper
    r14_center=(m,d['M'],d['S'])==(6,114,36)
    if r14_center:
        if not bound<arb(474)/1000:
            raise ValueError('R14 Schur contraction exceeds frozen ceiling')
    elif not bound<arb(1):
        raise ValueError('Schur contraction fails')
    tau=1/(1-(arb(far['rho_num'])/far['rho_den'])**(-2*m))
    _,_,c,B,_,_,_=center(a.center)
    hb=(arb(1)/2 if m==2 else arb(1))*c[1]**m/B
    alpha=hb*hb
    inverse_upper=(tau/(alpha*(1-bound))).upper()
    if not inverse_upper.is_finite():
        raise ValueError('nonfinite Schur inverse')
    if r14_center and not inverse_upper<arb(9):
        raise ValueError('R14 Schur inverse exceeds frozen ceiling')
    print('SCHUR_RECEIPTS_PASS_ONLY')
    print('m',m,'near columns',total,'near maximum',peak,
          'far upper',far_upper,'global K upper',bound,
          'Schur inverse norm upper',inverse_upper)


if __name__=='__main__':
    main()
