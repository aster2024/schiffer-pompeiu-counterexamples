#!/usr/bin/env python3
"""Fail-closed entrypoint for the R4 Berenstein interval proof package.

Run in this directory with python-flint 0.9.0. The default authenticates
the bounded-batch receipts and reruns the finite inverse, all analytic
far-tail bounds, nonlinear bounds, radii inequalities, and geometry.
`--recompute` also regenerates every batch receipt before printing PROVED.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import flint

ROOT=Path(__file__).resolve().parent
EXPECTED_MANIFEST_SHA256='24b1e46952290b5d634bfa5507b6f0d98bdc6d64a4b9781e1d5734e320b02482'
CENTER='finite_center_r4_M31_S24.json'
INVERSE='Ahat_arb_r4_M31_S24.json'
JAC_SPANS=[(0,104),(104,208),(208,312),(312,416)]
FINITE_SPANS=JAC_SPANS
SOURCE_SPANS=[(a,min(a+100,1336)) for a in range(0,1336,100)]
SHAPE_SPANS=[(0,23),(23,46),(46,58)]+[(a,min(a+10,134)) for a in range(58,134,10)]


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def dyadic(pair):
    m,e=pair
    assert type(m) is int and type(e) is int and abs(e)<1000000
    return Fraction(m*(2**e),1) if e>=0 else Fraction(m,2**(-e))


def run_stage(name,*args):
    env=dict(os.environ)
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1')
    proc=subprocess.run([sys.executable,'-B',str(ROOT/name),*map(str,args)],
                        cwd=ROOT,env=env,capture_output=True,text=True)
    if proc.returncode:
        raise RuntimeError(f'{name} failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}')
    return proc.stdout


def expected_files():
    names=[CENTER,INVERSE,'exact_finite_r4.py','verify_finite_inverse_r4.py',
           'verify_tail_source_batch_r4.py','verify_tail_shape_batch_r4.py',
           'verify_far_source_bound_r4.py','verify_far_shape_bound_r4.py',
           'verify_nonlinear_r4.py','verify_geometry_r4.py','verify_radii_r4.py',
           'verify_factor_exact_r4.py','verify_signchange_r4.py',
           'verify_ball_linearization.py']
    names +=[f'jac_batch_{a}_{b}.json' for a,b in JAC_SPANS]
    names +=[f'tailcol_{a}_{b}.json' for a,b in FINITE_SPANS]
    names +=[f'tail_source_{a}_{b}.json' for a,b in SOURCE_SPANS]
    names +=[f'tail_shape_{a}_{b}.json' for a,b in SHAPE_SPANS]
    return names


def authenticate_manifest():
    manifest=ROOT/'CERTIFICATE_MANIFEST.json'
    assert sha(manifest)==EXPECTED_MANIFEST_SHA256,'certificate manifest changed'
    with open(manifest) as f:data=json.load(f)
    assert data['version']==1 and set(data['files'])==set(expected_files())
    for name,expected in data['files'].items():
        assert sha(ROOT/name)==expected,f'file changed: {name}'
    return data


def verify_receipts():
    center_sha=sha(ROOT/CENTER);inverse_sha=sha(ROOT/INVERSE)
    groups=[('tailcol',FINITE_SPANS,416,Fraction(463,500),False),
            ('tail_source',SOURCE_SPANS,1336,Fraction(463,500),True),
            ('tail_shape',SHAPE_SPANS,134,Fraction(463,500),True)]
    maxima={}
    for prefix,spans,total,cap,has_inverse in groups:
        cursor=0;greatest=Fraction(0)
        for start,end in spans:
            assert start==cursor and 0<=start<end<=total
            cursor=end
            path=ROOT/f'{prefix}_{start}_{end}.json'
            with open(path) as f:d=json.load(f)
            assert d['center_sha256']==center_sha and d['bits']==128
            assert (d['start'],d['end'])==(start,end)
            assert d['dimension' if prefix=='tailcol' else 'mode_count']==total
            if has_inverse:assert d['inverse_sha256']==inverse_sha
            vals=[dyadic(pair) for pair in d['column_upper_dyadic']]
            assert len(vals)==end-start and all(v>=0 for v in vals)
            bound=dyadic(d['max_upper_dyadic'])
            assert bound==max(vals) and d['argmax']==start+vals.index(bound)
            assert bound<cap,f'{prefix} batch failed: {start}:{end}'
            greatest=max(greatest,bound)
        assert cursor==total
        maxima[prefix]=greatest
    return maxima


def recompute_receipts():
    with tempfile.TemporaryDirectory(prefix='r4_recompute_') as td:
        scratch=Path(td)
        for prefix,spans in [('jac_batch',JAC_SPANS),('tailcol',FINITE_SPANS),
                             ('tail_source',SOURCE_SPANS),('tail_shape',SHAPE_SPANS)]:
            for start,end in spans:
                name=f'{prefix}_{start}_{end}.json';out=scratch/name
                if prefix=='jac_batch':
                    run_stage('exact_finite_r4.py','--center',CENTER,'--bits','128',
                              '--batch',f'{start}:{end}','--output',out)
                elif prefix=='tailcol':
                    run_stage('exact_finite_r4.py','--center',CENTER,'--bits','128',
                              '--tail-column-batch',f'{start}:{end}','--output',out)
                elif prefix=='tail_source':
                    run_stage('verify_tail_source_batch_r4.py','--batch',
                              f'{start}:{end}','--output',out)
                else:
                    run_stage('verify_tail_shape_batch_r4.py','--batch',
                              f'{start}:{end}','--output',out)
                assert out.read_bytes()==(ROOT/name).read_bytes(),f'not reproducible: {name}'


def main():
    if sys.flags.optimize:
        raise RuntimeError('Python optimisation disables certificate assertions')
    p=argparse.ArgumentParser()
    p.add_argument('--recompute',action='store_true',help='regenerate all bounded batches')
    args=p.parse_args()
    assert flint.__version__=='0.9.0','python-flint 0.9.0 required'
    authenticate_manifest()
    if args.recompute:recompute_receipts()
    out=run_stage('verify_finite_inverse_r4.py','--check-arb-inverse',INVERSE,
                  *(f'jac_batch_{a}_{b}.json' for a,b in JAC_SPANS))
    assert 'CERTIFIED FINITE BLOCK ONLY' in out
    maxima=verify_receipts()
    for name,marker in [('verify_far_source_bound_r4.py','CERTIFIED FAR SOURCE BOUNDS ONLY'),
                        ('verify_far_shape_bound_r4.py','CERTIFIED FAR SHAPE BOUND ONLY'),
                        ('verify_factor_exact_r4.py','EXACT DIRICHLET FACTOR VERIFIED'),
                        ('verify_signchange_r4.py','CERTIFIED SIGN CHANGE'),
                        ('verify_ball_linearization.py','CERTIFIED BALL LINEARIZATION'),
                        ('verify_nonlinear_r4.py','C4 upper'),
                        ('verify_geometry_r4.py','CERTIFIED GEOMETRY MARGINS'),
                        ('verify_radii_r4.py','RADII INEQUALITIES PASS')]:
        assert marker in run_stage(name),f'{name} did not report its checked bounds'
    # The exact finite inverse defect is below 2.05e-34, so this bound plus
    # the finite-column maximum must still lie below the rational Z cap.
    assert maxima['tailcol']+Fraction(21,10**35)<Fraction(463,500)
    print('PROVED')


if __name__=='__main__':
    main()
