#!/usr/bin/env python3
"""Dimension-fourteen certificate entrypoint.

The SHA-pinned inputs and all inequalities of Sections 6 and 7 are checked
before PROVED. Use --recompute-all to regenerate every bounded column file.
"""
from proof_guard import require
import argparse
import glob
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def gate(script,*args,quiet=True,timeout_seconds=900):
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    cmd=[sys.executable,'-B',script,*map(str,args)]
    p=subprocess.run(cmd,cwd=ROOT,env=env,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                     timeout=timeout_seconds,check=False,preexec_fn=lambda:os.nice(19))
    if p.returncode:
        print(p.stdout,file=sys.stderr)
        raise RuntimeError(f'gate failed ({p.returncode}): {script}')
    if not quiet:
        print(p.stdout.rstrip())
    elif p.stdout.strip():
        print(script,':',p.stdout.strip().splitlines()[-1],flush=True)


def same_file(recomputed,original):
    require((sha(recomputed)==sha(original)), (recomputed,original))
    require(((ROOT/recomputed).read_bytes()==(ROOT/original).read_bytes()),
            ('recomputed bytes differ',recomputed,original))
    (ROOT/recomputed).unlink()


def receipt_spans(pattern):
    return sorted(Path(p).name for p in glob.glob(str(ROOT/pattern)))


def recompute_all():
    old='center_r14_M114_S36_arbrefined.json'
    for original in receipt_spans('r14_M114_jac_*.json'):
        data=json.loads((ROOT/original).read_text())
        temp='recompute_'+original
        gate('exact_finite_rank2.py','--center',old,'--bits',data['bits'],
             '--batch',f"{data['start']}:{data['end']}",'--output',temp,
             timeout_seconds=900)
        same_file(temp,original)
    for original in receipt_spans('r14_M114_schur_[0-9]*.json'):
        data=json.loads((ROOT/original).read_text())
        temp='recompute_'+original
        gate('cert_near_schur_rank2.py',old,'--J',data['J'],
             '--rho-num',data['rho_num'],'--rho-den',data['rho_den'],
             '--bits',data['bits'],'--start-index',data['start_index'],
             '--end-index',data['end_index'],'--output',temp,
             timeout_seconds=900)
        same_file(temp,original)
    for original in receipt_spans('r14_M114_cross_shape_*.json'):
        data=json.loads((ROOT/original).read_text())
        temp='recompute_'+original
        gate('cert_finite_cross_shape_rank2.py',old,'Ahat_r14_M114_S36.npy',
             '--bits',data['bits'],'--start-index',data['start_index'],
             '--end-index',data['end_index'],'--output',temp,
             timeout_seconds=900)
        same_file(temp,original)
    for original in receipt_spans('r14_M114_H_*.json'):
        data=json.loads((ROOT/original).read_text())
        temp='recompute_'+original
        gate('cert_H_near_rank2.py',old,'--J',data['J'],
             '--start-index',data['start_index'],'--end-index',data['end_index'],
             '--output',temp,timeout_seconds=900)
        same_file(temp,original)
    for original in receipt_spans('r14_M114_finite_coarse_*.json'):
        data=json.loads((ROOT/original).read_text())
        temp='recompute_'+original
        gate('cert_finite_input_coarse_rank2.py',old,
             '--start',data['start'],'--end',data['end'],'--output',temp,
             timeout_seconds=900)
        same_file(temp,original)


def main():
    if not __debug__:
        raise SystemExit('optimized Python disables proof assertions')
    ap=argparse.ArgumentParser()
    ap.add_argument('--recompute-all',action='store_true')
    a=ap.parse_args()
    manifest=json.loads((ROOT/'CERTIFICATE_MANIFEST_R14.json').read_text())
    require((manifest['kind']=='R14_CERTIFICATE_V2'), 'R14 certificate kind mismatch')
    for name,digest in manifest['files'].items():
        path=Path(name)
        require((not path.is_absolute() and '..' not in path.parts), 'verify_r14_package.py:103: proof gate failed')
        require((sha(name)==digest), name)
    print('PAYLOAD_SHA256_PASS',len(manifest['files']),'files',flush=True)
    if a.recompute_all:
        recompute_all()
        print('ALL_RECEIPTS_RECOMPUTED',flush=True)

    old='center_r14_M114_S36_arbrefined.json'
    new='center_r14_M186_S48_arbrefined2.json'
    inv='Ahat_r14_M114_S36.npy'
    gate('exact_finite_rank2.py','--center',old,'--bits',128,
         '--check-rho-num',102,'--check-rho-den',100,'--max-residual','1e-17')
    gate('verify_factor_dirichlet_rank2.py',old)
    gate('verify_source_tail_rank2.py',old,'--rho-num',102,'--rho-den',100)
    finite_receipts=receipt_spans('r14_M114_jac_*.json')
    require((len(finite_receipts)==8), 'verify_r14_package.py:118: proof gate failed')
    gate('verify_finite_block_rank2.py','--center',old,'--inverse',inv,
         '--batch-width',50,'--rho-num',102,'--rho-den',100,
         '--save-profile','recompute_profile.json',*finite_receipts)
    same_file('recompute_profile.json','r14_M114_finite_inverse_profile.json')
    gate('far_schur_bound_rank2.py',old,'--J',361,
         '--rho-num',102,'--rho-den',100,
         '--output','recompute_far.json')
    same_file('recompute_far.json','r14_M114_schur_far.json')
    gate('aggregate_schur_rank2.py',old,'--far','r14_M114_schur_far.json',
         'r14_M114_schur_0_5.json','r14_M114_schur_5_10.json',
         'r14_M114_schur_10_15.json','r14_M114_schur_15_20.json')
    gate('bound_finite_cross_source_rank2.py',old,
         'r14_M114_finite_inverse_profile.json')
    cross=receipt_spans('r14_M114_cross_shape_*.json')
    require((len(cross)==11), 'verify_r14_package.py:133: proof gate failed')
    gate('aggregate_finite_cross_shape_rank2.py',old,inv,*cross)
    gate('aggregate_H_rank2.py',old,'--far','r14_M114_schur_far.json',
         'r14_M114_H_0_5.json','r14_M114_H_5_10.json',
         'r14_M114_H_10_15.json','r14_M114_H_15_20.json')
    gate('bound_QD_source_rank2.py',old)
    gate('bound_preconditioner_norm_r14.py')
    gate('bound_high_source_defect_r14.py')
    gate('bound_high_shape_defect_r14.py')
    gate('aggregate_finite_input_coarse_r14.py')
    gate('verify_last_finite_column_r14.py')
    gate('nonlinear_majorant_rank2.py',old,'--rho-num',102,'--rho-den',100,
         '--bound-output','recompute_old_nonlinear.json')
    same_file('recompute_old_nonlinear.json','r14_M114_nonlinear_bounds.json')
    gate('verify_center_shift_r14.py')
    gate('verify_factor_dirichlet_rank2.py',new)
    gate('nonlinear_majorant_rank2.py',new,'--rho-num',102,'--rho-den',100,
         '--bound-output','recompute_new_nonlinear.json')
    same_file('recompute_new_nonlinear.json','r14_M186_nonlinear_bounds.json')
    gate('verify_radii_geometry_r14.py')
    print('PROVED',flush=True)


if __name__=='__main__':
    main()
