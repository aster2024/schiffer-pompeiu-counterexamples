#!/usr/bin/env python3
"""Compare typed interval values with the paper's rational ceilings.

Run under nice -n 19 next to the coefficient programs and certificates.
The comparisons check the numerical bounds of Sections 5 and 6."""
import argparse
import builtins
from fractions import Fraction
import json
import os
from pathlib import Path
import runpy
import sys
from flint import arb, ctx


def rational(value):
    q = Fraction(str(value))
    return arb(q.numerator)/q.denominator


def number(value):
    if isinstance(value,Fraction):
        return arb(value.numerator)/value.denominator
    return value if isinstance(value,arb) else arb(value)


def compare(value,cap,label,lower=False):
    x=number(value); c=rational(cap)
    if not x.is_finite() or not ((x>c) if lower else (x<c)):
        raise ValueError('numerical inequality failed: '+label+' '+str(x))
    print('INEQUALITY',label,'>',cap if lower else '',flush=True) if lower else print('INEQUALITY',label,'<',cap,flush=True)


def run(script,args=()):
    records=[]
    original_print=builtins.print
    original_argv=sys.argv
    def capture(*values,**kwargs):
        records.append(values)
        original_print(*values,**kwargs)
    builtins.print=capture
    sys.argv=[script,*map(str,args)]
    try:
        runpy.run_path(script,run_name='__main__')
    finally:
        builtins.print=original_print
        sys.argv=original_argv
    return records


def field(records,label):
    matches=[]
    for row in records:
        for i,v in enumerate(row[:-1]):
            if isinstance(v,str) and v==label:
                matches.append(row[i+1])
    if len(matches)!=1:
        raise ValueError('ambiguous printed field: '+label)
    return matches[0]


def caps(records,items):
    for label,cap in items:
        compare(field(records,label),cap,label)


def dyad(pair):
    a,e=pair
    if type(a) is not int or type(e) is not int:
        raise ValueError('invalid dyadic pair')
    return arb(a)*arb(2)**e


def r4():
    batches=['jac_batch_0_104.json','jac_batch_104_208.json',
             'jac_batch_208_312.json','jac_batch_312_416.json']
    out=run('verify_finite_inverse_r4.py',['--check-arb-inverse','Ahat_arb_r4_M31_S24.json',*batches])
    caps(out,[('normA upper','1665.434'),('defect upper','2.05e-34'),
              ('last boundary column weighted norm upper','54.778'),
              ('preconditioned finite residual upper','1.198868e-6'),
              ('preconditioned tail residual upper','1.838234e-6'),
              ('preconditioned total residual upper','3.038e-6')])
    # This is the tighter finite-defect ceiling used by the standalone driver.
    compare(field(out,'defect upper'),'21e-35','finite defect entering whole-column bound')
    for pattern,cap in [('tailcol_*.json','0.578'),('tail_source_*.json','0.679'),('tail_shape_*.json','0.720')]:
        paths=sorted(Path('.').glob(pattern))
        if not paths:
            raise ValueError('missing receipts')
        for path in paths:
            compare(dyad(json.loads(path.read_text())['max_upper_dyadic']),cap,path.name)
    caps(run('verify_far_source_bound_r4.py'),[
        ('P upper','9.09'),('Cq upper','86.059'),
        ('far angular m>=63 upper','0.885833'),('far radial s>=56 upper','0.006359')])
    caps(run('verify_far_shape_bound_r4.py'),[
        ('boundary upper','0.706765'),('interior upper','0.076108'),
        ('trace coupling upper','0.142951'),('finite cross upper','5.23e-10'),
        ('total upper','0.925823')])
    caps(run('verify_nonlinear_r4.py'),[
        ('P','9.089833'),('P2','80.886957'),('H0','7.543929'),
        ('H1','9.089833'),('H20','57.308957'),('Vnorm','36.702830'),
        ('normA upper','2444'),('C2 upper','268'),('C3 upper','0.0023'),('C4 upper','6e-8')])
    geometry=run('verify_geometry_r4.py')
    for name,cap in [('univalence derivative','5.39635'),('p modulus','5.66945'),
                     ('convexity criterion','0.10935'),('star-shaped criterion','0.84029'),
                     ('nonball c3','0.1531'),('q sign','24.345')]:
        rows=[row for row in geometry if row and row[0]==name]
        if len(rows)!=1 or rows[0][1]!='lower':
            raise ValueError('missing geometry output')
        compare(rows[0][2],cap,name,lower=True)
    run('verify_signchange_r4.py')
    run('verify_radii_r4.py')


def r14():
    old='center_r14_M114_S36_arbrefined.json'
    new='center_r14_M186_S48_arbrefined2.json'
    inverse='Ahat_r14_M114_S36.npy'
    jac=sorted(p.name for p in Path('.').glob('r14_M114_jac_*.json'))
    out=run('verify_finite_block_rank2.py',['--center',old,'--inverse',inverse,*jac])
    caps(out,[('normA upper','82859'),('defect upper','2.161e-9')])
    caps(run('verify_source_tail_rank2.py',[old]),[('delta upper','0.141426')])
    far=run('far_schur_bound_rank2.py',[old,'--J','361','--rho-num','102','--rho-den','100'])
    caps(far,[('Hbound','249.711'),('shapeK','0.275974'),('traceK','0.197588'),('far Schur K upper','0.473562')])
    near=sorted(p.name for p in Path('.').glob('r14_M114_schur_[0-9]*.json'))
    schur=run('aggregate_schur_rank2.py',[old,'--far','r14_M114_schur_far.json',*near])
    caps(schur,[('near maximum','0.059850'),('global K upper','0.474'),
                ('Schur inverse norm upper','9')])
    caps(run('bound_finite_cross_source_rank2.py',[old,'r14_M114_finite_inverse_profile.json']),[
        ('maximum upper','0.171168')])
    cross=sorted(p.name for p in Path('.').glob('r14_M114_cross_shape_*.json'))
    caps(run('aggregate_finite_cross_shape_rank2.py',[old,inverse,*cross]),[
        ('maximum upper','9.183045')])
    hs=sorted(p.name for p in Path('.').glob('r14_M114_H_*.json'))
    caps(run('aggregate_H_rank2.py',[old,'--far','r14_M114_schur_far.json',*hs]),[('H upper','250')])
    caps(run('bound_QD_source_rank2.py',[old]),[('upper','0.000112959')])
    caps(run('bound_preconditioner_norm_r14.py'),[
        ('H upper','1076'),('Q upper','0.133'),('A upper','94580')])
    caps(run('bound_high_source_defect_r14.py'),[('whole column upper','0.49')])
    caps(run('bound_high_shape_defect_r14.py'),[('whole column upper','0.018962')])
    caps(run('aggregate_finite_input_coarse_r14.py'),[('all-but-379 maximum','0.114825')])
    caps(run('verify_last_finite_column_r14.py'),[
        ('Schur residual','3.573e-15'),('whole column upper','0.160581')])
    caps(run('nonlinear_majorant_rank2.py',[old]),[('sum q*dq upper','37.149')])
    caps(run('verify_center_shift_r14.py'),[
        ('distance upper','6.440e-18'),('old nonlinear derivative sum upper','37.149'),
        ('new Z upper','0.490000000024')])
    for name in ['r14_M114_nonlinear_bounds.json','r14_M186_nonlinear_bounds.json']:
        data=json.loads(Path(name).read_text())
        compare(dyad(data['coeff_upper_dyadic'][0]),'18.551',name+' d2')
        compare(dyad(data['derivative_sum_upper_dyadic']),'37.149',name+' derivative sum')
    out=run('verify_radii_geometry_r14.py')
    caps(out,[('Y upper','1.063e-25'),('fixed-boundary-scale factor upper','1.000000000000000000000017'),
              ('map gap upper','-4.99989e-21'),('Lipschitz gap upper','-0.49999999999996')])
    for label,cap in [('convex numerator margin lower','25.4946'),('q lower','29.1875'),('h lower','0.98405')]:
        compare(field(out,label),cap,label,lower=True)
    signs=field(out,'sign values')
    compare(signs,'-0.0795470','negative center sign')
    rows=[row for row in out if any(isinstance(v,str) and v=='sign values' for v in row)]
    row=rows[0]; index=next(i for i,v in enumerate(row) if isinstance(v,str) and v=='sign values')
    compare(row[index+2],'1.4919379','positive center sign',lower=True)


def main():
    if sys.flags.optimize:
        raise SystemExit('optimized Python disables proof assertions')
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default=str(Path(__file__).resolve().parent))
    a=ap.parse_args()
    if os.getpriority(os.PRIO_PROCESS,0)!=19:
        raise SystemExit('run under nice -n 19')
    root=Path(a.root).resolve()
    os.chdir(root)
    sys.path.insert(0,str(root))
    ctx.threads=1
    if (root/'finite_center_r4_M31_S24.json').exists():
        r4()
    elif (root/'center_r14_M114_S36_arbrefined.json').exists():
        r14()
    else:
        raise SystemExit('certificate sources not found')
    print('ALL_PAPER_INEQUALITIES_CHECKED',flush=True)


if __name__=='__main__':
    main()
