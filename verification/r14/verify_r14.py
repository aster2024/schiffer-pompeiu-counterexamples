#!/usr/bin/env python3
"""Fail-closed Arb CAP verifier for the D12-symmetric Schiffer domain in R14."""

import hashlib
import json
import math
from pathlib import Path
from flint import arb, ctx
import rank2_cap_core as v

FROZEN_BUNDLE_SHA256='3d6973f5391d99816ab1621bc84b1be144dbbc15f730cd097cc31ee5922156d8'
FINITE_RANGES=[(0,20)]+[(a,a+40) for a in range(20,580,40)]
FINITE_RANGES += [(580,600)]+[(a,a+1) for a in range(600,610)]
G_RANGES=[(0,20)]+[(a,a+100) for a in range(20,1120,100)]
G_RANGES += [(1120,1215)]
SHAPE_RANGES=[(121,121,1),(133,145,2),(157,169,2)]
SHAPE_RANGES += [(j,j,1) for j in (181,193,205,217,229)]
FROZEN_FILES=[
    'rank2_cap_core.py','center_r14_M114_S60.json',
    'inverse_r14_M114_S60.npy','make_inverse_rank2.py',
    'finite_rank2.py','finite_batch_rank2.py','aggregate_finite_r14.py',
    'finite_bound_r14_M114_S60.json',
    'tail_inverse_rank2.py','tail_inverse_bound_r14_M114_S60.json',
    'g_tail_rank2.py','g_far_bound_r14_M114_S60.json',
    'shape_tail_rank2.py','far_shape_rank2.py',
    'far_shape_split_r14.py','far_shape_split_r14_j241.json',
    'check_identities_rank2.py','identities_r14_M114_S60.json',
]+[f'finite_r14_{a:03d}_{b:03d}.json' for a,b in FINITE_RANGES]
FROZEN_FILES += [f'g_near_r14_{a:04d}_{b:04d}.json' for a,b in G_RANGES]
FROZEN_FILES += [f'shape_near_r14_{a:03d}_{b:03d}.json' for a,b,_ in SHAPE_RANGES]


def digest():
    h=hashlib.sha256()
    for name in FROZEN_FILES:
        b=Path(name).read_bytes()
        h.update(name.encode()+b'\0')
        h.update(len(b).to_bytes(8,'big'))
        h.update(b)
    return h.hexdigest()


def need(ok,why):
    if not ok:
        raise RuntimeError(why)


def get(name):
    return json.loads(Path(name).read_text())


def interval(value):
    x=arb(value)
    need(x.is_finite() and x.upper().is_finite(),'nonfinite interval in receipt')
    return x


def main():
    ctx.prec=128
    ctx.threads=1
    need(len(FROZEN_BUNDLE_SHA256)==64 and digest()==FROZEN_BUNDLE_SHA256,
         'frozen R14 bundle missing or changed')
    d=get('center_r14_M114_S60.json')
    rows,g,c,p,V,F=v.make_algebra(d)
    need(d['m']==6 and d['step']==12 and d['M']==114 and d['S']==60
         and len(rows)==610,'unsupported R14 centre')
    rho=arb(11)/10
    finite=get('finite_bound_r14_M114_S60.json')
    need(finite['status']=='OPEN_COMPONENT' and finite['dimension']==14
         and finite['N']==610 and finite['bits']>=128
         and finite['batch_count']==len(FINITE_RANGES)
         and finite['centre_sha256']==hashlib.sha256(Path('center_r14_M114_S60.json').read_bytes()).hexdigest()
         and finite['inverse_sha256']==hashlib.sha256(Path('inverse_r14_M114_S60.npy').read_bytes()).hexdigest(),
         'finite receipt mismatch')
    cur=0
    defects=[]
    finite_z=[]
    for a,b in FINITE_RANGES:
        item=get(f'finite_r14_{a:03d}_{b:03d}.json')
        need(a==cur and item['kind']=='finite_batch' and item['dimension']==14
             and item['start']==a and item['stop']==b and item['N']==610,
             'finite-column coverage gap')
        cur=b
        defects.append(interval(item['inverse_defect_max']))
        finite_z.append(interval(item['Z_finite_max']))
    need(cur==610,'finite columns incomplete')
    Anfin=interval(finite['A_finite']).upper()
    inv_def=v.max_upper(defects)
    Y=interval(finite['Y']).upper()
    Zfin=v.max_upper(finite_z)
    need(inv_def<=interval(finite['inverse_defect']).upper()
         and Zfin<=interval(finite['Z_finite']).upper()
         and Anfin<arb(11767) and inv_def<arb(1),
         'finite inverse or aggregate failed')
    tail=get('tail_inverse_bound_r14_M114_S60.json')
    need(tail['status']=='OPEN_COMPONENT' and tail['dimension']==14
         and tail['M']==114 and tail['S']==60 and tail['count']==100
         and tail['far_start']==1326,'principal tail receipt mismatch')
    Atail=interval(tail['A_tail_upper']).upper()
    need(Atail<arb(1662)/1000 and Anfin>Atail,
         'principal tail inverse or full-norm comparison failed')
    gfar=get('g_far_bound_r14_M114_S60.json')
    need(gfar['kind']=='far' and gfar['dimension']==14
         and gfar['degree']==60 and gfar['near_count']==1215
         and gfar['n_far']==186 and gfar['s_far']==122,
         'g-far support mismatch')
    cur=0
    gnear=[]
    for a,b in G_RANGES:
        item=get(f'g_near_r14_{a:04d}_{b:04d}.json')
        need(a==cur and item['kind']=='near' and item['dimension']==14
             and item['start']==a and item['stop']==b
             and item['near_count']==1215,'g-near coverage gap')
        cur=b
        gnear.append(interval(item['max']))
    need(cur==1215,'g-near incomplete')
    Zgnear=v.max_upper(gnear)
    Zgfar=v.max_upper([interval(gfar['far_n']),interval(gfar['far_s'])])
    cur=121
    shnear=[]
    for a,b,count in SHAPE_RANGES:
        item=get(f'shape_near_r14_{a:03d}_{b:03d}.json')
        need(a==cur and item['kind']=='near' and item['dimension']==14
             and item['first']==a and item['last']==b and item['count']==count,
             'shape-near coverage gap')
        cur=b+12
        shnear.append(interval(item['max']))
    need(cur==241,'shape-near incomplete')
    Zshnear=v.max_upper(shnear)
    shfar=get('far_shape_split_r14_j241.json')
    need(shfar['status']=='OPEN_COMPONENT' and shfar['dimension']==14
         and shfar['jstar']==241 and shfar['M']==114 and shfar['L']==222
         and interval(shfar['Q_identity_error'])<arb('1e-20'),
         'far-shape split support or Q identity failed')
    Zshfar=interval(shfar['Z_shape_far']).upper()
    identities=get('identities_r14_M114_S60.json')
    need(identities['status']=='IDENTITIES_PASS' and identities['dimension']==14
         and identities['bits']>=256 and identities['K_cases']==12
         and identities['radial_power_cases']==3
         and interval(identities['ball_principal_max_error'])<arb('1e-50')
         and interval(identities['Q_connection_max_error'])<arb('1e-50'),
         'identity receipt failed')
    Zbound=v.max_upper([Zfin,Zgnear,Zgfar,Zshnear,Zshfar])

    # H=g+|psi'|^2(Kg+Im(psi^6)/B), degree eight.
    m=v.M_ROOT
    b=c[1]
    B=v.BOUNDARY_SCALE
    alpha=b**(m+1)/B
    P=sum((abs(x)*rho**a for a,x in p.items()),v.Z)
    C=sum((abs(x)*rho**j for j,x in c.items()),v.Z)
    Vn=v.norm(V,rho)
    theta0=1/(alpha*(m+1)*rho**(m-1))
    theta1=1/(alpha*rho**m)
    kap=arb(1)/((m+2)*(m+4))  # rigorous Arb enclosure of 1/80 (not a binary64 float)
    Anq=arb(11767)
    h={k:arb(math.comb(m,k))*C**(m-k)*theta0**k/B
       for k in range(1,m+1)}
    nonlinear={}
    for q in range(2,m+3):
        val=arb(0)
        if q in h:val+=P*P*h[q]
        if q-1 in h:val+=2*P*theta1*h[q-1]
        if q-2 in h:val+=theta1**2*h[q-2]
        if q==2:val+=2*P*theta1*kap+theta1**2*Vn
        if q==3:val+=theta1**2*kap
        nonlinear[q]=Anq*val
    Yq=arb(46)/1000000
    Zq=arb(61)/100
    rational={2:arb(211),3:arb(61)/10000,
              4:arb(2)/10**9,5:arb(3)/10**14,
              6:arb(2)/10**19,7:arb(7)/10**25,
              8:arb(2)/10**30}
    need(Y<Yq and Zbound<Zq
         and all(nonlinear[q]<rational[q] for q in rational),
         'R14 majorant exceeded its rational threshold')
    radius=arb(1)/6000
    polynomial=Yq+(Zq-1)*radius+sum((rational[q]*radius**q for q in rational),v.Z)
    derivative=Zq-1+sum((q*rational[q]*radius**(q-1) for q in rational),v.Z)
    need(polynomial<0 and derivative<0,'R14 radii polynomial failed')

    ext=arb(21)/20
    need(ext<rho,'extension radius exceeds coefficient radius')
    extension_margin=b-sum((j*abs(x)*ext**(j-1) for j,x in c.items() if j>1),v.Z)-radius*theta1
    derivative_lower=b-sum((j*abs(x) for j,x in c.items() if j>1),v.Z)-radius*theta1
    second_perturb=radius/(alpha*rho**(m-1)*arb.const_e()*rho.log())
    second_upper=sum((j*(j-1)*abs(x) for j,x in c.items() if j>1),v.Z)+second_perturb
    curvature_lower=1-second_upper/derivative_lower
    jfirst=1+v.STEP
    cfirst_lower=abs(c[jfirst])-radius/(alpha*(jfirst+m)*rho**(jfirst+m-1))
    need(b>0 and extension_margin>0 and derivative_lower>0
         and curvature_lower>0 and cfirst_lower>0,
         'R14 analytic convex non-disc geometry failed')
    out={'status':'PROVED','dimension':14,'bits':ctx.prec,
         'M':114,'S':60,'rho_weight':'11/10',
         'Y':Y.str(30),'Z':Zbound.str(30),
         'Z_finite':Zfin.str(30),'Z_g_near':Zgnear.str(30),
         'Z_g_far':Zgfar.str(30),'Z_shape_near':Zshnear.str(30),
         'Z_shape_far':Zshfar.str(30),
         'A_finite':Anfin.str(30),'A_tail':Atail.str(30),
         'inverse_defect':inv_def.str(30),
         'C2_to_C8':{str(q):nonlinear[q].str(30) for q in nonlinear},
         'radius':radius.str(30),'radii_polynomial':polynomial.str(30),
         'derivative_polynomial':derivative.str(30),
         'extension_margin':extension_margin.str(30),
         'curvature_lower':curvature_lower.str(30),
         'c13_lower':cfirst_lower.str(30),
         'frozen_bundle_sha256':FROZEN_BUNDLE_SHA256,
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path('certificate_r14.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PROVED')
    print('certificate: certificate_r14.json')


if __name__=='__main__':
    main()
