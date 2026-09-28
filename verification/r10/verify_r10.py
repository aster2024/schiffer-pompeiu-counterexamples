#!/usr/bin/env python3
"""Fail-closed Arb CAP verifier for the D8-symmetric Schiffer domain in R10."""

import hashlib
import json
import math
from pathlib import Path
from flint import arb, ctx
import rank2_cap_core as v

FROZEN_BUNDLE_SHA256='0c8b38d38b96138ae2d59cee3d35b3bf1eecdd9c737a0a6bf39c6e1e70b66884'
FROZEN_FILES=[
    'rank2_cap_core.py','center_r10_M52_S36.json','inverse_r10_M52_S36.npy',
    'finite_rank2.py','finite_bound_r10.json',
    'tail_inverse_rank2.py','tail_inverse_bound_r10.json',
    'g_tail_rank2.py','g_far_bound_r10.json',
    'shape_tail_rank2.py','shape_near_r10_057_169.json',
    'shape_near_r10_177_217.json','shape_near_r10_225_265.json',
    'shape_near_r10_273_297.json',
    'far_shape_rank2.py','far_shape_bound_r10.json',
    'check_identities_rank2.py','identities_r10.json',
]+[f'g_near_r10_{a:03d}_{min(a+20,358):03d}.json'
   for a in range(0,358,20)]


def digest():
    h=hashlib.sha256()
    for name in FROZEN_FILES:
        b=Path(name).read_bytes()
        h.update(name.encode()+b'\0')
        h.update(len(b).to_bytes(8,'big'))
        h.update(b)
    return h.hexdigest()


def need(ok,reason):
    if not ok:
        raise RuntimeError(reason)


def get(name):
    return json.loads(Path(name).read_text())


def interval(value):
    x=arb(value)
    need(x.is_finite() and x.upper().is_finite(),'nonfinite interval in receipt')
    return x


def upper_max(xs):
    return v.max_upper(xs)


def main():
    ctx.prec=128
    ctx.threads=1
    need(len(FROZEN_BUNDLE_SHA256)==64 and digest()==FROZEN_BUNDLE_SHA256,
         'frozen R10 bundle missing or changed')
    d=get('center_r10_M52_S36.json')
    rows,g,c,p,V,F=v.make_algebra(d)
    need(d['m']==4 and d['step']==8 and d['M']==52 and d['S']==36
         and len(rows)==259,'unsupported R10 centre')
    rho=arb(11)/10
    finite=get('finite_bound_r10.json')
    need(finite['status']=='OPEN_COMPONENT' and finite['dimension']==10
         and finite['N']==259 and finite['bits']>=128
         and finite['centre_sha256']==hashlib.sha256(Path('center_r10_M52_S36.json').read_bytes()).hexdigest()
         and finite['inverse_sha256']==hashlib.sha256(Path('inverse_r10_M52_S36.npy').read_bytes()).hexdigest(),
         'finite receipt mismatch')
    Anfin=interval(finite['A_finite']).upper()
    inv_def=interval(finite['inverse_defect']).upper()
    Y=interval(finite['Y']).upper()
    Zfin=interval(finite['Z_finite']).upper()
    need(Anfin<arb(10195) and inv_def<arb(1),'finite inverse failed')
    tail=get('tail_inverse_bound_r10.json')
    need(tail['status']=='OPEN_COMPONENT' and tail['dimension']==10
         and tail['M']==52 and tail['S']==36 and tail['count']==100
         and tail['far_start']==860,'tail inverse receipt mismatch')
    Atail=interval(tail['A_tail_upper']).upper()
    need(Atail<arb(2015)/1000 and Anfin>Atail,
         'principal tail inverse or full-norm comparison failed')
    gfar=get('g_far_bound_r10.json')
    need(gfar['kind']=='far' and gfar['dimension']==10
         and gfar['degree']==24 and gfar['near_count']==358
         and gfar['n_far']==84 and gfar['s_far']==62,
         'g-far support mismatch')
    cur=0
    gnear=[]
    for a in range(0,358,20):
        stop=min(a+20,358)
        item=get(f'g_near_r10_{a:03d}_{stop:03d}.json')
        need(item['kind']=='near' and item['dimension']==10
             and item['start']==cur and item['stop']==stop
             and item['near_count']==358,'g-near coverage gap')
        cur=stop
        gnear.append(interval(item['max']))
    need(cur==358,'g-near incomplete')
    Zgnear=upper_max(gnear)
    Zgfar=upper_max([interval(gfar['far_n']),interval(gfar['far_s'])])
    shapes=[(57,169,15),(177,217,6),(225,265,6),(273,297,4)]
    shape_bounds=[]
    for first,last,count in shapes:
        item=get(f'shape_near_r10_{first:03d}_{last:03d}.json')
        need(item['kind']=='near' and item['dimension']==10
             and item['first']==first and item['last']==last
             and item['count']==count,'shape-near coverage gap')
        shape_bounds.append(interval(item['max']))
    Zshnear=upper_max(shape_bounds)
    shfar=get('far_shape_bound_r10.json')
    need(shfar['status']=='OPEN_COMPONENT' and shfar['dimension']==10
         and shfar['jstar']==305 and shfar['M']==52
         and shfar['min_freq']>52
         and interval(shfar['Q_identity_error'])<arb('1e-20'),
         'shape-far support or Q identity failed')
    Zshfar=interval(shfar['Z_shape_far']).upper()
    identities=get('identities_r10.json')
    need(identities['status']=='IDENTITIES_PASS'
         and identities['dimension']==10 and identities['bits']>=256
         and identities['K_cases']==12 and identities['radial_power_cases']==3
         and interval(identities['ball_principal_max_error'])<arb('1e-50')
         and interval(identities['Q_connection_max_error'])<arb('1e-50'),
         'identity receipt failed')
    Zbound=upper_max([Zfin,Zgnear,Zgfar,Zshnear,Zshfar])

    # Bound all degrees 2,...,6 of g+|psi'|^2(Kg+Im(psi^4)/B).
    m=v.M_ROOT
    b=c[1]
    B=v.BOUNDARY_SCALE
    alpha=b**(m+1)/B
    P=sum((abs(x)*rho**a for a,x in p.items()),v.Z)
    C=sum((abs(x)*rho**j for j,x in c.items()),v.Z)
    Vn=v.norm(V,rho)
    theta0=1/(alpha*(m+1)*rho**(m-1))
    theta1=1/(alpha*rho**m)
    kap=arb(1)/((m+2)*(m+4))  # rigorous Arb enclosure of 1/48 (not a binary64 float)
    Anq=arb(10195)
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
    Yq=arb(1)/100000000
    Zq=arb(56)/100
    rational={2:arb(396),3:arb(72)/1000,
              4:arb(4)/10000000,5:arb(2)/100000000000,
              6:arb(4)/10000000000000000}
    need(Y<Yq and Zbound<Zq and all(nonlinear[q]<rational[q] for q in rational),
         'R10 majorant exceeded its rational threshold')
    radius=arb(1)/10000000
    polynomial=Yq+(Zq-1)*radius+sum((rational[q]*radius**q for q in rational),v.Z)
    derivative=Zq-1+sum((q*rational[q]*radius**(q-1) for q in rational),v.Z)
    need(polynomial<0 and derivative<0,'R10 radii polynomial failed')

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
         'R10 analytic convex non-disc geometry failed')
    out={'status':'PROVED','dimension':10,'bits':ctx.prec,
         'M':52,'S':36,'rho_weight':'11/10',
         'Y':Y.str(30),'Z':Zbound.str(30),
         'Z_finite':Zfin.str(30),'Z_g_near':Zgnear.str(30),
         'Z_g_far':Zgfar.str(30),'Z_shape_near':Zshnear.str(30),
         'Z_shape_far':Zshfar.str(30),
         'A_finite':Anfin.str(30),'A_tail':Atail.str(30),
         'inverse_defect':inv_def.str(30),
         'C2_to_C6':{str(q):nonlinear[q].str(30) for q in nonlinear},
         'radius':radius.str(30),'radii_polynomial':polynomial.str(30),
         'derivative_polynomial':derivative.str(30),
         'extension_margin':extension_margin.str(30),
         'curvature_lower':curvature_lower.str(30),
         'c9_lower':cfirst_lower.str(30),
         'frozen_bundle_sha256':FROZEN_BUNDLE_SHA256,
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path('certificate_r10.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PROVED')
    print('certificate: certificate_r10.json')


if __name__=='__main__':
    main()
