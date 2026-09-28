"""Independent Arb convexity certificate for the R14 Schiffer domain."""

import hashlib
import json
from pathlib import Path
from flint import arb, ctx


def need(ok,why):
    if not ok:
        raise RuntimeError(why)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dyadic(s):
    a,b=float.fromhex(s).as_integer_ratio()
    x=arb(a)/b
    need(x.is_exact(),'centre coefficient not exact dyadic')
    return x


def main():
    ctx.prec=128
    ctx.threads=1
    cert=json.loads(Path('certificate_r14.json').read_text())
    need(cert['status']=='PROVED' and cert['dimension']==14 and cert['bits']>=128
         and cert['verifier_sha256']==sha('verify_r14.py'),
         'main R14 certificate or verifier mismatch')
    d=json.loads(Path('center_r14_M114_S60.json').read_text())
    need(d['m']==6 and d['step']==12 and d['M']==114 and d['S']==60,
         'wrong R14 centre')
    js=list(range(1,1+12*len(d['c']),12))
    c={j:dyadic(s) for j,s in zip(js,d['c'])}
    b=c[1]
    B=dyadic(d['boundary_scale'])
    rho=arb(11)/10
    alpha=b**7/B
    radius=arb(1)/6000
    need(arb(cert['radius']).contains(radius),'main radius mismatch')
    theta1=1/(alpha*rho**6)
    first=sum((j*abs(x) for j,x in c.items() if j>1),arb(0))
    second=sum((j*(j-1)*abs(x) for j,x in c.items() if j>1),arb(0))
    derivative_lower=b-first-radius*theta1
    second_perturb=radius/(alpha*rho**5*arb.const_e()*rho.log())
    second_upper=second+second_perturb
    curvature_lower=1-second_upper/derivative_lower
    ext=arb(21)/20
    extension_margin=b-sum((j*abs(x)*ext**(j-1) for j,x in c.items() if j>1),arb(0))-radius*theta1
    omega13=alpha*(13+6)*rho**18
    c13_lower=abs(c[13])-radius/omega13
    need(derivative_lower>0 and curvature_lower>arb(89)/100
         and extension_margin>0 and c13_lower>0,
         'R14 convexity or non-disc inequality failed')
    out={'status':'PROVED_CONVEX','dimension':14,'bits':ctx.prec,
         'main_certificate_sha256':sha('certificate_r14.json'),
         'main_verifier_sha256':sha('verify_r14.py'),
         'centre_sha256':sha('center_r14_M114_S60.json'),
         'radius':radius.str(30),
         'extension_margin':extension_margin.str(30),
         'derivative_lower':derivative_lower.str(30),
         'second_upper':second_upper.str(30),
         'curvature_lower':curvature_lower.str(30),
         'c13_lower':c13_lower.str(30)}
    Path('convex_r14.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PROVED_CONVEX',curvature_lower.str(20))


if __name__=='__main__':
    main()
