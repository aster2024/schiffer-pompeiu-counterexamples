"""Independent Arb convexity check for the certified R8 conformal map."""

import argparse
import hashlib
import json
from pathlib import Path
from flint import arb, ctx


def need(test,why):
    if not test:
        raise RuntimeError(why)


def exact_dyadic(text):
    a,b=float.fromhex(text).as_integer_ratio()
    x=arb(a)/b
    need(x.is_exact(),'input centre is not an exact dyadic')
    return x


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--main-certificate',default='certificate_r8.json')
    p.add_argument('--output',default='convex_r8.json')
    p.add_argument('--bits',type=int,default=128)
    args=p.parse_args()
    need(args.bits>=128,'insufficient Arb precision')
    ctx.prec=args.bits
    ctx.threads=1
    cert=json.loads(Path(args.main_certificate).read_text())
    need(cert['status']=='PROVED' and cert['dimension']==8 and cert['bits']>=128,
         'main R8 existence certificate missing')
    need(cert['verifier_sha256']==sha('verify_r8.py'),
         'main verifier changed after certificate')
    data=json.loads(Path('center_r8_M45_S24.json').read_text())
    need(data['m']==3 and data['step']==6 and data['M']==45 and data['S']==24,
         'unexpected centre dimensions')
    js=list(range(1,44,6))
    need(len(js)==len(data['c'])==8,'shape coefficient count mismatch')
    c={j:exact_dyadic(s) for j,s in zip(js,data['c'])}
    rho=arb(11)/10
    b=c[1]
    alpha=b**4
    r=arb(1)/30000
    need(arb(cert['radius']).contains(r),'main certificate radius mismatch')
    omega=lambda j:alpha*(j+3)*rho**(j+2)
    theta1=1/(alpha*rho**3)
    first=sum((j*abs(x) for j,x in c.items() if j>=7),arb(0))
    second=sum((j*(j-1)*abs(x) for j,x in c.items() if j>=7),arb(0))
    derivative_lower=b-first-r*theta1
    second_perturb=r/(alpha*rho**2*arb.const_e()*rho.log())
    second_upper=second+second_perturb
    convex_lower=1-second_upper/derivative_lower
    ext=arb(21)/20
    extension_lower=b-sum((j*abs(x)*ext**(j-1) for j,x in c.items() if j>=7),arb(0))-r*theta1
    c7_lower=abs(c[7])-r/omega(7)
    need(derivative_lower>0 and convex_lower>0 and extension_lower>0
         and c7_lower>0,'convexity or non-disc inequality failed')
    need(convex_lower>arb(3)/4,'convexity margin below the released threshold')
    out={'status':'PROVED_CONVEX','bits':args.bits,
         'main_certificate_sha256':sha(args.main_certificate),
         'main_verifier_sha256':sha('verify_r8.py'),
         'centre_sha256':sha('center_r8_M45_S24.json'),
         'radius':r.str(30),'extension_margin':extension_lower.str(30),
         'derivative_lower':derivative_lower.str(30),
         'second_upper':second_upper.str(30),
         'curvature_lower':convex_lower.str(30),
         'c7_lower':c7_lower.str(30)}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print('PROVED_CONVEX',convex_lower.str(20))


if __name__=='__main__':
    main()
