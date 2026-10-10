"""Arb a posteriori bound for the sole R14 finite column 379 left by the
coarse estimate. Uses a frozen dyadic Schur seed plus the global K inverse.

This completes the finite-column class for the stated 13-step triangular
preconditioner, but by itself does not prove nonlinear existence.
"""
from proof_guard import require
import hashlib
import json
import numpy as np
from flint import arb,arb_mat,ctx
from exact_finite_rank2 import center,columns,linear_action,q_of,to_char
from exact_finite_r4 import ZERO,conv,scale,dyadic


def finite_upper(x):
    require((x.is_finite()), 'verify_last_finite_column_r14.py:16: proof gate failed')
    y=x.upper();require((y.is_finite()), 'verify_last_finite_column_r14.py:17: proof gate failed')
    return y


def main():
    ctx.prec=192;ctx.threads=1
    name='center_r14_M114_S36_arbrefined.json'
    seedname='r14_M114_lastshape_tail_seed.json'
    invname='Ahat_r14_M114_S36.npy'
    d,g,c,B,keys,ms,js=center(name)
    require(((d['m'],d['M'],d['S'])==(6,114,36)), 'verify_last_finite_column_r14.py:27: proof gate failed')
    m=d['m'];rho=arb(102)/100;b=c[1]
    with open(seedname) as f:seed=json.load(f)
    require((seed['center_sha256']==hashlib.sha256(open(name,'rb').read()).hexdigest()), 'verify_last_finite_column_r14.py:30: proof gate failed')
    i=379;require((seed['finite_index']==i and seed['m']==m), 'verify_last_finite_column_r14.py:31: proof gate failed')
    require((seed['jlist']==[js[-1]+2*m*(k+1) for k in range(8)]), 'verify_last_finite_column_r14.py:32: proof gate failed')
    chat={j:dyadic(s) for j,s in zip(seed['jlist'],seed['c_hex'])}
    F,G=columns(g,c,B,m,keys,js,i,i+1)[0]
    q0=q_of(g,m)
    QF=to_char(scale(conv(q0,q_of(F,m)),2/b),m)
    bprime={n:G.get(n,ZERO)-QF.get(n,ZERO) for n in set(G)|set(QF)}
    HF,HG=linear_action(g,c,B,m,{},chat)
    QH=to_char(scale(conv(q0,q_of(HF,m)),2/b),m)
    Cchat={n:HG.get(n,ZERO)-QH.get(n,ZERO) for n in set(HG)|set(QH)}
    R={n:bprime.get(n,ZERO)-Cchat.get(n,ZERO)
       for n in set(bprime)|set(Cchat) if n not in ms}
    Rnorm=sum((rho**(n-m)*abs(v) for n,v in R.items()),ZERO)
    hb=b**m/B;alpha=hb*hb
    c0={};run=arb(0)
    if bprime:
        for j in range(max(bprime)-m+1,js[-1],-2*m):
            run+=bprime.get(j+m-1,ZERO)/(alpha*(j+2*m))
            c0[j]=-run
    c0norm=sum(((j+2*m)*rho**(j-1)*abs(v) for j,v in c0.items()),ZERO)
    k=arb(474)/1000;Cinv=arb(9);steps=13
    cerr=Cinv*Rnorm+k**steps/(1-k)*c0norm
    ghat={key:F.get(key,ZERO)-HF.get(key,ZERO)
          for key in set(F)|set(HF) if key not in keys}
    JF,JG=linear_action(g,c,B,m,ghat,chat)
    y=arb_mat([[F.get(key,ZERO)-JF.get(key,ZERO)] for key in keys]+
              [[G.get(n,ZERO)-JG.get(n,ZERO)] for n in ms])
    aa=np.load(invname,mmap_mode='r')
    n=len(keys)+len(js);require((aa.shape==(n,n) and np.isfinite(aa).all()), 'verify_last_finite_column_r14.py:59: proof gate failed')
    A=arb_mat([[dyadic(float(aa[r,s]).hex()) for s in range(n)] for r in range(n)])
    x=A*y
    win=[rho**nn for nn,s in keys]+[(j+2*m)*rho**(j-1) for j in js]
    w=win[i]
    finite=sum((win[t]*abs(x[t,0]-(arb(1) if t==i else ZERO))
                for t in range(n)),ZERO)/w
    highg=sum((rho**nn*abs(v) for (nn,s),v in ghat.items()),ZERO)/w
    highc=sum(((j+2*m)*rho**(j-1)*abs(v) for j,v in chat.items()),ZERO)/w
    approx=finite+highg+highc
    difference=((1+arb(1)/5)*250+11)*cerr/w
    total=approx+difference
    require((finite_upper(total)<arb(1)/4), 'verify_last_finite_column_r14.py:71: proof gate failed')
    print('LAST_FINITE_COLUMN_BOUND_ONLY','index',i,
          'Schur residual',finite_upper(Rnorm),
          'c0 norm',finite_upper(c0norm),
          'tail inverse error',finite_upper(cerr/w),
          'finite defect',finite_upper(finite),
          'high source defect',finite_upper(highg),
          'high shape defect',finite_upper(highc),
          'remainder',finite_upper(difference),
          'whole column upper',finite_upper(total))


if __name__=='__main__':
    main()
