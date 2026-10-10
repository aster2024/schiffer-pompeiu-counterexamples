"""Exact finite coefficient algebra and Arb residual for I2(m) Berenstein.

This is a center diagnostic, not a proof of an infinite-dimensional zero.
"""
from proof_guard import require
import argparse
import hashlib
import json
from fractions import Fraction
from flint import arb,ctx
from exact_finite_r4 import ZERO,acc,add,scale,poly,abs2,conv,dyadic


def sine(f,m):
    out={}
    for (n,s),v in f.items():
        require((n>=m and (n-m)%(2*m)==0), 'exact_finite_rank2.py:16: proof gate failed')
        acc(out,(n,s),-v/2)
        acc(out,(-n,s),v/2)
    return out


def unsine(f,m):
    out={}
    for (n,s),v in f.items():
        require(((abs(n)-m)%(2*m)==0), 'exact_finite_rank2.py:25: proof gate failed')
        acc(out,(abs(n),s),-v if n>0 else v)
    return out


def KD(f,m):
    out={}
    for (n,s),v in f.items():
        require((n>=m and (n-m)%(2*m)==0 and s>=0), 'exact_finite_rank2.py:33: proof gate failed')
        if s==0:
            a=v/(4*(n+1)*(n+2))
            acc(out,(n,0),-a);acc(out,(n,1),a)
        else:
            d=n+2*s
            acc(out,(n,s-1),v/(4*d*(d+1)))
            acc(out,(n,s),-v/(2*d*(d+2)))
            acc(out,(n,s+1),v/(4*(d+1)*(d+2)))
    return out


def char(n,m):
    require((n>=m and (n-m)%(2*m)==0), 'exact_finite_rank2.py:46: proof gate failed')
    return {a:arb(1) for a in range(-(n-m),n-m+1,2*m)}


def power(c,k):
    out={0:arb(1)}
    for _ in range(k):
        out=conv(out,c)
    return out


def h_of(c,m,B):
    factor=(arb(1)/2 if m==2 else arb(1))/B
    return add(*(scale(char(n,m),factor*v) for n,v in power(c,m).items()))


def dh_of(c,j,m,B):
    factor=(arb(1)/2 if m==2 else arb(1))*m/B
    return add(*(scale(char(n+j,m),factor*v) for n,v in power(c,m-1).items()))


def q_of(g,m):
    return add(*(scale(char(n,m),v/(2*(n+1)))
                 for (n,s),v in g.items() if s==0))


def a2_of(c):
    p={j-1:j*v for j,v in c.items()}
    return conv(p,{-a:v for a,v in p.items()})


def to_char(f,m):
    require((all(a%(2*m)==0 for a in f)), 'exact_finite_rank2.py:78: proof gate failed')
    amax=max((abs(a) for a in f),default=0)
    return {a+m:f.get(a,ZERO)-f.get(a+2*m,ZERO)
            for a in range(0,amax+1,2*m)}


def center(path):
    with open(path) as f:d=json.load(f)
    m,M,S=d['m'],d['M'],d['S']
    ms=list(range(m,M+1,2*m))
    js=list(range(1,1+2*m*len(ms),2*m))
    keys=[(n,s) for n in ms for s in range(S+1)]
    if 'g_dyadic' in d:
        require((len(keys)==len(d['g_dyadic']) and len(js)==len(d['c_dyadic'])), 'exact_finite_rank2.py:91: proof gate failed')
        g={k:arb(man)*arb(2)**exp for k,(man,exp) in zip(keys,d['g_dyadic'])}
        c={j:arb(man)*arb(2)**exp for j,(man,exp) in zip(js,d['c_dyadic'])}
    else:
        require((len(keys)==len(d['g_hex']) and len(js)==len(d['c_hex'])), 'exact_finite_rank2.py:95: proof gate failed')
        g={k:dyadic(v) for k,v in zip(keys,d['g_hex'])}
        c={j:dyadic(v) for j,v in zip(js,d['c_hex'])}
    B=dyadic(d['B_hex'])
    require((B>0), 'exact_finite_rank2.py:99: proof gate failed')
    return d,g,c,B,keys,ms,js


def evaluate(g,c,B,m):
    p={j-1:j*v for j,v in c.items()}
    F=add(g,unsine(abs2(sine(KD(g,m),m),p),m))
    q=q_of(g,m);h=h_of(c,m,B);aa=a2_of(c)
    b=c[1]
    boundary=to_char(add(conv(q,q),scale(conv(conv(h,h),aa),-1)),m)
    return F,scale(boundary,1/b)


def columns(g,c,B,m,keys,js,start,end):
    p={j-1:j*v for j,v in c.items()}
    q=q_of(g,m);h=h_of(c,m,B);aa=a2_of(c)
    V=sine(KD(g,m),m)
    out=[]
    for idx in range(start,end):
        if idx<len(keys):
            key=keys[idx];e={key:arb(1)}
            Fi=add(e,unsine(abs2(sine(KD(e,m),m),p),m))
            if key[1]==0:
                dq=scale(char(key[0],m),arb(1)/(2*(key[0]+1)))
                Bi=to_char(scale(conv(q,dq),2),m)
            else:Bi={}
        else:
            j=js[idx-len(keys)]
            dp={j-1:arb(j)}
            mixed=add(poly(poly(V,p,True),dp),poly(poly(V,dp,True),p))
            Fi=unsine(mixed,m)
            dh=dh_of(c,j,m,B)
            da2=add(conv(dp,{-k:v for k,v in p.items()}),
                    conv(p,{-k:v for k,v in dp.items()}))
            Bi=to_char(add(scale(conv(conv(h,dh),aa),-2),
                            scale(conv(conv(h,h),da2),-1)),m)
        out.append((Fi,scale(Bi,1/c[1])))
    return out


def linear_action(g,c,B,m,dg,dc):
    """Exact derivative at the frozen center applied to arbitrary finite increments."""
    p={j-1:j*v for j,v in c.items()}
    dp={j-1:j*v for j,v in dc.items()}
    Fi=add(dg,unsine(abs2(sine(KD(dg,m),m),p),m))
    q=q_of(g,m);h=h_of(c,m,B);aa=a2_of(c)
    if dp:
        V=sine(KD(g,m),m)
        Fi=add(Fi,unsine(add(poly(poly(V,p,True),dp),
                              poly(poly(V,dp,True),p)),m))
        dh=add(*(scale(dh_of(c,j,m,B),v) for j,v in dc.items()))
        da2=add(conv(dp,{-k:v for k,v in p.items()}),
                conv(p,{-k:v for k,v in dp.items()}))
        Bs=add(scale(conv(conv(h,dh),aa),-2),
               scale(conv(conv(h,h),da2),-1))
    else:
        Bs={}
    Bq=scale(conv(q,q_of(dg,m)),2)
    Bi=to_char(add(Bq,Bs),m)
    return Fi,scale(Bi,1/c[1])


def norm(F,G,m,rho):
    return (sum((rho**n*abs(v) for (n,s),v in F.items()),ZERO)
            +sum((rho**(n-m)*abs(v) for n,v in G.items()),ZERO))


def finite_upper(v):
    require((v.is_finite()), 'non-finite Arb enclosure')
    y=v.upper()
    require((y.is_finite()), 'non-finite Arb upper endpoint')
    return y


def maximum_abs_upper(values):
    best=arb(0)
    for v in values:
        x=finite_upper(abs(v))
        if x>best:best=x
    return best


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--center',required=True)
    ap.add_argument('--bits',type=int,default=160)
    ap.add_argument('--batch',default='')
    ap.add_argument('--output',default='')
    ap.add_argument('--check-rho-num',type=int,default=0)
    ap.add_argument('--check-rho-den',type=int,default=1)
    ap.add_argument('--max-residual',default='')
    args=ap.parse_args()
    ctx.prec=args.bits;ctx.threads=1
    d,g,c,B,keys,ms,js=center(args.center)
    F,G=evaluate(g,c,B,d['m'])
    print('m',d['m'],'dimension',len(keys)+len(js),'supports',len(F),len(G))
    print('finite interior max',maximum_abs_upper(F.get(k,ZERO) for k in keys))
    print('finite boundary max',maximum_abs_upper(G.get(n,ZERO) for n in ms))
    for p,q in ((101,100),(102,100),(105,100),(11,10)):
        rho=arb(p)/q
        print('rho',p,'/',q,'full residual upper',finite_upper(norm(F,G,d['m'],rho)))
    if args.max_residual:
        require((args.check_rho_num>args.check_rho_den>0), 'exact_finite_rank2.py:201: proof gate failed')
        ceiling=Fraction(args.max_residual)
        require((ceiling>0), 'exact_finite_rank2.py:203: proof gate failed')
        check=norm(F,G,d['m'],arb(args.check_rho_num)/args.check_rho_den)
        require((finite_upper(check)<arb(ceiling.numerator)/ceiling.denominator), 'exact_finite_rank2.py:205: proof gate failed')
        print('RESIDUAL_THRESHOLD_PASS',args.check_rho_num,args.check_rho_den,
              args.max_residual)
    if args.batch:
        start,end=map(int,args.batch.split(':'))
        require((0<=start<end<=len(keys)+len(js) and args.output), 'exact_finite_rank2.py:210: proof gate failed')
        entries=[]
        for col,(Fi,Gi) in enumerate(columns(g,c,B,d['m'],keys,js,start,end),start):
            for row,k in enumerate(keys):
                v=Fi.get(k,ZERO)
                if not v.is_zero():
                    a,e=v.mid().man_exp();r,f=v.rad().man_exp()
                    entries.append([row,col,int(a),int(e),int(r),int(f)])
            for row,n in enumerate(ms,len(keys)):
                v=Gi.get(n,ZERO)
                if not v.is_zero():
                    a,e=v.mid().man_exp();r,f=v.rad().man_exp()
                    entries.append([row,col,int(a),int(e),int(r),int(f)])
        payload=dict(center=args.center,
                     center_sha256=hashlib.sha256(open(args.center,'rb').read()).hexdigest(),
                     bits=args.bits,start=start,end=end,dimension=len(keys)+len(js),
                     entries=entries)
        with open(args.output,'w') as f:json.dump(payload,f,separators=(',',':'))
        print('wrote batch',start,end,'entries',len(entries),args.output)


if __name__=='__main__':
    main()
