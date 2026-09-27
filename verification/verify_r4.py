#!/usr/bin/env python3
"""Standalone Arb certificate for the R4 Schiffer/Pompeiu construction.

Run: nice -n 19 python3 verify_r4.py --output certificate_r4.json
Requires python-flint (tested with 0.9.0); no external data or network.
The proof and all formulas are in FINAL.md. Partial modes never print PROVED.
"""

def require(condition, message="verification condition failed"):
    if not condition:
        raise RuntimeError(message)

from flint import arb, arb_mat

ZERO = arb(0)


def acc(f, k, v):
    if not v.is_zero():
        f[k] = f.get(k, ZERO) + v


def add(*fs):
    out = {}
    for f in fs:
        for k, v in f.items():
            acc(out, k, v)
    return out


def scale(f, a):
    return {k: a*v for k, v in f.items() if not v.is_zero()}


def z(f, conjugate=False):
    out = {}
    sign = -1 if conjugate else 1
    for (m, s), v in f.items():
        n = abs(m)
        d = n + 2*s + 1
        if sign*m >= 0:
            acc(out, (m+sign,s), v*(n+s+1)/d)
            if s:
                acc(out, (m+sign,s-1), v*s/d)
        else:
            acc(out, (m+sign,s), v*(n+s)/d)
            acc(out, (m+sign,s+1), v*(s+1)/d)
    return out


def poly(f, p, conjugate=False):
    """p maps nonnegative even exponents to Arb coefficients."""
    degree = max(p)
    require(degree % 2 == 0 and all((a % 2 == 0 for a in p)), 'Failed: degree % 2 == 0 and all((a % 2 == 0 for a in p))')
    out = scale(f,p[degree])
    for a in range(degree-2,-1,-2):
        out = z(z(out,conjugate),conjugate)
        if a in p:
            out = add(out,scale(f,p[a]))
    return out


def abs2(f,p):
    return poly(poly(f,p,True),p)


def sine(f):
    out = {}
    for (n,s),v in f.items():
        require(n > 0 and n % 2 == 1, 'Failed: n > 0 and n % 2 == 1')
        acc(out,(n,s),-v/2)
        acc(out,(-n,s),v/2)
    return out


def unsine(f):
    """Antisymmetric projection; exact for all real sine fields used here."""
    out = {}
    for (m,s),v in f.items():
        require(m % 2 != 0, (m, s))
        acc(out,(abs(m),s),-v if m>0 else v)
    return out


def K(f):
    out = {}
    for (n,s),v in f.items():
        require(n > 0 and n % 2 == 1 and (s >= 1), 'Failed: n > 0 and n % 2 == 1 and (s >= 1)')
        d=n+2*s
        acc(out,(n,s-1),v/(4*d*(d+1)))
        acc(out,(n,s),-v/(2*d*(d+2)))
        acc(out,(n,s+1),v/(4*(d+1)*(d+2)))
    return out


def norm(f,rho):
    return sum((abs(v)*rho**abs(m) for (m,s),v in f.items()),ZERO)


def upper_max(values):
    """Return a rigorous upper bound, never a midpoint comparison."""
    bounds=[]
    for v in values:
        if not v.is_finite():
            raise ArithmeticError('nonfinite interval in maximum')
        bounds.append(v.upper())
    return max(bounds,default=ZERO)


class Inverse:
    def __init__(self,M,S,rho,b,Ahat):
        self.M,self.S,self.rho,self.b,self.Ahat=M,S,rho,b,Ahat
        self.rows=[(n,s) for n in range(1,M+1,2) for s in range(1,S+1)]+[(j,0) for j in range(1,M+1,2)]
        self.pos={k:i for i,k in enumerate(self.rows)}
        self.N=len(self.rows)
        self.cw=[rho**n if s else self.omega(n) for n,s in self.rows]

    def omega(self,j):
        return self.b**2*(j+1)*self.rho**j

    def split(self,f):
        vf=[arb(0) for _ in self.rows]
        dg={}
        harmonics={}
        for (n,s),v in f.items():
            if (n,s) in self.pos:
                vf[self.pos[n,s]]=v
            elif s==0:
                require(n > self.M, 'Failed: n > self.M')
                harmonics[n]=v
            else:
                dg[n,s]=v
        eta={}
        run=arb(0)
        for n in range(max(harmonics,default=self.M),self.M,-2):
            run=run+harmonics.get(n,ZERO)/(self.b**2*(n+1))
            eta[n]=run
            k=(n-2,1)
            if k in self.pos:
                vf[self.pos[k]]+=self.b**2*run
            else:
                acc(dg,k,self.b**2*run)
            if n==self.M+2:
                vf[self.pos[self.M,0]]+=self.b**2*(self.M+1)*run
        return vf,eta,dg

    def tailnorm(self,eta,dg):
        return sum((self.omega(j)*abs(v) for j,v in eta.items()),ZERO)+norm(dg,self.rho)

    def apply_many(self,fs,subtract=None):
        splits=[self.split(f) for f in fs]
        B=arb_mat([[sp[0][i] for sp in splits] for i in range(self.N)])
        X=self.Ahat*B
        results=[]
        for q,(v,e,g) in enumerate(splits):
            if subtract is not None:
                typ,idx=subtract[q]
                if typ=='finite':
                    X[idx,q]-=1
                elif typ=='shape':
                    e[idx]=e.get(idx,ZERO)-1
                else:
                    raise ValueError(typ)
            results.append(sum((abs(X[i,q])*self.cw[i] for i in range(self.N)),ZERO)+self.tailnorm(e,g))
        return results

from flint import arb


def far_shape(g,c,p,rho,jstar,M,S,Antail,H):
    b=c[1];b2=b*b
    q1={};q0={}
    # Delta e_j = 2 Re [w^(j-1)(j Gamma1 + Gamma0)].
    # These q's are the real coefficients of 2 i Gamma_boundary.
    for a,pa in p.items():
        for j,cj in c.items():
            if a==0 and j==1:continue
            acc(q1,(j,a),pa*cj)
            acc(q1,(0,a+j),-pa*cj)
    for a,pa in p.items():
        for bb,pb in p.items():
            if a==0 and bb==0:continue
            acc(q0,(a+1,bb),pa*pb)
    keys=sorted(set(q1)|set(q0))
    kmin=min(a-bb for a,bb in keys);kmax=max(a-bb for a,bb in keys)
    require(kmin <= 1, 'Failed: kmin <= 1')
    require(jstar + kmin - 1 > M + 2, 'Failed: jstar + kmin - 1 > M + 2')
    require(jstar > 2 * M + max(p) + 2, 'Failed: jstar > 2 * M + max(p) + 2')
    Y={};eps={}
    for a,bb in keys:
        k=a-bb
        acc(Y,k,q1.get((a,bb),ZERO))
        acc(eps,k,(abs(q1.get((a,bb),ZERO))*(1+bb)+abs(q0.get((a,bb),ZERO)))/(jstar+1))
    ks=list(range(kmin,kmax+1,2))
    require(all((k % 2 == 1 for k in Y)), 'Failed: all((k % 2 == 1 for k in Y))')
    eta={};de={}
    run=arb(0);err=arb(0)
    denom=jstar+kmin
    require(denom > 0, 'Failed: denom > 0')
    for k in reversed(ks):
        yy=Y.get(k,ZERO);ee=eps.get(k,ZERO)
        run+=yy
        ratio=arb(abs(k-1))/denom
        err+=ee*(1+ratio)+abs(yy)*ratio
        eta[k]=run/b2
        de[k]=err/b2
    E={k:abs(eta[k])+de[k] for k in ks}
    tau=1/(1-rho**-2)
    window=sum(((1+arb(max(0,k-1))/(jstar+1))*rho**(k-1)*E[k] for k in ks),ZERO)
    low=E[kmin]*rho**(kmin-3)*tau
    gcoupling=(sum((rho**(k-3)*E[k] for k in ks),ZERO)+E[kmin]*rho**(kmin-5)*tau)/(jstar+1)
    finite=H*E[kmin]/((jstar+1)*rho**jstar)
    nonharm=sum(((abs(q1.get((a,bb),ZERO))+abs(q0.get((a,bb),ZERO))/jstar)*bb*rho**(a-bb-1)/(b2*(jstar+a)) for a,bb in keys),ZERO)
    # K g = (1-|w|^2)^2 Q, using a positive connection of Jacobi alpha=2 to alpha=0.
    Q={}
    for n in range(1,M+1,2):
        C1=[{0:arb(1)}];C2=[{0:arb(1)}]
        for k in range(1,S):
            one={i:v*(k+n)/(k+n+1) for i,v in C1[-1].items()}
            one[k]=arb(2*k+n+1)/(k+n+1)
            two={i:v*(2*k+n+2)/(k+n+2) for i,v in one.items()}
            for i,v in C2[-1].items():acc(two,i,v*(k+n)/(k+n+2))
            C1.append(one);C2.append(two)
        for s in range(1,S+1):
            for t,v in C2[s-1].items():
                acc(Q,(n,t),g.get((n,s),ZERO)*v/(4*s*(s+1)))
    L=M+max(p)
    shifted=poly(sine(Q),p,True)
    for _ in range(L):shifted=z(shifted)
    require(all((m >= 0 for (m, t) in shifted)), 'Failed: all((m >= 0 for (m, t) in shifted))')
    require(jstar > L, 'Failed: jstar > L')
    wsum=sum((abs(v)*(2*t+1)*(2*t+3)*rho**(m-1-L) for (m,t),v in shifted.items()),ZERO)
    W=8*wsum/(b2*(jstar-L)**2)
    Wpre=Antail*W
    total=window+low+gcoupling+finite+nonharm+Wpre
    limit=sum((rho**(k-1)*abs(eta[k]) for k in ks),ZERO)+abs(eta[kmin])*rho**(kmin-3)*tau
    return dict(Zshape_far=total.upper(),shape_limit_bound=limit,far_window=window,far_low=low,far_gcoupling=gcoupling,
                far_finite=finite,far_nonharmonic=nonharm,far_W_unpreconditioned=W,far_W_preconditioned=Wpre,
                far_wsum=wsum,far_shift=L,far_kmin=kmin,far_kmax=kmax)

CENTRE = {'M': 41, 'S': 24, 'rho': [21, 20]}
CENTRE['g'] = """
-0x1.63ddd929de7bap+9 -0x1.324988ffa0fbcp+11 -0x1.9e929ff68ff96p+11 0x1.39d28d78b33cep+10
0x1.8dc317df3eb1cp+12 -0x1.b76635ae30c0dp+12 0x1.c5ec2d3473ceap+11 -0x1.291670bcd763dp+10
0x1.14f2ddf9f003bp+8 -0x1.88a37b24f59fdp+5 0x1.ba7179885cfa8p+2 -0x1.9984dc61736dfp-1
0x1.404510d188266p-4 -0x1.bc6e3b89946ddp-8 0x1.3b4fc8b4fe067p-11 -0x1.1efa72ff16054p-14
0x1.47903d371a65dp-17 -0x1.70a71e2bb6af5p-20 0x1.682f7fe737707p-23 -0x1.2c774c83f04f0p-26
0x1.b9a38abaf517ep-30 -0x1.332fd8306061cp-33 0x1.c7a8f791e559bp-37 -0x1.7c1268d0a25c3p-40
-0x1.1e1b5641d5afep+2 0x1.45a014419c875p+3 0x1.2523bb57e7d05p+5 -0x1.d1cc9dc6769b1p+0
-0x1.124304864ff1ap+5 0x1.a5adc308a019cp+2 0x1.01722d5c13242p+4 -0x1.ac4d17816835ap+3
0x1.5c6a299e6961ap+2 -0x1.72b7793f6e2bep+0 0x1.20d0942d0465ep-2 -0x1.64f29aa00ecfbp-5
0x1.79566165a260ep-8 -0x1.6b27328dbfdebp-11 0x1.48c4794edb3c6p-14 -0x1.1f2303a92c89ap-17
0x1.fe288f01be21cp-21 -0x1.de23ca2b6dc37p-24 0x1.c5010ebc0a370p-27 -0x1.98148ae8c5836p-30
0x1.5a5a22a248536p-33 -0x1.1d3680cd2d6d0p-36 0x1.85e4c2eff6545p-40 -0x1.69b8ff4805f38p-43
-0x1.9224021f154dfp+0 0x1.4381cfe891c1bp+4 0x1.1a2143eee9161p+4 -0x1.2fb12085ed60dp+2
-0x1.701405c06cad6p+4 0x1.9c6268d46f4d8p+3 0x1.20cf98f4235c1p+0 -0x1.b99c42cbdd1cap+1
0x1.a87fdb7208c0fp+0 -0x1.e1ffe87af4abfp-2 0x1.845732d66bc2ap-4 -0x1.e06df1f30494dp-7
0x1.de12769d16833p-10 -0x1.8b70e82f92a92p-13 0x1.1991f311d2f76p-16 -0x1.74f0ebaeeee18p-20
0x1.111d7d506956cp-23 -0x1.fc12b23c48498p-27 0x1.0957d06235038p-29 -0x1.03459b9704705p-32
0x1.bbd07cbaf1191p-36 -0x1.453ba8210fa93p-39 0x1.6381da96c9556p-42 -0x1.26546eeae90a4p-45
-0x1.2eca198ca5d44p+7 -0x1.0b0eba932384bp+8 0x1.9b379cfd5bd27p+8 -0x1.7e26baaa1f507p+7
0x1.519c22f683abap+5 -0x1.eb449abd9e8fcp+2 0x1.b72f0156c4038p+1 -0x1.9c478ab601a8cp+0
0x1.00d6550769b21p-1 -0x1.c1390ad694bdcp-4 0x1.286988923bca3p-6 -0x1.3723597e039abp-9
0x1.0eacce301ca5dp-12 -0x1.95183f1ccb581p-16 0x1.124ea5b275a4dp-19 -0x1.7610770d31571p-23
0x1.2d6c84c77bd1bp-26 -0x1.2485d60266ad1p-29 0x1.2157af5ee64afp-32 -0x1.002add1d51027p-35
0x1.6f3ffbf1ba1a3p-39 -0x1.511c4984e6574p-42 -0x1.1841d04eb2665p-42 0x1.3bb4952935d50p-42
0x1.29a085c959fd5p-1 -0x1.638138e5c8c1dp+2 0x1.49549ef7c937fp+3 -0x1.d0bf5dab0e05fp+2
0x1.1f70362fcab8bp+1 -0x1.cc0a460694df9p-3 -0x1.8a37d42d168afp-7 -0x1.677b67ec54052p-6
0x1.25fa18170636bp-6 -0x1.98b936402c848p-8 0x1.6ce18f2d4008bp-10 -0x1.dedfefbf896d8p-13
0x1.f35b30c86792ap-16 -0x1.b74c176615c89p-19 0x1.5603d0a09d1b6p-22 -0x1.ebe16d3e700f6p-26
0x1.5f2127e5f8c12p-29 -0x1.15ceffecb98e6p-32 0x1.f302712768429p-36 -0x1.92ff4736a239dp-39
0x1.3cbbfda875735p-42 0x1.a2f657cbb0e74p-45 -0x1.7ceac1df63412p-42 -0x1.41ee319ae5122p-44
-0x1.749f7043ea35ap+4 0x1.803d40e5c9335p+3 0x1.48c7ea8ffd9f7p+1 -0x1.ad9850798a6abp+1
0x1.e03a774717080p-1 -0x1.557d2fa9db646p-4 0x1.64fb69895144bp-8 -0x1.89e0a8dcb679cp-7
0x1.b1657aabc0bbdp-8 -0x1.00d91a82d64f9p-9 0x1.a0f6ee68d830bp-12 -0x1.00552b81acc80p-14
0x1.f9154207ae18fp-18 -0x1.9da8cec4d4341p-21 0x1.219f0adf7710bp-24 -0x1.69340db484882p-28
0x1.bf60cc22c8b5cp-32 -0x1.4bdf89fb47d8fp-35 0x1.281ace5720fc8p-38 -0x1.174a3d44256e2p-41
-0x1.c31fa11851426p-44 0x1.dba3c9fe44250p-45 0x1.abe013813d2f4p-42 -0x1.c57265e2c5a7ep-42
-0x1.82cbfaed341c6p+3 0x1.5e9234f517050p+2 0x1.1777c7d39c94bp+1 -0x1.1927dfe90a031p+1
0x1.5f6446370442ap-1 -0x1.cfd2e980fe299p-4 0x1.e6fdf876e9f41p-7 -0x1.0d0f88f19c3cfp-8
0x1.81304415b92d0p-10 -0x1.8eb4497307799p-12 0x1.29d56cfa7149dp-14 -0x1.558ff25480461p-17
0x1.3c82fe086d98fp-20 -0x1.ee919e92a0c6cp-24 0x1.529ca1e737673p-27 -0x1.aaba763b4f627p-31
0x1.11f8225ec2b26p-34 -0x1.abb58d1403406p-38 0x1.7b608c79e58d5p-41 -0x1.3191871243e6bp-44
0x1.227e5a4815e32p-43 -0x1.89305a2643038p-45 -0x1.97c0b1df29bdep-44 0x1.25f70effcc87cp-43
-0x1.2661de02c66c3p+0 0x1.8690e37ca850fp-1 -0x1.daa6df3ed0b58p-6 -0x1.14024cf09e18fp-3
0x1.e5e78d80d8b3cp-5 -0x1.6e6af9da9ebfap-7 0x1.b7d77a50b150ap-11 -0x1.0fbfe3d2e368ap-14
0x1.2e9e1d66f51ebp-14 -0x1.07910e0cd56cep-15 0x1.035c1be8045fdp-17 -0x1.6383956feda7bp-20
0x1.7674ad4f5dd4ap-23 -0x1.4299cef52307ap-26 0x1.dc4eb086f09fep-30 -0x1.393aa50c04474p-33
0x1.81ab684804f1dp-37 -0x1.9e1b7163372b4p-41 0x1.30215e26cbe3cp-43 -0x1.2275a066a78c4p-44
-0x1.a9a70376983b4p-43 -0x1.c8581c9cc8d50p-46 -0x1.c93a680a47fb8p-45 0x1.bb73c505908fap-43
-0x1.26bf4d3e62e72p+0 0x1.894699ffa168fp-1 -0x1.6c585f924db5ep-4 -0x1.127930b43c36dp-4
0x1.eeb4817fca1fap-6 -0x1.6f5931c66cf78p-8 0x1.1f42bfaf14af8p-11 -0x1.1f021ce7c2a1ep-14
0x1.ebcbbdc87bec4p-16 -0x1.4574f30f74e90p-17 0x1.1a0f15ccb9ee9p-19 -0x1.643fba1489dfdp-22
0x1.60d6210d3f256p-25 -0x1.1fa08ddbb6a36p-28 0x1.8f43cae76ab5fp-32 -0x1.ebcc3c02a626ap-36
0x1.1fda5bffd0e7ap-39 -0x1.a0f8675118685p-43 0x1.04a9d9a4a1cfap-43 -0x1.48b78a3a31340p-43
0x1.720e84a8f1960p-45 0x1.051aefaf5cf84p-41 0x1.b55c59d112038p-47 0x1.ecb352dca7f81p-43
-0x1.795d15713dcbfp-2 0x1.fecbaad45df93p-3 -0x1.0da5d186d1226p-5 -0x1.3ed0544c77a6ep-6
0x1.3787ef8825b59p-7 -0x1.0116fb3b84002p-9 0x1.edf3cb9578776p-13 -0x1.a44ff2b5b329ap-16
0x1.8a0da4e50f590p-18 -0x1.c86573fab1d5dp-20 0x1.838b446c57247p-22 -0x1.e4f5b4b162cc5p-25
0x1.daa294d1c9ef5p-28 -0x1.7e455654f09f5p-31 0x1.094b96bd7e42cp-34 -0x1.31b7c78e292cap-38
0x1.802ba4965b254p-44 0x1.7bf87c9240a88p-43 -0x1.dc0a6a0c7aabap-43 0x1.45ac0d9ed03fap-44
-0x1.b76d9f7fea7c8p-44 -0x1.58c07787c63f8p-44 -0x1.aaf972f2305adp-42 0x1.97f396e153900p-44
-0x1.29bb09b2a35fep-4 0x1.c0bc814634432p-5 -0x1.a71eddac15a65p-7 -0x1.f88a8bd724078p-11
0x1.3cb19ef88bc95p-10 -0x1.36fa9477cf285p-12 0x1.2d944e37e502ap-15 -0x1.496f15bc50b0dp-19
0x1.ec987069258c7p-22 -0x1.a3e0356a202c7p-23 0x1.bab89278ac028p-25 -0x1.394cfa5f68b33p-27
0x1.4c333c37f4758p-30 -0x1.1b28e0b0a48a2p-33 0x1.92f75efad912dp-37 -0x1.4848d7ca08732p-40
-0x1.8c2c1378e0722p-43 0x1.e77a10b03a385p-43 0x1.9494771a29b60p-45 -0x1.6f1ce3767bcacp-44
-0x1.6514e50fce448p-43 0x1.31f1950e5cc9ap-42 0x1.9b5e13a520b6fp-43 -0x1.fce175e178ee6p-44
-0x1.4465bebc292a9p-5 0x1.e559955551574p-6 -0x1.d2534a2995226p-8 -0x1.c2fef6de51645p-13
0x1.0084834c7d2f3p-11 -0x1.fda5586c95936p-14 0x1.056307b874ea3p-16 -0x1.62ac8598b14bbp-20
0x1.9d2afdf214f59p-23 -0x1.e6743287f18bbp-25 0x1.c22e06f2c3d7cp-27 -0x1.294d6082f08dep-29
0x1.2c0956e50130ap-32 -0x1.eaeba93303ddbp-36 0x1.7473a23aec3b7p-39 0x1.e76d3da32b7c8p-44
-0x1.b01bd0346dea3p-43 -0x1.8ab02ebf4ae8fp-43 0x1.a0c8bdb6cb454p-42 -0x1.7a5794ef54eb0p-45
-0x1.440ed1b12badep-43 -0x1.ecf55d93d4641p-45 -0x1.7c16d64d147fdp-44 0x1.1040a452ed6c6p-42
-0x1.58865181e0e14p-7 0x1.067aeacd06b6bp-7 -0x1.0d942d4f90feep-9 0x1.3a38f598becc4p-16
0x1.dbbcaa2a209a0p-14 -0x1.04e1f3f6a25e9p-15 0x1.250d21fa088f6p-18 -0x1.a4607fda39da2p-22
0x1.5f9abcdd845c7p-25 -0x1.4a24539e041fcp-27 0x1.3313bc1a72589p-29 -0x1.a04c341039546p-32
0x1.aa7f2f0b31033p-35 -0x1.626f1b695a50cp-38 0x1.cf23dd852d2bep-42 0x1.26c9f7de103e8p-43
0x1.00fba31e091e0p-43 -0x1.d0c76c202aab8p-47 -0x1.498b976bbebc8p-42 0x1.71d954a2541fbp-43
-0x1.1d9c016fd058ap-43 -0x1.2013083c24accp-43 -0x1.cf6c5e03012e0p-44 -0x1.79dbe2c560939p-46
-0x1.70d80c66e0b7cp-9 0x1.23927c201224ep-9 -0x1.53ab90b3acc4bp-11 0x1.bb88261f7ab02p-15
0x1.23cb1a7b99d6fp-16 -0x1.9d6d6e400c67dp-18 0x1.f256fce1648c5p-21 -0x1.5856ab559fe07p-24
0x1.b40d80fbe3a4cp-28 -0x1.90faf78b2fda0p-30 0x1.a75d4c7f4d931p-32 -0x1.35ee200dfea5fp-34
0x1.4c11d392416edp-37 -0x1.2dfda5641211cp-40 0x1.40b73a39e99a8p-42 0x1.65316295ba660p-46
-0x1.34ac74b233d81p-43 0x1.604121a73c3c3p-42 -0x1.a35f2a6d26e90p-45 0x1.ceb62513ada84p-44
0x1.005cd7e1baeb8p-43 -0x1.f461d17f65329p-43 0x1.5665496d48adep-44 -0x1.7d17ae3d0fb10p-48
-0x1.27b748920a54fp-10 0x1.d2eaa3239adbcp-11 -0x1.11aeace317faep-12 0x1.966924af2dc53p-16
0x1.66edb10bc3d48p-18 -0x1.14e61c5983f4ap-19 0x1.590687050b633p-22 -0x1.016663113e0bbp-25
0x1.4cf1071ee0f74p-29 -0x1.c0a7fee2a2d59p-32 0x1.971d8aa4e60e0p-34 -0x1.1816b6fb2a21ep-36
0x1.4d9ec66b10750p-39 -0x1.b4f2bef9683aap-42 -0x1.1942b2ee95e70p-42 0x1.2703e4c0a59ccp-42
-0x1.650442d01f92bp-45 -0x1.6dbc1684200dep-43 0x1.fc0bfd4cee348p-45 0x1.a24a8238eec00p-47
0x1.e8b3e8f8571b2p-44 0x1.c01a681a26948p-43 0x1.252cfde085930p-44 -0x1.8a758830a1e98p-44
-0x1.32b026ed75d31p-12 0x1.eb8b22de19449p-13 -0x1.2ce536b517681p-14 0x1.1294e300bb8a9p-17
0x1.e886870e9e3c9p-21 -0x1.f669a59647cc4p-22 0x1.5375900294e57p-24 -0x1.0bbe87336eab3p-27
0x1.3f3650688fb75p-31 -0x1.47c8ddd933d4dp-34 0x1.1cde4b5c73b96p-36 -0x1.bc97f0b881838p-39
0x1.62a4d3a89e2a1p-42 0x1.62a1c4e6a13f4p-46 0x1.28bed043baadep-43 0x1.7486115aa7202p-44
-0x1.5fe3197f7237ep-43 -0x1.3b1156633982cp-45 0x1.e6c544a169ec8p-44 -0x1.165a2644ce528p-45
0x1.5376f27441977p-43 -0x1.e7e5fa9b2ad52p-42 -0x1.aeeb012b03ee0p-47 -0x1.2958903d774f8p-45
-0x1.7069eb570b4c8p-14 0x1.2b81bac95d59ap-14 -0x1.7eacc5f51e681p-16 0x1.a2389b400d65fp-19
0x1.51ce2d4696699p-24 -0x1.bf72c2cac6b35p-24 0x1.4b39c4d263819p-26 -0x1.0d9538413e9bcp-29
0x1.309e445ea72c7p-33 -0x1.032c6fac5d048p-36 0x1.8929ea3d94a32p-39 -0x1.87b04544585cep-42
-0x1.2d52911415134p-44 0x1.0cc0745f261e4p-44 -0x1.3b5969dc0621cp-42 -0x1.5afd0341a8a80p-44
0x1.4c57e3509677bp-42 -0x1.45bed5b8af5b8p-45 -0x1.0faa686af7909p-43 0x1.fd228313655a6p-45
0x1.ea484a9539027p-42 0x1.342751d8e54c1p-43 -0x1.1657214403732p-43 0x1.e4303fff74dfep-45
-0x1.f0fefa7229572p-16 0x1.94bcd17bfc97ep-16 -0x1.0482a192bbfd1p-17 0x1.297fc612a3a38p-20
0x1.8f2ddda121c9fp-29 -0x1.fd9da7ddd0a9cp-26 0x1.8d110953a8675p-28 -0x1.51c058afbe8a8p-31
0x1.8e7a546c61dc4p-35 -0x1.1bf87be63adcdp-38 0x1.fb74bf6cb5c52p-41 -0x1.e10d80889fc1ap-44
0x1.4e0d31c5b6ceap-42 -0x1.207ef0711dab8p-42 -0x1.e04122daf54d0p-43 0x1.b97bf98b01a0fp-42
0x1.aa3ea68761fc0p-47 0x1.9c70f1cb7b822p-45 -0x1.928e3d9769060p-44 -0x1.949ad612b3a95p-43
-0x1.53335028d0b89p-43 0x1.1493626e358e0p-45 -0x1.a92b5301ecb80p-49 0x1.c6888c1d49acdp-43
-0x1.08f984881b871p-17 0x1.b414adbd7fb89p-18 -0x1.2070a9191aca7p-19 0x1.65f664f13fcabp-22
-0x1.3f4df7c77a98ap-27 -0x1.ae99c260077e0p-28 0x1.7a1a537567683p-30 -0x1.5274079deece3p-33
0x1.8556366ecefccp-37 -0x1.3eef055912458p-40 0x1.2b570547b02bep-41 -0x1.baf86d64c08e6p-42
-0x1.2e27ecf666934p-43 0x1.282aea12e8b00p-42 0x1.c2df7be4aaaf1p-43 -0x1.13877e40c235ap-41
-0x1.8c05e0a0c9150p-44 0x1.fc8162f6a34e8p-46 -0x1.279bf474233f4p-43 -0x1.60d22e039a2b6p-44
0x1.091b9804a1a6fp-43 0x1.23b4931011f77p-42 0x1.7c7fd19db8ea8p-45 -0x1.ce9558c9d76b6p-43
-0x1.459f9823eabd6p-19 0x1.0da6d96906552p-19 -0x1.6b2ea3c972d9fp-21 0x1.ddad169819100p-24
-0x1.804493f4e4a4ap-28 -0x1.7b8231e93def9p-30 0x1.81ac77774f98dp-32 -0x1.6ca6b821c9d75p-35
0x1.017105a3b20b9p-38 -0x1.f47d38eb0b78cp-45 -0x1.43128bd35256dp-42 0x1.637251ae17805p-42
0x1.103c2effc9025p-42 -0x1.312d38ebb4d0bp-42 -0x1.4dfd3d47ea216p-42 0x1.baccc3b78aefdp-43
0x1.f1c8091f0b3fcp-43 -0x1.0034a8495cbb7p-43 0x1.16c7c6b550a69p-43 0x1.2c1651fccd0dbp-43
-0x1.7c86546a18032p-45 -0x1.702d8c135710ap-43 0x1.5fb2d8224eafep-43 -0x1.7da79b30ced13p-43
-0x1.9033769158565p-21 0x1.4c224202d3d7dp-21 -0x1.c2ec3b06d743ep-23 0x1.30d2f818a0effp-25
-0x1.29bbe1e129141p-29 -0x1.6b93b96532095p-32 0x1.a52995500ac95p-34 -0x1.99d37a1ec65f9p-37
0x1.951774f9d407bp-41 -0x1.67636b3990a48p-44 0x1.03d03671620c7p-43 0x1.71b4e9eb97c12p-44
-0x1.3e26df2b454f3p-42 0x1.d07ad572d3a08p-43 0x1.f3f3be1f52a37p-44 -0x1.bb6283b7b3e8cp-44
-0x1.6bb751d37deb5p-42 0x1.7b94f13390e82p-42 -0x1.bb077a401c976p-45 -0x1.c2fe61ac86918p-45
0x1.d6d6b61c3cdf4p-44 0x1.7e3d664180860p-47 -0x1.a7aede1e4d4a0p-44 0x1.48ce6230643e8p-44
""".split()
CENTRE['c'] = """
0x1.dae1a1e0bb837p+3 0x1.7589bfe6d254bp-4 0x1.81f65da739513p-4 0x1.8aaf75e3b7500p-4
0x1.41c6ef6ce8a22p-8 0x1.6579de4b2afe9p-8 0x1.cd5f18327b429p-10 0x1.71528c206f046p-13
0x1.c719084ad329dp-14 0x1.c404a1255983cp-16 0x1.43b14c14bea05p-18 0x1.12caa5bbcb15ep-19
0x1.efe0af2707320p-22 0x1.d6abe46bb6255p-24 0x1.3f4166576214cp-25 0x1.24be82f894c89p-27
0x1.3933505b371a4p-29 0x1.7346a7d82f022p-31 0x1.6097beb9ab9bep-33 0x1.740da5729e3b1p-35
0x1.4dfa474b1759bp-37
""".split()

import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
if hasattr(os,'getpriority'):
    increase=19-os.getpriority(os.PRIO_PROCESS,0)
    if increase>0:os.nice(increase)
import argparse
import json
import time
from pathlib import Path
from flint import arb, arb_mat, ctx

START=time.monotonic()


def log(*s):
    print('[%8.2fs]'%(time.monotonic()-START),*s,flush=True)


def dyadic(h):
    n,d=float.fromhex(h).as_integer_ratio()
    x=arb(n)/d
    require(x.is_exact(), 'Failed: x.is_exact()')
    return x


def checkpoint(values,path):
    Path(path).write_text(json.dumps({k:(v.str(30) if isinstance(v,arb) else v) for k,v in values.items()},indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['base','all'],default='all')
    parser.add_argument('--bits',type=int,default=128)
    parser.add_argument('--jstar',type=int,default=1001)
    parser.add_argument('--output',default='certificate_r4.json')
    args=parser.parse_args()
    require(args.bits >= 96, 'Failed: args.bits >= 96')
    ctx.prec=args.bits
    ctx.threads=1
    data=CENTRE
    M,S=data['M'],data['S']
    require(M == 41 and S == 24 and (len(data['g']) == 504) and (len(data['c']) == 21), "Failed: M == 41 and S == 24 and (len(data['g']) == 504) and (len(data['c']) == 21)")
    rho=arb(data['rho'][0])/data['rho'][1]
    rows=[(n,s) for n in range(1,M+1,2) for s in range(1,S+1)]+[(n,0) for n in range(1,M+1,2)]
    nG=(M+1)//2*S
    g={k:dyadic(v) for k,v in zip(rows[:nG],data['g'])}
    c={j:dyadic(v) for j,v in zip(range(1,M+1,2),data['c'])}
    p={j-1:v*j for j,v in c.items()}
    b=c[1]
    V=add(K(g),{(j,0):v for j,v in c.items()})
    Vy=sine(V)
    P=sum((abs(v)*rho**a for a,v in p.items()),ZERO)
    Vnorm=norm(V,rho)
    result={'bits':ctx.prec,'M':M,'S':S,'rho':rho,'b':b,'Pprime':P,'Vnorm':Vnorm}
    log('centre',b,'Pprime',P,'Vnorm',Vnorm)
    F0=add(g,unsine(abs2(Vy,p)))
    log('residual fully expanded',len(F0),'coefficients')
    cols=[]
    for i,k in enumerate(rows[:nG]):
        e={k:arb(1)}
        cols.append(add(e,unsine(abs2(sine(K(e)),p))))
        if i%40==39:log('finite g columns',i+1,nG)
    P1=poly(Vy,p,True)
    Q=abs2({(1,0):arb(1)},p)
    Aj,Bj=P1,Q
    for j in range(1,M+1,2):
        cols.append(scale(unsine(add(scale(Aj,arb(j)),scale(Bj,-arb(1)/2))),arb(2)))
        Aj=z(z(Aj));Bj=z(z(Bj))
    log('all finite columns',len(cols))
    N=len(rows)
    Mf=arb_mat([[col.get(k,ZERO) for col in cols] for k in rows])
    log('Arb inverse begin',N)
    enclosed_inverse=Mf.inv()
    # Ahat is exactly the dyadic midpoint of the rigorous enclosure, not an interval operator.
    Ahat=arb_mat([[enclosed_inverse[i,j].mid() for j in range(N)] for i in range(N)])
    inv=Inverse(M,S,rho,b,Ahat)
    cw=inv.cw
    E=arb_mat(N,N)
    for i in range(N):E[i,i]=1
    E=E-Ahat*Mf
    result['inverse_defect']=upper_max(sum((abs(E[i,j])*cw[i] for i in range(N)),ZERO)/cw[j] for j in range(N))
    require(result['inverse_defect'] < 1, "Failed: result['inverse_defect'] < 1")
    Anfin=upper_max(sum((abs(Ahat[i,j])*cw[i] for i in range(N)),ZERO)/rho**rows[j][0] for j in range(N))
    # Exact finite correction caused by a harmonic tail residual.
    H=sum((abs(Ahat[i,inv.pos[M,1]]+(M+1)*Ahat[i,inv.pos[M,0]])*cw[i] for i in range(N)),ZERO)
    hvals=[];hs=arb(0);gs=arb(0)
    NH=M+402
    for n in range(M+2,NH+1,2):
        hs+=(n+1)*rho**n
        if n>M+2:gs+=rho**(n-2)
        hvals.append((H+hs+gs)/((n+1)*rho**n))
    tau=1/(1-rho**-2)
    beyond=tau+tau/((NH+3)*rho**2)+H/((NH+3)*rho**(NH+2))
    Antail=upper_max(hvals+[beyond,arb(1)])
    An=upper_max([Anfin,Antail])
    result.update(normA=An,normR=Anfin,normA_tail=Antail,H=H)
    log('A inverse defect',result['inverse_defect'],'normA',An,'tail',Antail)
    result['Y']=inv.apply_many([F0])[0].upper()
    log('Y',result['Y'])
    zfin=[]
    for q in range(0,N,32):
        fs=cols[q:q+32]
        vv=inv.apply_many(fs,[('finite',i) for i in range(q,q+len(fs))])
        zfin.extend(v/cw[i] for i,v in enumerate(vv,start=q))
    result['Zfin']=upper_max(zfin)
    result['Zfin_argmax']=max(range(N),key=lambda i:zfin[i].upper())
    log('Zfin',result['Zfin'],rows[result['Zfin_argmax']])
    checkpoint(result,str(Path(args.output).with_suffix('.base.json')))
    if args.stage=='base':
        log('BASE COMPLETE; no theorem certification in partial mode')
        return
    # Boundary g-tail rectangle. The omitted polynomial factors are bounded in norm.
    degree=max(p);ps={a:v for a,v in p.items() if a<=16}
    Ps=sum((abs(v)*rho**a for a,v in ps.items()),ZERO)
    Pt=sum((abs(v)*rho**a for a,v in p.items() if a>16),ZERO)
    delta_a=2*Ps*Pt+Pt**2
    nmax=M+2*degree+2;smax=S+degree+3
    gb=[(n,s) for n in range(1,nmax+1,2) for s in range(1,smax+1) if not(n<=M and s<=S)]
    Zgb=arb(0);gbarg=None
    for q in range(0,len(gb),32):
        ks=gb[q:q+32]
        fs=[unsine(abs2(sine(K({k:arb(1)})),ps)) for k in ks]
        vals=inv.apply_many(fs)
        for (n,s),v in zip(ks,vals):
            d=n+2*s
            bound=v/rho**n+An*delta_a/(d*(d+2))
            require(bound.is_finite(), 'Failed: bound.is_finite()')
            if bound.upper()>Zgb:
                Zgb=bound.upper();gbarg=(n,s)
        if q%256==0:log('boundary g',q,len(gb),'max',Zgb,'arg',gbarg)
    result.update(Zgb=Zgb,Zgb_argmax=gbarg,Zgb_count=len(gb),g_short_degree=16,g_remainder_norm=delta_a)
    Ds=1+2*(smax+1);Dn=nmax+4
    result['Zg_far']=upper_max([P**2/(Ds*(Ds+2)),Antail*P**2/(Dn*(Dn+2))])
    log('g tails',Zgb,result['Zg_far'])
    # Consecutive shape columns, complete supports, no numerical truncation.
    Zsh=arb(0);sharg=None
    pending=[];pjs=[]
    for j in range(M+2,args.jstar+1,2):
        f=scale(unsine(add(scale(Aj,arb(j)),scale(Bj,-arb(1)/2))),arb(2))
        pending.append(f);pjs.append(j)
        if len(pending)==32 or j+2>args.jstar:
            vals=inv.apply_many(pending,[('shape',i) for i in pjs])
            for jj,val in zip(pjs,vals):
                zz=val/inv.omega(jj)
                require(zz.is_finite(), 'Failed: zz.is_finite()')
                if zz.upper()>Zsh:Zsh=zz.upper();sharg=jj
            pending=[];pjs=[]
            log('shape columns through',j,'max',Zsh)
        Aj=z(z(Aj));Bj=z(z(Bj))
    result.update(Zsh=Zsh,Zsh_argmax=sharg,jstar=args.jstar)
    result.update(far_shape(g,c,p,rho,args.jstar,M,S,Antail,H))
    log('far shape bound',result['Zshape_far'])
    Z=upper_max([result[k] for k in ['Zfin','Zgb','Zg_far','Zsh','Zshape_far']])
    theta1=1/(b**2*rho);theta0=1/(2*b**2);kap=arb(1)/15
    C2=An*upper_max([P*theta1*kap,2*P*theta1*theta0+theta1**2*Vnorm])
    C3=An*theta1**2*upper_max([kap,theta0])
    result.update(Z=Z,C2=C2,C3=C3)
    # Exact rational majorants are accepted only after all interval comparisons.
    Yq=arb(7)/1000000;Zq=arb(37)/100;C2q=arb(4);C3q=arb(1)/1000
    r=arb(1)/50000
    require(result['Y'] < Yq and Z < Zq and (C2 < C2q) and (C3 < C3q), "Failed: result['Y'] < Yq and Z < Zq and (C2 < C2q) and (C3 < C3q)")
    radius_poly=Yq+(Zq-1)*r+C2q*r*r+C3q*r*r*r
    derivative_poly=Zq-1+2*C2q*r+3*C3q*r*r
    require(radius_poly < 0 and derivative_poly < 0, 'Failed: radius_poly < 0 and derivative_poly < 0')
    dc1=r/inv.omega(1)
    c1lo=b-dc1
    # Sup_j j/(b^2(j+1)rho^j) <= 1/(b^2 rho), valid also for all omitted coefficients.
    derivative_tail=sum((j*abs(v) for j,v in c.items() if j>=3),ZERO)+r/(b**2*rho)
    map_tail=sum((abs(v) for j,v in c.items() if j>=3),ZERO)+r/(4*b**2*rho**3)
    sig1=derivative_tail/c1lo;sig0=map_tail/c1lo
    star=(sig1+sig0)/(1-sig0)
    c3lo=abs(c[3])-r/inv.omega(3)
    extension_radius=arb(41)/40
    extension_margin=c1lo-sum((j*abs(v)*extension_radius**(j-1) for j,v in c.items() if j>=3),ZERO)-r/(b**2*rho)
    require(1 < extension_radius and extension_radius < rho and (extension_margin > 0), 'Failed: 1 < extension_radius and extension_radius < rho and (extension_margin > 0)')
    require(c1lo > 0 and sig0 < 1 and (sig1 < 1) and (star < 1) and (c3lo > 0), 'Failed: c1lo > 0 and sig0 < 1 and (sig1 < 1) and (star < 1) and (c3lo > 0)')
    result.update(radius=r,radii_polynomial=radius_poly,derivative_polynomial=derivative_poly,
                  c1_lower=c1lo,derivative_defect=sig1,star_defect=star,abs_c3_lower=c3lo,
                  extension_radius=extension_radius,extension_derivative_margin=extension_margin,
                  status='PROVED')
    checkpoint(result,args.output)
    log('radii polynomial',radius_poly,'derivative',derivative_poly)
    log('geometry: derivative defect',sig1,'star defect',star,'c3 lower',c3lo)
    print('PROVED',flush=True)


if __name__=='__main__':
    main()

