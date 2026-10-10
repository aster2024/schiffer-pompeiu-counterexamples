"""Exact fmpq check that the finite Dirichlet field factors by 1-r^2."""
from proof_guard import require
import argparse
import hashlib
import json
from flint import fmpq


def rat_pair(man,exp):
    return fmpq(int(man)*(2**int(exp))) if exp>=0 else fmpq(int(man),2**(-int(exp)))


def rat_hex(value):
    num,den=float.fromhex(value).as_integer_ratio()
    return fmpq(num,den)


def up(n,s):
    return fmpq((s+1)*(n+s+1),(n+2*s+1)*(n+2*s+2))


def down(n,s):
    return fmpq(0) if s==0 else fmpq(s*(n+s),(n+2*s)*(n+2*s+1))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('center')
    a=ap.parse_args()
    raw=open(a.center,'rb').read();d=json.loads(raw)
    m,M,S=d['m'],d['M'],d['S']
    ms=list(range(m,M+1,2*m))
    keys=[(n,s) for n in ms for s in range(S+1)]
    if 'g_dyadic' in d:g={k:rat_pair(*v) for k,v in zip(keys,d['g_dyadic'])}
    else:g={k:rat_hex(v) for k,v in zip(keys,d['g_hex'])}
    require((len(g)==len(keys)), 'verify_factor_dirichlet_rank2.py:35: proof gate failed')
    checked=0
    for n in ms:
        V=[fmpq(0) for _ in range(S+2)]
        for s in range(S+1):
            v=g[n,s]
            if s==0:
                z=v/fmpq(4*(n+1)*(n+2))
                V[0]-=z;V[1]+=z
            else:
                D=n+2*s
                V[s-1]+=v/fmpq(4*D*(D+1))
                V[s]-=v/fmpq(2*D*(D+2))
                V[s+1]+=v/fmpq(4*(D+1)*(D+2))
        q=[fmpq(0) for _ in range(S+3)]
        for s in range(S,-1,-1):
            q[s]=((up(n,s+1)+down(n,s+1))*q[s+1]
                  -down(n,s+2)*q[s+2]-V[s+1])/up(n,s)
        for s in range(S+2):
            val=(up(n,s)+down(n,s))*q[s]
            if s:val-=up(n,s-1)*q[s-1]
            val-=down(n,s+1)*q[s+1]
            require((val==V[s]), (n,s))
            checked+=1
    print('EXACT_DIRICHLET_FACTOR_ONLY','m',m,'M',M,'S',S,
          'checked coefficients',checked,
          'center SHA256',hashlib.sha256(raw).hexdigest())


if __name__=='__main__':
    main()
