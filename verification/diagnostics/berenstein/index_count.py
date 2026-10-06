import numpy as np
from scipy.special import jn_zeros, jv
from geom import cu, maps, rho_c, m
N=1<<16
s=2*np.pi*np.arange(N)/N
w,psi,d,dd=maps(cu,s)
Rin = rho_c*np.abs(psi).min(); Rout = rho_c*np.abs(psi).max()
print('R_in = %.9f, R_out = %.9f (physical units; centre of the D_77 symmetry at origin)'%(Rin,Rout))
def count(R):
    tot=0; per=[]
    for n in range(0, int(R)+2):
        nt = int(R/3)+5
        z = jn_zeros(n, nt) if n>0 else jn_zeros(0, nt)
        assert z[-1] > R, (n, z[-1])
        c = int(np.sum(z<=R))
        per.append(c)
        tot += c if n==0 else 2*c
    return tot, per
for R in [Rin, rho_c, Rout]:
    tot, per = count(R)
    print('R=%.9f  #Dirichlet eigenvalues of the disc of radius R that are <= 1: %d' % (R, tot))
# also Neumann-type (for reference), counted via zeros of J_n'
from scipy.special import jnp_zeros
def countN(R):
    tot=0
    for n in range(0,int(R)+2):
        nt=int(R/3)+5
        z = jnp_zeros(n, nt)
        assert z[-1]>R
        c=int(np.sum(z<=R))
        tot += (c if n==0 else 2*c)
    return tot+1   # mu_1=0
for R in [Rin, rho_c, Rout]:
    print('Neumann: R=%.9f count(<=1, incl. mu_1=0) = %d'%(R,countN(R)))
# Weyl check
print('Weyl: A/(4pi) - P/(4pi)sqrt(1) =', np.pi*rho_c**2/(4*np.pi) - 2*np.pi*rho_c/(4*np.pi))
