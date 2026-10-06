# Geometry of the frozen centre (72 modes), unit-conformal-radius normalisation
import json, numpy as np, mpmath as mp
mp.mp.dps = 60
from pathlib import Path
D = json.loads((Path(__file__).resolve().parents[2] / 'berenstein' / 'centre_frozen.json').read_text())
m = D['m']
cmp_ = [mp.mpf(x) for x in D['p']]
rho_c = float(cmp_[0])
cu = np.array([float(x/cmp_[0]) for x in cmp_])      # unit conformal radius coefficients
cj = np.array([float(x) for x in cmp_])               # physical (eigenvalue 1) coefficients
def maps(coef, s):
    w = np.exp(1j*s)
    psi = np.zeros_like(w); dpsi = np.zeros_like(w); ddpsi = np.zeros_like(w)
    winv = np.exp(-1j*s)
    for j,c in enumerate(coef):
        e = m*j+1
        u = np.exp(1j*(m*j)*s)          # w^(m j)
        psi += c*w*u
        dpsi += c*e*u
        ddpsi += c*e*(e-1)*u*winv       # w^(m j - 1)
    return w, psi, dpsi, ddpsi
if __name__ == '__main__':
    print('rho_c', rho_c, 'alpha_u = rho_c^2 =', rho_c**2)
    print('c_j physical, j=0..12:', cj[:13])
    print('c_j/c_0:', cu[:8])
    N = 1<<16
    s = 2*np.pi*np.arange(N)/N
    w, psi, dpsi, ddpsi = maps(cu, s)
    cf = 1 + w*ddpsi/dpsi
    kap = cf.real/np.abs(dpsi)           # curvature (unit conformal radius)
    print('convexity factor min/max', cf.real.min(), cf.real.max())
    print('curvature min/max (unit scale)', kap.min(), kap.max(), ' radius of curvature', 1/kap.max(), 1/kap.min())
    print('|psi| min/max', np.abs(psi).min(), np.abs(psi).max())
    phi = s + np.angle(dpsi)             # outward normal angle
    h = (psi*np.exp(-1j*phi)).real; hp = (psi*np.exp(-1j*phi)).imag
    print('support function min/max', h.min(), h.max(), ' max|h\'|', np.abs(hp).max())
    # width in direction theta: h(theta)+h(theta+pi); need pairing -> do by interpolation
    phiu = np.unwrap(phi)
    print('phi range', phiu[0], phiu[-1]-phiu[0])
