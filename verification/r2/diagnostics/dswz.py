#!/usr/bin/env python3
"""Stationary-phase quantities for the comparison with Dai-Sun-Wei-Zhang
(arXiv:2511.19819v2), evaluated for the forty-mode centre normalized to
conformal radius one (binary64, trapezoid rule; no enclosures).

Reads ../center_j40_frozen.json.  Notation as in Section 8 of the paper:
h support function in the normal angle, varrho = h + h'' radius of curvature,
  B    = 2 max sqrt(varrho) sinh|h'|                    (two-point amplitude defect),
  T(s) = sqrt(2 pi (alpha + s^2)) max 2 sqrt(varrho) sinh(s|h'|)   (dilate s*Omega),
  R_D  = I_D - leading two-point term at t = 0,  I_D = int exp(i lambda xi.x) ds.
None of these numbers enters the existence proof.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
centre = json.loads((HERE.parent / "center_j40_frozen.json").read_text())
m = int(centre["m"])
c = np.array([float(x) for x in centre["p"]])
rho = c[0]
alpha = rho * rho
c = c / rho                                  # conformal radius one
deg = m * np.arange(len(c)) + 1


def curve(s):
    """psi, psi', w*psi'' at w = exp(i s) for the normalized centre."""
    s = np.asarray(s, float)
    z = np.zeros(s.shape, complex)
    p = np.zeros(s.shape, complex)
    q = np.zeros(s.shape, complex)
    w = np.exp(1j * s)
    for j in range(len(c)):
        wd = np.exp(1j * deg[j] * s)
        z += c[j] * wd
        p += deg[j] * c[j] * wd / w
        q += deg[j] * (deg[j] - 1) * c[j] * wd / w
    return w, z, p, q


def support(theta):
    """h, h', varrho at normal angle theta (Newton iteration for the boundary point)."""
    theta = np.asarray(theta, float)
    s = theta.copy()
    for _ in range(60):
        w, z, p, q = curve(s)
        s = s - np.angle(w * p * np.exp(-1j * theta)) / np.real(1 + q / p)
    w, z, p, q = curve(s)
    assert np.max(np.abs(np.angle(w * p * np.exp(-1j * theta)))) < 1e-11
    zn = z * np.exp(-1j * theta)
    return zn.real, zn.imag, np.abs(p) / np.real(1 + q / p)


print("alpha = rho^2 = %.4f (eigenvalue at conformal radius one)" % alpha)
th = np.linspace(0.0, 2 * np.pi / m, 20001)
h, hp, vr = support(th)
B = (2 * np.sqrt(vr) * np.sinh(np.abs(hp))).max()
lam1 = np.sqrt(alpha + 1)
print("sup|h'| = %.6e   varrho in [%.6f, %.6f]" % (np.abs(hp).max(), vr.min(), vr.max()))
print("B = 2 max sqrt(varrho) sinh|h'| = %.10f" % B)
print("sqrt(2 pi) * sqrt(alpha+1) * B   = %.4f" % (np.sqrt(2 * np.pi) * lam1 * B))


def T(s):
    return np.sqrt(2 * np.pi * (alpha + s * s)) * (2 * np.sqrt(vr) * np.sinh(s * np.abs(hp))).max()


smax = np.sqrt(alpha / 8)
print("dilates s*Omega have eigenvalue alpha/s^2 >= 8 iff s <= %.3f" % smax)
for s in (1.0, 10.0, 100.0, 200.0, smax):
    print("  s = %8.3f   T(s) = %10.3f" % (s, T(s)))


def G(lam, nphase, nq):
    """lam^{3/2} max_theta |R_D(lam, 0, theta)| over nphase normal angles in one period."""
    sq = 2 * np.pi * np.arange(nq) / nq
    _, zq, pq, _ = curve(sq)
    speed = np.abs(pq) * (2 * np.pi / nq)
    thetas = (np.arange(nphase) + 0.5) * (2 * np.pi / m) / nphase
    hh, _, rr = support(thetas)
    best = 0.0
    for k, theta in enumerate(thetas):
        I = np.sum(np.exp(1j * lam * np.real(zq * np.exp(-1j * theta))) * speed)
        # central symmetry: h and varrho are pi-periodic
        lead = np.sqrt(2 * np.pi / lam) * 2 * np.sqrt(rr[k]) * np.cos(lam * hh[k] - np.pi / 4)
        best = max(best, abs(I - lead))
    return best, lam ** 1.5 * best


for lam, nphase, nq in ((rho, 720, 1 << 17), (1.6e5, 480, 1 << 21)):
    r, g = G(lam, nphase, nq)
    print("lambda = %12.4f  (%d normal angles, %d quadrature nodes):  max|R_D| = %.5f   lambda^1.5 max|R_D| = %.1f"
          % (lam, nphase, nq, r, g))
    sys.stdout.flush()

# C^5 size of the support function: Fourier series of h in the normal angle over one period
M = 4096
phi = 2 * np.pi / m * np.arange(M) / M
hg, _, _ = support(phi)
H = np.fft.rfft(hg) / M
k = m * np.arange(len(H))
for K in (42, 60):
    spec = H * (1j * k) ** 5
    spec[K + 1:] = 0
    print("sup|h^(5)| with %d harmonics: %.3e" % (K, np.abs(np.fft.irfft(spec * M, M)).max()))
