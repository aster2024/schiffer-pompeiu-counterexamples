#!/usr/bin/env python3
"""Descriptive geometry of the forty-mode centre psi_0 (binary64 / mpmath, no enclosures).

Reads ../center_j40_frozen.json.  Prints the radius range, the curvature range,
the convexity factor, the width ratio, the relative Hausdorff distance to
centred discs, and the two numbers used in the comparison with Kobayashi's
perturbation theorem (Section 8 of the paper).  None of these numbers enters
the existence proof.
"""
import json
from pathlib import Path

import numpy as np
from mpmath import mp, besselj, findroot, pi as mp_pi, e as mp_e, mpf

HERE = Path(__file__).resolve().parent
centre = json.loads((HERE.parent / "center_j40_frozen.json").read_text())
m = int(centre["m"])
c = np.array([float(x) for x in centre["p"]])
deg = m * np.arange(len(c)) + 1
rho = c[0]
print("m =", m, " modes =", len(c) - 1)
print("rho = c_0 = %.13f   rho^2 = %.6f" % (rho, rho * rho))

N = 1_000_000                       # points on one symmetry period
t = (np.arange(N) + 0.5) * (2 * np.pi / m) / N
z = np.zeros(N, complex)
p = np.zeros(N, complex)
q = np.zeros(N, complex)            # w * psi''(w)
for j in range(len(c)):
    wd = np.exp(1j * deg[j] * t)    # w^(jm+1)
    z += c[j] * wd
    p += deg[j] * c[j] * wd * np.exp(-1j * t)
    q += deg[j] * (deg[j] - 1) * c[j] * wd * np.exp(-1j * t)
radius = np.abs(z)
factor = np.real(1 + q / p)
kappa = factor / np.abs(p)
print("|psi_0| on the circle: min %.10f  max %.10f" % (radius.min(), radius.max()))
print("peak-to-trough radial variation / rho = %.5e" % ((radius.max() - radius.min()) / rho))
print("Re(1 + w psi''/psi') on the circle: min %.7f  max %.7f" % (factor.min(), factor.max()))
print("curvature * rho: min %.6f  max %.6f" % (kappa.min() * rho, kappa.max() * rho))
print("|psi_0'| on the circle: min %.4f  max %.4f" % (np.abs(p).min(), np.abs(p).max()))

# support function h in the normal angle: outer normal direction = arg(w psi'(w))
normal = np.exp(1j * t) * p / np.abs(p)
h = np.real(z * np.conj(normal))
print("support function: min %.10f  max %.10f" % (h.min(), h.max()))
print("minimum width / maximum width = min h / max h = %.10f   (central symmetry: width = 2h)" % (h.min() / h.max()))
dH = (h.max() - h.min()) / (h.max() + h.min())
print("relative Hausdorff distance to centred discs, (max h - min h)/(max h + min h) = %.6e" % dH)

# Kobayashi comparison (see Section 8).  Normalize to conformal radius one: the radial
# function g of the normalized centre satisfies sup|g - 1| >= dH, and the zero circle of
# the Fourier transform has radius rho.  The first smallness condition of
# [Kobayashi 1993, Proposition 4.4] requires
#     t_0 B(g) exp(||g||_inf) < eps(2,R) delta(2,R),      R >= rho,
# whose left side is >= 2 pi e sup|g(t_0,.) - 1| and whose right side is
# <= 2 pi |J_0(j_{1,k})| / j_{1,k} for every zero j_{1,k} <= R of J_1.
mp.dps = 30
gmax = max(radius.max() / rho - 1, 1 - radius.min() / rho)
print("normalized radial function: sup|g - 1| = %.6e" % gmax)
print("2 pi e * dH                = %.4e   (lower bound for the left side, any dilate)" % float(2 * mp_pi * mp_e * dH))
zero = findroot(lambda x: besselj(1, x), mpf(rho))        # zero of J_1 nearest to rho
zeros = [zero - mp_pi, zero]                               # refine the neighbouring zero below
zeros[0] = findroot(lambda x: besselj(1, x), zeros[0])
for j1 in zeros:
    k = int(round(float(j1 / mp_pi - mpf(1) / 4)))        # McMahon: j_{1,k} ~ (k + 1/4) pi
    print("zero j_{1,%d} of J_1 = %s  (%s rho):  2 pi |J_0(j)| / j = %.4e"
          % (k, mp.nstr(j1, 16), "below" if j1 < mpf(centre["p"][0]) else "above",
             float(2 * mp_pi * abs(besselj(0, j1)) / j1)))
