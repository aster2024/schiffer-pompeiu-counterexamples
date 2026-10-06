#!/usr/bin/env python3
"""Exact rational evaluation of the algebraic constants of the two decimal inputs.

Reads center_j40_frozen.json and reference.json from this directory and evaluates,
with fractions.Fraction (no rounding), the finite sums that are quoted with
explicit bounds in Sections 3-7 of the paper.  Every assertion below is an
exact comparison of rational numbers.  The enclosures that involve Bessel
functions (spectral matrices, boundary traces, the chain of Table 2) are
produced by verify_planar_convex_v2.py, not by this script.
"""
from fractions import Fraction as F
import json
from pathlib import Path

here = Path(__file__).resolve().parent
centre = json.loads((here / "center_j40_frozen.json").read_text())
reference = json.loads((here / "reference.json").read_text())
m = centre["m"]
assert m == reference["m"] == 182
c = [F(x) for x in centre["p"]]            # psi_0 = sum c_j w^(jm+1)
c_ref = [F(x) for x in reference["p"]]
assert len(c) == 41 and len(c_ref) == 9
p = [(j * m + 1) * x for j, x in enumerate(c)]          # p_0 = psi_0'
p_ref = [(j * m + 1) * x for j, x in enumerate(c_ref)]


def show(name, value, digits=13):
    print("%-34s %.*g" % (name, digits, float(value)))


# Section 3: eight-mode reference, constants of the complement gap
rho_ref = p_ref[0]
s = sum(abs(x) for x in p_ref[1:])
v_star = 2 * rho_ref * s + s * s
g_star = 160 * (2 * rho_ref - 160)
gamma = g_star - v_star
show("rho_ref", rho_ref, 15)
show("v_*", v_star)
show("g_*", g_star)
show("gamma = g_* - v_*", gamma)
assert v_star < F("5004.03") and gamma > F("230022.34")
assert rho_ref - 160 > 1 and 1092 > rho_ref + 160

# Sections 4 and 6: forty-mode centre
rho = c[0]
assert F("814.457") < rho < F("814.458")
P = sum(abs(x) * 2 ** j for j, x in enumerate(p))
P1 = sum(abs(x) * (j * m) * 2 ** j for j, x in enumerate(p))
V_star = P * P - rho * rho
h = 1
while h * m <= rho:
    h += 1
d_h = (h * m) ** 2 - rho * rho
beta = V_star / d_h
show("rho", rho, 15)
show("P", P)
show("P_1", P1)
show("V_* = P^2 - rho^2", V_star)
print("%-34s %d   W_l = %d" % ("h", h, 1 + 2 * sum(2 ** j for j in range(1, h))))
show("d_h", d_h)
show("beta = V_*/d_h", beta)
assert h == 5
assert P < F("820.7003") and P1 < F("1176.5538") and V_star < F("10208.097")
assert d_h > F("164759.13") and beta < F("0.061958")
E0_upper = F("9.643e-61")                  # enclosed bound for E_0 (Table 2 of the paper)
d_star = rho * rho - V_star - E0_upper
j_star = rho - (P - rho)
show("d_* (using E_0 < 9.643e-61)", d_star)
show("j_*", j_star)
assert d_star > F("6.5313e5") and j_star > F("808.21")

# Section 3.4: transfer from the reference to the centre
unweighted = sum(abs(x) for x in p)
unweighted_ref = sum(abs(x) for x in p_ref)
change = sum(abs((p[j] if j < len(p) else 0) - (p_ref[j] if j < len(p_ref) else 0))
             for j in range(max(len(p), len(p_ref))))
delta_q = (unweighted + unweighted_ref) * change
show("delta_q", delta_q)
assert delta_q < F("8.376e-14")

# Section 5: strip constants
S = sum(abs(x) * 4 ** j for j, x in enumerate(c) if j)
show("S", S)
assert S < F("0.07") and rho > 2 * S

# Section 7: geometry on the ball of radius r_* = 1e-32
r = F("1e-32")
first = sum(abs(x) for x in p[1:])
second = sum(abs(x) * (j * m) for j, x in enumerate(p) if j)
L = rho - first - r
Theta = second + F(m + 1, 2) * r
show("L_*", L, 16)
show("Theta_*", Theta, 16)
show("1 - Theta_*/L_*", 1 - Theta / L, 16)
show("additive convexity margin", rho - first - second - F(m + 1, 2) * r, 16)
show("|c_1| - r_*/(2(m+1))", abs(c[1]) - r / (2 * (m + 1)), 16)
show("rho - (P - rho) - r_*", rho - (P - rho) - r, 16)
assert L > 811 and 1 - Theta / L > F("0.2998879")
assert rho - first - second - F(m + 1, 2) * r > F("243.3264")
assert abs(c[1]) - r / (2 * (m + 1)) > F("0.01645")
assert rho - (P - rho) - r > 808 and P + r < 821
print("all exact comparisons hold")
