"""Descriptive centre plot; no interval-certified claim is inferred from it."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent
centre = json.loads((BASE.parent.parent / "anc/r2/center_j40_frozen.json").read_text())
m = int(centre["m"])
c = np.array([float(x) for x in centre["p"]])
degree = m * np.arange(len(c)) + 1
rho = c[0]

def boundary(theta):
    w = np.exp(1j * theta)
    z = np.sum(c[:, None] * w[None, :] ** degree[:, None], axis=0)
    p = np.sum((degree * c)[:, None] *
               w[None, :] ** (degree[:, None] - 1), axis=0)
    wpprime = np.sum((degree * (degree - 1) * c)[:, None] *
                    w[None, :] ** (degree[:, None] - 1), axis=0)
    factor = np.real(1 + wpprime / p)
    curvature = factor / np.abs(p)
    return z, factor, curvature

theta = np.linspace(0, 2*np.pi, 32769)
z, _, _ = boundary(theta)
angle = np.angle(z)
radius = np.abs(z)
schematic = (1 + 2000*(radius/rho - 1)) * np.exp(1j*angle)
phi = np.linspace(0, 2*np.pi, 2049)
_, factor, curvature = boundary(phi/m)
plt.rcParams.update({"font.family": "serif", "font.size": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout="constrained",
                         gridspec_kw={"width_ratios": [1, 1.35]})
axes[0].plot(schematic.real, schematic.imag, color="#173e69", lw=.65)
axes[0].plot(np.cos(theta), np.sin(theta), color="#b5b5b5", lw=.7)
axes[0].set_aspect("equal")
axes[0].set_axis_off()
axes[0].set_title(r"182 ripples; radial deviation $\times 2000$", fontsize=10)
axes[1].plot(phi/(2*np.pi), rho*curvature, color="#173e69", lw=1.5)
axes[1].axhline(1, color="#999999", linestyle="--", lw=.7)
axes[1].set(xlabel=r"Fraction of one symmetry period",
            ylabel=r"Curvature ratio $\rho\kappa$", xlim=(0, 1), ylim=(.2, 1.85))
axes[1].set_xticks([0, .25, .5, .75, 1])
axes[1].grid(axis="y", color="#dddddd", lw=.5)
fig.savefig(BASE/"domain.pdf")
fig.savefig(BASE/"domain.png", dpi=180)
diagnostic = {"scope": "binary64 descriptive centre evaluation, not a certificate",
              "rho": rho, "radial_min_sample": float(radius.min()),
              "radial_max_sample": float(radius.max()),
              "relative_peak_to_trough": float(np.ptp(radius)/rho),
              "curvature_ratio_min_sample": float((rho*curvature).min()),
              "curvature_ratio_max_sample": float((rho*curvature).max()),
              "convexity_factor_min_sample": float(factor.min()),
              "schematic_amplification": 2000}
(BASE/"DIAGNOSTIC.json").write_text(json.dumps(diagnostic, indent=2)+"\n")
