#!/usr/bin/env python3
"""Arb certificate of convexity for the exact R6 Schiffer meridian domain.

Adapted from the R^4 script paper/anc/verify_convex.py. Run after
verify_r6.py --stage final, which proves that the exact zero x*=(g*,c*) lies in
the X-ball of radius r=1/5000 about the frozen dyadic centre x°. Here the shape
coefficients are indexed by j = 1 (mod 4) and weighted by
    omega_j = alpha (j+2) rho^(j+1),  alpha = b^3/2,  b = c°_1,  rho = 11/10,
so every shape vector c in the certified ball satisfies
    sum_j omega_j |c_j - c°_j| <= r            (the infinite tail included).
The script checks a global, unsampled lower bound for Re(1 + w psi''/psi') on
the closed unit disc that holds simultaneously for every such c. No boundary
grid and no truncation of the tail enters the certificate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from flint import arb, ctx

# verify_r6.py with the fail-closed max_upper, frozen for this release.
EXPECTED_VERIFIER_SHA256 = "67a565438b6a7e459145a18b40c8c12d21134b097280575fce4bc9db8b6941c9"
EXPECTED_CENTRE_SHA256 = "83018f0975e0cd61f4168e2de8bb1946032265869b5bfa168773c5920198aaba"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def exact_dyadic(hex_string):
    numerator, denominator = float.fromhex(hex_string).as_integer_ratio()
    value = arb(numerator) / denominator
    require(value.is_exact(), "dyadic centre coefficient lost exactness")
    return value


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(x, what):
    require(x.is_finite() and x.upper().is_finite() and x.lower().is_finite(),
            "nonfinite interval: " + what)
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=128)
    parser.add_argument("--main-certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.bits >= 96, "at least 96 bits are required")
    ctx.prec = args.bits
    ctx.threads = 1

    here = Path(__file__).resolve().parent
    source_path = here / "verify_r6.py"
    data_path = here / "center_r6.json"
    cert_path = args.main_certificate or here / "certificate_r6.json"
    output_path = args.output or here / "convex_r6.json"

    require(sha256(source_path) == EXPECTED_VERIFIER_SHA256,
            "main verifier source is not the released fail-closed version")
    require(sha256(data_path) == EXPECTED_CENTRE_SHA256, "centre data changed")

    cert = json.loads(cert_path.read_text())
    require(cert.get("status") == "PROVED" and cert.get("bits", 0) >= 96,
            "main interval certificate has not proved the radius bound")
    require(cert.get("M") == 58 and cert.get("S") == 30,
            "main certificate uses a different centre dimension")
    require(cert.get("rho_weight") == "11/10", "main certificate uses a different weight")
    require(cert.get("verifier_sha256") == EXPECTED_VERIFIER_SHA256,
            "main certificate was produced by a different verifier")
    require(cert.get("sha256", {}).get("center_r6.json") == EXPECTED_CENTRE_SHA256,
            "main certificate refers to a different centre")
    radius = arb(1) / 5000
    require(arb(cert["radius"]).contains(radius), "main certificate uses a different radius")

    data = json.loads(data_path.read_text())
    require(data["M"] == 58 and data["S"] == 30, "unexpected centre dimensions")
    js = list(range(1, data["M"], 4))
    require(len(data["c"]) == len(js) == 15 and len(data["g"]) == 450,
            "unexpected centre coefficient count")
    require(all(j % 4 == 1 for j in js), "shape index set is not j = 1 mod 4")
    c = {j: exact_dyadic(h) for j, h in zip(js, data["c"])}

    rho = arb(11) / 10
    b = c[1]
    require(b > 0 and rho > 1, "invalid centre or weight")
    alpha = b**3 / 2

    def omega(j):
        return alpha * (j + 2) * rho**(j + 1)

    # Perturbation delta c = c - c° with sum_j omega_j |delta c_j| <= r
    # (j = 1, 5, 9, ...; infinitely many j). For |w| <= 1:
    #  |delta c_1| + sum_{j>=5} j |delta c_j| <= r * sup_j kappa_j,
    #    kappa_1 = 1/omega_1 = 1/(3 alpha rho^2),
    #    kappa_j = j/omega_j < 1/(alpha rho^(j+1)) <= 1/(alpha rho^6)   (j >= 5);
    #  sum_{j>=5} j(j-1) |delta c_j| <= r * sup_{j>=5} j(j-1)/omega_j
    #    <= r * sup_{x>0} x rho^(-x-1)/alpha = r/(alpha rho e log rho),
    #  because j(j-1)/(j+2) < j and x rho^(-x) attains its maximum 1/(e log rho)
    #  at x = 1/log rho.  These bounds cover the entire infinite tail.
    kappa1 = 1 / omega(1)
    kappa_rest = 1 / (alpha * rho**6)
    first_perturbation = finite(radius * max(kappa1.upper(), kappa_rest.upper()), "E1")
    second_perturbation = finite(radius / (alpha * rho * arb.const_e() * rho.log()), "E2")

    first_centre = finite(sum((j * abs(v) for j, v in c.items() if j >= 5), arb(0)), "A1")
    second_centre = finite(sum((j * (j - 1) * abs(v) for j, v in c.items() if j >= 5), arb(0)), "A2")
    derivative_lower = finite(b - first_centre - first_perturbation, "L")
    second_upper = finite(second_centre + second_perturbation, "U")
    curvature_lower = finite(1 - second_upper / derivative_lower, "1-U/L")

    require(derivative_lower > 0, "psi' may vanish on the closed disk")
    require(curvature_lower > 0, "Re(1+w psi''/psi') is not certified positive")
    # A rational bound with 8 decimals, itself certified: q < lower endpoint.
    k = int(float(curvature_lower.mid()) * 10**8) - 1
    require(0 < k < 10**8, "unexpected size of the convexity margin")
    rational_lower = arb(k) / 10**8
    require(curvature_lower > rational_lower, "rational convexity bound not certified")
    rational_text = f"{k // 10**8}.{k % 10**8:08d}"
    result = {
        "status": "PROVED_CONVEX",
        "bits": args.bits,
        "verifier_sha256": sha256(source_path),
        "centre_sha256": sha256(data_path),
        "main_certificate_sha256": sha256(cert_path),
        "radius": str(radius),
        "rho": str(rho),
        "centre_b": str(b),
        "alpha": str(alpha),
        "first_centre_absolute_sum": str(first_centre),
        "second_centre_absolute_sum": str(second_centre),
        "first_perturbation_bound": str(first_perturbation),
        "second_perturbation_bound": str(second_perturbation),
        "psi_prime_modulus_lower": str(derivative_lower),
        "psi_double_prime_modulus_upper": str(second_upper),
        "real_1_plus_w_psi2_over_psi1_lower": str(curvature_lower),
        "certified_rational_lower_bound": rational_text,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    print("PROVED_CONVEX  Re(1+w psi''/psi') >", rational_text,
          "on |w|<=1 for every shape in the certified ball", flush=True)
    print("certificate:", output_path, flush=True)


if __name__ == "__main__":
    main()
