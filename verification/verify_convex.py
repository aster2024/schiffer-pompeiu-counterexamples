#!/usr/bin/env python3
"""Arb certificate of convexity for the exact Schiffer domain.

Run after verify_r4.py. The latter proves that the exact coefficient vector
lies in the X-ball of radius 1/50000 about the dyadic centre. This script
checks a global, unsampled lower bound for Re(1+w*psi''/psi') on |w| <= 1.
The perturbation estimates include every odd shape coefficient, including
the infinite tail. No numerical boundary grid enters the certificate.
"""

import argparse
import hashlib
import json
import runpy
from pathlib import Path

from flint import arb, ctx

EXPECTED_SOURCE_SHA256 = "32c4334a215df722651b610f3c2922b3472b65edab2be21d32df7739195f85ed"


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
    source_path = here / "verify_r4.py"
    data_path = here / "data" / "centre.json"
    cert_path = args.main_certificate or here / "logs" / "certificate_r4_128.json"
    output_path = args.output or here / "logs" / "convex_128.json"
    require(sha256(source_path) == EXPECTED_SOURCE_SHA256,
            "main verifier source is not the audited current version")
    source_data = runpy.run_path(str(source_path))["CENTRE"]
    data = json.loads(data_path.read_text())
    require(data == source_data, "centre.json differs from embedded source data")
    require(data["M"] == 41 and data["S"] == 24, "unexpected centre dimensions")
    require(len(data["c"]) == 21 and len(data["g"]) == 504,
            "unexpected centre coefficient count")

    cert = json.loads(cert_path.read_text())
    require(cert.get("status") == "PROVED" and cert.get("bits", 0) >= 96,
            "main interval certificate has not proved the radius bound")
    require(cert.get("M") == 41 and cert.get("S") == 24,
            "main certificate uses a different centre dimension")
    radius = arb(1) / 50000
    require(arb(cert["radius"]).contains(radius),
            "main certificate uses a different radius")

    require(data["rho"] == [21, 20], "unexpected coefficient weight")
    rho = arb(data["rho"][0]) / data["rho"][1]
    c = {j: exact_dyadic(h) for j, h in zip(range(1, 42, 2), data["c"])}
    b = c[1]
    require(b > 0 and rho > 1, "invalid coefficient weight")

    # From ||c*-c°||_X <= r, with omega_j=b^2(j+1)rho^j:
    # sum j|delta c_j| <= r/(b^2 rho).
    # For the second derivative, j(j-1)/(j+1) <= j and
    # sup_{j odd >= 1} j rho^{-j} <= sum_{j >= 1} j rho^{-j}
    # = rho/(rho-1)^2. The last bound covers the entire infinite tail.
    first_centre = sum((j * abs(v) for j, v in c.items() if j >= 3), arb(0))
    second_centre = sum((j * (j - 1) * abs(v)
                         for j, v in c.items() if j >= 3), arb(0))
    first_perturbation = radius / (b * b * rho)
    second_perturbation = (radius / (b * b)) * rho / (rho - 1)**2
    derivative_lower = b - first_centre - first_perturbation
    second_upper = second_centre + second_perturbation
    curvature_lower = 1 - second_upper / derivative_lower

    require(derivative_lower > 0, "psi' may vanish on the closed disk")
    require(curvature_lower > 0,
            "Re(1+w psi''/psi') is not certified positive")
    result = {
        "status": "PROVED_CONVEX",
        "bits": args.bits,
        "source_sha256": sha256(source_path),
        "data_sha256": sha256(data_path),
        "main_certificate_sha256": sha256(cert_path),
        "radius": str(radius),
        "rho": str(rho),
        "centre_b": str(b),
        "first_centre_absolute_sum": str(first_centre),
        "second_centre_absolute_sum": str(second_centre),
        "first_perturbation_bound": str(first_perturbation),
        "second_perturbation_bound": str(second_perturbation),
        "psi_prime_modulus_lower": str(derivative_lower),
        "psi_double_prime_modulus_upper": str(second_upper),
        "real_1_plus_w_psi2_over_psi1_lower": str(curvature_lower),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    print("PROVED_CONVEX", curvature_lower, flush=True)


if __name__ == "__main__":
    main()
