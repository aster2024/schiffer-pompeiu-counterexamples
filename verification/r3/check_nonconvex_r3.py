#!/usr/bin/env python3
"""Exact check of the non-convexity of the three-dimensional domain at its poles.

This script is separate from verify_r3.py and does not modify it. It reads the
frozen dyadic centre c° embedded in verify_r3.py (after checking the SHA-256 of
the verifier, of its payload and of the centre file) and verifies, in exact
rational arithmetic, the inequalities used in the paper for Theorem 1.1(iv):

  S1 = sum_j j c°_j        = psi°'(1),
  S2 = sum_j j^2 c°_j      = psi°'(1) + psi°''(1),

and, for every c with sum_j j R^j |c_j - c°_j| <= E = r/(3b)
(R = 21/20, r = 3e-7, b the fixed dyadic scale of the paper),

  psi_c'(1)               >= S1 - E/R   > 0,
  psi_c'(1) + psi_c''(1)  <= S2 + K E   < 0,     K = 8 >= max_j j R^(-j).

Hence 1 + psi_c''(1)/psi_c'(1) < 0 on the whole existence ball.
Only the Python standard library is used. Exit code 0 and the last line
NONCONVEX_AT_POLES_VERIFIED on success; exit code 2 otherwise.
"""
import base64
import hashlib
import json
import re
import sys
import zlib
from fractions import Fraction
from pathlib import Path

VERIFIER_SHA256 = "c2802ee1202ac2dd899f01c4e5794f1a7e5f78872f53becd1728f23e19180cdb"
PAYLOAD_SHA256 = "0ba7a364b8bed66c4e50531f871762b8a185968a9bb38b59350c4cc9ad749a4d"
CENTRE_NAME = "refined_M180_step1.json"
CENTRE_SHA256 = "1c2e959714280ca4bbd2e60c5a7bc378486f5c96205d858e0de0c561ef26a344"


def need(ok, message):
    if not ok:
        print("NOT_VERIFIED", message)
        sys.exit(2)


def main():
    here = Path(__file__).resolve().parent
    raw = (here / "verify_r3.py").read_bytes()
    need(hashlib.sha256(raw).hexdigest() == VERIFIER_SHA256, "verify_r3.py checksum")
    match = re.search(r'PACKED = """(.*?)"""', raw.decode("ascii"), re.S)
    need(match is not None, "embedded payload not found")
    data = zlib.decompress(base64.b85decode("".join(match.group(1).split())))
    need(hashlib.sha256(data).hexdigest() == PAYLOAD_SHA256, "payload checksum")
    payload = json.loads(data)
    centre = base64.b64decode(payload["files"][CENTRE_NAME])
    need(hashlib.sha256(centre).hexdigest() == CENTRE_SHA256, "centre checksum")
    need(payload["sha256"][CENTRE_NAME] == CENTRE_SHA256, "centre checksum in payload")
    coeffs = json.loads(centre)["c"]

    c = {}
    for key, (mantissa, exponent) in coeffs.items():
        j = int(key)
        need(j % 2 == 1 and 1 <= j <= 181, "unexpected shape index %r" % key)
        c[j] = Fraction(int(mantissa)) * Fraction(2) ** int(exponent)
    need(sorted(c) == list(range(1, 182, 2)), "shape indices must be 1, 3, ..., 181")

    R = Fraction(21, 20)
    r = Fraction(3, 10**7)
    b = Fraction(float.fromhex("0x1.85067a9afca2dp+5"))
    E = r / (3 * b)
    K = Fraction(8)
    # j R^(-j) <= K for all odd j >= 1: check j <= 21 directly; for j >= 21 the
    # ratio of consecutive odd terms is (j+2)/(j R^2) <= 23/(21 R^2) < 1.
    need(all(Fraction(j) / R**j <= K for j in range(1, 23, 2)), "j R^-j <= K for j <= 21")
    need(Fraction(23, 21) < R * R, "monotonicity of j R^-j for j >= 21")

    S0 = sum(c.values())
    S1 = sum(j * v for j, v in c.items())
    S2 = sum(j * j * v for j, v in c.items())
    lower_d1 = S1 - E / R
    upper_d12 = S2 + K * E
    need(lower_d1 > 0, "psi'(1) > 0 on the ball")
    need(upper_d12 < 0, "psi'(1) + psi''(1) < 0 on the ball")
    ratio_upper = upper_d12 / (S1 + E / R)  # 1 + psi''/psi' = (psi'+psi'')/psi' < this
    need(ratio_upper < 0, "negative curvature quantity")

    def dec(x):
        return format(float(x), ".15g")

    print("psi_centre(1)            =", dec(S0))
    print("psi_centre'(1)           =", dec(S1))
    print("psi_centre'(1)+psi''(1)  =", dec(S2))
    print("1+psi''(1)/psi'(1) at c° =", dec(1 + (S2 - S1) / S1))
    print("E = r/(3b)               =", dec(E))
    print("ball: psi'(1)           >=", dec(lower_d1))
    print("ball: psi'(1)+psi''(1)  <=", dec(upper_d12))
    print("ball: 1+psi''/psi'(1)   < ", dec(ratio_upper))
    # exact two-sided enclosures used in the paper
    need(Fraction(4281903, 100000) < S1 < Fraction(4281904, 100000), "S1 enclosure")
    need(Fraction(-7758015, 1000000) < S2 < Fraction(-7758014, 1000000), "S2 enclosure")
    need(upper_d12 < Fraction(-7758, 1000), "S2 + K E < -7.758")
    need(ratio_upper < Fraction(-181, 1000), "ratio < -0.181")
    print("NONCONVEX_AT_POLES_VERIFIED")


if __name__ == "__main__":
    main()
