"""Proof gates that remain active under optimized Python."""


class ProofGateError(ValueError):
    pass


def require(condition, message="proof gate failed"):
    if not condition:
        raise ProofGateError(message)
