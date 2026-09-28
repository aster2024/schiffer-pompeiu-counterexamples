"""Conformal fixed-disc formulation for the I2(m) Weyl Schiffer problem.

Planar problem, wavenumber normalised to 1:
    F Weyl anti-invariant, (Delta+1)F=0, F-pi=|grad(F-pi)|=0 on dD.

Fixed-disc form: psi(w)=sum c_j w^j, j=1 mod m, real c_j;
    W = F o psi - Im(psi**m), g = Laplacian W, W = K g,
    H(g, psi) = g + |psi'|^2 (K g + Im(psi**m)) = 0.

Basis: S_{n,s} = r^n P_s^{(0,n)}(2 r^2 - 1) sin(n theta), n=0 mod m, s >= 0.
  <S_{m,s}, S_{m,s}>_{L^2(disc)} = pi / (2 (m + 2 s + 1)).
  K S_{m,s} = S_{m,s-1}/(4D(D+1)) - S_{m,s}/(2D(D+2)) + S_{m,s+1}/(4(D+1)(D+2)),  D = m + 2 s, s >= 1.
Shape: psi(w) = sum_{j=1 mod m} c_j w^j, c_j real.
"""
import numpy as np
from scipy import special as sps


def jacobi_table(m, smax, x):
    """P_s^{(0,m)}(x), s = 0..smax (rows), via the standard three-term recurrence (alpha=0, beta=m)."""
    a, b = 0.0, float(m)
    P = np.zeros((smax + 1, len(x)))
    P[0] = 1.0
    if smax >= 1:
        P[1] = (a + 1) + (a + b + 2) * (x - 1) / 2
    for n in range(2, smax + 1):
        c1 = 2 * n * (n + a + b) * (2 * n + a + b - 2)
        c2 = (2 * n + a + b - 1) * ((2 * n + a + b) * (2 * n + a + b - 2) * x + a * a - b * b)
        c3 = 2 * (n + a - 1) * (n + b - 1) * (2 * n + a + b)
        P[n] = (c2 * P[n - 1] - c3 * P[n - 2]) / c1
    return P


class Disc:
    """Quadrature grid on a pi/m sector and basis tables."""

    def __init__(self, Weyl_m, M, S, Nx=None, Nq=None, extra_s=1, step=None,
                 boundary_scale=1.0):
        if step is None:
            step = Weyl_m
        if Weyl_m not in (3, 4, 6) or step not in (Weyl_m, 2*Weyl_m) or M % step != Weyl_m % step:
            raise ValueError("M must match the first sine frequency modulo step")
        self.Weyl_m, self.step, self.M, self.S = Weyl_m, step, M, S
        self.boundary_scale = float(boundary_scale)
        if not np.isfinite(self.boundary_scale) or self.boundary_scale<=0:
            raise ValueError("boundary scale must be positive and finite")
        self.ms = np.arange(Weyl_m, M + 1, step)
        self.js = np.arange(1, 1 + step*len(self.ms), step)
        Jc = int(self.js[-1])
        self.Jc = Jc
        dF = (M + 2 * S + 2) + 2 * (Jc - 1) + 2
        if Nq is None:
            Nq = (dF + M) // Weyl_m + 8
        if Nx is None:
            Nx = ((dF + M) // 2 + S) // 2 + 8
        self.Nx, self.Nq = Nx, Nq
        x, wx = np.polynomial.legendre.leggauss(Nx)
        self.x, self.wx = x, wx
        self.r = np.sqrt((1 + x) / 2)
        self.th = (np.arange(Nq) + 0.5) * (np.pi / Weyl_m) / Nq
        R, TH = np.meshgrid(self.r, self.th, indexing='ij')
        self.R, self.TH = R.ravel(), TH.ravel()
        self.W = self.R * np.exp(1j * self.TH)
        # Full disc = 2m copies of this sector; r dr=dx/4.
        self.qw = (np.pi / (2 * Nq)) * np.repeat(wx, Nq)
        self.S_ext = S + extra_s
        # basis tables: dict m -> array (S_ext+1, grid)
        self.tab = {}
        for m in self.ms:
            P = jacobi_table(m, self.S_ext, x)                  # (s, Nx)
            rad = (self.r[None, :] ** m) * P                     # (s, Nx)
            ang = np.sin(m * self.th)                            # (Nq,)
            self.tab[m] = (rad[:, :, None] * ang[None, None, :]).reshape(self.S_ext + 1, -1)
        # index maps
        self.gidx = [(m, s) for m in self.ms for s in range(1, S + 1)]
        self.hidx = [(m, 0) for m in self.ms]
        self.ng, self.nc = len(self.gidx), len(self.js)
        self.rows = self.gidx + self.hidx
        self.norm2 = {}

    # --- coefficient <-> grid ---------------------------------------------------------------
    def eval_S(self, coef):
        """coef: dict or array over (m, s) with s in 0..S_ext -> grid values."""
        out = np.zeros(self.R.shape)
        for m in self.ms:
            out += coef[m] @ self.tab[m]
        return out

    def project(self, f, smax=None):
        """L2 projection of grid function f (odd-odd sector) onto S_{m,s}, s=0..smax -> dict m -> array."""
        if smax is None:
            smax = self.S
        out = {}
        fw = f * self.qw
        for m in self.ms:
            n2 = np.pi / (2 * (m + 2 * np.arange(smax + 1) + 1))
            out[m] = (self.tab[m][:smax + 1] @ fw) / n2
        return out

    def K_coef(self, gc):
        """gc: dict m -> array over s=1..S (index s-1). Returns dict m -> array over s=0..S+1."""
        out = {}
        for m in self.ms:
            v = np.zeros(self.S + 2)
            for s in range(1, self.S + 1):
                D = m + 2 * s
                a = gc[m][s - 1]
                v[s - 1] += a / (4 * D * (D + 1))
                v[s] += -a / (2 * D * (D + 2))
                v[s + 1] += a / (4 * (D + 1) * (D + 2))
            out[m] = v
        return out

    # --- unknown vector ------------------------------------------------------------------------
    def unpack(self, z):
        gz = z[:self.ng].reshape(len(self.ms), self.S)
        gc = {m: gz[i] for i, m in enumerate(self.ms)}
        c = z[self.ng:]
        return gc, c

    def pack(self, gc, c):
        return np.concatenate([np.concatenate([gc[m] for m in self.ms]), c])

    def psi_vals(self, c):
        W = self.W
        psi = np.zeros(W.shape, complex)
        dpsi = np.zeros(W.shape, complex)
        for j, cj in zip(self.js, c):
            psi += cj * W ** j
            dpsi += cj * j * W ** (j - 1)
        return psi, dpsi

    def residual_grid(self, z):
        gc, c = self.unpack(z)
        gfull = {m: np.concatenate([[0.0], gc[m], np.zeros(self.S_ext - self.S)]) for m in self.ms}
        g = self.eval_S(gfull)
        Kc = self.K_coef(gc)
        Kg = self.eval_S({m: np.concatenate([Kc[m], np.zeros(self.S_ext + 1 - len(Kc[m]))]) for m in self.ms})
        psi, dpsi = self.psi_vals(c)
        a2 = np.abs(dpsi) ** 2
        V = Kg + (psi ** self.Weyl_m).imag / self.boundary_scale
        return g + a2 * V, dict(g=g, Kg=Kg, psi=psi, dpsi=dpsi, a2=a2, V=V)

    def residual(self, z):
        Fg, aux = self.residual_grid(z)
        pr = self.project(Fg)
        res = np.concatenate([pr[m][1:self.S + 1] for m in self.ms] + [np.array([pr[m][0] for m in self.ms])])
        return res, aux

    def jacobian(self, z):
        gc, c = self.unpack(z)
        Fg, aux = self.residual_grid(z)
        a2, V, dpsi, psi = aux['a2'], aux['V'], aux['dpsi'], aux['psi']
        cols = []
        # g columns
        for m in self.ms:
            for s in range(1, self.S + 1):
                D = m + 2 * s
                KS = (self.tab[m][s - 1] / (4 * D * (D + 1)) - self.tab[m][s] / (2 * D * (D + 2))
                      + self.tab[m][s + 1] / (4 * (D + 1) * (D + 2)))
                cols.append(self.tab[m][s] + a2 * KS)
        W = self.W
        for j in self.js:
            dW = j * W ** (j - 1)
            cols.append(2 * np.real(np.conj(dpsi) * dW) * V
                        + a2 * self.Weyl_m * np.imag(psi ** (self.Weyl_m - 1) * W ** j)
                        / self.boundary_scale)
        C = np.array(cols).T   # grid x ncols
        # project all columns
        Cw = C * self.qw[:, None]
        rows = []
        for m in self.ms:
            n2 = np.pi / (2 * (m + 2 * np.arange(self.S + 1) + 1))
            rows.append((self.tab[m][1:self.S + 1] @ Cw) / n2[1:, None])
        hr = []
        for m in self.ms:
            n2 = np.pi / (2 * (m + 1))
            hr.append((self.tab[m][0] @ Cw) / n2)
        return np.vstack(rows + [np.array(hr)])


def ball_rho(m, k):
    """k-th zero of J_(m+1), radial Neumann wavenumber in R^(2m+2)."""
    return sps.jn_zeros(m + 1, k)[-1]
