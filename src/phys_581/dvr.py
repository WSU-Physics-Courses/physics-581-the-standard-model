import itertools
import math
from warnings import warn

import numpy as np
import scipy as sp


from . import bessel

_TINY = np.finfo(float).tiny

__all__ = [
    "roots_genlaguerre_radau",
    "LaguerreDVR",
    "LaguerreRadauDVR",
    "SphericalDVRBasis",
]


def roots_genlaguerre_radau(n, alpha):
    """Return the generalize Laguerre-Radau quadrature points and weights (x, w).

    This includes the fixed endpoint `x[0]=0`.
    """
    x, _w = sp.special.roots_genlaguerre(n, alpha + 1)
    Lx = sp.special.eval_genlaguerre(n, alpha, x)
    if int(alpha) == alpha:
        # Specialize for integer alpha using exact factorial function.
        alpha, n = int(alpha), int(n)
        gamma = lambda n: sp.special.factorial(n - 1, exact=True)
    else:
        gamma = sp.special.gamma

    comb = sp.special.comb

    w = np.concatenate(
        [
            [gamma(alpha + 1) / comb(n + alpha + 1, n)],
            gamma(alpha + 1) / (n + alpha + 1) * comb(n + alpha, n) / Lx**2,
        ]
    )

    x = np.concatenate([[0], x])
    return (x, w)


class LaguerreDVR:
    N = 10
    alpha = 0

    def __init__(self, **kw):
        for key in kw:
            if not hasattr(self, key):
                raise ValueError(f"Unknown {key=}")
            setattr(self, key, kw[key])
        self.init()

    def init(self):
        if int(self.alpha) == self.alpha:
            # Specialize for integer alpha using exact factorial function.
            self.alpha = int(self.alpha)
            self._gamma = lambda n: sp.special.factorial(n - 1, exact=True)
        else:
            self._gamma = sp.special.gamma

        N, a = self.N, self.alpha
        r, w = sp.special.roots_genlaguerre(N, a)
        if np.any(w == 0):
            warn(f"{sum(w==0)} zero weights (out of {N=})")

        n = np.arange(N)

        # Compute the weights by evaluating the polynomials at the abscissa.
        # Our P's are real, so no conj() here.
        C = np.sqrt(w[:, None]) * self.P(r=r[:, None], n=n[None, :])

        self.r, self.w, self.C = r, w, C

        # Compute the kinetic energy (negative laplacian / 2... no mass)
        du_w = self._du_w(r=r[:, None], m=n[None, :])
        self.T = np.einsum("i,im,in->mn", w, du_w, du_w) / 2

    def get_w(self, r):
        """Return the integration weight."""
        return r**self.alpha * np.exp(-r)

    def P(self, r, n, d=0):
        """Return the orthonormal function `P_n(r)`."""
        gamma = self._gamma
        normalization = np.sqrt(
            np.asarray(gamma(n + self.alpha + 1) / gamma(n + 1), dtype=float)
        )
        if d == 0:
            res = sp.special.eval_genlaguerre(n, self.alpha, r)
        elif d == 1:
            res = -sp.special.eval_genlaguerre(n - 1, self.alpha + 1, r)
        elif d == 2:
            res = sp.special.eval_genlaguerre(n - 2, self.alpha + 2, r)
        else:
            raise NotImplementedError(f"{d=}")
        return res / normalization

    def F(self, r, m, d=0):
        """Return F_m(r).

        Note: these do not have any factors of the weight function.
        """
        # Could do this with broadcasting, but it is not faster and is more complicated.
        return sum(self.C[m, n] * self.P(r, n, d=d) for n in np.arange(self.N))

    def _du_w(self, r, m):
        """Return `(du_n/dr)/sqrt(w(r))`."""
        Fm = self.F(r, m, d=0)
        dFm = self.F(r, m, d=1)
        if self.alpha == 0:
            return -0.5 * Fm + dFm

        # Needs careful computation of Fm/r...
        return (self.alpha / r - 1) / 2 * Fm + dFm

    def u(self, r, m):
        """Return the coordinate basis functions $u_n(r) = \sqrt{w(r)}F_n(r)$."""
        return np.sqrt(self.get_w(r)) * self.F(r, m)


class LaguerreRadauDVR:
    N = 10
    alpha = 0

    def __init__(self, **kw):
        for key in kw:
            if not hasattr(self, key):
                raise ValueError(f"Unknown {key=}")
            setattr(self, key, kw[key])
        self.init()

    def init(self):
        if int(self.alpha) == self.alpha:
            # Specialize for integer alpha using exact factorial function.
            self.alpha = int(self.alpha)
            self._gamma = lambda n: sp.special.factorial(n - 1, exact=True)
        else:
            self._gamma = sp.special.gamma

        N, a = self.N, self.alpha
        r, w = roots_genlaguerre_radau(N, a)
        if np.any(w == 0):
            warn(f"{sum(w==0)} zero weights (out of {N=})")

        # Skip the constant P_0(x)
        n = np.arange(N + 1)

        # Compute the weights by evaluating the polynomials at the abscissa.
        # Our P's are real, so no conj() here.
        C = np.sqrt(w[:, None]) * self.P(r=r[:, None], n=n[None, :])

        self.r, self.w, self.C = r, w, C

        # Compute the kinetic energy (negative laplacian / 2... no mass)
        du_w = self._du_w(r=r[:, None], m=n[None, :])
        self.T = np.einsum("i,im,in->mn", w, du_w, du_w) / 2

    def get_w(self, r):
        """Return the integration weight."""
        return r**self.alpha * np.exp(-r)

    def P(self, r, n, d=0):
        """Return the orthonormal function `P_n(r)`."""
        gamma = self._gamma
        normalization = np.sqrt(
            np.asarray(gamma(n + self.alpha + 1) / gamma(n + 1), dtype=float)
        )
        if d == 0:
            res = sp.special.eval_genlaguerre(n, self.alpha, r)
        elif d == 1:
            res = -sp.special.eval_genlaguerre(n - 1, self.alpha + 1, r)
        elif d == 2:
            res = sp.special.eval_genlaguerre(n - 2, self.alpha + 2, r)
        else:
            raise NotImplementedError(f"{d=}")
        return res / normalization

    def F(self, r, m, d=0):
        """Return F_m(r).

        Note: these do not have any factors of the weight function.
        """
        # Could do this with broadcasting, but it is not faster and is more complicated.
        return sum(self.C[m, n] * self.P(r, n, d=d) for n in np.arange(self.N + 1))

    def _du_w(self, r, m):
        """Return `(du_n/dr)/sqrt(w(r))`."""
        Fm = self.F(r, m, d=0)
        dFm = self.F(r, m, d=1)
        if self.alpha == 0:
            return -0.5 * Fm + dFm

        # Needs careful computation of Fm/r...
        return (self.alpha / r - 1) / 2 * Fm + dFm

    def u(self, r, m):
        """Return the coordinate basis functions $u_n(r) = \sqrt{w(r)}F_n(r)$."""
        return np.sqrt(self.get_w(r)) * self.F(r, m)


class SphericalDVRBasis:
    """Spherical DVR basis.

    This represents d-dimensional problems with rotational symmetry.

    Parameters
    ----------
    d : int
        Dimension of space.  I.e. `d=2` for cylindrical coordinates, `d=3` for spherical.
    k_max : float
        Cutoff (wavenumber in momentum space).
    R : float
        Approximate radius of basis.  The maximum number of lattice points will be
        roughly `R * k_max / np.pi`.
    """

    d = 3
    k_max = 1.0
    R = 10.0

    def __init__(self, **kw):
        for key in kw:
            if not hasattr(self, key):
                raise ValueError(f"Unknown key `{key}`")
            setattr(self, key, kw[key])

    def get_weights(self, l=0, volume=False):
        """Return the integration weights.

        Arguments
        ---------
        l : int
            Angular momentum quantum number.
        volume : bool
            If `True` then `sum(f*w)` will be the spherical integral (including the
            factors of the area `S` of the sphere and `r**(d-1)`), otherwise it will
            be simply the integral of the radial function from 0.
        """
        nu = self.get_nu(l=l)
        rn = self.get_rn(l=l)

        if volume:
            d_2 = self.d / 2.0
            S = 2 * np.pi**d_2 / sp.special.gamma(d_2)
            metric = r ** (self.d - 1) * S
        else:
            metric = 1
        ws = [metric / self._F(_n, _rn, nu=nu) ** 2 for _n, _rn in enumerate(rn)]
        return ws
        # Quadrature weights
        w = 2.0 / (self.k_max * z * bessel.J(nu=nu, d=1)(z) ** 2)

    def get_N(self, nu=None, l=0):
        """Return the number of abscissa below R."""
        # For large n, the roots of the bessel function are approximately
        # z[n] ~ (n + 0.75 + 0.5*nu)*pi, so R = r_max ~ (N - 0.25 + 0.5*nu)*pi/k_max
        if nu is None:
            nu = self.get_nu(l=l)
        N = int(np.ceil(0.25 - 0.5 * nu + self.k_max * self.R / np.pi))
        return N

    def get_nu(self, l=0):
        """Return `nu = l + d/2 - 1` for the centrifugal term.

        Arguments
        ---------
        l : int
           Angular quantum number.
        """
        nu = l + self.d / 2 - 1
        return nu

    def get_K(self, l=0):
        r"""Return `K`, the DVR kinetic energy matrix for the radial function.

        This term effects the $-d^2/dr^2 + (\nu^2 - 1/4)/r^2$ term.

        Returns
        -------
        K : array
           Operates on radial wavefunctions
        """
        nu = self.get_nu(l=l)
        rn = self.get_rn(nu=nu)
        zn = self.k_max * rn
        n = np.arange(len(zn))
        i1 = (slice(None), None)
        i2 = (None, slice(None))

        # Quadrature weights
        w = 2.0 / (self.k_max * zn * bessel.J(nu=nu, d=1)(zn) ** 2)

        # DVR kinetic term for radial function:
        K = np.ma.divide(
            (-1.0) ** (n[i1] - n[i2]) * 8.0 * zn[i1] * zn[i2],
            (zn[i1] ** 2 - zn[i2] ** 2) ** 2,
        ).filled(0)
        K[n, n] = 1.0 / 3.0 * (1.0 + 2.0 * (nu**2 - 1.0) / zn**2)
        K *= self.k_max**2
        return K

    def get_rn(self, nu=None, l=0):
        r"""Return the abscissa."""
        if nu is None:
            nu = self.get_nu(l=l)
        N = self.get_N(nu=nu)
        return bessel.j_root(nu=nu, N=N) / self.k_max

    def _F(self, n, r, nu, d=0):
        r"""Return the dth derivative of the n'th basis function."""
        rn = self.get_rn(nu=nu)
        zn = self.k_max * rn
        z = self.k_max * r
        H = bessel.J_sqrt_pole(nu=nu, zn=zn, d=0)
        coeff = math.sqrt(2.0 * self.k_max) * (-1.0) ** (n + 1) / (1.0 + r / rn)
        if 0 == d:
            return coeff * H(z)
        elif 1 == d:
            dH = bessel.J_sqrt_pole(nu=nu, zn=zn, d=1)
            return coeff * (dH(z) - H(z) / (z + zn)) * self.k_max
        else:
            raise NotImplementedError

    def get_F(self, r, l=0):
        """Return a function that can extrapolate a radial
        wavefunction to a new set of abscissa `r`."""
        nu = self.nu(l=l)
        rn = self.get_rn(nu=nu)
        n = np.arange(rn.size)[:, None]

        # Here is the transform matrix
        _F = self._F(n, r, nu=nu) / self._F(n, rn.T, nu=nu)

        def F(u):
            return np.dot(u, _F)

        return F
