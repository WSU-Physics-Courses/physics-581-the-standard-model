"""Tools for solving the radial Schrodinger equation to high precision,"""
from functools import partial

import numpy as np

import scipy.optimize
import scipy.integrate
import scipy as sp


class CoulombSEQ:
    """Helper to solve the radial Schrodinger equation with a Coulomb tail.

    Assumes that the potential behaves like `V(r) = -alpha/r` at long distances, and
    that `V(r)` is not more singular than the centrifugal term `1/r^2` at short
    distance.

    Attributes
    ----------
    dimension : int
        Dimension of space.  I.e `dimension = 2` for cylindrical coordinates (without z).
    hbar, m : float
        Physical constants
    alpha : float
        Coefficient for the tail of the potential (Units of energy * distance).
    """

    ######################################################################
    # Main interface and user-define methods
    dimension = 3
    hbar = 1.0
    m = 1.0
    alpha = 1.0

    def V(self, r):
        """Potential.  (Overload to change.)

        The default is the pure Coulomb potential.  Custom potentials must include this
        tail.
        """
        return -self.alpha / r

    def compute_E(self, E, tol=1e-8, l=0, lam=0.999, **kw):
        """Return the energy eigenvalue E by shooting.

        Arguments
        ---------
        E : float
            Initial guess.
        l : int
            Angular momentum quantum number.
        tol : float
            Desired relative tolerance.
        lam : float
            Parameter for bracketing the energy.  Initial search is in the range
            `(lam*E, E/lam)`
        kw : dict
            Additional arguments passed down to `solve_ivp`.
        """
        f_ = partial(self.objective, tol=tol, **kw)
        f0 = f1 = 1
        while f0 * f1 > 0:
            lam *= lam
            E0 = lam * E
            E1 = E / lam
            f0, f1 = map(f_, [E0, E1])
        return sp.optimize.brentq(f_, E0, E1, xtol=tol)

    ######################################################################
    # Helpers etc.  Should not need modification in general.
    def __init__(self, **kw):
        for key in kw:
            if not hasattr(self, key):
                raise ValueError(f"Unknown key {key}")
            setattr(self, key, kw[key])

    def get_a(self, E):
        """Return a, the length scale associated with the exponential decay."""
        assert E < 0
        return self.hbar / np.sqrt(-2 * self.m * E)

    def get_r_v(self, E):
        """Return `r_v`, an approximate location of the turning point."""
        assert E < 0
        r_v = self.alpha / (-E)

        # Could actually find this by looking at the potential.
        return r_v

    def compute_du_dr(self, r, udu, E, l=0):
        """Return `(du, ddu)` for solve_ivp.

        Integrates the radial equation at the specified energy and angular momentum.

        Arguments
        ---------
        r : array-like
            Radius.  `r.shape == (n,)`
        udu: array-like
            Current solution and derivative, `u, du = udu`.  `udu.shape == (2, n)`.
        E : float
            Energy for state.
        l : int
            Angular momentum quantum number.
        """

        nu = l + self.dimension / 2 - 1
        u, du = udu
        V = self.V(r=r)
        ddu = (2 * self.m * (V - E) / self.hbar ** 2 + (nu ** 2 - 0.25) / r ** 2) * u
        return (du, ddu)

    def get_r_u_du_backwards(
        self, E, u0=None, R=None, R_max=None, tol=1e-8, max_step=None, l=0
    ):
        """Return `(rs, us, dus)` for the radial equation.

        Solves the radial equation by integrating from `R_max` backwards to `R` in order
        to kill any exponentially growing solutions.

        Arguments
        ---------
        E, l : float, int
            Energy, and angular momentum.
        u0 : float | None
            Initial size of the radial wavefunction `u0 = u(R_max)`.  If not provided,
            then we estimate by first performing a quick low-tolerance integration to
            `r_v`.
        R, R_max : float | None
            Integrate from `R_max` down to `R`.  If not provided, we estimate based on the
            length scales of the problem.
        tol : float
            Tolerance for final integration.  Assumes everything is roughly order unity,
            so this should be interpreted as a relative tolerance.
        max_step : float | None
            Passed to `solve_ivp`.  Estimated if not provided.

        Returns
        -------
        rs : array-like
            Radial points (abscissa) for solution.
        us : array-like
            Radial wavefunction at these points `u(r)`
        dus : array-like
            Derivative of the radial wavefunction at these points.
        """
        a = self.get_a(E=E)
        r_v = self.get_r_v(E=E)

        if R is None:
            R = np.sqrt(tol) * a

        r_v = max(r_v, R)

        if R_max is None:
            assert tol < 1
            R_max = r_v + a * max(1, -np.log(tol))

        if u0 is None:
            # Choose a reasonable initial condition
            if R_max <= r_v:
                u0 = 1
            else:
                u0 = np.sqrt(tol)
                y0 = (u0, -a * u0)
                res = sp.integrate.solve_ivp(
                    partial(self.compute_du_dr, E=E, l=l),
                    t_span=(R_max, r_v),
                    y0=y0,
                    max_step=(R_max - r_v) / 10,
                    atol=tol,
                    rtol=1e-3,
                )
                us_, du_ = res.y
                u0 = us_[0] / us_[-1]

        y0 = (u0, -u0 * a)
        R_span = (R_max, R)
        if max_step is None:
            max_step = abs(R_max - R) / 10
        res = sp.integrate.solve_ivp(
            partial(self.compute_du_dr, E=E, l=l),
            t_span=R_span,
            y0=y0,
            max_step=max_step,
            atol=tol,
            rtol=tol,
        )
        rs = res.t
        us, dus = res.y
        return rs, us, dus

    def objective(self, E, r0_a=None, R_max_a=None, N=10, deg=6, tol=1e-8, **kw):
        """Short-distance objective function for finding E.

        Arguments
        ---------
        tol : float
            Tolerance for final integration.  Assumes everything is roughly order unity,
            so this should be interpreted as a relative tolerance.
        r0_a : float | None
            Dimensionless quantity `r0/a` specifying how far down to integrate.
            If not provided, then estimate from tolerance and length scales.
        R_max_a : float | None
            Dimensionless quantity `R_max/a` specifying R_max.
        N, deg : int
           How many points and the degree to use when extrapolating from `r0` to 0.  We
           use `np.polyfit`.
        """
        a = self.get_a(E=E)
        R_max = None if R_max_a is None else a * R_max_a
        r0 = None if r0_a is None else a * r0_a
        rs_, us_, dus_ = self.get_r_u_du_backwards(E, R=r0, R_max=R_max, tol=tol, **kw)
        P = np.polyfit(rs_[-N:], us_[-N:], deg=deg)
        u0 = np.polyval(P, 0)
        return u0 / abs(us_).max()
