---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.2
kernelspec:
  display_name: Python 3 (phys-581)
  language: python
  name: phys-581
---

```{code-cell} ipython3
:tags: [hide-cell]

import mmf_setup;mmf_setup.nbinit()
import logging;logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
```

How to Solve the (Radial) Schrödinger Equation II
=================================================

:::{margin}
Here we will consider only the spherically symmetric ($S$-wave) orbits and phase-shifts
$l=m=0$ with constant spherical harmonics $Y_{0}^{0}(\theta, \phi)$. If you are not
familiar with these, you should consider them as "Fourier" components on the sphere.
This is most obvious for $l=0$ where $Y_{0}^{m}(\theta, \phi) \propto e^{\I m\phi}$ are
plane waves.  The $\theta$ dependence is more complicated due to the coordinate system.
:::
## Hydrogen Atom

For testing we use the exact energies for the non-relativistic [hydrogen atom][]:
\begin{gather*}
  V(r) = \frac{-\alpha}{r}, \qquad
  E_{l,n} = \frac{-m\alpha^2/2\hbar^2}{2(1+l+n)^2}.
\end{gather*}
The eigenstates can also be expressed analytically:
\begin{gather*}
  \psi_{n,l,m}(r, \theta, \phi) \propto e^{r/na}\left(\frac{2r}{na}\right)^{l}
  L_{n-l-1}^{2l+1}\left(\frac{2r}{na}\right)Y_{l}^{m}(\theta, \phi), \qquad
  a = \frac{\hbar^2}{m \alpha},
\end{gather*}
where $m$ is the reduced mass $m = m_em_p/(m_e+m_p)$, $a$ is the [reduced Bohr
radius][], $L_{n-l-1}^{2l+1}(\rho)$ is a [generalized Laguerre polynomial][] of degree 
$n-l-1$, and $Y_{l}^{m}(\theta, \phi)$ is a [spherical harmonic][] function of degree 
$l$ and order $m$.

These exact solutions might be used to deal with the singular properties of the Coulomb
potential, but we postpone discussing this for now in favour of more general
techniques.


## To Do

* Fix asymptotic behaviour at $r=0$ for $l>0$.
* Improve estimate of $u_0$ in `get_r_u_du_backwards`: this is causing major performance
  issues for large $n$.  If it is set to a reasonable value (so that the maximum
  solution is order 1) then with a good `max_step` we get nice convergence.

## Simple Solution: Shooting

The simple solution is to use {meth}`scipy.integrate.solve_ivp` to shoot a solution to
some large radius $R$:

\begin{gather*}
  u(0) = 0, \qquad u'(0) = 1, \qquad u(R)\approx 0.
\end{gather*}

A couple of issues need to be dealt with:

:::{margin}
The long-range behaviour of the Coulomb potential leads to it being classified as a non-local potential.  Specifically, this leads to a infinite total cross-section for [Rutherford scattering], and non-extensivity.
*(If you have not done so, you should calculate the energy of a uniformly charged sphere of charge $Q$ and radius $R$.  If the system were extensive, then holding the charge densty $\rho = 3Q/4\pi R^3$ fixed, the energy would scale as the volume $E = f(\rho) R^{3}$.  Instead, it grows much more quickly.  Basically: with local interactions, particles only feel the effects of nearby particles.  The diverging cross-section means that charged particles feel the effects of all other particles.)*
:::
1. The integrand is singular at $r=0$, so we must start from a small non-zero value
   $r_0$, thus we replace our boundary conditions with:
   

   \begin{gather*}
     u(r_0) = r_0, \qquad u'(r_0) = 1, \qquad u(R)\approx 0.
   \end{gather*}
  
2. The potential extends to $r\rightarrow \infty$ without falling off very quickly $V(r) \propto r^{-1}$.  Thus, determining what radius $R$ is "large" is a bit subtle.

:::{note}
The results of this discussion are codified in the `phys_581.seq.SEQ` class which we use elsewhere.  These notes document some of the explorations that led to that code.
:::

```{code-cell} ipython3
from scipy.integrate import solve_ivp
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import brentq
from functools import partial

from phys_581 import dvr

d = 3
hbar = m = alpha = 1.0
r0 = 0.2

rng = np.random.default_rng(seed=2)

# 3rd order random polynomial
P = (rng.random(3) - 0.5)*2

def V(r):
    return np.polyval(P, (r/r0)**2) * np.exp(-(r/r0)**2/2) - alpha/np.sqrt(r**2 + r0**2)

def V(r):
    """Hydrogen atom for testing."""
    return -alpha/r

def compute_du_dr(r, u_, E, l=0):
    nu = l + d/2 - 1
    u, du = u_
    ddu = (2*m*(V(r) - E) / hbar**2 + (nu**2-0.25)/r**2)*u
    return (du, ddu)

def get_r_u_du(E, r0=1e-12, R=30.0, R_max=None, tol=1e-8, max_step=0.1, l=0):
    y0 = (r0, 1)
    R_span = (r0, R)
    res = solve_ivp(partial(compute_du_dr, E=E, l=l), t_span=R_span, y0=y0,
                    max_step=max_step, atol=tol, rtol=tol)
    rs = res.t
    us, dus = res.y
    return rs, us, dus
```

```{code-cell} ipython3
Es = -m*alpha**2/2/hbar**2/ np.arange(1, 20)**2

fig, ax = plt.subplots()
lss = ["-", "--", "-.", ":"]
for _nE, E in enumerate(Es[:3]):
    for _nt, log10_tol in enumerate([-8, -10, -12, -13]):
        r, u, du = get_r_u_du(E=E, tol=10**log10_tol, R=40)
        label = f"tol = $10^{{{log10_tol}}}$" if _nE == 0 else None
        ax.plot(r, u, ls=lss[-1-_nt], c=f"C{_nE}", label=label)
ax.set(ylim=(-1, 1), ylabel='$u(r)$', xlabel="$r$")
ax.legend();
```

We see that this works, but because of the exponentially growing solution for larger
$r$, one needs quite hight tolerances on the integrator.  We see another problem, that
as the binding energy gets smaller, the radius of the wavefunction gets larger.  This will
introduce errors in the energy.

These are two different types of errors: the first is a result of the integrator not
taking a small enough step size.  This is sometimes called an "ultraviolet" (UV) error
because it results from short-wavelength physics.  The second error is sometimes called
an "infrared" (IR) error because it results from long-wavelength physics (our box is not
large enough).  To get an accurate result, one must understand and mitigate both types
of error, each of which increases the computational cost of the calculation.

Let's complete the simple solution by shooting to find the eigenvalues.  We will polish
the solution using {py:func}`scipy.optimize.brentq`.  This requires a bracket, so we
start with a guess for $E$ then expand our search for a bracket.

```{code-cell} ipython3
import scipy.optimize

def f0(E, **kw):
    """Simple objective function u(R)=0"""
    r, u, du = get_r_u_du(E, **kw)
    return u[-1]

def get_E(E, f=f0, tol=1e-8, lam=0.9, **kw):
    """Shoot for the best E."""
    f_ = partial(f, tol=tol, **kw)
    f0 = f1 = 1
    while f0*f1 > 0:
        lam *= lam
        E0 = lam*E
        E1 = E/lam
        f0, f1 = map(f_, [E0, E1])
    return scipy.optimize.brentq(f_, E0, E1, xtol=tol)

get_E(Es[1], f=f0, R=20)
```

We can now explore the errors:

```{code-cell} ipython3
from tqdm import tqdm
from IPython.display import display, clear_output

tols = 10**(np.linspace(-4, -12, 10))
Rs = np.linspace(10.0, 60.0, 10)
R = 30.0
tol = 1e-9

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
ax = axs[0]
for _nE, E in enumerate(Es[:3]):
    Es_ = np.array([get_E(E, R=R, tol=_tol) for _tol in tqdm(tols)])
    ax.loglog(tols, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")
ax.axvline([tol], c='y', ls=':')
ax.set(xlabel="tol", ylim=(1e-12, 1e-2), ylabel="rel. error", title="UV Convergence")
ax.legend();

ax = axs[1]
for _nE, E in enumerate(Es[:3]):
    Es_ = np.array([get_E(E, R=_R, tol=tol) for _R in tqdm(Rs)])
    ax.semilogy(Rs, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")
ax.axvline([R], c='y', ls=':')
ax.set(xlabel="R", ylim=(1e-12, 1e-2), ylabel="rel. error", title="IR Convergence")
ax.legend();
plt.suptitle("Convergence with the $u(R_\max) = 0$ boundary condition")
clear_output()
display(fig)
plt.close('all')
```

On the left, we show the UV convergence for a fixed $R_\max$, and on the right, we show
$IR$ convergence for a fixed tolerance.  We see that the $n=0$ state is approaching
convergence, but that the $n=1$ and $n=2$ states have difficulty.  The presense of
convergence plateaus in one plot indicates that we have reached convergence in the
other.  On the right plot, we see that choosing an integration tolerance of $10^{-9}$
allows us to reach this level of precision for all states, but that the $n=1$ state
requires at least $R_\max > 40$ while the $n=2$ state requires $R_\max > 60$.

On the left, we see that a box size of $R_\max = 30$ allows us to achive high
convergence for the $n=0$ state, but that the $n=1$ and $n=2$ states are not converged.

## Asymptotic Behaviour

The IR convergence is best understood by considering the asymptotic behaviour.  For
large $r$ we have:

\begin{gather*}
  \lim_{r\rightarrow \infty} \frac{-\hbar^2}{2m} u''(r) \approx E u(r) + \frac{\alpha}{r}u(r).
\end{gather*}

This introduces two length scales $a$ and $r_v$:

\begin{gather*}
  a = \frac{\hbar}{\sqrt{-2m E}}, \qquad
  r_v = \frac{\alpha}{-E},\qquad
  \frac{u''(r)}{u(r)} \approx \frac{1}{a^2}\left(1  - \frac{r_v}{r}\right).
\end{gather*}

:::{margin}
For the ground state of the hydrogen atom, $a$ is the [Bohr radius].
:::
The meaning of these is that wavefunction exponential decays with length scale $a$

\begin{gather*}
  \lim_{r\rightarrow \infty} u(r) \propto e^{-r/a}.
\end{gather*}

with a turning point at $r_v$: $u''(r_v)\approx 0$.

:::{admonition} Notes
:class: dropdown
Quantitatively, we can include the centrifugal term when solving for $r_v$:

\begin{gather*}
  a^2r^2\frac{u''(r)}{u(r)} \approx \nu^2 - 1/4 - \frac{2mE}{\hbar^2}r_v - \frac{2m\alpha}{\hbar^2} r
  \approx (\nu^2 - 1/4)a^2 + r^2 + \frac{\alpha}{E} r,\\
  r_v \approx \frac{\alpha}{-E}\frac{1 + \sqrt{1 + (1 - 4\nu^2)\frac{a^2E^2}{\alpha^2}}}{2}.
\end{gather*}

This will give a better approximation of the turning point for large angular momenta $l$.
We also get a qualitative limit on $E$ for large $l$:

\begin{gather*}
  \frac{-\alpha}{a\sqrt{4\nu^2 - 1}} < E < 0
\end{gather*}

:::

This can be implemented as the boundary condition:

\begin{gather*}
  au'(R) = -u(R).
\end{gather*}

```{code-cell} ipython3
def f1(E, **kw):
    """Better objective function encoding long-distance behaviour."""
    a = hbar / np.sqrt(-2*m*E)
    r, u, du = get_r_u_du(E, **kw)
    return a*du[-1] + u[-1]

get_E(Es[1], f=f0, R=20), get_E(Es[1], f=f1, R=20)
```

```{code-cell} ipython3
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
ax = axs[0]
for _nE, E in enumerate(Es[:3]):
    Es_ = np.array([get_E(E, f=f1, R=R, tol=_tol) for _tol in tqdm(tols)])
    ax.loglog(tols, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")
ax.axvline([tol], c='y', ls=':')
ax.set(xlabel="tol", ylim=(1e-12, 1e-2), ylabel="rel. error", title="UV Convergence")
ax.legend();

ax = axs[1]
for _nE, E in enumerate(Es[:3]):
    Es_ = np.array([get_E(E, f=f1, R=_R, tol=tol) for _R in tqdm(Rs)])
    ax.semilogy(Rs, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")
ax.axvline([R], c='y', ls=':')
ax.set(xlabel="R", ylim=(1e-12, 1e-2), ylabel="rel. error", title="IR Convergence")
ax.legend();
plt.suptitle("Convergence with the $u(R_\max) + au'(R)= 0$ boundary condition")
clear_output()
display(fig)
plt.close('all')
```

This gives slightly better behaviour, but it is not dramatic.

## A Solution: Integrate Backwards

The real problem is that, at large $r$, we have an exponentially diverging solution in
addition to the desired exponentially decaying solution.  A good strategy is to
integrate backwards from a large $R_\max$ to $R$: this will quickly dampen any of the
exponentially diverging solution, allowing us to match the boundary condition with high
accuracy.

What initial conditions should we use?  A good strategy is to use the asymptotic form, now including the Coulomb piece to guess.  We do a preliminary integration from $R_\max$ to the turning $r_v$ to get the scales, then repeat with high accuracy and an initial value so that $u(r_v) \approx 1$.

How large must $R_\max$ be?  If we want an accuracy of $\epsilon$, then a good choice is

\begin{gather*}
   \epsilon \approx e^{-{R_\max - r_v}r/a}, \qquad
   R_\max > r_v - a \ln \epsilon.
\end{gather*}

```{code-cell} ipython3
def get_r_u_du_backwards(E, u0=None, R=30.0, R_max=None, tol=1e-8, max_step=None, l=0):
    a = np.sqrt(-2*m*E) / hbar
    r_v = max(alpha/(-E), R)   # Could be improved
    if R_max is None:
        R_max = r_v - a * min(-1, np.log(tol))

    if u0 is None:
        # Choose a reasonable initial condition
        if R_max <= r_v:
            u0 = 1
        else:
            u0 = np.sqrt(tol)
            y0 = (u0, -a*u0)        
            res = solve_ivp(partial(compute_du_dr, E=E, l=l), 
                            t_span=(R_max, r_v), 
                            y0=y0, 
                            max_step=(R_max-r_v)/10,
                            atol=tol, rtol=1e-3)
            us_, du_ = res.y
            u0 = us_[0] / us_[-1]
    y0 = (u0, -a * u0)
    R_span = (R_max, R)
    if max_step is None:
        max_step = abs(R_max - R)/10
    res = solve_ivp(partial(compute_du_dr, E=E, l=l), t_span=R_span, y0=y0,
                    max_step=max_step, atol=tol, rtol=tol)
    rs = res.t
    us, dus = res.y
    return rs, us, dus
```

In principle, one can now integrate both forward an backwards, meeting in the middle.

```{code-cell} ipython3
E = Es[2]
R = 10.0
rs, us, dus = get_r_u_du(E, R=R)
rs_, us_, dus_ = get_r_u_du_backwards(E, R=R)
plt.plot(rs, us)
plt.plot(rs_, us_*us[-1]/us_[-1])
```

We now show that it is sufficient to run the backwards iteration.  As running to $r=0$
is slightly problematic, we choose a small $r_0$ and then do a polynomial extrapolation
of the last few points to extrapolate to $r=0$, requiring that $\lim_{r\rightarrow 0} u(r) = 0$.

Notes:
* Using the simple condition $r_0 u'(r_0) - u(r_0) = 0$ corresponds to a linear
  extrapolation of the last two points, but does not give high precision, even if $r_0$
  is carefully chosen.  (This is the usual issue of balancing truncation and roundoff
  errors).
* The parameters in the default version have been tweaked a bit by hand exploring up to
  the $n=7$ hydrodgen solution.  This could be done a bit better.

```{code-cell} ipython3
def f_(E, ar0=1e-2, aR_max=None, N=10, order=6, **kw):
    a = hbar/np.sqrt(-2*m*E)
    if aR_max is None:
        R_max = None
    else:
        R_max = aR_max * a
    r0 = ar0 * a
    rs_, us_, dus_ = get_r_u_du_backwards(E, R=r0, R_max=R_max, **kw)
    return np.polyval(np.polyfit(rs_[-N:], us_[-N:], deg=order), 0)/abs(us_).max()
    return (r0*dus_[-1] - us_[-1])/np.sqrt((r0*dus_[-1])**2 + (us_[-1])**2)
```

```python
print(np.abs(list(map(partial(f_, aR_max=20), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=40), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=50), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=60), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=100), Es))).max())
```

```{code-cell} ipython3
print(np.abs(list(map(partial(f_, aR_max=20), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=40), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=50), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=60), Es))).max())
print(np.abs(list(map(partial(f_, aR_max=100), Es))).max())
```

```{code-cell} ipython3
ns_ = np.linspace(0.8, 20, 200)
Es_ = -m*alpha**2/2/hbar**2/ ns_**2
fs_ = np.array(list(map(partial(f_), tqdm(Es_))))
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(ns_, fs_)
ax.grid(True)
ax.set(xlabel="$n$", ylabel="boundary condition", xticks=np.arange(20));
```

From the previous plot, we see that the normalized boundary condition is very regular,
lending itself to nice root-finding.

```{code-cell} ipython3
Es = -m*alpha**2/2/hbar**2/ np.arange(1, 10)**2
#get_E(Es[-1], f=f_, tol=1e-9)/Es[-1] - 1

# Broken if tol too large... check
get_E(Es[-1], f=f_, tol=1e-8)/Es[-1] - 1
```

```python
E = Es[-1]
a = np.sqrt(-2*m*E)/hbar
r, u, du = get_r_u_du_backwards(E=E, R=0.1/a, R_max=60/a, u0=1e-12, tol=1e-10, max_step=)
plt.plot(r, u, '+')
axis = plt.axis()
plt.plot(r, 1e10*np.exp(-a*r))
plt.axis(axis)
```

```python
from scipy.integrate import solve_bvp
def compute_du_dr_E(r, u_, E_):
    E, = E_
    du, ddu = compute_du_dr(r, u_, E=E)
    return du, ddu

def bc(u0_, uR_, E_):
    E, = E_
    a = np.sqrt(-2*m*E)/hbar
    u0, du0 = u0_
    uR, duR = uR_
    return (u0, du0 - 1, a*duR + uR)

E = Es[0]
a = np.sqrt(-2*m*E)/hbar
r0 = 1e-3/a
R_max = 60.0/a
r = np.linspace(r0, R_max)
u = r*np.exp(-a*r)
du = (1-a*r)*np.exp(-a*r)
sol = solve_bvp(compute_du_dr_E, bc, x=r, y=(u, du), p=[E], tol=1e-12, bc_tol=1e-12)
r = sol.x
u, du = sol.y
E_ = sol.p[0]
plt.plot(r, u), E/E_
#bc((u[0], du[0]), (u[-1], du[-1]), [E])
```

```python
from tqdm import tqdm

Es = -m*alpha**2/2/hbar**2/ np.arange(1, 5)**2

tols = 10**(np.linspace(-4, -12, 10))
aRs = np.linspace(10.0, 60.0, 10)
aR = 30.0
tol = 1e-6

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
ax = axs[0]
for _nE in reversed([0, 1, 2, len(Es)-1]):
    E = Es[_nE]
    Es_ = np.array([get_E(E, f=f_, aR_max=aR, tol=_tol) for _tol in tqdm(reversed(tols))])
    axs[0].loglog(tols, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")

    Es_ = np.array([get_E(E, f=f_, aR_max=_aR, tol=tol) for _aR in tqdm(reversed(aRs))])
    axs[1].semilogy(aRs, abs(Es_/E - 1), '+-', label=f"$n={_nE}$", c=f"C{_nE}")

axs[0].axvline([tol], c='y', ls=':')
axs[0].set(xlabel="tol", ylim=(1e-12, 1e-2), ylabel="rel. error", title="UV Convergence")
axs[0].legend();

axs[1].axvline([R], c='y', ls=':')
axs[1].set(xlabel="R", ylim=(1e-12, 1e-2), ylabel="rel. error", title="IR Convergence")
axs[1].legend();

plt.suptitle("Convergence with the $u(R_\max) + au'(R)= 0$ boundary condition")
clear_output()
display(fig)
plt.close('all')
```

## How to Solve the Schrödinger Equation (Old)

We will consider here only spherically symmetric problems for simplicity.  We start with
some numerical code to solve the radial Schrödinger equation.  We do this in two ways:
one with finite differences (easy to code, but not very accurate), and once with a
Bessel-function discrete variable representation (DVR) basis as described in
{cite:p}`LC:2002`.  The latter is a very accurate spectral method, (exponentially
accurate for analytic potentials) but might seem like a bit of a black box until you
study DVR bases a bit more.

### Radial Equation

:::{margin}
With the usual spherical coordinates $\phi$ and $\theta$, and
angular momentum quantum numbers $l$ and $m$, we have 

| $d$ | $\Omega$         | $\lambda$ |
|-----|------------------|-----------|
| $2$ | $\phi$           | $\abs{m}$ |
| $3$ | $(\phi, \theta)$ | $l$       |

for $d=2$ dimensions (i.e. cylindrical coordinates, but no $z$ dependence) and $d=3$
dimensions respectively.  The total orbital [angular momentum operator] is related to the
[Laplace-Beltrami operator] by

\begin{gather*}
  \op{L}^2 = \hbar^2\op{\Delta}_{S^{d-1}}
\end{gather*}

which has eigenvalues $\hbar^2 m^2$ in $d=2$ (i.e. $\op{L}^2 \equiv \op{L}_z^2$ and
$\hbar^2 l (l - 1)$ in $d=3$.
:::
Our aim is to satisfy the Schrödinger equation for central potentials in $d$-dimensions,
which we can do in the usual way by expressing the wavefunction $\Psi(r, \Omega) =
\psi(r)Y_{\lambda}(\Omega)$ in terms of the radial wavefunction $u(r)$ and
the appropriate generalized spherical harmonics $Y_{\lambda}(\Omega)$:

\begin{gather*}
  \left(\frac{-\hbar^2\nabla^2}{2m} + V(r) - E\right)\Psi(r, \Omega) = 0, \\
  \left[
    \frac{\hbar^2}{2m}
    \underbrace{
      \left(-\diff[2]{}{r} + \frac{\nu^2 - 1/4}{r^2}\right)
    }_{\op{K}} +
  V(r)\right]u(r) = E u(r),\\
  u(r) = r^{(d-1)/2}\psi(r), \qquad
  \nu = \lambda + \frac{d}{2} - 1.  
\end{gather*}

:::::{admonition} Details
:class: dropdown

Here $\Omega$ is the generalized solid angle.  Rotational invariance implies that
angular momentum is a good quantum number, so all eigenstates can be factored $\Psi(r, \Omega) =
\psi(r)Y_{\lambda}(\Omega)$ where $Y_{\lambda}(\Omega)$ is the generalized spherical
harmonic on the ($d-1$)-dimensional sphere and $\lambda \in \{0, 1, \dots\}$ is the
generalized angular momentum:

\begin{gather*}
  \nabla^2 = \frac{1}{r^{d-1}}\diff{}{r}\left(r^{d-1} \diff{}{r}\right)
  + \frac{1}{r^2}\Delta_{S^{d-1}}, \\
  \Delta_{S^{d-1}}Y_{\lambda}(\Omega) =
  \lambda (\lambda + d - 2)Y_{\lambda}(\Omega),
\end{gather*}

where $\Delta_{S^{d-1}}$ is the [Laplace-Beltrami operator].  Introducing the radial
wavefunction $u(r)$, this becomes quadratic:

\begin{gather*}
  u(r) = r^{(d-1)/2}\psi(r), \qquad
  \nu = \lambda + \frac{d}{2} - 1,\\
  r^{(d-1)/2}\nabla^2 \psi(r) = 
  \left(\diff[2]{}{r} - \frac{\nu^2 - 1/4}{r^2}\right)u(r).
\end{gather*}

This follows after a little algebra from 

\begin{gather*}
  \frac{1}{r^{d-1}}\diff{}{r}\left(r^{d-1} \diff{}{r}\right)r^{(1-d)/2}u(r)\\
  =
  r^{(1-d)/2}
  \left(
    \diff[2]{}{r}
    -
    \frac{(d-3)(d-1)}{4}
  \right)u(r).
\end{gather*}
:::::

This can be solved quite simply -- but not very accurately -- with finite differences.
Highly accurate solutions can be obtained by shooting, but this can be inefficient.  I
highly recommend that you stop and implement your own solution to this problem for
various functions $V(r)$ and parameters $\nu$.  As a check, you should be able to find
the eigenstates for [hydrogenic atoms] with $V(r) \propto 1/r$.

```{code-cell} ipython3
from phys_581 import seq
from importlib import reload
reload(seq)

class SEQ(seq.CoulombSEQ):
    """Classic Coulomb potential."""
    hbar = m = 1
    alpha = 1.0

    def V(self, r):
        return -self.alpha / r
    
    def get_E0(self, ns, l=0):
        """"Return the exact energies for hydrogen."""
        E_Ry = -self.m * self.alpha**2 / 2 / self.hbar**2
        return E_Ry / (ns + l + 1)**2

    def plot(self, N=None, Es=None,l=0, fig=None):
        """Plot N wavefunctions"""
        if Es is None:
            Es = self.get_E0(ns=np.arange(N), l=l)
    
        if fig is None:
            fig, ax = plt.subplots(figsize=(8,5))
        else:
            ax = fig.gca()
        
        for n, E in enumerate(Es):
            r, u, du = self.get_r_u_du_backwards(E=E)
            ax.plot(r/(1+n)/self.get_a(E=E), u, label=f"n={n}, l={l}")
        ax.set(xlabel="r/a(1+n)", ylabel="u(r)", xlim=(-0.1, 3))
        ax.legend()
        return fig

s = SEQ()
l = 0
s.plot(5, l=l);
```

```{code-cell} ipython3
class SEQ1(SEQ):
    """Truncated Coulomb potential."""
    r0 = 0.5
    V0 = None
    
    def V(self, r):
        V0 = self.V0
        if V0 is None:
            V0 = super().V(self.r0)
        return np.where(r<self.r0, V0, super().V(r))
    
    def plot(self, N, l=0, **kw):
        E0s = super().get_E0(ns=np.arange(N), l=l)
        Es = [self.compute_E(_E, l=l) for _E in E0s]
        super().plot(Es=Es, **kw)
    
s1 = SEQ1()
s1.plot(5, l=l);
```

### DVR Basis for the Radial Equation

:::{margin}
See {cite:p}`LC:2002` for details.  We follow most of the notations, except we use
$u(r) \equiv \phi(r)$ for the radial wavefunction, and $k_\max \equiv K$ for the maximum
momentum.
:::
Here we will present without proof a spectral method based on the Bessel-function discrete
variable representation (DVR).  Here we introduce a basis $\ket{F_{\nu,n}}$ obtained by
projecting onto a space with wave-vectors less than $\abs{k} < k_\max$:

\begin{gather*}
  \op{P} = \int_{\abs{\vect{k}}<k_\max}\!\!\!\!\!\d^{d} \vect{k}\; \ket{\vect{k}}\bra{\vect{k}}.
\end{gather*}

In this basis, the following representation for the operator $\op{K}$ is exact:

:::{margin}
In $d=3$ dimensions, we must consider fractional $\nu$, so the [numerical routine `jn_zeros` in 
SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jn_zeros.html)
does not suffice.  We provide our own implementation in the accompanying code.  One
could use the spherical [Bessel function]s of the first kind instead:

\begin{gather*}
  j_{n}(z) = \sqrt{\frac{\pi}{2z}}J_{n + 1/2}(z).
\end{gather*}
:::
\begin{gather*}
  \mat{K}^{(\lambda)}_{m,n} 
  = \braket{F_{\nu, m}|\left(-\diff[2]{}{r} + \frac{\nu^2 - 1/4}{r^2}\right)|F_{\nu, n}}\\
  = k_{\max}
  \begin{cases}
    \frac{1}{3}\left(1 + \frac{2(\nu^2 - 1)}{z_{\nu,n}^2}\right) & m=n,\\
    (-1)^{n-m}\frac{8 z_{\nu,m}z_{\nu,n}}{(z_{\nu, m}^2 - z_{\nu, n}^2)^2}, & m \neq n,
  \end{cases}\\
  J_{\nu}(z_{\nu, n}) = 0,
\end{gather*}

where $z_{\nu, n}$ are the roots of the [Bessel function]s of the first kind.

Furthermore, the basis is quasi-local, so that the potential operator can be expressed
as a diagonal matrix

\begin{gather*}
  \braket{F_{\nu, m}|V(\op{r})|F_{\nu, n}} \approx
  V(r_{\nu, n}) \delta_{mn}, \qquad
  r_{\nu, n} = \frac{z_{\nu, n}}{k_\max}.
\end{gather*}

This is not exact, but provides exponential accuracy for analytic potentials.

:::::{admonition} Basis functions
:class: dropdown

We will not need them for our work here, but for reference, the basis functions are:

\begin{gather*}
  \braket{F_{\nu, n}|F_{\nu, m}} = \int_0^{\infty}\d{r}\; F_{\nu, n}(r)F_{\nu, m}(r) =
  \delta_{mn},\\
  F_{\nu, n}(r) =
  (-1)^{n+1} \frac{r_{\nu, n} \sqrt{2r}}
                     {r^2 - r_{\nu, n}^2}J_\nu(k_\max r).
\end{gather*}

:::{note}
One must be careful evaluating these near the roots where the denominator vanishes.  In
our code, we do this with a careful application of Taylor series and L'Hopital's rule.
:::

Note that these are zero at all abscissa $r_{\nu, m}$ except for $m = n$:

\begin{gather*}
  F_{\nu, n}(r_{\nu, m}) = \delta_{n,m}
  \underbrace{
    \sqrt{k_\max}
    (-1)^{n+1} \sqrt{\frac{z_{\nu, n}}{2}}
    J_\nu'(z_{\nu, n})
  }_{w_n}.
\end{gather*}

This is a key property of DVR bases which are closely related to the classical
[orthogonal polynomials] where a careful choice of both abscissa and weights leads to
twice the expected accuracy when integrating.  It also allows us to express a function
in the basis by simply computing its value at the abscissa:

\begin{gather*}
  f(r) = \sum_{n} f_n F_{\nu, n}(r), \qquad
  f(r_{\nu, m}) = \sum_{n} f_n F_{\nu, n}(r_{\nu, m}) = f_m w_{m},\\
  f_n = \braket{F_{\nu, n}|f} = \frac{f(r_{\nu, n})}{w_n}.
\end{gather*}

Finally, we note that from these properties, the
coefficients $w_n$ act as integration weights for functions expressed in the basis:

\begin{gather*}
  \int_0^{r} f(r) \d{r} =   
  \sum_{n} f_{n} \int_0^{r} F_{\nu, n}(r) \d{r} =   
\end{gather*}
:::::

### Demonstration

As a quick demonstration, we find the eigenstates of the spherical harmonic oscillator,
and the eigenstates of [hydrogenic atoms]:

\begin{gather*}
  V_{HO}(r) = \frac{m\omega^2r^2}{2}, \qquad 
  E_{l,n} = \hbar\omega\left(2n + l + \frac{d}{2}\right), \qquad
  a_{ho} = \frac{\hbar}{\sqrt{m\omega}},\\
  V_{H}(r) = \frac{-\alpha}{r}, \qquad 
  E_{l,n} = \underbrace{
    -\frac{m \alpha^2}{2\hbar^2}}_{-13.6\;\mathrm{eV}}
  \frac{1}{(l+n+1)^2},\qquad
  a_{h} = \frac{\hbar^2 }{m \alpha}.
\end{gather*}

The numerical value is given for hydrogen where $\alpha = e^2/4\pi\epsilon_0 \approx
14.4$eV Å.

```{code-cell} ipython3
from phys_581 import bessel, dvr

N = 10
d = 3
hbar = m = w = 1.0
a_ho = hbar / np.sqrt(m*w)
R = N*a_ho
k_max = N/a_ho

def V(r):
    return m * (w*r)**2 / 2
    
def get_E(l, N=N):
    n = np.arange(N)
    return hbar * w * (2*n + l + d/2)

for d in [2, 3, 4]:
  basis = dvr.SphericalDVRBasis(R=R, d=d, k_max=k_max)
  for l in [0, 1, 2, 3, 4]:
    r = basis.get_rn(l=l)
    H = hbar**2 / 2/ m * basis.get_K(l) + np.diag(V(r))
    assert np.allclose(np.linalg.eigvalsh(H)[:N], get_E(l=l))
```

```{code-cell} ipython3
from phys_581 import bessel, dvr

N = 5
d = 3
hbar = m = e = alpha = 1.0

a_h = hbar**2 / m / alpha
R = 10*N*a_h
k_max = 10*N/a_h

def V(r):
    return -alpha / r
    
def get_E(l, N=N):
    n = np.arange(N)
    return -m * alpha**2 / 2 / hbar**2 / (1 + n + l)**2

basis = dvr.SphericalDVRBasis(R=R, d=d, k_max=k_max)
for l in [0, 1, 2, 3, 4]:
    r = basis.get_rn(l=l)
    H = hbar**2 / 2/ m * basis.get_K(l) + np.diag(V(r))
    print(d, l)
    print(np.linalg.eigvalsh(H)[:N])
    print(get_E(l=l))
    #assert np.allclose(np.linalg.eigvalsh(H)[:N], get_E(l=l))
```

[Bessel function]: <https://en.wikipedia.org/wiki/Bessel_function>
[Bohr radius]: <https://en.wikipedia.org/wiki/Bohr_radius>
[Jacobi elliptic functions]: <https://en.wikipedia.org/wiki/Jacobi_elliptic_functions>
[Jupyter Book with Sphinx]: <https://jupyterbook.org/sphinx/index.html>
[Jupyter Book]: <https://jupyterbook.org>
[Jupyter]: <https://jupyter.org> "Jupyter"
[Jupytext]: <https://jupytext.readthedocs.io> "Jupyter Notebooks as Markdown Documents, Julia, Python or R Scripts"
[Laplace-Beltrami operator]: <https://en.wikipedia.org/wiki/Laplace%E2%80%93Beltrami_operator>
[Liouville's Theorem]: <https://en.wikipedia.org/wiki/Liouville%27s_theorem_(Hamiltonian)>
[Manim Community]: <https://www.manim.community/>
[Markdown]: <https://daringfireball.net/projects/markdown/>
[MyST Cheatsheet]: <https://jupyterbook.org/reference/cheatsheet.html>
[MyST]: <https://myst-parser.readthedocs.io/en/latest/> "MyST - Markedly Structured Text"
[MySt-NB]: <https://myst-nb.readthedocs.io>
[Rutherford scattering]: <https://en.wikipedia.org/wiki/Rutherford_scattering>
[Sphinx]: <https://www.sphinx-doc.org/>
[angular momentum operator]: <https://en.wikipedia.org/wiki/Angular_momentum_operator>
[glue]: <https://myst-nb.readthedocs.io/en/latest/use/glue.html>
[hydrogenic atoms]: <https://en.wikipedia.org/wiki/Hydrogen-like_atom>
[orthogonal polynomials]: <https://en.wikipedia.org/wiki/Orthogonal_polynomials>

```{code-cell} ipython3
k = 2*np.pi * np.fft.fftfreq(N, dx)
K = (hbar*k)**2/2/m
A = np.fft.ifft(K*np.fft.fft(np.eye(N), axis=1), axis=1)
V = np.diag(Vx)
C = None

E, psi = sp.linalg.eig(A+V, C)
inds = np.argsort(abs(E))
E, psi = E[inds], psi[:, inds]
E[:10].real
```

```{code-cell} ipython3

```


## See Also

* [Physics 555: Bound States in the 1D Schrödinger Equation](
  https://physics-555-quantum-technologies.readthedocs.io/en/latest/Notes/Shooting.html)
  

[hydrogen atom]: <https://en.wikipedia.org/wiki/Hydrogen_atom>
[reduced Bohr radius]: <https://en.wikipedia.org/wiki/Bohr_radius#Reduced_Bohr_radius>
[generalized Laguerre polynomial]: <https://en.wikipedia.org/wiki/Laguerre_polynomial#Generalized_Laguerre_polynomials>
[spherical harmonic]: <https://en.wikipedia.org/wiki/Spherical_harmonics>
[finite difference methods]: <https://en.wikipedia.org/wiki/Finite_difference_method>
[Dirichlet boundary conditions]: <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>
