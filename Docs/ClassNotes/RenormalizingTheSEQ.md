---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.8
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

(sec:RG-SEQ)=
How to Renormalize The Schrödinger Equation
===========================================

Here we work through the example discussed in {cite:p}`Lepage:1997`, which uses the
example of bound states in a spherically symmetric potential, with Coulomb-like
behaviour for large radii $r\rightarrow \infty$, but with "unknown" corrections at short
distances $r\rightarrow 0$:

\begin{gather*}
  \left(\frac{-\hbar^2\nabla^2}{2m} + V(r) - E\right)\Psi(r, \Omega) = 0, \qquad
  \lim_{r\rightarrow \infty} V(r) \rightarrow \frac{-\alpha}{r}.
\end{gather*}


## How to Solve the Schrödinger Equation

To easily work through the discussion in {cite:p}`Lepage:1997`, we must be able to
formulate and solve the Schrödinger equation fr spherically symmetric potentials.

### Radial Equation
Spherical symmetry allows use to express this as a simple 1D boundary value problem
(BVP) for the radial equation:

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
Highly accurate solutions can be obtained by shooting, and we present details in
{ref}`sec:RadialSEQ` about how to do this.

## The Essence

The essential idea is low-energy properties, like
the bound state energies $E_n$, should not be highly sensitive to details about
the nature of the potential $V(r)$ at short distances $r \ll \lambda$ 

## Numerical Approach

The technical problem is to compute properties such as the bound state energies $E_n$ of
the potential $V(r)$. {cite:p}`Lepage:1997` presents a nice derivation about how to do
this using perturbation theory with various contact interaction terms, and I urge you to
follow and reproduce this discussion.

Here we will take an alternative numerical approach, defining our own effective
potential as follows:

\begin{gather*}
  V(r) = P(r/z)f_a(r) - \frac{\alpha}{r}\bigl(1-f_a(r)\bigr), \qquad
  f_a(r) = e^{-r^2/2a^2}, \qquad
  P(z) = \sum_{n} c_n \frac{z^n}{n!}.
\end{gather*}

This satisfies the criteria laid out in {cite:p}`Lepage:1997`:

1. We incorporate the correct long-range behavior through the cutoff function $1-f_a(r)$
   which goes to zero exponentially fast for $r>a$.
2. We have introduced an ultraviolet cutoff $a$ into our theory which softens the
   potential at short distances.
3. We have added "local" corrections via the parameters $c_n$.  The locality is provided
   by the cutoff factor $f_a(r)$.

```{code-cell} ipython3
import warnings;warnings.simplefilter("error")
from functools import partial
from scipy.optimize import root
from tqdm import tqdm
from phys_581 import seq

class SEQ(seq.CoulombSEQ):
    c = [-2.0]
    a = 1.0
    
    E_tol = 1e-4
    
    def f(self, r):
        return np.exp(-(r/self.a)**2/2)
        
    def V(self, r):
        f = self.f(r)
        return np.polyval(self.c, r/self.a)*f + (1-f)*super().V(r)
    
    def fit(self, Es, cs=None, **kw):
        """Fit the specified energies."""
        if cs is None:
            cs = self.c
        cs = list(cs)[:len(Es)]
        cs = cs + [0.0]*(len(Es) - len(cs))
        self.c = cs
        
        def objective(cs, Es):
            self.c[:len(Es)] = cs
            Es_ = np.array([self.compute_E(_E, tol=self.E_tol, lam=0.999) for _E in Es])
            return Es_/Es - 1
        
        for n in tqdm(range(1, len(Es)+1)):
            f = partial(objective, Es=Es[:n])
            res = root(f, self.c[:n], **kw)
        self.res = res
        if not res.success:
            raise Exception(res.message)
        self.c[:len(Es)] = res.x

        
s = SEQ()

Es = np.array([
    -1.28711542, 
    -0.183325753,
    -0.0703755485,
    -0.0371495726,
    -0.0229268241,
    -0.0155492598,
    -0.00534541931,
    -0.00129205010])
```

```{code-cell} ipython3
#s.compute_E(-1)
s.fit(Es[:2])
```

```{code-cell} ipython3
%time s.fit(Es[:2])
```

## Figure 1

Here we reproduce Fig. 1 from {cite:p}`Lepage:1997` using our potential:

```{code-cell} ipython3
:tags: [hide-cell]

s0 = seq.CoulombSEQ()
Es0 = np.array([s0.compute_E(_E, lam=0.999) for _E in tqdm(Es)])

fig, ax = plt.subplots()
ax.loglog(-Es, abs(Es0/Es - 1), ':o', label=r"$-\alpha/r$")
for a in [1.0, 0.5, 0.25, 0.125]:
    s1 = SEQ(a=a)
    s1.fit(Es[:1])
    Es1 = np.array([s1.compute_E(_E, lam=0.999) for _E in tqdm(Es)])
    ax.loglog(-Es, abs(Es1/Es - 1), ':s', label=f"$a={a:.2f}$, $c_0={s1.c[0]:.4f}$")
ax.legend(loc='upper left')
ax.set(ylim=(1e-4, 1), xlabel="$-E$", ylabel=r"$|\Delta E/E|$")
```

```{code-cell} ipython3
s2 = SEQ(a=1, E_tol=1e-4)
s2.fit(Es[:2])
Es2 = np.array([s2.compute_E(_E, lam=0.999) for _E in tqdm(Es)])
```

```{code-cell} ipython3
s2.compute_E(Es[0]) - Es[0]
```

```{code-cell} ipython3
plt.loglog(abs(Es), abs(Es2/Es - 1))
```

[Manim Community]: <https://www.manim.community/>
[Jacobi elliptic functions]: <https://en.wikipedia.org/wiki/Jacobi_elliptic_functions>
[glue]: <https://myst-nb.readthedocs.io/en/latest/use/glue.html>
[MyST]: <https://myst-parser.readthedocs.io/en/latest/> "MyST - Markedly Structured Text"
[Sphinx]: <https://www.sphinx-doc.org/>
[Markdown]: <https://daringfireball.net/projects/markdown/>
[MyST Cheatsheet]: <https://jupyterbook.org/reference/cheatsheet.html>
[Jupyter Book]: <https://jupyterbook.org>
[Jupyter Book with Sphinx]: <https://jupyterbook.org/sphinx/index.html>
[Jupyter]: <https://jupyter.org> "Jupyter"
[Jupytext]: <https://jupytext.readthedocs.io> "Jupyter Notebooks as Markdown Documents, Julia, Python or R Scripts"
[MySt-NB]: <https://myst-nb.readthedocs.io>
[Liouville's Theorem]: <https://en.wikipedia.org/wiki/Liouville%27s_theorem_(Hamiltonian)>
[Laplace-Beltrami operator]: <https://en.wikipedia.org/wiki/Laplace%E2%80%93Beltrami_operator>
[angular momentum operator]: <https://en.wikipedia.org/wiki/Angular_momentum_operator>
[hydrogenic atoms]: <https://en.wikipedia.org/wiki/Hydrogen-like_atom>
[Bessel function]: <https://en.wikipedia.org/wiki/Bessel_function>
[orthogonal polynomials]: <https://en.wikipedia.org/wiki/Orthogonal_polynomials>
