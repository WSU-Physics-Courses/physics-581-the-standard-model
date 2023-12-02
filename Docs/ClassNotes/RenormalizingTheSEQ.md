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

To do this, I suggest you complete the following tasks:

1. Get a method to numerically solve the radial Schrödinger equation.  You will use this
   to compute your "experimental" data.  Note: this is a little challenging due to the
   nature of the Coulomb potential being singular at the origin, and long-ranged
   (falling off slowly $V(r) \sim 1/r$ at long distances).  The following notes have
   many hints and suggestions, but please try to use them only as needed.

   ```{toctree}
   ---
   maxdepth: 2
   titlesonly:
   glob:
   ---
   RG_SEQ/*
   ```

2. Perform the perturbative analysis discussed in {cite:p}`Lepage:1997` but following
   with your own data.  If there is anything you don't know how to do, please ask on the
   [Hypothes.is document](https://hyp.is/ITZMcI5tEe6EBmPaQth57Q/emailwsu-my.sharepoint.com/personal/m_forbes_wsu_edu/Documents/Courses/Physics%20581%20Standard%20Model/1997_Lepage-HowToRenormalizeTheSchrodingerEquation.pdf?CT=1701230953632&OR=ItemsView) so I can help fill in background in class or online.


## How to Solve the Schrödinger Equation

To easily work through the discussion in {cite:p}`Lepage:1997`, we must be able to
formulate and solve the Schrödinger equation for spherically symmetric potentials.

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

The essential idea is low-energy properties, like the bound state energies $E_n$, should
not be highly sensitive to details about the nature of the potential $V(r)$ at short
distances $r \ll 1/\Lambda$ where $\Lambda$ is a large momentum scale.

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
