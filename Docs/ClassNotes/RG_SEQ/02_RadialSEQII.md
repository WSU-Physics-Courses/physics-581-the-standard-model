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

Here we continue our previous explorations of numerical methods for solving the radial
Schrödinger Equation, trying to overcome some of the shortcomings seen there.  We will
focus on problems with a Coulomb potential that has both the mild singularity at $r=0$
and the long tail $V \sim 1/r$.

For testing we use the exact energies for the non-relativistic [hydrogen atom][]:
\begin{gather*}
  V(r) = \frac{-\alpha}{r}, \qquad
  E_{l,n} = \frac{-m\alpha^2/2\hbar^2}{2(1+l+n)^2}.
\end{gather*}
The eigenstates can also be expressed analytically:
\begin{gather*}
  \psi_{n,l,m}(r, \theta, \phi) \propto e^{-r/na}\left(\frac{2r}{na}\right)^{l}
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

## See Also

* [Physics 555: Bound States in the 1D Schrödinger Equation](
  https://physics-555-quantum-technologies.readthedocs.io/en/latest/Notes/Shooting.html)
  

[hydrogen atom]: <https://en.wikipedia.org/wiki/Hydrogen_atom>
[reduced Bohr radius]: <https://en.wikipedia.org/wiki/Bohr_radius#Reduced_Bohr_radius>
[generalized Laguerre polynomial]: <https://en.wikipedia.org/wiki/Laguerre_polynomial#Generalized_Laguerre_polynomials>
[spherical harmonic]: <https://en.wikipedia.org/wiki/Spherical_harmonics>
[finite difference methods]: <https://en.wikipedia.org/wiki/Finite_difference_method>
[Dirichlet boundary conditions]: <https://en.wikipedia.org/wiki/Dirichlet_boundary_condition>
