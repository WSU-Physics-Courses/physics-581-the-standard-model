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

Diagonalizing the Radial Schrödinger Equation
=============================================

Here we continue our previous explorations of numerical methods for solving the radial
Schrödinger Equation, trying to overcome some of the shortcomings seen there.  We will
focus on problems with a Coulomb potential that has both the mild singularity at $r=0$
and the long tail $V \sim 1/r$:
\begin{gather*}
  \Biggl(
    \frac{-\hbar^2}{2m}\biggl(\diff[2]{}{r} + \frac{l(l+1)}{r^2}\biggr) + V(r) - E
  \Biggr)u(r) = 0.
\end{gather*}

:::{margin}
Note that we can choose units where $\alpha = \hbar = m = 1$ so that
\begin{gather*}
  1 = \underbrace{a = \frac{\hbar^2}{m\alpha}}_{\text{length}}
    = \underbrace{\frac{\alpha}{l}}_{\text{energy}}.
\end{gather*}
There are no dimensionless parameters in this theory.
:::
## Coulomb Potential

For testing we use the exact energies for the non-relativistic [hydrogen atom][]:
\begin{gather*}
  V(r) = \frac{-\alpha}{r}, \qquad
  E_{ln} = \frac{-m\alpha^2/\hbar^2}{2(1+l+n)^2} = \frac{-\hbar^2}{2ma^2(1+l+n)^2},\\
  u_{ln}(r) = \sqrt{\left(\frac{2}{na}\right)^3\frac{(n-l-1)!}{2n(n+l)!}}
              re^{-r/na}
              \left(\frac{2r}{na}\right)^{l}
              L_{n-l-1}^{2l+1}\left(\frac{2r}{na}\right), \\
  a = \frac{\hbar^2}{m \alpha},\qquad
  L_{n}^{\alpha}(x) = \frac{x^{-\alpha}}{n!}\left(\diff{}{x} - 1\right)^{n}x^{n+\alpha}.
\end{gather*}
where $m$ is the reduced mass $m = m_em_p/(m_e+m_p)$, $a$ is the [reduced Bohr
radius][], $L_{n-l-1}^{2l+1}(\rho)$ is a [generalized Laguerre polynomial][] of degree 
$n-l-1$, and $Y_{l}^{m}(\theta, \phi)$ is a [spherical harmonic][] function of degree 
$l$ and order $m$. The full 3D eigenstates are:
\begin{gather*}
  \psi_{n,l,m}(r, \theta, \phi) u_{nl}(r)Y_{l}^{m}(\theta, \phi).
\end{gather*}
The states in these formulae are classified by the following quantum numbers:
* $n \in \{1, 2, 3, \dots\}$: [principal quantum number][],
* $l \in \{0, 1, 2, \dots, n-1\}$: [azimuthal quantum number][].  These are sometimes
  referred to by orbital terminology S ($l=0$), P ($l=1$), D ($l=2$) etc.  Here we focus
  on the S-wave properties $l=0$ which are spherically symmetric.
* $m \in \{-l, -l+1, \dots, l-1, l\}$: [magnetic quantum number][],

Note: Orthogonality between different values of $l$ is enforced by the spherical
harmonics.  The radial wavefunctions are orthogonal in the following sense:
\begin{gather*}
  \braket{u_{nl}|u_{ml}} = \int_0^{\infty} u_{nl}^{*}(r)u_{ml}(r)\d{r} = \delta_{mn}.
\end{gather*}
We now check these properties numerically.


```{code-cell}
:tags: [hide-input]

from scipy.special import genlaguerre as L, factorial
from scipy.integrate import quad

class Coulomb:
    m = alpha = hbar = 1
    a = hbar**2/m/alpha

    def get_u(self, r, n, l=0):
        a = self.a
        A = np.sqrt((2/n/a)**3*factorial(n-l-1)/2/n/factorial(n+l))
        rho = 2*r/n/a
        return A*r*np.exp(-rho/2)*rho**l*L(n-l-1, 2*l+1)(rho)

    def get_E(self, n, l=0):
        return -self.hbar**2/2/self.m/self.a**2/(1+l+n)**2

coulomb = Coulomb()

# Check orthonormality

N = 10
def f(r):
    return (coulomb.get_u(r, n, l).conj() * coulomb.get_u(r, m, l))

for n in range(1, N+1):
    for m in range(1, n+1):
        for l in range(min(m, n)):
            res = quad(f, 0, np.inf)[0]
            if m == n:
                assert np.allclose(res, 1)
            else:
                assert np.allclose(res, 0)

r = np.linspace(0, 60, 1000)[1:]
fig, ax = plt.subplots()
for l in range(3):
    for n in range(l+1, 5):
        u_r = coloumb.get_u(r, n=n, l=l)
        ax.plot(r, u_r, ls=['-', '--', ':', '-.'][l], c=f"C{n}",
                label=f"{n=}, {l=}")
ax.set(xlabel="$r$", ylabel="$u(r)$") 
ax.legend();
```
These exact solutions might be used to deal with the singular properties of the Coulomb
potential, but we postpone discussing this for now in favour of more general
techniques.

## Diagonalization

Our diagonalization approach did not work very well before, most likely because of the
singularities and long-range tail of the Coulomb potential.  One approach for dealing
with this is to use the exact solutions above as a basis, and to expand any other
terms.  This can work quite well, but computing the matrix elements 

```{code-cell} ipython3
import scipy.linalg
import scipy as sp

# Hydrogen atom
N = 64*8
R = 20.0
a = R/N
n = np.arange(N)
r = (n + 1)*a

N_states = 30
hbar = m = alpha = 1

Vr = -alpha/r
l = 0
En = -m*alpha**2/hbar**2/2/(1+l+n)**2

D2 = (np.diag(np.ones(N-1), k=1) 
      -2*np.diag(np.ones(N)) 
      + np.diag(np.ones(N-1), k=-1))

I = (np.diag(np.ones(N-1), k=1) 
     + 10*np.diag(np.ones(N)) 
     + np.diag(np.ones(N-1), k=-1))/12

V = np.diag(Vr)

fig, axs = plt.subplots(2, 1)

ax = axs[0]

# Finite Difference
H = -hbar**2/2/m * D2/a**2 + V
E, psi = sp.linalg.eig(H)
inds = np.argsort(E.real)
E, psi = E[inds], psi[:, inds]
ax.semilogy(n[:N_states], abs(E/En-1)[:N_states], label="Finite Difference")
print(En[:4])
print(E[:4].real)

axs[1].plot(r, abs(psi[:,0])**2)
axs[1].plot(r, abs(psi[:,1])**2)

# Numerov
H = -hbar**2/2/m * D2/a**2 + I @ V
E, psi = sp.linalg.eig(H, I)
inds = np.argsort(E.real)
E, psi = E[inds], psi[:, inds]
axs[1].plot(r, abs(psi[:,0])**2)
axs[1].plot(r, abs(psi[:,1])**2)
ax.semilogy(n[:N_states], abs(E/En-1)[:N_states], label="Numerov");
print(E[:4].real)
```

We now see that, although we qualitatively correct solutions, we only have 3 bound
states (because our box gets too small), and about 4 places of accuracy, even with a lot
of points.


## Change of Basis

:::{margin}
I think of the basis functions $\phi_{n}(r)$ as a "unitary matrix" $[\mat{\phi}]_{rn}$ whose
first index $r$ is continuous and whose second index $n$ is discrete.  This is unitary
in the sense of
\begin{gather*}
  \mat{\phi}^\dagger \mat{\phi} = \sum_{r}[\mat{\phi}^\dagger]_{mr}[\mat{\phi}]_{rn}\\
  = \int \d{r}\;w(r) \phi^*_{m}(r)\phi_{n}(r) = \delta_{mn}.
\end{gather*}
I.e.
\begin{gather*}
  \sum_{r} \equiv \int \d{r}\; w(r).
\end{gather*}
If the basis is complete, then we also have
\begin{gather*}
  \mat{\phi}\mat{\phi}^\dagger = \sum_{n}\phi_{n}(r)\phi^*_{n}(r') = \delta(r-r').
\end{gather*}
If the basis is finite, however, then this might only be true on some set of discrete
abscissa $r_n$. Now we write:
\begin{gather*}
  \op{H}u(r) = u(r)E, \\
  u(r) = \sum_{n}s_n\phi_{n}(r) = \mat{\phi}\ket{s},\\
  \underbrace{\mat{\phi}^\dagger \op{H}\mat{\phi}}_{\mat{H}}\ket{s} 
  = \mat{\phi}^\dagger\mat{\phi}\ket{s}E = \ket{s}E,\\
  \mat{H}\mat{S} = \mat{S}\mat{E}.
\end{gather*}
:::
A general strategy is to express the problem in an orthonormal basis $\ket{\phi_n}$ of
functions $\phi_n(r) = \braket{r|\phi_n}$,
\begin{gather*}
  \braket{\phi_m|\phi_n} = \int\d{r}\; w(r) \phi_m^*(r)\phi_n(r) = \delta_{mn},
\end{gather*}
using an inner product with a weight function $w(r)$.  E.g., in spherical coordinates we
might take $w(r) \propto r^2$.  One now simply computes the matrix elements of the
Hamiltonian and diagonalize:
\begin{gather*}
  \braket{\phi_m|\op{H}|\phi_n} = [\mat{H}]_{mn}
  = \int\d{r}\; w(r) \phi^*_m(r)
  \Biggl(
    \frac{-\hbar^2}{2m}\Bigl(\phi''_n(r) + \frac{l(l+1)}{r^2}\phi_n(r)\Bigr)
    + V(r)\phi_n(r)
  \Biggr),\\
  \mat{H} = \mat{S}\mat{E}\mat{S}^{-1},\qquad
  u_n(r) = \sum_{m}\phi_m(r)\mat{S}_{mn},
\end{gather*}
where $\mat{E} = \diag(E_0, E_1, \dots)$.  Depending on the complexity of the basis
functions, the terms in this integral may be easier or harder to compute.

An obvious application is to use the exact solutions to the Coulomb potential:
\begin{gather*}
  \op{H} = \op{H}_{C} + V_{s}(\op{r}).
\end{gather*}
Here we demonstrate this technique by computing the lowest 3 eigenvalues for
\begin{gather*}
  V(r) = \underbrace{\frac{-\alpha}{r}}_{V_C} + \underbrace{\frac{\alpha}{a}e^{-r^2/2a^2}}_{V_s}.
\end{gather*}

```{code-cell}
coulomb = Coulomb()

l = 0

def get_V(r):
    alpha, a = coulomb.a, coulomb.alpha
    return alpha/a * np.exp(-(r/a)**2/2)

N = 40
V = np.zeros((N, N))
ns = np.arange(1, N+1)

# These loops are quite slow...
for n in ns:
    for m in range(n, N+1):
        def f(r):
            un = coulomb.get_u(r, n, l=l)
            um = coulomb.get_u(r, m, l=l)
            return um.conj()*get_V(r)*un
        V[n-1, m-1] = V[m-1, n-1] = quad(f, 0, np.inf)[0]

H0 = np.diag(coulomb.get_E(n=ns, l=l))
H = H0 + V
E_ = np.linalg.eigvalsh(H)

Es = []
N0 = 3
Ns = 2**np.arange(1+int(np.log2(N0)), 1+int(np.log2(N)))
Es = np.array([np.linalg.eigvalsh(H[:N, :N])[:N0] for N in Ns])
fig, ax = plt.subplots()
ax.loglog(Ns, abs(Es - E_[:N0]), '-+')
ax.set(xlabel="$N$", ylabel="$E$ error", xticks=Ns);
```

We see that this basis does a reasonable job, but the integrals are quite expensive to
compute.  For a one-off computation, this should work to generate some data for
analysis.  *(Note: we have only crudely estimated the "correct" answers here by
computing in a larger basis, so take the observed scaling with a grain of salt.)*

## DVR Basis for the Radial Equation

Another approach for diagonalization is to use a discrete variable representation
([DVR][]) basis (sometimes called [pseudo-spectral methods][DVR]).  These can be
constructed for any [orthogonal polynomials][] (see {cite}`Baye:1986si`) but some care
is needed when applying these.

The idea is define a set of $N$ basis functions $\phi_{n}(x)$ that are orthogonorma on
some range $[a, b]$ with a weight function $w(x)$:
\begin{gather*}
  \int_{a}^{b}\d{x}\; w(x)\phi_m^*(x)\phi_n(x) = \delta_{mn}.
\end{gather*}
Additionally, we need a [quadrature rule][numerical integration] with $N$ points $x_i$
and weights $w_i$ such that





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
[numerical integration]: <https://en.wikipedia.org/wiki/Numerical_integration>
[DVR]: <https://en.wikipedia.org/wiki/Pseudo-spectral_method>
[orthogonal polynomials]: <https://en.wikipedia.org/wiki/Orthogonal_polynomials>
[principal quantum number]: <https://en.wikipedia.org/wiki/Principal_quantum_number>
[azimuthal quantum number]: <https://en.wikipedia.org/wiki/Azimuthal_quantum_number>
[magnetic quantum number]: <https://en.wikipedia.org/wiki/Magnetic_quantum_number>
