---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.0
kernelspec:
  display_name: Python 3 (phys-581)
  language: python
  name: phys-581
---

```{code-cell}
:tags: [hide-cell]

import mmf_setup;mmf_setup.nbinit()
import logging;logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
```

(sec:DiracEq)=
# Dirac Equation

## Rotations
:::{margin}
Here $\varepsilon_{ijk}$ is the [Levi-Civita symbol][] which effects the cross product:
\begin{gather*}
  [\vect{\theta}\times\vect{B}]_i = \varepsilon_{iaj}\theta_aB_j = [\mat{\vect{\theta}\times}]_{ij}B_j.
\end{gather*}
:::
We start with some group theory.  Active rotations about the axis $\vect{\theta}$ of
magnitude $\theta = \abs{\theta}$ in 3D can be effected by the following linear
transformations, which form a faithful representation of the SO(3) group:
\begin{gather*}
  \mat{R}_{\vect{\theta}} = e^{\mat{\vect{\theta}\times}}, \qquad
  [\mat{\vect{\theta}\times}]_{ij} = \varepsilon_{iaj}\theta_a.
\end{gather*}
From this, we can determine the elements of the corresponding Lie Algebra by expanding
about the origin to linear order in $\theta$:
\begin{gather*}
  \mat{R}_{\vect{\theta}} = e^{\mat{\vect{\theta}\times}} =
  \mat{1} + \mat{\vect{\theta}\times} + O(\theta^2).
\end{gather*}
:::{margin}
I remember these as follows: $\mat{T}_{k}$ has zeros on both rows and columns
$i=j=k$.  The other $2\times 2$ block is anti-symmetric with $\mat{T}_{y}$ having
the opposite sign as in the formula for a determinant. To get the sign, act on a vector:
Under active rotations,
$\mat{T}_{z}\uvect{x} = \uvect{y}$ while $\mat{T}_{z}\uvect{y} = -\uvect{y}$:
<!-- \begin{gather*} -->
<!--   \mat{R}_{\theta\uvect{z}}\begin{pmatrix} -->
<!--     1 & 0\\ -->
<!--     0 & 1\\ -->
<!--     0 & 0 -->
<!--   \end{pmatrix} -->
<!--   = -->
<!--   \mat{R}_{\theta\uvect{z}}\begin{pmatrix} -->
<!--     \cos\theta & -\sin\theta\\ -->
<!--     \sin\theta & \cos\theta\\ -->
<!--     0 & 0 -->
<!--   \end{pmatrix}. -->
<!-- \end{gather*} -->
:::
The three generators are $[\mat{T}_{k}]_{ij} = \varepsilon_{ikj}$:
\begin{gather*}
  \mat{T}_{x} = \begin{pmatrix}
    0 & 0 & 0\\
    0 & 0 & -1\\
    0 & 1 & 0
  \end{pmatrix}, \qquad
  \mat{T}_{y} = \begin{pmatrix}
    0 & 0 & 1\\
    0 & 0 & 0\\
    -1 & 0 & 0
  \end{pmatrix}, \qquad
  \mat{T}_{z} = \begin{pmatrix}
    0 & -1 & 0\\
    1 & 0 & 0\\
    0 & 0 & 0
  \end{pmatrix}.
\end{gather*}
Rotations are obtained from the algebra by exponentiating:
\begin{gather*}
  \mat{R}_{\vect{\theta}} = e^{\vect{\theta}\cdot\vect{\mat{T}}}.
\end{gather*}
For rotations in 3D, this notation make sense, but in particle physics a
different convention is used that more closely connects with quantum mechanics.
The generators $\mat{T}_a = - \mat{T}_a^T$ are anti-symmetric, but in physics, these
will correspond to observables like angular momentum which need to be Hermitian.

Thus, the convention in particle physics is to include a factor of $\I$:
\begin{gather*}
  \mat{L}_{a} = \I\mat{T}_{a},\qquad
  \mat{R}_{\vect{\theta}} 
  = e^{\vect{\theta}\cdot\vect{\mat{T}}}
  = e^{\vect{\theta}\cdot\vect{\mat{L}}/\I}.
\end{gather*}
The resulting generators are now Hermitian:
\begin{gather*}
  \mat{L}_{x} = \I\begin{pmatrix}
    0 & 0 & 0\\
    0 & 0 & -1\\
    0 & 1 & 0
  \end{pmatrix}, \qquad
  \mat{L}_{y} = \I\begin{pmatrix}
    0 & 0 & 1\\
    0 & 0 & 0\\
    -1 & 0 & 0
  \end{pmatrix}, \qquad
  \mat{L}_{z} = \I\begin{pmatrix}
    0 & -1 & 0\\
    1 & 0 & 0\\
    0 & 0 & 0
  \end{pmatrix}.
\end{gather*}
With an additional factor of $\hbar$ (which we set to one in this course), these form a
representation of [angular momentum][].

The group elements are formed from the algebra by exponentiating, which looks like this:
\begin{gather*}
  \mat{R}_{\vect{θ}} 
  = e^{\vect{θ}\cdot\vect{\mat{T}}} 
  = e^{\vect{θ}\cdot\vect{\mat{L}}/\I}
\end{gather*}
From now on, we will follow the physics convention.  Note that the generators can be
obtained by expanding the rotation to linear order:
\begin{gather*}
  \mat{R}_{\vect{\theta}}
  = \mat{1} + \vect{\theta}\cdot\frac{\mat{L}}{\I} + O(\theta^2).
\end{gather*}
For example, a rotation about the $z$ axis has the form:
\begin{gather*}
  \mat{R}_{\theta \uvect{z}} 
  = \begin{pmatrix}
    \cos\theta & -\sin\theta & 0\\
    \sin \theta & \cos \theta & 0\\
    0 & 0 & 1
  \end{pmatrix}
  = \mat{1}
  + \theta
    \underbrace{
    \begin{pmatrix}
      0 & -1 & 0\\
      1\theta & 0 & 0\\
      0 & 0 & 0
    \end{pmatrix}
  }_{\mat{L}_x/\I}
  + O(\theta)^2.
\end{gather*}





### Lie Algebras
These matrices form a basis for the [Lie algebra][], which can be defined in terms of
the [structure constants][]:
\begin{gather*}
  [\mat{L}_{a}, \mat{L}_{b}] = \I f_{abc}\mat{L}_{c}.
\end{gather*}
For $\mathfrak{so}(3)$, we have $f_{abc} = \epsilon_{abc}$.

:::{margin}
The raising and lowing of indices here has no specific significance, but is common so
that summed indices are of mixed height.  This will be important later when we use the
relativistic metric.
:::
:::{admonition} Convention in Mathematics
In terms of the anti-symmetric matrices $\mat{T}_{a}$, we have
\begin{gather*}
  [\mat{T}_{a}, \mat{T}_{b}] = c_{ab}^{}{c} \mat{T}_{c}.
\end{gather*}
For rotations in 3D, this notation make sense.

Each has their use.  For example, since the structure constants are real, if one has a
complex representation $\mat{T}_{a}$ *(not purely real)*, then one can
form another complex representation by conjugating $\overline{\mat{T}}_{a}$.  In the
second formulation, this requires an additional sign, $-\overline{\mat{L}}_{a}$, which
inequivalent only if $\mat{L}_{a}$ is not purely imaginary.
:::

### Trivial Representation
We start with the [trivial representation][], which always exists
\begin{gather*}
  \mat{L}_k = \mat{0}.
\end{gather*}
This is he only one-dimensional representation, since finite numbers compute.  It is thus referred to by the number $\mathbf{1}$.

### Adjoint Representation
:::{margin}
Note that the [adjoint representation][] is sometimes written:
\begin{gather*}
  [\mat{L}_{a}]_{bc} = \I f_{acb}
\end{gather*}
with the indices reversed.  This simply means that if $\mat{L}_{a}$ is a representation,
then so is $-\mat{L}_{a}^{T}$, which follows from taking the transpose of the structure
equation:
\begin{gather*}
  [\mat{L}_{a}, \mat{L}_{b}]^T 
  = \I f_{abc} \mat{L}_{c}^T\\
  = [\mat{L}_{b}^T, \mat{L}_{a}^T]
  = -[\mat{L}_{a}^T, \mat{L}_{b}^T]\\
  [-\mat{L}_{a}^T, -\mat{L}_{b}^T] = \I f_{abc} (-\mat{L}_{c}^T).
\end{gather*}
:::
The [structure constants][] also form a matrix representation $\mat{L}_{a}$ called the
[adjoint representation][] (we follow the form in {cite-p}`Georgi:2019`):
\begin{gather*}
  [\mat{L}_{a}]_{bc} = -\I f_{abc} = -\I \epsilon_{abc}.
\end{gather*}
It is purely imaginary and has the same dimension as the number of generators.  For
$\mathfrak{so}(3)$ this is $3$ and so this representation is often referred to
by the number $\mathbf{3}$.

::::{doit} Do It!
Show that the [adjoint representation][] defined by the structure factors indeed
satisfies the algebra.  
*Hint: use the [Jacobi identity][]: $\bigl[\mat{A}, [\mat{B}, \mat{C}]\bigr]
  +
  \bigl[\mat{B}, [\mat{C}, \mat{A}]\bigr]
  +
  \bigl[\mat{C}, [\mat{A}, \mat{B}]\bigr]
  = \mat{0}$.*
::::
::::{solution}
First we compute:
\begin{gather*}
  \bigl[\mat{L}_{a}, [\mat{L}_{b}, \mat{L}_{c}]\bigr]
  = \I f_{bcd} [\mat{L}_{a}, \mat{L}_{d}] 
  = - f_{bcd}f_{ade}\mat{L}_{e}.
\end{gather*}
Summing over permutations of $\{abc, bca, cab\}$ we can express [Jacobi identity][] as
\begin{gather*}
  f_{bcd}f_{ade} + f_{abd}f_{cde} + f_{cad}f_{bde} = 0.
\end{gather*}
Now using the [adjoint representation][] $[\mat{L}_{a}]_{bc} = -\I f_{abc}$ and the
anti-symmetry $f_{abc} = -f_{bac}$ implied by the commutator, we can rearrange this as
\begin{gather*}
  \underbrace{
    \underbrace{(-\I f_{bcd})}_{[\mat{L}_{b}]_{cd}}
    \underbrace{(-\I f_{ade})}_{[\mat{L}_{a}]_{de}}
    - 
    \underbrace{(-\I f_{acd})}_{[\mat{L}_{a}]_{cd}}
    \underbrace{(-\I f_{bde})}_{[\mat{L}_{b}]_{de}}
  }_{-[\mat{L}_{b}, \mat{L}_{a}]_{ce}}
  = -\I f_{abd} \underbrace{(-\I f_{dce})}_{[\mat{L}_{d}]_{ce}},\\
  [\mat{L}_{a}, \mat{L}_{b}] = \I f_{abd}\mat{L}_{d}.
\end{gather*}
This proves that the [adjoint representation][] satisfies the algebra.
::::

```{code-cell}
:tags: [hide-cell]

from scipy.linalg import expm
eps = np.zeros((3, 3, 3)) # Levi-Civita symbol.
eps[0,1,2] = eps[1,2,0] = eps[2,0,1] = 1
eps[2,1,0] = eps[1,0,2] = eps[0,2,1] = -1
Tx, Ty, Tz = T = np.einsum('iaj->aij', eps)
Lx, Ly, Lz = L = 1j*T

def com(A, B):
    return A@B - B@A

assert np.allclose(com(Tx, Ty), Tz)
assert np.allclose(com(Ty, Tz), Tx)
assert np.allclose(com(Tz, Tx), Ty)
assert np.allclose(com(Lx, Ly), 1j*Lz)
assert np.allclose(com(Ly, Lz), 1j*Lx)
assert np.allclose(com(Lz, Lx), 1j*Ly)
print(T.astype(int))
```

### Pauli Matrices

The Pauli matrices provide another two-dimensional representation of
$\mathfrak{so}(3)$.  To see this, we note:
\begin{gather*}
  \mat{\sigma}_{x} = \begin{pmatrix}
    0 & 1\\
    1 & 0
  \end{pmatrix}, \qquad
  \mat{\sigma}_{y} = \begin{pmatrix}
    0 & -\I\\
    \I & 0
  \end{pmatrix}, \qquad
  \mat{\sigma}_{z} = \begin{pmatrix}
    1 & 0\\
    0 & -1
  \end{pmatrix}
\end{gather*}
from which the following properties can be deduced:
\begin{gather*}
  \mat{σ}_{a}\mat{σ}_{b} = \delta_{ab}\mat{1} + \I\varepsilon_{abc}\mat{\sigma}_c,\qquad
  \Tr\mat{σ}_{a} = 0,\\
  [\mat{\sigma}_{a}, \mat{\sigma}_{b}] = 2\I\varepsilon_{abc}\mat{\sigma}_c, \qquad
  \{\mat{\sigma}_{a}, \mat{\sigma}_{b}\} = 2\delta_{ab}\mat{1}.
\end{gather*}
Scaling appropriately, we thus see that the following set of matrices form a
2-dimensional representation of $\mathfrak{s0}(3)$:
\begin{gather*}
  \left[\frac{\mat{\sigma}_{a}}{2}, \frac{\mat{\sigma}_{b}}{2}\right] = \I \varepsilon_{abc}\frac{\mat{\sigma}_c}{2}, \qquad
  \mat{L}_{a} = \frac{\mat{\sigma}_{a}}{2}.
\end{gather*}
This is a complex representation referred to by the number $\textbf{2}$.

As mentioned above, since the [structure constants][] are real, complex representations
like this always appear in pairs:
\begin{gather*}
  [-\overline{\mat{L}}_{a}, -\overline{\mat{L}}_{b}] = \I f_{abc}(-\overline{\mat{L}}_{c}).
\end{gather*}
For $\mathfrak{so}(3)$, this representation is called $\mathbf{\bar{2}}$:
\begin{gather*}
  \mat{L}_{a} = \frac{-\overline{\mat{\sigma}}_{a}}{2}.
\end{gather*}

:::{doit} Show that $\mathbf{2}$ and $\mathbf{\bar{2}}$ are equivalent.

I.e. find a matrix $\mat{S}$ such that
\begin{gather*}
  \mat{S}\mat{\sigma}_{a}\mat{S}^{-1} = -\overline{\mat{\sigma}_{a}}.
\end{gather*}
:::

:::{solution}
\begin{gather*}
  \mat{S} = \begin{pmatrix}
    0 & 1\\
    -1 & 0
  \end{pmatrix}, \qquad
  \mat{S}^{-1} = \mat{S}^T = -\mat{S},\qquad
  [\mat{S}]_{ij} = \varepsilon_{ij}.
\end{gather*}
This is the 2-dimensional [Levi-Civita symbol][].
:::
```{code-cell}
:tags: [hide-input]

σ = np.array([
    [[0, 1],
     [1, 0]],
    [[0, -1j],
     [1j, 0]],
    [[1, 0],
     [0, -1]]])
σbar = -σ.conj() 
assert np.allclose(σ[0]@σ[1], 1j*σ[2])
S = np.array([
    [0, 1],
    [-1, 0]])
Sinv = S.T
assert np.allclose(S @ Sinv, np.eye(2))
assert np.allclose(S @ σ @ S.T, σbar)
```

### Extended Pauli Matrices

Note that, if we include the identity $\mat{σ}_{0} = \mat{1}$, then [Pauli matrices][]
form a basis for Hermitian matrices:
\begin{gather*}
  \mat{A} = \mat{A}^\dagger = a^{\mu}\mat{\sigma}_\mu.
\end{gather*}


Exponentiating the $2$ and $\bar{2}$ representations, we 

# Lorentz Group

In addition to rotations $\vect{\theta}$, the Lorentz group has boosts.  A boost along
the $x$ axis with [rapidity][] $\eta$ has the form
\begin{gather*}
  \mat{B}_{\eta\hat{x}} = \I \begin{pmatrix}
    \cosh\eta & \sinh \eta & 0 & 0\\
    \sinh \eta & \cosh\eta & 0 & 0\\
    0 & 0 & 1 & 0\\
    0 & 0 & 0 & 1
  \end{pmatrix}
  = e^{\eta\mat{K}_x/\I}
\end{gather*}
:::{margin}
I.e.:
\begin{gather*}
  \mat{J}_{a} = \begin{pmatrix}
  1 \\
  & \mat{L}_{a}
  \end{pmatrix}.
\end{gather*}
:::
The corresponding generators are
\begin{gather*}
  \mat{J}_{x} = \I\begin{pmatrix}
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & -1\\
    0 & 0 & 1 & 0
  \end{pmatrix}, \qquad
  \mat{J}_{y} = \I\begin{pmatrix}
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 1\\
    0 & 0 & 0 & 0\\
    0 & -1 & 0 & 0
  \end{pmatrix}, \qquad
  \mat{J}_{z} = \I\begin{pmatrix}
    0 & 0 & 0 & 0\\
    0 & 0 & -1 & 0\\
    0 & 1 & 0 & 0\\
    0 & 0 & 0 & 0
  \end{pmatrix},\\
  \mat{K}_{x} = \I 
  \begin{pmatrix}
    0 & 1 & 0 & 0\\
    1 & 0 & 0 & 0\\
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0
  \end{pmatrix}, \qquad
  \mat{K}_{y} = \I 
  \begin{pmatrix}
    0 & 0 & 1 & 0\\
    0 & 0 & 0 & 0\\
    1 & 0 & 0 & 0\\
    0 & 0 & 0 & 0
  \end{pmatrix}, \qquad
  \mat{K}_{z} = \I 
  \begin{pmatrix}
    0 & 0 & 0 & 1\\
    0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0\\
    1 & 0 & 0 & 0
  \end{pmatrix},\\
\end{gather*}
Together, these generate the proper orthochronos Lorentz group:
\begin{gather*}
  \mat{\Lambda}(\vect{\theta}, \vect{\eta}) 
  = e^{(\vect{\theta}\cdot\vect{\mat{J}} + \vect{\eta}\cdot\vect{\mat{K}})/\I},\qquad
  \vect{\eta} = \uvect{v}\tanh^{-1}\frac{v}{c},
\end{gather*}
where the algebra is defined by
\begin{gather*}
  [\mat{J}_{i}, \mat{J}_{j}] = \I \varepsilon_{ijk}\mat{J}_{k}, \qquad
  [\mat{K}_{i}, \mat{J}_{j}] = \I \varepsilon_{ijk}\mat{K}_{k}, \qquad
  [\mat{K}_{i}, \mat{K}_{j}] = -\I \varepsilon_{ijk}\mat{J}_{k}.
\end{gather*}
The defining representation given above defines the 4×4 space-time transformation, but
quantum fields can transform under different representations.  To classify these, note
that these six generators can be rearranged to form a computing
$\mathfrak{su}(2)_L\times\mathfrak{su}(2)_R$ algebra:
\begin{gather*}
  \vect{\mat{J}}^{\pm} = \tfrac{1}{2}(\vect{\mat{J}} \pm \I \vect{\mat{K}}), \qquad
  [\mat{J}^{\pm}_{i}, \mat{J}^{\pm}_{j}] = \I\varepsilon_{ijk}\mat{J}^{\pm}_{k}, \qquad
  [\mat{J}^{+}_{i}, \mat{J}^{-}_{j}] = \mat{0}.
\end{gather*}
:::{margin}
Note that the dimension of the $\mathfrak{su}(2)$ representation with Casimir $j$ is $d
= 1+2j$.  Thus, $j=0$ is 1 dimensional (scalar), $j=\tfrac{1}{2}$ is 2-dimensional
(spinor), $j=1$ is 3-dimensional (vector) etc.
:::
Thus, we can label the representations by the pair $(j_+, j_-)$ where $j_{\pm} = 0,
\tfrac{1}{2}, 1, \cdots$ are defined by the Casimirs
\begin{gather*}
  \vect{J}_{\pm}^2 = j_{\pm}(j_{\pm}+1).
\end{gather*}
The spin-1/2 representations are $\vect{\mat{J}} = \vect{\mat{\sigma}}/2$ and
$\vect{\mat{K}} = \pm \vect{\mat{\sigma}}/2\I$:
\begin{gather*}
  \newcommand{\spin}{{\small\frac{1}{2}}}
  \mat{\Lambda}_{(\spin,0)} 
  = e^{(\I\vect{\theta} + \vect{\eta})\cdot\vect{\mat{\sigma}}/2}, \qquad
  \mat{\Lambda}_{(0, \spin)} 
  = e^{(\I\vect{\theta} - \vect{\eta})\cdot\vect{\mat{\sigma}}/2}
  = \mat{\Lambda}_{(\spin,0)}^{-1\dagger}, \\
  \mat{\Lambda}_{(\spin^{*},0)} 
  = e^{(-\I\vect{\theta} - \vect{\eta})\cdot\vect{\mat{\sigma}^*}/2}
  = \mat{\Lambda}_{(\spin,0)}^{-1T}, \qquad
  \mat{\Lambda}_{(0, \spin^*)} 
  = e^{(-\I\vect{\theta} + \vect{\eta})\cdot\vect{\mat{\sigma}}/2}
  = \mat{\Lambda}_{(\spin,0)}^{*}
  = \mat{\Lambda}_{({\spin},0)}^{\dagger T}.
\end{gather*}




## Covariant Formulation

To formulate things covariantly, we define the following four-vectors (now choosing
units such that $c = \hbar = 1$):
\begin{gather*}
  A^{\mu} = (A^0, \vect{A}), \qquad x^{\mu} = (t, \vect{x}), \qquad
  p^{\mu} (E, \vect{p}).
\end{gather*}
Note that the contravariant four-vectors, with a raised index, have the corresponding
physical quantities.  We now introduce the metric
\begin{gather*}
  g_{\mu\nu} = g^{\mu\nu} = \diag(1, -1, -1, -1).
\end{gather*}
We follow the conventions in {cite-p}`Langacker:2017`:

:::{figure} _images/LangackerTable1.2.png
:name: tab:LangackerTable1.2

Table 1.2 from {cite-p}`Langacker:2017`.  Note how the metric $g_{\mu\nu}$ is used to
raise and lower indices, changing the sign of the spatial part in the process.
:::

Using these expressions, we can express the Lorentz transform more covariantly by
packaging the angles and rapidities into 4×4 tensor $\omega_{\rho \sigma}$:
\begin{gather*}
  \theta^{i} = \tfrac{1}{2}\varepsilon^{ijk} \omega_{jk}, \qquad
  \eta^{i} = \omega^{0i},\\
  \omega_{\mu\nu} = \begin{pmatrix}
    0 & \vect{\eta}\\
    \vect{\eta}^T & \vect{\theta}\cdot\vect{\mat{T}}
  \end{pmatrix}
  = \begin{pmatrix}
    0 & \eta_x & \eta_y & \eta_z\\
    \eta_x & 0 & \theta_z & -\theta_y\\
    \eta_y & -\theta_z & 0 & \theta_x\\
    \eta_z & \theta_y & -\theta_x & 0\\
  \end{pmatrix}.
\end{gather*}
We similarly package the matrices $\vect{\mat{J}}$ and $\vect{\mat{K}}$ into a 4×4 of
matrices $\mat{M}_{\rho\sigma}$:
\begin{gather*}
  \mat{J}^{i} = \tfrac{1}{2}\varepsilon^{ijk}\mat{M}_{jk}, \qquad
  \mat{K}^{i} = \mat{M}^{0i}, \\
  \tfrac{1}{2}\omega^{\rho\sigma}\mat{M}_{\rho\sigma} 
  = \vect{\theta}\cdot\vect{\mat{J}} + \vect{\eta}\cdot\vect{\mat{K}}.
\end{gather*}
Including the generators of translations $\mat{P}_{\mu}$ we have the full Poincaré
group, which has the algebra
\begin{gather*}
  [\mat{P}_{\mu}, \mat{P}_{\nu}] = 0, \qquad
  [\mat{M}_{\mu\nu}, \mat{P}_{\rho}] 
  = -\I(g_{\mu\rho}\mat{P}_{\nu} - g_{\nu\sigma}\mat{P}_{\rho}),\\
  [\mat{M}_{\mu\nu}, \mat{M}_{\rho\sigma}] =  
  \I(g_{\nu\rho}\mat{M}_{\mu\sigma} - g_{\nu\sigma}\mat{M}_{\mu\rho} 
  - g_{\mu\rho}\mat{M}_{\nu\sigma} + g_{\mu\sigma}\mat{M}_{\nu\rho}).
\end{gather*}







```{code-cell}

#:tags: [hide-cell]
J = np.zeros((3, 4, 4), dtype=complex)
J[:, 1:, 1:] = L
K = np.zeros((3, 4, 4), dtype=complex)
for i in [0, 1, 2]:
    K[i, 0, i+1] = K[i, i+1, 0] = 1j
Jx, Jy, Jz = J
Kx, Ky, Kz = K

eta = 0.1
c, s = np.cosh(eta), np.sinh(eta)
assert np.allclose(
    expm(eta*K[0]/1j),
    np.array([
        [c, s, 0, 0],
        [s, c, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]]))

def com_(A, B):
    """Return the full set of commutators."""
    return (np.einsum('iab,jbc->ijac', A, B) 
            - np.einsum('ibc,jab->ijac', A, B))

def eps_(A):
    """Return eps_{ijk}A_k."""
    return np.einsum('ijk,kab->ijab', eps, A)
    
assert np.allclose(com_(J, J), 1j*eps_(J))
assert np.allclose(com_(K, J), 1j*eps_(K))
assert np.allclose(com_(K, K), -1j*eps_(J))
```




[moments]: <https://en.wikipedia.org/wiki/Moment_(physics)>
[Levi-Civita symbol]: <https://en.wikipedia.org/wiki/Levi-Civita_symbol>
[structure constants]: <https://en.wikipedia.org/wiki/Structure_constants>
[adjoint representation]: <https://en.wikipedia.org/wiki/Adjoint_representation#Structure_constants>
[Jacobi identity]: <https://en.wikipedia.org/wiki/Jacobi_identity>
[Pauli matrices]: <https://en.wikipedia.org/wiki/Pauli_matrices>
[trivial representation]: <https://en.wikipedia.org/wiki/Trivial_representation>
[Lie algebra]: <https://en.wikipedia.org/wiki/Lie_algebra>
[angular momentum]: <https://en.wikipedia.org/wiki/Angular_momentum>



