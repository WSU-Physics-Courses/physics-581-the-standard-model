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

```{code-cell} ipython3
:tags: [hide-cell]

import mmf_setup;mmf_setup.nbinit()
import logging;logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
```

(sec:Electromagnetism)=
# Electromagnetism (Classical)

Here we follow the natural approach used to specify the Standard Model in order to
see the consequences that lead to Maxwell's equations for electromagnetic fields.

The approach we shall use is the following:

1. Establish an appropriate field with appropriate symmetry properties.
2. Write the most general Lagrangian density we can consider that satisfies the
   symmetries and limit this to the lowest non-trivial order.
3. Derive the equations of motions and see the consequencies.

Here is how it will play out:

*  We will consider a vector field $A^{\mu} = (\phi, \vect{A})$ that will ultimately
   contain the scalar and vector potentials of classical E&M.  This vector field will
   transform under the same [adjoint representation][] of the Lorentz group as other
   four-vectors:
   \begin{gather*}
     A^\mu \rightarrow \Lambda^{\mu}{}_{\nu}A^{\nu}.
   \end{gather*}
   After writing down the lowest order Lagrangian, we shall see that this gives rise to
   negative energy modes.  Note: these are not negative energy in the same sense as
   anti-particles, but fundamental -- quantizing without any further restrictions would
   lead to a Hamiltonian that is unbounded from below.
*  Currently, the only know solution, is to impose an additional symmetry called
   gauge-symmetry that renders these negative-energy modes unphysical.
   
## Lorentz Invariant Terms

:::{margin}
In terms of indices, we have
\begin{gather*}
  A_{\mu}A^{\mu} = A^{\mu}g_{\mu\nu}A^{\nu},\\
  A^{\mu} \rightarrow \Lambda^{\mu}{}_{\nu}A^{\nu},\\
  A_{\mu}A^{\mu} \rightarrow \Lambda_{\mu}{}_{\sigma}\Lambda^{\mu}{}_{\nu}A^{\nu}A^{\sigma},\\
 \Lambda_{\mu}{}_{\sigma}\Lambda^{\mu}{}_{\nu} = 
 \Lambda^{\alpha}{}_{\sigma}g_{\alpha\mu}\Lambda^{\mu}{}_{\nu} = g_{\sigma\nu},\\
 [\mat{\Lambda}^T]_{\sigma}{}^{\alpha}[\mat{g}]_{\alpha\mu}[\mat{\Lambda}]^{\mu}{}_{\nu} 
 = [\mat{g}]_{\sigma\nu}.
\end{gather*}
*Note that the metric $\mat{g}$ is not a matrix like $\mat{\Lambda}$: it makes vectors to
co-vectors, whereas matrices map vectors to vectors.*
:::
It is easy to construct some Lorentz-invariant terms by using the property of the
[adjoint representation][] that it leaves the metric invariant:
\begin{gather*}
  \mat{\Lambda}^T \mat{g} \mat{\Lambda} = \mat{g}.
\end{gather*}
The metric $\mat{g}$ is built into the raising and lowering of indices so that invariant
quantities can be formed by simply ensuring that all indices are contracted.  This
allows us to form Lorentz-invariant quantities (scalars) like these:
:::{margin}
It turns out that one can form additional Lorentz-invariant quantities by using the
conjugate representations and the 4D [Levi-Civita symbol][].  We will come back to this
later.

Note also that
\begin{multline*}
  \tfrac{1}{2}\partial^{2}A^{2}
  = \\
  (\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu}) \\
  + A^{\nu}\partial^{2}A_{\nu}.
\end{multline*}
:::

\begin{gather*}
  \partial_{\mu}A^{\mu}, \qquad
  A_{\mu}A^{\mu} \equiv A^2, \\
  \partial^2 A^2, \qquad
  (\partial_{\mu}A^{\mu})^2, \qquad
  (\partial_{\mu}A_{\nu})(\partial^{\nu}A^{\mu}), \qquad
  (\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu}).
\end{gather*}
The first term is the only way to make a scalar with a single $A^{\mu}$, but on its own,
is a total derivative: hence it plays no role in the local equations of motion.  It can
be used in combination with other terms, however.
The second term is valid and will turn out to give a mass to the corresponding
particle.

The remaining terms on the next line are the only other terms that can appear if we
restrict ourselves to at most quadratic order in $A^{\mu}$ and quadratic order in
derivatives $\partial^{\mu}$ (for reasons of renormalizability as we shall discuss
later).  If these appear on their own, then we can further reduce this to two
independent terms since the first is also a total derivative, and the second two are
equivalent after integrating by parts.

:::{margin}
Here we follow the argument given in {cite:p}`t-Hooft:2016` §4, but with our metric conventions.
:::
Thus, to this lowest order, we have the following Lagrangian density:
\begin{gather*}
  \mathcal{L}(A^{\mu}, \partial_{\mu}A^{\nu}) = 
  -\tfrac{1}{2}\alpha (\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu})
  + \tfrac{1}{2}\beta (\partial_{\mu}A^{\mu})^2
  + \tfrac{1}{2} m^2 A^{\mu}A_{\mu}
  - A^{\mu}J_{\mu}.
\end{gather*}
Following the standard procedure, we find the equations of motion
\begin{gather*}
  \partial_{\mu} \pdiff{L}{(\partial_{\mu}A_{\nu})} = \pdiff{L}{A_{\nu}},\qquad
  (\beta \partial^{\nu}\partial_{\mu} - \alpha \delta^{\nu}_{\mu} \partial^2)A^{\mu}
  = m^2 A^{\nu} - J^{\nu}.
\end{gather*}
:::{margin}
\begin{gather*}
  \dot{A}^{\mu} \equiv \partial_0 A^{\mu}.
\end{gather*}
:::
We can also compute the canonical momentum $\pi_{\mu}$ and form the Hamiltonian density
\begin{gather*}
  \pi_{\mu} = \pdiff{L}{\dot{A}^{\mu}} = 
  \beta\delta^{0}_{\mu}\partial_{\nu}A^{\nu}
  -\alpha\dot{A}_{\mu},\qquad
  \mathcal{H} = \pi_{\mu}\dot{A}^{\mu} - \mathcal{L}.
\end{gather*}
Now consider the UV limit of high energy and momentum.  In this limit, the mass and
current terms become negligible and we can take $m, J^{\mu} \rightarrow 0$.  Now
consider the form of the Hamiltonian in two limits:

1. First consider when all space-like derivatives $\partial_{i} = \nabla^i = \vect{\nabla}$ and
   spacelike components $A^{i} = \vect{A}$ are small compared with $\dot{A}^{0}$.
   In this limit
   \begin{gather*}
     \mathcal{H} \rightarrow \tfrac{1}{2}(\beta - \alpha)(\dot{A}^0)^2.
   \end{gather*}
   Positivity of the Hamiltonian density in this limit thus requires $\beta \geq \alpha$.
   :::{margin}
   The first term here from $-\mathcal{L}$ is positive due to the metric.  In
   these expressions, there are no more signs from the metric.
   :::
2. Now consider the opposite case where $A^{0}$ and all time-like derivatives are negligible:
   \begin{gather*}
     \mathcal{H} \rightarrow 
     \tfrac{1}{2}\alpha (\nabla^{i}A^{j})^2
     -
     \tfrac{1}{2}\beta(\vect{\nabla}\cdot \vect{A})^2.
   \end{gather*}
   Positivity requires $\alpha \geq \beta$.
Thus, we conclude that $\alpha = \beta$.
::::{doit} Do It!
Show that the second condition implies $\alpha \geq \beta$.
::::
::::{solution}
:show:

We write this as a quadratic form:
\begin{gather*}
  \mathcal{H} \rightarrow \tfrac{1}{2}(\nabla^{i}A^{j})M_{ij;ab}(\nabla^{a}A^{b}), 
  \qquad
  M_{ij;ab} = \alpha\delta_{ai}\delta_{bj} - \beta\delta_{ij}\delta_{ab}.
\end{gather*}
This matrix $\mat{M}$ has an eigenvector $x_{ab} = \delta_{ab}$ with eigenvalue $\alpha -
3\beta$.  Positivity thus requires $\alpha \geq 3\beta$.
\begin{gather*}
  x_{ab} = \delta_{ab}, \\
  (\alpha-3\beta)\delta_{ij}
\end{gather*}



\begin{gather*}
  a A_{i,j}A_{i,j} - b A_{i, i}A_{j,j}
\end{gather*}
::::

```{code-cell} ipython3
a = 1.23
b = 0.0000456

I = np.eye(3)

M = (a*np.einsum('ai,bj->ijab', I, I) - b*np.einsum('ab,ij->ijab', I, I)).reshape(9, 9)
np.linalg.eigvalsh(M)
```

::::{admonition} Details
Since we had some difficulties, here is a careful reckoning.  We start with the standard form
\begin{gather*}
  \mathcal{L} = -\frac{1}{4}F_{\mu\nu}F^{\mu\nu} = \tfrac{1}{2}(E^2-B^2)
\end{gather*}
where
\begin{gather*}
  \mat{F}^{\mu\nu}
  =
  \begin{pmatrix}
    0 & -\vect{E}\\
    \vect{E}^{T} & \mat{\vect{B}\times}
  \end{pmatrix}
  =
  \begin{pmatrix}
    0 & -E_x & -E_y & -E_z\\
    E_x & 0 & -B_z & B_y\\
    E_y & B_z & 0 & -B_x\\
    E_z & -B_y & B_x & 0
  \end{pmatrix},\\
  \mat{F}^{\mu}{}_{\nu}
  =
  \begin{pmatrix}
    0 & \vect{E}\\
    \vect{E}^{T} & -\mat{\vect{B}\times}
  \end{pmatrix},\qquad
  \mat{F}_{\mu\nu}
  =
  \begin{pmatrix}
    0 & \vect{E}\\
    -\vect{E}^{T} & \mat{\vect{B}\times}
  \end{pmatrix}.
\end{gather*}
The dual field-strength tensor is also important, and can be formed from the
4-dimensional [Levi-Civita symbol][]:
\begin{gather*}
  \varepsilon_{0123} = -\varepsilon^{0123} = 1, \qquad
  \tilde{F}_{\mu\nu} = \tfrac{1}{2}\varepsilon_{\mu\nu\rho\sigma}F^{\rho\sigma},\\
  \mat{\tilde{F}}_{\mu\nu} = 
  \begin{pmatrix}
    0 & -\vect{B}\\
    \vect{B}^{T} & -\mat{\vect{E}\times}
  \end{pmatrix}=
  \begin{pmatrix}
    0 & -B_x & -B_y & -B_z \\
    B_x & 0 & E_z & - E_y\\
    B_y & -E_z &  0 & E_x\\
    B_z & E_y & -E_x & 0
  \end{pmatrix},
  
\end{gather*}


Note that the dual tensor $\tilde{F}_{\mu\nu}$ is obtained from $F^{\mu\nu}$ by replacing
\begin{gather*}
  \vect{E} \mapsto \vect{B}, \qquad
  \vect{B} \mapsto -\vect{E}.
\end{gather*}
Contracting, we have the following, which gives us the expression for $\mathcal{L}$ above:
\begin{align*}
  F_{\mu\nu}F^{\mu\nu} &= 2(B^2 - E^2), &
  \tilde{F}_{\mu\nu}\tilde{F}^{\mu\nu} &= 2(E^2 - B^2), &
  \tilde{F}_{\mu\nu}F^{\mu\nu} &= 2\vect{B}\cdot\vect{E}.
\end{align*}
\begin{gather*}
  F^{\mu\rho}F_{\rho\nu} =
  2\begin{pmatrix}
    E^2 & \vect{B}\times\vect{E}\\
    (\vect{B}\times\vect{E})^T &
    \vect{B}\vect{B}^T - B^2\mat{1}
  \end{pmatrix}
\end{gather*}





We consider terms in the
Lagrangian one at a time, including their contribution to the momentum:
\begin{gather*}
  (\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu})
  = \dot{A}_{\nu}\dot{A}^{\nu} - \vect{\nabla}A_{\nu}\cdot\vect{\nabla}A^{\nu}
  = \dot{A}_0^2 - \dot{\vect{A}}^2 - (\vect{\nabla}A_0)^2 + (\nabla_{i}A_{j})(\nabla_{i}A_{j}),\\
  (\partial_{\mu}A_{\nu})(\partial^{\nu}A^{\mu})
  = \dot{A}^{\nu}\partial_{\nu}A_0 + \vect{\nabla}A^{\nu}\cdot\partial_{\nu}\vect{A} 
  = \dot{A}_0^2 + 2\dot{\vect{A}}\cdot\vect{\nabla}A_0 + (\nabla_{i}A_{j})(\nabla_{j}A_{i}),\\
  (\partial_{\mu}A^{\mu})^2
  = (\dot{A}_0 + \vect{\nabla}\cdot\vect{A})^2
  = \dot{A}_0^2 + 2\dot{A}_0\vect{\nabla}\cdot\vect{A} + (\vect{\nabla}\cdot\vect{A})^2.
\end{gather*}
The last two are equivalent under integration by parts twice.  Here are the
corresponding contributions to the momentum $\pi^{\mu}$:
\begin{gather*}
  \tfrac{1}{2}\pdiff{(\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu})}{\dot{A}_{\mu}}
  = \dot{A}^{\mu} = \begin{pmatrix}
    \dot{A}_0\\
    \dot{\vect{A}}
  \end{pmatrix},\\
  \tfrac{1}{2}\pdiff{(\partial_{\mu}A_{\nu})(\partial^{\nu}A^{\mu})}{\dot{A}_{\mu}}
  = \partial^{\mu}A_0
  = \begin{pmatrix}
    \dot{A}_0\\
    -\vect{\nabla}A_0
  \end{pmatrix},\\
  \tfrac{1}{2}\pdiff{(\partial_{\mu}A^{\mu})^2}{\dot{A}_{\mu}}
  = g^{0\mu}\partial_{\nu}A^{\nu}
  = \begin{pmatrix}
    \dot{A}_0 + \vect{\nabla}\cdot\vect{A}\\
    \vect{0}
  \end{pmatrix}.
\end{gather*}
The corresponding contributions to the Hamiltonian density are:
\begin{gather*}
  \dot{A}^{\mu}\dot{A}_{\mu} 
  - \tfrac{1}{2}(\partial_{\mu}A_{\nu})(\partial^{\mu}A^{\nu}) 
  = \tfrac{1}{2}\dot{A}_0^2 - \tfrac{1}{2}\dot{\vect{A}}^2
  + \tfrac{1}{2}[(\vect{\nabla}A_0)^2 - (\nabla_{i}A_{j})(\nabla_{i}A_{j})],\\
  \dot{A}^{\mu}\partial_{\mu}A_0 
  - \tfrac{1}{2}(\partial_{\mu}A_{\nu})(\partial^{\nu}A^{\mu})
  = \tfrac{1}{2}\dot{A}_0^2 - \tfrac{1}{2}(\nabla_{i}A_{j})(\nabla_{j}A_{i}),\\
  \dot{A}_{0}\partial_{\nu}A^{\nu}
  -\tfrac{1}{2}(\partial_{\mu}A^{\mu})^2
  = \tfrac{1}{2}\dot{A}_0^2 - \tfrac{1}{2}(\vect{\nabla}\cdot\vect{A})^2.
\end{gather*}

Thus, if we start with the Lagrangian density
\begin{gather*}
  \mathcal{L} = \tfrac{1}{2}\Bigl(
    -(\partial_{\mu} A_\nu)(\partial^{\mu} A^\nu)
    +
    (\partial_{\mu} A_\nu)(\partial^{\nu} A^\mu)
  \Bigr),
\end{gather*}
then our Hamiltonian density will have the form
\begin{gather*}
  \mathcal{H} = \tfrac{1}{2}\Bigl(
    \dot{\vect{A}}^2
    -(\vect{\nabla}A_0)^2 
    + 
    \underbrace{(\nabla_{i}A_{j})(\nabla_{i}A_{j})
    - (\nabla_{i}A_{j})(\nabla_{j}A_{i})}_{B^2} 
    \Bigr),\\
    = \tfrac{1}{2}\Bigl(
    E^2 + B^2 
    -2\vect{\nabla}A_0\cdot(\vect{\nabla}A_0 - \dot{\vect{A}})
    \Bigr),\\
    = \tfrac{1}{2}\Bigl(
    E^2 + B^2 
    +2\vect{\nabla}A_0\cdot\vect{E}
    \Bigr),\\
\end{gather*}
Expressing this in terms of the electric and magnetic fields, we have
\begin{gather*}
  \vect{E} = -\vect{\nabla}A_0 - \dot{\vect{A}}, \qquad
  \vect{B} = \vect{\nabla}\times \vect{A},\\
  E^2 = \dot{\vect{A}}^2 + (\vect{\nabla}A_0)^2 + 2\dot{\vect{A}}\cdot\vect{\nabla}A_0,\\
  B^2 = \underbrace{\varepsilon_{ija}\varepsilon_{kla}}_{\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk}}
        (\nabla_iA_j)(\nabla_k A_l)
      = (\nabla_iA_j)(\nabla_i A_j)
        -
        (\nabla_iA_j)(\nabla_j A_i)
\end{gather*}

As a check, the stress-energy tensor is
\begin{gather*}
  \theta^{\mu}_{\nu} = \pdiff{\mathcal{L}}{(\partial_{\mu}A_{\alpha})}\partial_{\nu}A_{\alpha} -
  \delta^{\mu}_{\nu}\mathcal{L}\\
  =
  -(\partial^{\mu} A^{\alpha})(\partial_{\nu}A_{\alpha})
  +
  (\partial^{\alpha} A^{\mu})(\partial_{\nu}A_{\alpha})
  +
  \frac{\delta^{\mu}_{\nu}}{2}
  \Bigl(
   (\partial_{\alpha} A_{\beta})(\partial^{\alpha} A^{\beta})
    -
    (\partial_{\alpha} A_\beta)(\partial^{\beta} A^{\alpha})
  \Bigr).
\end{gather*}
The Hamiltonian density is
\begin{gather*}
  \mathcal{H} = \theta^{0}_{0} 
  =
  -(\partial^{0} A^{\alpha})(\partial_{0}A_{\alpha})
  +
  (\partial^{\alpha} A^{0})(\partial_{0}A_{\alpha})
  +
  \frac{1}{2}
  \Bigl(
   (\partial_{\alpha} A_{\beta})(\partial^{\alpha} A^{\beta})
    -
    (\partial_{\alpha} A_\beta)(\partial^{\beta} A^{\alpha})
  \Bigr)\\
  =
  \tfrac{1}{2}
  \Bigl(
   \dot{\vect{A}}^2 - (\vect{\nabla}A_{0})^2 
   +(\nabla_{i} A^{j})(\nabla_{i} A^{j})
   -(\nabla_{i} A^{j})(\nabla_{j} A^{i})
  \Bigr).
\end{gather*}

According to Marggiore, we can add to this a total derivative to get an "improved"
energy-momentum tensor that is gauge invariant:
\begin{gather*}
  T^{\mu}_{\nu} = \theta^{\mu}_{\nu} + \partial_{\rho}(F^{\mu \rho}A_{\nu}).
\end{gather*}
This adds the following piece to $\mathcal{H}$:
\begin{gather*}
  \partial_{\rho}\Bigl((\partial_0A^{\rho} - \partial^{\rho}A^{0})A_{0}\Bigr)
  =
  (\dot{\vect{A}} + \vect{\nabla}A_{0})\cdot\vect{\nabla}A_{0}
  +
  (\partial_{\rho}\dot{A}^{\rho} - \partial^2A_0)A_{0}\\
  =
  (\dot{\vect{A}} + \vect{\nabla}A_{0})\cdot\vect{\nabla}A_{0}
  +
  A_{0}(\vect{\nabla}\cdot\dot{\vect{A}} + \nabla^2 A_0)
  =
  (\dot{\vect{A}} + \vect{\nabla}A_{0})\cdot\vect{\nabla}A_{0}
  -\dot{\vect{A}}\cdot\vect{\nabla}A_{0} - (\vect{\nabla}A_0)^2
\end{gather*}
This gives:
\begin{gather*}
  \mathcal{H} = T^{0}_{0}
  =
  \tfrac{1}{2}
  \Bigl(
   \dot{\vect{A}}^2 +
   (\vect{\nabla}A_{0})^2 
   + 2\dot{\vect{A}}\cdot\vect{\nabla}A_{0} 
   + 2(\partial_{\rho}\dot{A}^{\rho} - \partial^2A_0)A_{0}
   +(\nabla_{i} A^{j})(\nabla_{i} A^{j})
   -(\nabla_{i} A^{j})(\nabla_{j} A^{i})
  \Bigr).
\end{gather*}

::::
 

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
[adjoint representation][] (we follow the form in {cite:p}`Georgi:2019`):
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

```{code-cell} ipython3
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

```{code-cell} ipython3
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



[moments]: <https://en.wikipedia.org/wiki/Moment_(physics)>
[Levi-Civita symbol]: <https://en.wikipedia.org/wiki/Levi-Civita_symbol>
[structure constants]: <https://en.wikipedia.org/wiki/Structure_constants>
[adjoint representation]: <https://en.wikipedia.org/wiki/Adjoint_representation#Structure_constants>
[Jacobi identity]: <https://en.wikipedia.org/wiki/Jacobi_identity>
[Pauli matrices]: <https://en.wikipedia.org/wiki/Pauli_matrices>
[trivial representation]: <https://en.wikipedia.org/wiki/Trivial_representation>
[Lie algebra]: <https://en.wikipedia.org/wiki/Lie_algebra>
[angular momentum]: <https://en.wikipedia.org/wiki/Angular_momentum>
