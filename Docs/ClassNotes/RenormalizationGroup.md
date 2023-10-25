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

(sec:RG)=
Renormalization Group Techniques
================================

Here we present a practical introduction to the key concepts [renormalization group] (RG)
techniques in physics by way of a few examples, followed by a bit of philosophy.


## The Essence

:::{margin}
changing scales
:::
The essence of RG techniques is to understand how the structure of a physical theory
changes as we consider a system at **different scales**.  I.e. at a microscopic scale,
water concerns the dynamics of a bunch of particles, but if we scale out, we can
describe many behaviors using the equations of fluid dynamics, for example, as embodied
in the [Navier-Stokes equations].

:::{margin}
RG flow
:::
One possibility is that the theory looks very different at different scales.  Even the
appropriate degrees of freedom might be different -- water molecules in one case, fluid
density, velocity, and temperature in another.  How a theory changes as we change the
scale is called **RG flow**.

:::{margin}
running coupling
:::
If the structure of the theory remains the same, but the parameters values change, then
these are said to **run**.  In field-theories, the parameters "couple" different fields
together, and so are called "couplings": in this context the term **[running coupling]** is
often used.

:::{margin}
fixed points
:::
Sometimes, one finds regions of a theory where the couplings run very slowly, and the
same theory applies over a wide range of scales.  In the extreme case, we may be able to
change the scale while rescaling the parameters in such a way that the theory does not
change.  If we can find a point where this behavior is seen, then we have found a
**fixed point** of the RG flow.  At or near a fixed point, the theory will exhibit
[self-similarity].  This is common near second-order phase transitions, an in these
cases RG techniques really shine, providing accurate quantitative techniques for
calculating properties of the system such as [critical exponents].

:::{margin}
asymptotic freedom
:::
Taken to an extreme, one might demand that a fundamental theory be well defined at all
scales.  This was a common strategy in the [reductionist](https://en.wikipedia.org/wiki/Reductionism) approach towards a [theory of
everything (ToE)](https://en.wikipedia.org/wiki/Theory_of_everything), and an ideal
realized in [QCD] -- the theory of quarks and gluons that describes the strong
interactions.  [QCD] exhibits the a property called **[asymptotic freedom]**  whereby
the coupling constants get weaker as one goes to smaller and smaller length scales
(equivalently, to higher and higher energies).  Asymptotic freedom causes the theory to
become non-interacting at small distances (high energies), allowing one to consider the
theory at arbitrarily small scales.  Asymptotically free theories like [QCD] can in
principle be complete, providing a complete description of everything. 

:::{margin}
[Standard Model] of particle physics
:::
This reductionist ideal, however, is thwarted by [QED] which is not asymptotically
free.  Instead, the QED coupling constant (the fine-structure constant $\alpha$) becomes
larger as one moves to smaller scales.  The coupling constant diverges at the [Landau
pole], rendering the theory invalid.  Thus, QED -- one of the most accurate physical theory to
date -- is known to be incomplete, and the quest for [physics beyond the Standard Model]
continues, with billions of dollars being directed to smash particles together with
higher and higher energies to try to figure out what the heck is wrong.[^1]

[^1]: As of 2022, the current answer is: not much.  We still have very little evidence
      guiding us towards the new physics needed to complete the [Standard Model].


:::{margin}
relevant, irrelevant, and marginal terms
:::
The view we shall take here, as clearly described in {cite:p}`Huang:2013`, is that
simple physical theories generically arise from RG flow through "theory space" which
gets "stuck" near (but not at) fixed-points.  Near these fixed-points, most of the
terms in the theory becomes small (**irrelevant**) and the theory takes on a simple form,
with only a handful of **relevant** or **marginal** terms remaining.  Simplicity in
physics, thus emerges naturally and somewhat generically from RG flow: it does not
necessarily imply that there is a simple underlying theory.  

:::{margin}
quantum computing
:::
Conversely, this picture also suggests that strange behavior may exist if we finely tune
parameters or the scales.  For example, I conjecture that the advantages offered by
[quantum computing] and technologies will arise because of such unusual behaviours,
however, realizing these advantages will be difficult precisely because of the required
fine-tuning.

:::{margin}
applications
:::
Renormalization group techniques are not limited to particle physics.  They apply to
many different systems where the behaviour changes slowly from one scale to another.  In
these lectures, we will consider a few representative example, but keep in mind that many
applications exists to, percolation theory *(i.e. how forest fires spread, see
{cite:p}`Creswick:1992` for examples)*, traffic flow *(see [Kshitij
Jerath](http://selforganizing.systems/)s work, or his [iSciMath
seminar](https://discourse.iscimath.org/t/mar-17-renormalization-group-theory-the-systems-engineering-perspective/491)
-- request an account if you would like access and email [Michael](mailto:m.forbes+iscimath@wsu.edu))*,
fluid dynamics *(see [Jorge Noronha's WSU
Colloquium](https://www.youtube.com/watch?v=49AOXU06oJk) or [Pavel Kovtun's Aspen seminar](https://www.youtube.com/watch?v=XaF_22S4scs))*, statistical mechanics, condensed matter, and many other.

The key RG process involves two steps:
:::{margin}
semigroup
:::
1. **Coarse Graining:** The first aspect of an RG analysis is to consider a system at
   different scales. Often this involves some type of averaging, or "decimation" that
   implies loss of information (in which cased the RG should be considered as a
   [semigroup] with no exact inverse). 
2. **Rescaling:** The second step is to rescale the theory so that it "looks like" the
   original theory.  Most of the quantitative power of RG comes when the theory flows to
   a fixed point where the structure of the theory is almost invariant under an
   appropriate coarse-graining + rescaling.

## Random Walks

:::{margin}
When I last checked, it was available for only [$2,242.19 on
Amazon!](https://www.amazon.com/Introduction-Renormalization-Group-Methods-Physics/dp/0486793451).
Creswick has given us permission to distribute electronic copies, and we include one in
the project [Resources folder on CoCalc].  Ask if you need permission to ask.
:::
The first example we will consider ({ref}`sec:RG-RandomWalks`) is that of [random
walk]s.  Much of the content comes from an excellent, but generally hard-to-obtain book
{cite:p}`Creswick:1992`, but the salient points can be found in the following notes by
John McGreevy {cite:p}`McGreevy:2018`: 

* [Physics 217, The Renormalization Group, Fall 2018 (UCSD)](https://mcgreevy.physics.ucsd.edu/f18/)
* [McGreevy's Lecture Notes (PDF)]

:::{margin}
fixed points, coarse graining, and rescaling
:::
We start with a gaussian random walk, which provides an example of a fixed point.
Coarse graining amounts to taking $N$ steps, for which the probability distribution 
is again a gaussian (products of gaussians are gaussians), but with a different
variance.  The rescaling step brings the variance back to the starting value,
demonstrating exact self-similarity in this case.

:::{margin}
RG flow, universaity
:::
Moving away from gaussian distribution, we discuss RG flow by looking at the [cumulant]s
as parameters, and show that a wide range of different distributions flow towards a
gaussian distribution -- a manifestation of universality reminiscent of the [central
limit theorem].

## How to Renormalize the Schrödinger Equation

:::{margin}
effective field theory
:::
The second example we will consider in {ref}`sec:RG-SEQ` works through parts of 
{cite:p}`Lepage:1997`.  Here the central techniques of 
modern [effective field theory] (especially in the context of nuclear physics) are
demonstrated by considering the radial Schrödinger equation.  The analogies made in
{cite:p}`Lepage:1997` are accurate and deep: although understanding everything will
likely be challenging, it is highly worth working through the examples and trying to
understand as much as possible.

:::{margin}
regularization
:::
This example explicitly demonstrates the nature of many of the divergences seen in field
theory, and provides a constructive approach to "regularizing" these divergences to
obtain accurate physical results by developing and fitting an appropriate low-energy effective
theory.

Discussion points:

* RG Limit cycles and the $1/r^2$ potential.  (Start with $1/r$ and show that a
  dimensionful scale emerges.  $1/r^2$ has no scale.  What does this mean?  How does a
  scale for $E_0$ get set?


## Ising Model

Time-permitting, we will work though some examples related to the [Ising model].
Exactly solvable examples are available in 1D and 2D, allowing one to test the RG
approach for studying critical phenomena and phase transitions.  This discussion follows
[McGreevy's Lecture Notes (PDF)] {cite:p}`McGreevy:2018` and the discussion in
{cite:p}`Creswick:1992`.

Discussion points:

* The coarse-graining procedure for the Ising model s not what one might naturally
  expect.  This is so that the structure of the theory remains fixed.  Can one do an RG
  analysis with the natural coarse graining?
* Good example of mean field theory.

## Philosophy

Finally, we will end with a discussion of renormalization group techniques based on
{cite:p}`Huang:2013`, putting forth the idea that our relatively simple physical theories result
from RF flow through the "space of theories" which linger near fixed-points where
physics is well approximated by a simple theory over a large range of scales.

Discussion points:

* Asymptotic freedom and charge screening: charged gluons is key.
* Why quantum gravity is not renormalizable?  At each level in perturbation theory one
  needs to 
  * [Shomer, A: "A pedagogical explanation for the non-renormalizability of
    gravity"](https://arxiv.org/abs/0709.3555v2)
  * ['t Hooft: "Dimensional Reduction in Quantum Gravity"](https://arxiv.org/abs/gr-qc/9310026)
* Renormalizability:
  * For a good formal discussion, see Chapter 4 {cite:p}`Coleman:1988` which is based on
    Hepp's theorem (stated without proof).

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

[universality class]: <https://en.wikipedia.org/wiki/Universality_class>
[critical exponents]: <https://www.google.com/search?client=firefox-b-1-d&q=critical+exponents>
[running coupling]: <https://en.wikipedia.org/wiki/Coupling_constant#Running_coupling>
[physics beyond the Standard Model]: <https://en.wikipedia.org/wiki/Physics_beyond_the_Standard_Model>
[Landau pole]: <https://en.wikipedia.org/wiki/Landau_pole>
[QED]: <https://en.wikipedia.org/wiki/Quantum_electrodynamics>
[QCD]: <https://en.wikipedia.org/wiki/Quantum_chromodynamics>
[asymptotic freedom]: <https://en.wikipedia.org/wiki/Asymptotic_freedom>
[Standard Model]: <https://en.wikipedia.org/wiki/Standard_Model>
[quantum computing]: <https://en.wikipedia.org/wiki/Quantum_computing>
[renormalization group]: <https://en.wikipedia.org/wiki/Renormalization_group>
[random walk]: <https://en.wikipedia.org/wiki/Random_walk>
[semigroup]: <https://en.wikipedia.org/wiki/Semigroup>
[Navier-Stokes equations]: <https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations>
[self-similarity]: <https://en.wikipedia.org/wiki/Self-similarity>
[cumulant]: <https://en.wikipedia.org/wiki/Cumulant>
[central limit theorem]: <https://en.wikipedia.org/wiki/Central_limit_theorem>
[effective field theory]: <https://en.wikipedia.org/wiki/Effective_field_theory>
[Ising model]: <https://en.wikipedia.org/wiki/Ising_model>
[Resources folder on CoCalc]: <https://cocalc.com/projects/629fc51f-1039-41e6-a64f-b4301f3ee35b/files/Resources>

[McGreevy's Lecture Notes (PDF)]: <https://mcgreevy.physics.ucsd.edu/f18/2018F-217-lectures.pdf>
