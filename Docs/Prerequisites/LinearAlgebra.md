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
:tags: [remove-cell]

import mmf_setup;mmf_setup.nbinit()
import logging;logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
```

(sec:linear-algebra)=
Linear Algebra
==============

Although this course does not assume you have a strong physics background, a solid
grounding in linear algebra is essential.  These notes cover the important material you
will be expected to understand deeply before participating in the course.  If anything
is unfamiliar, please refresh your understanding prior to the course, or ask one of the
instructors for a review.

Please see {ref}`sec:linear_algebra_resources` for additional resources, including and
excellent set of videos.

## Terminology

We will freely use the following terminology throughout these notes and in the course.
You will be expected to understand what these terms mean, and the implications.  For
example, you should know that if the **linear operator** $\op{H}$ is **Hermitian**, then its
**eigenvalues** are real, and its **eigenvectors** can be chosen to form a **complete**,
**orthonormal basis**.  In equations:

\begin{gather*}
  \op{H}\ket{n} = \ket{n}\lambda_n, \qquad
  \lambda_n \in \mathbb{R}, \qquad
  \braket{m|n} = \delta_{mn}.
\end{gather*}

Some of these terms will be described here, but you may need to refer elsewhere to
complete your understanding.  We will restrict ourselves to finite-dimensional vector
spaces in this course, so some of the terminology will be redundant (i.e. for finite
dimensional spaces, Hermitian and self-adjoint mean the same thing)

* Set.
* Vector, column vector.
* Co-vector, row vector.
* Dyad.
* Vector space (complex).
* Subspace.
* Rank.
* Linear and anti-linear operator.
* Inner product $\braket{a|b}$.
* Basis (pl. bases).
* Components.
* Complete.
* Orthogonal.
* Orthonormal.
* Transpose $\mat{A}^T$ and Hermitian conjugate $\mat{A}^\dagger = (\mat{A}^{T})^*$.
* Symmetric and anti-symmetric, Hermitian, self-adjoint.
* Orthogonal matrix, unitary matrix.
* Singular value decomposition (SVD): $\mat{A} = \mat{U}\mat{D}\mat{V}^\dagger$.
* QR decomposition $\mat{A} = \mat{Q}\mat{R}$. I.e. Gram-Schmidt orthogonalization.
* LR decomposition $\mat{A} = \mat{L}\mat{U}$. I.e. Gaussian elimination.
* Trace and determinant.
* Eigenvectors and eigenvalues.
* Tensor product.
* Cartesian product.
* Partial trace.

### Notation

We will use the following "bra-ket" or "braket" notation due to Dirac, which is
ubiquitous and very useful in quantum mechanics.

```{margin}
Pronounce $\ket{\psi}$ "ket psi".
If you must think of $\ket{\psi}$ as a "list of numbers, think of $\ket{\psi}$ as
**column vector**: 

\begin{gather*}
  \ket{\psi} \= \begin{pmatrix}
    1\\
    3 - 4\I
  \end{pmatrix}.
\end{gather*}
```
* $\ket{\psi}$: An **abstract vector** or simply, a **vector**.  This represents a
  physical object: think of a physical arrow in space, in contradistinction to a "list
  of numbers".  Although one can treat a "list of numbers" as an abstract **column
  vector**, this is better thought of as the components of an abstract vector in a
  particular **basis**. This distinction causes quite a bit of confusion, and is the
  focus of the section {ref}`sec:components` below.
  
```{margin}
Pronounce $\bra{\psi}$ as "bra psi", and think of it as s **row vector**, i.e.:
\begin{gather*}
  \bra{\psi} \= \begin{pmatrix}
    1 & 3 + 4\I
  \end{pmatrix}.
\end{gather*}
Don't forget to take the complex conjugate!
```
* $\bra{\psi}$: An **abstract co-vector** or simply, a **co-vector**.  This is a linear
  operator that takes vectors into numbers, and is the **Hermitian conjugate** of the
  vector $\ket{\psi}$:
  
  \begin{gather*}
    \bra{\psi} = \ket{\psi}^\dagger.
  \end{gather*}

  Although one can think of $\bra{\psi}$ as a row vector, it is better to think
  $\bra{f}$ as a linear operator that takes vectors into numbers.
  
```{margin}
Pronounce $\braket{f|g}$ "the braket of $f$ and $g$":
\begin{gather*}
  \braket{\psi|\psi} \equiv \bra{\psi}\!\ket{\psi} \\
  =
  \begin{pmatrix}
    1 & 3 + 4\I
  \end{pmatrix}
  \begin{pmatrix}
    1\\
    3 - 4\I
  \end{pmatrix}\\
  = 26.
\end{gather*}
```
* $\braket{f|g}$: **Inner product**.  Much of the power of Dirac's notation comes from
  the expression of the inner product as an operator product of bras and kets.  This
  inner product can be thought of as the dot product of a row vector and column vector.
* $\norm{\ket{\psi}} = \sqrt{\braket{\psi|\psi}}$: **Norm** or **$L_2$ norm**.  We will
  only make use of the $L_2$ norm in this course:
  
  \begin{gather*}
    \norm{\ket{\psi}} = \sqrt{\sum_{n} \abs{\psi_n}^2}.
  \end{gather*}

* $\{\ket{n}\}$: **Orthonormal basis**.  More carefully, we should specify the number of
  basis vectors such as:
  
  \begin{gather*}
    \{\ket{n} \;|\; n=0, 1, \dots, N-1\} = 
    \{\ket{0}, \ket{1}, \dots, \ket{N-1}\}.
  \end{gather*}

  We will assume that all bases have been properly orthogonalized an normalized so that
  
  \begin{gather*}
    \braket{m|n} = \delta_{nm} = \begin{cases}
      1 & m = n,\\
      0 & \text{otherwise},
    \end{cases}
  \end{gather*}
  ```{margin}
  The Kronecker-delta is the index form of the identity matrix:
  \begin{gather*}
    \delta_{mn} = [\Bbb{1}]_{mn}
  \end{gather*}
  ```
  where $\delta_{mn}$ is the [Kronecker-delta]. This assumption greatly simplifies the
  notation of things like the components of a vector (see below).  When working with
  concrete vectors, then you may need to first orthonormalize your bases using the
  Gram-Schmidt procedure, or (equivalently) the QR decomposition.
  
* $\psi_n = \braket{n|\psi}$: The **components** of the vector $\ket{\psi}$ in the
  **orthonormal basis**
  $\{\ket{n}\}$.  Here $\bigl\{\psi_n \; |\;  n\in\{0, 1, \dots, N-1\}\bigr\}$ is the
  "list of numbers" version of a vector you might be used to from elementary linear
  algebra courses.

<!-- See https://docs.withorbit.com/ -->
<orbit-reviewarea>
  <orbit-prompt
    question="Is the ket $\ket{\psi}$ best thought of as a row-vector or a column vector?"
    answer="A column vector:
    $$\ket{\psi} = \begin{pmatrix}5 \\ 2+3i \end{pmatrix}$$"
  ></orbit-prompt>
  <orbit-prompt
    question="What is the bra $\bra{\psi}$ corresponding to the following ket?
    $$\ket{\psi} = \begin{pmatrix}5 \\ 2+3i \end{pmatrix}$$"
    answer="$$\bra{\psi} = \begin{pmatrix}5 & 2-3i \end{pmatrix}$$"
  ></orbit-prompt>
  <orbit-prompt
    cloze="Two vectors for which the braket $\langle a|b\rangle = 0$ are said to be {orthogonal}."
  ></orbit-prompt>
  <orbit-prompt
    cloze="If, in addition, $\langle a|a\rangle = 1$ and
    $\langle b|b\rangle = 1$, they are said to be {orthonormal}."
  ></orbit-prompt>
</orbit-reviewarea>

(sec:components)=
## Vectors and Components

:::{margin}
  This is (almost) the last time we will use the traditional vector notation $\vect{v}$
  and $\uvect{i}$ etc.  From now on, we will favour their **ket** equivalents:
  \begin{gather*}
    \vect{v} \equiv \ket{v}, \quad
    \uvect{i} \equiv \ket{i},
  \end{gather*}
  etc.  The rare exception will be when the vector represents an index on something
  else, like a matrix.  Thus, we may occasionally use $\vect{\mat{\sigma}} =
  (\mat{\sigma}^x, \mat{\sigma}^y, \mat{\sigma}^z)$ as the "vector" of [Pauli matrices]
  so we can form expressions like:
  \begin{gather*}
    \vect{B}\cdot\vect{\mat{\sigma}} = 
    B_x\mat{\sigma}^x + B_y\mat{\sigma}^y + B_z\mat{\sigma}^z,
  \end{gather*}
  but this will be rare.
:::
Students sometimes have difficulty working with quantum mechanics because of confusion
between a **vector** and the **components of a vector**.  To be clear, when we write a
vector $\vect{v} \equiv \ket{v}$ as:

\begin{gather*}
  \ket{v} \=
  \begin{pmatrix}
    1\\
    -2\\
    3
  \end{pmatrix}
\end{gather*}

we really have in mind a **standard orthonormal basis** $\{\uvect{i}=\ket{i}, \uvect{j}=\ket{j},
\uvect{k}=\ket{k}\}$

\begin{gather*}
  \ket{i} \= \begin{pmatrix}
    1\\
    0\\
    0
  \end{pmatrix}, \qquad
  \ket{j} \= \begin{pmatrix}
    0\\
    1\\
    0
  \end{pmatrix}, \qquad
  \ket{k} \= \begin{pmatrix}
    0\\
    0\\
    1
  \end{pmatrix}.
\end{gather*}

such that

\begin{gather*}
  \ket{v} = 1\ket{i} - 2\ket{j} + 3\ket{k}.
\end{gather*}

```{note}
In my work, I try to write this as

\begin{gather*}
  \ket{v} = \ket{i}1 - \ket{j}2 + \ket{k}3.
\end{gather*}

Although not common, this helps me keep track of things.  The reason is that I think of
$\ket{i}$ etc. as column vectors -- matrices with shape $N\times 1$ -- and I think of
the coefficients (numbers) as $1\times 1$ matrices.  Thus, I can think of $\ket{i}1$ as
matrix multiplication between a $N\times 1$ matrix and a $1 \times 1$ matrix, which 
makes sense in terms of the dimensions, and gives an $N \times 1$ matrix as a result.

Similarly, I think of the long bar on the left side of the ket as corresponding to the
first index of the matrix that is "long" in the sense that it takes $N$ values, while
the point on the right corresponds second index which takes only $1$ value.

Thus, when I look at $\braket{f|g}$, the points on both sides tell me that this is a
$1\times 1$ matrix, or a number, while the dyad $\ket{f}\bra{g}$ has long bars on both
sides, and so is clearly an $N \times N$ matrix.

Following this idea, I write the eigenvalue problem as:

\begin{gather*}
  \underbrace{\op{H}}_{N\times N}\;\underbrace{\ket{\psi}}_{N\times 1} =
  \underbrace{\ket{\psi}}_{N\times 1}\;\underbrace{\lambda}_{1\times 1}.
\end{gather*}

I think of the matrix $\mat{H}$ moving "through" the ket from left to right, and turning
into the eigenvalue $\lambda$.

Thus, when I diagonalize a matrix $\mat{H} = \mat{U}\mat{D}\mat{U}^\dagger$, I remember
that this is really just a collection of eigenvalue problems where the columns of
$\mat{U}$ are the eigenvectors:

\begin{gather*}
  \mat{H}\mat{U} = \mat{U}\mat{D},\\
  \mat{H}
  \overbrace{\begin{pmatrix}
    &\\
    \ket{u_0} & 
    \ket{u_1} & 
    \dots\\
    &
  \end{pmatrix}}^{\mat{U}}
  =\\
  \begin{pmatrix}
    &\\
    \mat{H}\ket{u_0} & 
    \mat{H}\ket{u_1} & 
    \dots\\
    &
  \end{pmatrix}
  =
  \begin{pmatrix}
    &\\
    \ket{u_0}\lambda_0 & 
    \ket{u_1}\lambda_1 & 
    \dots\\
    &
  \end{pmatrix}=\\
  = 
  \underbrace{\begin{pmatrix}
    &\\
    \ket{u_0} & 
    \ket{u_1} & 
    \dots\\
    &
  \end{pmatrix}}_{\mat{U}}
  \underbrace{\begin{pmatrix}
    \lambda_0\\
    & \lambda_1\\
    & & \ddots
  \end{pmatrix}}_{\mat{D}}
\end{gather*}
```

Mathematicians will insist that vectors $\ket{v}$ be treated as abstract objects, and
that the "list of numbers" only makes sense after one has specified a basis.  Hence,
expressions like the one above should be regarded with suspicion (hence the $\=$
notation).

:::{margin}
For example, in quantum mechanics, the relationship between an abstract quantum state
$\ket{\psi}$ and the **wavefunction** $\psi(x)$ is exactly such an expression:

\begin{gather*}
  \psi(x) = \braket{x|\psi}.
\end{gather*}

The "numbers" represented by the wavefunction $\psi(x)$ are just the coefficients of
$\ket{\psi}$ in the position basis $\{\ket{x}\}$.
:::
It is worthwhile developing a nose for such suspect relations, and keeping
this in mind can help clarify what is going on:

\begin{gather*}
  \ket{v} \=
  \begin{pmatrix}
    1\\
    -2\\
    3
  \end{pmatrix}
  =
  \begin{pmatrix}
    \braket{i|v}\\
    \braket{j|v}\\
    \braket{k|v}
  \end{pmatrix}.
\end{gather*}

For the purposes of this course, however, we shall always **assume the existence of a
standard orthonormal basis**, for example, with the $z$-axis pointing "up", the $x$-axis pointing
"to the right" along, and the $y$-axis pointing "into the page" (or "into the board" in
the classroom).

Thus, we shall freely write

\begin{gather*}
  \ket{v} =
  \begin{pmatrix}
    1\\
    -2\\
    3
  \end{pmatrix}
\end{gather*}

without qualification $\=$, with the knowledge that this description depends on our
implicit choice of the **standard orthonormal basis**.

:::{admonition} Discussion: Abstract vectors $\ket{\psi}$ vs. lists of numbers -- which is more fundamental?
:class: dropdown

Physicists often regard the quantum state $\ket{\psi}$ as an abstract quantity
representing physical reality, with the list of numbers having secondary meaning
dependent on ones (arbitrary) choice of basis.

This position, however, is not completely tenable.  As both Heisenberg and Einstein
point out **Fred: can you provide the references?**, what really matters **are the
numbers** which result from **measurement**.  The physical theory simply provides some
mathematical framework which describes how these numbers relate to each other in
different circumstances.  In this view, it is the numbers themselves that play the
fundamental role, and hence, it is reasonable to consider the numbers themselves as the
representation of the physics, with the choice of basis explicitly defined by the nature
of the experiment.

In quantum mechanics, measurement is subtle, and we will learn that the components of a
wavefunction cannot actually be measured, nevertheless, there is a reasonable argument
that the abstract vector $\ket{\psi}$ is no more fundamental that the components in a
particular well-defined basis.  We will use this argument to favour the simplified
notation.
:::

<!-- See https://docs.withorbit.com/ -->
<orbit-reviewarea>
  <orbit-prompt
    question="What assumption are we making implicitly when we write 
    $$\ket{v} = \begin{pmatrix} 1\\ -2\\ 3\end{pmatrix}?$$"
    answer="A choice of a standard orthonormal basis."
  ></orbit-prompt>
</orbit-reviewarea>

## Special Matrices

In quantum mechanics, two general classes of matrices play a central role: **hermitian
matrices** and **unitary matrices**.  These are the complex generalizations of
**symmetric** and **orthogonal** matrices respectively.

### Hermitian Matrices

**Hermitian** matrices are square matrices that are the same as their **hermitian
conjugate**:

\begin{gather*}
  \mat{H} = \mat{H}^\dagger = (\mat{H}^{T})^*.
\end{gather*}

E.g., for a 3×3 matrix:

\begin{gather*}
  \mat{H} = \begin{pmatrix}
    a & b & c\\
    d & e & f\\
    g & h & k
  \end{pmatrix}, \qquad
  \mat{H}^\dagger = \begin{pmatrix}
    a^* & d^* & g^*\\
    b^* & e^* & h^*\\
    c^* & f^* & k^*
  \end{pmatrix},  
\end{gather*}

thus, we can write the most general 3×3 hermitian matrix as:

\begin{gather*}
  \mat{H} = \begin{pmatrix}
    a & b+c\I & d+e\I\\
    b-c\I & f & g+h\I\\
    d-e\I & g-h\I & k
  \end{pmatrix},
\end{gather*}

where $a, b, \dots, k$ are real.  In particular, note that the diagonals are real.

For quantum mechanics, there are two key properties of hermitian matrices:

:::{margin}
If the eigenvalues are distinct $\lambda_1 \neq \lambda_2 \neq \cdots$, then one can
simply say that the eigenvectors are orthogonal. The subtlety is when there are
**degenerate eigenvalues** $\lambda_1 = \lambda_2$.  In this case, one must be careful
to chose the eigenvectors $\ket{1}$ and $\ket{2}$ to be orthogonal, since any linear
combination of corresponding eigenvectors is also an eigenvector, e.g. $\ket{1} + \ket{2}$.
:::
1. Their eigenvalues are real: $\lambda_i \in \mathbb{R}$.
2. Their eigenvectors can be chosen to form a complete orthonormal basis.

:::{admonition} Do it!  Prove these statements.
:class: dropdown

First suppose that $\ket{1}$ is an eigenvector:

\begin{gather*}
  \mat{H}\ket{1} = \ket{1}\lambda_1, \qquad
  (\mat{H}\ket{1})^\dagger = \bra{1}\mat{H}^\dagger = \lambda_1^*\bra{1}.
\end{gather*}

Hence, if $\mat{H} = \mat{H}^{\dagger}$ is hermitian, then

\begin{gather*}
  \braket{1|\mat{H}|1} = \braket{1|1}\lambda_1 = \braket{1|\mat{H}^\dagger|1} = \braket{1|1}\lambda_1^*.
\end{gather*}

Thus $\lambda_1 = \lambda_1^*$ is real.  This holds for any eigenvalue of $\mat{H}$

\begin{gather*}
  (\mat{H}\ket{2})^\dagger = \bra{2}\mat{H}^\dagger = \lambda_2^*\bra{2}.
\end{gather*}

Now we present a simplified proof in the case that the eigenvectors are distinct.  Suppose that

\begin{gather*}
  \mat{H}\ket{1} = \ket{1}\lambda_1, \qquad
  \mat{H}\ket{2} = \ket{2}\lambda_2, \qquad
  \lambda_1 \neq \lambda_2.
\end{gather*}

Then, since $\lambda_2^* = \lambda_2$ is real:

\begin{gather*}
  \braket{2|\mat{H}|1} = \braket{2|1}\lambda_1
  =\braket{2|\mat{H}^\dagger|1} \\
  = \braket{1|\mat{H}|2}^* 
  = \braket{1|2}^*\lambda_2^*
  = \braket{2|1}\lambda_2.
\end{gather*}

Since $\lambda_1 \neq \lambda_2$, the only way for this to be true is for $\braket{1|2}
= 0$.  I.e. the eigenvectors $\ket{1}$ and $\ket{2}$ are orthogonal.
:::

In quantum mechanics, physical quantities are associated with matrices as follows:

> Physical systems are represented by normalized state vectors $\ket{\psi}$ with
> $\braket{\psi|\psi} = 1$.  Physical quantities are represented by hermitian matrices,
> e.g. $\mat{H} = \mat{H}^\dagger$.  The result of measuring $\mat{H}$ for a state
> $\ket{\psi}$ is one of the eigenvalues $\lambda_n$ of $\mat{H}$ with probability $p_n
> = \braket{n|\psi}$.  After obtaining a measurement of $\lambda_n$, the physical system
> will be left in the corresponding eigenstate $\ket{n}$.

By restricting physical quantities to be associated with hermitian matrices, we ensure
that measured values are real, and that all events are properly accounted for with unit
probability $1 = p_1 + p_2 + \cdots$ of measuring at least one of the outcomes.  *(Note:
there are extensions of quantum mechanics that consider non-hermitian matrices that have
become quite popular recently.  This could form the basis for an interesting project if
someone is interested.  These generally represent "open" systems where particles can
leave, so the total probability of measuring one of the outcomes may be less than 1 for example.)*

### Positive Operator

A [positive operator][] (or positive semi-definite matrix) is an operator (matrix)
$\mat{A}$ for which

\begin{gather*}
  \braket{\psi|\mat{A}|\psi} \geq 0  
\end{gather*}

is real and non-negative for all vectors $\ket{psi}$ in the space.  If the vector space
is complex, then $\mat{A} = \mat{A}^\dagger$ is hermitian.

:::{admonition} Do it!  Prove this statement.
:class: dropdown

This is Exercise (2.24) in {cite:p}`Nielsen:2010`.  As a start, consider the following matrix:

\begin{gather*}
  \mat{A} = \begin{pmatrix}
    1 & 1\\
    0 & 1
  \end{pmatrix}.
\end{gather*}

1. Show that the two eigenvalues are both $1$.
2. Show that $\mat{A}$ is positive over all real vectors by considering

   \begin{gather*}
     \ket{\psi} = \begin{pmatrix}
       \cos(\theta)\\
       \sin(\theta)
     \end{pmatrix}
   \end{gather*}

3. Find a complex vector $\ket{\psi}$ where the expectation value
   $\braket{\psi|\mat{A}|\psi}$ is not real.
4. Follow the hint in Exercise (2.24) that every matrix $\mat{A} = \mat{B} + \I\mat{C}$
   can be expressed in this way where $\mat{B}=\mat{B}^\dagger$ and
   $\mat{C}=\mat{C}^\dagger$ are hermitian.
:::

### Unitary Matrices

**Unitary** matrices are square matrices whose hermitian conjugate is their inverse:

\begin{gather*}
  \mat{U}^\dagger = \mat{U}^{-1}, \qquad
  \mat{U}^\dagger \mat{U} = \mat{1}.
\end{gather*}

The most important properties are:

1. The columns of a unitary matrix form a complete orthonormal basis.
2. The eigenvalues of a unitary matrix are phases $e^{\I \theta}$.
3. The eigenvectors of a unitary matrix can be chosen to form a complete orthonormal
basis.
4. $\mat{U}\mat{U}^\dagger = \mat{1}$.

:::{admonition} Do it!  Prove these statements.
:class: dropdown

The first property is easy to show by partitioning the matrix.  Let $\ket{1}$ be the
first column, $\ket{2}$ the second, etc.  Then $\bra{1}$ forms the first row of
$\mat{U}^\dagger$, $\bra{2}$ the second row, etc.:

\begin{gather*}
  \mat{U} = \begin{pmatrix}
    \ket{1} & \ket{2} & \cdots & \ket{N}
  \end{pmatrix}, \qquad
  \mat{U}^\dagger = \begin{pmatrix}
    \bra{1} \\ \bra{2} \\ \vdots \\
    \bra{N}
  \end{pmatrix}.
\end{gather*}

Now simply compute the product:

\begin{gather*}
  \mat{U}^\dagger \mat{U} = 
  \begin{pmatrix}
    \braket{1|1} & \braket{1|2} & \cdots & \braket{1|N}\\
    \braket{2|1} & \braket{2|2} & \cdots & \braket{2|N}\\
    \vdots & & \ddots & \vdots\\
    \braket{N|1} & \braket{N|2} & \cdots & \braket{N|N}.
  \end{pmatrix}
  = 
  \mat{1}.
\end{gather*}

Hence $\braket{n|m} = \delta_{m,n}$ and the columns $\ket{m}$ form a complete
orthonormal basis.

Now assume $\mat{U}\ket{1} = \ket{1}\lambda$.  Then

\begin{gather*}
  \abs{\lambda}^2\braket{1|1} = 
  \overbrace{(\bra{1}\mat{U}^\dagger)}^{\lambda^*\bra{1}}\overbrace{\mat{U}\ket{1}}^{\ket{1}\lambda}
  =\\
  = \braket{1|\mat{U}^\dagger\mat{U}|1}
  = \braket{1|\mat{1}|1} = \braket{1|1} = 1.
\end{gather*}

Hence $\abs{\lambda} = 1$ and so $\lambda = e^{\I\theta}$ for some real angle $\theta$.

The other properties are easier to prove after we relate hermitian and unitary matrices.
:::

Unitary matrices play a central role in quantum mechanics because they preserve the norm
of a state $\ket{\psi}$, which is associated with conservation of probability.  I.e. if
$\ket{\phi} = \mat{U}\ket{\psi}$ where $\mat{U}$ is unitary, then

\begin{gather*}
    \braket{\phi|\phi} =
    \braket{\psi|\mat{U}^\dagger\mat{U}\psi} = 
    \braket{\psi|\mat{1}\psi} = 
    \braket{\psi|\psi}.
\end{gather*}

The **propagator** in quantum mechanics -- the linear operator that takes one state into
another state in the future -- is a unitary matrix $\mat{U}$.  In general, we will
represent the action of a quantum computer (without measurement) as a large unitary
matrix $\mat{U}$.

(sec:Projections)=
### Projection Matrices

[Projection][] matrices also play an important role, providing a way of decomposing vectors
into their components.  Their key property is that they are [idempotent][]:

\begin{gather*}
  \mat{P}^2 = \mat{P}.
\end{gather*}

I.e. once you project a vector into a subspace, projecting it again does not change anything.

For example, suppose


\begin{gather*}
  \newcommand{\proj}[2][]{\mat{P}^{#1}_{\ket{#2}}}
  \ket{\psi} = \ket{a}a + \ket{b}b
\end{gather*}


where $\ket{a}$ and $\ket{b}$ are linearly independent unit vectors.  Projection
matrices $\proj{a}$ and $\proj{b}$ that act as

\begin{gather*}
  \proj{a}\ket{\psi} = \ket{a}a, \qquad
  \proj{a}\ket{\psi} = \ket{b}b,
\end{gather*}

are useful:

\begin{gather*}
  \proj{a} = \frac{\ket{a}\!\bra{a}\proj[\perp]{b}}{\braket{a|\proj[\perp]{b}|a}}, \qquad
  \proj{b} = \frac{\ket{b}\!\bra{b}\proj[\perp]{a}}{\braket{b|\proj[\perp]{a}|b}},
\end{gather*}

where $\proj[\perp]{a}$ is the projector into the subspace orthogonal to $\ket{a}$ etc.

:::{admonition} Do it!  Derive projectors $\proj{a}$ and $\proj{b}$.
:class: dropdown

First, lets take $\ket{a}$ and $\ket{b}$ to be unit vectors $\braket{a|a} = \braket{b|b}
= 1$.  If this is not the case, then redefine and the coefficients as

\begin{gather*}
  \ket{a} \rightarrow \frac{\ket{a}}{\sqrt{\braket{a|a}}}, \qquad
  a \rightarrow a\sqrt{\braket{a|a}},
\end{gather*}

etc.

The key properties that we want are:

\begin{gather*}
  \proj{a}\ket{a} = \ket{a}, \qquad
  \proj{b}\ket{b} = \ket{b}, \\
  \proj{a}\ket{b} = \proj{b}\ket{a} = 0.
\end{gather*}

The second property can be ensured by using the projectors $\proj[\perp]{a}$ and
$\proj[\perp]{b}$ into 
the subspaces that are orthogonal to $\ket{b}$ and $\ket{a}$ respectively.  Similarly,
we need the "output" of the operators to be along the appropriate direction.  Thus, we
can write

\begin{gather*}
  \proj[\perp]{a} = \mat{1} - \ket{a}\!\bra{a}, \qquad
  \proj[\perp]{b} = \mat{1} - \ket{b}\!\bra{b},\\
  \proj{a} = \ket{a}\bra{A}\proj[\perp]{b},\qquad
  \proj{b} = \ket{b}\bra{B}\proj[\perp]{a},
\end{gather*}

where we have yet to specify $\bra{A}$ and $\bra{B}$.

To have $\proj{a}\ket{a} = \ket{a}\braket{A|\proj[\perp]{b}|a} = \ket{a}$, we need

\begin{gather*}
  1 = \braket{A|\proj[\perp]{b}|a} = \braket{B|\proj[\perp]{a}|b}.
\end{gather*}

This is enough to make the projectors idempotent, i.e.:

\begin{gather*}
  \proj[2]{a} = \ket{a}\braket{A|\proj[\perp]{b}|a}\bra{A}\proj[\perp]{b}
              = \ket{a}\bra{A}\proj[\perp]{b} = \proj{a}.
\end{gather*}

We can thus simply take (using the fact that $\proj[\perp]{b}$ is idempotent and hermitian)

\begin{gather*}
  \ket{A} = \frac{\proj[\perp]{b}\ket{a}}{\braket{a|\proj[\perp]{b}|a}}, \qquad
  \ket{B} = \frac{\proj[\perp]{a}\ket{b}}{\braket{b|\proj[\perp]{a}|b}}.
\end{gather*}

In higher dimensional spaces, one can add to $\ket{A}$ and $\ket{B}$ anything in the
subspace orthogonal to $\span\{\ket{a}, \ket{b}\}$.

Putting this all together, we have explicitly

\begin{gather*}
  \proj{a} = \frac{\ket{a}\!\bra{a}\proj[\perp]{b}}{\braket{a|\proj[\perp]{b}|a}}
           = \frac{\ket{a}\!\bra{a}\bigl(\bra{a} - \braket{a|b}\bra{b}\bigr)}
                  {1 - \abs{\braket{a|b}}^2}.
\end{gather*}

The first form generalizes to any number of projectors by generalizing $\proj[\perp]{b}$
to be the projection into the orthogonal subspace to all other vectors.

Note that if $\braket{a|b} = 0$, then this reduces to the familiar form for orthogonal
vectors $\proj{a} = \ket{a}\!\bra{a}$.
:::

We usually work with orthogonal vectors $\ket{a|b} = 0$ or orthonormal vectors
$\braket{a|a} = 1$, which greatly simplifies the formula:

\begin{gather*}
  \proj{a} = \frac{\ket{a}\!\bra{a}}{\braket{a|a}} = \ket{a}\!\bra{a} \text{ if }
  \braket{a|a} = 1.
\end{gather*}

::::{note}

If the vectors are not orthogonal, $\braket{a|b} \neq 0$, then one must be careful
because several properties enjoyed by orthogonal projects do not hold:

1. If the vectors are not orthogonal, then the projectors depend on the full set of
   vectors.  I.e. our formula for $\proj{a}$ depends on $\ket{b}$.  If the vectors are
   orthogonal, then this dependence vanishes.
2. If the vectors are not orthogonal, then the projectors are not hermitian, and
   therefore **not** positive operators on complex vector spaces.

For these reasons, it is usually preferably to change your metric so that your vectors
become orthonormal.

:::{admonition} To do.
Demonstrate how these projectors can be simply constructed using a different metric.
:::

::::

```{code-cell}
from qiskit.visualization import array_to_latex

def normalize(psi):
    return np.divide(psi, np.linalg.norm(psi))

a = np.array([1, 0])
b = np.array([1, 1])

psi = 1.2*a - 3.4*b
bhat = normalize(b)

I = np.eye(2)

P_perp_a = I - np.outer(a, a)
P_perp_b = I - np.outer(bhat, bhat)
P_a = (np.outer(a, a) @ P_perp_b) / (a @ P_perp_b @ a)
P_b = (np.outer(b, b) @ P_perp_a) / (b @ P_perp_a @ b)

assert np.allclose(P_a, P_a @ P_a)
assert np.allclose(P_b, P_b @ P_b)
assert np.allclose(P_a @ psi, 1.2*a)
assert np.allclose(P_b @ psi, -3.4*b)
display(array_to_latex(P_a, prefix="\proj{a} = "))
display(array_to_latex(P_b, prefix="\proj{b} = "))
```

(sec:MatrixExponential)=
### Matrix Exponential

Hermitian matrices and unitary matrices are closely related through the **matrix
exponential** and **matrix logarithm**.  I.e. if $\mat{H} = \mat{H}^{\dagger}$ is a
hermitian matrix, then the matrix exponential of $\I \mat{H}$ is unitary:

\begin{gather*}
  \mat{U} = e^{\I\mat{H}}.
\end{gather*}

Stated another way, the matrix exponential of an **anti-hermitian** matrix is unitary.

:::{note}
In quantum mechanics, a special hermitian matrix $\mat{H}$ called the **Hamiltonian**
represents the energy of the system and defines the **Schrödinger equation**

\begin{gather*}
  \I\hbar \pdiff{}{t}\ket{\psi(t)} = \mat{H}\ket{\psi(t)}.
\end{gather*}

If the Hamiltonian $\mat{H}$ does not depend on time, then the matrix exponential of this with an
additional numerical factor of $t/\I \hbar$ gives the explicit solution:

\begin{gather*}
  \mat{U}_t = e^{\mat{H}t/\I\hbar}, \qquad
  \ket{\psi(t)} = \mat{U}_t\ket{\psi(0)}.
\end{gather*}

The unitary matrix $\mat{U}_t$ is called the **propagator** and encodes the complete
solution to the time-independent the Schrödinger equation where the Hamiltonian
$\mat{H}$ does not depend on time.

In quantum computing, we generally work directly with the unitary propagator
$\mat{U}_t$, which represents the action of our computer.  The problem for physicists is
to find the appropriate physical system to get a Hamiltonian $\mat{H}$ that
exponentiates to the appropriate $\mat{U}_t$.
:::

:::{margin}
See {cite:p}`Moler:2003` for details about how matrix exponentials can be computed -- a
core problem in physics.  Quantum computers are expected to play a major role in solving
some of these problems since this is the essence of **quantum simulation**.
:::

Although somewhat complicated in general to do efficiently, the problem of computing the
matrix exponential of a hermitian matrix is quite straight forward.  The idea works for
any analytic function $f(\mat{H})$ whose Taylor series converges:

\begin{gather*}
  f(x) = f_0 + f_1x + f_2\frac{x^2}{2} + \cdots + f_n\frac{x^n}{n!} + \cdots\\
  f(\mat{H}) = f_0\mat{1} + f_1\mat{H}+ f_2\frac{\mat{H}^2}{2} + \cdots + f_n\frac{\mat{H}^n}{n!} + \cdots.
\end{gather*}

This can in principle be applied to any matrix, but we can get a simple form for
hermitian matrices by using the fact that the eigenvectors of the hermitian matrix
$\mat{H}\ket{n} = \ket{n}\lambda_n$ can be chosen to form a complete orthonormal basis
$\ket{n}$, we can write these as the columns of a unitary matrix $\mat{V}$:

\begin{gather*} 
  \mat{V} = 
  \begin{pmatrix} 
    \ket{1} & \ket{2} & \cdots & \ket{N}
  \end{pmatrix}.
\end{gather*}

Using this, we can write:

\begin{gather*}
  \mat{H}\mat{V} =
  \mat{V}\mat{D}, \qquad
  \mat{H} =
  \mat{V}\mat{D}\mat{V}^\dagger,
\end{gather*}

:::{hint}
:class: dropdown

\begin{gather*}
  \mat{H}\mat{V} =
  \begin{pmatrix}
    \mathstrut\\ 
    \mat{H}\ket{1} & \mat{H}\ket{2} & \cdots & \mat{H}\ket{N}\\
    \mathstrut
  \end{pmatrix} =\\
  =
  \begin{pmatrix}
    \mathstrut\\
    \ket{1}\lambda_1 & \ket{2}\lambda_2 & \cdots & \ket{N}\lambda_N\\
    \mathstrut
  \end{pmatrix} = \mat{V}\mat{D}.
\end{gather*}
:::
where $\mat{D}$ is diagonal

:::{margin}
One can apply this approach further to non-hermitian matrices $\mat{A}$ as long as are **normal**
meaning that they **commute** with their hermitian conjugate $\mat{A}\mat{A}^\dagger =
\mat{A}^\dagger\mat{A}$.  This is enough to ensure that they have a complete (though not
necessarily orthonormal) set of eigenvectors and we can write:

\begin{gather*}
  \mat{A} = \mat{S}\mat{D}\mat{S}^{-1}.
\end{gather*}
:::
\begin{gather*}
  \mat{D} = 
  \begin{pmatrix} 
    \lambda_1 \\
    & \lambda_2 \\
    & & \cdots \\
    & & & \lambda_N
  \end{pmatrix}.
\end{gather*}

This is called **diagonalization** or we say that we are **diagonalizing** $\mat{H}$.
We can thus express any of the powers $\mat{H}^n$ as:
\begin{gather*}
  \mat{H}^n = (\mat{V}\mat{D}\mat{V}^\dagger)
  (\mat{V}\mat{D}\mat{V}^\dagger)\cdots
  (\mat{V}\mat{D}\mat{V}^\dagger)
  = \mat{V}\mat{D}\mat{D}\cdots \mat{D}\mat{V}^\dagger = \\
  = \mat{V}\mat{D}^n\mat{V}^\dagger.
\end{gather*}

Hence, we can write

\begin{gather*}
  f(\mat{H}) = \mat{V}f(\mat{D})\mat{V}^\dagger
  = \mat{V}
  \begin{pmatrix} 
    f(\lambda_1) \\
    & f(\lambda_2) \\
    & & \cdots \\
    & & & f(\lambda_N)
  \end{pmatrix}
  \mat{V}^\dagger
\end{gather*}

by factoring the remaining factors of $\mat{V}$ and $\mat{V}^\dagger$ on the left and
right of each term in the series.

We now see that the $f(\mat{H})$ and $\mat{H}$ have the same eigenvectors $\ket{n}$ and
that the eigenvalues of $f(\mat{H})$ are just the $f(\lambda_n)$.  Thus, if $\mat{H} =
\mat{H}^\dagger$ is hermitian with real eigenvalues $\lambda_n$, then the matrix
exponential of $\I\mat{H}$ has eigenvalues that are phases $e^{\I\lambda_n}$.

#### Review
<orbit-reviewarea>
  <orbit-prompt
    cloze="The eigenvalues of a hermitian matrix are {real}."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The eigenvectors of a hermitian matrix can be
    chosen to form a {complete orthonormal} basis."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The hermitian conjugate of a unitary matrix is its {inverse}."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The eigenvalues of a unitary matrix are {phases $\exp(i\theta)$}."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The eigenvectors of a unitary matrix can be chosen to form a {complete
  orthonormal} basis."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The matrix exponential can be defined and computed from the T{aylor series}."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The matrix exponential can be defined and computed from its {Taylor} series."
  ></orbit-prompt>
  <orbit-prompt
    cloze="The matrix exponential of an {anti-hermitian} matrix is a unitary matrix."
  ></orbit-prompt>
  <orbit-prompt
    cloze="In quantum mechanics, the hermitian matrix $H$ on the right-hand-side of the Schrödinger equation is called the {Hamiltonian}."
  ></orbit-prompt>
  <orbit-prompt
    question="In quantum mechanics, the unitary propagator $U_t$ is the matrix exponential of the hermitian Hamiltonian $H$ multiplied by what factor?"
    answer="$$\frac{t}{i \hbar}, \qquad U_t = e^{H t/i\hbar}$$."
  ></orbit-prompt>
</orbit-reviewarea>

## Matrix Decompositions

(sec:QR-Decomposition)=
### QR Decomposition

Given a set of vectors $\ket{a_0}$, $\ket{a_1}$, etc., one can follow the Gram-Schmidt
procedure to produce an orthonormal basis.  To simplify the notation, we define the
"normalization" operator $\newcommand{\N}{\mathcal{N}}\N(\ket{a})$

\begin{gather*}
  \N(\ket{a}) = \frac{\ket{a}}{\norm{\ket{a}}} = \frac{\ket{a}}{\sqrt{\braket{a|a}}}.
\end{gather*}

We can now express the Gram-Schmidt procedure as a sequence of steps projecting the
remaining basis vectors into an orthogonal subspace, then 

\begin{align*}
  \ket{0} &= \N(\ket{a_0}) = \frac{\ket{a_0}}{\sqrt{\braket{a_0|a_0}}} \\
  \ket{1} &= \N\Bigl((\op{1} - \ket{0}\bra{0})\ket{a_1}\Bigr)
          = \frac{\ket{a_1} - \ket{0}\braket{0|a_1}}{\norm{\ket{a_1} - \ket{0}\braket{0|a_1}}}\\  
  \ket{2} &= \N\Bigl((\op{1} - \ket{0}\bra{0} - \ket{1}\bra{1})\ket{a_2}\Bigr)\\
          &\;\;\vdots
\end{align*}

The QR decomposition numerically encodes the results of this process in terms of a
unitary matrix $\mat{Q}$ and an upper triangular matrix $\mat{R}$:

\begin{gather*}
  \mat{A} = \mat{Q}\mat{R}\\
  \underbrace{
    \begin{pmatrix}
    &\\
    \ket{a_0} & 
    \ket{a_1} & 
    \dots\\
    &
  \end{pmatrix}}_{\mat{A}}
  =
  \underbrace{
    \begin{pmatrix}
    &\\
    \ket{0} & 
    \ket{1} & 
    \dots\\
    &
  \end{pmatrix}
  }_{\mat{Q}}
  \underbrace{
  \begin{pmatrix}
    R_{00} & R_{01} & R_{02} & \dots\\
    & R_{11} & R_{12} & \dots\\
    & & R_{22} & \dots\\
    & & & \ddots
  \end{pmatrix}}_{\mat{R}}\\
  \begin{aligned}
    \ket{a_0} &= \ket{0}R_{00} = \sqrt{\braket{a_0|a_0}}\\
    \ket{a_1} &= \ket{0}R_{01} + \ket{1}R_{11}\\
    \ket{a_2} &= \ket{0}R_{02} + \ket{1}R_{12} + \ket{2}R_{22}\\
    &\;\;\vdots
  \end{aligned}
\end{gather*}

Recall that the columns of a **unitary matrix** $\mat{Q}$ form an orthonormal basis.

### Singular Value Decomposition

:::{margin}
Typically, the singular values are sorted in decreasing order $\sigma_1 \geq \sigma_2
\geq \sigma_3 \cdots$.  The **rank** of the matrix is the number of positive singular
values (or the number of singular values greater than some threshold $\sigma_{k} > \epsilon$).
:::
Any matrix $\mat{A}$ can be decomposed into its **singular value decomposition** or
**SVD**:

\begin{gather*}
  \mat{A} = \mat{U}\mat{D}\mat{V}^{\dagger},
\end{gather*}

where $\mat{U}$ and $\mat{V}$ are unitary, and $\mat{D}$ is diagonal with non-negative
entries $\sigma_n$ called the **singular values**.

Note: $\mat{A}$ need not be square:

\begin{gather*}
  \underbrace{\mat{A}}_{2\times 3}
  =
  \underbrace{
    \begin{pmatrix}
      &\\
      \ket{u_1} & 
      \ket{u_2} \\
      &
    \end{pmatrix}
  }_{\mat{U}}
  \underbrace{
    \begin{pmatrix}
      \sigma_1 & 0 & 0\\
      0 & \sigma_2 & 0\\
    \end{pmatrix}
  }_{\mat{D}}
  \underbrace{
    \begin{pmatrix}
      & \bra{v_1} &\\
      & \bra{v_2} &\\
      & \bra{v_3} &
    \end{pmatrix}
  }_{\mat{V}^\dagger},\\
  \underbrace{\mat{A}}_{3\times 2}
  =
  \underbrace{
    \begin{pmatrix}
      &\\
      \ket{u_1} & 
      \ket{u_2} &
      \ket{u_3} \\
      &
    \end{pmatrix}
  }_{\mat{U}}
  \underbrace{
    \begin{pmatrix}
      \sigma_1 & 0\\
      0 & \sigma_2\\
      0 & 0
    \end{pmatrix}
  }_{\mat{D}}
  \underbrace{
    \begin{pmatrix}
      & \bra{v_1} &\\
      & \bra{v_2} &\\
    \end{pmatrix}
  }_{\mat{V}^\dagger}
\end{gather*}

The SVD can be computed efficiently by finding the eigenvalues and eigenvectors of the
two hermitian matrices $\mat{A}^{\dagger}\mat{A}$ and $\mat{A}\mat{A}^\dagger$
respectively.

## Direct Sum and Tensor Product

Consider two vectors $\ket{a}$ (with $N_a$ components) and $\ket{b}$ (with $N_b$
components).  We can form two new vectors from these:

1. The **direct sum** $\ket{a} \oplus \ket{b}$ which has $N_a + N_b$ components.
2. The **tensor product** $\ket{a} \otimes \ket{b}$ which has $N_a N_b$ components.

The first can be visualized as stacking the vectors on top of each other:

\begin{gather*}
    \ket{a} \oplus \ket{b} = \begin{pmatrix}
      \ket{a}\\
      \ket{b}
    \end{pmatrix}.
\end{gather*}

The direct sum is useful when we decompose a vector space into subspaces, for example,
when using projection matrices.  The corresponding matrices are block diagonal:

\begin{gather*}
  \mat{A} \oplus \mat{B} = \begin{pmatrix}
    \mat{A} \\
    & \mat{B}
  \end{pmatrix}
\end{gather*}

The **tensor product** is a little more tricky to visualize:

\begin{gather*}
  \ket{a} = \begin{pmatrix}
    a_1\\
    a_2\\
    \vdots\\
    a_{N_a}
  \end{pmatrix}, \qquad
  \ket{b} = \begin{pmatrix}
    b_1\\
    b_2\\
    \vdots\\
    b_{N_b}
  \end{pmatrix}, \\
  \ket{a}\otimes\ket{b} = 
  \begin{pmatrix}
    a_1\ket{b}\\
    a_2\ket{b}\\
    \vdots\\
    a_{N_a}\ket{b}
  \end{pmatrix}
  =
  \begin{pmatrix}
    a_1b_1\\
    a_1b_2\\
    \vdots\\
    a_2b_{N_b}\\
    a_2b_1\\
    a_2b_2\\
    \vdots\\
    a_2b_{N_b}\\
    \vdots\\
    a_{N_a}b_1\\
    a_{N_a}b_2\\
    \vdots\\
    a_{N_a}b_{N_b}
  \end{pmatrix}.
\end{gather*}

It is easier to think of $\ket{a}\otimes\ket{b}$ in terms of grouped indices:

\begin{gather*}
  [\ket{a}\otimes\ket{b}]_{ij} = a_{i}b_{j}.
\end{gather*}

We can generalize to matrices, first in index form:

\begin{gather*}
  [\mat{A}\otimes\mat{B}]_{ik,jl} = A_{ij}B_{kl}.
\end{gather*}

Note the grouping carefully.  Visually, this corresponds to

\begin{gather*}
  \mat{A} = \begin{pmatrix}
    a_{11} & a_{12} & \cdots\\
    a_{21} & a_{22} & \cdots\\
    \vdots & \vdots & \ddots\\
  \end{pmatrix}, \qquad
  \mat{B} = \begin{pmatrix}
    b_{11} & b_{12} & \cdots\\
    b_{21} & b_{22} & \cdots\\
    \vdots & \vdots & \ddots\\
  \end{pmatrix}, \\
  \mat{A}\otimes\mat{B}
  =
  \begin{pmatrix}
    a_{11}\mat{B} & a_{12}\mat{B} & \cdots\\
    a_{21}\mat{B} & a_{22}\mat{B} & \cdots\\
    \vdots & \vdots & \ddots
  \end{pmatrix}.
\end{gather*}

These satisfy the following properties:

\begin{gather*}
  (\mat{A}\ket{a})\otimes(\mat{B}\ket{b}) = 
  (\mat{A}\otimes \mat{B})(\ket{a}\otimes\ket{b}),\\
  \Tr \mat{A}\otimes\mat{B} = 
  (\Tr \mat{A})(\Tr\mat{B}).
\end{gather*}

We also talk about the **tensor product** of vector spaces, so that if $\ket{a} \in \mathcal{H}_a$
and $\ket{b} \in \mathcal{H}_b$, then $\ket{a}\otimes\ket{b} \in \mathcal{H}_a \otimes
\mathcal{H}_b$.  If we have an orthonormal basis $\{\ket{a_i}\}$ for $\mathcal{H}_a$ and
an orthonormal basis
$\{\ket{b_j}\}$ for $\mathcal{H}_b$, then the basis $\ket{ij} =
\ket{a_i}\otimes\ket{b_j}$ is an orthonormal basis for $\mathcal{H}_a\otimes
\mathcal{H}_b$.

It is useful to be able to work with the tensor product numerically.  The function
{py:func}`numpy.kron`, named for [Kronecker product][], directly implements the tensor
product, but I find it clearer to explicitly work with indices using
{py:func}`numpy.einsum`.

:::{margin}
We switch here to indexing from $0$ to $N-1$ to match Python.
:::
First note that we can "combine indices" using {py:func}`numpy.reshape`.  Thus, if we
have a matrix $A_{i,j}:$

\begin{gather*}
  \mat{A} = \begin{pmatrix}
    A_{00} & A_{01} & A_{02}\\
    A_{10} & A_{11} & A_{12}
  \end{pmatrix},
\end{gather*}

this will "ravel" to a vector

\begin{gather*}
  \begin{pmatrix}
    A_{00} \\
    A_{01} \\
    A_{02} \\
    A_{10} \\
    A_{11} \\
    A_{12}
  \end{pmatrix}.
\end{gather*}

although with a single index (so it appears flat):

```{code-cell}
A = np.array([
    [1000, 1001, 1002],
    [1010, 1011, 1012]])
display(A)
display(A.ravel())
```

Here is the tensor product demonstrating the index ordering:

\begin{gather*}
  [\mat{A}\otimes \mat{B}]_{ik,jl} = A_{ij}B_{kl}
\end{gather*}

```{code-cell}
A = np.array([
    [0, 1],
    [2, 3]])
B = np.array([
    [10, 20],
    [30, 40]])
display(np.kron(A, B))
display(np.einsum('ij,kl->ikjl', A, B).reshape((4,4)))
```

### Schmidt Decomposition

The **[Schmidt decomposition][]** is an important decomposition of a tensor product space
which states that for any vector $\ket{w} \in \mathcal{H}_a \otimes \mathcal{H}_b$ we
can find a set of orthonormal bases $\{\ket{u_i}\}$ and $\{\ket{v_j}\}$ for $\mathcal{H}_a$ and
an $\mathcal{H}_b$ respectively such that

\begin{gather*}
  \ket{w} = \sum_{i=1}^{\min(m,n)} \alpha_i \ket{u}_{i}\otimes\ket{v_{i}}
\end{gather*}

where $\min(m, n)$ is the dimension of the smaller space $m = \dim\mathcal{H}_a$ or $n =
\dim\mathcal{H}_b$ respectively.

:::{admonition} Do it!  Prove this using the SVD.
:class: dropdown

First express $\ket{w}$ in some orthonormal basis

\begin{gather*}
  \ket{w} = \sum_{ij} \beta_{ij}\ket{a_i}\otimes\ket{b_j}.
\end{gather*}

Use the fact that any matrix has a SVD to factor the matrix of coefficients $\beta_{ij}$
and use this to explicitly construct the bases $\{\ket{u_i}\}$ and $\{\ket{v_j}\}$ and
coefficients $\alpha_i = \sigma_i$ that provides the **Schmidt decomposition**.
:::


[Pauli matrices]: https://en.wikipedia.org/wiki/Pauli_matrices
[direct sum]: <https://en.wikipedia.org/wiki/Direct_sum>
[tensor product]: <https://en.wikipedia.org/wiki/Tensor_product>
[Schmidt decomposition]: <https://en.wikipedia.org/wiki/Schmidt_decomposition>
[Kronecker product]: <https://en.wikipedia.org/wiki/Kronecker_product>
[positive operator]: <https://en.wikipedia.org/wiki/Positive_operator_(Hilbert_space)>
[Projection]: <https://en.wikipedia.org/wiki/Projection_(linear_algebra)>
