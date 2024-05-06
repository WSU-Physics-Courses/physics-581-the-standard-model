---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
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

(sec:demonstration)=
JupyterBook Demonstration
=========================

This document demonstrates some of the features provided by the documentation.  It can
serve as a starting point for your own documentation.

This documentation is formatted using an extension to [Markdown][] called [MyST][] which
allows full interoperability with [Sphinx][], the system used to build the documentation.
Here we demonstrate some features.  For more details, please see:

* [Jupyter Book][]
* [Jupyter Book with Sphinx][]
* [MyST Cheatsheet][]

## Math

You can use LaTeX math with [MathJaX][].

### Example
Here we implement a simple example of a pendulum of mass $m$ hanging down distance $r$
from a pivot point in a gravitational field $g>0$ with coordinate $\theta$ so that the
mass is at $(x, z) = (r\sin\theta, -r\cos\theta)$ and $\theta=0$ is the downward
equilibrium position: 

\begin{gather*}
  L(\theta, \dot{\theta}, t) = \frac{m}{2}r^2\dot{\theta}^2 + mgr\cos\theta\\
  p_{\theta} = \pdiff{L}{\dot{\theta}} = mr^2 \dot{\theta}, \qquad
  \dot{\theta}(\theta, p_{\theta}, t) = \frac{p_{\theta}}{mr^2},\\
  H(\theta, p_{\theta}) = p_{\theta}\dot{\theta} - L = \frac{p_{\theta}^2}{2mr^2} - mgr\cos\theta,\\
  \vect{y} = \begin{pmatrix}
    \theta\\
    p_{\theta}
  \end{pmatrix},\qquad
  \dot{\vect{y}} = 
  \begin{pmatrix}
    0 & 1\\
    -1 & 0
  \end{pmatrix}
  \cdot
  \begin{pmatrix}
    \pdiff{H}{\theta}\\
    \pdiff{H}{p_{\theta}}
  \end{pmatrix}
  =
  \begin{pmatrix}
    p_{\theta}/mr^2\\
    -mgr\sin\theta
  \end{pmatrix}.
\end{gather*}

```{note}
On [CoCalc][], only a subset of the math will render in the live Markdown editor.  This
uses [KaTeX][], which is much faster than [MathJaX][], but more limited.  The final
documentation will use [MathJaX][], and loads the macros defined in
`Docs/_static/math_defs.tex`.
```

### SI Units

In LaTeX, I like to use the [siunitx][] package.  This is not currently supported in
MathJax (they [could not get funding](https://github.com/mathjax/MathJax/issues/447)),
but can be mocked up with [some defines](https://github.com/KaTeX/KaTeX/issues/817).
The idea is that we use the macros for HTML, but substitute the appropriate
`\usepackage{siunitx}` in the LaTeX files. *(This second step is currently not implemented)*

:::{admonition] To Do: Add a mechanism for different MathJaX vs LaTeX defines.
:::

[math]: <https://jupyterbook.org/content/math.html>
[KaTeX]: <https://katex.org/>
[siunitx]: <https://ctan.org/pkg/siunitx>

## Jupyter Notebooks

This document is actually a [Jupyter][] notebook, synchronized with [Jupytext][].  The top
of the document contains information about the kernel that should be used etc.  These
are parsed using [MyST-NB][] and the output of cells will be displayed.  For example, here
is a plot demonstrating [Liouville's Theorem][].

```{code-cell}
:tags: [hide-input, full-width]

plt.rcParams['figure.dpi'] = 300
from scipy.integrate import solve_ivp

alpha = 0.0
m = r = g = 1.0
w = g/r
T = 2*np.pi / w   # Period of small oscillations.


def f(t, y):
    # We will simultaneously evolve N points.
    N = len(y)//2
    theta, p_theta = np.reshape(y, (2, N))
    dy = np.array([p_theta/m/r**2*np.exp(-alpha*t),
                   -m*g*r*np.sin(theta)*np.exp(alpha*t)])
    return dy.ravel()


# Start with a circle of points in phase space centered here
def plot_set(y, dy=0.1, T=T, phase_space=True, N=10, Nt=5, c='C0', 
             Ncirc=1000, max_step=0.01, _alpha=0.7, 
             fig=None, ax=None):
    """Plot the phase flow of a circle centered about y0.
    
    Arguments
    ---------
    y : (theta0, ptheta_0)
        Center of initial circle.
    dy : float
        Radius of initial circle.
    T : float
        Time to evolve to.
    phase_space : bool
        If `True`, plot in phase space $(q, p)$, otherwise plot in
        the "physical" phase space $(q, P)$ where $P = pe^{-\alpha t}$.
    N : int
        Number of points to show on circle and along path.
    Nt : int
        Number of images along trajectory to show.
    c : color
        Color.
    alpha : float
        Transparency of regions.
    max_step : float
        Maximum spacing for times dt.
    Ncirc : int
        Minimum number of points in circle.
    """
    global alpha
    skip = int(np.ceil(Ncirc // N))
    N_ = N * skip
    th = np.linspace(0, 2*np.pi, N_ + 1)[:-1]
    z = dy * np.exp(1j*th) + np.asarray(y).view(dtype=complex)
    y0 = np.ravel([z.real, z.imag])

    skip_t = int(np.ceil(T / max_step / Nt))
    Nt_ = Nt * skip_t + 1
    t_eval = np.linspace(0, T, Nt_)
    res = solve_ivp(f, t_span=(0, T), y0=y0, t_eval=t_eval)
    assert Nt_ == len(res.t)
    thetas, p_thetas = res.y.reshape(2, N_, Nt_)
    ylabel = r"$p_{\theta}$"
    if phase_space:
        ylabel = r"$P_{\theta}=p_{\theta}e^{-\alpha t}$"
        p_thetas *= np.exp(-alpha * res.t)
    
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(thetas[::skip].T, p_thetas[::skip].T, "-k", lw=0.1)
    for n in range(Nt+1):
        tind = n*skip_t
        ax.plot(thetas[::skip, tind], p_thetas[::skip, tind], '.', ms=0.5, c=c)
        ax.fill(thetas[:, tind], p_thetas[:, tind], c=c, alpha=_alpha)
    ax.set(xlabel=r"$\theta$", ylabel=ylabel, aspect=1)
    return fig, ax


fig, ax = plt.subplots(figsize=(10,5))


for n, y in enumerate(np.linspace(0.25, 1.75, 6)):
    plot_set(y=(y, y), c=f"C{n}", ax=ax)
```

:::{sidebar} Phase flow for a pendulum.

The small colored circular region is evolved forward in
time according to the Hamiltonian equations.  The trajectories of 10 equally spaced
points are shown at 6 different times equally spaced from $t=0$ to $T =
2\pi\sqrt{r/g}$, the period of oscillation for small amplitude modes.  At small energies
(blue circle), the period is almost independent of the amplitude, and the circle stays
together.  As the energy increases (orange, green, and red), the period starts
to depend more sensitively on the amplitude, and the initial circular region starts to
shear.  The purple region divides the two qualitatively different regions of rotation
and libration, and gets stretched as some points oscillate and others orbit. 
The areas remain constant, despite this stretching, as a demonstration of
Liouville's theorem.
:::

`````{admonition} Details
We start a code-block with:

````markdown
```{code-cell}
:tags: [hide-input, full-width]
...
```
````
This indicates that it is a code cell, to be run with python 3, and has some
[tags](https://myst-nb.readthedocs.io/en/latest/use/hiding.html) to make the cell and
output full width, but hiding the code.  (There is a link to click to show the code.)
`````

If you need more control, you can [glue][] the output to a variable, then load it as a
proper figure, insert it in a table, etc.

## Manim


We provide experimental support for [Manim Community][] -- the community edition of the animation
software used by the [3Blue1Brown](https://www.3blue1brown.com/) project.  Here we
demonstrate an animation of the [Jacobi elliptic functions][].

```{code-cell}
:tags: [hide-cell]

import manim.utils.ipython_magic
```

```{code-cell}
:tags: [hide-input]

%%manim -v WARNING --progress_bar None -qm JacobiEllipse
from scipy.special import ellipj, ellipkinc
from manim import *

config.media_width = "100%"
config.media_embed = True

my_tex_template = TexTemplate()
with open("_static/math_defs.tex") as f:
    my_tex_template.add_to_preamble(f.read())

def get_scd(phi, k):
    m = k**2
    u = ellipkinc(phi, m)
    s, c, d, phi_ = ellipj(u, m)
    return s, c, d

def get_xyr(phi, k):
    s, c, d = get_scd(phi, k)
    r = 1/d
    x, y = r*c, r*s
    return x, y, r

class JacobiEllipse(Scene):
    def construct(self):
        config["tex_template"] = my_tex_template
        config["media_width"] = "100%"
        class colors:
            ellipse = YELLOW
            x = BLUE
            y = GREEN
            r = RED
            
        phi = ValueTracker(0.01)
        k = 0.8
        axes = Axes(x_range=[0, 2*np.pi, 1], 
                    y_range=[-2, 2, 1],
                    x_length=6, 
                    y_length=4,
                    axis_config=dict(include_tip=False),
                    x_axis_config=dict(numbers_to_exclude=[0]),
                   ).add_coordinates()
        
        
        plane = PolarPlane(radius_max=2).add_coordinates()

        x = lambda phi: get_xyr(phi, k)[0]
        y = lambda phi: get_xyr(phi, k)[1]
        r = lambda phi: get_xyr(phi, k)[2]
        
        g_ellipse = plane.plot_polar_graph(r, [0, 2*np.pi], color=colors.ellipse)
        
        points_colors = [(x, colors.x), (y, colors.y), (r, colors.r)]
        x_graph, y_graph, r_graph = [axes.plot(_x, color=_c) for _x, _c in points_colors]
        graphs = VGroup(x_graph, y_graph, r_graph)
        
        dot = always_redraw(lambda:
            Dot(plane.polar_to_point(r(phi.get_value()), phi.get_value()), 
                fill_color=colors.ellipse, 
                fill_opacity=0.8))

        @always_redraw
        def lines():
            c2p = plane.coords_to_point
            phi_ = phi.get_value()
            x_, y_ = x(phi_), y(phi_)
            return VGroup(
                Line(c2p(0, 0), c2p(x_, 0), color=colors.x),
                Line(c2p(0, y_), c2p(x_, y_), color=colors.x),
                Line(c2p(0, 0), c2p(0, y_), color=colors.y),
                Line(c2p(x_, 0), c2p(x_, y_), color=colors.y),
                Line(axes.c2p(phi_, y_), c2p(x_, y_), color=colors.y),
                Line(c2p(0, 0), c2p(x_, y_), color=colors.r),
            ).set_opacity(0.8)

        dots = always_redraw(lambda:
            VGroup(*(Dot(axes.c2p(phi.get_value(), _x(phi.get_value())), fill_color=_c, fill_opacity=1)
                     for _x, _c in points_colors)))

        a_group = VGroup(axes, dots, graphs)
        p_group = VGroup(plane, g_ellipse, dot, lines)
        a_group.shift(RIGHT*2)
        p_group.shift(LEFT*4)
    
        labels = VGroup(
            axes.get_graph_label(
                x_graph, label=r"x=\cn(u, k)", color=colors.x,
                x_val=2.5, direction=DOWN).shift(0.2*LEFT),
            axes.get_graph_label(
                y_graph, label=r"y=\sn(u, k)", color=colors.y,
                x_val=4.5, direction=DR),
            axes.get_graph_label(
                r_graph, label=r"r=1/\dn(u, k)", color=colors.r,
                x_val=2, direction=UR),
        )
        #labels.next_to(a_group, RIGHT)
        self.add(p_group, a_group, labels)
        self.play(phi.animate.set_value(2*np.pi), run_time=5, rate_func=linear)
```

```{warning}
We need a Manim >= 0.15.0 to get this to work, resolving [issue
2441](https://github.com/ManimCommunity/manim/issues/2441) by embeding the video in the
output.
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
