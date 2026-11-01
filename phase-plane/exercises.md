---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/phase-plane-exercises.pdf
    id: phase-plane-exercises-pdf
downloads:
  - id: phase-plane-exercises-pdf
    title: Download PDF
---

# Exercises: The Phase Plane

The exercises follow the pages of this part, and the last section adds two
classical models that are not in the notes. Some are calculations by hand, and
a few ask for a short computation in Python (or any other language).

## Linear systems

:::{exercise} Three matrices
:label: exr-pp-matrices

For each matrix, find the trace and the determinant and classify the origin
with {prf:ref}`prop-pp-trace-det`. Then find the solution with the given
initial value, using {prf:ref}`prop-pp-eigen-solutions`,
{prf:ref}`prop-pp-complex` or {prf:ref}`prop-pp-repeated`.

$$
A = \begin{pmatrix} 1 & 2 \\ 3 & 2 \end{pmatrix},\ \mathbf u(0) = \begin{pmatrix} 3 \\ 2 \end{pmatrix};
\qquad
B = \begin{pmatrix} -1 & -4 \\ 1 & -1 \end{pmatrix},\ \mathbf u(0) = \begin{pmatrix} 0 \\ 1 \end{pmatrix};
\qquad
C = \begin{pmatrix} -3 & 1 \\ -1 & -1 \end{pmatrix},\ \mathbf u(0) = \begin{pmatrix} 1 \\ 0 \end{pmatrix}.
$$

For $B$, does the flow turn clockwise or counterclockwise?
:::

:::{exercise} A path through the trace-determinant plane
:label: exr-pp-path

Let

$$
J(k) = \begin{pmatrix} -1 & k \\ -1 & -2 \end{pmatrix}.
$$

1. For which $k$ is the origin a saddle, a stable node, or a stable spiral?
2. Mark the path that $(T, \Delta)$ follows as $k$ runs from $-5$ to $5$.
3. {prf:ref}`prop-pp-no-spirals` says that cross effects of the same sign rule
   out spirals. Check this here. Do cross effects of opposite signs guarantee a
   spiral?
:::

:::{exercise} A tethered cell
:label: exr-pp-spring

A cell is held back by an elastic tether that pulls it toward its anchor with
the force $-kx$. With Stokes drag, Newton's second law {eq}`eq-eom-full`
becomes $m\ddot x + \gamma\dot x + kx = 0$.

1. With $v = \dot x$, write it as a linear system for $(x, v)$ and find its
   trace and determinant.
2. When does the cell spiral into its rest position, and when does it creep in
   without oscillating? Which case applies to a crawling cell, with the
   numbers of [](../modeling/why-inertia-is-irrelevant.md)?
3. Show that for small $m$ the eigenvalues are approximately $-k/\gamma$ and
   $-\gamma/m$. Which one is the motion of the force balance
   $\gamma\dot x + kx = 0$, and what does the other one describe? Compare with
   {prf:ref}`prop-pp-fast-slow`.
:::

## Drug and bacteria in the plane

:::{exercise} Without collective resistance
:label: exr-pp-no-collective

Suppose the cells do not destroy the drug, so that $N_c \to \infty$ in
{eq}`eq-pp-model` and $\dot d = \delta(1 - d)$.

1. Draw the nullclines and find the steady states.
2. Compute the Jacobian at each steady state and classify them.
3. Is there an inoculum effect? Which bifurcation is left, and which one has
   disappeared?
:::

:::{exercise} The separatrix at the threshold
:label: exr-pp-separatrix

Take $r = 1$/day, $K = 10$, $N_c = 1$, $a = 2$/day and the TB clearance rate
$\delta = 17$/day in {eq}`eq-pp-model`. The threshold of the reduced model is
$N_- \approx 1.298$ (Inoculum and dose, in
[](../one-dimensional/exercises.md)).

1. Find the saddle $(N_-, \bar d(N_-))$, its Jacobian, and its eigenvalues.
   Compare the positive eigenvalue with the slope $f'(N_-) \approx 0.362$ of
   the reduced equation.
2. Find the eigenvector of the negative eigenvalue, and the slope
   $\dd d/\dd N$ of the separatrix at the saddle.
3. Use this slope to estimate how much the threshold inoculum rises when the
   cells meet the full dose, $d(0) = 1$, instead of the quasi-steady level.
   Compare with {prf:ref}`prop-qss-threshold-shift`.
:::

## Quasi-steady state

:::{exercise} Messenger RNA and protein
:label: exr-pp-mrna-protein

An mRNA $m$ is made at rate $\alpha$ and degraded at rate $\delta_m$; each
mRNA is translated into protein $p$ at rate $\beta$, and the protein is lost at
rate $\delta_p$:

$$
\dot m = \alpha - \delta_m m, \qquad \dot p = \beta m - \delta_p p .
$$

In bacteria mRNA lasts minutes, while proteins last as long as the cell cycle,
so $\ve = \delta_p/\delta_m$ is small.

1. Identify the fast and the slow variable, the critical manifold and the
   reduced equation. Is the critical manifold attracting?
2. Solve both the full model and the reduced one with $m(0) = p(0) = 0$.
3. Show that the difference between the two protein curves is at most
   $p^\ast\ve/(1 - \ve)$, with $p^\ast = \alpha\beta/(\delta_m\delta_p)$.
   Which way does the full model err, and why?
:::

:::{exercise} A gene that switches itself on
:label: exr-pp-switch

A transcription factor activates its own gene. In scaled units, with $u$ the
protein and $m$ its mRNA,

$$
\ve\,\dot m = s(u) - m, \qquad \dot u = m - b\,u,
\qquad s(u) = \frac{u^2}{1 + u^2} .
$$

1. Find the critical manifold, show that it attracts, and write the reduced
   equation $\dot u = F(u)$.
2. Find the steady states of the reduced equation for $b = 0.4$ and classify
   them. For which $b$ is the gene a switch, with two stable states?
3. Show that the full model has the same steady states, and that a steady
   state is a saddle exactly when $F' > 0$ there, for every $\ve > 0$.
   Compare with {prf:ref}`prop-qss-inoc-exact`.
:::

## FitzHugh–Nagumo

:::{exercise} A symmetry
:label: exr-pp-fhn-symmetry

Show that if $(w(t), v(t))$ solves {eq}`eq-fhn` with parameter $a$, then
$(-w(t), -v(t))$ solves it with parameter $-a$. What does an excitable system
at $a = -1.1$ look like, and which way does it respond to kicks?
:::

:::{exercise} The threshold for small ε
:label: exr-pp-fhn-threshold

1. Evaluate the threshold {eq}`eq-fhn-threshold` at $a = 1.1$ and $a = 1.5$.
2. Write a program that kicks the activator from rest by $\Delta v$, at
   $a = 1.1$, and decides whether a pulse follows (for example, whether $v$
   exceeds $0.5$ within $4$ time units). Find the threshold kick by bisection
   for $\ve = 0.05$, $0.02$, $0.01$ and $0.005$, and compare with part 1. Why
   is the threshold larger for $\ve > 0$?
:::

:::{exercise} The inhibitor decays
:label: exr-pp-fhn-decay

FitzHugh's model has a decay term in the inhibitor equation
({prf:ref}`rmk-fhn-original`):

$$
\ve\,\dot v = c(v) - w, \qquad \dot w = v + a - b\,w, \qquad b > 0 .
$$

1. Show that a large box $|w| \le W$, $|v| \le V$ traps every trajectory. Why
   did {eq}`eq-fhn`, with $b = 0$, need the hexagon of
   {prf:ref}`prop-fhn-trap` instead?
2. Show that at a steady state $T = c'(v^\ast)/\ve - b$ and
   $\Delta = \big(1 - b\,c'(v^\ast)\big)/\ve$. Conclude that a steady state is a
   saddle exactly when the S is steeper than the $w$-nullcline there, and that
   saddles need $b > 1$.
3. Take $a = 0$ and $b = 2$. Find the three steady states and classify them.
   What kind of system is this?
:::

## Two classical models

:::{exercise} Lotka-Volterra
:label: exr-pp-lotka-volterra

The prey $x$ and the predators $y$ obey

$$
\dot x = x(1 - y), \qquad \dot y = y(x - 1) .
$$

1. Find the steady states and draw the nullclines with the directions of the
   arrows on them. Which way does the flow turn?
2. Compute the Jacobian and classify the steady states. What does
   linearization fail to decide?
3. Show that $H(x, y) = x - \ln x + y - \ln y$ is constant along solutions
   with $x, y > 0$. (Hint: separate variables in $\dd y/\dd x$.) What does
   this say about the trajectories near $(1, 1)$?
:::

:::{exercise} Predators and crowded prey
:label: exr-pp-predator-prey

Give the prey of the Lotka-Volterra model above logistic growth:

$$
\dot u = u(1 - u) - \alpha uv, \qquad \dot v = -\delta v + \beta uv,
\qquad \alpha, \beta, \delta > 0, \quad \delta < \beta .
$$

1. Find the nullclines and the three steady states, and the directions of the
   arrows on the nullclines. Where does the assumption $\delta < \beta$ enter?
2. Classify the steady states. When is the coexistence state a spiral?
3. Let $W = \beta u + \alpha v$. Show that
   $\dot W \le \beta(1 + \delta)^2/4 - \delta W$, and conclude that the
   triangle $u, v \ge 0$, $W \le M$ traps every trajectory if $M$ is large
   enough.
:::

