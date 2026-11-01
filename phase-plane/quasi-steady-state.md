---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/quasi-steady-state.pdf
    id: phase-plane-quasi-steady-state-pdf
downloads:
  - id: phase-plane-quasi-steady-state-pdf
    title: Download PDF
---

# Quasi-Steady State

:::{tip} Big Idea
Many models contain a variable that is much faster than the others. If it
relaxes to a stable value, it can be replaced by that value, as a function of
the slow variables, and the model loses a dimension. This is the
**quasi-steady state** (QSS) approximation, and it is the singular limit of
[](../overdamped/singular-limit.md) again, now with two variables. The recipe
is to scale until a small $\ve$ multiplies the fast derivative, check that the
fast variable relaxes, reduce, and then ask what the reduction gets wrong:
usually only the first moments, and sometimes that is exactly what matters.
This page gives the method; the two pages after it apply it to an enzyme and
to the drug of [](drug-in-the-plane.md).
:::

## Two variables, two clocks

In [](../overdamped/scaling-the-equation.md) the equation of motion of a cell
was scaled to $\ve\,\dd\tilde v/\dd\tilde t = \tilde F - \tilde v$
({eq}`eq-eom-scaled`), with $\ve = \tau_i/T$ the ratio of the cell's inertial
time to the time of the question. In [](../overdamped/singular-limit.md),
setting $\ve = 0$ replaced the differential equation by the algebraic one
$\tilde v = \tilde F$: the velocity followed the force instantaneously, and the
initial velocity was lost in a layer of width $\ve$. There the slow quantity,
the force, was prescribed. Now let it have an equation of its own.

:::{prf:definition} Fast-slow system
:label: def-qss-fast-slow

A **fast-slow system** is

$$
\frac{\dd x}{\dd\tau} = f(x, y),
\qquad
\ve\,\frac{\dd y}{\dd\tau} = g(x, y),
\qquad 0 < \ve \ll 1 ,
$$ (eq-qss-fast-slow)

with $x$ the **slow** variable and $y$ the **fast** variable. Time $\tau$ is
measured on the slow scale, and $\ve$ is the ratio of the fast time scale to
the slow one.
:::

$\ve$ is never given; it is found by scaling, as in
[](../overdamped/scaling-the-equation.md). Choose units for each variable and
for time so that $f$ and $g$ are of order one; $\ve$ is whatever dimensionless
group is left in front of a derivative. If no small group appears there, there
is no separation of time scales and the QSS approximation has no
justification. Setting $\ve = 0$ gives the same objects as in
[](../overdamped/singular-limit.md).

:::{prf:definition} Critical manifold, reduced problem, layer problem
:label: def-qss-reduced

The **critical manifold** is the set where the fast equation is at rest,
$g(x, y) = 0$. Where it can be solved as $y = h(x)$, the fast variable is
**slaved** to the slow one.

The **reduced problem** (the outer problem) is the slow equation on the
critical manifold,

$$
\frac{\dd x}{\dd\tau} = f\big(x, h(x)\big) .
$$ (eq-qss-reduced)

The **layer problem** (the inner problem) is the fast equation on the fast
time $T = \tau/\ve$, with $\ve = 0$:

$$
\frac{\dd y}{\dd T} = g(x_0, y), \qquad x = x_0 \text{ frozen} .
$$ (eq-qss-layer)
:::

The reduced problem is one dimension smaller than the full one and cannot take
an initial value for $y$; that value is lost in the layer, as the initial
velocity was. The layer problem says where $y$ goes while $x$ has not had time
to move, and for the reduction to make sense $y$ must go to the critical
manifold. The layer problem is a flow on a line with a steady state at
$h(x_0)$, so by the slope test ({prf:ref}`prop-ir-slope-test`) this holds when

$$
\frac{\partial g}{\partial y}\big(x, h(x)\big) < 0 .
$$ (eq-qss-attracting)

A critical manifold with this property is **attracting**. Where the slope is
positive the manifold repels, and the reduction describes nothing that can be
observed.

| | singular limit of a cell | fast-slow system |
|---|---|---|
| equation | $\ve\,\dd\tilde v/\dd\tilde t = \tilde F(\tilde t) - \tilde v$ | $\dd x/\dd\tau = f$, $\ve\,\dd y/\dd\tau = g(x, y)$ |
| slow quantity | the force $\tilde F$, prescribed | $x$, with its own equation |
| outer solution | $\tilde v = \tilde F$ ({prf:ref}`def-outer-inner`) | $y = h(x)$ and the reduced problem |
| inner solution | on $\sigma = \tilde t/\ve$, force frozen ({prf:ref}`def-outer-inner`) | on $T = \tau/\ve$, $x$ frozen |
| attracting | $\partial(\tilde F - \tilde v)/\partial\tilde v = -1 < 0$ | {eq}`eq-qss-attracting` |
| error | at most $\ve M$ ({prf:ref}`thm-composite-uniform`) | $O(\ve)$ (Tikhonov, below) |

:::{prf:theorem} Tikhonov
:label: thm-qss-tikhonov

Let the critical manifold $y = h(x)$ be attracting, in the sense of
{eq}`eq-qss-attracting`, along a solution $\bar x(\tau)$ of the reduced
problem {eq}`eq-qss-reduced` on a finite interval $0 \le \tau \le \tau_1$, with
$\bar x(0) = x(0)$. Let $y(0)$ lie in the basin of attraction of $h(x(0))$ for
the layer problem {eq}`eq-qss-layer`. Then, as $\ve \to 0$,

$$
x(\tau) \to \bar x(\tau) \quad\text{for } 0 \le \tau \le \tau_1,
\qquad
y(\tau) \to h\big(\bar x(\tau)\big) \quad\text{for } 0 < \tau_0 \le \tau \le \tau_1 ,
$$

uniformly, for every fixed $\tau_0 > 0$. For smooth $f$ and $g$ the error is
$O(\ve)$.
:::

We state this without proof {cite:p}`tikhonov1952`. It is the nonlinear
version of {prf:ref}`thm-composite-uniform`, as noted in
{prf:ref}`rmk-reading-bound`: the slow variable is right from the start, the
fast variable after the layer. For the drug model,
{prf:ref}`prop-pp-catch-up` of [](drug-in-the-plane.md) proves the part
that matters by hand.

:::{warning} Finite time
Tikhonov's theorem holds on a fixed finite interval $[0, \tau_1]$. It does
not say which steady state a solution reaches as $\tau \to \infty$. Near a
threshold a solution can hesitate for a long time before it decides, and that
long-time question needs a separate argument. The drug example of
[](fast-drug.md) is such a case.
:::

:::{important} The QSS recipe
1. **Scale.** Find units in which a small $\ve$ multiplies one derivative.
   Write $\ve$ in the original parameters and estimate it.
2. **Slave.** Solve $g(x, y) = 0$ for $y = h(x)$, the critical manifold.
3. **Check.** Verify $\partial g/\partial y < 0$ on it. Without this the
   reduction is wrong, not just inaccurate.
4. **Reduce.** Substitute $y = h(x)$ into the slow equation.
5. **Ask what was lost.** The initial value of $y$, and $O(\ve)$ corrections.
   Decide whether the question at hand can see them.
:::

## Two examples

- [](michaelis-menten.md): an enzyme converts its substrate. The small
  parameter is the ratio of enzyme to substrate, and the reduction gives the
  Michaelis-Menten law. Step 5 finds that only the first instant is lost.
- [](fast-drug.md): the drug and the bacteria of [](drug-in-the-plane.md). The
  small parameter is the growth rate of the cells over the clearance rate of
  the drug. Here step 5 matters: the steady states are exact, the dynamics are
  slowed slightly, and the threshold inoculum depends on the initial drug
  level, which the reduced model cannot express.

In both examples the critical manifold is attracting everywhere, so the recipe works along every solution. When the
critical manifold folds back on itself, condition {eq}`eq-qss-attracting`
fails at the folds, and the fast variable jumps. That is the situation of
[](fitzhugh-nagumo.md).
