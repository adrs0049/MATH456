---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/scaling-the-equation.pdf
    id: overdamped-scaling-the-equation-pdf
downloads:
  - id: overdamped-scaling-the-equation-pdf
    title: Download PDF
---

# Nondimensionalization

:::{tip} Big Idea
Measure every variable in a unit taken from the problem instead of in seconds
and metres. If the units are well chosen, every rescaled variable is of order
one, and all the information about how big a term is moves into the constant
in front of it. Four constants collapse into one number, $\ve = \tau_i/T$, and
it stands in front of the inertial term.
:::

## Reference units

:::{prf:definition} Reference unit, dimensionless variable
:label: def-reference-units

A **reference unit** (or **scale**) for a variable is a fixed constant carrying
that variable's units, chosen from the problem. Writing
$$
x = X\,\tilde x
\qquad\text{equivalently}\qquad
\tilde x = \frac{x}{X}
$$
replaces $x$ by the dimensionless $\tilde x$, which counts how many $X$'s.
Dimensionless variables carry a tilde.
:::

:::{prf:definition} Order one
:label: def-order-one

A dimensionless quantity is **of order one**, written $O(1)$, if it is neither
large nor small compared with $1$: as a working rule, between about $0.1$ and
$10$. A scale is **well chosen** when the corresponding dimensionless variable,
and its derivatives, are $O(1)$ over the range we care about.
:::

There are three reasons to scale.

1. **Only dimensionless numbers can be called small.** The cell's mass is
   $10^{-12}$ kg, which sounds tiny, but that is a statement about the
   kilogram: the same mass is $10^{9}$ femtograms. "The mass is negligible" is
   not a meaningful sentence; "$\ve = 10^{-9}$" is.
2. **Parameters collapse.** Each scaled variable removes one parameter. Four
   mechanical constants become one below.
3. **Size moves into the prefactors.** If every scale is well chosen, every
   tilde variable and derivative is $O(1)$, so terms can be ranked by their
   prefactors. A term with prefactor $10^{-9}$ multiplying something $O(1)$
   contributes $10^{-9}$.

:::{warning} The catch is the "if"
A badly chosen scale leaves a tilde variable that is not $O(1)$, and reading
the prefactors then gives the wrong answer. The model usually signals this: a
term expected to be $O(1)$ comes out enormous.
:::

## Choosing the scales

Two scales are measured, one follows from them, and one is a decision.

| Variable | Scale | Value | Where it comes from |
|---|---|---|---|
| length | $L$ | $10^{-5}$ m | measured: the cell's diameter |
| velocity | $V = v_0$ | $10^{-8}$ m/s | measured: the speed cells crawl at |
| time | $T = L/V$ | $10^{3}$ s | derived: {prf:ref}`def-migration-time` |
| force | $\gamma V$ | $10^{-14}$ N | chosen: the drag on a cell moving at $V$ |

The time scale is the question's, not the equation's; the equation's own time,
$\tau_i$, will appear as the coefficient left over. The force scale is a
genuine choice, and a bet: it makes the drag on a cell moving at its observed
speed of order one, which amounts to betting that drag is the force that
matters. If the bet were wrong, a term expected to be $O(1)$ would come out
large.

## The substitution

Three replacements:

$$
t = T\,\tilde t, \qquad
v(t) = V\,\tilde v(\tilde t), \qquad
F(t) = \gamma V\,\tilde F(\tilde t) .
$$ (eq-three-replacements)

Each says: old quantity $=$ fixed unit $\times$ pure number. The units live in
$T$, $V$, $\gamma V$; the tildes are numbers. By the chain rule, with
$\tilde t(t) = t/T$,

$$
\frac{\dd v}{\dd t}
= V\,\frac{\dd\tilde v}{\dd\tilde t}\,\frac{\dd\tilde t}{\dd t}
= \frac{V}{T}\,\frac{\dd\tilde v}{\dd\tilde t} .
$$

:::{warning} Substitute in the denominator too
The derivative contains $t$ as well as $v$, so it needs both replacements.
Replacing $v$ and forgetting $\dd t$ gives $V\,\dd\tilde v/\dd\tilde t$,
which has the units of a velocity; the correct expression has units
$[V/T] = \mathrm{m/s^2}$, an acceleration. Each derivative with respect to $t$
brings one factor of $1/T$.
:::

Replace each term of {eq}`eq-eom-v`:

$$
\underbrace{m\cdot\frac{V}{T}\,\frac{\dd\tilde v}{\dd\tilde t}}_{\text{inertia}}
\;=\;
\underbrace{\gamma V\,\tilde F(\tilde t)}_{\text{applied force}}
\;-\;
\underbrace{\gamma V\,\tilde v(\tilde t)}_{\text{drag}} .
$$

Every term is a force, so divide by the force scale $\gamma V$:

$$
\underbrace{\frac{m}{\gamma T}}_{\text{dimensionless}}\,
\frac{\dd\tilde v}{\dd\tilde t} = \tilde F - \tilde v .
$$ (eq-scaled-before-naming)

The two force terms now have coefficient $1$, $V$ has dropped out, and the only
prefactor that can be large or small is $m/(\gamma T)$, in front of the
inertial term.

The initial condition must be scaled with the equation. The instant $t = 0$ is
the instant $\tilde t = 0$, so $v(0) = v_0$ becomes
$\tilde v(0) = v_0/V =: \tilde v_0$, a pure number. With $V = v_0$ it equals
$1$, but we keep the symbol, since the reduction is about what happens to it.

:::{prf:definition} The dimensionless equation of motion
:label: def-epsilon-inertia

With $t = T\tilde t$, $v = v_0\tilde v$, $F = \gamma v_0\tilde F$, the model
{eq}`eq-eom-v` together with its initial condition reads

$$
\ve\,\frac{\dd\tilde v}{\dd\tilde t} = \tilde F(\tilde t) - \tilde v,
\qquad
\tilde v(0) = \tilde v_0,
\qquad
\ve = \frac{m}{\gamma T} = \frac{\tau_i}{T} = \frac{m v_0}{\gamma L} .
$$ (eq-eom-scaled)

The single parameter $\ve$ is the ratio of the inertial time of
{prf:ref}`def-tau-i`, which belongs to the cell, to the migration time of
{prf:ref}`def-migration-time`, which belongs to the question.
For a crawling cell

$$
\ve = \frac{\tau_i}{T}
\approx \frac{10^{-6}\ \mathrm{s}}{10^{3}\ \mathrm{s}} = 10^{-9} .
$$ (eq-epsilon-number)
:::

:::{prf:remark} What the scaling bought
:label: rmk-scaling-bought

Before scaling, "inertia is negligible" was a comparison between two numbers
carried around outside the model. After scaling it is a statement about the
model itself: in {eq}`eq-eom-scaled` every tilde quantity is $O(1)$ by
construction, one term carries a prefactor of $10^{-9}$, and the other two
carry a prefactor of $1$. The small parameter is now something we can take a
limit in, and taking that limit is the reduction ([](singular-limit.md)).
:::
