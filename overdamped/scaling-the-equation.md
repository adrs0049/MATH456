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
In [](../modeling/why-inertia-is-irrelevant.md) inertia looked negligible
because $\tau_i/T \approx 10^{-9}$, a comparison of two numbers made outside
the equation. To delete the inertial term and ask how much error that commits,
the ratio has to be inside the equation. Measure every variable in a unit taken from the problem
instead of in seconds and metres. If the units are well chosen, every rescaled
variable is of order one, and all the information about how big a term is
moves into the constant in front of it. Four constants collapse into one
number, $\ve = \tau_i/T$, and it stands in front of the inertial term.
:::

## The question

[](../modeling/why-inertia-is-irrelevant.md) ended with an observation:
except for a few microseconds after each change in force, a cell's velocity is
$v \approx F/\gamma$. Turning that into an equation of motion means deleting
the term $m\ddot{\mathbf x}$ from Newton's second law, which leaves the
overdamped force balance $\gamma\dot{\mathbf x} = \mathbf F$. The reduced
equation cannot be exactly right: a cell that loses its traction while moving
glides to a stop under Newton's law, but stops dead under the force balance.
So the question of this part is:

> How much error do we make by using the overdamped equation in place of
> Newton's law?

The first step is to put the equation in a form where "small" means something.

## The equation

In one dimension, with $F$ the net force the cell generates for itself, Newton's
second law with Stokes drag {eq}`eq-stokes-drag` is the following.

:::{prf:definition} Equation of motion
:label: def-eom

$$
m\,\ddot x = F(t) - \gamma\,\dot x
\qquad\Longleftrightarrow\qquad
m\,\dot v = F(t) - \gamma\,v, \quad v = \dot x .
$$ (eq-eom-full)

Second order in $x$, so it needs $x(0)$ **and** $\dot x(0)$. Since $x$ does
not appear on the right it is first order in $v$, and $x$ follows by one more
integration.
:::

We work with the equation for the velocity,

$$
m\,\frac{\dd v}{\dd t} = F(t) - \gamma\,v ,
\qquad v(0) = v_0 .
$$ (eq-eom-v)

No value of $F$ is needed anywhere in this part: the whole argument is about
time scales.

## Units

| Symbol | Meaning | SI units | At cell scale |
|---|---|---|---|
| $t$ | time | $\mathrm{s}$ | s, min |
| $x$ | position | $\mathrm{m}$ | $\mu$m |
| $v = \dot x$ | velocity | $\mathrm{m\,s^{-1}}$ | $\mu$m/min |
| $a = \ddot x$ | acceleration | $\mathrm{m\,s^{-2}}$ | never measured |
| $m$ | mass | $\mathrm{kg}$ | pg |
| $F$ | force | $\mathrm{N} = \mathrm{kg\,m\,s^{-2}}$ | pN |
| $\gamma$ | drag coefficient | $\mathrm{N\,s\,m^{-1}} = \mathrm{kg\,s^{-1}}$ | pN s $\mu\mathrm{m}^{-1}$ |
| $\eta$ | viscosity | $\mathrm{Pa\,s} = \mathrm{kg\,m^{-1}s^{-1}}$ | water: $10^{-3}$ |

Every term of {eq}`eq-eom-full` is a force:

$$
\underbrace{\mathrm{kg}\cdot\mathrm{m\,s^{-2}}}_{m\ddot x}
\;=\;
\underbrace{\mathrm{N}}_{F}
\;=\;
\underbrace{\mathrm{kg\,s^{-1}}\cdot\mathrm{m\,s^{-1}}}_{\gamma\dot x}
\;=\; \mathrm{kg\,m\,s^{-2}} . \qquad\checkmark
$$ (eq-dimensional-check)

Both spellings of the units of $\gamma$ are useful. As $\mathrm{N\,s\,m^{-1}}$
it is a force per unit velocity, which is what $\gamma v = F$ says; as
$\mathrm{kg\,s^{-1}}$ the kilograms cancel against $m$ and leave seconds.

Nothing in the equation is a time except $t$ itself. The units say which
combinations of the constants, and of the cell's own length $L$ and speed
$v_0$, are times.

:::{prf:example} Which combinations are times
:label: ex-combinations-are-times
:class: dropdown

**From the equation.** The only combination of $m$ and $\gamma$ with units of
time is

$$
\frac{m}{\gamma}
= \frac{\mathrm{kg}}{\mathrm{kg\,s^{-1}}}
= \mathrm{s} ,
$$

the inertial stopping time $\tau_i$ of {prf:ref}`def-tau-i`.

**From the biology.** A length the cell has and a speed at which it moves give

$$
\frac{L}{v_0}
= \frac{\mathrm{m}}{\mathrm{m\,s^{-1}}}
= \mathrm{s} ,
$$

the migration time $T$ of {prf:ref}`def-migration-time`.

These are the two times compared in
[](../modeling/why-inertia-is-irrelevant.md), and the units alone would have
found them before anything was solved.
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
:class: dropdown

Before scaling, "inertia is negligible" was a comparison between two numbers
carried around outside the model. After scaling it is a statement about the
model itself: in {eq}`eq-eom-scaled` every tilde quantity is $O(1)$ by
construction, one term carries a prefactor of $10^{-9}$, and the other two
carry a prefactor of $1$. The small parameter is now something we can take a
limit in, and taking that limit is the reduction ([](singular-limit.md)).
:::
