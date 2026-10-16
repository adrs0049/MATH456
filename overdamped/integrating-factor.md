---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/integrating-factor.pdf
    id: overdamped-integrating-factor-pdf
downloads:
  - id: overdamped-integrating-factor-pdf
    title: Download PDF
---

# Background: The Integrating Factor

:::{tip} Big Idea
A linear first-order equation $y' + p\,y = q$ can always be solved in closed
form: multiplying by $\mu = e^{\int p}$ turns the left side into the exact
derivative $(\mu y)'$, and one integration finishes the job. For the equation
of motion the result writes the velocity as a forgotten initial condition plus
a weighted average of the recent force.
:::

Every equation solved exactly in [](index.md) is linear and first order. This
page derives the one formula that solves all of them, and then shows where each
of those exact solutions comes from.

:::{prf:proposition} Linear first-order equations
:label: prop-integrating-factor

Let $p$ and $q$ be continuous on $[0, \infty)$ and put
$P(t) = \int_0^t p(s)\,\dd s$. The initial value problem

$$
\frac{\dd y}{\dd t} + p(t)\,y = q(t),
\qquad y(0) = y_0,
$$ (eq-if-ivp)

has exactly one solution, namely

$$
y(t) = e^{-P(t)}\Bigl(y_0 + \int_0^t e^{P(s)}\,q(s)\,\dd s\Bigr) .
$$ (eq-if-solution)
:::

:::{prf:proof}
:class: dropdown

Put $\mu(t) = e^{P(t)}$, the **integrating factor**. It is positive and
satisfies $\mu' = p\,\mu$ with $\mu(0) = 1$.

*Every solution has the form {eq}`eq-if-solution`.* If $y$ solves
{eq}`eq-if-ivp`, the product rule gives

$$
\frac{\dd}{\dd t}(\mu y) = \mu y' + \mu' y = \mu\,(y' + p y) = \mu q .
$$

The left side is now an exact derivative. Integrating from $0$ to $t$ and
using $\mu(0)y(0) = y_0$,

$$
\mu(t)\,y(t) - y_0 = \int_0^t \mu(s)\,q(s)\,\dd s ,
$$

and dividing by $\mu(t) > 0$ gives {eq}`eq-if-solution`. So there is at most
one solution.

*The formula is a solution.* Differentiating {eq}`eq-if-solution` with the
product rule and the fundamental theorem of calculus,

$$
y' = -p\,e^{-P}\Bigl(y_0 + \int_0^t e^{P(s)}q(s)\,\dd s\Bigr)
+ e^{-P}e^{P}q = -p\,y + q ,
$$

and at $t = 0$ the integral vanishes, so $y(0) = y_0$.
:::

:::{prf:remark} The method in three steps
:label: rmk-if-steps
:class: dropdown

As a recipe:

1. Write the equation in the standard form {eq}`eq-if-ivp` and form
   $\mu = \exp\bigl(\int p\,\dd t\bigr)$.
2. Multiply through by $\mu$. The left side becomes $(\mu y)'$, because
   $\mu$ was built to satisfy $\mu' = p\mu$.
3. Integrate both sides and divide by $\mu$. The arbitrary constant enters
   here, once, and is fixed by the initial condition.
:::

## The equation of motion

In scaled variables the equation of motion is
$\ve\,\dd\tilde v/\dd\tilde t = \tilde F(\tilde t) - \tilde v$ with
$\tilde v(0) = \tilde v_0$. In the standard form {eq}`eq-if-ivp`,

$$
\frac{\dd\tilde v}{\dd\tilde t} + \frac{1}{\ve}\,\tilde v
= \frac{1}{\ve}\,\tilde F(\tilde t),
\qquad p = \frac{1}{\ve}, \quad P(\tilde t) = \frac{\tilde t}{\ve},
$$

and {eq}`eq-if-solution` gives, for any continuous force,

$$
\tilde v(\tilde t) = \tilde v_0\,e^{-\tilde t/\ve}
+ \frac{1}{\ve}\int_0^{\tilde t} e^{-(\tilde t - s)/\ve}\,\tilde F(s)\,\dd s .
$$ (eq-if-motion)

The first term is the initial velocity, forgotten on the time $\ve$. The
second is a weighted average of the force over the recent past, with weights
$\ve^{-1}e^{-(\tilde t - s)/\ve}$ that integrate to $1 - e^{-\tilde t/\ve}$
and are concentrated on the last few multiples of $\ve$. Every exact solution
of the equation of motion is a case of {eq}`eq-if-motion`.

**Constant force.** The integral is
$\tilde F\bigl(1 - e^{-\tilde t/\ve}\bigr)$, so

$$
\tilde v(\tilde t) = \tilde F + (\tilde v_0 - \tilde F)\,e^{-\tilde t/\ve} ,
$$

which is {eq}`eq-general-solution` of [](singular-limit.md). In dimensional variables, with
$p = 1/\tau_i$ and $q = F/m$, the same computation gives
{eq}`eq-layer-solution`.

**Ramped force.** For $\tilde F(s) = 1 + s$ and $\tilde v_0 = 0$, integrate
by parts with $u = 1 + s$ and
$\dd g = \ve^{-1}e^{-(\tilde t - s)/\ve}\,\dd s$, so that
$g = e^{-(\tilde t - s)/\ve}$:

$$
\frac{1}{\ve}\int_0^{\tilde t} e^{-(\tilde t - s)/\ve}(1 + s)\,\dd s
= \Bigl[(1 + s)\,e^{-(\tilde t - s)/\ve}\Bigr]_0^{\tilde t}
- \int_0^{\tilde t} e^{-(\tilde t - s)/\ve}\,\dd s
= 1 + \tilde t - e^{-\tilde t/\ve} - \ve\bigl(1 - e^{-\tilde t/\ve}\bigr) .
$$

This is $1 + \tilde t - \ve + (\ve - 1)e^{-\tilde t/\ve}$, which is
{eq}`eq-ramp-exact`.

**The error of the composite.** In the proof of
{prf:ref}`thm-composite-uniform` ([](singular-limit.md)) the error $w = \tilde v - \tilde v_\mathrm{unif}$
satisfies $w' + w/\ve = -\tilde F'$ with $w(0) = 0$. Here $p = 1/\ve$,
$q = -\tilde F'$ and $y_0 = 0$, and {eq}`eq-if-solution` gives

$$
w(\tilde t) = -\int_0^{\tilde t} e^{-(\tilde t - s)/\ve}\,\tilde F'(s)\,\dd s ,
$$

the formula used there.
