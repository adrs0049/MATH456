---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/singular-limit.pdf
    id: overdamped-singular-limit-pdf
downloads:
  - id: overdamped-singular-limit-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# The Singular Limit

:::{tip} Big Idea
Setting $\ve = 0$ in $\ve\,\dd\tilde v/\dd\tilde t = \tilde F - \tilde v$
deletes the derivative, and with it the room for the initial condition. The
condition is not lost: it lives in a window one $\tau_i$ wide, which the time
scale of our question cannot resolve but the cell's own time scale resolves in
full. One equation, one solution, two time scales, and a rule (matching) for
joining the two descriptions.
:::

## Set $\ve = 0$

The full problem is an initial value problem, one first-order equation and one
condition to fix the one constant it produces:

$$
\ve\,\frac{\dd\tilde v}{\dd\tilde t} = \tilde F - \tilde v,
\qquad\qquad
\tilde v(0) = \tilde v_0 .
$$ (eq-ivp)

Setting $\ve = 0$ deletes the derivative:

$$
\tilde v = \tilde F,
\qquad\text{that is,}\qquad
\gamma v = F .
$$ (eq-outer)

| | order | arbitrary constants | conditions to impose | balanced |
|---|---|---|---|---|
| full, $\ve \ne 0$ | $1$ | $1$ | $1$ | yes |
| reduced, $\ve = 0$ | $0$ | $0$ | $1$ | no |

At $\tilde t = 0$ the reduced equation insists that $\tilde v(0) = \tilde F(0)$.
If the cell started at any other speed, the reduced model does not approximate
the initial condition badly; it contradicts it.

:::{prf:definition} Singular limit
:label: def-singular-limit

Let a model depend on a small parameter $\ve$, and let the **reduced problem**
be what is left after setting $\ve = 0$. The limit is **singular** if setting
$\ve = 0$ lowers the order, so the reduced problem cannot carry all of the
data. The reduced solution then fails in a thin region called a **boundary
layer** (an **initial layer** when the thin variable is time), whose thickness
is set by $\ve$, and a second time scale is needed to describe it. A limit
that keeps the order is **regular**: the reduced solution is then uniformly
close to the true one, and no second time scale is needed.
:::

Our equation of motion is singular: the order went from $1$ to $0$. Yet the
reduced equation came from a correct one by deleting a term with prefactor
$10^{-9}$, so it can only be wrong on a set small enough that nine orders of
magnitude went unnoticed. Where that set is, and how small, follows from the
two time scales in the equation.

## Two time scales

$\ve = \tau_i/T$ is a ratio of two times. $\tau_i = m/\gamma$ belongs to the
cell; $T = L/v_0$ belongs to the question, and writing the equation in units of
$T$, as in {eq}`eq-ivp`, is what compressed the layer to a point. Measuring
time instead in the cell's own unit, $\sigma = t/\tau_i = \tilde t/\ve$, the
chain rule gives $\dd/\dd\tilde t = \ve^{-1}\,\dd/\dd\sigma$, and {eq}`eq-ivp`
becomes

$$
\frac{\dd\tilde v}{\dd\sigma} = \tilde F - \tilde v,
\qquad \tilde v(0) = \tilde v_0 ,
$$ (eq-inner-equation)

with no $\ve$ in it. This is the unit the equation picks for itself.

{eq}`eq-ivp` and {eq}`eq-inner-equation` are the same equation in two units of
time, and they differ only in what they resolve. On $\sigma$ the transient
occupies an interval of order one and is resolved in full, while the migration
has not begun. On $\tilde t$ the transient is compressed into a span $\ve$ of
the axis and vanishes as $\ve \to 0$. The shape of the transient does not
depend on $\ve$ at all; only the fraction of the $\tilde t$ axis it occupies
does.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER, C_INNER, C_COMP, C_LAYER = "#111111", "#0072B2", "#D55E00", "#009E73", "#F0C808"
V0, FBAR = 0.0, 1.0            # start at rest; force balance value F/gamma = 1

EPS = 0.005                      # exaggerated tau_i / T; a real cell has 1e-9
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
sig = np.linspace(0, 6, 600)
ax[0].plot(sig, FBAR + (V0 - FBAR) * np.exp(-sig), color=C_EXACT, lw=1.8)
ax[0].axhline(FBAR, color=C_OUTER, ls="--", lw=1.4)
ax[0].set_xlabel(r"time of the cell  $\sigma = t/\tau_i$"); ax[0].set_ylabel(r"$\tilde v$")
ax[0].set_title("the layer, resolved"); ax[0].set_ylim(-0.08, 1.18)
ax[0].annotate("inner problem\n" + r"$d\tilde v/d\sigma = \tilde F - \tilde v$",
               xy=(2.2, 0.55), fontsize=8, color=C_EXACT)
ax[0].annotate(r"$\tilde v = \tilde F$", xy=(4.6, 1.05), fontsize=8, color=C_OUTER)
t = np.linspace(0, 1, 20000)
ax[1].plot(t, FBAR + (V0 - FBAR) * np.exp(-t / EPS), color=C_EXACT, lw=1.8, label="exact")
ax[1].axhline(FBAR, color=C_OUTER, ls="--", lw=1.4, label=r"outer, $\tilde v = \tilde F$")
ax[1].axvspan(0, EPS, color=C_LAYER, alpha=0.9, lw=0)
ax[1].set_xlabel(r"time of the question  $\tilde t = t/T$")
ax[1].set_title("the layer is invisible"); ax[1].set_ylim(-0.08, 1.18)
ax[1].annotate("the whole of the left panel\nis this yellow sliver", xy=(EPS, 0.5),
               xytext=(0.22, 0.38), fontsize=8, color="#8a6d00",
               arrowprops=dict(arrowstyle="->", color="#8a6d00", lw=1.0))
ax[1].legend(loc="lower right", frameon=False)
fig.tight_layout();
```

One solution on two time axes, starting from rest with $\tilde F = 1$ and an
exaggerated $\ve = 0.005$. **Left:** on the cell's time $\sigma$ the transient
fills the frame. **Right:** on the time of the question $\tilde t$ the same
curve is a jump at the origin followed by the outer solution
$\tilde v = \tilde F$; the whole left panel lives inside the yellow sliver of
width $\ve$. For a real cell the sliver is a million times thinner than for
$\ve = 10^{-3}$.

Neither view is wrong; each answers a different question. Ours is where the
cell is after a body length, and on that time scale $\tilde v = \tilde F$ from
the start. Whether the ignored microsecond matters for that question is settled
in [](force-sets-velocity.md).

## Inner and outer solutions

The equation is linear, so it can be solved exactly and every approximation
checked against the exact solution. For constant $\tilde F$,

$$
\tilde v(\tilde t) = \tilde F + C\,e^{-\tilde t/\ve},
\qquad
C = \tilde v_0 - \tilde F \ \text{ fixed by the initial condition,}
$$ (eq-general-solution)

by the integrating factor ({prf:ref}`prop-integrating-factor`; the steps are in
[](integrating-factor.md)). Letting $\ve \to 0$ with $\tilde t > 0$ fixed,
$e^{-\tilde t/\ve} \to 0$, so the second term disappears and with it every
trace of $C$.

```{code-cell} python
:tags: [hide-input]
fig, ax = plt.subplots(figsize=(5.4, 2.7))
t = np.linspace(0, 1, 4000)
for eps, col in ((0.2, "#9ecae1"), (0.05, "#4292c6"), (0.01, C_OUTER)):
    ax.plot(t, np.exp(-t / eps), color=col, lw=1.9, label=rf"$\varepsilon = {eps}$")
ax.plot([0], [1], "o", color=C_INNER, ms=6, zorder=5)
ax.annotate("every curve starts at $1$,\nfor every $\\varepsilon$", xy=(0.0, 1.0),
            xytext=(0.17, 0.86), fontsize=8, color=C_INNER,
            arrowprops=dict(arrowstyle="->", color=C_INNER, lw=1.0))
ax.annotate("but for any fixed $\\tilde t > 0$\nit tends to $0$", xy=(0.45, 0.0),
            xytext=(0.45, 0.30), fontsize=8, color="#555555",
            arrowprops=dict(arrowstyle="->", color="#777777", lw=1.0))
ax.set_xlabel(r"time of the question  $\tilde t$"); ax.set_ylabel(r"$e^{-\tilde t/\varepsilon}$")
ax.set_xlim(-0.02, 1.0); ax.set_ylim(-0.03, 1.08)
ax.legend(loc="upper right", frameon=False)
ax.set_title("the term that carried the free constant");
```

The function $e^{-\tilde t/\ve}$ that multiplies the free constant $C$, for
three values of $\ve$. Every curve equals $1$ at $\tilde t = 0$, but at any
fixed $\tilde t > 0$ it tends to $0$ as $\ve \to 0$.

:::{prf:remark} The constant did not vanish; its carrier did
:label: rmk-carrier

Nothing removed $C$ from the problem. What the limit removed is the function it
multiplies: $e^{-\tilde t/\ve} \to 0$ for every $\tilde t > 0$, so $C$
multiplies nothing and becomes invisible. Only at $\tilde t = 0$ does the
exponential equal $1$, for every $\ve$. The initial condition is compressed
into a window of width $\ve$, which the time scale of the question cannot
resolve and the cell's time scale $\sigma$ resolves in full.
:::

The exact solution therefore has two different limits as $\ve \to 0$,
depending on what is held fixed while $\ve$ shrinks.

:::{prf:definition} Outer solution
:label: def-outer

Let a problem with a small parameter $\ve$ have a thin initial layer. Its
**outer solution** $u_\mathrm{out}$ is obtained by setting $\ve = 0$ in the
equation written on the time scale of the question. It has no free constant,
follows the forcing, and is valid for times long compared with the layer.
:::

:::{prf:example} The outer solution
:label: ex-outer-solution

Setting $\ve = 0$ in {eq}`eq-ivp` leaves $0 = \tilde F - \tilde v$, so

$$
\tilde v_\mathrm{out}(\tilde t) = \tilde F .
$$ (eq-toy-outer)

There is no free constant, and $\tilde v_\mathrm{out}(0) = \tilde v_0$ only if
$\tilde v_0 = \tilde F$. This is the limit of the exact solution
{eq}`eq-general-solution` as $\ve \to 0$ with $\tilde t > 0$ held fixed.
:::

:::{prf:definition} Inner solution
:label: def-inner

The **inner solution** $u_\mathrm{in}$ is obtained by first rescaling time by
the thickness of the layer, so that $\ve$ leaves the equation, and then
setting $\ve = 0$. It keeps the order of the equation, carries the initial
condition, and is valid inside the layer.
:::

:::{prf:example} The inner solution
:label: ex-inner-solution

With $\sigma = \tilde t/\ve$ the equation is {eq}`eq-inner-equation`. It
contains no $\ve$, so setting $\ve = 0$ changes nothing here; when the force
varies it freezes the force at $\tilde F(0)$, as in
{prf:ref}`ex-ramped-force`. The equation is first order, takes the initial
condition, and has the solution

$$
\tilde v_\mathrm{in}(\sigma) = \tilde F + (\tilde v_0 - \tilde F)\,e^{-\sigma} .
$$ (eq-toy-inner)

This is the limit of {eq}`eq-general-solution` as $\ve \to 0$ with
$\sigma = \tilde t/\ve$ held fixed.
:::

The interval $0 \le \tilde t \lesssim \ve$ on which the two limits differ is
the initial layer. Its thickness $\ve$ is set by the equation, not chosen; in
dimensional time it is one $\tau_i$, a microsecond.

:::{prf:definition} Matching and the composite solution
:label: def-matching

The outer and inner solutions both hold in the **overlap region**, where the
time is long compared with the layer and short compared with the slow scale.
Requiring them to agree there,

$$
\lim_{\sigma\to\infty} u_\mathrm{in}(\sigma)
= \lim_{t\to 0^+} u_\mathrm{out}(t)
\;\equiv\; u_\mathrm{common} ,
$$ (eq-matching)

is the **matching condition**, and $u_\mathrm{common}$ is the **common
limit**. It is counted twice if the two solutions are added, so the
**composite** (or **uniform**) solution subtracts it once,

$$
u_\mathrm{unif} = u_\mathrm{in} + u_\mathrm{out} - u_\mathrm{common} .
$$ (eq-composite)
:::

Where each one is right:

| | $\tilde t \sim \ve$ | $\ve \ll \tilde t \ll 1$ | $\tilde t \sim 1$ |
|---|---|---|---|
| $\tilde v_\mathrm{in}$ | right | right, $= \tilde F(0)$ | wrong: stuck at $\tilde F(0)$ |
| $\tilde v_\mathrm{out}$ | wrong: no transient | right, $= \tilde F(0)$ | right |

In the overlap region both equal $\tilde F(0)$; if the two limits in
{eq}`eq-matching` disagreed, the layer would have been scaled wrongly. Inside
the layer $\tilde v_\mathrm{out} \approx \tilde v_\mathrm{common}$ and the
composite reduces to $\tilde v_\mathrm{in}$; outside it
$\tilde v_\mathrm{in} \approx \tilde v_\mathrm{common}$ and the composite
reduces to $\tilde v_\mathrm{out}$.

:::{prf:example} The composite solution
:label: ex-composite-solution

From {prf:ref}`ex-inner-solution`, $\tilde v_\mathrm{in} \to \tilde F$ as
$\sigma \to \infty$, and from {prf:ref}`ex-outer-solution`,
$\tilde v_\mathrm{out} = \tilde F$ for all $\tilde t$. The matching condition
holds with $\tilde v_\mathrm{common} = \tilde F$, and

$$
\tilde v_\mathrm{unif}
= \tilde F + (\tilde v_0 - \tilde F)\,e^{-\tilde t/\ve},
\qquad\text{that is}\qquad
v(t) = \frac{F}{\gamma} + \Bigl(v_0 - \frac{F}{\gamma}\Bigr)e^{-t/\tau_i} .
$$ (eq-cell-composite)
:::

:::{prf:proposition} The initial layer
:label: prop-initial-layer

For constant $F$, the solution of
$m\,\dd v/\dd t = F - \gamma v$ with $v(0) = v_0$ is

$$
v(t) = \underbrace{\frac{F}{\gamma}}_{\text{outer}}
     + \underbrace{\Bigl(v_0 - \frac{F}{\gamma}\Bigr)e^{-t/\tau_i}}_{\text{layer}} ,
$$ (eq-layer-solution)

which is exactly the composite {eq}`eq-cell-composite`.
:::

:::{prf:proof}
:class: dropdown

In standard form the equation is $\dd v/\dd t + v/\tau_i = F/m$ with
$\tau_i = m/\gamma$. By {prf:ref}`prop-integrating-factor` with
$p = 1/\tau_i$, $P(t) = t/\tau_i$ and $q = F/m$,

$$
v(t) = e^{-t/\tau_i}\Bigl(v_0 + \frac{F}{m}\int_0^t e^{s/\tau_i}\,\dd s\Bigr)
= e^{-t/\tau_i}\Bigl(v_0 + \frac{F}{\gamma}\bigl(e^{t/\tau_i} - 1\bigr)\Bigr),
$$

which is {eq}`eq-layer-solution`. With velocity unit $V$, the scaled variables
are $v = V\tilde v$, $v_0 = V\tilde v_0$, $F/\gamma = V\tilde F$ and
$t/\tau_i = \tilde t/\ve$, so dividing by $V$ gives
{eq}`eq-cell-composite`.
:::

For a constant force the composite is not an approximation: it reproduces the
exact solution, because the outer solution is itself the common limit. When the
force varies this is no longer true, and the next section measures the error.

:::{prf:remark} What to do about a singular limit
:label: rmk-singular-recipe

Singular limits are not to be avoided; almost every useful reduction in
mathematical biology is one. The method has four steps:

1. set $\ve = 0$ and get the **outer** problem;
2. rescale time and get the **inner** problem, which says how thick the layer
   is and what happens inside it;
3. **match** the two and form the composite;
4. estimate what the layer contributes to a named observable, and only then
   decide whether to keep it.

Step 4 is the one most often skipped; for the cell it is carried out in
[](force-sets-velocity.md).
:::

## A force that varies on the slow time scale

With $\tilde F$ constant, the outer solution is the common limit and the
composite collapses onto the inner solution, so the three curves coincide and
there is nothing to see. A cell's traction changes over minutes rather than
over $\tau_i$, so the interesting case is a force that varies on the slow time
scale. Take the simplest one.

:::{prf:example} A ramped force
:label: ex-ramped-force

Consider

$$
\ve\,\frac{\dd\tilde v}{\dd\tilde t} = \tilde F(\tilde t) - \tilde v,
\qquad \tilde F(\tilde t) = 1 + \tilde t,
\qquad \tilde v(0) = 0 .
$$ (eq-ramp-problem)

**Outer.** $\tilde v_\mathrm{out}(\tilde t) = \tilde F(\tilde t) = 1 + \tilde t$,
which equals $1$ at $\tilde t = 0$. The initial condition is not met.

**Inner.** With $\tilde t = \ve\sigma$ the forcing is

$$
\tilde F(\ve\sigma) = 1 + \ve\sigma = \tilde F(0) + O(\ve) ,
$$ (eq-frozen-force)

so the inner problem is

$$
\frac{\dd\tilde v}{\dd\sigma} = 1 - \tilde v + O(\ve),
\qquad \tilde v(0) = 0
\qquad\Longrightarrow\qquad
\tilde v_\mathrm{in}(\sigma) = 1 - e^{-\sigma} .
$$ (eq-ramp-inner)

To leading order the force is frozen at $\tilde F(0)$: the layer is too brief
for the ramp to act, and the inner solution never learns that the force grows.

**Match.** $\tilde v_\mathrm{in} \to 1$ as $\sigma\to\infty$ and
$\tilde v_\mathrm{out} \to 1$ as $\tilde t \to 0$, so
$\tilde v_\mathrm{common} = \tilde F(0) = 1$ and

$$
\tilde v_\mathrm{unif}
= \underbrace{(1 - e^{-\tilde t/\ve})}_{\text{inner}}
+ \underbrace{(1 + \tilde t)}_{\text{outer}}
- \underbrace{1}_{\text{common}}
= 1 + \tilde t - e^{-\tilde t/\ve} .
$$ (eq-ramp-composite)

**Check.** With the integrating factor $e^{\tilde t/\ve}$
({prf:ref}`prop-integrating-factor`) and one integration by parts
([](integrating-factor.md)),

$$
\tilde v_\mathrm{exact}
= 1 + \tilde t - \ve + (\ve - 1)\,e^{-\tilde t/\ve} ,
$$ (eq-ramp-exact)

so that

$$
\tilde v_\mathrm{exact} - \tilde v_\mathrm{unif}
= -\ve\bigl(1 - e^{-\tilde t/\ve}\bigr),
\qquad
\bigl|\tilde v_\mathrm{exact} - \tilde v_\mathrm{unif}\bigr| \le \ve
\quad\text{for all }\tilde t .
$$ (eq-ramp-error)

The error is uniformly $O(\ve)$: small everywhere, not only outside the layer.
:::

```{code-cell} python
:tags: [hide-input]
eps = 0.04
fig, ax = plt.subplots(figsize=(6.2, 3.1))
t = np.linspace(0, 0.6, 4000); sig = t / eps
exact = 1 + t - eps + (eps - 1) * np.exp(-sig)
outer = 1 + t
inner = 1 - np.exp(-sig)
comp = inner + outer - 1.0                          # common limit is F(0) = 1
ax.axvspan(4 * eps, 0.28, color=C_LAYER, alpha=0.28, lw=0)
ax.plot(t, exact, color=C_EXACT, lw=3.2, alpha=0.3, label="exact")
ax.plot(t, outer, color=C_OUTER, ls="--", lw=1.5, label=r"outer  $\tilde v_{\rm out} = \tilde F(\tilde t)$")
ax.plot(t, inner, color=C_INNER, ls=":", lw=2.0,
        label=r"inner  $\tilde v_{\rm in} = \tilde F(0) + (\tilde v_0 - \tilde F(0))e^{-\sigma}$")
ax.plot(t, comp, color=C_COMP, lw=1.3, label=r"composite  $\tilde v_{\rm in} + \tilde v_{\rm out} - \tilde v_{\rm common}$")
ax.set_xlabel(r"time of the question  $\tilde t = t/T$"); ax.set_ylabel(r"$\tilde v$")
ax.set_ylim(-0.08, 1.80)
ax.annotate("overlap region\n" + r"$\tau_i \ll t \ll T$", xy=(0.205, 0.45), fontsize=8,
            color="#8a6d00", ha="center")
ax.annotate("outer misses\nthe layer", xy=(0.004, 1.0), xytext=(0.045, 1.58), fontsize=8,
            color=C_OUTER, arrowprops=dict(arrowstyle="->", color=C_OUTER, lw=1.0))
ax.annotate("inner misses\nthe ramp", xy=(0.50, 1.0), xytext=(0.30, 1.12), fontsize=8,
            color=C_INNER, arrowprops=dict(arrowstyle="->", color=C_INNER, lw=1.0))
ax.legend(loc="lower right", frameon=False)
ax.set_title(r"a ramped force, $\tilde F(\tilde t) = 1 + \tilde t$,  with $\varepsilon = 0.04$");
```

The ramped force of {prf:ref}`ex-ramped-force` with $\ve = 0.04$. The outer
solution has the late behavior right and the initial condition wrong; the
inner solution has the initial condition right and never learns that the force
is growing; both are correct in the overlap region (shaded); and the composite
agrees with the exact solution to $O(\ve)$ everywhere.

| | What it knows about $\tilde F$ | Valid when |
|---|---|---|
| inner | only $\tilde F(0)$: the force is frozen | $t \sim \tau_i$ |
| outer | the whole of $\tilde F(\tilde t)$, instantaneously | $t \gg \tau_i$ |
| composite | both | everywhere, to $O(\ve)$ |

The reason is {eq}`eq-frozen-force`. On the cell's time scale
$\tilde t = \ve\sigma$, so
$\tilde F(\tilde t) = \tilde F(0) + \ve\sigma\tilde F'(0) + O(\ve^2)$: the
variation of the force is itself $O(\ve)$ and invisible at leading order. The
same bound holds for every force that varies on the slow time scale.

:::{prf:theorem} The composite is uniformly accurate
:label: thm-composite-uniform

Let $\tilde F$ be continuously differentiable on $[0, \infty)$ with
$|\tilde F'(\tilde t)| \le M$ for all $\tilde t \ge 0$. Let $\tilde v$ solve

$$
\ve\,\frac{\dd\tilde v}{\dd\tilde t} = \tilde F(\tilde t) - \tilde v,
\qquad \tilde v(0) = \tilde v_0 ,
$$

and let

$$
\tilde v_\mathrm{unif}(\tilde t)
= \tilde F(\tilde t) + \bigl(\tilde v_0 - \tilde F(0)\bigr)\,e^{-\tilde t/\ve}
$$ (eq-general-composite)

be the composite of the outer solution $\tilde F(\tilde t)$ and the inner
solution $\tilde F(0) + (\tilde v_0 - \tilde F(0))e^{-\sigma}$, whose common
limit is $\tilde F(0)$. Then for all $\tilde t \ge 0$,

$$
\bigl|\tilde v(\tilde t) - \tilde v_\mathrm{unif}(\tilde t)\bigr|
\le \ve M\bigl(1 - e^{-\tilde t/\ve}\bigr) \le \ve M .
$$ (eq-uniform-bound)
:::

:::{prf:proof}
:class: dropdown

Let $w = \tilde v - \tilde v_\mathrm{unif}$. At $\tilde t = 0$,
$\tilde v_\mathrm{unif}(0) = \tilde F(0) + \tilde v_0 - \tilde F(0) = \tilde v_0$,
so $w(0) = 0$.

Differentiating {eq}`eq-general-composite`,

$$
\ve\,\frac{\dd\tilde v_\mathrm{unif}}{\dd\tilde t}
= \ve\tilde F'(\tilde t) - \bigl(\tilde v_0 - \tilde F(0)\bigr)e^{-\tilde t/\ve}
= \ve\tilde F'(\tilde t) + \tilde F(\tilde t) - \tilde v_\mathrm{unif}(\tilde t) .
$$

So the composite satisfies the equation of motion up to the residual
$\ve\tilde F'$. Subtracting this from the equation for $\tilde v$,

$$
\ve\,\frac{\dd w}{\dd\tilde t} = -w - \ve\tilde F'(\tilde t),
\qquad\text{that is}\qquad
\frac{\dd w}{\dd\tilde t} + \frac{w}{\ve} = -\tilde F'(\tilde t) .
$$

By {prf:ref}`prop-integrating-factor` with $p = 1/\ve$ and $w(0) = 0$,

$$
w(\tilde t) = -\int_0^{\tilde t} e^{-(\tilde t - s)/\ve}\,\tilde F'(s)\,\dd s .
$$

Bound the integrand by $M e^{-(\tilde t - s)/\ve}$ and integrate:

$$
|w(\tilde t)| \le M\int_0^{\tilde t} e^{-(\tilde t - s)/\ve}\,\dd s
= \ve M\bigl(1 - e^{-\tilde t/\ve}\bigr) .
$$
:::

:::{prf:remark} Reading the bound
:label: rmk-reading-bound

Four consequences of {eq}`eq-uniform-bound`.

1. For a constant force $M = 0$, and the composite is exact, as in
   {prf:ref}`prop-initial-layer`.
2. For the ramp $\tilde F = 1 + \tilde t$, $M = 1$, and
   {eq}`eq-ramp-error` shows that the error is exactly
   $-\ve(1 - e^{-\tilde t/\ve})$. The bound {eq}`eq-uniform-bound` is attained,
   so it cannot be improved.
3. In dimensional variables $\ve M$ times the velocity unit is
   $\tau_i \max_t |\dd(F/\gamma)/\dd t|$: the velocity lags the force balance
   by one inertial time. A force that changes over minutes produces a lag that
   no measurement could detect.
4. The proof used that the equation is linear. For a nonlinear equation
   $\ve\,\dd\tilde v/\dd\tilde t = f(\tilde v, \tilde t)$ the same construction
   gives an $O(\ve)$ approximation provided the outer solution is a stable
   root of $f = 0$, that is $\partial f/\partial\tilde v < 0$ there
   (Tikhonov's theorem). We state this without proof.
:::
