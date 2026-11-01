---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/fast-drug.pdf
    id: phase-plane-fast-drug-pdf
downloads:
  - id: phase-plane-fast-drug-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# The Drug Is Fast

:::{tip} Big Idea
Meropenem is cleared in about an hour and TB divides about once a day, so the
drug is the fast variable. Following the recipe recovers the phase line of
[](../one-dimensional/saddle-node.md), and step 5 says exactly what it costs.
The steady states and the bifurcation diagram are exact for every $\ve$, the
dynamics are slowed by a few percent, and the threshold inoculum moves with
the initial drug level, which the reduced model cannot express.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_OUTER, C_INNER, C_COMP, C_GREY = "#0072B2", "#D55E00", "#009E73", "#999999"
```

The model of [](drug-in-the-plane.md), {eq}`eq-pp-model`, follows TB cells $N$
that destroy the drug meropenem, $d = D/D_0$, which kills them.
[](../one-dimensional/saddle-node.md) replaced the drug by its steady state
$\bar d(N) = 1/(1 + N/N_c)$, because the drug is cleared in about an hour
while TB needs about a day to divide, and found a bistable phase line with a
threshold inoculum $N_-$ and a saddle-node bifurcation.
[](drug-in-the-plane.md) showed geometrically why this works: the arrows are
nearly vertical, the distance to the $d$-nullcline decays at least at the
clearance rate ({prf:ref}`prop-pp-catch-up`), and the slow eigenvalue is the
slope of the phase line ({prf:ref}`prop-pp-fast-slow`). Here we follow the
recipe of [](quasi-steady-state.md) instead, and compute what the reduction
gets wrong.

**Scale.** Measure cells against the carrying capacity and time in cell
generations, $n = N/K$, $\tau = rt$, and keep $d = D/D_0$. Then
{eq}`eq-pp-model` becomes

$$
\frac{\dd n}{\dd\tau} = n(1 - n) - \alpha\,d\,n,
\qquad
\ve\,\frac{\dd d}{\dd\tau} = 1 - d - \kappa\,n\,d,
$$ (eq-qss-inoc)

with

$$
\ve = \frac{r}{\delta}, \qquad
\alpha = \frac{a}{r} = \frac{kD_0}{r}, \qquad
\kappa = \frac{K}{N_c} .
$$

The small parameter is the growth rate of the cells divided by the clearance
rate of the drug. TB grows at $r = 1$/day and meropenem's half-life of an hour
gives $\delta \approx 17$/day, so $\ve \approx 0.06$. For the *E. coli* of
[](../one-dimensional/populations-and-fractions.md), which doubles every
20 minutes, $r \approx 2$/h and $\ve \approx 3$: the cells would be faster
than the drug, and the reduction would have no justification.

**Slave.** The critical manifold is $d = h(n) = 1/(1 + \kappa n)$, which is
{eq}`eq-bf-drug-qss` in scaled form.

**Check.** $\partial g/\partial d = -(1 + \kappa n) < 0$ for every $n \ge 0$.
The manifold is attracting everywhere, and more cells make it more strongly
attracting, because they destroy drug as well as clearance removes it.

**Reduce.**

$$
\frac{\dd n}{\dd\tau} = F(n) = n(1 - n) - \frac{\alpha\,n}{1 + \kappa n},
$$ (eq-qss-inoc-reduced)

which is {eq}`eq-bf-inoculum` in scaled form.

For TB the recipe could stop here. To see what step 5 is for, treat $\ve$ as a
dial, from TB at $\ve \approx 0.06$ to *E. coli* at $\ve \approx 3$, and ask
what the reduction gets wrong. The answer is much less than the size of $\ve$
suggests.

## What does not depend on ε

:::{prf:proposition} The bifurcation diagram is exact
:label: prop-qss-inoc-exact

For every $\ve > 0$, the steady states of {eq}`eq-qss-inoc` are the points
$\big(n^\ast, h(n^\ast)\big)$ with $n^\ast$ a steady state of the reduced
equation {eq}`eq-qss-inoc-reduced`, and

- if $F'(n^\ast) > 0$, the steady state is a saddle;
- if $F'(n^\ast) < 0$, it is stable.

The bifurcation diagram of [](../one-dimensional/saddle-node.md), with its
transcritical bifurcation at $a = r$ and saddle-node at $a = a_\mathrm{sn}$,
therefore holds for the full model, whatever the value of $\ve$.
:::

:::{prf:proof}
:class: dropdown

At a steady state $\dd d/\dd\tau = 0$, which is $d = h(n)$, and then
$\dd n/\dd\tau = F(n) = 0$. Stability is {prf:ref}`prop-pp-agree` in scaled
variables. Directly: the Jacobian of {eq}`eq-qss-inoc` at $n^\ast = 0$ is
triangular with eigenvalues $1 - \alpha = F'(0)$ and $-1/\ve$. At
$n^\ast > 0$, using $1 - n^\ast = \alpha d^\ast$ and
$d^\ast = 1/(1 + \kappa n^\ast)$,

$$
\operatorname{tr} J = -n^\ast - \frac{1 + \kappa n^\ast}{\ve} < 0,
\qquad
\det J = \frac{n^\ast(1 + \kappa n^\ast)}{\ve}
\Big[1 - \frac{\alpha\kappa}{(1 + \kappa n^\ast)^2}\Big]
= -\frac{1 + \kappa n^\ast}{\ve}\,F'(n^\ast),
$$

since $F'(n^\ast) = -n^\ast\big[1 - \alpha\kappa/(1 + \kappa n^\ast)^2\big]$
by the same relations. Conclude with {prf:ref}`prop-pp-trace-det`.
:::

The steady states do not depend on the approximation, since at a steady state
the drug is at rest anyway, and neither does their stability.

## The first correction slows time

Tikhonov's theorem says the error is $O(\ve)$; here it can be computed. The
drug does not sit exactly on $h(n)$: as $n$ changes, $h(n)$ moves and the drug
chases it. This is the lag that {prf:ref}`prop-pp-catch-up` bounds by
$M/\delta$; to first order it can be found exactly.

:::{prf:proposition} The lag slows, it does not steer
:label: prop-qss-inoc-slow

The system {eq}`eq-qss-inoc` has an invariant curve
$d = h(n) + \ve h_1(n) + O(\ve^2)$ with

$$
h_1(n) = \frac{\kappa\,F(n)}{(1 + \kappa n)^3},
$$

and on it

$$
\frac{\dd n}{\dd\tau}
= F(n)\,\Big[1 - \ve\,\frac{\alpha\kappa\,n}{(1 + \kappa n)^3}\Big] + O(\ve^2) .
$$ (eq-qss-inoc-corrected)

To first order in $\ve$, the drug's lag multiplies the reduced right-hand side
by a factor between $1 - \tfrac{4}{27}\,\ve\alpha$ and $1$. The phase line,
the steady states and their stability are unchanged, and the dynamics are
slowed by at most the fraction $\tfrac{4}{27}\,\ve\alpha$.
:::

:::{prf:proof}
:class: dropdown

On a curve $d = h(n) + \ve h_1(n) + \dots$ the chain rule gives
$\ve\,\dd d/\dd\tau = \ve\,h'(n)\,\dd n/\dd\tau = \ve\,h'(n)F(n) + O(\ve^2)$.
The drug equation requires this to equal
$1 - d(1 + \kappa n) = -\ve h_1(n)(1 + \kappa n) + O(\ve^2)$, since
$h(n)(1 + \kappa n) = 1$. Matching the terms of order $\ve$,

$$
h_1(n) = -\frac{h'(n)\,F(n)}{1 + \kappa n} = \frac{\kappa\,F(n)}{(1 + \kappa n)^3},
$$

because $h'(n) = -\kappa/(1 + \kappa n)^2$. Substituting into the cell
equation,
$\dd n/\dd\tau = F(n) - \ve\,\alpha n\,h_1(n) + O(\ve^2)$, which is
{eq}`eq-qss-inoc-corrected`. Finally $x/(1 + x)^3$ has its maximum $4/27$ on
$x \ge 0$ at $x = 1/2$, so with $x = \kappa n$ the bracket lies between
$1 - \tfrac{4}{27}\ve\alpha$ and $1$.
:::

When the population grows ($F > 0$) the drug lies above its QSS value: it has
not caught up with the extra destruction. When the population shrinks it lies
below. With $\alpha = 2$ the slowdown is at most 2% for TB, and at most 15% at
$\ve = 0.5$, although that $\ve$ is no longer small.

## The layer moves the threshold

What the reduction cannot see is the initial drug level. The reduced model
starts the drug at $h(n_0)$, its QSS value. A real experiment starts it
wherever the protocol puts it:

- cells introduced into a culture or a patient that already holds the drug:
  $d(0) = 1$, more drug than QSS;
- an infusion started on an established infection: $d(0) = 0$, less.

[](drug-in-the-plane.md) showed the consequence in the plane: the separatrix
of the saddle is tilted, so the threshold inoculum depends on $d(0)$, and
for a slow drug by a factor of seven ({prf:ref}`ex-pp-protocol`). The layer
problem gives the size of the tilt.

:::{prf:proposition} The threshold inoculum depends on the protocol
:label: prop-qss-threshold-shift

To first order in $\ve$, an inoculum starting with drug level $D(0)$ is
cleared when $N(0) < N_\mathrm{th}$, where

$$
N_\mathrm{th} \approx N_-\,\Big[1 + \frac{k\,\big(D(0) - D_-\big)}{\delta + \gamma N_-}\Big]
= N_-\,\Big[1 + \frac{a\,\big(d(0) - \bar d(N_-)\big)}{\delta\,(1 + N_-/N_c)}\Big]
$$ (eq-qss-threshold-shift)

and $D_- = D_0/(1 + N_-/N_c)$ is the QSS drug level at the threshold.
:::

:::{prf:proof}
:class: dropdown

This is a first-order calculation, not a rigorous bound. Inside the layer
$n \approx n_0$ is frozen and, by {eq}`eq-qss-layer`, the drug relaxes as

$$
d(T) = h(n_0) + \big(d(0) - h(n_0)\big)\,e^{-(1 + \kappa n_0)T},
\qquad T = \tau/\ve .
$$

The cell equation is $\dd n/\dd\tau = F(n) - \alpha n\,\big(d - h(n)\big)$, and
near the threshold $F(n) \approx 0$, so the layer changes $n$ by

$$
\Delta n \approx -\alpha n_0 \int_0^\infty \big(d - h(n_0)\big)\,\dd\tau
= -\ve\,\alpha n_0\,\frac{d(0) - h(n_0)}{1 + \kappa n_0} .
$$

After the layer the solution is close to the critical manifold, where the
reduced phase line applies, so it is cleared when $n_0 + \Delta n < n_-$. To
first order this is $n_0 < n_-\big[1 + \ve\alpha\,(d(0) - h(n_-))/(1 + \kappa n_-)\big]$.
In the original units $\ve\alpha = a/\delta$, $1 + \kappa n_- = 1 + N_-/N_c$
and $h(n_-) = \bar d(N_-)$, which is the second form; the first follows from
$a\,d = kD$ and $\delta(1 + N/N_c) = \delta + \gamma N$.
:::

The correction is a product: the extra kill rate $k\,(D(0) - D_-)$, times the
time $1/(\delta + \gamma N_-)$ for which the drug is away from its QSS level.
Cells that meet a full dose before they have destroyed any of it need a larger
inoculum to survive; cells that get a head start before the drug arrives need
a smaller one. The growth rate $r$ does not appear: near the threshold growth
and killing balance, and what matters is how long the drug is out of balance
compared with how fast it kills. For TB at $a = 2$/day,
{eq}`eq-qss-threshold-shift` gives $1.27$ with no drug at the start and $1.34$
with the full dose, the thresholds read off the separatrix in
[](drug-in-the-plane.md) to two decimals.

```{code-cell} python
:tags: [hide-input]
r, Kc, Nc, a = 1.0, 10.0, 1.0, 2.0          # 1/day, 1e8 cells/mL, 1e8 cells/mL, 1/day
alpha, kappa = a / r, Kc / Nc
h = lambda n: 1 / (1 + kappa * n)
nm = (9 - np.sqrt(41)) / 20                 # the saddle n_- = N_-/K


def full(t, y, eps):
    n, d = y
    return [n * (1 - n) - alpha * d * n, (1 - d - kappa * n * d) / eps]


def separatrix(eps):
    """Stable manifold of the saddle (n_-, h(n_-)), integrated backward."""
    J = np.array([[-nm, -alpha * nm],
                  [-kappa * h(nm) / eps, -(1 + kappa * nm) / eps]])
    w, V = np.linalg.eig(J)
    v = V[:, np.argmin(w.real)].real   # numpy >= 2.5 returns complex arrays

    def leave(t, y):
        return min(y[1] + 0.05, 1.05 - y[1], y[0])
    leave.terminal = True
    back = lambda t, y: [-z for z in full(t, y, eps)]
    halves = [solve_ivp(back, (0, 50), [nm, h(nm)] + 1e-7 * s * v, events=leave,
                        rtol=1e-10, atol=1e-12, method="LSODA").y for s in (1, -1)]
    return np.hstack([halves[0][:, ::-1], halves[1]])


def threshold(eps, d0):
    n, d = separatrix(eps)
    i = np.argsort(d)
    return np.interp(d0, d[i], n[i])


fig, ax = plt.subplots(figsize=(4.8, 3.1))
E = np.linspace(0.02, 3.0, 40)
for d0, col, lab in ((1.0, C_INNER, "drug already at $D_0$"),
                     (h(nm), C_OUTER, "drug at its QSS level"),
                     (0.0, C_COMP, "no drug yet")):
    ax.plot(E, [10 * threshold(e, d0) for e in E], color=col, lw=1.6, label=lab)
    first = nm * (1 + E * alpha * (d0 - h(nm)) / (1 + kappa * nm))
    ax.plot(E, 10 * first, color=col, ls=":", lw=1.2)
for e, lab in ((1 / 16.6, "TB"), (2.9, "E. coli")):
    ax.axvline(e, color=C_GREY, lw=0.8, ls="--")
    ax.text(e + 0.04, 0.15, lab, fontsize=8, color="#555555",
            style="normal" if lab == "TB" else "italic")
ax.set_xlim(0, 3); ax.set_ylim(0, 3.6)
ax.set_xlabel(r"$\varepsilon = r/\delta$")
ax.set_ylabel(r"threshold inoculum ($10^8$ cells/mL)")
ax.legend(loc="upper left", frameon=False, fontsize=8)
ax.set_title("the threshold depends on the protocol")
fig.tight_layout();
```

The threshold inoculum of the full model {eq}`eq-qss-inoc` against $\ve$, with
$\alpha = a/r = 2$ and $\kappa = 10$, for three starting drug levels. Solid:
the full model, where the separatrix crosses the line $d = d(0)$. Dotted: the
first-order formula {eq}`eq-qss-threshold-shift`. At TB's $\ve \approx 0.06$
the three thresholds differ by a few percent; at *E. coli*'s $\ve \approx 3$
they differ by a factor of six.

:::{prf:remark} Why the QSS start gives the exact threshold
:label: rmk-qss-separatrix
:class: dropdown

The blue line in the figure is flat: a population that starts with the drug at
its QSS level has threshold exactly $N_-$, for every $\ve$. For
$n_- < n < n_+$ the critical manifold $d = h(n)$ lies below the curve
$d = (1 - n)/\alpha$ on which $\dd n/\dd\tau = 0$, and between the two curves
$n$ grows and $d$ falls. Neither curve can be crossed outward: on the first
$d$ is at rest while $n$ grows, which carries the solution above $h(n)$, and
on the second $n$ is at rest while $d$ falls. So a solution that starts on
the critical manifold with $n_- < n_0 < n_+$ stays between the curves, with
$n$ increasing, and reaches $n_+$. For $n_0 < n_-$ the two curves swap
places, $n$ falls and $d$ rises, and the solution is cleared. In other words,
the separatrix crosses the critical manifold exactly at the saddle.
:::

The answer to "what does the reduction cost?" is therefore specific. The
steady states, their stability and the bifurcation diagram cost nothing. The
dynamics are slowed, by at most 2% for TB. The threshold inoculum is shifted by
{eq}`eq-qss-threshold-shift`, which depends on an initial condition that the
reduced model cannot even express.

