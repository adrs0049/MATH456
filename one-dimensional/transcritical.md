---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/transcritical.pdf
    id: one-dimensional-transcritical-pdf
downloads:
  - id: one-dimensional-transcritical-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Transcritical Bifurcation: Turning the Dose

:::{tip} Big Idea
Some parameters are in our hands: an experimenter chooses the dose of a drug.
As a parameter changes, the steady states move, and at special values they
appear, disappear, or exchange stability. Such a change is a **bifurcation**,
and at a bifurcation the biology changes its outcome. Turning the dose in the
shared-good model and drawing all its phase lines in one picture, the
**bifurcation diagram**, reveals a **transcritical** bifurcation: it sets the
smallest dose at which resistance persists, and near it the culture takes a
long time to settle.
:::

## The shared-good model at three doses

Return to resistance as a shared good ({prf:ref}`def-pf-shared` in
[](populations-and-fractions.md)):

$$
\dot p = f(p; a) = p\,(1-p)\,\big[\,(a - c) - a\,p\,\big] .
$$ (eq-bf-shared)

It has two parameters. The cost $c$ is set by the cell. The dose, through the
kill rate $a = kD_0$, is set by us. The semicolon in $f(p; a)$ separates the
variable $p$, which changes in time, from the parameter $a$, which we hold
fixed during an experiment and change between experiments. Hold the cost
fixed and turn the dose. There are three cases, according to how the dose
compares with the cost.

- **Below the cost, $a < c$.** The bracket $(a - c) - ap$ is negative for
  every $p$ in $[0, 1]$, so $f < 0$ and every arrow points left. The formula
  $p^\ast = 1 - c/a$ gives a negative number, not a fraction. The resistant
  cells are lost: the dose is too small to pay for the enzyme.
- **At the cost, $a = c$.** The bracket is $-ap$, and
  $f(p) = -c\,p^2(1-p) < 0$ on $(0, 1)$. The arrows still point left, and
  $p = 0$ is still stable, but $f'(0) = 0$: the slope test says nothing.
- **Above the cost, $a > c$.** As in {prf:ref}`ex-pf-shared`, a stable mixture
  $p^\ast = 1 - c/a$ lies inside $(0, 1)$, and both ends are unstable.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_INNER, C_COMP = "#111111", "#D55E00", "#009E73"
C_GREY = "#999999"


def arrow(ax, xy0, xy1, color=C_INNER, lw=1.5, ms=11):
    ax.add_patch(FancyArrowPatch(xy0, xy1, arrowstyle="-|>", mutation_scale=ms,
                                 color=color, lw=lw, zorder=3))


def flow_arrows(ax, f, zeros, lo, hi, at=None, gap=0.03, **kw):
    """Arrows of x' = f(x) between the steady states on [lo, hi]: horizontal
    along y = 0, or vertical at x = at (a phase line inside a bifurcation
    diagram)."""
    pts = sorted(set([lo, hi, *zeros]))
    g = gap * (hi - lo)
    for x0, x1 in zip(pts[:-1], pts[1:]):
        if x1 - x0 < 3 * g:
            continue
        u, v = x0 + g, x1 - g
        if f(0.5 * (x0 + x1)) < 0:
            u, v = v, u
        if at is None:
            arrow(ax, (u, 0), (v, 0), **kw)
        else:
            arrow(ax, (at, u), (at, v), **kw)


def mark_states(ax, f, zeros, lo, hi, ms=7):
    """Filled: stable. Open: unstable. Only the sides inside [lo, hi] count."""
    eps = 1e-4 * (hi - lo)
    for z in zeros:
        inward = []
        if z - eps > lo:
            inward.append(f(z - eps) > 0)
        if z + eps < hi:
            inward.append(f(z + eps) < 0)
        if all(inward):
            ax.plot([z], [0], "o", color=C_COMP, ms=ms, zorder=4)
        else:
            ax.plot([z], [0], "o", mfc="white", mec=C_EXACT, ms=ms, zorder=4)


def phase_line(ax, f, zeros, lo, hi):
    ax.plot([lo, hi], [0, 0], color=C_EXACT, lw=1.0)
    flow_arrows(ax, f, zeros, lo, hi)
    mark_states(ax, f, zeros, lo, hi)
    ax.set_ylim(-1, 1); ax.set_yticks([])
    ax.spines["left"].set_visible(False); ax.spines["bottom"].set_visible(False)


c = 0.1
f_shared = lambda p, a: p * (1 - p) * ((a - c) - a * p)
shared_zeros = lambda a: [0, 1] + ([1 - c / a] if a > c else [])

fig, axs = plt.subplots(2, 3, figsize=(7.2, 3.1), sharex=True, sharey="row",
                        gridspec_kw={"height_ratios": [3, 1]})
pp = np.linspace(0, 1, 400)
for j, (a, title) in enumerate([(0.05, r"$a = 0.05 < c$"),
                                (0.1, r"$a = c = 0.1$"),
                                (0.5, r"$a = 0.5 > c$")]):
    f = lambda p, a=a: f_shared(p, a)
    ax = axs[0, j]
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.plot(pp, f(pp), color=C_EXACT, lw=1.8)
    ax.set_title(title)
    phase_line(axs[1, j], f, shared_zeros(a), 0, 1)
    axs[1, j].set_xlabel("resistant fraction $p$")
axs[0, 0].set_ylabel(r"$\dot p$ (1/h)")
axs[0, 0].set_ylim(-0.025, 0.06)
fig.tight_layout();
```

The shared-good model {eq}`eq-bf-shared` at three doses, with cost
$c = 0.1$/h as in {prf:ref}`ex-pf-resistance`; all figures on this page use
this cost. Top: the right-hand side. Bottom: the phase line, with stable steady states
filled and unstable ones open. Below the cost (left) resistance is lost. At the
cost (middle) it is still lost, but the graph touches zero at $p = 0$ with zero
slope. Above the cost (right) a stable mixture has appeared.

Between the first and the third phase line, the stability moved from $p = 0$
to a new steady state that was not there before. The rest of this page finds
where exactly this happens, and how.

## The bifurcation diagram

Instead of drawing phase lines one dose at a time, draw them all in one
picture. Turn each phase line so that it stands vertically, and place it at its
value of $a$. The steady states of neighbouring phase lines join into curves.

:::{prf:definition} Bifurcation and bifurcation diagram
:label: def-bf-bifurcation

Let $\dot x = f(x; \mu)$ depend on a parameter $\mu$. The **bifurcation
diagram** plots the steady states $x^\ast$ against $\mu$, stable states as
solid curves and unstable states as dashed curves. A **bifurcation** happens at
a value $\mu_c$ where the phase line changes its picture: steady states appear,
disappear, or change stability. $\mu_c$ is the **bifurcation value**.
:::

The steady states of {eq}`eq-bf-shared` lie on three curves,

$$
p = 0, \qquad p = 1, \qquad p^\ast(a) = 1 - \frac{c}{a},
$$ (eq-bf-branches)

and the slopes there are

$$
f'(0) = a - c, \qquad
f'(1) = c, \qquad
f'(p^\ast) = -a\,p^\ast(1 - p^\ast) = -\frac{c\,(a - c)}{a} .
$$ (eq-bf-slopes)

Reading the signs: $p = 1$ is always unstable. $p = 0$ is stable for $a < c$
and unstable for $a > c$. $p^\ast$ is the other way round: unstable for
$a < c$, where it is negative and not a fraction, and stable for $a > c$. Both
change at $a = c$, which is exactly where $p^\ast$ passes through $0$.

### The break-even dose

We found $p^\ast$ by solving the bracket for the fraction. Solve it for the
dose instead: the bracket $(a - c) - ap = a(1-p) - c$ vanishes when

$$
a = H(p) = \frac{c}{1 - p} .
$$ (eq-bf-breakeven)

$H(p)$ is the **break-even dose**: the dose at which resistant cells at
fraction $p$ exactly pay for their enzyme. It rises with $p$: the more
resistant cells there are, the less drug they leave, so the larger the dose
must be for resistance still to pay.

:::{prf:remark} The graphical method
:label: rem-bf-graphical
:class: dropdown

Factoring $1 - p$ out of the bracket,

$$
f(p; a) = p\,(1-p)^2\,\big[\,a - H(p)\,\big], \qquad 0 \le p < 1,
$$

so for $0 < p < 1$ the resistant fraction grows exactly when the dose lies
above the break-even dose. Draw the graph of $H$ and a horizontal line at
height $a$: $p$ grows wherever the line lies above the graph. Swap the axes,
and the graph of $H$ becomes the branch $p^\ast = 1 - c/a$, while the line at
$a$ becomes the vertical phase line.
:::

```{code-cell} python
:tags: [hide-input]
a = 0.3
ps = 1 - c / a
H = lambda p: c / (1 - p)
P = np.linspace(-0.6, 1 - c / 0.62, 400)
grows = (P >= 0) & (H(P) < a)
states = ((0, False), (ps, True), (1, False))
dot = lambda stable: (dict(color=C_COMP) if stable
                      else dict(mfc="white", mec=C_EXACT))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))

ax = axs[0]
ax.axvspan(-0.6, 0, color="0.93", zorder=0)
ax.plot(P, H(P), color=C_EXACT, lw=1.8)
for x in (0, 1):
    ax.axvline(x, color=C_EXACT, lw=1.2)
ax.axhline(a, color=C_INNER, lw=1.0)
ax.fill_between(P, H(P), a, where=grows, color=C_COMP, alpha=0.18, lw=0)
for x0, x1 in ((0, ps), (ps, 1)):
    u, v = x0 + 0.04, x1 - 0.04
    if f_shared((x0 + x1) / 2, a) < 0:
        u, v = v, u
    arrow(ax, (u, a), (v, a), ms=10)
for z, stable in states:
    ax.plot([z], [a], "o", ms=6, zorder=4, **dot(stable))
ax.text(0.7, 0.2, "$H(p)$", fontsize=9)
ax.text(-0.57, 0.03, "not a fraction", fontsize=7.5, color="#777777")
ax.set_xlim(-0.6, 1.1); ax.set_ylim(0, 0.6)
ax.set_xlabel("resistant fraction $p$"); ax.set_ylabel("$H(p)$ (1/h)")
ax.set_title("graphical method: $a$ is a horizontal line")

ax = axs[1]
ax.axhspan(-0.6, 0, color="0.93", zorder=0)
ax.plot(H(P), P, color=C_EXACT, lw=1.8)
for y in (0, 1):
    ax.axhline(y, color=C_EXACT, lw=1.2)
ax.axvline(a, color=C_INNER, lw=1.0)
ax.fill_betweenx(P, H(P), a, where=grows, color=C_COMP, alpha=0.18, lw=0)
flow_arrows(ax, lambda p: f_shared(p, a), [0, ps, 1], 0, 1, at=a, ms=10)
for z, stable in states:
    ax.plot([a], [z], "o", ms=6, zorder=4, **dot(stable))
ax.text(0.42, 0.62, r"$p^\ast = 1 - c/a$", fontsize=8)
ax.text(0.33, -0.3, "$p < 0$: not a fraction", fontsize=7.5, color="#777777")
ax.set_xlim(0, 0.6); ax.set_ylim(-0.6, 1.1)
ax.set_xlabel("dose $a$ (1/h)"); ax.set_ylabel("$p$")
ax.set_title("axes swapped: $a$ is a vertical phase line")
fig.tight_layout();
```

The same picture twice, for cost $c = 0.1$/h and dose $a = 0.3$/h. Left: the
break-even dose $H(p)$ and the line at height $a$. Where the line lies above
$H$ (shaded) the resistant fraction grows; the arrows of the phase line are
drawn along the line, with filled circles for stable steady states and open
ones for unstable. The vertical lines $p = 0$ and $p = 1$ are steady at every
dose. Right: the axes swapped. The graph of $H$ is now the branch
$p^\ast = 1 - c/a$, the line at $a$ stands up as the phase line at that dose,
and $p = 0$ and $p = 1$ lie flat. Repeating the right panel for every dose
gives the left diagram below.

The picture also gives the stability. At $p^\ast$ the line crosses a *rising*
graph, with growth to its left and decline to its right, so $p^\ast$ is stable
wherever it is a fraction. In formulas,
$f'(p^\ast) = -p^\ast(1-p^\ast)^2\,H'(p^\ast)$, and $H' > 0$.

```{code-cell} python
:tags: [hide-input]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
A = np.linspace(0, 0.6, 400)
below, above = A <= c, A >= c
for ax in axs:
    ax.axhspan(-0.6, 0, color="0.93", zorder=0)
    ax.set_xlim(0, 0.6); ax.set_ylim(-0.6, 1.1)
    ax.set_xlabel("dose $a$ (1/h)")
    for a in (0.05, 0.3, 0.5):
        ax.plot([a, a], [0, 1], color="0.85", lw=0.8, zorder=1)

ax = axs[0]
ax.plot(A[below], 0 * A[below], color=C_EXACT, lw=2)
ax.plot(A[above], 0 * A[above], color=C_EXACT, lw=1.2, ls="--")
ax.plot(A, 1 + 0 * A, color=C_EXACT, lw=1.2, ls="--")
ax.plot(A[above], 1 - c / A[above], color=C_EXACT, lw=2)
A_neg = np.linspace(c / 1.6, c, 100)
ax.plot(A_neg, 1 - c / A_neg, color=C_EXACT, lw=1.2, ls="--")
for a in (0.05, 0.3, 0.5):
    flow_arrows(ax, lambda p, a=a: f_shared(p, a), shared_zeros(a), 0, 1,
                at=a, ms=9, lw=1.2)
ax.plot([c], [0], "o", color=C_INNER, ms=5, zorder=5)
ax.annotate("transcritical\nat $a = c$", xy=(c, 0), xytext=(0.13, -0.45),
            fontsize=8, color=C_INNER)
ax.text(0.33, -0.3, "$p < 0$: not a fraction", fontsize=7.5, color="#777777")
ax.text(0.36, 0.55, r"$p^\ast = 1 - c/a$", fontsize=8)
ax.set_ylabel("steady states $p$")
ax.set_title(r"shared good: $p^\ast$ rises from $0$")

ax = axs[1]
ax.plot(A[below], 0 * A[below], color=C_EXACT, lw=2)
ax.plot(A[above], 0 * A[above], color=C_EXACT, lw=1.2, ls="--")
ax.plot(A[below], 1 + 0 * A[below], color=C_EXACT, lw=1.2, ls="--")
ax.plot(A[above], 1 + 0 * A[above], color=C_EXACT, lw=2)
ax.plot([c, c], [0, 1], color=C_INNER, lw=2, alpha=0.6)
for a in (0.05, 0.3, 0.5):
    flow_arrows(ax, lambda p, a=a: (a - c) * p * (1 - p), [0, 1], 0, 1,
                at=a, ms=9, lw=1.2)
ax.annotate("every $p$ is\nsteady at $a = c$", xy=(c, 0.5), xytext=(0.12, 0.4),
            fontsize=8, color=C_INNER)
ax.set_title("no sharing: the outcome jumps")
fig.tight_layout();
```

Bifurcation diagrams in the dose $a$, for cost $c = 0.1$/h. Solid: stable.
Dashed: unstable. The vertical arrows are the phase lines at $a = 0.05$, $0.3$
and $0.5$/h. Left: the shared-good model {eq}`eq-bf-shared`. The branch
$p^\ast = 1 - c/a$ comes up through $p = 0$ at $a = c$, and the two exchange
stability. Right: the model without sharing, $\dot p = (a - c)\,p(1-p)$. The
outcome jumps from $p = 0$ to $p = 1$ as the dose passes $c$.

## The transcritical bifurcation

:::{prf:proposition} A transcritical bifurcation at the cost
:label: prop-bf-transcritical

Fix $c > 0$ in {eq}`eq-bf-shared`. As the dose $a$ increases through $c$, the
steady state $p^\ast = 1 - c/a$ passes through $p = 0$, and the two exchange
stability. For $a < c$, $p = 0$ is stable and $p^\ast < 0$ is unstable. For
$a > c$, $p = 0$ is unstable and $p^\ast \in (0, 1)$ is stable. The steady
resistant fraction therefore grows continuously from $0$, and just above the
bifurcation value

$$
p^\ast = \frac{a - c}{a} \approx \frac{a - c}{c} .
$$
:::

:::{prf:proof}
:class: dropdown

$p^\ast = 1 - c/a$ is negative for $a < c$, zero at $a = c$, and in $(0, 1)$
for $a > c$. The stability follows from the signs of the slopes
{eq}`eq-bf-slopes` and the slope test ({prf:ref}`prop-ir-slope-test`):
$f'(0) = a - c$ changes sign from negative to positive at $a = c$, and
$f'(p^\ast) = -c(a - c)/a$ changes sign from positive to negative. The
approximation uses $a \approx c$ in the denominator.
:::

This is called a **transcritical** bifurcation: two branches of steady states
cross at the critical value, and each passes its stability to the other.

As biology: below $a = c$, the drug kills sensitive cells more slowly than the
enzyme slows resistant ones, and resistance is lost from any mixture. Above
$a = c$, resistance persists, at a level that starts at zero and rises with the
dose. At no dose do the sensitive cells vanish.

### Without sharing

Without sharing, for example when resistant cells pump the drug out of
themselves without lowering the level in the culture, the difference of rates
is the constant $a - c$, and $\dot p = (a - c)\,p(1-p)$. The right panel above
shows its bifurcation diagram. Below $a = c$ the resistant cells are lost;
above it they take over completely. At exactly $a = c$ the right-hand side is
zero for every $p$, so every mixture is a steady state.

:::{warning} A jump is fragile
A whole segment of steady states at a single parameter value needs the rates
of the two types to be *exactly* equal for every $p$. The smallest dependence of
the rates on $p$ breaks it. With a little sharing, $D = D_0(1 - \ve p)$ with
small $\ve > 0$, the transcritical bifurcation returns, and the interior steady
state crosses $(0, 1)$ over the window of doses $c < a < c/(1 - \ve)$, which
shrinks to a point as $\ve \to 0$.
:::

### Resistance at low doses

:::{note} Biology background: the MIC
The **minimum inhibitory concentration** (MIC) of an antibiotic is the
smallest dose that stops the growth of a sensitive strain. Doctors choose doses
above it. In our model sensitive cells grow at $r - a$, so the MIC is the dose
with $a = r$.
:::

Resistance persists above $a = c$. The cost of resistance is small compared
with the growth rate, $c \ll r$, so this threshold lies far below the MIC:
doses too small to treat anything still keep resistant cells in the
population. The threshold $a = c$ is the same with or without sharing; sharing
decides only what happens above it. Gullberg and colleagues competed resistant and sensitive bacteria at
very low drug concentrations and found that resistant cells were enriched at
concentrations up to several hundred times below the MIC
{cite:p}`gullberg2011`.

## Critical slowing down

At $a = c$ the slope at $p = 0$ is zero, and near $p = 0$

$$
\dot p = -c\,p^2(1 - p) \approx -c\,p^2 .
$$ (eq-bf-at-critical)

The phase line says $p = 0$ is stable, but the decay is slow: separating
variables in $\dot p = -c\,p^2$ gives $p(t) = p_0/(1 + c\,p_0\,t)$, a power of
$t$ instead of an exponential.

Near the bifurcation value the approach is exponential again, but slow. The
time scale is $1/|f'|$ at the stable steady state, from {eq}`eq-bf-slopes`:

$$
\frac{1}{c - a} \quad (a < c), \qquad
\frac{a}{c\,(a - c)} \quad (a > c),
$$

and both grow without bound as $a \to c$.

:::{prf:remark} Critical slowing down
:label: rem-bf-critical-slowing
:class: dropdown

Near a bifurcation the system takes a long time to settle, because the arrows
on the phase line near the steady state are short. This is **critical slowing
down**. In an experiment, a culture near $a = c$ that is followed for a day has
not reached its steady state, and the fraction measured at the end is not
$p^\ast$.
:::

[](saddle-node.md) turns the same dose in a model that counts cells instead of
fractions, and finds a second kind of bifurcation.
