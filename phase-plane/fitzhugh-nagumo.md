---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/fitzhugh-nagumo.pdf
    id: phase-plane-fitzhugh-nagumo-pdf
downloads:
  - id: phase-plane-fitzhugh-nagumo-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Pulses and Oscillations: FitzHugh–Nagumo

:::{tip} Big Idea
A drug and cells that harm each other can only switch. An **activator** that
amplifies itself and makes its own **inhibitor** can also oscillate. With the
inhibitor frozen, the activator is a bistable switch, and its bifurcation
diagram is an S with two folds. Let the inhibitor move, and it drives the
activator round the hysteresis loop of the S by itself. A trapping region and
the Poincaré–Bendixson theorem make this a closed orbit: an **oscillation**.
Move the steady state just past a fold and the system rests instead, but it is
**excitable**: it ignores small kicks and answers every large one with the
same large pulse.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER, C_INNER, C_COMP = "#111111", "#0072B2", "#D55E00", "#009E73"
C_GREY = "#999999"

c = lambda v: v - v**3 / 3            # the cubic: fast nullcline is w = c(v)


def fhn(t, y, eps, a):
    v, w = y
    return [(v - v**3 / 3 - w) / eps, v + a]


def solve(y0, T, eps, a, n=2000, **kw):
    t = np.linspace(0, T, n)
    s = solve_ivp(fhn, (0, T), y0, args=(eps, a), t_eval=t, rtol=1e-9,
                  atol=1e-11, method="LSODA", **kw)
    return s.t, s.y


def mark(ax, x, y, stable, ms=6):
    style = dict(color=C_COMP) if stable else dict(mfc="white", mec=C_EXACT)
    ax.plot([x], [y], "o", ms=ms, zorder=6, **style)


VV = np.linspace(-2.3, 2.3, 400)
mid, outer_l, outer_r = np.abs(VV) <= 1, VV <= -1, VV >= 1


def cubic(ax, lw=1.8, color=C_OUTER):
    """The fast nullcline, drawn with w across and v up: stable outer branches
    solid, unstable middle branch dashed."""
    ax.plot(c(VV[outer_l]), VV[outer_l], color=color, lw=lw)
    ax.plot(c(VV[outer_r]), VV[outer_r], color=color, lw=lw)
    ax.plot(c(VV[mid]), VV[mid], color=color, lw=lw * 0.7, ls="--")


def phase_plane(ax, F, G, xlim, ylim, trajs=(), t_max=20.0, n_grid=200,
                density=1.2, flow_color="thistle", null_colors=(C_INNER, C_OUTER),
                lw_null=2.0, lw_traj=1.0, null=True, fixed_points=True):
    """Phase portrait of x' = F(x, y), y' = G(x, y), after
    A. Buttenschoen, NumericalPDEs/de (2019): streamlines whose width grows
    with the speed, nullclines as zero contours (x-nullcline in the first
    colour), trajectories from the points in trajs (squares), and steady
    states found by fsolve from a grid of guesses and classified by the
    eigenvalues of a finite-difference Jacobian (filled: stable, open:
    unstable). F and G must accept numpy arrays. Returns the steady states."""
    X, Y = np.meshgrid(np.linspace(*xlim, n_grid), np.linspace(*ylim, n_grid))
    U, V = F(X, Y), G(X, Y)
    speed = np.hypot(U, V)
    ax.streamplot(X, Y, U, V, density=density, color=flow_color, arrowsize=1.0,
                  linewidth=0.4 + 1.6 * np.sqrt(speed / np.nanmax(speed)), zorder=1)
    if null:
        for Z, col in zip((U, V), null_colors):
            ax.contour(X, Y, Z, [0], colors=[col], linewidths=lw_null, zorder=3)
    rhs = lambda t, z: [F(*z), G(*z)]
    for z0 in trajs:
        s = solve_ivp(rhs, (0, t_max), z0, rtol=1e-9, atol=1e-11,
                      method="LSODA", dense_output=True)
        t = np.concatenate([np.linspace(0, 0.02 * t_max, 400),
                            np.linspace(0.02 * t_max, t_max, 1600)])
        ax.plot(*s.sol(t), color=C_EXACT, lw=lw_traj, zorder=4)
        ax.plot(*z0, "s", color=C_EXACT, ms=3, zorder=4, clip_on=False)
    fps = []
    if fixed_points:
        f2 = lambda z: [F(*z), G(*z)]
        for gx in np.linspace(*xlim, 12):
            for gy in np.linspace(*ylim, 12):
                z, _, ier, _ = fsolve(f2, [gx, gy], full_output=True)
                inside = (xlim[0] - 1e-9 <= z[0] <= xlim[1] + 1e-9
                          and ylim[0] - 1e-9 <= z[1] <= ylim[1] + 1e-9)
                if ier == 1 and inside and not any(np.allclose(z, p, atol=1e-6) for p in fps):
                    fps.append(z)
        e = 1e-6
        for x0, y0 in fps:
            J = np.array([[F(x0 + e, y0) - F(x0 - e, y0), F(x0, y0 + e) - F(x0, y0 - e)],
                          [G(x0 + e, y0) - G(x0 - e, y0), G(x0, y0 + e) - G(x0, y0 - e)]]) / (2 * e)
            stable = np.all(np.linalg.eigvals(J).real < 0)
            style = dict(color=C_COMP) if stable else dict(mfc="white", mec=C_EXACT)
            ax.plot([x0], [y0], "o", ms=6, zorder=6, clip_on=False, **style)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    return fps
```

(fhn-oscillate)=
## What it takes to oscillate

By {prf:ref}`prop-pp-no-spirals` in [](drug-in-the-plane.md), cross effects
of the same sign give real eigenvalues, so nothing circles. Circling needs
cross effects of *opposite* signs, $x$ promoting $y$ and $y$ suppressing $x$:
a **negative feedback loop**.

| | $\partial\dot x/\partial y$ | $\partial\dot y/\partial x$ | what it can do | example |
|---|---|---|---|---|
| mutual harm | $-$ | $-$ | switch: nodes and saddles | drug and cells |
| mutual help | $+$ | $+$ | switch: nodes and saddles | two genes that activate each other |
| activator and inhibitor | $-$ | $+$ | spirals, oscillations | FitzHugh–Nagumo |

For the spiral to wind outward the trace must also be positive
({prf:ref}`prop-pp-trace-det`), which needs a positive diagonal entry: a
variable that amplifies itself. The recipe for an oscillator is
**self-amplification plus delayed negative feedback**.

:::{note} Biology background: excitable cells
A nerve cell at rest keeps a voltage across its membrane. A small push to that
voltage dies away. A push past a **threshold** opens channels that let charge
in, which pushes the voltage further and opens more channels: a
self-amplifying runaway. The voltage spikes to the same size however far past
threshold the push was (**all or none**). More slowly, other channels let
charge out and the voltage returns to rest; for a while afterwards the cell
cannot fire again (it is **refractory**). Hodgkin and Huxley described this
with four equations in 1952. FitzHugh {cite:p}`fitzhugh1961` reduced them to
two, a fast voltage-like variable and a slow recovery variable, and Nagumo and
colleagues built the same equations as a circuit {cite:p}`nagumo1962`.
:::

## The model

:::{prf:definition} FitzHugh–Nagumo model
:label: def-fhn

$$
\ve\,\dot v = v - \frac{v^3}{3} - w,
\qquad
\dot w = v + a,
$$ (eq-fhn)

with $0 < \ve \ll 1$. Time is measured on the slow scale; $v$ is the fast
**activator**, $w$ the slow **inhibitor**.
:::

| term | in words | sign in the Jacobian |
|---|---|---|
| $v$ in $\ve\dot v$ | the activator amplifies itself | $+$ on the diagonal |
| $-v^3/3$ | the amplification saturates | keeps $v$ bounded |
| $-w$ | the inhibitor suppresses the activator | $\partial\dot v/\partial w = -1/\ve < 0$ |
| $v + a$ | the activator drives the inhibitor up; below $v = -a$ the inhibitor decays | $\partial\dot w/\partial v = 1 > 0$ |
| $\ve$ small | the activator is fast, the inhibitor slow | the delay in the feedback |

The cross effects have opposite signs at every point, and the self-term is
positive for $|v| < 1$.

:::{warning} A caricature, not a mechanism
FitzHugh–Nagumo is not derived from a list of molecules, as the drug model
was. It is the simplest pair of equations with the right *shape*: a cubic
fast nullcline and a slow negative feedback. What carries over to real
systems is the geometry, not the numbers.
:::

:::{prf:remark} FitzHugh's version
:label: rmk-fhn-original
:class: dropdown

FitzHugh wrote $\dot v = v - v^3/3 - w + I$, $\dot w = \phi\,(v + a - b\,w)$,
with an applied current $I$ and a small rate $\phi$. Measuring time in units
of $1/\phi$ gives {eq}`eq-fhn` with $\ve = \phi$ and $b = 0$. We drop the
decay term $-bw$ to keep the algebra short; then $I$ only shifts $w$, and the
parameter that matters is $a$.
:::

The small $\ve$ is used twice. First we freeze the slow inhibitor and study
the activator alone: it is a switch. Then we let the inhibitor move and watch
it flip the switch back and forth.

## The inhibitor frozen: a switch

Since $\ve$ is small, the activator moves much faster than the inhibitor.
Freeze the inhibitor at a value $w$. The activator then obeys a phase line
with one parameter,

$$
\ve\,\dot v = c(v) - w, \qquad c(v) = v - \frac{v^3}{3} .
$$ (eq-fhn-fast)

**The activator grows exactly when $c(v) > w$.** This is the graphical method
of [](../one-dimensional/saddle-node.md), with $c$ in place of the tolerated
dose $h$: draw the graph of $c$ and a horizontal line at height $w$. The
activator grows where the graph lies above the line, and the steady states
are the crossings. At a crossing $v^\ast$ the slope of the right-hand side is

$$
\lambda = \frac{c'(v^\ast)}{\ve} = \frac{1 - (v^\ast)^2}{\ve} ,
$$ (eq-fhn-fast-rate)

so by the slope test ({prf:ref}`prop-ir-slope-test`) a crossing where $c$
falls is stable and one where $c$ rises is unstable. The factor $1/\ve$ makes
both fast: the activator settles long before the inhibitor has moved.

The cubic is odd. It falls on $(-\infty, -1)$, rises on $(-1, 1)$ and falls
again on $(1, \infty)$, with a minimum $c(-1) = -2/3$ and a maximum
$c(1) = 2/3$, and it takes these values once more at $v = \pm 2$. A
horizontal line therefore meets it three times if $|w| < 2/3$ and once if
$|w| > 2/3$.

:::{prf:proposition} The fast equation is bistable
:label: prop-fhn-fast-bistable

Fix $w$ in the fast equation {eq}`eq-fhn-fast`.

1. If $|w| < 2/3$, there are three steady states, one in each of
   $(-2, -1)$, $(-1, 1)$ and $(1, 2)$. The outer two are stable and the middle
   one is unstable. Initial values below the middle state tend to the lower
   state; initial values above it tend to the upper state.
2. If $|w| > 2/3$, there is one steady state, with $v^\ast < -2$ if $w > 2/3$
   and $v^\ast > 2$ if $w < -2/3$. It is stable, and every initial value tends
   to it.

There are saddle-node bifurcations at $w = 2/3$, where the upper and middle
states meet at $v = 1$, and at $w = -2/3$, where the lower and middle states
meet at $v = -1$.
:::

:::{prf:proof}
:class: dropdown

On an interval where $c$ is monotone, it takes every value between its end
values exactly once. On $(-\infty, -1)$, $c$ falls from $+\infty$ to $-2/3$; on
$(-1, 1)$ it rises from $-2/3$ to $2/3$; on $(1, \infty)$ it falls from $2/3$
to $-\infty$. For $|w| < 2/3$ each interval contains one solution of
$c(v) = w$, and since $c(-2) = 2/3 > w$ and $c(2) = -2/3 < w$, the outer two
lie in $(-2, -1)$ and $(1, 2)$. For $w > 2/3$ only the first interval reaches
height $w$, and since $c(-2) = 2/3 < w$, its solution has $v^\ast < -2$; the
case $w < -2/3$ is its mirror image, since $c$ is odd. At $w = 2/3$ the middle
and the right interval both end at the maximum $c(1) = 2/3$, so their two
solutions meet at $v = 1$; at $w = -2/3$ the left and the middle interval meet
at $v = -1$.

Stability is {eq}`eq-fhn-fast-rate`: $c' < 0$ for $|v| > 1$ and $c' > 0$ for
$|v| < 1$. For the basins, $c(v) - w$ is positive below the lowest state
(since $c(v) \to +\infty$ as $v \to -\infty$) and changes sign at each simple
crossing, so the arrows point toward the outer states and away from the
middle one. A solution that starts between two neighbouring steady states is
monotone and cannot cross either, so it converges to the one its arrow
points to.
:::

```{code-cell} python
:tags: [hide-input]
w0 = 0.2
roots0 = sorted(x.real for x in np.roots([-1 / 3, 0, 1, -w0]) if abs(x.imag) < 1e-9)


def phase_arrows(ax, w, roots, vertical, color=C_INNER, lo=-2.4, hi=2.4, g=0.12):
    """Arrows of eps v' = c(v) - w between the steady states, drawn along the
    line at height w (vertical=False) or up the line at w (vertical=True)."""
    pts = [lo] + list(roots) + [hi]
    for x0, x1 in zip(pts[:-1], pts[1:]):
        if x1 - x0 < 3 * g:
            continue
        u, v = (x0 + g, x1 - g) if c(0.5 * (x0 + x1)) > w else (x1 - g, x0 + g)
        xy, xyt = ((w, v), (w, u)) if vertical else ((v, w), (u, w))
        ax.annotate("", xy=xy, xytext=xyt,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1))


fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3))

ax = axs[0]
ax.plot(VV[outer_l], c(VV[outer_l]), color=C_OUTER, lw=1.8)
ax.plot(VV[outer_r], c(VV[outer_r]), color=C_OUTER, lw=1.8)
ax.plot(VV[mid], c(VV[mid]), color=C_OUTER, lw=1.26, ls="--")
ax.fill_between(VV, w0, c(VV), where=c(VV) > w0, color=C_COMP, alpha=0.15, lw=0)
for y in (2 / 3, -2 / 3):
    ax.axhline(y, color=C_GREY, lw=0.7, ls=":")
ax.axhline(w0, color=C_INNER, lw=1.0)
phase_arrows(ax, w0, roots0, vertical=False)
for x in roots0:
    mark(ax, x, w0, abs(x) > 1, ms=5)
ax.text(-0.6, w0 + 0.08, f"$w = {w0}$", fontsize=8, color=C_INNER, ha="center")
ax.text(-2.35, 0.73, "$2/3$", fontsize=8, color=C_GREY)
ax.text(-2.35, -0.6, "$-2/3$", fontsize=8, color=C_GREY)
ax.set_xlim(-2.4, 2.4); ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("activator $v$"); ax.set_ylabel("$c(v)$")
ax.set_title("graphical method: $w$ is a horizontal line", fontsize=8.5)

ax = axs[1]
cubic(ax)
ax.fill_betweenx(VV, w0, c(VV), where=c(VV) > w0, color=C_COMP, alpha=0.15, lw=0)
ax.axvline(w0, color=C_INNER, lw=1.0)
phase_arrows(ax, w0, roots0, vertical=True)
for x in roots0:
    mark(ax, w0, x, abs(x) > 1, ms=5)
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-2.4, 2.4)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
ax.set_title("axes swapped: $w$ is a vertical phase line", fontsize=8.5)
fig.tight_layout();
```

The same picture twice, at $w = 0.2$. Left: the graph of $c$ and the line at
height $w$. Where the graph lies above the line (shaded) the activator grows;
the arrows of the phase line are drawn along the line, with filled circles for
stable steady states and an open one for the unstable state. Dotted: the
special heights $\pm 2/3$. Right: the axes swapped. The graph of $c$ lies on
its side, the line at $w$ stands up, and its arrows are the phase line at that
$w$. Repeating the right panel for every $w$ gives the bifurcation diagram
below: the curve is the same for every $w$, and only the vertical line moves.

### The bifurcation diagram is an S

The steady states satisfy $w = c(v)$, which gives the parameter as a function
of the state. There is no need to solve the cubic for $v$: letting $v$ run and
plotting the points $(c(v), v)$, with $w$ across and $v$ up, traces out the
diagram. It is the graph of $c$ turned on its side, an S. The outer branches
lean back, where $c' < 0$, and are stable; the middle branch leans forward and
is unstable. The turning points of $c$ become the rightmost and leftmost points
of the S, $(w, v) = (\pm 2/3, \pm 1)$, where $c(v) = w$ and $c'(v) = 0$
({prf:ref}`rem-bf-slope-fails`). There a stable and the unstable branch meet
and disappear: two saddle-node bifurcations, or **folds**.

```{code-cell} python
:tags: [hide-input]
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
for ax in axs:
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-2.4, 2.4)
    ax.set_xlabel("inhibitor $w$ (frozen: a parameter)")

ax = axs[0]
ax.fill_betweenx(VV, -1.5, c(VV), color=C_COMP, alpha=0.12, lw=0, zorder=0)
cubic(ax)
for w in (-0.9, w0, 0.9):
    roots = sorted(x.real for x in np.roots([-1 / 3, 0, 1, -w]) if abs(x.imag) < 1e-9)
    ax.plot([w, w], [-2.4, 2.4], color="0.85", lw=0.8, zorder=1)
    phase_arrows(ax, w, roots, vertical=True)
    for x in roots:
        mark(ax, w, x, abs(x) > 1, ms=5)
ax.plot([2 / 3, -2 / 3], [1, -1], "o", color=C_INNER, ms=5, zorder=7)
ax.annotate("fold at $(2/3, 1)$", xy=(2 / 3, 1), xytext=(0.75, 1.6), fontsize=8,
            color=C_INNER, arrowprops=dict(arrowstyle="-", color=C_INNER, lw=0.6))
ax.annotate("fold at $(-2/3, -1)$", xy=(-2 / 3, -1), xytext=(-1.45, -1.75), fontsize=8,
            color=C_INNER, arrowprops=dict(arrowstyle="-", color=C_INNER, lw=0.6))
ax.set_ylabel("activator $v$")
ax.set_title("steady states and phase lines", fontsize=8.5)

ax = axs[1]
cubic(ax, color=C_GREY)
up, lo_ = VV[VV >= 1], VV[VV <= -1]
ax.plot(c(up), up, color=C_INNER, lw=2.2)
ax.plot(c(lo_), lo_, color=C_INNER, lw=2.2)
for xy0, xy1 in (((c(1.85), 1.85), (c(1.6), 1.6)),     # w up along the upper branch
                 ((2 / 3, 0.85), (2 / 3, -1.85)),        # jump down
                 ((c(-1.6), -1.6), (c(-1.35), -1.35)),   # w down along the lower branch
                 ((-2 / 3, -0.85), (-2 / 3, 1.85))):     # jump up
    ax.annotate("", xy=xy1, xytext=xy0,
                arrowprops=dict(arrowstyle="-|>", color=C_INNER, lw=1.5))
ax.text(-0.45, 1.35, "$w$ up", fontsize=8, color=C_INNER)
ax.text(0.75, -0.4, "jump\ndown", fontsize=8, color=C_INNER)
ax.text(0.05, -1.3, "$w$ down", fontsize=8, color=C_INNER)
ax.text(-1.05, 0.2, "jump\nup", fontsize=8, color=C_INNER)
ax.annotate("", xy=(-2 / 3, 2.32), xytext=(2 / 3, 2.32),
            arrowprops=dict(arrowstyle="<->", color="#777777", lw=0.9))
ax.text(0.72, 2.25, "bistable", fontsize=8, color="#777777")
ax.set_title("hysteresis", fontsize=8.5)
fig.tight_layout();
```

Bifurcation diagram of the fast equation {eq}`eq-fhn-fast` in the inhibitor
$w$. Left: solid branches are stable, the dashed middle branch unstable, and
the activator grows in the shaded region to the left of the S; the vertical
arrows are the phase lines at $w = -0.9$, $0.2$ and $0.9$. The middle branch
meets the upper one at the fold $(2/3, 1)$ and the lower one at the fold
$(-2/3, -1)$. Right: raising $w$ slowly, the activator follows the upper
branch until the fold at $w = 2/3$ and drops to $v = -2$; lowering $w$ again,
it follows the lower branch until the fold at $w = -2/3$ and jumps up to
$v = 2$. Between the folds (double arrow) both outer branches are stable.

Between the folds, which branch holds the activator depends on its history:
the hysteresis loop of [](../one-dimensional/saddle-node.md), in the right
panel. There an experimenter had to turn the dose to go round the loop. In
FitzHugh–Nagumo nobody has to. By $\dot w = v + a$, the inhibitor rises while
$v > -a$ and falls while $v < -a$. If $|a| < 1$, it rises on the whole upper
branch, where $v > 1$, and falls on the whole lower branch, where $v < -1$:
exactly the motion that drives the activator round the loop. Whether the loop
really closes into an oscillation is what the rest of the page settles.

The S is also the $v$-nullcline of the full model {eq}`eq-fhn`, the critical
manifold of [](quasi-steady-state.md). **The fast nullcline is the bifurcation
diagram of the fast equation, with the slow variable as the parameter.** This
holds for every fast/slow system ({prf:ref}`rem-bf-nullcline`); in the drug
model the diagram had no folds.

:::{prf:remark} Where the critical manifold attracts
:label: rmk-fhn-attracting
:class: dropdown

The pull toward the critical manifold $w = c(v)$ is
$\partial(\ve\dot v)/\partial v = 1 - v^2$. It is negative on the outer
branches $|v| > 1$, which attract in the sense of {eq}`eq-qss-attracting`,
positive on the middle branch, which repels, and **zero at the folds**
$v = \pm 1$. For the drug it was $-\delta(1 + N/N_c)$, never zero. So
{prf:ref}`prop-pp-catch-up` and Tikhonov's theorem {prf:ref}`thm-qss-tikhonov`
apply on each outer branch, away from the folds, and break down at them.
:::

## The phase plane

Now let the inhibitor move. We draw the phase plane with $w$ across and $v$
up, as in the bifurcation diagram, so that the fast motion is vertical. The
nullclines ({prf:ref}`def-pp-nullclines`) are the S, $w = c(v)$, on which the
arrows are horizontal, and the line $v = -a$, on which they are vertical. They
cross once, at the steady state $(w^\ast, v^\ast) = (c(-a), -a)$. It lies on
the middle branch of the S if $|a| < 1$ and on an outer branch if $|a| > 1$.

**Directions.** Left of the S the activator rises, and right of it it falls.
Above the line $v = -a$ the inhibitor grows and the state moves right; below
it the state moves left. So the state moves right above the steady state,
down to its right, left below it and up to its left: **the flow turns
clockwise**. This is the rotation that opposite-sign cross effects allow
([](#fhn-oscillate)); in the drug model the arrows sent the state toward a
steady state or away from it, but never around it.

**Speeds.** The vertical velocity carries the factor $1/\ve$, the horizontal
one does not. Away from the S the arrows are nearly vertical, and only within
a distance of order $\ve$ of the S do the two velocities become comparable.

```{code-cell} python
:tags: [hide-input]
eps, a = 0.05, 0.5
gw, gv = np.meshgrid(np.linspace(-1.3, 1.3, 14), np.linspace(-2.4, 2.4, 15))
U, V = gv + a, (c(gv) - gw) / eps                   # (w', v')
L = np.hypot(U / 2.8, V / 5.0)                       # length in axis fractions
fig, ax = plt.subplots(figsize=(5.4, 3.6))
ax.quiver(gw, gv, 0.045 * U / L, 0.045 * V / L * 5 / 2.8, color="0.65",
          angles="xy", scale_units="xy", scale=1, width=0.004)
cubic(ax)
ax.axhline(-a, color=C_INNER, lw=1.8)
mark(ax, c(-a), -a, False)
ax.plot([2 / 3, -2 / 3], [1, -1], "o", color=C_OUTER, ms=4, zorder=6)
for w, v, s in ((-1.05, 1.6, "$\\dot v > 0,\\ \\dot w > 0$"),
                (0.55, 1.95, "$\\dot v < 0,\\ \\dot w > 0$"),
                (0.55, -1.15, "$\\dot v < 0,\\ \\dot w < 0$"),
                (-1.25, -2.2, "$\\dot v > 0,\\ \\dot w < 0$")):
    ax.text(w, v, s, fontsize=7.5, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
ax.text(1.0, -a + 0.1, "$\\dot w = 0$", fontsize=8, color=C_INNER)
ax.text(0.74, 0.62, "$\\dot v = 0$", fontsize=8, color=C_OUTER,
        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-2.5, 2.5)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
fig.tight_layout();
```

Nullclines and directions of {eq}`eq-fhn` for $a = 0.5$, $\ve = 0.05$. Blue:
the $v$-nullcline, the S (dashed between the folds, marked by dots). Orange:
the $w$-nullcline $v = -a$. Grey: the direction of the arrow at a grid of
points, all drawn the same length. Away from the S they are nearly vertical;
near it they turn horizontal. The signs label the four regions; around the
steady state (open circle) the arrows turn clockwise.

So a trajectory drops or rises almost vertically onto the S, then moves
slowly along it. Where the steady state sits on the S decides what happens
next.

## The steady state

:::{prf:proposition} Stability of the steady state
:label: prop-fhn-stability

At $(w^\ast, v^\ast) = (c(-a), -a)$ the Jacobian of {eq}`eq-fhn`, ordered as
$(v, w)$, is

$$
J = \begin{pmatrix} (1 - a^2)/\ve & -1/\ve \\ 1 & 0 \end{pmatrix},
\qquad
T = \frac{1 - a^2}{\ve}, \quad \Delta = \frac1\ve > 0 .
$$ (eq-fhn-jacobian)

The steady state is never a saddle, and

- for $|a| > 1$ it is stable; it lies on an outer, stable branch of the S;
- for $|a| < 1$ it is unstable; it lies on the middle branch;
- it is a spiral when $(1 - a^2)^2 < 4\ve$, that is, close to the folds
  (for $\ve = 0.05$: $0.74 < |a| < 1.20$);
- at $a = \pm 1$ it sits exactly on a fold, $T = 0$, and the eigenvalues are
  $\pm i/\sqrt\ve$.
:::

:::{prf:proof}
:class: dropdown
The partial derivatives of $\dot v = (v - v^3/3 - w)/\ve$ are
$(1 - v^2)/\ve$ and $-1/\ve$; those of $\dot w = v + a$ are $1$ and $0$.
At $v^\ast = -a$ this gives {eq}`eq-fhn-jacobian`. Since $\Delta > 0$, the
classification of {prf:ref}`prop-pp-trace-det` gives stability for $T < 0$,
instability for $T > 0$, and complex eigenvalues when $T^2 < 4\Delta$, which
is $(1 - a^2)^2/\ve^2 < 4/\ve$. At $T = 0$ the eigenvalues solve
$\lambda^2 + 1/\ve = 0$.
:::

The steady state is unstable exactly when it sits on the middle branch, where
the frozen switch is unstable too.

## An oscillation

For $|a| < 1$ trajectories leave the steady state, and there is no other one
to go to. Two steps show that they settle on a closed orbit: a region they
cannot leave, and the Poincaré–Bendixson theorem.

### A trapping region

The arrows turn clockwise, but that alone does not keep a trajectory from
spiralling out to infinity. We need a bounded region that the flow never
leaves. A box does not work: the inhibitor has no decay, so on the right side
of any box, above the line $v = -a$, the inhibitor grows and the flow leaves.
But there the state is also right of the S, where $v$ falls at a speed of
order $1/\ve$, so a side that leans slightly to the left is crossed inward.
The same holds, mirrored, on the left.

:::{prf:proposition} A trapping hexagon
:label: prop-fhn-trap

Let $|a| \le 2$ and $\ve > 0$. There are hexagons $H$ of the shape in the
figure below, as large as we like, such that no trajectory of {eq}`eq-fhn`
that starts in $H$ leaves it. In particular every trajectory stays bounded.
:::

:::{prf:proof}
:class: dropdown

Choose $V \ge 3$, $W$ and $s > 0$ with

$$
W - s > 1, \qquad W < \frac{V^3}{3} - V, \qquad
s\,\Big(W - s - \frac23\Big) \ge \ve\,(V + |a|)^2 ,
$$ (eq-fhn-trap)

and let $H$ be the hexagon with corners, in the $(w, v)$ plane,

$$
(-W, V),\quad (W - s, V),\quad (W, -a),\quad (W, -V),\quad
(-W + s, -V),\quad (-W, -a) .
$$

We check that on each side the arrows point into $H$ or along the side.

- **Top**, $v = V$, $-W \le w \le W - s$. Here $c(V) = -(V^3/3 - V) < -W \le w$,
  so $\dot v = (c(V) - w)/\ve < 0$: the arrows point down.
- **Bottom**, $v = -V$, $-W + s \le w \le W$. Here $w \le W < c(-V)$, so
  $\dot v > 0$: up.
- **Right, lower part**, $w = W$, $-V \le v \le -a$: $\dot w = v + a \le 0$,
  so the arrows point left or along the side.
- **Left, upper part**, $w = -W$, $-a \le v \le V$: $\dot w \ge 0$, right.
- **Right, slanted part**, from $(W - s, V)$ to $(W, -a)$. On it
  $w \ge W - s > 1$, while $c(v) \le 2/3$ for $v \ge -2$, so
  $\dot v \le -(W - s - 2/3)/\ve$; and $0 \le \dot w = v + a \le V + |a|$.
  A normal pointing out of $H$ is $(V + a, s)$, in components $(w, v)$. The
  velocity $(\dot w, \dot v)$ has with it the product
  $(V + a)\,\dot w + s\,\dot v \le (V + |a|)^2 - s\,(W - s - 2/3)/\ve \le 0$
  by {eq}`eq-fhn-trap`: the arrows point in.
- **Left, slanted part**, from $(-W + s, -V)$ to $(-W, -a)$: the same with all
  signs reversed, using $c(v) \ge -2/3$ for $v \le 2$.

For large hexagons, take $s = (W - 2/3)/2$. Then the last condition in
{eq}`eq-fhn-trap` reads $(W - 2/3)^2 \ge 4\ve\,(V + |a|)^2$, which holds once
$W$ is large compared with $\sqrt\ve\,V$, and the second condition allows
$W$ to grow like $V^3$.
:::

For $\ve = 0.05$ and $|a| \le 1.5$, the numbers $V = 3$, $W = 5$, $s = 1/4$
satisfy {eq}`eq-fhn-trap`. The hexagon is far larger than any orbit; its job
is only to show that nothing escapes.

```{code-cell} python
:tags: [hide-input]
eps, a = 0.05, 0.5
Vh, Wh, sh = 3.0, 5.0, 0.25
hexagon = np.array([(-Wh, Vh), (Wh - sh, Vh), (Wh, -a), (Wh, -Vh),
                    (-Wh + sh, -Vh), (-Wh, -a), (-Wh, Vh)])
fig, ax = plt.subplots(figsize=(6.2, 3.4))
phase_plane(ax, lambda w, v: v + a, lambda w, v: (c(v) - w) / eps,
            (-5.6, 5.6), (-3.5, 3.5), trajs=[(-4.9, 2.5), (4.95, -2.5), (4.9, 1.0),
                                              (-4.95, -1.5), (0.0, 2.95), (0.0, -2.95)],
            t_max=12, density=1.0, null=False)
cubic(ax, lw=1.4)
ax.plot(*hexagon.T, color=C_INNER, lw=1.6, zorder=5)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
fig.tight_layout();
```

The trapping hexagon of {prf:ref}`prop-fhn-trap` (orange) for $a = 0.5$,
$\ve = 0.05$, with $V = 3$, $W = 5$, $s = 1/4$. The two slanted sides, above
the line $v = -a$ on the right and below it on the left, lean only slightly.
Trajectories from points on the boundary (squares) all move inward, drop onto
the S and settle on a loop around it, the closed orbit of
{prf:ref}`prop-fhn-closed-orbit`.

### Poincaré–Bendixson

In the plane, a trajectory trapped in a bounded region that never comes to
rest cannot cross its own earlier laps, so the laps line up like the grooves
of a record and close in on a loop.

:::{prf:theorem} Poincaré–Bendixson
:label: thm-fhn-pb

Let $R$ be a closed, bounded region of the plane that contains no steady
state. If a trajectory starts in $R$ and stays in $R$ for all later times,
then it is a closed orbit or it approaches one.
:::

For a proof see Strogatz §7.3 {cite:p}`strogatz2015`. The hexagon traps
trajectories; to keep the steady state out of the region, we cut a small hole
around it.

:::{prf:proposition} FitzHugh–Nagumo oscillates
:label: prop-fhn-closed-orbit

If $|a| < 1$, then for every $\ve > 0$ the model {eq}`eq-fhn` has a closed
orbit, and every trajectory except the steady state approaches a closed orbit.
:::

:::{prf:proof}
:class: dropdown

For $|a| < 1$ both eigenvalues of {eq}`eq-fhn-jacobian` have positive real
part ({prf:ref}`prop-fhn-stability`), whatever $\ve$ is. So there is a small
ellipse $E$ around the steady state that every trajectory crosses outward: for
the linearized flow, take a level curve of a quadratic form that grows along
it, and near the steady state the nonlinear terms are too small to change the
direction of crossing. Let $R$ be a trapping hexagon $H$ of
{prf:ref}`prop-fhn-trap` with the inside of $E$ removed, $H$ chosen to contain
both the ellipse and the starting point. Trajectories that start in $R$
cannot leave through the hexagon, nor enter the ellipse, so they stay in $R$.
The only steady state is inside $E$, so $R$ contains none, and
{prf:ref}`thm-fhn-pb` applies.
:::

In the trapping-hexagon figure above, the trajectories end on a closed orbit
that runs round the hysteresis loop of the fast equation: right along the
upper branch of the S to its fold, a drop to the lower branch, left along it
to the lower fold, and a jump back up. The same trip, made once, is the pulse
of an excitable system.

## Excitability: just past the fold

Take $a$ slightly larger than $1$. The steady state $v^\ast = -a$ lies on the
lower branch just past the fold and is stable ({prf:ref}`prop-fhn-stability`).
Nothing oscillates. Now kick the activator. A small kick that leaves $v$ below
the middle branch dies away: the state falls back to the lower branch and to
rest. A large kick that puts $v$ above the middle branch is carried up to the
upper branch, and the state makes the trip of the closed orbit once: crawl
right to the upper fold, jump down to $v = -2$, crawl left back to rest. One
full pulse.

```{code-cell} python
:tags: [hide-input]
a = 1.3                                   # exaggerated, so the gap shows
ws = c(-a)                                # rest: (w*, v*) = (c(-a), -a)
root = np.sqrt(12 - 3 * a**2)
vm, vu = (a - root) / 2, (a + root) / 2   # middle and upper states at w*


def arrow(ax, p, q, fast, color=C_EXACT, lw=1.4):
    """Arrow from p to q: two heads for a fast move, one for a slow one."""
    head = dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=0, shrinkB=0,
                mutation_scale=11)
    ax.annotate("", xy=q, xytext=p, arrowprops=head)
    if fast:
        m = (0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]))
        ax.annotate("", xy=m, xytext=p, arrowprops=head)


def along(ax, v0, v1, color=C_EXACT, lw=2.2):
    """Slow crawl along the S from v0 to v1, with one head in the middle."""
    v = np.linspace(v0, v1, 200)
    ax.plot(c(v), v, color=color, lw=lw, solid_capstyle="round", zorder=4)
    k = len(v) // 2
    arrow(ax, (c(v[k - 1]), v[k - 1]), (c(v[k]), v[k]), False, color, lw=1.4)


fig, ax = plt.subplots(figsize=(6.0, 3.8))
cubic(ax, color="0.65", lw=1.4)
ax.axhline(-a, color=C_INNER, lw=1.0, alpha=0.7)

# small kick: stays below the threshold and falls back
ax.annotate("", xy=(ws, vm - 0.22), xytext=(ws, -a),
            arrowprops=dict(arrowstyle="-|>", color=C_OUTER, lw=1.1, ls=":",
                            shrinkA=0, shrinkB=0))
ax.annotate("", xy=(ws - 0.02, -a + 0.04), xytext=(ws - 0.02, vm - 0.22),
            arrowprops=dict(arrowstyle="-|>", color=C_OUTER, lw=1.2, shrinkA=0,
                            shrinkB=0, connectionstyle="arc3,rad=0.6"))

# large kick: crosses the threshold, then the pulse
kick = vm + 0.3
ax.plot([ws + 0.035, ws + 0.035], [-a, kick], color=C_EXACT, lw=1.2, ls=":")
arrow(ax, (ws + 0.035, kick), (ws + 0.035, vu - 0.04), True)     # jump up
along(ax, vu, 1.0)                                                # crawl right
arrow(ax, (2 / 3, 1.0), (2 / 3, -1.96), True)                     # drop
along(ax, -2.0, -a - 0.02)                                        # crawl back

mark(ax, ws, -a, True, ms=7)
ax.plot([ws], [vm], "o", ms=6, mfc="white", mec=C_EXACT, zorder=6)
ax.plot([2 / 3], [1], "o", color="0.65", ms=4, zorder=5)

txt = dict(fontsize=8, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
ax.annotate("threshold:\nthe middle branch", xy=(ws, vm), xytext=(-1.25, 0.2),
            arrowprops=dict(arrowstyle="-", color="#555555", lw=0.6), **txt)
ax.text(ws - 0.3, -a - 0.3, "rest", **txt)
ax.text(-1.25, -0.95, "small kick:\nfalls back", color=C_OUTER, **txt)
ax.text(-0.42, 0.45, "large kick:\njump up", **txt)
ax.text(0.05, 2.15, "crawl right (slow)", **txt)
ax.text(0.72, -0.45, "drop\n(fast)", **txt)
ax.text(-0.15, -2.35, "crawl back: refractory", **txt)
ax.text(1.0, -a + 0.08, "$\\dot w = 0$", fontsize=8, color=C_INNER)
ax.set_xlim(-1.3, 1.3); ax.set_ylim(-2.45, 2.45)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
fig.tight_layout();
```

A sketch of excitability for small $\ve$, with the steady state (filled) just
past the lower fold, drawn at $a = 1.3$ so that the gap shows. Double arrows
are fast, single arrows slow. A small kick (blue) stays below the middle
branch, the threshold, and falls straight back to rest. A large kick (dotted)
crosses it: the activator jumps up to the upper branch, crawls right to the
upper fold, drops to the lower branch and crawls back to rest. On the way
back the threshold is further away than at rest: the system is refractory.

:::{prf:definition} Excitable
:label: def-fhn-excitable

A system with a stable steady state is **excitable** if small perturbations
of the steady state decay directly, while every perturbation beyond a
**threshold** produces a large excursion of the same stereotyped shape before
the state returns to rest (**all or none**).
:::

Why does the threshold sit at the middle branch, and why does a kick beyond it
always make the whole trip? A kick changes $v$ at once and leaves $w$ at its
rest value $w^\ast = c(-a)$. Right after the kick the fast equation
{eq}`eq-fhn-fast` at this frozen $w^\ast$ takes over, and
{prf:ref}`prop-fhn-fast-bistable` decides: below its middle state $v$ returns
to the lower state, which is the rest state; above it, $v$ rises to the upper
branch. The next proposition says how large the threshold is and how fast
the decision is made.

:::{prf:proposition} The threshold just past the fold
:label: prop-fhn-threshold

Let $1 < a < 2$. With $w$ frozen at $w^\ast = c(-a)$, the fast equation has
the rest state $-a$, the middle state $v_m = \big(a - \sqrt{12 - 3a^2}\big)/2$
and an upper state. So as $\ve \to 0$ a kick $\Delta v$ from rest fires
exactly when

$$
\Delta v > \Delta v_\mathrm{th}(a) = v_m + a = \frac{3a - \sqrt{12 - 3a^2}}{2} .
$$ (eq-fhn-threshold)

The threshold rises from $0$ at $a = 1$ to $3$ at $a = 2$. Just past the
fold, $a = 1 + \delta$ with $\delta$ small,

$$
\Delta v_\mathrm{th} \approx 2\delta ,
$$

and the rest state attracts and the middle state repels at the rate
$2\delta/\ve$, to leading order.
:::

:::{prf:proof}
:class: dropdown

For any $u$, $c(v) - c(u) = (v - u)\big(1 - (v^2 + uv + u^2)/3\big)$. With
$u = -a$, the other solutions of $c(v) = c(-a)$ solve $v^2 - av + a^2 - 3 = 0$,
that is $v = \big(a \pm \sqrt{12 - 3a^2}\big)/2$, real for $a \le 2$. The minus
sign gives a value between $-1$ (at $a = 1$) and $1$ (at $a = 2$), the middle
state; the plus sign gives one above $1$. At $a = 1 + \delta$,
$12 - 3a^2 = 9 - 6\delta - 3\delta^2$, so
$\sqrt{12 - 3a^2} = 3 - \delta + O(\delta^2)$ and
$v_m = -1 + \delta + O(\delta^2)$. By {eq}`eq-fhn-fast-rate` the rates are
$(1 - a^2)/\ve = -(2\delta + \delta^2)/\ve$ at the rest state and
$(1 - v_m^2)/\ve = (2\delta + O(\delta^2))/\ve$ at the middle state.
:::

At $a = 1.1$ the threshold is $0.20$. A kick beyond it can only end in a full
pulse: on the upper branch $\dot w = v + a > 0$ and there is no steady state,
so the state must crawl to the upper fold and come round. Nothing on this trip
depends on the size of the kick. The response is set by the S, not by the
kick, and that is the all-or-none response.

The proposition also says why excitability lives **just past the fold**.
First, the threshold grows with the distance from the fold, to $3$ at $a = 2$;
beyond $a = 2$, $w^\ast > 2/3$, the fast equation at $w^\ast$ has only the
lower state, and no kick of the activator alone produces a pulse. Second, the
decision between returning and firing is made at the rate $2(a - 1)/\ve$. The
threshold is sharp only when this is fast compared with the inhibitor, that
is when $a - 1$ is large compared with $\ve$. Closer to the fold a kicked state
hesitates near the middle branch long enough for $w$ to move, and the
threshold blurs.

```{code-cell} python
:tags: [hide-input]
eps, a = 0.05, 1.1
v0, w0 = -a, c(-a)


def kicked(kicks, T=10.0, n=5000):
    """Kick v by dv at the given times; return t, v, w."""
    ts, vs, ws = [], [], []
    y, t0 = [v0, w0], 0.0
    for tk, dv in list(kicks) + [(T, 0.0)]:
        if tk > t0:
            t, (v, w) = solve(y, tk - t0, eps, a, n=max(2, int(n * (tk - t0) / T)))
            ts.append(t + t0); vs.append(v); ws.append(w)
            y = [v[-1], w[-1]]
        y = [y[0] + dv, y[1]]
        t0 = tk
    return np.concatenate(ts), np.concatenate(vs), np.concatenate(ws)


def fires(dv):
    t, v, w = kicked([(0.0, dv)], T=4.0, n=2000)
    return v.max() > 0.5


lo, hi = 0.0, 1.0
for _ in range(30):
    m = 0.5 * (lo + hi)
    lo, hi = (lo, m) if fires(m) else (m, hi)
thr = hi

fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1), gridspec_kw={"width_ratios": [1, 1.3]})
ax = axs[0]
cubic(ax)
ax.axhline(-a, color=C_INNER, lw=1.0, alpha=0.7)
for dv, col in ((thr - 0.03, C_OUTER), (thr + 0.03, C_EXACT)):
    t, v, w = kicked([(0.0, dv)], T=6.0)
    ax.plot(w, v, color=col, lw=1.1)
    ax.plot([w0, w0], [v0, v0 + dv], color=col, lw=1.1, ls=":")
mark(ax, w0, v0, True)
ax.set_xlim(-1.0, 1.2); ax.set_ylim(-2.4, 2.4)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
ax.set_title(rf"$a = {a}$: excitable")

ax = axs[1]
for dv, col, lab in ((thr - 0.03, C_OUTER, "kick below threshold"),
                     (thr + 0.03, C_EXACT, "kick above threshold")):
    t, v, w = kicked([(0.5, dv)], T=10.0)
    ax.plot(t, v, color=col, lw=1.2, label=lab)
t, v, w = kicked([(0.5, 0.8), (2.6, 0.8), (5.6, 0.8)], T=10.0)
ax.plot(t, v, color=C_INNER, lw=1.0, alpha=0.8, label="three kicks of 0.8")
for tk in (0.5, 2.6, 5.6):
    ax.annotate("", xy=(tk, -1.6), xytext=(tk, -2.3),
                arrowprops=dict(arrowstyle="-|>", color=C_INNER, lw=0.8))
ax.set_xlabel("time"); ax.set_ylabel("activator $v$"); ax.set_xlim(0, 10)
ax.set_ylim(-2.4, 2.4)
ax.legend(frameon=False, fontsize=7.5, loc="upper right")
fig.tight_layout();
print(f"threshold kick at eps = {eps}: {thr:.3f}")
```

Excitability at $a = 1.1$, $\ve = 0.05$. Left: kicks from rest (dotted)
just below and just above threshold; the first returns at once, the second
makes a full excursion around the S. Right: the same kicks in time, and, in
orange, three equal kicks of size $0.8$ (arrows). The second comes during the
pulse, while the inhibitor is still high, and does nothing: the system is
**refractory**.

For $\ve = 0.05$ the threshold kick is about $0.33$, larger than the
$\ve \to 0$ value $0.2$, because $w$ moves while $v$ moves: at $a = 1.1$ the
middle branch repels at a rate of only about $2(a - 1)/\ve = 4$. The
threshold is not perfectly sharp: kicks in a very narrow band give
intermediate responses.

:::{prf:remark} Refractoriness from the phase plane
:label: rmk-fhn-refractory
:class: dropdown

During the crawl back along the lower branch, $w$ is larger than at rest.
The threshold is the middle branch at the current $w$, which rises as $w$
increases, so a larger kick is needed. Right after the upper fold, at
$w = 2/3$, the lower and middle branches are at $v = -2$ and $v = 1$: a kick
of $3$ is needed. The refractory period is the time it takes $w$ to come back
down.
:::

## Summary

FitzHugh–Nagumo has three behaviours, set by where the slow nullcline
$v = -a$ crosses the S: on the middle branch it oscillates, just past a fold
it is excitable, and far past a fold it is simply stable. The first two are
side by side below, in the plane and in time.

```{code-cell} python
:tags: [hide-input]
eps, T = 0.05, 8.0
C_TRAJ = (C_EXACT, "#CC79A7")             # one colour per starting point
cases = ((0.5, [(-0.38, -0.3), (-1.2, 2.3)], r"$a = 0.5$: steady state on the middle branch"),
         (1.1, [(0.6, 2.2), (-0.62, -0.55)], r"$a = 1.1$: just past the lower fold"))
fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6), sharey=True,
                        gridspec_kw={"height_ratios": [1.3, 1]})
for k, (a, starts, title) in enumerate(cases):
    F = lambda w, v: v + a
    G = lambda w, v: (c(v) - w) / eps
    ax, axt = axs[0, k], axs[1, k]
    phase_plane(ax, F, G, (-1.4, 1.4), (-2.5, 2.5))
    t = np.linspace(0, T, 4000)
    for z0, col in zip(starts, C_TRAJ):
        s = solve_ivp(lambda t, z: [F(*z), G(*z)], (0, T), z0, rtol=1e-9,
                      atol=1e-11, method="LSODA", dense_output=True)
        w, v = s.sol(t)
        ax.plot(w, v, color=col, lw=1.1, zorder=4)
        ax.plot(*z0, "s", color=col, ms=4, zorder=5)
        axt.plot(t, v, color=col, lw=1.1)
    axt.axhline(-a, color=C_INNER, lw=0.8, alpha=0.7)
    ax.set_xlabel("inhibitor $w$"); ax.set_title(title, fontsize=8.5)
    axt.set_xlabel("time"); axt.set_xlim(0, T)
axs[0, 0].set_ylabel("activator $v$"); axs[1, 0].set_ylabel("activator $v$")
axs[0, 0].plot([], [], color=C_OUTER, lw=2, label=r"$\dot v = 0$")
axs[0, 0].plot([], [], color=C_INNER, lw=2, label=r"$\dot w = 0$")
axs[0, 0].legend(loc="lower left", frameon=False, fontsize=7.5)
fig.tight_layout();
```

Top: phase planes of {eq}`eq-fhn` with $\ve = 0.05$. Blue: the S. Orange: the
$w$-nullcline $v = -a$. Pale lines: streamlines, thicker where the flow is
faster. Squares mark the two starting points, black and pink. Bottom: the
activator $v$ against time for the same two trajectories, in the same
colours, with the level $v = -a$ of the steady state in orange. Left,
$a = 0.5$: the steady state is on the middle branch and unstable. The black
trajectory starts next to it, the pink one far away, and both end on the same
closed orbit ({prf:ref}`prop-fhn-closed-orbit`): pulses that repeat for ever.
Right, $a = 1.1$: the steady state is just past the lower fold and stable.
Both trajectories make one excursion around the S, a single pulse, and then
rest at $v = -a$ ({prf:ref}`prop-fhn-threshold`).

### Two kinds of threshold

| | drug and cells | FitzHugh–Nagumo, excitable |
|---|---|---|
| what decides | which side of the separatrix you start on | whether the kick crosses the middle branch of the S |
| steady states | two stable, one saddle | one, stable |
| the threshold is | the stable manifold of a saddle | a repelling branch of the fast nullcline |
| after the decision | stays at a different steady state forever | returns to the *same* rest state after a pulse |
| sharpness | exactly sharp | sharp up to an exponentially small band |
| memory | permanent: a switch | none after the pulse; brief refractoriness |

The switch remembers which side it was put on. The excitable system forgets,
and only reports, with a stereotyped pulse, that something large enough
happened.
