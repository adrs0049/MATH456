---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/conserved-quantities.pdf
    id: phase-plane-conserved-quantities-pdf
downloads:
  - id: phase-plane-conserved-quantities-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Conserved Quantities and Limit Cycles

:::{tip} Big Idea
After a perturbation, a predator-prey cycle continues at a new amplitude
indefinitely, while a firing nerve returns to the same rhythm within one
pulse. Which of the two occurs can be read off from the equations without
solving them. A **conserved quantity** confines every trajectory to one of
its level curves, so a perturbation moves the system to another level curve
permanently: a center. The **divergence** measures how the flow changes
area. A flow that preserves area has no attractors, a flow that contracts
area everywhere has no closed orbit (Bendixson, Dulac), and an attracting
**limit cycle** requires a region where the flow expands area and a region
where it contracts area. FitzHugh–Nagumo has such a limit cycle.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER, C_INNER, C_COMP = "#111111", "#0072B2", "#D55E00", "#009E73"
C_GREY, C_PINK = "#999999", "#CC79A7"

c = lambda v: v - v**3 / 3            # FitzHugh-Nagumo: fast nullcline w = c(v)
VV = np.linspace(-2.3, 2.3, 400)


def fhn(a, eps):
    """FitzHugh-Nagumo in the order (w, v), as in fitzhugh-nagumo.md."""
    return lambda t, z: [z[1] + a, (c(z[1]) - z[0]) / eps]


def lv(t, z):
    """Lotka-Volterra, prey x and predators y."""
    x, y = z
    return [x * (1 - y), y * (x - 1)]


def run(f, z0, T, n=4000, t0=0.0):
    s = solve_ivp(f, (t0, t0 + T), z0, rtol=1e-10, atol=1e-12, method="LSODA",
                  dense_output=True)
    t = np.linspace(t0, t0 + T, n)
    return t, s.sol(t)


def flow(ax, F, G, xlim, ylim, density=0.9):
    X, Y = np.meshgrid(np.linspace(*xlim, 160), np.linspace(*ylim, 160))
    U, V = F(X, Y), G(X, Y)
    speed = np.hypot(U, V)
    ax.streamplot(X, Y, U, V, density=density, color="thistle", arrowsize=0.9,
                  linewidth=0.4 + 1.4 * np.sqrt(speed / np.nanmax(speed)), zorder=1)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)


def mark(ax, x, y, stable, ms=5):
    style = dict(color=C_COMP) if stable else dict(mfc="white", mec=C_EXACT)
    ax.plot([x], [y], "o", ms=ms, zorder=6, **style)


def cubic(ax, lw=1.4, color=C_OUTER):
    """The S, with w across and v up: outer branches solid, middle dashed."""
    for part, ls in ((VV <= -1, "-"), (VV >= 1, "-"), (np.abs(VV) <= 1, "--")):
        ax.plot(c(VV[part]), VV[part], color=color, lw=lw, ls=ls)
```

## Two oscillators, one kick

Take two oscillators from biology: predators and their prey, and a nerve cell
firing. The first is Lotka–Volterra. In scaled units, prey $x$ grow and are
eaten, and predators $y$ die without prey and grow by eating them:

$$
\dot x = x(1 - y), \qquad \dot y = y(x - 1)
$$

({prf:ref}`ex-cq-lotka-volterra` below). The second is FitzHugh–Nagumo
{eq}`eq-fhn` with $a = 0.5$, where every trajectory except the steady state
ends on one closed orbit ([](fitzhugh-nagumo.md)). Both oscillate. The
figure below perturbs each of them once.

```{code-cell} python
:tags: [hide-input]
eps, a = 0.05, 0.5
fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.0), gridspec_kw={"height_ratios": [1.3, 1]})
# Lotka-Volterra: halve the prey at t = 8
tA, ZA = run(lv, [2.2, 1.0], 8.0, n=3000)
tB, ZB = run(lv, [ZA[0, -1] * 0.5, ZA[1, -1]], 16.0, n=6000, t0=8.0)
ax, axt = axs[0, 0], axs[1, 0]
flow(ax, lambda x, y: x * (1 - y), lambda x, y: y * (x - 1), (0, 3.6), (0, 3.2), density=0.7)
mark(ax, 1, 1, True)
ax.plot(*ZA, color=C_EXACT, lw=1.2, zorder=4)
ax.plot(*ZB, color=C_PINK, lw=1.2, zorder=4)
ax.annotate("", xy=ZB[:, 0], xytext=ZA[:, -1],
            arrowprops=dict(arrowstyle="-|>", color=C_INNER, lw=1.4), zorder=5)
ax.set_xlabel("prey $x$"); ax.set_ylabel("predators $y$")
ax.set_title("Lotka–Volterra: a center", fontsize=9)
axt.plot(tA, ZA[0], color=C_EXACT, lw=1.1); axt.plot(tB, ZB[0], color=C_PINK, lw=1.1)
axt.axvline(8.0, color=C_INNER, lw=0.8, ls="--")
axt.set_xlabel("time"); axt.set_ylabel("prey $x$"); axt.set_xlim(0, 24)
# FitzHugh-Nagumo: kick the activator up at t = 4.4, late on the lower branch
F = fhn(a, eps)
tA, ZA = run(F, [0.0, 2.0], 4.4, n=4400)
tB, ZB = run(F, [ZA[0, -1], ZA[1, -1] + 1.3], 5.6, n=5600, t0=4.4)
ax, axt = axs[0, 1], axs[1, 1]
flow(ax, lambda w, v: v + a, lambda w, v: (c(v) - w) / eps, (-1.3, 1.3), (-2.5, 2.5), density=0.7)
cubic(ax)
mark(ax, c(-a), -a, False)
ax.plot(*ZA[:, tA > 1.5], color=C_EXACT, lw=1.2, zorder=4)
ax.plot(*ZB, color=C_PINK, lw=1.0, zorder=4)
ax.annotate("", xy=ZB[:, 0], xytext=ZA[:, -1],
            arrowprops=dict(arrowstyle="-|>", color=C_INNER, lw=1.4), zorder=5)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
ax.set_title("FitzHugh–Nagumo: a limit cycle", fontsize=9)
tC, ZC = run(F, ZA[:, -1], 5.6, n=5600, t0=4.4)            # the same run, not kicked
axt.plot(tC, ZC[1], color=C_GREY, lw=1.0, ls=":")
axt.plot(tA, ZA[1], color=C_EXACT, lw=1.1); axt.plot(tB, ZB[1], color=C_PINK, lw=1.1)
axt.axvline(4.4, color=C_INNER, lw=0.8, ls="--")
axt.set_xlabel("time"); axt.set_ylabel("activator $v$"); axt.set_xlim(0, 10)
fig.tight_layout();
```

The kick test. Left: Lotka–Volterra. At $t = 8$ (orange arrow and dashed
line) half of the prey are removed. The population moves to a larger orbit
(pink) and cycles on it indefinitely, with a larger amplitude and a longer
period. Right: FitzHugh–Nagumo with $\ve = 0.05$, $a = 0.5$. At $t = 4.4$,
late in the slow return along the lower branch, the activator is pushed up by
$1.3$, across the middle branch of the S. The next pulse fires immediately
rather than about half a time unit later, and the trajectory (pink) is back on
the same orbit within one lap. Compared with the run that was not kicked
(grey dotted), the later pulses have the same shape and are only shifted
earlier.

The prey cycle retains the effect of the perturbation indefinitely; the
nerve loses it within one pulse. This page identifies the properties of the
equations that determine which of the two occurs. Three tools are used, each
built on the previous one: a conserved quantity, the rate at which the flow
changes area, and a criterion derived from that rate which excludes closed
orbits.

## Conserved quantities

:::{prf:definition} Conserved quantity
:label: def-cq-conserved

A **conserved quantity** of $\dot x = F(x, y)$, $\dot y = G(x, y)$ on a region
$R$ is a continuously differentiable function $H(x, y)$ that is not constant
on any open patch of $R$ and satisfies

$$
\frac{\dd}{\dd t}H\big(x(t), y(t)\big)
= \frac{\partial H}{\partial x}\,F + \frac{\partial H}{\partial y}\,G = 0
\quad\text{in } R .
$$ (eq-cq-conserved)

Every trajectory then lies on one level curve $H = \text{const}$.
:::

The condition "not constant on any patch" excludes the constant function,
which every system conserves.

:::{prf:example} The harmonic oscillator
:label: ex-cq-harmonic
:class: dropdown

A mass on a spring without friction, $\ddot x = -\omega^2x$, is the system
$\dot x = y$, $\dot y = -\omega^2x$ for position $x$ and velocity $y$. Its
energy

$$
H = \tfrac12\,y^2 + \tfrac12\,\omega^2x^2
$$ (eq-cq-energy)

is conserved:

$$
\dot H = y\,\dot y + \omega^2x\,\dot x = -\omega^2xy + \omega^2xy = 0 .
$$

The level curves are ellipses, and every trajectory traverses its ellipse
with the same period $2\pi/\omega$. This is the center (g) of the gallery in
[](linear-systems.md). Its eigenvalues are $\pm i\omega$, the turning of
{prf:ref}`ex-pp-spiral` without stretching.
:::

:::{prf:example} Lotka–Volterra
:label: ex-cq-lotka-volterra
:class: dropdown

In scaled units, prey $x$ grow and are eaten, and predators $y$ die without
prey and grow by eating them:

$$
\dot x = x(1 - y), \qquad \dot y = y(x - 1) .
$$ (eq-cq-lv)

The steady states are $(0, 0)$, a saddle, and $(1, 1)$. Along a trajectory
with $x, y > 0$, $\dd y/\dd x = y(x - 1)/\big(x(1 - y)\big)$ separates into
$(1 - y)\,\dd y/y = (x - 1)\,\dd x/x$, and integrating both sides shows that

$$
H(x, y) = x - \ln x + y - \ln y
$$ (eq-cq-lv-H)

is constant. Directly,

$$
\dot H = \Big(1 - \frac1x\Big)\,x(1 - y) + \Big(1 - \frac1y\Big)\,y(x - 1)
= (x - 1)(1 - y) + (y - 1)(x - 1) = 0 .
$$

At $(1, 1)$ the Jacobian is $\begin{pmatrix} 0 & -1 \\ 1 & 0\end{pmatrix}$,
with $T = 0$ and $\lambda = \pm i$: a linear center. That is a borderline
case, and by the warning at the end of [](linear-systems.md) the nonlinear
terms could turn it into a slow spiral either way. The conserved quantity
decides: see {prf:ref}`ex-cq-lv-closed`.
:::

Where the linearization is a center, a conserved quantity can decide whether
the nonlinear system is one.

:::{prf:proposition} A minimum of a conserved quantity is a center
:label: prop-cq-center

Let $H$ be a conserved quantity and $\mathbf x^\ast$ a steady state, with no
other steady state within a distance $r_0$ of it. Suppose that $H$ increases
along every ray out of $\mathbf x^\ast$: for every unit vector $\mathbf e$,

$$
\frac{\dd}{\dd r}\,H(\mathbf x^\ast + r\,\mathbf e) > 0
\qquad\text{for } 0 < r \le r_0 .
$$ (eq-cq-bowl)

Then every trajectory that starts close enough to $\mathbf x^\ast$, other
than $\mathbf x^\ast$ itself, is a closed orbit that goes once around
$\mathbf x^\ast$. The same holds if $H$ decreases along every ray.

Condition {eq}`eq-cq-bowl` holds, for some $r_0 > 0$, whenever $H$ is twice
continuously differentiable and has a strict minimum at $\mathbf x^\ast$ by
the second-derivative test: $\nabla H(\mathbf x^\ast) = 0$, $H_{xx} > 0$ and
$H_{xx}H_{yy} - H_{xy}^2 > 0$ at $\mathbf x^\ast$.
:::

Condition {eq}`eq-cq-bowl` says that the graph of $H$ is a bowl with its
bottom at $\mathbf x^\ast$. Since $H$ is constant along trajectories, they
can move neither up nor down the bowl; they circle its bottom and return to
their starting point after one turn.

:::{prf:proof}
:class: dropdown

*The second-derivative test.* Write $H''$ for the matrix of second
derivatives at $\mathbf x^\ast$. The two conditions on it say that it is
positive definite, so $\mathbf e^{\mathsf T}H''\,\mathbf e \ge \lambda > 0$
for every unit vector $\mathbf e$, with $\lambda$ its smaller eigenvalue.
Since $\nabla H(\mathbf x^\ast) = 0$, the derivative along the ray is
$r\,\mathbf e^{\mathsf T}H''\,\mathbf e + o(r) \ge r\lambda + o(r)$, which is
positive for $0 < r \le r_0$ once $r_0$ is small enough.

Put $\mathbf x^\ast$ at the origin and use polar coordinates $(r, \theta)$.
Let $D$ be the disc $r \le r_0$, and split the vector field at each point
$\mathbf u \ne 0$ into a radial part $f_r$ and an angular part $f_\theta$,
along $\mathbf e_r$ and $\mathbf e_\theta$. The trajectory turns at the rate
$\dot\theta = f_\theta/r$.

*Trajectories turn.* In polar coordinates the gradient of $H$ is
$\partial_rH\,\mathbf e_r + r^{-1}\partial_\theta H\,\mathbf e_\theta$, so
conservation {eq}`eq-cq-conserved` reads

$$
\frac{\partial H}{\partial r}\,f_r + \frac1r\,\frac{\partial H}{\partial\theta}\,f_\theta = 0 .
$$

If $f_\theta = 0$ at some $\mathbf u \ne 0$ in $D$, then $f_r = 0$ too,
because $\partial H/\partial r > 0$ by {eq}`eq-cq-bowl`; so $\mathbf u$ would
be a steady state, and there is none. Hence $f_\theta \ne 0$ on the punctured
disc. It is continuous there, so it has one sign; say $f_\theta > 0$, so that
every trajectory in $D$ turns counterclockwise.

*An invariant annulus.* By {eq}`eq-cq-bowl`, $H > H(0)$ on the circle
$r = r_0$; let $m > H(0)$ be its smallest value there. Start at $\mathbf p$
with $0 < |\mathbf p| < r_0$ and $h = H(\mathbf p) < m$, which holds for
$\mathbf p$ close enough to $0$. Along the trajectory $H = h$, so it never
reaches the circle, where $H \ge m$. By continuity there is $\delta > 0$ with
$H < h$ on the disc $r < \delta$, so it never enters that disc either. The
trajectory stays in the annulus $\delta \le r \le r_0$, which is closed and
bounded, and on it the continuous positive function $f_\theta/r$ has a
positive minimum $c$. So $\dot\theta \ge c$.

*The orbit closes.* Following the angle continuously along the trajectory,
$\theta(t) \ge \theta(0) + ct$, so at some time $P$ it reaches
$\theta(0) + 2\pi$: the trajectory is back on the ray through $\mathbf p$.
On that ray $H$ is strictly increasing in $r$, so $\mathbf p$ is the only
point of the ray where $H = h$. Hence $\mathbf u(P) = \mathbf p$, and the
trajectory is a closed orbit, which goes once around the origin. For a
decreasing $H$, apply the argument to $-H$.
:::

:::{prf:example} Every Lotka–Volterra orbit is closed
:label: ex-cq-lv-closed
:class: dropdown

For Lotka–Volterra {eq}`eq-cq-lv`, along a ray from $(1, 1)$, with
$x = 1 + re_1$ and $y = 1 + re_2$,

$$
\frac{\dd}{\dd r}\,H = e_1\Big(1 - \frac1x\Big) + e_2\Big(1 - \frac1y\Big)
= r\,\Big(\frac{e_1^2}{x} + \frac{e_2^2}{y}\Big) > 0 ,
$$

and this holds as far as the ray stays in the open quadrant, not only near
$(1, 1)$. The proof of {prf:ref}`prop-cq-center` therefore applies to
every orbit. $H$ grows without bound
towards the axes and towards infinity, so each level set $H = h$ is closed
and bounded and keeps away from $(1, 1)$, which plays the role of the
annulus. So every trajectory with $x, y > 0$ other than $(1, 1)$ is a closed
orbit.
:::

For Lotka–Volterra the bowl condition holds on the whole open quadrant, so
every trajectory with $x, y > 0$ other than $(1, 1)$ is a closed orbit
({prf:ref}`ex-cq-lv-closed`).

```{code-cell} python
:tags: [hide-input]
starts = (1.3, 2.2, 3.5)
cols = (C_COMP, C_OUTER, C_INNER)
fig, (ax, axt) = plt.subplots(1, 2, figsize=(7.4, 3.0),
                              gridspec_kw={"width_ratios": [1, 1.4]})
flow(ax, lambda x, y: x * (1 - y), lambda x, y: y * (x - 1), (0, 4.2), (0, 3.2))
ax.axhline(1, color=C_GREY, lw=0.8, ls=":"); ax.axvline(1, color=C_GREY, lw=0.8, ls=":")
mark(ax, 1, 1, True); mark(ax, 0, 0, False)
for x0, col in zip(starts, cols):
    ev = lambda t, z: z[1] - 1
    ev.direction = 1
    s = solve_ivp(lv, (0, 30), [x0, 1], rtol=1e-11, atol=1e-12, events=ev,
                  dense_output=True)
    P = s.t_events[0][1] - s.t_events[0][0]          # one full lap
    t, (x, y) = run(lv, [x0, 1], 14.0, n=3000)
    ax.plot(*run(lv, [x0, 1], P, n=800)[1], color=col, lw=1.2, zorder=4)
    axt.plot(t, x, color=col, lw=1.2, label=f"$x(0) = {x0}$: period ${P:.2f}$")
ax.set_xlabel("prey $x$"); ax.set_ylabel("predators $y$")
ax.set_title(r"level curves of $H$", fontsize=9)
axt.axhline(1, color=C_GREY, lw=0.8, ls=":")
axt.set_xlabel("time"); axt.set_ylabel("prey $x$"); axt.set_xlim(0, 14)
axt.legend(frameon=False, fontsize=7.5, loc="upper right")
axt.set_title("each orbit keeps its own amplitude and period", fontsize=9)
fig.tight_layout();
```

Lotka–Volterra {eq}`eq-cq-lv`. Left: three closed orbits, level curves of
{eq}`eq-cq-lv-H`, all starting at $y = 1$, and the streamlines (thicker where
the flow is faster). Dotted: the nullclines $x = 1$ and $y = 1$. Right: the
prey on the three orbits. Small orbits are nearly sinusoidal with period
close to $2\pi$, the period of the linearization; larger ones have sharper
peaks and longer periods.

A center therefore has three features, two of which the figure shows.
There is a continuum of closed orbits, one for each value of $H$, and the
initial condition determines which one the system is on. The period depends
on the orbit; only for the harmonic oscillator, which is linear, is it the
same for every amplitude. The third feature is that no orbit attracts its
neighbours. Stating this precisely requires two notions that earlier pages
used without definition.

:::{prf:definition} Invariant set, attractor
:label: def-cq-attractor

A set $S$ of the plane is **forward invariant** if every trajectory that
starts in $S$ stays in $S$ for all $t \ge 0$. It is **invariant** if every
trajectory that starts in $S$ stays in $S$ for all $t$, forward and backward
in time.

A closed and bounded set $A$ is an **attractor** if

1. $A$ is invariant;
2. $A$ attracts an open set $U$ that contains it: every trajectory that
   starts in $U$ approaches $A$, in the sense that its distance to $A$ tends
   to $0$ as $t \to \infty$;
3. $A$ is minimal: no smaller closed set inside $A$ has properties 1 and 2.

A closed orbit that is an attractor is a **stable limit cycle**. On this page,
as in most of biology, *limit cycle* means a stable one.
:::

```{code-cell} python
:tags: [hide-input]
mu = 0.4                                      # r' = mu r (1 - r^2), theta' = 1
cyc = lambda t, z: [mu * z[0] * (1 - z[0]**2 - z[1]**2) - z[1],
                    mu * z[1] * (1 - z[0]**2 - z[1]**2) + z[0]]
TH = np.linspace(0, 2 * np.pi, 300)
ring = lambda r: (r * np.cos(TH), r * np.sin(TH))
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
# forward invariant, not invariant: the disc r <= 1.6
ax = axs[0]
ax.fill(*ring(1.6), color=C_COMP, alpha=0.12, lw=0)
ax.plot(*ring(1.6), color=C_COMP, lw=1.0)
p = [1.5, 0.0]
ax.plot(*run(cyc, p, 12.0, n=2000)[1], color=C_EXACT, lw=1.1, zorder=4)
back = solve_ivp(lambda t, z: [-f for f in cyc(t, z)], (0, 3), p, rtol=1e-9,
                 events=lambda t, z: z[0]**2 + z[1]**2 - 2.2**2, max_step=0.01)
ax.plot(*back.y, color=C_EXACT, lw=1.1, ls="--", zorder=4)
ax.plot(*p, "o", ms=3.5, color=C_EXACT, zorder=5)
ax.text(1.05, -1.75, "$S$", color=C_COMP, fontsize=10)
ax.set_title("forward invariant, not invariant", fontsize=9)
# attractor: the unit circle, attracting the annulus U
ax = axs[1]
ax.fill(*ring(1.6), color=C_OUTER, alpha=0.10, lw=0)
ax.fill(*ring(0.35), color="white", lw=0)
ax.plot(*ring(1.6), color=C_OUTER, lw=0.8, ls=":")
ax.plot(*ring(0.35), color=C_OUTER, lw=0.8, ls=":")
for z0 in ([1.55, 0.0], [0.4, 0.0]):
    ax.plot(*run(cyc, z0, 12.0, n=2000)[1], color=C_GREY, lw=0.9, zorder=3)
ax.plot(*ring(1.0), color=C_EXACT, lw=1.8, zorder=4)
mark(ax, 0, 0, False, ms=4)
ax.text(0.72, 0.72, "$A$", fontsize=10)
ax.text(1.05, -1.75, "$U$", color=C_OUTER, fontsize=10)
ax.set_title("attractor: a stable limit cycle", fontsize=9)
# invariant, not an attractor: closed orbits of Lotka-Volterra, on log scales
ax = axs[2]
Uu, Ss = np.meshgrid(np.linspace(-1.9, 1.9, 300), np.linspace(-1.9, 1.9, 300))
Hl = np.exp(Uu) - Uu + np.exp(Ss) - Ss
ax.contourf(Uu, Ss, Hl, [2.25, 2.75], colors=[C_OUTER], alpha=0.10)
ax.contour(Uu, Ss, Hl, [2.25, 2.75], colors=[C_OUTER], linewidths=0.8, linestyles=":")
ax.contour(Uu, Ss, Hl, [2.33, 2.67], colors=[C_GREY], linewidths=0.9)
ax.contour(Uu, Ss, Hl, [2.5], colors=[C_EXACT], linewidths=1.8)
mark(ax, 0, 0, True, ms=4)
ax.text(0.55, 0.75, r"$\gamma$", fontsize=10)
ax.text(0.75, -1.45, "$U$", color=C_OUTER, fontsize=10)
ax.set_title("invariant, not an attractor", fontsize=9)
for ax, (lo, hi) in zip(axs, ((-2.0, 2.0), (-2.0, 2.0), (-1.8, 1.4))):
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout();
```

Left: in a flow with a limit cycle, the disc $S$ (green) is forward
invariant: every arrow on its edge points inwards, and the trajectory from
the dot stays inside. It is not invariant: followed backward in time
(dashed), the same trajectory leaves. Middle: the circle $A$ is an
attractor. It is invariant, and every trajectory from the annulus $U$ (blue)
winds onto it, from inside and from outside. The closed disc bounded by $A$
is invariant too and also attracts $U$, but it is not minimal, since $A$ is
inside it. Right: closed orbits of Lotka–Volterra, on log scales. The orbit
$\gamma$ is invariant, but the trajectories that start in the band $U$
around it (grey) stay on their own orbits and do not approach it.

The box of {prf:ref}`prop-pp-box` and the hexagon of {prf:ref}`prop-fhn-trap`
are forward invariant, which is the property used there. They are not invariant:
followed backward in time, trajectories leave them. Every trajectory is
invariant, and so is every steady state, every closed orbit, and every level
set of a conserved quantity. A stable node or spiral is an attractor, with
$U$ a small disc around it, and so is the closed orbit of FitzHugh–Nagumo, as
the kick test above and {prf:ref}`rmk-cq-floquet` below show. A closed orbit of a
center is invariant, but it fails property 2: nearby trajectories remain on
their own orbits and do not approach it.

:::{prf:proposition} A conserved quantity rules out attractors
:label: prop-cq-no-attractor

If a system has a conserved quantity on a region $R$, then no steady state
and no closed orbit in $R$ is an attractor whose open set $U$ lies in $R$.
:::

:::{prf:proof}
:class: dropdown

Suppose a steady state $\mathbf x^\ast$ is an attractor, with open set $U$.
Every trajectory that starts in $U$ tends to $\mathbf x^\ast$. $H$ is constant
along it and continuous, so its value at the start is the limit
$H(\mathbf x^\ast)$. Then $H = H(\mathbf x^\ast)$ on all of $U$, and $H$ is
constant on an open patch, which a conserved quantity is not.

For a closed orbit $\gamma$ the argument is the same. $H$ has one value $h$ on
$\gamma$, since $\gamma$ is a trajectory. A trajectory $\mathbf u(t)$ that
starts in $U$ comes within distance $1/n$ of some point $\mathbf q_n$ of
$\gamma$ at some time $t_n$, for every $n$. Since $\gamma$ is closed and
bounded, a subsequence of the $\mathbf q_n$ converges to a point $\mathbf q$ of
$\gamma$, and along it $\mathbf u(t_n) \to \mathbf q$. So
$H(\mathbf u(0)) = H(\mathbf u(t_n)) \to H(\mathbf q) = h$, and $H = h$ on all
of $U$.
:::

Hence a kick moves the system to a neighbouring level set of $H$, where it
remains; this explains the left half of the kick test. A Lotka–Volterra
population that loses part of its prey, for example in a harsh winter,
continues to cycle with a different amplitude.

## How the flow changes area

Finding $H$ relied on a separation of variables, and most systems have no
conserved quantity. A quantity that can always be computed is the rate at
which the flow changes area. Consider a patch of initial conditions carried
by the flow.

```{code-cell} python
:tags: [hide-input]
n = 3000
th = np.linspace(0, 2 * np.pi, n, endpoint=False)
area = lambda x, y: 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def carry(rhs, x0, y0, times, **kw):
    """Carry the closed curve (x0, y0), the edge of a patch, along the flow."""
    s = solve_ivp(lambda t, z: np.concatenate(rhs(z[:n], z[n:])), (0, times[-1]),
                  np.concatenate([x0, y0]), t_eval=times, rtol=1e-10, atol=1e-12, **kw)
    return [(s.y[:n, k], s.y[n:, k]) for k in range(len(times))]


eps, a = 0.05, 0.5
shades = (C_EXACT, "#555555", C_INNER, C_COMP)
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2))
# Lotka-Volterra on log scales, u = ln x, s = ln y
ax = axs[0]
U = np.linspace(-1.6, 1.6, 300)
for h in (2.1, 2.4, 3.0):                      # a few orbits, H = e^u - u + e^s - s
    Uu, Ss = np.meshgrid(U, U)
    ax.contour(Uu, Ss, np.exp(Uu) - Uu + np.exp(Ss) - Ss, [h], colors=[C_GREY],
               linewidths=0.6)
times = [0, 1.0, 2.5, 4.5]
for (x, y), t, col in zip(carry(lambda u, s: (1 - np.exp(s), np.exp(u) - 1),
                                0.55 + 0.25 * np.cos(th), 0.25 * np.sin(th), times),
                          times, shades):
    ax.fill(x, y, color=col, alpha=0.25, lw=0)
    ax.plot(np.append(x, x[0]), np.append(y, y[0]), color=col, lw=0.9)
    ax.plot([], [], color=col, lw=4, alpha=0.5, label=f"$t = {t:g}$: area ${area(x, y):.3f}$")
ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal")
ax.set_xlabel(r"$\ln x$ (prey)"); ax.set_ylabel(r"$\ln y$ (predators)")
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax.set_title("Lotka–Volterra, on log scales", fontsize=9)
# FitzHugh-Nagumo
ax = axs[1]
cubic(ax)
times = [0, 0.03, 0.08, 0.2]
for (x, y), t, col in zip(carry(lambda w, v: (v + a, (c(v) - w) / eps),
                                0.3 + 0.25 * np.cos(th), 1.4 + 0.25 * np.sin(th), times,
                                method="LSODA"),
                          times, shades):
    ax.fill(x, y, color=col, alpha=0.25, lw=0)
    ax.plot(np.append(x, x[0]), np.append(y, y[0]), color=col, lw=0.9)
    ax.plot([], [], color=col, lw=4, alpha=0.5, label=f"$t = {t:g}$: area ${area(x, y):.3f}$")
ax.set_xlim(-1.2, 1.2); ax.set_ylim(-2.4, 2.4)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
ax.legend(frameon=False, fontsize=7, loc="lower left")
ax.set_title(r"FitzHugh–Nagumo, $\varepsilon = 0.05$, $a = 0.5$", fontsize=9)
fig.tight_layout();
```

A disc of starting points carried by the flow, with the area of the patch at
each time. Left: Lotka–Volterra {eq}`eq-cq-lv`, drawn on log scales, with a
few closed orbits in grey. The patch is sheared and bent as it circulates,
but its area does not change. Right: FitzHugh–Nagumo with $\ve = 0.05$,
$a = 0.5$. Within $0.2$ time units the patch is flattened onto the outer
branch of the S, and its area falls by a factor of about $30$.

A patch whose area shrinks is being compressed onto a smaller set, such as
an attractor; a patch whose area is preserved cannot be. To make this
precise, we compute the rate at which the flow changes area.

**Divergence.** Take a tiny rectangle of starting points with corners
$(x, y)$ and $(x + \Delta x, y + \Delta y)$, carried by the flow
$(\dot x, \dot y) = (F, G)$. Its right edge moves horizontally at speed
$F(x + \Delta x, y)$ and its left edge at $F(x, y)$, so
$\frac{\dd}{\dd t}\Delta x \approx F_x\,\Delta x$, and likewise
$\frac{\dd}{\dd t}\Delta y \approx G_y\,\Delta y$. By the product rule, the
area $A = \Delta x\,\Delta y$ changes at the rate

$$
\frac{\dd A}{\dd t} \approx (\nabla\cdot\mathbf f)\,A ,
\qquad
\nabla\cdot\mathbf f = \frac{\partial F}{\partial x} + \frac{\partial G}{\partial y} .
$$ (eq-cq-area)

```{code-cell} python
:tags: [hide-input]
Fx, Gy, dt = 0.8, -0.5, 0.5                  # F = F(x), G = G(y): no shear
F = lambda x: 0.3 + Fx * x
G = lambda y: 0.6 + Gy * y
x0, y0, dx, dy = 0.0, 0.0, 1.0, 1.0
X0, X1 = x0 + F(x0) * dt, x0 + dx + F(x0 + dx) * dt
Y0, Y1 = y0 + G(y0) * dt, y0 + dy + G(y0 + dy) * dt

fig, ax = plt.subplots(figsize=(5.6, 3.1))
ax.fill([x0, x0 + dx, x0 + dx, x0], [y0, y0, y0 + dy, y0 + dy], color=C_GREY, alpha=0.18, lw=0)
ax.plot([x0, x0 + dx, x0 + dx, x0, x0], [y0, y0, y0 + dy, y0 + dy, y0], color=C_EXACT, lw=1.3)
ax.fill([X0, X1, X1, X0], [Y0, Y0, Y1, Y1], color=C_INNER, alpha=0.15, lw=0)
ax.plot([X0, X1, X1, X0, X0], [Y0, Y0, Y1, Y1, Y0], color=C_INNER, lw=1.3, ls="--")
arrow = dict(arrowstyle="-|>", color=C_OUTER, lw=1.3, shrinkA=0, shrinkB=0)
ym = y0 + dy / 2
ax.annotate("", xy=(x0 + F(x0) * dt, ym), xytext=(x0, ym), arrowprops=arrow)
ax.annotate("", xy=(x0 + dx + F(x0 + dx) * dt, ym), xytext=(x0 + dx, ym), arrowprops=arrow)
xm = x0 + dx / 2
ax.annotate("", xy=(xm, y0 + G(y0) * dt), xytext=(xm, y0), arrowprops=arrow)
ax.annotate("", xy=(xm, y0 + dy + G(y0 + dy) * dt), xytext=(xm, y0 + dy), arrowprops=arrow)
ax.text(x0 - 0.04, ym, "$F(x, y)$", color=C_OUTER, fontsize=8, ha="right", va="center")
ax.text(x0 + dx + 0.02, ym + 0.05, r"$F(x + \Delta x, y)$", color=C_OUTER, fontsize=8)
ax.text(xm + 0.03, y0 + 0.08, "$G(x, y)$", color=C_OUTER, fontsize=8)
ax.text(xm, Y1 + 0.1, r"$G(x, y + \Delta y)$", color=C_OUTER, fontsize=8, ha="center")
# side lengths
ax.annotate("", xy=(x0 + dx, -0.22), xytext=(x0, -0.22),
            arrowprops=dict(arrowstyle="<->", color=C_EXACT, lw=0.8))
ax.text(xm, -0.36, r"$\Delta x$", ha="center", fontsize=8.5)
ax.annotate("", xy=(X1, -0.5), xytext=(X0, -0.5),
            arrowprops=dict(arrowstyle="<->", color=C_INNER, lw=0.8))
ax.text((X0 + X1) / 2, -0.64, r"$\Delta x\,(1 + F_x\,\delta t)$", ha="center",
        fontsize=8.5, color=C_INNER)
ax.annotate("", xy=(-0.75, y0 + dy), xytext=(-0.75, y0),
            arrowprops=dict(arrowstyle="<->", color=C_EXACT, lw=0.8))
ax.text(-0.8, ym, r"$\Delta y$", ha="right", va="center", fontsize=8.5)
ax.annotate("", xy=(X1 + 0.12, Y1), xytext=(X1 + 0.12, Y0),
            arrowprops=dict(arrowstyle="<->", color=C_INNER, lw=0.8))
ax.text(X1 + 0.17, (Y0 + Y1) / 2, r"$\Delta y\,(1 + G_y\,\delta t)$", va="center",
        fontsize=8.5, color=C_INNER)
ax.text(x0 + 0.04, y0 + dy - 0.1, "time $t$", fontsize=8)
ax.text(X1 - 0.04, Y0 + 0.05, r"time $t + \delta t$", fontsize=8, color=C_INNER, ha="right")
ax.set_xlim(-1.1, 2.6); ax.set_ylim(-0.75, 1.3); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout();
```

A tiny rectangle of starting points (grey) and where the flow has carried it
a short time $\delta t$ later (orange, dashed). Its right edge moves faster
than its left edge, by $F_x\,\Delta x$, so the rectangle widens by the factor
$1 + F_x\,\delta t$. Its top edge moves slower than its bottom edge, so it
gets shorter by the factor $1 + G_y\,\delta t$, here less than $1$. The area
changes by the product, $1 + (F_x + G_y)\,\delta t$ to first order. In this
sketch the field has no shear; a general field also tilts the rectangle into
a parallelogram, which does not change its area to first order.

The **divergence** $\nabla\cdot\mathbf f$ is the rate at which the flow
stretches area, per unit area. A larger patch is a union of tiny rectangles,
so its area changes at the rate $\iint_\text{patch}\nabla\cdot\mathbf f\,\dd A$.
At a steady state the divergence is the trace of the Jacobian, the sum of the
eigenvalues.

The flows of frictionless mechanics preserve area.

:::{prf:definition} Hamiltonian system
:label: def-cq-hamiltonian

A system is **Hamiltonian**, with **Hamiltonian** $H(x, y)$, if

$$
\dot x = \frac{\partial H}{\partial y},
\qquad
\dot y = -\frac{\partial H}{\partial x} .
$$ (eq-cq-hamiltonian)
:::

Then

$$
\dot H = \frac{\partial H}{\partial x}\,\dot x + \frac{\partial H}{\partial y}\,\dot y
= \frac{\partial H}{\partial x}\,\frac{\partial H}{\partial y}
- \frac{\partial H}{\partial y}\,\frac{\partial H}{\partial x} = 0 :
$$

the Hamiltonian is a conserved quantity. The harmonic oscillator of
{prf:ref}`ex-cq-harmonic` is Hamiltonian with its energy {eq}`eq-cq-energy`,
since $\partial H/\partial y = y$ and $-\partial H/\partial x = -\omega^2x$. This is
the form of every frictionless mechanical system, with $x$ a position, $y$ a
momentum and $H$ the energy.

For a Hamiltonian system the mixed partial derivatives cancel,

$$
\nabla\cdot\mathbf f
= \frac{\partial}{\partial x}\,\frac{\partial H}{\partial y}
- \frac{\partial}{\partial y}\,\frac{\partial H}{\partial x} = 0 ,
$$

so **a Hamiltonian flow preserves area** (Liouville's theorem). A patch of
initial conditions may be sheared and bent, but its area is constant, so it
cannot be compressed onto a point or a curve. This is
{prf:ref}`prop-cq-no-attractor` again, from the point of view of area.

Lotka–Volterra is not Hamiltonian as it stands: its divergence is

$$
\frac{\partial}{\partial x}\,x(1 - y) + \frac{\partial}{\partial y}\,y(x - 1)
= (1 - y) + (x - 1) = x - y .
$$

But dividing the vector field by $xy > 0$ gives

$$
\Big(\frac1y - 1,\ 1 - \frac1x\Big)
= \Big(-\frac{\partial H}{\partial y},\ \frac{\partial H}{\partial x}\Big)
$$

with $H$ from {eq}`eq-cq-lv-H`, which is Hamiltonian with $-H$ in place of
$H$. Dividing by a positive function changes the speed along the
trajectories and not the trajectories themselves. Lotka–Volterra is
therefore a Hamiltonian system up to a change of speed, and the rescaled flow
preserves area.

The original flow preserves a weighted area. The computation behind
{eq}`eq-cq-area`, carried out with area weighted by $B(x, y) > 0$, shows that a patch's weighted area $\iint B\,\dd A$ changes
at the rate $\iint\nabla\cdot(B\mathbf f)\,\dd A$. With $B = 1/(xy)$ the
rescaled field above has divergence zero, so Lotka–Volterra preserves
$\iint \dd x\,\dd y/(xy)$. With $u = \ln x$ and $s = \ln y$, that is
$\iint \dd u\,\dd s$: the ordinary area on log scales, which is what the
figure at the start of this section showed.

:::{prf:remark} Every conserved quantity comes from a rescaling
:label: rmk-cq-every-conserved
:class: dropdown

The rescaling of Lotka–Volterra is an instance of a general fact in the
plane. If $H$ is conserved, {eq}`eq-cq-conserved` says that $(F, G)$ is
perpendicular to $\nabla H = (H_x, H_y)$, so it is parallel to
$(H_y, -H_x)$: $(F, G) = \mu\,(H_y, -H_x)$ for some function $\mu$, wherever
$\nabla H \ne 0$. Where $\mu > 0$, dividing by $\mu$ makes the system
Hamiltonian. A system with a conserved quantity is therefore Hamiltonian up
to the speed at which its trajectories are traversed.
:::

:::{prf:remark} Conserved quantities and symmetries
:label: rmk-cq-noether
:class: dropdown

In mechanics, conserved quantities are tied to symmetries. **Noether's
theorem** states that for a system with a Lagrangian or Hamiltonian structure,
every continuous symmetry gives a conserved quantity: invariance under
shifts in space gives conservation of momentum, invariance under rotations
gives angular momentum, and invariance under shifts in time gives energy.
The Hamiltonian $H$ of {prf:ref}`def-cq-hamiltonian` is of the last kind: it
is conserved because the equations do not depend on $t$.

A general system $\dot x = F$, $\dot y = G$ has no such structure, and
Noether's theorem says nothing about it. In the plane, though,
{prf:ref}`rmk-cq-every-conserved` shows that a conserved quantity supplies
the structure itself, after a change of speed, so every conserved quantity
is the energy of a rescaled system and its symmetry is the shift in time.
On log scales Lotka–Volterra is Hamiltonian without a change of speed: with
$u = \ln x$ and $s = \ln y$, it reads $\dot u = 1 - e^s$, $\dot s = e^u - 1$, which is
Hamiltonian {eq}`eq-cq-hamiltonian` with $-H$, written in $u$ and $s$ as
$-(e^u - u + e^s - s)$. With one position and one momentum, the energy is
all that Noether's theorem can give. It says more for systems with several,
where a symmetry such as a rotation gives a second conserved quantity next
to the energy.
:::

## The criteria of Bendixson and Dulac

The same computation of area yields a criterion that excludes closed orbits.

```{code-cell} python
:tags: [hide-input]
th = np.linspace(0, 2 * np.pi, 400)
rr = 1 + 0.22 * np.cos(2 * th) + 0.08 * np.sin(3 * th)
X, Y = 1.25 * rr * np.cos(th), rr * np.sin(th)
dX, dY = np.gradient(X), np.gradient(Y)                 # counterclockwise tangent
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.fill(X, Y, color=C_COMP, alpha=0.12, lw=0)
ax.plot(X, Y, color=C_EXACT, lw=1.5)
for i in range(15, 400, 50):
    t = np.array([dX[i], dY[i]]) / np.hypot(dX[i], dY[i])
    ax.annotate("", xy=(X[i] + 0.35 * t[0], Y[i] + 0.35 * t[1]), xytext=(X[i], Y[i]),
                arrowprops=dict(arrowstyle="-|>", color=C_OUTER, lw=1.2, shrinkA=0, shrinkB=0))
i = 65
t = np.array([dX[i], dY[i]]) / np.hypot(dX[i], dY[i])
n = np.array([t[1], -t[0]])                               # outward normal
ax.annotate("", xy=(X[i] + 0.45 * n[0], Y[i] + 0.45 * n[1]), xytext=(X[i], Y[i]),
            arrowprops=dict(arrowstyle="-|>", color=C_INNER, lw=1.2, shrinkA=0, shrinkB=0))
ax.text(X[i] + 0.5 * n[0], Y[i] + 0.5 * n[1], r"$\mathbf{n}$", color=C_INNER, fontsize=9)
ax.text(X[i] + 0.45 * t[0] - 0.05, Y[i] + 0.4 * t[1] + 0.06, r"$\mathbf{f}$", color=C_OUTER, fontsize=9)
ax.text(1.55, -1.05, r"$\mathbf{f}\cdot\mathbf{n} = 0$ on $C$", color=C_EXACT, fontsize=8.5)
ax.text(-0.15, -0.1, "$A$", fontsize=12, color=C_COMP)
ax.text(-1.55, 0.95, "$C$", fontsize=11)
ax.set_xlim(-2.0, 2.5); ax.set_ylim(-1.4, 1.5); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout();
```

A closed orbit $C$ and the region $A$ inside it. The vector field
$\mathbf f$ (blue) is tangent to $C$ at every point, so it is perpendicular to
the outward normal $\mathbf n$ (orange): nothing flows into $A$ or out of it
across $C$. The flow moves the points of $A$ around inside $A$, and the area
of $A$ cannot change.

:::{prf:theorem} Dulac's criterion
:label: thm-cq-dulac

Let $(\dot x, \dot y) = (F, G)$ with $F$, $G$ continuously differentiable on a
simply connected region $R$ (one without holes), and let $B(x, y) > 0$ be
continuously differentiable on $R$. If

$$
\nabla\cdot(B\mathbf f) = \frac{\partial (BF)}{\partial x} + \frac{\partial (BG)}{\partial y}
$$

is strictly negative throughout $R$, or strictly positive throughout $R$, then
no closed orbit lies in $R$. With $B = 1$ this is Bendixson's criterion.
:::

:::{prf:proof}
:class: dropdown

Suppose a closed orbit $C$ lies in $R$. It is a simple closed curve, because a
trajectory cannot cross itself; let $A$ be the region it encloses. Because $R$
has no holes, $A$ lies in $R$ too. The divergence theorem in the plane
(Green's theorem) gives

$$
\iint_A \nabla\cdot(B\mathbf f)\,\dd A = \oint_C B\,\mathbf f\cdot\mathbf n\;\dd s ,
$$

with $\mathbf n$ the outward unit normal to $C$ and $s$ arc length. Since $C$
is a trajectory, $\mathbf f$ is tangent to $C$ and $\mathbf f\cdot\mathbf n = 0$,
so the right-hand side is zero. The left-hand side is the integral of a
continuous function of one strict sign over a region of positive area, so it
is not zero. This contradiction proves the claim.
:::

The statements that the area inside $C$ is constant and that there is no
flux across $C$ are equivalent. Multiplying the vector field by $B > 0$ changes the speed along the trajectories
but not the trajectories themselves, as for Lotka–Volterra above, so it does
not change the closed orbits. Equivalently, $B$ is a weight on area, and
Dulac's criterion tests whether the weighted area inside a closed orbit can
be constant. The freedom in choosing $B$ makes it more widely applicable than
Bendixson's criterion.

:::{prf:example} Friction destroys a center
:label: ex-cq-friction
:class: dropdown

Add friction to the spring of {prf:ref}`ex-cq-harmonic`:
$\dot x = y$, $\dot y = -\omega^2x - ky$ with $k > 0$. The divergence is

$$
\nabla\cdot\mathbf f = \frac{\partial y}{\partial x}
+ \frac{\partial}{\partial y}\big(-\omega^2x - ky\big) = -k < 0
$$

everywhere, so there is no closed orbit for any $k > 0$. All closed orbits
of the center disappear, and the energy {eq}`eq-cq-energy` decreases:

$$
\dot H = y\,\dot y + \omega^2x\,\dot x = y\,(-\omega^2x - ky) + \omega^2xy = -ky^2 \le 0 .
$$
:::

Centers are therefore not robust: arbitrarily small dissipation destroys
them.

:::{prf:example} The drug and the bacteria cannot oscillate
:label: ex-cq-drug-dulac
:class: dropdown

For the model {eq}`eq-pp-model` of [](drug-in-the-plane.md), with
$F = rN(1 - N/K) - aNd$ and $G = \delta\,(1 - d - Nd/N_c)$, take
$B = 1/(Nd)$ on the open quadrant $N, d > 0$. Then
$BF = r(1 - N/K)/d - a$ and $BG = \delta\big(1/(Nd) - 1/N - 1/N_c\big)$, and

$$
\frac{\partial (BF)}{\partial N} + \frac{\partial (BG)}{\partial d}
= -\frac{r}{Kd} - \frac{\delta}{Nd^2} < 0 .
$$

So no closed orbit lies in the open quadrant, for any dose and any clearance
rate. Nor can a closed orbit meet the axes: $N = 0$ is made of trajectories, and
every arrow on $d = 0$ points up. With $B = 1$ the divergence is
$r(1 - 2N/K) - ad - \delta(1 + N/N_c)$, which is positive near the origin when
$\delta < r$, and Bendixson alone would not decide. The factor $1/(Nd)$ is
the change to log scales, on which both per-capita rates fall as their own
variable rises.
:::

:::{prf:remark} Dulac's criterion cannot be reversed
:label: rmk-cq-dulac-converse
:class: dropdown

Dulac's criterion says: if $\nabla\cdot(B\mathbf f)$ has one strict sign on
$R$, then $R$ contains no closed orbit. Its contrapositive is logically
equivalent, and so also true: if $R$ contains a closed orbit, then for
every $B > 0$ the divergence $\nabla\cdot(B\mathbf f)$ changes sign in $R$
or vanishes somewhere. [Pump and drain](#cq-pump-drain) below uses it in
this form.

The converse is false: if $R$ contains no closed orbit, there need not be
any $B$ that makes the divergence one-signed. Take the saddle $\dot x = x$,
$\dot y = -y$, which has no closed orbit at all. At the origin $\mathbf f = 0$,
so $\nabla\cdot(B\mathbf f) = B\,(\nabla\cdot\mathbf f) + \nabla B\cdot\mathbf f = 0$
there, whatever $B$ is. So a divergence that changes sign or vanishes, for
every choice of $B$, still does not give a closed orbit. In terms of the
proof: a closed orbit forces the integral of $\nabla\cdot(B\mathbf f)$ over
the region inside it to vanish, but that integral also vanishes over many
regions whose edges are not trajectories.
:::

Dulac's criterion can only exclude closed orbits
({prf:ref}`rmk-cq-dulac-converse`). To show that a closed orbit exists takes
a different kind of argument, one that traps a trajectory: the
Poincaré–Bendixson theorem {prf:ref}`thm-fhn-pb`, by which a trapping region
that contains no steady state contains a closed orbit. That is how
[](fitzhugh-nagumo.md) found its orbit. The full theorem is stated at the end
of this page ({prf:ref}`thm-cq-pb`).

(cq-pump-drain)=
## Pump and drain: what a limit cycle needs

By the proof of Dulac's criterion, the region $A$ inside any closed orbit
satisfies $\iint_A\nabla\cdot(B\mathbf f)\,\dd A = 0$ for every $B > 0$, so
$\nabla\cdot(B\mathbf f)$ either takes both signs in $A$ or vanishes
identically there. For a system with a conserved quantity the second
alternative holds for a suitable $B$, such as $B = 1/(xy)$ for
Lotka–Volterra. For an isolated attracting orbit the first holds: inside it
there is a region where the flow expands area, a **pump**, and a region
where it contracts area, a **drain**.

For FitzHugh–Nagumo {eq}`eq-fhn`, with $(\dot w, \dot v) = \big(v + a,\ (c(v) - w)/\ve\big)$,

$$
\nabla\cdot\mathbf f = \frac{\partial (v + a)}{\partial w}
+ \frac{\partial}{\partial v}\,\frac{c(v) - w}{\ve} = \frac{1 - v^2}{\ve} .
$$ (eq-cq-fhn-div)

It is positive in the strip $|v| < 1$, where the activator amplifies itself,
and negative outside it, where the amplification saturates. Each of the
three regions has no holes, so no closed orbit lies inside any one of them:
the closed orbit of {prf:ref}`prop-fhn-closed-orbit` must cross both lines
$v = \pm1$.

```{code-cell} python
:tags: [hide-input]
eps, a = 0.05, 0.5
t, (w, v) = run(fhn(a, eps), [0.0, 2.0], 20.0, n=40000)
lap = t > 15
fig, (ax, axt) = plt.subplots(1, 2, figsize=(7.4, 3.0),
                              gridspec_kw={"width_ratios": [1, 1.3]})
ax.axhspan(-1, 1, color=C_INNER, alpha=0.10, lw=0)
for lo, hi in ((1, 3), (-3, -1)):
    ax.axhspan(lo, hi, color=C_OUTER, alpha=0.08, lw=0)
cubic(ax)
ax.plot(w[lap], v[lap], color=C_EXACT, lw=1.3, zorder=4)
mark(ax, c(-a), -a, False)
ax.text(-1.15, 0.05, "pump\n$\\nabla\\cdot\\mathbf{f} > 0$", fontsize=7.5, color=C_INNER)
ax.text(-1.15, 1.6, "drain", fontsize=7.5, color=C_OUTER)
ax.text(0.6, -1.9, "drain", fontsize=7.5, color=C_OUTER)
ax.set_xlim(-1.2, 1.2); ax.set_ylim(-2.4, 2.4)
ax.set_xlabel("inhibitor $w$"); ax.set_ylabel("activator $v$")
ax.set_title("the closed orbit crosses both regions", fontsize=9)

div = (1 - v**2) / eps
cum = np.concatenate([[0], np.cumsum(0.5 * (div[1:] + div[:-1]) * np.diff(t))])
tl, cl = t[lap] - t[lap][0], cum[lap] - cum[lap][0]
axt.plot(tl, cl, color=C_EXACT, lw=1.3)
axt.set_xlabel("time"); axt.set_ylabel(r"$\int_0^t \nabla\cdot\mathbf{f}\,\mathrm{d}s$")
axt.set_xlim(0, tl[-1])
axt.set_title("accumulated divergence along the orbit", fontsize=9)
fig.tight_layout();
```

FitzHugh–Nagumo with $\ve = 0.05$ and $a = 0.5$. Left: the closed orbit
(black) around the S (blue; dashed where repelling), with the strip
$|v| < 1$, where the divergence {eq}`eq-cq-fhn-div` is positive, shaded
orange, and the drain outside it shaded blue. Right: the divergence
integrated along the orbit, over about two laps. It rises briefly during each
jump across the strip and falls steadily while the orbit moves slowly along the
outer branches of the S.

:::{prf:remark} Net drain makes a limit cycle attract
:label: rmk-cq-floquet
:class: dropdown

For a closed orbit $\gamma$ of period $P$, a small displacement across
$\gamma$ is multiplied, after one lap, by

$$
\exp\Big(\int_0^P \nabla\cdot\mathbf f\big(\gamma(t)\big)\,\dd t\Big) .
$$

If the integral is negative, every nearby trajectory is drawn onto $\gamma$:
it is an attractor, a stable limit cycle ({prf:ref}`def-cq-attractor`). We
state this without proof {cite:p}`perko2001`. For FitzHugh–Nagumo with $\ve = 0.05$ and $a = 0.5$ the
orbit spends about $14\%$ of its time in the pump strip, where the integral
gains about $4$, and the rest on the outer branches, where it loses about
$81$. The factor per lap is $e^{-77}$, so the effect of a kick has decayed by the
time the orbit returns to the S. For the rescaled Lotka–Volterra flow the
divergence is zero, so the factor is $1$ on every orbit: perturbations
neither decay nor grow.
:::

This explains the right half of the kick test. The nerve returns to its
cycle because, on every lap, its flow contracts area far more than it
expands it; the prey cycle does not, because its flow preserves area on log
scales.

:::{prf:remark} The Hopf bifurcation
:label: rmk-cq-hopf
:class: dropdown

Close to the parameter value at which it is created, a limit cycle resembles
a harmonic oscillator. At $a = \pm1$ the steady state of FitzHugh–Nagumo sits
at a fold, and by {prf:ref}`prop-fhn-stability` its eigenvalues are
$\pm i/\sqrt\ve$: the linearization is a harmonic oscillator, with frequency
$\omega = 1/\sqrt\ve$. For $|a|$ slightly below $1$ the trace
$(1 - a^2)/\ve$ is small and positive, and the steady state is a slowly
growing spiral. The nonlinear terms stop the growth. Differentiating the fast
equation of {eq}`eq-fhn` and using $\dot w = v + a$ gives

$$
\ve\,\ddot v + (v^2 - 1)\,\dot v + (v + a) = 0 ,
$$ (eq-cq-vdp)

a mass $\ve$ on a spring with a friction coefficient $v^2 - 1$ that is
negative for $|v| < 1$, where it pumps energy in, and positive for
$|v| > 1$, where it drains it: the two signs of {eq}`eq-cq-fhn-div`. The
oscillation grows until it reaches the drain and settles where the two
balance over a lap, on a small, nearly elliptical limit cycle with period
close to $2\pi\sqrt\ve$. (With $a = 0$, {eq}`eq-cq-vdp` is the **van der
Pol oscillator**, from which FitzHugh built his model.) The birth of a limit
cycle from a steady state, as a pair of complex eigenvalues crosses the
imaginary axis, is a **Hopf bifurcation** (Strogatz §8.2
{cite:p}`strogatz2015`). For small $\ve$ the small cycle exists only in a
narrow interval, about $0.993 < a < 1$ for $\ve = 0.05$, and then grows into
the relaxation oscillation of [](fitzhugh-nagumo.md), which is far from
harmonic.
:::

## Summary

Centers and limit cycles are both periodic, but they differ in every other
respect. A center comes from a conserved quantity. Its closed orbits form
a continuum, one for each level of $H$, the amplitude and the period are set
by the initial condition, and a kick moves it to another orbit permanently.
Its flow keeps area, after a change of speed, so it attracts nothing, and
arbitrarily small friction destroys all of its orbits: crowding of the
prey turns the Lotka–Volterra center into a stable spiral (see the
[exercises](exercises.md)). A limit cycle is isolated. Its amplitude and
period are set by the parameters, a kick only shifts its phase, and around
the orbit its flow pumps area in some places and drains it in others,
draining more on balance. A rhythm whose amplitude is restored after any
perturbation, such as a nerve that fires a train of identical pulses, must
be a limit cycle.

The Poincaré–Bendixson theorem classifies the long-time behaviour of a
bounded trajectory in the plane. The **limit set** of a trajectory
$\mathbf u(t)$ is the set of points $\mathbf q$ with
$\mathbf u(t_n) \to \mathbf q$ for some sequence of times $t_n \to \infty$.

:::{prf:theorem} Poincaré–Bendixson
:label: thm-cq-pb

Let $F$ and $G$ be continuously differentiable, and let $\mathbf u(t)$ be a
trajectory of $\dot x = F$, $\dot y = G$ that stays in a closed, bounded
region $R$ for all $t \ge 0$, where $R$ contains only finitely many steady
states. Then the limit set of $\mathbf u(t)$ is one of the following:

1. a steady state;
2. a closed orbit;
3. finitely many steady states together with trajectories that connect them,
   each of which tends to one of the steady states as $t \to \infty$ and to
   one as $t \to -\infty$.
:::

For a proof see Perko §3.7 {cite:p}`perko2001`.
