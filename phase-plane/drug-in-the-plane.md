---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/drug-in-the-plane.pdf
    id: phase-plane-drug-in-the-plane-pdf
downloads:
  - id: phase-plane-drug-in-the-plane-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Drug and Bacteria in the Plane

:::{tip} Big Idea
In [](../one-dimensional/saddle-node.md) the drug was replaced by its steady
state, leaving one equation for the cells. Here we keep both variables, so the
state is a point in a **phase plane**, organized by nullclines and classified
by the trace and determinant of the Jacobian. The plane reproduces every
verdict of the phase line, and it shows why the drug could be eliminated: the
drug is pulled toward its nullcline at a rate never below the clearance rate.
What the reduction could not see is that the threshold inoculum is a curve,
the separatrix of a saddle, and that this curve is tilted.
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

r, K, Nc = 1.0, 10.0, 1.0           # 1/day, 1e8 cells/mL, 1e8 cells/mL
dbar = lambda N: 1 / (1 + N / Nc)   # drug nullcline (the QSS drug level)
NN = np.linspace(0, K, 400)


def field(N, d, a, delta):
    return r * N * (1 - N / K) - a * N * d, delta * (1 - d - N * d / Nc)


def rhs(t, y, a, delta):
    return field(y[0], y[1], a, delta)


def interior_states(a):
    z = np.roots([-r / (K * Nc), r * (1 / Nc - 1 / K), r - a])
    return sorted(x.real for x in z if abs(x.imag) < 1e-12 and x.real > 0)


def jacobian(N, d, a, delta):
    return np.array([[r * (1 - 2 * N / K) - a * d, -a * N],
                     [-delta * d / Nc, -delta * (1 + N / Nc)]])


def separatrix(a, delta):
    """Stable manifold of the saddle (N_-, dbar(N_-)), integrated backward
    until it leaves the box 0 <= d <= 1."""
    Nm = interior_states(a)[0]
    w, V = np.linalg.eig(jacobian(Nm, dbar(Nm), a, delta))
    v = V[:, np.argmin(w.real)].real   # numpy >= 2.5 returns complex arrays

    def leave(t, y):
        return min(y[1] + 0.02, 1.02 - y[1], y[0])
    leave.terminal = True
    back = lambda t, y: [-z for z in field(y[0], y[1], a, delta)]
    halves = [solve_ivp(back, (0, 60), [Nm, dbar(Nm)] + 1e-7 * s * v,
                        events=leave, rtol=1e-10, atol=1e-12, method="LSODA").y
              for s in (1, -1)]
    return np.hstack([halves[0][:, ::-1], halves[1]])


def mark(ax, x, y, stable, ms=6):
    style = dict(color=C_COMP) if stable else dict(mfc="white", mec=C_EXACT)
    ax.plot([x], [y], "o", ms=ms, zorder=6, clip_on=False, **style)


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

## Keeping the drug

[](../one-dimensional/saddle-node.md) started from two balance equations, the
logistic cells with a kill term and the drug equation {eq}`eq-bf-drug`, and
then let the drug sit at its steady state {eq}`eq-bf-drug-qss`. Here we do not
take that step.

:::{note} Biology background: half-life and clearance
A drug in the blood is removed by the kidneys and the liver. If each molecule
has the same chance per unit time $\delta$ of being removed, the drug level
decays like $e^{-\delta t}$ and halves in the **half-life**
$t_{1/2} = \ln 2/\delta$. Meropenem, the β-lactam used against TB in
[](../one-dimensional/saddle-node.md), has a half-life of about one hour, so
$\delta \approx 0.69$ per hour $\approx 17$ per day. Other antibiotics stay
in the body much longer: a half-life of two days gives $\delta \approx 1/3$
per day, slower than a TB cell divides. Below, $\delta$ is a dial between
these two cases.
:::

Measure the drug in units of its drug-free level $D_0 = \sigma/\delta$, so
that $d = D/D_0$ lies between $0$ and $1$. Dividing {eq}`eq-bf-drug` by $D_0$
and using $N_c = \delta/\gamma$ gives the following model.

:::{prf:definition} Collective resistance with the drug kept
:label: def-pp-model

$$
\dot N = rN\Big(1 - \frac{N}{K}\Big) - a\,N\,d,
\qquad
\dot d = \delta\Big(1 - d - \frac{N\,d}{N_c}\Big),
$$ (eq-pp-model)

with $N \ge 0$ in $10^8$ cells per mL, $0 \le d \le 1$, time in days, and
$a = kD_0$ the kill rate at the full dose.
:::

All parameters except $\delta$ are those of
[](../one-dimensional/saddle-node.md): $r = 1$/day, $K = 10$, $N_c = 1$, and
the dose $a = 2$/day. The clearance rate $\delta$ multiplies every term of the
drug equation, including the supply, because the unit $D_0$ itself depends on
$\delta$: in these units the drug-free level is approached in a time
$1/\delta$.

## The vector field

The state of the infection at time $t$ is the point $(N(t), d(t))$ in the
plane. At every point, {eq}`eq-pp-model` gives a velocity $(\dot N, \dot d)$,
an arrow that says where the state moves next. The arrows together form the
**vector field**, and a solution is a curve that follows the arrows, a
**trajectory**.

:::{prf:example} Three arrows
:label: ex-pp-arrows

Take $a = 2$, $\delta = 2$. At $(1, 1)$: $\dot N = 0.9 - 2 = -1.1$ and
$\dot d = 2(1 - 1 - 1) = -2$, down and left. At $(5, 0.1)$: $\dot N = 1.5$,
$\dot d = 0.8$, up and right. At $(8, 0.5)$: $\dot N = -6.4$, $\dot d = -7$,
steeply down and left.
:::

```{code-cell} python
:tags: [hide-input]
a, delta = 2.0, 2.0
F = lambda N, d: r * N * (1 - N / K) - a * N * d
G = lambda N, d: delta * (1 - d - N * d / Nc)
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9), sharey=True)
ax = axs[0]
gN, gd = np.meshgrid(np.linspace(0.5, 9.5, 10), np.linspace(0.05, 0.95, 10))


def unit(u, v, L):
    """Arrow of length L (in axis fractions) along (u, v) on the box [0,K]x[0,1]."""
    n = np.hypot(u / K, v)
    return L * K * (u / K) / n, L * v / n


ax.quiver(gN, gd, *unit(F(gN, gd), G(gN, gd), 0.06), color="0.6", angles="xy",
          scale_units="xy", scale=1, width=0.005)
for N0, d0 in ((1, 1), (5, 0.1), (8, 0.5)):
    ax.quiver([N0], [d0], *unit(F(N0, d0), G(N0, d0), 0.12), color=C_INNER,
              angles="xy", scale_units="xy", scale=1, width=0.009, zorder=5)
    ax.plot([N0], [d0], "o", color=C_INNER, ms=3, zorder=6, clip_on=False)
ax.set_xlim(0, K); ax.set_ylim(0, 1)
ax.set_xlabel(r"$N$ ($10^8$ cells/mL)"); ax.set_ylabel(r"drug $d = D/D_0$")
ax.set_title("the vector field (directions)")
ax = axs[1]
phase_plane(ax, F, G, (0, K), (0, 1), trajs=[(1, 1), (5, 0.1), (8, 0.5)],
            null=False, fixed_points=False, t_max=8)
ax.set_xlabel(r"$N$ ($10^8$ cells/mL)")
ax.set_title("trajectories follow the arrows")
fig.tight_layout();
```

The vector field of {eq}`eq-pp-model` with $a = 2$/day, $\delta = 2$/day.
Left: the direction of the arrow at a grid of points (grey), and the three
arrows of {prf:ref}`ex-pp-arrows` (orange). Right: streamlines of the same
field, thicker where the state moves faster, and the trajectories from the
three points (squares).

Two facts carry over from the phase line. First, **trajectories do not
cross**: through each point passes exactly one trajectory, because the initial
value problem has exactly one solution. On the line this meant that a solution
cannot pass a steady state; in the plane, trajectories can squeeze past each
other but never cut through. Second, the biologically meaningful states form
a trap.

:::{prf:proposition} The box is invariant
:label: prop-pp-box

A solution of {eq}`eq-pp-model` that starts in the box
$0 \le N \le K$, $0 \le d \le 1$ stays in it for all $t \ge 0$.
:::

:::{prf:proof}
:class: dropdown

On each edge the arrows point inward or along the edge: at $N = 0$,
$\dot N = 0$; at $N = K$, $\dot N = -aKd \le 0$; at $d = 0$,
$\dot d = \delta > 0$; at $d = 1$, $\dot d = -\delta N/N_c \le 0$. So no
solution can leave through an edge.
:::

Cells and drug stay nonnegative, the drug never exceeds its drug-free level,
and the cells never exceed the carrying capacity.

## Nullclines

The phase line was organized by its steady states. The phase plane is
organized by the curves on which one of the two velocities vanishes.

:::{prf:definition} Nullclines
:label: def-pp-nullclines

The **$N$-nullcline** is the set where $\dot N = 0$, and the
**$d$-nullcline** the set where $\dot d = 0$. On the $N$-nullcline the arrows
are vertical; on the $d$-nullcline they are horizontal. A steady state is a
point where both vanish: a crossing of the two nullclines.
:::

For {eq}`eq-pp-model` the $N$-nullcline is
$N\big(r(1 - N/K) - a\,d\big) = 0$. It has two pieces: the axis $N = 0$, and
the straight line

$$
d = \frac{r}{a}\Big(1 - \frac{N}{K}\Big),
$$ (eq-pp-nnull)

from height $r/a$ on the $d$ axis down to $(K, 0)$. Below the line the drug is
too weak and the cells grow; above it they shrink. The $d$-nullcline is

$$
d = \bar d(N) = \frac{1}{1 + N/N_c},
$$ (eq-pp-dnull)

the quasi-steady drug level {eq}`eq-bf-drug-qss`. Above it the drug falls,
below it the drug rises.

:::{prf:remark} What the clearance rate changes
:label: rmk-pp-delta
:class: dropdown

Neither nullcline depends on $\delta$, because $\delta$ multiplies the whole
drug equation and so cannot move the set where $\dot d = 0$. The steady states,
which sit on the nullclines, therefore do not depend on $\delta$ either. What
$\delta$ changes is the length of the vertical component of every arrow, and so
the direction of the arrows away from the nullclines and the paths that
trajectories take.
:::

**Steady states.** The axis $N = 0$ meets the $d$-nullcline at $(0, 1)$: no
cells, full drug. The line {eq}`eq-pp-nnull` meets it where

$$
\frac{r}{a}\Big(1 - \frac{N}{K}\Big) = \frac{1}{1 + N/N_c}
\quad\Longleftrightarrow\quad
r\Big(1 - \frac{N}{K}\Big)\Big(1 + \frac{N}{N_c}\Big) = a
\quad\Longleftrightarrow\quad
h(N) = a ,
$$

with $h$ the tolerated dose {eq}`eq-bf-tolerance`. The graphical method of
[](../one-dimensional/saddle-node.md), a horizontal line at height $a$ crossing
the hump of $h$, is this picture in disguise. Both bifurcations are visible:

- **Transcritical at $a = r$.** The line {eq}`eq-pp-nnull` starts at height
  $r/a$ on the $d$ axis. At $a = r$ it passes through $(0, 1)$, and a crossing
  moves from $N > 0$ onto the axis.
- **Saddle-node at $a = a_\mathrm{sn}$.** Raising $a$ tilts the line down,
  pivoting about $(K, 0)$. At $a_\mathrm{sn}$ it touches the hyperbola
  $\bar d$, and the two crossings $N_\pm$ merge and disappear.

```{code-cell} python
:tags: [hide-input]
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.8), sharey=True)
for ax, (a, title) in zip(axs, [(0.5, r"$a = 0.5 < r$"),
                                (2.0, r"$r < a = 2 < a_\mathrm{sn}$"),
                                (3.5, r"$a = 3.5 > a_\mathrm{sn}$")]):
    phase_plane(ax, lambda N, d: r * N * (1 - N / K) - a * N * d,
                lambda N, d: 2.0 * (1 - d - N * d / Nc), (0, K), (0, 1))
    ax.plot([0, 0], [0, 1], color=C_INNER, lw=2.5, zorder=3, clip_on=False)
    ax.set_title(title); ax.set_xlabel(r"$N$ ($10^8$ cells/mL)")
axs[0].plot([], [], color=C_OUTER, lw=2, label=r"$\dot d = 0$")
axs[0].plot([], [], color=C_INNER, lw=2, label=r"$\dot N = 0$")
axs[0].set_ylabel(r"drug $d = D/D_0$")
axs[0].legend(loc="upper right", frameon=False, fontsize=7.5)
fig.tight_layout();
```

Nullclines of {eq}`eq-pp-model` for three doses. Blue: the $d$-nullcline
$\bar d(N)$. Orange: the $N$-nullcline, the axis $N = 0$ and the line
{eq}`eq-pp-nnull`. Pale lines are streamlines with $\delta = 2$/day. Filled
circles are stable steady states, open circles unstable (classified
numerically here, and by hand in {prf:ref}`prop-pp-agree`). Left: one positive
steady state. Middle: two, the threshold $N_-$ and the established infection
$N_+$. Right: the line has tilted below the hyperbola and only the cleared
state $(0, 1)$ is left.

In the middle panel the nullclines cut the box into regions in each of which
$\dot N$ and $\dot d$ have fixed signs. For instance, a trajectory starting at
$(3, 1)$ first loses drug, then the cells grow, and it settles at $N_+$; one
starting at $(0.5, 1)$ is killed while the drug stays high, and goes to
$(0, 1)$.

## Linearization in the plane

On the phase line the slope $f'(N^\ast)$ decided stability
({prf:ref}`prop-ir-slope-test`). In the plane, write $N = N^\ast + x$,
$d = d^\ast + y$ with $x, y$ small, and expand the right-hand sides
$F(N, d)$, $G(N, d)$ of {eq}`eq-pp-model` to first order:

$$
\begin{pmatrix}\dot x\\ \dot y\end{pmatrix}
\approx J \begin{pmatrix}x\\ y\end{pmatrix},
\qquad
J = \begin{pmatrix} F_N & F_d \\ G_N & G_d \end{pmatrix},
$$ (eq-pp-linear)

with every partial derivative taken at the steady state. $J$ is the
**Jacobian matrix**. Its entries say how each rate responds to each variable:
$F_d$ is how much faster the cells die per unit of extra drug, $G_N$ how much
faster the drug falls per extra cell.

This is the linear system of [](linear-systems.md). Its solutions are
combinations of
$e^{\lambda t}\mathbf v$ with $J\mathbf v = \lambda\mathbf v$
({prf:ref}`prop-pp-eigen-solutions`): along the eigenvector $\mathbf v$ the
deviation grows or shrinks like $e^{\lambda t}$, as on a phase line with slope
$\lambda$. The type of the steady state is read off from the trace
$T = F_N + G_d$ and the determinant $\Delta = F_N G_d - F_d G_N$
({prf:ref}`prop-pp-trace-det`).

Differentiating {eq}`eq-pp-model`, the Jacobian of the model is

$$
J = \begin{pmatrix}
r\big(1 - \tfrac{2N}{K}\big) - a\,d & -a\,N \\[3pt]
-\delta\,d/N_c & -\delta\,\big(1 + N/N_c\big)
\end{pmatrix}.
$$ (eq-pp-jacobian)

:::{prf:proposition} The plane agrees with the phase line
:label: prop-pp-agree

For every clearance rate $\delta > 0$:

1. $(0, 1)$ is a stable node if $a > r$ and a saddle if $a < r$; its
   eigenvalues are $r - a$ and $-\delta$.
2. At a positive steady state $(N^\ast, \bar d(N^\ast))$,

   $$
   T = -\frac{rN^\ast}{K} - \delta\Big(1 + \frac{N^\ast}{N_c}\Big) < 0,
   \qquad
   \Delta = -\,\delta\,N^\ast\,h'(N^\ast) .
   $$ (eq-pp-det)

3. Hence a positive steady state on the rising side of $h$
   ($h'(N^\ast) > 0$, the threshold $N_-$) is a saddle, and one on the falling
   side ($h'(N^\ast) < 0$, the established state $N_+$) is a stable node.

These are exactly the verdicts of the phase line of
[](../one-dimensional/saddle-node.md).
:::

:::{prf:proof}
:class: dropdown

At $(0, 1)$ the top-right entry $-aN$ of {eq}`eq-pp-jacobian` vanishes, so $J$
is triangular and its eigenvalues are the diagonal entries $r - a$ and
$-\delta$. The first is the slope $f'(0) = r - a$ of the phase line; the
second has eigenvector $(0, 1)$, the drug direction. Both are real, so the
point is a node when both are negative and a saddle when $r - a > 0$.

At a positive steady state, $r(1 - N^\ast/K) = a\,d^\ast$, so the top-left
entry of $J$ is $-rN^\ast/K$, which gives $T$. For the determinant, use
$a\,d^\ast/N_c = r(1 - N^\ast/K)/N_c$:

$$
\Delta
= \frac{rN^\ast}{K}\,\delta\Big(1 + \frac{N^\ast}{N_c}\Big)
- \frac{a\,\delta\,N^\ast d^\ast}{N_c}
= \delta N^\ast\Big[\frac{r}{K}\Big(1 + \frac{N^\ast}{N_c}\Big)
- \frac{r}{N_c}\Big(1 - \frac{N^\ast}{K}\Big)\Big]
= -\,\delta\,N^\ast\,h'(N^\ast),
$$

since the bracket is $-h'(N^\ast)$ for
$h(N) = r(1 - N/K)(1 + N/N_c)$. Now apply {prf:ref}`prop-pp-trace-det`:
$\Delta < 0$ gives a saddle; $\Delta > 0$ with $T < 0$ gives a stable steady
state, which is a node and not a spiral by {prf:ref}`prop-pp-no-spirals`.
:::

:::{prf:remark} The determinant at a bifurcation
:label: rmk-pp-det-bifurcation
:class: dropdown

At $a_\mathrm{sn}$ the two crossings merge at the top of the hump, where
$h' = 0$, so $\Delta = 0$ and one eigenvalue is zero. This is the planar form
of "the slope test fails at a bifurcation". At the transcritical bifurcation
$a = r$, the eigenvalue $r - a$ of $(0, 1)$ is zero.
:::

### No spirals

:::{prf:proposition} Same-sign cross effects give real eigenvalues
:label: prop-pp-no-spirals

If the off-diagonal entries of a $2\times 2$ matrix $J$ satisfy
$F_d\,G_N > 0$, then both eigenvalues of $J$ are real. In
{eq}`eq-pp-jacobian`, $F_d = -aN < 0$ and $G_N = -\delta d/N_c < 0$, so every
positive steady state of {eq}`eq-pp-model` is a saddle or a node, never a
spiral.
:::

:::{prf:proof}
:class: dropdown

$T^2 - 4\Delta = (F_N + G_d)^2 - 4F_NG_d + 4F_dG_N = (F_N - G_d)^2 + 4F_dG_N > 0$.
:::

In words: the drug harms the cells, and the cells harm the drug. The two
**cross effects** have the same sign, and solutions near a stable state slide
in along nearly straight paths without circling it. To circle, one variable
must promote the other while the other suppresses the first.

## The saddle and its separatrix

At $a = 2$ the saddle sits at $N_- \approx 1.30$. Along one of its
eigenvectors solutions leave it, toward the cleared state on one side and
toward $N_+$ on the other. Along the other eigenvector solutions arrive.

:::{prf:definition} Stable manifold and separatrix
:label: def-pp-separatrix

The **stable manifold** of a saddle is the set of starting points whose
trajectories tend to the saddle. It consists of the saddle and two
trajectories that arrive at it, one from each side, tangent to the
eigenvector of the negative eigenvalue. When it divides the plane into
regions with different fates, it is called a **separatrix**.
:::

Every start on one side of the separatrix is cleared; every start on the other
side settles at $N_+$. A start exactly on it approaches the saddle and
hesitates there forever, the planar version of starting exactly at $N_-$.

So **the threshold is a curve**. On the phase line it was one number, $N_-$.
In the plane, whether an inoculum is cleared depends on the starting drug level
as well as on the starting number of cells. Two protocols differ exactly in
$d(0)$:

- cells arrive in a patient who already has the full dose: $d(0) = 1$;
- treatment starts on an infection that is already there: $d(0) = 0$.

```{code-cell} python
:tags: [hide-input]
a = 2.0
Nm, Np = interior_states(a)
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
for ax, (delta, title) in zip(axs, [(17.0, r"TB and meropenem, $\delta = 17$/day"),
                                    (1 / 3, r"a slow drug, $\delta = 1/3$/day")]):
    ax.plot(NN, dbar(NN), color=C_OUTER, lw=1.4, ls="--", label=r"$\dot d = 0$")
    ax.plot(NN, (r / a) * (1 - NN / K), color=C_INNER, lw=1.0, alpha=0.6,
            label=r"$\dot N = 0$")
    sep = separatrix(a, delta)
    ax.fill_betweenx(np.linspace(0, 1, 200), 0,
                     np.interp(np.linspace(0, 1, 200), *sep[::-1][:, np.argsort(sep[1])]),
                     color=C_INNER, alpha=0.08, lw=0)
    ax.plot(sep[0], sep[1], color=C_INNER, lw=2, label="separatrix")
    for N0 in (0.4, 1.0, 2.0, 4.0, 7.0, 9.5):
        for d0 in (0.0, 1.0):
            sol = solve_ivp(rhs, (0, 40), [N0, d0], args=(a, delta), rtol=1e-9,
                            atol=1e-12, method="LSODA", dense_output=True)
            t = np.concatenate([np.linspace(0, 6 / delta, 200), np.linspace(6 / delta, 40, 400)])
            n, d = sol.sol(t)
            ax.plot(n, d, color=C_EXACT, lw=0.8)
            ax.plot([N0], [d0], "s", color=C_EXACT, ms=2.5)
            i = 40
            ax.annotate("", xy=(n[i + 3], d[i + 3]), xytext=(n[i], d[i]),
                        arrowprops=dict(arrowstyle="-|>", color=C_EXACT, lw=0.8))
    mark(ax, 0, 1, True); mark(ax, Nm, dbar(Nm), False); mark(ax, Np, dbar(Np), True)
    ax.set_xlim(-0.1, K); ax.set_ylim(-0.02, 1.04)
    ax.set_xlabel(r"$N$ ($10^8$ cells/mL)"); ax.set_title(title)
axs[0].set_ylabel(r"drug $d = D/D_0$")
axs[0].text(0.15, 0.06, "cleared", color=C_INNER, fontsize=8)
axs[1].legend(loc="upper right", frameon=False, fontsize=7.5)
fig.tight_layout();
```

Phase portraits of {eq}`eq-pp-model` at $a = 2$/day. Trajectories start at
the squares, with no drug ($d = 0$) or the full dose ($d = 1$). Inocula left
of the separatrix (orange, shaded) are cleared; the rest settle at $N_+$.
Left: for TB the drug is fast. Trajectories move almost vertically onto the
$d$-nullcline (dashed) and then slide along it, and the separatrix is almost
the vertical line $N = N_- \approx 1.3$. Right: for a slow drug the cells move
before the drug has settled. The separatrix tilts, and the threshold inoculum
runs from about $0.46$ (no drug yet) to about $3.2$ (full dose already
there), a factor of seven.

:::{prf:example} Same cells, same dose, opposite outcomes
:label: ex-pp-protocol

With the slow drug, let $10^8$ cells per mL ($N = 1$) meet a patient who
already holds the full dose. At $d = 1$ the separatrix sits at
$N \approx 3.2$, so the inoculum lies to its left and is cleared. If instead
the same inoculum is there first and treatment starts afterwards, then at
$d = 0$ the separatrix sits at $N \approx 0.46$, the inoculum lies to its
right, and the infection establishes itself. For TB the two thresholds are
$1.27$ and $1.34$, both close to $N_- = 1.30$, and the protocol barely
matters.
:::

::::{prf:remark} Computing the separatrix
:label: rmk-pp-computing-separatrix
:class: dropdown

The separatrix in the figure above is computed in two steps. Linearization
gives its direction at the saddle, and the flow, run backward in time,
extends it from there.

:::{prf:theorem} Stable manifold theorem
:label: thm-pp-stable-manifold

Let $\mathbf x^\ast$ be a saddle of $\dot{\mathbf x} = \mathbf f(\mathbf x)$,
and let the Jacobian at $\mathbf x^\ast$ have eigenvalues
$\lambda_s < 0 < \lambda_u$ with eigenvectors $\mathbf v_s$, $\mathbf v_u$.
The stable manifold of $\mathbf x^\ast$ is a smooth curve through
$\mathbf x^\ast$, tangent there to $\mathbf v_s$. The **unstable manifold**,
the set of points whose trajectories tend to $\mathbf x^\ast$ backward in time,
is a smooth curve through $\mathbf x^\ast$ tangent to $\mathbf v_u$.
:::

The proof is beyond these notes {cite:p}`strogatz2015`. Only $\mathbf v_s$
is tangent to the separatrix. The unstable eigenvector $\mathbf v_u$ is
tangent to the unstable manifold: the two trajectories that leave the saddle,
one toward the cleared state and one toward $N_+$.

**Step 1: the local picture.** Near the saddle, the separatrix is
approximately the line through $\mathbf x^\ast$ along $\mathbf v_s$,

$$
\mathbf x \approx \mathbf x^\ast \pm \varepsilon\,\mathbf v_s ,
$$ (eq-pp-sep-start)

with an error of order $\varepsilon^2$, since the curve bends away from its
tangent. This already shows the tilt. At $a = 2$ and for TB ($\delta = 17$),
the eigenvalues at $N_-$ are $0.357$ and $-39.6$, and
$\mathbf v_s \approx (0.066, 0.998)$ points almost along the drug axis: the
separatrix leaves the saddle almost vertically. For the slow drug
($\delta = 1/3$) they are $0.243$ and $-1.14$, and
$\mathbf v_s \approx (0.93, 0.36)$ has slope $0.39$.

**Step 2: follow the flow backward.** On the separatrix, trajectories run
into the saddle. To trace the separatrix out of the saddle, reverse time:
solve $\dot{\mathbf x} = -\mathbf f(\mathbf x)$ from the two starting points
{eq}`eq-pp-sep-start`, one for each sign, until the solution leaves the region
of interest. Each run draws one half of the separatrix. The figure uses
$\varepsilon = 10^{-7}$.

Backward integration is the right direction, not only the convenient one.
Near the saddle, a point at distance $\eta$ from the stable manifold moves
away from it like $\eta\,e^{\lambda_u t}$. This is why starts near the
threshold hesitate and why the smallest error decides their fate. Reversing
time turns the factor into $e^{-\lambda_u t}$: backward in time, the
separatrix *attracts* nearby solutions, so the $O(\varepsilon^2)$ error of the
start and the errors of the numerical method shrink instead of growing. The
price of a small $\varepsilon$ is time: the backward run needs a time of
about $\ln(1/\varepsilon)/|\lambda_s|$ to get away from the saddle, about two
weeks for the slow drug.

**A check: bisection.** For a fixed starting drug level $d_0$, the threshold
inoculum can also be found by simulating forward. Pick an interval of $N_0$
whose left end is cleared and whose right end settles at $N_+$, start at its
midpoint, and keep the half on which the outcome changes. Twenty halvings
narrow the interval by a factor of about a million. Each run only has to
decide which side it started on, so this needs no linearization, but it
gives one point of the curve, while one backward run gives the whole curve.
The two methods agree:

| | $\delta = 17$, $d_0 = 0$ | $\delta = 17$, $d_0 = 1$ | $\delta = 1/3$, $d_0 = 0$ | $\delta = 1/3$, $d_0 = 1$ |
|---|---|---|---|---|
| backward from the saddle | $1.2699$ | $1.3359$ | $0.4564$ | $3.1645$ |
| bisection | $1.2699$ | $1.3359$ | $0.4563$ | $3.1644$ |
::::

## Why the reduction works

The left panel above shows the two motions that the reduction assumes: a fast
vertical one onto the $d$-nullcline, then a slow one along it. The drug
equation carries the factor $\delta$, so for TB the arrows are nearly vertical
except near the $d$-nullcline. This section makes that precise with the
integrating factor of [](../overdamped/integrating-factor.md).

Let $z(t) = d(t) - \bar d(N(t))$ be the vertical distance from the state to
the $d$-nullcline. Since
$\dot d = -\delta(1 + N/N_c)\big(d - \bar d(N)\big)$, the chain rule gives

$$
\dot z = -\,\delta\Big(1 + \frac{N}{N_c}\Big)\,z \;-\; \bar d'(N)\,\dot N .
$$ (eq-pp-z)

The first term pulls $z$ to zero at a rate never below $\delta$. The second is
the speed of the moving target, and it is bounded.

:::{prf:proposition} The drug catches up
:label: prop-pp-catch-up

For a solution of {eq}`eq-pp-model` that starts in the box
$0 \le N \le K$, $0 \le d \le 1$,

$$
\big|d(t) - \bar d(N(t))\big|
\;\le\; \big|z(0)\big|\,e^{-\delta t} + \frac{M}{\delta},
\qquad M = \tfrac14\max(r, a) .
$$ (eq-pp-catch-up)
:::

:::{prf:proof}
:class: dropdown

By {prf:ref}`prop-pp-box` the solution stays in the box. There

$$
\big|\bar d'(N)\,\dot N\big|
= \frac{N/N_c}{(1 + N/N_c)^2}\,\Big|r\Big(1 - \frac NK\Big) - a\,d\Big|
\le \frac14\,\max(r, a) = M ,
$$

because $x/(1 + x)^2 \le 1/4$ for $x \ge 0$, and
$r(1 - N/K) - a\,d$ lies between $-a$ and $r$.

Write {eq}`eq-pp-z` as $\dot z + p(t)\,z = q(t)$ with
$p(t) = \delta(1 + N(t)/N_c) \ge \delta$ and $q(t) = -\bar d'(N)\dot N$,
$|q| \le M$; both are continuous because the solution is. By
{prf:ref}`prop-integrating-factor`, with $P(t) = \int_0^t p$,

$$
z(t) = z(0)\,e^{-P(t)} + \int_0^t e^{-(P(t) - P(s))}\,q(s)\,\dd s .
$$

Since $p \ge \delta$, $P(t) - P(s) \ge \delta\,(t - s)$, so

$$
|z(t)| \le |z(0)|\,e^{-\delta t} + M\int_0^t e^{-\delta(t - s)}\,\dd s
\le |z(0)|\,e^{-\delta t} + \frac{M}{\delta} .
$$
:::

The bound has a transient and a lag. For TB the distance decays on the time
scale $1/\delta \approx 1.4$ hours, and afterwards the drug trails its target
by at most $M/\delta \approx 0.03$ of the full dose ($a = 2$, $M = 0.5$/day).
The cells feel the error through their per-capita growth rate,

$$
\frac{\dot N}{N} = \frac{f(N)}{N} - a\,z ,
$$ (eq-pp-percapita)

with $f$ the reduced right-hand side {eq}`eq-bf-inoculum`. The transient
changes $\ln N$ by at most $a|z(0)|/\delta \approx 0.12$, a kick of about 12%
in the number of cells, and the lag changes the per-capita rate by at most
$aM/\delta \approx 0.06$/day. Both shrink like $1/\delta$, that is, with the
small parameter $\ve = r/\delta$ of {prf:ref}`rem-bf-qss-small-parameter`. A
12% kick matters only near the threshold, where it produces the tilt of the
separatrix. For the slow drug
$M/\delta = 1.5$, and the bound says nothing.

```{code-cell} python
:tags: [hide-input]
a, delta = 2.0, 17.0
M = max(r, a) / 4
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9))
ax = axs[0]
t = np.linspace(0, 1.0, 600)
for N0, d0, col in ((0.5, 0.0, C_OUTER), (3.0, 1.0, C_INNER), (8.0, 0.0, C_COMP)):
    sol = solve_ivp(rhs, (0, 1.0), [N0, d0], args=(a, delta), t_eval=t,
                    rtol=1e-10, atol=1e-13, method="LSODA")
    z = sol.y[1] - dbar(sol.y[0])
    ax.semilogy(24 * t, np.abs(z), color=col, lw=1.5,
                label=rf"$N_0 = {N0:g}$, $d_0 = {d0:g}$")
    ax.semilogy(24 * t, abs(d0 - dbar(N0)) * np.exp(-delta * t) + M / delta,
                color=col, lw=0.9, ls=":")
ax.axhline(M / delta, color=C_GREY, lw=0.8, ls="--")
ax.text(23.5, M / delta * 1.25, r"$M/\delta$", ha="right", fontsize=8, color="#666666")
ax.set_ylim(1e-5, 1.5); ax.set_xlim(0, 24)
ax.set_xlabel("time (hours)"); ax.set_ylabel(r"$|z| = |d - \bar d(N)|$")
ax.legend(frameon=False, fontsize=7, loc="lower left")
ax.set_title("the drug catches up within hours")

ax = axs[1]
f_red = lambda t, N: r * N * (1 - N / K) - a * N * dbar(N)
T = np.linspace(0, 15, 500)
for N0 in (0.8, 1.6, 4.0, 10.0):
    red = solve_ivp(f_red, (0, 15), [N0], t_eval=T, rtol=1e-10)
    ax.plot(T, red.y[0], color=C_OUTER, lw=2.6, alpha=0.35)
    for d0, ls in ((0.0, "-"), (1.0, "--")):
        sol = solve_ivp(rhs, (0, 15), [N0, d0], args=(a, delta), t_eval=T,
                        rtol=1e-10, atol=1e-13, method="LSODA")
        ax.plot(T, sol.y[0], color=C_EXACT, lw=0.9, ls=ls)
ax.plot([], [], color=C_OUTER, lw=2.6, alpha=0.35, label="reduced (phase line)")
ax.plot([], [], color=C_EXACT, lw=0.9, label="full, no drug at start")
ax.plot([], [], color=C_EXACT, lw=0.9, ls="--", label="full, full dose at start")
ax.set_xlabel("time (days)"); ax.set_ylabel(r"$N$ ($10^8$ cells/mL)")
ax.set_xlim(0, 15); ax.set_ylim(0, 10.5)
ax.legend(frameon=False, fontsize=7, loc="center right")
ax.set_title(r"TB: full and reduced models agree")
fig.tight_layout();
```

Left: the distance $|z|$ from the $d$-nullcline for three starts, with
$\delta = 17$/day and $a = 2$/day (solid), and the bound {eq}`eq-pp-catch-up`
(dotted). It falls by orders of magnitude within hours and stays below the lag
$M/\delta$ (dashed). Right: the cells in the reduced model of
[](../one-dimensional/saddle-node.md) (thick, pale) and in the full model with
either starting drug level (thin). They differ visibly only at $N_0 = 1.6$,
the start closest to the threshold.

### The eigenvalues say the same

:::{prf:proposition} Fast and slow eigenvalues
:label: prop-pp-fast-slow

Let $N^\ast > 0$ be a steady state of the reduced equation {eq}`eq-bf-inoculum`
and $f'(N^\ast)$ its slope. The eigenvalues of {eq}`eq-pp-jacobian` at
$(N^\ast, \bar d(N^\ast))$ satisfy, as $\delta \to \infty$,

$$
\lambda_\mathrm{fast} = -\,\delta\Big(1 + \frac{N^\ast}{N_c}\Big) + O(1),
\qquad
\lambda_\mathrm{slow} = f'(N^\ast) + O(1/\delta) .
$$
:::

:::{prf:proof}
:class: dropdown

By {eq}`eq-pp-det` and the formula
$f'(N^\ast) = N^\ast h'(N^\ast)/(1 + N^\ast/N_c)$ of {eq}`eq-bf-slope-h`,
with $c = 1 + N^\ast/N_c$,

$$
\lambda_\mathrm{fast}\,\lambda_\mathrm{slow} = \Delta = -\,\delta\,c\,f'(N^\ast),
\qquad
\lambda_\mathrm{fast} + \lambda_\mathrm{slow} = T = -\,\delta\,c - \frac{rN^\ast}{K} .
$$

The eigenvalues are real ({prf:ref}`prop-pp-no-spirals`). Since $T^2$ is of
order $\delta^2$ while $\Delta$ is of order $\delta$, expanding
$\sqrt{T^2 - 4\Delta} = |T|\big(1 - 2\Delta/T^2 + O(\delta^{-2})\big)$ in
{prf:ref}`prop-pp-trace-det` gives
$\lambda_\mathrm{fast} = T + O(1) = -\delta c + O(1)$ and
$\lambda_\mathrm{slow} = \Delta/T + O(1/\delta)$. Finally
$\Delta/T = f'(N^\ast)\big/\big(1 + rN^\ast/(K\delta c)\big) = f'(N^\ast) + O(1/\delta)$.
:::

For TB the slow eigenvalue is the slope of the phase line, and the slow
eigenvector is tangent to the $d$-nullcline. For the slow drug the two rates
are no longer far apart, and the agreement is gone.

Everything above rests on one fact:
$\partial\dot d/\partial d = -\delta(1 + N/N_c) \le -\delta$ for every $N$.
The drug is pulled toward its nullcline at a rate that is never small,
anywhere along it.
