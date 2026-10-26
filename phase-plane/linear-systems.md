---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/linear-systems.pdf
    id: phase-plane-linear-systems-pdf
downloads:
  - id: phase-plane-linear-systems-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Linear Systems in the Plane

:::{tip} Big Idea
Near a steady state, a model with two variables behaves like a linear system
$\dot{\mathbf u} = J\mathbf u$. Look for solutions that keep their direction
and only change their size: the direction must be an eigenvector, and along it
the motion is a phase line, $e^{\lambda t}$. Adding such solutions gives every
solution. Two numbers, the trace and the determinant of $J$, decide what the
picture looks like: a saddle, a node, or a spiral.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.linalg import expm
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER, C_INNER, C_COMP = "#111111", "#0072B2", "#D55E00", "#009E73"
C_GREY = "#999999"

r, K, Nc = 1.0, 10.0, 1.0           # 1/day, 1e8 cells/mL, 1e8 cells/mL
dbar = lambda N: 1 / (1 + N / Nc)   # drug nullcline (the QSS drug level)
NN = np.linspace(0, K, 400)


def interior_states(a):
    z = np.roots([-r / (K * Nc), r * (1 / Nc - 1 / K), r - a])
    return sorted(x.real for x in z if abs(x.imag) < 1e-12 and x.real > 0)


def jacobian(N, d, a, delta):
    return np.array([[r * (1 - 2 * N / K) - a * d, -a * N],
                     [-delta * d / Nc, -delta * (1 + N / Nc)]])


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

On the phase line, near a steady state $N^\ast$ the deviation
$x = N - N^\ast$ obeys $\dot x \approx f'(N^\ast)\,x$
({prf:ref}`prop-ir-slope-test`). In the plane, a model
$\dot{\mathbf x} = \mathbf f(\mathbf x)$ with two variables is linearized at
a steady state $\mathbf x^\ast$ in the same way: the deviation
$\mathbf u = \mathbf x - \mathbf x^\ast$ obeys, to first order,

$$
\dot{\mathbf u} = J\,\mathbf u ,
\qquad
J = \begin{pmatrix} J_{11} & J_{12} \\ J_{21} & J_{22} \end{pmatrix},
$$ (eq-pp-linear-system)

where $J$ is the **Jacobian matrix** of $\mathbf f$ at $\mathbf x^\ast$, the
matrix of partial derivatives $J_{ij} = \partial f_i/\partial x_j$. This page
solves {eq}`eq-pp-linear-system` and sorts its solutions into types;
[](drug-in-the-plane.md) applies the result to a model.

## Separating time from direction

On the phase line, $\dot x = \lambda x$ is solved by $x(0)\,e^{\lambda t}$. The
two equations of {eq}`eq-pp-linear-system` are coupled, so they cannot be
solved one at a time. Instead, separate time from direction, as in separation
of variables: look for a solution that keeps a fixed direction $\mathbf v$ and
changes only its size,

$$
\mathbf u(t) = s(t)\,\mathbf v .
$$ (eq-pp-separation)

Substituting into $\dot{\mathbf u} = J\mathbf u$ gives

$$
\dot s(t)\,\mathbf v = s(t)\,J\mathbf v .
$$

The left side points along $\mathbf v$, so the right side must too. If
$\mathbf v$ is an **eigenvector**, $J\mathbf v = \lambda\mathbf v$, the vector
equation becomes the scalar equation $\dot s = \lambda s$: along $\mathbf v$
the motion is a phase line with slope $\lambda$, and
$s(t) = s(0)\,e^{\lambda t}$.

:::{prf:proposition} Superposition principle
:label: prop-pp-superposition

If $\mathbf u_1(t)$ and $\mathbf u_2(t)$ solve $\dot{\mathbf u} = J\mathbf u$,
then so does $c_1\mathbf u_1(t) + c_2\mathbf u_2(t)$ for all constants $c_1$,
$c_2$, real or complex.
:::

:::{prf:proof}
:class: dropdown

Matrix multiplication is linear, so
$\frac{d}{dt}(c_1\mathbf u_1 + c_2\mathbf u_2)
= c_1J\mathbf u_1 + c_2J\mathbf u_2 = J(c_1\mathbf u_1 + c_2\mathbf u_2)$.
:::

:::{prf:proposition} Eigenvector solutions
:label: prop-pp-eigen-solutions

If $J\mathbf v = \lambda\mathbf v$ with $\mathbf v \ne 0$, then
$\mathbf u(t) = e^{\lambda t}\mathbf v$ solves $\dot{\mathbf u} = J\mathbf u$.
If $J$ has eigenvalues $\lambda_1 \ne \lambda_2$ with eigenvectors
$\mathbf v_1$, $\mathbf v_2$, every solution is

$$
\mathbf u(t) = c_1\,e^{\lambda_1 t}\,\mathbf v_1 + c_2\,e^{\lambda_2 t}\,\mathbf v_2,
$$ (eq-pp-general)

with $c_1$, $c_2$ fixed by $\mathbf u(0) = c_1\mathbf v_1 + c_2\mathbf v_2$.
:::

:::{prf:proof}
:class: dropdown

$\frac{d}{dt}\big(e^{\lambda t}\mathbf v\big) = \lambda e^{\lambda t}\mathbf v
= e^{\lambda t}J\mathbf v = J\big(e^{\lambda t}\mathbf v\big)$. By
{prf:ref}`prop-pp-superposition`, {eq}`eq-pp-general` then solves the system
for all $c_1$, $c_2$. The eigenvectors are not parallel: if
$\mathbf v_2 = k\mathbf v_1$, then $J\mathbf v_2 = \lambda_1\mathbf v_2$ and
$\lambda_1 = \lambda_2$. So every initial value can be written as
$c_1\mathbf v_1 + c_2\mathbf v_2$, and since a solution is determined by its
initial value, every solution is of the form {eq}`eq-pp-general`.
:::

In coordinates along the eigenvectors the system splits into two phase lines,
$\dot c_1 = \lambda_1 c_1$ and $\dot c_2 = \lambda_2 c_2$: these are the
variables that separate. The eigenvalues are the roots of a quadratic
({prf:ref}`prop-pp-trace-det` below), so two cases are left: a complex pair,
and a repeated root.

## Complex eigenvalues

Since $J$ is real, complex eigenvalues come in a conjugate pair
$\lambda = \alpha \pm i\omega$, and so do their eigenvectors,
$\mathbf v = \mathbf a \pm i\mathbf b$ with real vectors $\mathbf a$,
$\mathbf b$. The computation behind {prf:ref}`prop-pp-eigen-solutions` works
just as well with complex numbers, but $e^{\lambda t}\mathbf v$ is then a
complex solution. Its real and imaginary parts are real solutions.

:::{prf:proposition} Complex eigenvalues
:label: prop-pp-complex

Let $J\mathbf v = \lambda\mathbf v$ with $\lambda = \alpha + i\omega$,
$\omega \ne 0$, and $\mathbf v = \mathbf a + i\mathbf b$. Then

$$
\mathbf u_1(t) = e^{\alpha t}\big(\mathbf a\cos\omega t - \mathbf b\sin\omega t\big),
\qquad
\mathbf u_2(t) = e^{\alpha t}\big(\mathbf a\sin\omega t + \mathbf b\cos\omega t\big)
$$ (eq-pp-complex)

are real solutions, and every real solution is $c_1\mathbf u_1 + c_2\mathbf u_2$
with $c_1\mathbf a + c_2\mathbf b = \mathbf u(0)$.
:::

:::{prf:proof}
:class: dropdown

By Euler's formula, $e^{\lambda t} = e^{\alpha t}(\cos\omega t + i\sin\omega t)$.
Multiplying by $\mathbf a + i\mathbf b$ and collecting terms,
$e^{\lambda t}\mathbf v = \mathbf u_1 + i\,\mathbf u_2$. This solves the
system, so $\dot{\mathbf u}_1 + i\dot{\mathbf u}_2 = J\mathbf u_1 + iJ\mathbf u_2$.
Since $J$ is real, $J\mathbf u_1$ and $J\mathbf u_2$ are real, and comparing
real and imaginary parts gives $\dot{\mathbf u}_1 = J\mathbf u_1$ and
$\dot{\mathbf u}_2 = J\mathbf u_2$.

At $t = 0$, $\mathbf u_1 = \mathbf a$ and $\mathbf u_2 = \mathbf b$. These are
not parallel: otherwise $\mathbf v$ would be a complex multiple of a real
vector $\mathbf w \ne 0$, and $J\mathbf w = \lambda\mathbf w$ is impossible
because $J\mathbf w$ is real and $\lambda\mathbf w$ is not. So every initial
value is $c_1\mathbf a + c_2\mathbf b$, and
{prf:ref}`prop-pp-superposition` and uniqueness finish the proof as before.
:::

The brackets in {eq}`eq-pp-complex` are periodic with period $2\pi/\omega$:
on their own they go around an ellipse in the plane. The factor $e^{\alpha t}$
shrinks the ellipse if $\alpha < 0$ and grows it if $\alpha > 0$, so
trajectories are **spirals**. If $\alpha = 0$ the ellipse neither shrinks
nor grows, and the orbits are closed: a **center**. The imaginary part
$\omega$ sets the rotation, and the real part $\alpha$ alone decides
stability.

:::{prf:example} A pure spiral
:label: ex-pp-spiral

Take

$$
J = \begin{pmatrix} \alpha & -\omega \\ \omega & \alpha \end{pmatrix},
\qquad \omega > 0 .
$$

Then $J\mathbf v = (\alpha + i\omega)\mathbf v$ for $\mathbf v = (1, -i)$, so
$\mathbf a = (1, 0)$, $\mathbf b = (0, -1)$, and {eq}`eq-pp-complex` gives

$$
\mathbf u_1(t) = e^{\alpha t}\begin{pmatrix}\cos\omega t\\ \sin\omega t\end{pmatrix},
\qquad
\mathbf u_2(t) = e^{\alpha t}\begin{pmatrix}\sin\omega t\\ -\cos\omega t\end{pmatrix}.
$$

The solution $\mathbf u_1$ is at distance $e^{\alpha t}$ from the origin and
at angle $\omega t$: it turns counterclockwise at the constant rate $\omega$
while its distance shrinks or grows at the rate $\alpha$.
:::

## Repeated eigenvalues

The two roots coincide when $T^2 = 4\Delta$, in the notation of
{prf:ref}`prop-pp-trace-det`. If $J = \lambda I$, every vector is an
eigenvector and every solution $e^{\lambda t}\mathbf u(0)$ moves straight in
or out along its own line, a **star**. Otherwise all
eigenvectors are multiples of one vector $\mathbf v$, and
{prf:ref}`prop-pp-eigen-solutions` gives only one solution. A triangular
example, solved with the integrating factor, shows what the second looks
like.

:::{prf:example} A repeated eigenvalue
:label: ex-pp-repeated

Take

$$
\dot x = \lambda x + y,
\qquad
\dot y = \lambda y,
\qquad
J = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}.
$$

The only eigenvectors are multiples of $\mathbf v = (1, 0)$. The second
equation gives $y = y_0\,e^{\lambda t}$, and the first becomes
$\dot x = \lambda x + y_0\,e^{\lambda t}$: $x$ is forced at its own rate.
With the integrating factor $e^{-\lambda t}$,
$\frac{d}{dt}\big(e^{-\lambda t}x\big) = y_0$, so
$x = (x_0 + y_0\,t)\,e^{\lambda t}$ and

$$
\mathbf u(t) = x_0\,e^{\lambda t}\begin{pmatrix}1\\0\end{pmatrix}
+ y_0\,e^{\lambda t}\begin{pmatrix}t\\1\end{pmatrix}.
$$

The second solution is $e^{\lambda t}(t\,\mathbf v + \mathbf w)$ with
$\mathbf w = (0, 1)$.
:::

:::{prf:proposition} Repeated eigenvalue
:label: prop-pp-repeated

Let $\lambda$ be a double root, $J \ne \lambda I$, and $J\mathbf v = \lambda\mathbf v$
with $\mathbf v \ne 0$. There is a vector $\mathbf w$ with
$(J - \lambda I)\,\mathbf w = \mathbf v$, and every solution is

$$
\mathbf u(t) = c_1\,e^{\lambda t}\,\mathbf v
+ c_2\,e^{\lambda t}\big(t\,\mathbf v + \mathbf w\big),
$$ (eq-pp-repeated)

with $c_1\mathbf v + c_2\mathbf w = \mathbf u(0)$.
:::

:::{prf:proof}
:class: dropdown

Following {prf:ref}`ex-pp-repeated`, try
$\mathbf u = e^{\lambda t}(t\,\mathbf v + \mathbf w)$. Then

$$
\dot{\mathbf u} = e^{\lambda t}\big(\lambda t\,\mathbf v + \lambda\mathbf w + \mathbf v\big),
\qquad
J\mathbf u = e^{\lambda t}\big(\lambda t\,\mathbf v + J\mathbf w\big),
$$

which agree exactly when $J\mathbf w = \lambda\mathbf w + \mathbf v$.

Such a $\mathbf w$ exists. Let $M = J - \lambda I$. A double root means
$T = 2\lambda$ and $\Delta = \lambda^2$, so $M$ has trace $T - 2\lambda = 0$
and determinant $\Delta - \lambda T + \lambda^2 = 0$. For every $2\times 2$
matrix, $M^2 = (\operatorname{tr} M)\,M - (\det M)\,I$, as multiplying out
shows, so here $M^2 = 0$. Since $M \ne 0$, some $\mathbf z$ has
$M\mathbf z \ne 0$, and $M(M\mathbf z) = 0$ says $M\mathbf z$ is an
eigenvector, $M\mathbf z = k\mathbf v$ with $k \ne 0$. Take
$\mathbf w = \mathbf z/k$.

Finally $\mathbf v$ and $\mathbf w$ are not parallel, since $M\mathbf v = 0$
and $M\mathbf w = \mathbf v \ne 0$. So every initial value is
$c_1\mathbf v + c_2\mathbf w$, and {prf:ref}`prop-pp-superposition` and
uniqueness finish the proof.
:::

For large $|t|$ the term $t\,\mathbf v$ outgrows $\mathbf w$ in
{eq}`eq-pp-repeated`, so every trajectory ends up parallel to the one
eigenvector line, a **degenerate node**. The factor $t$ only delays the
decay: if $\lambda < 0$, then $t\,e^{\lambda t} \to 0$, and the steady state
is still stable. A repeated root is the borderline between a node
($T^2 > 4\Delta$) and a spiral ($T^2 < 4\Delta$): trajectories start to turn,
but not all the way round.

## Trace and determinant

:::{prf:proposition} Trace and determinant
:label: prop-pp-trace-det

Let $T = \operatorname{tr} J = J_{11} + J_{22}$ and
$\Delta = \det J = J_{11}J_{22} - J_{12}J_{21}$. The eigenvalues of $J$ are the roots of
$\lambda^2 - T\lambda + \Delta = 0$,

$$
\lambda_{1,2} = \frac{T \pm \sqrt{T^2 - 4\Delta}}{2},
\qquad
\lambda_1 + \lambda_2 = T, \quad \lambda_1\lambda_2 = \Delta .
$$

Consequently:

| condition | eigenvalues | the steady state is |
|---|---|---|
| $\Delta < 0$ | real, opposite signs | a **saddle**: unstable |
| $\Delta > 0$, $T < 0$, $T^2 > 4\Delta$ | real, both negative | a **stable node** |
| $\Delta > 0$, $T < 0$, $T^2 < 4\Delta$ | complex, negative real part | a **stable spiral** |
| $\Delta > 0$, $T > 0$, $T^2 > 4\Delta$ | real, both positive | an **unstable node** |
| $\Delta > 0$, $T > 0$, $T^2 < 4\Delta$ | complex, positive real part | an **unstable spiral** |

The borderlines are: $T^2 = 4\Delta$, a repeated eigenvalue (star or
degenerate node); $\Delta > 0$, $T = 0$, purely imaginary eigenvalues (a
center); and $\Delta = 0$, a zero eigenvalue.
:::

:::{prf:proof}
:class: dropdown

$J - \lambda I$ has determinant
$(J_{11} - \lambda)(J_{22} - \lambda) - J_{12}J_{21} = \lambda^2 - T\lambda + \Delta$,
and $\lambda$ is an eigenvalue exactly when this vanishes. The quadratic
formula gives the roots, and comparing coefficients gives their sum and
product. If $\Delta < 0$ the product is negative, so the roots are real with
opposite signs. If $\Delta > 0$ and the roots are real, they have the same
sign, which is the sign of their sum $T$; if they are complex they are
$T/2 \pm i\sqrt{4\Delta - T^2}/2$, with real part $T/2$.
:::

A **saddle** has one direction along which solutions approach and one along
which they leave: it is an unstable point of a phase line with one extra,
stable direction. The table is in Strogatz §5.2 and §6.3
{cite:p}`strogatz2015`.

Here is every case, each drawn from the exact solution of the linear system.

```{code-cell} python
:tags: [hide-input]
from scipy.linalg import expm

LIN = [("saddle", [[1, 1], [4, -2]]),
       ("stable node", [[-1, 1], [0, -3]]),
       ("unstable node", [[1, -1], [0, 3]]),
       ("degenerate node", [[-1, 1], [0, -1]]),
       ("stable spiral", [[-0.5, -3], [1, -0.5]]),
       ("unstable spiral", [[0.3, -3], [1, 0.3]]),
       ("center", [[0, -3], [1, 0]]),
       ("star", [[-1, 0], [0, -1]])]

fig, axs = plt.subplots(2, 4, figsize=(7.4, 4.0))
tt = np.linspace(-8, 8, 1601)
for k, (ax, (name, A)) in enumerate(zip(axs.flat, LIN)):
    A = np.array(A, float)
    phase_plane(ax, lambda x, y: A[0, 0] * x + A[0, 1] * y,
                lambda x, y: A[1, 0] * x + A[1, 1] * y, (-1, 1), (-1, 1),
                null=False, fixed_points=False, density=0.7, n_grid=81)
    lam, V = np.linalg.eig(A)
    if np.all(np.isreal(lam)) and name != "star":        # eigenvector lines
        for j in range(2):
            v = V[:, j].real
            ax.plot([-2 * v[0], 2 * v[0]], [-2 * v[1], 2 * v[1]], color=C_OUTER,
                    lw=1.3, ls="--", zorder=3)
    if name == "stable spiral":
        starts = [(0.9, 0), (-0.9, 0)]
    elif name == "unstable spiral":
        starts = [(0.15, 0), (-0.15, 0)]
    elif name == "center":
        starts = [(0.3, 0), (0.6, 0), (0.9, 0)]
    else:
        starts = [0.75 * np.array([np.cos(th), np.sin(th)])
                  for th in np.linspace(0, 2 * np.pi, 6, endpoint=False) + 0.3]
    for z0 in starts:                                     # x(t) = e^{At} x(0)
        Z = np.array([expm(A * s) @ np.array(z0, float) for s in tt]).T
        Z[:, np.abs(Z).max(axis=0) > 3] = np.nan
        ax.plot(*Z, color=C_EXACT, lw=0.8, zorder=4)
        i = np.nanargmin(np.abs(np.hypot(*Z) - 0.6))
        ax.annotate("", xy=Z[:, i + 15], xytext=Z[:, i],
                    arrowprops=dict(arrowstyle="-|>", color=C_EXACT, lw=0.8))
    if name == "center":
        ax.plot([0], [0], "o", color=C_GREY, ms=5, zorder=6)
    else:
        mark(ax, 0, 0, np.all(lam.real < 0), ms=5)
    ax.set_title(f"({'abcdefgh'[k]}) {name}", fontsize=8.5)
    ax.text(0.03, 0.03, rf"$T = {np.trace(A):g}$, $\Delta = {np.linalg.det(A):.3g}$",
            transform=ax.transAxes, fontsize=7,
            bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
fig.tight_layout();
```

Phase portraits of $\dot{\mathbf x} = J\mathbf x$ on the square
$-1 \le x, y \le 1$, with trace and determinant in each corner. Black: exact
solutions $e^{Jt}\mathbf x(0)$, forward and backward in time. Dashed blue:
eigenvector lines, which are themselves trajectories. The matrices are

$$
\text{(a)}\begin{pmatrix} 1 & 1 \\ 4 & -2 \end{pmatrix}\!,\;
\text{(b)}\begin{pmatrix} -1 & 1 \\ 0 & -3 \end{pmatrix}\!,\;
\text{(c)}\begin{pmatrix} 1 & -1 \\ 0 & 3 \end{pmatrix}\!,\;
\text{(d)}\begin{pmatrix} -1 & 1 \\ 0 & -1 \end{pmatrix}\!,\;
\text{(e)}\begin{pmatrix} -0.5 & -3 \\ 1 & -0.5 \end{pmatrix}\!,\;
\text{(f)}\begin{pmatrix} 0.3 & -3 \\ 1 & 0.3 \end{pmatrix}\!,\;
\text{(g)}\begin{pmatrix} 0 & -3 \\ 1 & 0 \end{pmatrix}\!,\;
\text{(h)}\begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\!.
$$

Things to notice:

- **(a)** The eigenvalues are $2$ and $-3$, with eigenvectors $(1, 1)$ and
  $(1, -4)$. They need not be perpendicular. Trajectories come in along the
  steep line and leave along the line $y = x$.
- **(b)** The eigenvalues are $-1$, eigenvector $(1, 0)$, and $-3$,
  eigenvector $(1, -2)$. Trajectories first move roughly parallel to the
  fast direction $(1, -2)$, then turn and arrive tangent to the slow
  direction, the $x$ axis.
- **(c)** is (b) with time reversed: $J$ is replaced by $-J$, the same
  curves with the arrows flipped.
- **(e), (f), (g)** The three matrices differ only on the diagonal, that is,
  in the trace. All three turn counterclockwise: at $(1, 0)$ the velocity is the first
  column $(J_{11}, J_{21})$ of $J$, which points up because $J_{21} = 1 > 0$. In (g),
  $\dot x = -3y$ and $\dot y = x$ give $\frac{\dd}{\dd t}(x^2 + 3y^2) = 0$, so
  the orbits are ellipses $\sqrt 3$ times as wide as they are tall. The sign of $T$ decides
  whether the turning goes inward, outward, or neither.
- **(d)** and **(h)** have the same trace and determinant, $\lambda = -1$
  twice, but different pictures: one eigenvector line in (d), every line in
  (h), as in {prf:ref}`prop-pp-repeated`. Trace and determinant fix the
  eigenvalues, not the eigenvectors.

**The trace–determinant plane.** {prf:ref}`prop-pp-trace-det` sorts every
matrix by the point $(T, \Delta)$. Plot that point, and the type of the
steady state can be read off from the region it lands in.

```{code-cell} python
:tags: [hide-input]
fig, ax = plt.subplots(figsize=(5.6, 3.4))
TT = np.linspace(-4.5, 4.5, 400)
ax.fill_between(TT, -2, 0, color=C_INNER, alpha=0.07, lw=0)
ax.plot(TT, TT ** 2 / 4, color=C_OUTER, lw=1.5)
ax.axhline(0, color=C_EXACT, lw=0.8); ax.axvline(0, color=C_EXACT, lw=0.8)
ax.plot([0, 0], [0, 4.5], color=C_COMP, lw=2.5)
for x, y, s in ((-2.3, 4.3, "stable spirals"), (0.6, 4.3, "unstable spirals"),
                (-4.4, 0.5, "stable nodes"), (3.0, 0.5, "unstable nodes"),
                (-0.5, -1.3, "saddles")):
    ax.text(x, y, s, fontsize=8, color="#555555", ha="left", va="center")
ax.text(0.1, 1.9, "centers", fontsize=8, color=C_COMP)
ax.text(3.15, 2.0, r"$T^2 = 4\Delta$", fontsize=8, color=C_OUTER)
pts = {}
for k, (name, A) in enumerate(LIN):
    A = np.array(A, float)
    pts.setdefault((round(np.trace(A), 3), round(np.linalg.det(A), 3)), []).append("abcdefgh"[k])
for (Tr, Dt), labels in pts.items():
    ax.plot([Tr], [Dt], "o", color=C_EXACT, ms=4, zorder=5)
    ax.annotate(", ".join(labels), (Tr, Dt), xytext=(4, 4), textcoords="offset points",
                fontsize=8)
states = [((0.0, 1.0), "cleared", (-40, 2))] + [
    ((n, dbar(n)), lab, (6, -10)) for n, lab in zip(interior_states(2.0), (r"$N_-$", r"$N_+$"))]
for (N, d), lab, off in states:
    Jd = jacobian(N, d, 2.0, 1 / 3)
    ax.plot([np.trace(Jd)], [np.linalg.det(Jd)], "*", color=C_INNER, ms=9, zorder=6)
    ax.annotate(lab, (np.trace(Jd), np.linalg.det(Jd)), xytext=off,
                textcoords="offset points", fontsize=8, color=C_INNER)
ax.set_xlim(-4.5, 4.5); ax.set_ylim(-2, 5)
ax.set_xlabel(r"trace $T$"); ax.set_ylabel(r"determinant $\Delta$")
fig.tight_layout();
```

The trace–determinant plane. Below the axis $\Delta = 0$ (shaded): saddles.
Above it, the parabola $T^2 = 4\Delta$ separates nodes (outside) from
spirals (inside), and the axis $T = 0$ separates stable (left) from unstable
(right). Black dots: the gallery, except the saddle (a) at
$(T, \Delta) = (-1, -6)$, below the frame; (d) and (h) sit at the same point on the
parabola, and the center (g) on the green half-axis. Orange stars: the three
steady states of the drug model {eq}`eq-pp-model` of [](drug-in-the-plane.md),
with $a = 2$/day and the
slow drug, $\delta = 1/3$/day. None of them is inside the
parabola.

:::{warning} Borderline cases do not survive linearization
For a nonlinear system, the Jacobian gives the right picture near a steady
state whenever $(T, \Delta)$ lies strictly inside one of the open regions:
saddles, nodes and spirals keep their type, because a small nonlinear term
cannot push the point across a boundary. On a boundary it can. A linear
center may be a slow stable spiral, a slow unstable spiral, or a true center
of the nonlinear system; a degenerate node may become a node or a spiral;
and $\Delta = 0$ says nothing at all. That is where the saddle-node
bifurcation happens. Strogatz §6.3 {cite:p}`strogatz2015`.
:::
