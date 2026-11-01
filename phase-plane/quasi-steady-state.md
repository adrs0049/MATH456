---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/quasi-steady-state.pdf
    id: phase-plane-quasi-steady-state-pdf
downloads:
  - id: phase-plane-quasi-steady-state-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Quasi-Steady State

:::{tip} Big Idea
Many models contain a variable that is much faster than the others. If it
relaxes to a stable value, it can be replaced by that value, as a function of
the slow variables, and the model loses a dimension. This is the
**quasi-steady state** (QSS) approximation, and it is the singular limit of
[](../overdamped/singular-limit.md) again, now with two variables. The recipe
is to scale until a small $\ve$ multiplies the fast derivative, check that the
fast variable relaxes, reduce, and then ask what the reduction gets wrong:
usually only the first moments, and sometimes that is exactly what matters.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER, C_INNER, C_COMP = "#111111", "#0072B2", "#D55E00", "#009E73"
C_GREY = "#999999"
```

## Two variables, two clocks

In [](../overdamped/scaling-the-equation.md) the equation of motion of a cell
was scaled to $\ve\,\dd\tilde v/\dd\tilde t = \tilde F - \tilde v$
({eq}`eq-eom-scaled`), with $\ve = \tau_i/T$ the ratio of the cell's inertial
time to the time of the question. In [](../overdamped/singular-limit.md),
setting $\ve = 0$ replaced the differential equation by the algebraic one
$\tilde v = \tilde F$: the velocity followed the force instantaneously, and the
initial velocity was lost in a layer of width $\ve$. There the slow quantity,
the force, was prescribed. Now let it have an equation of its own.

:::{prf:definition} Fast-slow system
:label: def-qss-fast-slow

A **fast-slow system** is

$$
\frac{\dd x}{\dd\tau} = f(x, y),
\qquad
\ve\,\frac{\dd y}{\dd\tau} = g(x, y),
\qquad 0 < \ve \ll 1 ,
$$ (eq-qss-fast-slow)

with $x$ the **slow** variable and $y$ the **fast** variable. Time $\tau$ is
measured on the slow scale, and $\ve$ is the ratio of the fast time scale to
the slow one.
:::

$\ve$ is never given; it is found by scaling, as in
[](../overdamped/scaling-the-equation.md). Choose units for each variable and
for time so that $f$ and $g$ are of order one; $\ve$ is whatever dimensionless
group is left in front of a derivative. If no small group appears there, there
is no separation of time scales and the QSS approximation has no
justification. Setting $\ve = 0$ gives the same objects as in
[](../overdamped/singular-limit.md).

:::{prf:definition} Critical manifold, reduced problem, layer problem
:label: def-qss-reduced

The **critical manifold** is the set where the fast equation is at rest,
$g(x, y) = 0$. Where it can be solved as $y = h(x)$, the fast variable is
**slaved** to the slow one.

The **reduced problem** (the outer problem) is the slow equation on the
critical manifold,

$$
\frac{\dd x}{\dd\tau} = f\big(x, h(x)\big) .
$$ (eq-qss-reduced)

The **layer problem** (the inner problem) is the fast equation on the fast
time $T = \tau/\ve$, with $\ve = 0$:

$$
\frac{\dd y}{\dd T} = g(x_0, y), \qquad x = x_0 \text{ frozen} .
$$ (eq-qss-layer)
:::

The reduced problem is one dimension smaller than the full one and cannot take
an initial value for $y$; that value is lost in the layer, as the initial
velocity was. The layer problem says where $y$ goes while $x$ has not had time
to move, and for the reduction to make sense $y$ must go to the critical
manifold. The layer problem is a flow on a line with a steady state at
$h(x_0)$, so by the slope test ({prf:ref}`prop-ir-slope-test`) this holds when

$$
\frac{\partial g}{\partial y}\big(x, h(x)\big) < 0 .
$$ (eq-qss-attracting)

A critical manifold with this property is **attracting**. Where the slope is
positive the manifold repels, and the reduction describes nothing that can be
observed.

| | singular limit of a cell | fast-slow system |
|---|---|---|
| equation | $\ve\,\dd\tilde v/\dd\tilde t = \tilde F(\tilde t) - \tilde v$ | $\dd x/\dd\tau = f$, $\ve\,\dd y/\dd\tau = g(x, y)$ |
| slow quantity | the force $\tilde F$, prescribed | $x$, with its own equation |
| outer solution | $\tilde v = \tilde F$ ({prf:ref}`def-outer-inner`) | $y = h(x)$ and the reduced problem |
| inner solution | on $\sigma = \tilde t/\ve$, force frozen ({prf:ref}`def-outer-inner`) | on $T = \tau/\ve$, $x$ frozen |
| attracting | $\partial(\tilde F - \tilde v)/\partial\tilde v = -1 < 0$ | {eq}`eq-qss-attracting` |
| error | at most $\ve M$ ({prf:ref}`thm-composite-uniform`) | $O(\ve)$ (Tikhonov, below) |

:::{prf:theorem} Tikhonov
:label: thm-qss-tikhonov

Let the critical manifold $y = h(x)$ be attracting, in the sense of
{eq}`eq-qss-attracting`, along a solution $\bar x(\tau)$ of the reduced
problem {eq}`eq-qss-reduced` on a finite interval $0 \le \tau \le \tau_1$, with
$\bar x(0) = x(0)$. Let $y(0)$ lie in the basin of attraction of $h(x(0))$ for
the layer problem {eq}`eq-qss-layer`. Then, as $\ve \to 0$,

$$
x(\tau) \to \bar x(\tau) \quad\text{for } 0 \le \tau \le \tau_1,
\qquad
y(\tau) \to h\big(\bar x(\tau)\big) \quad\text{for } 0 < \tau_0 \le \tau \le \tau_1 ,
$$

uniformly, for every fixed $\tau_0 > 0$. For smooth $f$ and $g$ the error is
$O(\ve)$.
:::

We state this without proof {cite:p}`tikhonov1952`. It is the nonlinear
version of {prf:ref}`thm-composite-uniform`, as noted in
{prf:ref}`rmk-reading-bound`: the slow variable is right from the start, the
fast variable after the layer. For the drug model,
{prf:ref}`prop-pp-catch-up` of [](drug-in-the-plane.md) proves the part
that matters by hand.

:::{warning} Finite time
Tikhonov's theorem holds on a fixed finite interval $[0, \tau_1]$. It does
not say which steady state a solution reaches as $\tau \to \infty$. Near a
threshold a solution can hesitate for a long time before it decides, and that
long-time question needs a separate argument. The drug example below is such
a case.
:::

:::{important} The QSS recipe
1. **Scale.** Find units in which a small $\ve$ multiplies one derivative.
   Write $\ve$ in the original parameters and estimate it.
2. **Slave.** Solve $g(x, y) = 0$ for $y = h(x)$, the critical manifold.
3. **Check.** Verify $\partial g/\partial y < 0$ on it. Without this the
   reduction is wrong, not just inaccurate.
4. **Reduce.** Substitute $y = h(x)$ into the slow equation.
5. **Ask what was lost.** The initial value of $y$, and $O(\ve)$ corrections.
   Decide whether the question at hand can see them.
:::

## First example: Michaelis-Menten

:::{note} Biology background: enzymes
An **enzyme** is a protein that speeds up a chemical reaction without being
used up. The molecule it acts on is its **substrate**. The enzyme binds the
substrate, forming a **complex**, and the complex either falls apart again or
converts the substrate into **product** and releases the enzyme unchanged.
Cells contain far fewer enzyme molecules than substrate molecules, and each
enzyme is used again and again.
:::

With $S$ the substrate, $E$ the free enzyme, $C$ the complex and $P$ the
product, the reactions are

$$
S + E \;\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}\; C
\;\overset{k_2}{\longrightarrow}\; E + P .
$$

Mass action gives one equation per species. The enzyme is never used up, so
$E + C = E_0$, the total amount of enzyme, and $E = E_0 - C$. That leaves

$$
\frac{\dd S}{\dd t} = -k_1(E_0 - C)\,S + k_{-1}C,
\qquad
\frac{\dd C}{\dd t} = k_1(E_0 - C)\,S - (k_{-1} + k_2)\,C ,
$$ (eq-qss-mm-dim)

with $S(0) = S_0$ and $C(0) = 0$: the reaction starts when substrate meets
free enzyme.

**Scale.** Measure substrate against its initial amount, complex against the
total enzyme, and time in units of $1/(k_1E_0)$:
$u = S/S_0$, $v = C/E_0$, $\tau = k_1E_0\,t$. Then {eq}`eq-qss-mm-dim` becomes

$$
\frac{\dd u}{\dd\tau} = -u + (u + K - \lambda)\,v,
\qquad
\ve\,\frac{\dd v}{\dd\tau} = u - (u + K)\,v,
$$ (eq-qss-mm)

with

$$
\ve = \frac{E_0}{S_0}, \qquad
K = \frac{k_{-1} + k_2}{k_1S_0}, \qquad
\lambda = \frac{k_2}{k_1S_0} .
$$

The small parameter is the ratio of enzyme to substrate. It is small because
the enzyme is recycled: there is much less of it than of the substrate it
processes.

:::{prf:proposition} The Michaelis-Menten law
:label: prop-qss-mm

The critical manifold of {eq}`eq-qss-mm` is $v = h(u) = u/(u + K)$, and it is
attracting for every $u \ge 0$. The reduced problem is

$$
\frac{\dd u}{\dd\tau} = -\frac{\lambda\,u}{u + K},
$$

which in the original units says that product forms at the rate

$$
\frac{\dd P}{\dd t} = k_2C = \frac{V_\mathrm{max}\,S}{K_M + S},
\qquad V_\mathrm{max} = k_2E_0,\quad K_M = \frac{k_{-1} + k_2}{k_1} .
$$ (eq-qss-mm-law)
:::

:::{prf:proof}
:class: dropdown

Setting $\ve = 0$ gives $u = (u + K)v$, so $v = u/(u + K)$. With
$g(u, v) = u - (u + K)v$, $\partial g/\partial v = -(u + K) \le -K < 0$, so the
manifold attracts at rate at least $K$ on the fast time. Substituting into the
slow equation,

$$
\frac{\dd u}{\dd\tau} = -u + (u + K - \lambda)\,\frac{u}{u + K}
= -\frac{\lambda\,u}{u + K} .
$$

In the original units, $C = E_0 v = E_0 S/(S + K_M)$ since
$u/(u + K) = S/(S + K S_0) = S/(S + K_M)$, and $\dd P/\dd t = k_2 C$.
:::

This is the **Michaelis-Menten law** {cite:p}`michaelis1913`, derived by this
argument by Briggs and Haldane {cite:p}`briggs1925`. The rate is proportional
to $S$ when substrate is scarce and saturates at $V_\mathrm{max}$ when every
enzyme is busy.

**What was lost.** The reduced problem forces $C = E_0S_0/(K_M + S_0)$ at
$t = 0$, but the complex starts at $C = 0$. On the fast time
$T = \tau/\ve = k_1S_0\,t$, with $u = 1$ frozen, the layer problem is
$\dd v/\dd T = 1 - (1 + K)\,v$, so

$$
v = \frac{1}{1 + K}\Big(1 - e^{-(1 + K)T}\Big),
\qquad\text{that is}\qquad
C = \frac{E_0S_0}{K_M + S_0}\Big(1 - e^{-k_1(S_0 + K_M)\,t}\Big).
$$

The complex builds up in a time $1/\big(k_1(S_0 + K_M)\big)$, during which
almost no substrate is used. Unless the experiment measures that first
instant, the reduction loses nothing.

```{code-cell} python
:tags: [hide-input]
K, lam = 1.0, 0.5


def mm(t, y, eps):
    u, v = y
    return [-u + (u + K - lam) * v, (u - (u + K) * v) / eps]


fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
ax = axs[0]
U = np.linspace(0, 1.05, 200)
ax.plot(U, U / (U + K), color=C_OUTER, lw=1.6, ls="--",
        label=r"critical manifold $v = u/(u+K)$")
eps = 0.05
for u0 in (1.0, 0.6, 0.3):
    sol = solve_ivp(mm, (0, 20), [u0, 0], args=(eps,), rtol=1e-9, atol=1e-12,
                    method="LSODA", dense_output=True)
    T = np.concatenate([np.linspace(0, 5 * eps, 200), np.linspace(5 * eps, 20, 400)])
    u, v = sol.sol(T)
    ax.plot(u, v, color=C_EXACT, lw=1.2)
    ax.annotate("", xy=(u[60], v[60]), xytext=(u[20], v[20]),
                arrowprops=dict(arrowstyle="-|>", color=C_EXACT, lw=1))
ax.set_xlabel("substrate $u = S/S_0$"); ax.set_ylabel("complex $v = C/E_0$")
ax.set_xlim(0, 1.05); ax.set_ylim(0, 0.6)
ax.legend(loc="upper left", frameon=False, fontsize=8)
ax.set_title(r"phase plane, $\varepsilon = 0.05$")

ax = axs[1]
T = np.linspace(0, 6, 600)
for eps, col in ((0.3, "#9ecae1"), (0.05, C_EXACT)):
    sol = solve_ivp(mm, (0, 6), [1, 0], args=(eps,), t_eval=T, rtol=1e-9,
                    atol=1e-12, method="LSODA")
    ax.plot(T, sol.y[0], color=col, lw=1.6, label=rf"full, $\varepsilon = {eps}$")
red = solve_ivp(lambda t, u: -lam * u / (u + K), (0, 6), [1], t_eval=T, rtol=1e-9)
ax.plot(T, red.y[0], color=C_OUTER, ls="--", lw=1.6, label="reduced")
ax.set_xlabel(r"time $\tau$"); ax.set_ylabel("substrate $u$")
ax.legend(frameon=False, fontsize=8)
ax.set_title("the reduction tracks the substrate")
fig.tight_layout();
```

The Michaelis-Menten system {eq}`eq-qss-mm` with $K = 1$, $\lambda = 0.5$.
Left: three solutions starting with no complex. Each jumps almost vertically
onto the critical manifold (the layer, in which $u$ hardly moves) and then
slides along it (the reduced problem). Right: the substrate. For
$\ve = 0.05$ the full and reduced solutions are hard to tell apart; for
$\ve = 0.3$ the layer uses up a visible amount of substrate.

:::{prf:remark} Which ratio has to be small
:label: rmk-qss-segel
:class: dropdown

The scaling above gives $\ve = E_0/S_0$, and $E_0 \ll S_0$ is the textbook
condition. Segel and Slemrod {cite:p}`segel1989` compared the time the complex
needs to build up with the time the substrate needs to be used up, and found
the sharper condition

$$
\frac{E_0}{S_0 + K_M} \ll 1 .
$$

When $K_M$ is large the reduction works even if there is as much enzyme as
substrate. Different scalings can give different $\ve$, and the right one
compares the two time scales the question actually involves.
:::

## Second example: the drug is fast

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
recipe instead, and compute what the reduction gets wrong.

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

### What does not depend on ε

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

### The first correction slows time

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

### The layer moves the threshold

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

In both examples, Michaelis-Menten and the drug, the critical manifold is
attracting everywhere, so the recipe works along every solution. When the
critical manifold folds back on itself, condition {eq}`eq-qss-attracting`
fails at the folds, and the fast variable jumps. That is the situation of
[](fitzhugh-nagumo.md).
