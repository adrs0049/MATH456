---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/michaelis-menten.pdf
    id: phase-plane-michaelis-menten-pdf
downloads:
  - id: phase-plane-michaelis-menten-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Michaelis-Menten Kinetics

:::{tip} Big Idea
An enzyme is recycled, so a cell needs far less of it than of the substrate it
converts. That ratio is the small parameter: the complex of enzyme and
substrate is fast, reaches its quasi-steady state almost at once, and the
substrate is then used up at the saturating Michaelis-Menten rate. The only
thing the reduction loses is the brief moment in which the complex forms.
:::

```{code-cell} python
:tags: [hide-input, remove-output]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_OUTER = "#111111", "#0072B2"
```

This page applies the recipe of [](quasi-steady-state.md) to the simplest
enzyme reaction.

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

