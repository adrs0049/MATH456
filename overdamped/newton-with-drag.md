---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/newton-with-drag.pdf
    id: overdamped-newton-with-drag-pdf
downloads:
  - id: overdamped-newton-with-drag-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Newton's Second Law with Drag

:::{tip} Big Idea
Write down Newton's second law for a cell in a viscous medium, with every term
present. Its units reveal two time scales before anything is solved: the
inertial stopping time $\tau_i = m/\gamma$, which belongs to the equation, and
the migration time $T = L/v_0$, which belongs to the question. For a cell the
first is a microsecond and the second is about twenty minutes. Their ratio,
$10^{-9}$, drives everything that follows.
:::

## The model

For a particle in vacuum, $F = ma$. A cell moves through a viscous medium
which resists it:

$$
m\,a = F - F_\mathrm{drag} .
$$ (eq-newton-drag)

Here $F$ is whatever force the cell generates for itself, $F_\mathrm{drag}$ is
the resistance of the medium, and $ma$ is the inertial force. At
$\mathrm{Re} \approx 10^{-8}$ ({prf:ref}`ex-reynolds`) the flow around the
cell is a Stokes flow, and the drag is linear in the velocity,

$$
F_\mathrm{drag} = \gamma\,v, \qquad \gamma = 6\pi\eta a ,
$$ (eq-stokes-drag)

for a sphere of radius $a$ in a medium of viscosity $\eta$.

:::{prf:definition} Equation of motion
:label: def-eom

$$
m\,\ddot x = F(t) - \gamma\,\dot x
\qquad\Longleftrightarrow\qquad
m\,\dot v = F(t) - \gamma\,v, \quad v = \dot x .
$$ (eq-eom-full)

Second order in $x$, so it needs $x(0)$ **and** $\dot x(0)$. Since $x$ does
not appear on the right it is first order in $v$, and $x$ follows by one more
integration.
:::

We work with the equation for the velocity,

$$
m\,\frac{\dd v}{\dd t} = F(t) - \gamma\,v ,
\qquad v(0) = v_0 .
$$ (eq-eom-v)

No value of $F$ is needed anywhere in this part: the whole argument is about
time scales.

## Units

| Symbol | Meaning | SI units | At cell scale |
|---|---|---|---|
| $t$ | time | $\mathrm{s}$ | s, min |
| $x$ | position | $\mathrm{m}$ | $\mu$m |
| $v = \dot x$ | velocity | $\mathrm{m\,s^{-1}}$ | $\mu$m/min |
| $a = \ddot x$ | acceleration | $\mathrm{m\,s^{-2}}$ | never measured |
| $m$ | mass | $\mathrm{kg}$ | pg |
| $F$ | force | $\mathrm{N} = \mathrm{kg\,m\,s^{-2}}$ | pN |
| $\gamma$ | drag coefficient | $\mathrm{N\,s\,m^{-1}} = \mathrm{kg\,s^{-1}}$ | pN s $\mu\mathrm{m}^{-1}$ |
| $\eta$ | viscosity | $\mathrm{Pa\,s} = \mathrm{kg\,m^{-1}s^{-1}}$ | water: $10^{-3}$ |

Every term of {eq}`eq-eom-full` is a force:

$$
\underbrace{\mathrm{kg}\cdot\mathrm{m\,s^{-2}}}_{m\ddot x}
\;=\;
\underbrace{\mathrm{N}}_{F}
\;=\;
\underbrace{\mathrm{kg\,s^{-1}}\cdot\mathrm{m\,s^{-1}}}_{\gamma\dot x}
\;=\; \mathrm{kg\,m\,s^{-2}} . \qquad\checkmark
$$ (eq-dimensional-check)

Both spellings of the units of $\gamma$ are useful. As $\mathrm{N\,s\,m^{-1}}$
it is a force per unit velocity, which is what $\gamma v = F$ says; as
$\mathrm{kg\,s^{-1}}$ the kilograms cancel against $m$ and leave seconds.

Nothing in the equation is a time except $t$ itself. The units say which
combinations of the constants, and of the cell's own length $L$ and speed
$v_0$, are times.

:::{prf:example} Which combinations are times
:label: ex-combinations-are-times

**From the equation.** The only combination of $m$ and $\gamma$ with units of
time is

$$
\frac{m}{\gamma}
= \frac{\mathrm{kg}}{\mathrm{kg\,s^{-1}}}
= \mathrm{s} .
$$

**From the biology.** A length the cell has and a speed at which it moves give

$$
\frac{L}{v_0}
= \frac{\mathrm{m}}{\mathrm{m\,s^{-1}}}
= \mathrm{s} .
$$

Two time scales are visible from the units alone, before anything is solved.
Everything that follows is about their ratio.
:::

## Switch the forcing off

Let the cell stop generating traction at $t = 0$, so $F(t) = 0$ for $t > 0$.
Then {eq}`eq-eom-v` is separable:

$$
m\,\frac{\dd v}{\dd t} = -\gamma\,v
\qquad\Longrightarrow\qquad
\frac{\dd v}{v} = -\frac{\gamma}{m}\,\dd t
\qquad\Longrightarrow\qquad
v(t) = v(0)\,e^{-\gamma t/m} .
$$ (eq-free-decay)

The combination $m/\gamma$ flagged by the units now appears as the decay time
of the velocity: the dynamics picked it out.

:::{prf:definition} Inertial stopping time
:label: def-tau-i

$$
\tau_i = \frac{m}{\gamma},
\qquad
v(t) = v(0)\,e^{-t/\tau_i} .
$$ (eq-tau-i)

The time over which a cell that stops pushing loses its velocity.
:::

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_INNER, C_COMP = "#111111", "#D55E00", "#009E73"

fig, ax = plt.subplots(figsize=(5.4, 2.9))
x = np.linspace(0, 5, 800)
ax.plot(x, np.exp(-x), color=C_EXACT, lw=2.0, zorder=3)
xt = np.linspace(0, 1.25, 50)                       # tangent at the origin
ax.plot(xt, 1 - xt, color="#888888", ls=(0, (4, 3)), lw=1.1)
ax.annotate("initial slope meets\nzero at $t=\\tau_i$", xy=(1.0, 0.0), xytext=(2.25, 0.20),
            fontsize=8, color="#666666", arrowprops=dict(arrowstyle="->", color="#888888"))
half = np.log(2.0)                                  # half-life
ax.plot([0, half], [0.5, 0.5], color=C_INNER, ls=":", lw=1.3)
ax.plot([half, half], [0, 0.5], color=C_INNER, ls=":", lw=1.3)
ax.plot([half], [0.5], "o", color=C_INNER, ms=5, zorder=4)
ax.annotate(r"half-life  $t_{1/2}=\tau_i\ln 2 \approx 0.69\,\tau_i$", xy=(half, 0.5),
            xytext=(1.15, 0.70), fontsize=8, color=C_INNER,
            arrowprops=dict(arrowstyle="->", color=C_INNER))
ax.plot([1], [np.exp(-1)], "o", color=C_COMP, ms=5, zorder=4)   # the 1/e point
ax.annotate(r"$1/e \approx 0.37$ at $t=\tau_i$", xy=(1.0, np.exp(-1)), xytext=(1.9, 0.47),
            fontsize=8, color=C_COMP, arrowprops=dict(arrowstyle="->", color=C_COMP))
for k, lab in ((3, r"$5\%$"), (5, r"$0.7\%$")):
    ax.plot([k], [np.exp(-k)], "o", color="#999999", ms=3.5, zorder=4)
    ax.annotate(lab, xy=(k, np.exp(-k)), xytext=(k - 0.1, np.exp(-k) + 0.07),
                fontsize=7.5, color="#777777")
ax.set_xlabel(r"$t/\tau_i$"); ax.set_ylabel(r"$v(t)\,/\,v(0)$")
ax.set_xlim(-0.08, 5.1); ax.set_ylim(-0.03, 1.06)
ax.set_title(r"free decay after the forcing is switched off:  $v = v(0)\,e^{-t/\tau_i}$");
```

The free decay {eq}`eq-free-decay`, with three equivalent ways of reading
$\tau_i$ off it: the **half-life** $t_{1/2} = \tau_i\ln 2 \approx 0.69\,\tau_i$;
the time to fall to $1/e \approx 37\%$ of the starting speed, which is $\tau_i$
itself; and the **initial slope** extended until it meets zero, also at
$\tau_i$. After $3\tau_i$ the cell keeps $5\%$ of its speed; after $7\tau_i$,
one part in a thousand. Since the axis is in units of $\tau_i$, nothing in
this picture is large or small yet; that judgement needs a number for $\tau_i$
and something to compare it with.

## Put numbers in

A cell is mostly water. Take its radius as $a = 5\ \mu$m and the medium ten
times more viscous than water.

| Symbol | Meaning | Value |
|---|---|---|
| $a$ | cell radius | $5\ \mu\mathrm{m} = 5\times10^{-6}$ m |
| $m$ | cell mass, $\tfrac43\pi a^3\rho$ | $\approx 10^{-12}$ kg |
| $\eta$ | medium viscosity | $10^{-2}$ Pa s |
| $\gamma = 6\pi\eta a$ | drag coefficient | $\approx 10^{-6}$ kg/s $= 10^{-6}$ N s/m |
| $L$ | cell diameter | $10\ \mu\mathrm{m} = 10^{-5}$ m |
| $v_0$ | migration speed | $0.1$ to $1\ \mu$m/min $\approx 10^{-8}$ m/s |

$$
\tau_i = \frac{m}{\gamma}
\approx \frac{10^{-12}\ \mathrm{kg}}{10^{-6}\ \mathrm{kg/s}}
= 10^{-6}\ \mathrm{s} .
$$ (eq-tau-i-number)

A microsecond must be compared with the time scale on which the behavior we
want to explain happens, not with the speed of a camera. That time scale has to
be constructed, and it is the second time of
{prf:ref}`ex-combinations-are-times`. Neither $L$ nor $v_0$ comes out of the
equation: $L$ is the cell's diameter, the natural yardstick for a question
about where a cell goes, and $v_0$ is the observed speed of its centroid,
$v_0 \sim 0.1$ to $1\ \mu\mathrm{m/min} \approx 10^{-8}$ m/s in order of
magnitude.

:::{prf:definition} Migration time
:label: def-migration-time

$$
T = \frac{L}{v_0}
\approx \frac{10^{-5}\ \mathrm{m}}{10^{-8}\ \mathrm{m/s}}
= 10^{3}\ \mathrm{s} \approx 17\ \mathrm{min} ,
$$ (eq-migration-time)

the time a cell takes to move one body length. It is built from two measured
quantities rather than measured itself.

$T$ is the time scale of the question: we want to predict where the cell is,
in body lengths, and that choice of observable is what put $L$ and $v_0$ into
it. $\tau_i$ belongs to the cell; $T$ belongs to us.
:::

A microsecond against tens of minutes:

$$
\frac{\tau_i}{T} \approx \frac{10^{-6}\ \mathrm{s}}{10^{3}\ \mathrm{s}}
= 10^{-9} .
$$

Tens of minutes is the time scale of everything we observe about a crawling
cell: crossing its own width, changing speed, turning. So far, though, two
numbers computed separately have been compared, which is not yet grounds for
deleting a term from the equation. [](scaling-the-equation.md) rewrites the
equation of motion so that the ratio $\tau_i/T$ appears in it explicitly.
