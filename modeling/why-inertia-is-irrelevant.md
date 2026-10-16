---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/why-inertia-is-irrelevant.pdf
    id: modeling-why-inertia-is-irrelevant-pdf
downloads:
  - id: modeling-why-inertia-is-irrelevant-pdf
    title: Download PDF
kernelspec:
  name: python3
  display_name: Python 3
---

# Doing Some Modeling: Why Inertia Is Irrelevant

:::{tip} Big Idea
A modeling argument follows the same moves almost every time: sharpen the
question, start from something trusted, get the numbers, solve the simplest
case, compare the answer with the scales of the question, and read off what it
means. Here they are applied once, to a crawling cell that stops pulling.
Newton's second law with drag says it coasts for about a microsecond and
travels about ten femtometres. Both numbers are $10^{-9}$ of the cell's own
scales, the twenty minutes it takes to move one body length and the body
length itself. So except for a few microseconds after each change in force, a
cell's velocity is proportional to the force on it: $v \approx F/\gamma$.
:::

This page asks one question about a cell and answers it with numbers. The six
steps below are not specific to cells; they are the moves of almost every
modeling argument.

## 1. Sharpen the question

How a cell gets around a corner ([](what-is-a-model.md)) cannot be attacked
yet, because before asking how a cell *steers* we have to know what sets its
motion at all. A swimmer who stops stroking glides on for a body length or
two: the velocity carries memory of the force that produced it. Does a cell's?
So the question is the one underneath:

> If a crawling cell stops pulling, how long does it keep moving, and how far
> does it get?

## 2. Start from something trusted

For a particle in vacuum, $F = ma$. A cell moves through a viscous medium,
which resists it. Newton's second law for the cell center, with an active
force, a contact force, and viscous drag, is

$$
m\ddot{\mathbf x} = \mathbf F_\mathrm{active} + \mathbf F_\mathrm{contact} -
    \gamma\dot{\mathbf x}.
$$ (eq-newton-cell)

The drag is linear in the velocity, with drag coefficient $\gamma$; step 5
shows why that is the right form for a cell.

## 3. Get the numbers

Before solving anything, estimate everything. For a crawling cell:

| Quantity | Symbol | Typical value |
|---|---|---|
| Cell diameter | $L$ | $10$ to $20\ \mu\mathrm{m}$ |
| Migration speed | $v_0$ | $0.1$ to $1\ \mu\mathrm{m}/\mathrm{min}$ |
| Single-motor stall force | $F_s$ | $\sim 5\ \mathrm{pN}$ |
| Traction force of a cell | $F_\mathrm{tr}$ | $10$ to $100\ \mathrm{nN}$ |
| Cytoplasmic viscosity | $\eta$ | $10$ to $1000 \times$ water |

Nothing in this table is exact, and nothing needs to be. Each entry is an
order of magnitude taken from observation, not from a model; what a model does
is compare them. Two of the entries are already worth comparing. A single
myosin or actin-polymerization motor stalls at a few piconewtons, and a cell
exerts tens of nanonewtons of traction, so the number of motors acting
together is

$$
\frac{F_\mathrm{tr}}{F_s}
\approx \frac{10^{-8}\ \mathrm{N}}{5\times10^{-12}\ \mathrm{N}}
\approx 10^{4} .
$$

Ten thousand motors act together in a cell, but only about ten protrusions.
The first number is large enough that averaging over motors gives a meaningful
deterministic force-velocity relation for a protrusion. The second is not: ten
is too few to average over, so protrusions must be kept as discrete objects.

The two coefficients in {eq}`eq-newton-cell` follow from the table. A cell is
mostly water, so a cell of diameter $10\ \mu$m has mass
$m \approx 10^{-12}$ kg. For the drag, use Stokes drag on a sphere of radius
$a = 5\ \mu$m in a medium ten times as viscous as water,
$\eta = 10^{-2}$ Pa$\cdot$s:

$$
\gamma = 6\pi\eta a \approx 6\pi \cdot 10^{-2} \cdot 5\times10^{-6}
       \approx 10^{-6}\ \mathrm{kg/s}.
$$ (eq-stokes-drag)

## 4. Solve the simplest case

Let a cell moving at the migration speed $v_0$ stop pulling at $t = 0$. In one
dimension, with $v = \dot x$, {eq}`eq-newton-cell` becomes
$m\,\dd v/\dd t = -\gamma v$ with $v(0) = v_0$, which is separable:

$$
\frac{\dd v}{v} = -\frac{\gamma}{m}\,\dd t
\qquad\Longrightarrow\qquad
v(t) = v_0\,e^{-\gamma t/m} .
$$ (eq-free-decay)

The constants $m$ and $\gamma$ have combined into a time.

:::{prf:definition} Inertial stopping time
:label: def-tau-i

$$
\tau_i = \frac{m}{\gamma},
\qquad
v(t) = v_0\,e^{-t/\tau_i} .
$$ (eq-tau-i)

The time over which a cell that stops pushing loses its velocity: its
**coasting time**. For a crawling cell

$$
\tau_i \approx \frac{10^{-12}\ \mathrm{kg}}{10^{-6}\ \mathrm{kg/s}}
= 10^{-6}\ \mathrm{s} .
$$ (eq-tau-i-number)
:::

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
C_EXACT, C_INNER, C_COMP, C_OUTER = "#111111", "#D55E00", "#009E73", "#0072B2"

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
ax.set_xlabel(r"$t/\tau_i$"); ax.set_ylabel(r"$v(t)\,/\,v_0$")
ax.set_xlim(-0.08, 5.1); ax.set_ylim(-0.03, 1.06)
ax.set_title(r"free decay after the forcing is switched off:  $v = v_0\,e^{-t/\tau_i}$");
```

The free decay {eq}`eq-free-decay`, with three equivalent ways of reading
$\tau_i$ off it: the **half-life** $t_{1/2} = \tau_i\ln 2 \approx 0.69\,\tau_i$;
the time to fall to $1/e \approx 37\%$ of the starting speed, which is $\tau_i$
itself; and the **initial slope** extended until it meets zero, also at
$\tau_i$. After $3\tau_i$ the cell keeps $5\%$ of its speed; after $7\tau_i$,
one part in a thousand.

The cell is not stopped at $\tau_i$, only slowed by a factor $e$; strictly it
never stops. The distance it covers is nevertheless finite, because the speed
decays exponentially. Integrating {eq}`eq-tau-i` from the moment the cell
stops pulling,

$$
x(t) - x(0) = \int_0^t v_0\,e^{-s/\tau_i}\,\dd s
= v_0\,\tau_i\left(1 - e^{-t/\tau_i}\right) .
$$ (eq-coasting-position)

:::{prf:example} Coasting distance of a crawling cell
:label: ex-coasting-cell

As $t \to \infty$, {eq}`eq-coasting-position` tends to the **coasting
distance**

$$
x_\mathrm{coast} = v_0\,\tau_i
\approx \bigl(10^{-8}\ \mathrm{m/s}\bigr)\bigl(10^{-6}\ \mathrm{s}\bigr)
= 10^{-14}\ \mathrm{m} ,
$$ (eq-coasting-distance)

ten femtometres, about the diameter of a heavy atomic nucleus and $10^{-4}$ of
the diameter of a hydrogen atom. Here $v_0 \approx 0.6\ \mu$m/min
$\approx 10^{-8}$ m/s.
:::

```{code-cell} python
:tags: [hide-input]
V0, TAU = 1e-8, 1e-6                     # m/s, s
t = np.linspace(0, 6e-6, 600)
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
ax[0].plot(t * 1e6, V0 * np.exp(-t / TAU) * 6e7, color=C_EXACT, lw=1.8)   # m/s -> um/min
ax[0].set_xlabel(r"time after traction stops ($\mu$s)")
ax[0].set_ylabel(r"speed ($\mu$m/min)")
ax[0].set_title("coasting time: about a microsecond")
ax[0].set_ylim(0, 0.66)
ax[1].plot(t * 1e6, V0 * TAU * (1 - np.exp(-t / TAU)) * 1e15, color=C_COMP, lw=1.8)
ax[1].axhline(V0 * TAU * 1e15, color=C_COMP, ls="--", lw=1.0)
ax[1].annotate(r"$x_{\rm coast} = v_0\tau_i \approx 10$ fm", xy=(2.2, 10.4),
               fontsize=8, color=C_COMP)
ax[1].set_xlabel(r"time after traction stops ($\mu$s)")
ax[1].set_ylabel("distance coasted (fm)")
ax[1].set_title("coasting distance: ten femtometres")
ax[1].set_ylim(0, 12.5)
fig.tight_layout();
```

The same coast in physical units. **Left:** a cell moving at
$0.6\ \mu$m/min loses its speed within a few microseconds. **Right:** the
distance it covers meanwhile, {eq}`eq-coasting-position`, levels off at
$x_\mathrm{coast} = v_0\tau_i \approx 10$ fm.

This answers the question of step 1: the cell keeps moving for about a
microsecond and gets about ten femtometres.

## 5. Compare with the scales of the question

A microsecond and ten femtometres are small only compared with something. They
must be compared with the scales on which the behavior we want to explain
happens, not with the resolution of a microscope. The length is the cell's own
diameter $L$, and the time is built from a length and a speed in the table.

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

Measured against these, the coasting time and the coasting distance are the
same number:

$$
\frac{\tau_i}{T} \approx \frac{10^{-6}\ \mathrm{s}}{10^{3}\ \mathrm{s}}
= 10^{-9},
\qquad
\frac{x_\mathrm{coast}}{L} = \frac{v_0\tau_i}{L} = \frac{\tau_i}{L/v_0}
= \frac{\tau_i}{T} .
$$ (eq-inertia-ratio)

This is not a coincidence of the numbers: the cell coasts for a fraction
$\tau_i/T$ of the time it takes to move one body length, at no more than its
usual speed, so it covers at most that fraction of a body length. The same
ratio appears term by term in {eq}`eq-newton-cell`: a cell moving at $v_0$ and
changing its speed over $T$ has inertia term $m v_0/T$ and drag term
$\gamma v_0$, whose ratio is again $\tau_i/T$.

```{code-cell} python
:tags: [hide-input]
C_GREY = "#666666"
fig, ax = plt.subplots(2, 1, figsize=(6.6, 3.4))
rows = [
    (ax[0], "time (s)", (1e-7, 1e4),
     [(1e-6, r"coasting time $\tau_i$", C_INNER, 0.5),
      (60, "one frame of\na cell movie", C_GREY, 0.5),
      (1e3, r"migration time $T$", C_OUTER, 1.5)]),
    (ax[1], "length (m)", (1e-15, 1e1),
     [(1e-14, r"cell coasts, $v_0\tau_i$", C_INNER, 0.5),
      (1e-11, "bacterium coasts", C_GREY, 1.5),
      (1e-10, "hydrogen atom", C_GREY, 0.5),
      (1e-5, r"cell diameter $L$", C_OUTER, 0.5),
      (2, "swimmer glides", C_GREY, 1.5)]),
]
for a, lab, (lo, hi), marks in rows:
    a.set_xscale("log"); a.set_xlim(lo, hi); a.set_ylim(-1.3, 2.3)
    a.set_yticks([]); a.spines["left"].set_visible(False)
    a.axhline(0, color="#999999", lw=0.8, zorder=1)
    a.set_xlabel(lab)
    for xm, text, c, h in marks:
        a.plot([xm, xm], [0, h - 0.1], color=c, lw=0.7, zorder=2)
        a.plot([xm], [0], "o", color=c, ms=6, zorder=3)
        a.text(xm, h, text, ha="center", va="bottom", fontsize=7.5, color=c)
    b0, b1 = marks[0][0], [m for m in marks if m[2] == C_OUTER][0][0]
    a.annotate("", xy=(b1, -0.7), xytext=(b0, -0.7),
               arrowprops=dict(arrowstyle="<->", color="#444444", lw=0.9))
    a.text(np.sqrt(b0 * b1), -0.7, " nine orders of magnitude ", ha="center",
           va="center", fontsize=7.5, color="#444444", backgroundcolor="white")
fig.tight_layout();
```

The coast against the scales of the question, on logarithmic axes. **Top:**
the coasting time $\tau_i$ lies nine orders of magnitude below the migration
time $T$, and eight below a single frame of a typical cell movie. **Bottom:**
the coasting distance lies nine orders of magnitude below the cell diameter
$L$ and four below the size of an atom. A bacterium coasts about a tenth of an
ångström {cite}`purcell1977`; a swimmer glides a metre or two.

The fluid gives the same verdict. The **Reynolds number**

$$
\mathrm{Re} = \frac{\rho u L}{\eta}
$$

compares inertia with viscosity in the flow around a body of size $L$ moving
at speed $u$ through a fluid of density $\rho$.

:::{prf:example} Reynolds number of a crawling cell
:label: ex-reynolds

With $\rho \approx 10^3\ \mathrm{kg/m^3}$, $u \approx 10^{-8}\ \mathrm{m/s}$
($0.6\ \mu$m/min), $L \approx 10^{-5}$ m and $\eta \approx 10^{-2}$ Pa$\cdot$s
(ten times water),
$$
\mathrm{Re} \approx \frac{10^3 \cdot 10^{-8} \cdot 10^{-5}}{10^{-2}} = 10^{-8}.
$$
:::

At such a small Reynolds number the flow around the cell is a Stokes flow, and
that is why the drag in {eq}`eq-newton-cell` is linear in the velocity, with
the coefficient {eq}`eq-stokes-drag`. Purcell's swimming bacterium has
$\mathrm{Re} \approx 10^{-4}$ {cite}`purcell1977`; a crawling cell is four
orders of magnitude further into the viscous world.

## 6. Read off what it means

Coasting is the case of zero force. The same calculation works with the force
on. Let the net force $F = F_\mathrm{active} + F_\mathrm{contact}$ be constant
from $t = 0$, and let the cell start at any speed $v(0)$. In one dimension
{eq}`eq-newton-cell` is $m\,\dd v/\dd t = F - \gamma v$, which is linear,
and the integrating factor $e^{t/\tau_i}$ ([](../overdamped/integrating-factor.md))
gives

$$
v(t) = \frac{F}{\gamma} + \Bigl(v(0) - \frac{F}{\gamma}\Bigr)e^{-t/\tau_i} .
$$ (eq-relax-to-force)

Whatever speed the cell had, it forgets it on the coasting time $\tau_i$ and
settles at $F/\gamma$; {eq}`eq-tau-i` is the case $F = 0$. A force that is
switched on, off, or changed starts a new transient of the same microsecond
length. Protrusion forces change over seconds to minutes, so the cell spends
all but a few microseconds of every second with its transient long finished.

A cell therefore has, for every practical purpose, no momentum. A cell that
stops pulling stops, to within ten femtometres. This matches what is known: a
cell whose motors are disrupted stops essentially instantly, and nobody has
filmed a cell coasting. Whenever a cell track looks as though it has inertia,
and it often does, the memory is in an internal state of the cell, not in its
mass.

The conclusion of the whole page is a single observation. **Except for a few
microseconds after each change in the force, the velocity of a cell is
proportional to the force acting on it now:**

$$
v \approx \frac{F}{\gamma} .
$$ (eq-v-proportional-f)

Force sets velocity, not acceleration. Turning this observation into an
equation of motion means deleting the term $m\ddot{\mathbf x}$, which changes
the order of the equation, and [](../overdamped/index.md) does this carefully
and computes how much error it commits.
