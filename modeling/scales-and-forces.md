---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/scales-and-forces.pdf
    id: modeling-scales-and-forces-pdf
downloads:
  - id: modeling-scales-and-forces-pdf
    title: Download PDF
---

# Sizes, Speeds, and Forces

:::{tip} Big Idea
Before any equation, get the numbers. Every reduction is justified by a ratio
of two entries in one table: two times, two forces, or two lengths. When two
entries of the same kind are far apart, something can be dropped; when they
are close, nothing can. The first such ratio is the Reynolds number, and for a
crawling cell it is $10^{-8}$.
:::

## The table

| Quantity | Symbol | Typical value |
|---|---|---|
| Cell diameter | $L$ | $10$ to $20\ \mu\mathrm{m}$ |
| Migration speed | $u_0$ | $0.1$ to $1\ \mu\mathrm{m}/\mathrm{min}$ |
| Protrusion length | $\ell$ | $1$ to $10\ \mu\mathrm{m}$ |
| Protrusion lifetime | $T$ | $1$ to $10\ \mathrm{min}$ |
| Single-motor stall force | $F_s$ | $\sim 5\ \mathrm{pN}$ |
| Traction force of a cell | $F$ | $10$ to $100\ \mathrm{nN}$ |
| Cytoplasmic viscosity | $\eta$ | $10$ to $1000 \times$ water |
| Drag coefficient | $\gamma$ | see below |

Nothing in this table is exact, and nothing needs to be. Each entry is an
order of magnitude taken from observation, not from a model; what a model does
is compare them.

Two entries are times, and a third can be built from a length and a speed:

$$
\tau_\mathrm{mig} = \frac{L}{u_0}
\approx \frac{10\ \mu\mathrm{m}}{0.5\ \mu\mathrm{m/min}} = 20\ \mathrm{min} ,
$$

the time for a cell to move its own length. Compared with the protrusion
lifetime $T \approx 1$ to $10$ min, this says a cell replaces its entire set of
protrusions a few times while translating one body length. These two times are
close, so neither process can be treated as fast relative to the other.

:::{warning} Two uses of $T$
In this table $T$ is the protrusion lifetime. In
[](../overdamped/index.md), $T = L/v_0$ is the migration time, the quantity
called $\tau_\mathrm{mig}$ here, and the migration speed is written $v_0$.
:::

## Life at low Reynolds number

The Reynolds number compares inertia with viscosity:

$$
\mathrm{Re} = \frac{\rho u L}{\eta} .
$$

:::{prf:example} Reynolds number of a crawling cell
:label: ex-reynolds

With $\rho \approx 10^3\ \mathrm{kg/m^3}$, $u \approx 10^{-8}\ \mathrm{m/s}$
($0.6\ \mu$m/min), $L \approx 10^{-5}$ m and $\eta \approx 10^{-2}$ Pa$\cdot$s
(ten times water),
$$
\mathrm{Re} \approx \frac{10^3 \cdot 10^{-8} \cdot 10^{-5}}{10^{-2}} = 10^{-8}.
$$
:::

For comparison, Purcell's swimming bacterium {cite}`purcell1977` has
$\mathrm{Re} \approx 10^{-4}$ and coasts about a tenth of an ångström after it
stops swimming. A crawling cell is four orders of magnitude further into the
viscous world.

:::{prf:proposition} Inertia is irrelevant
:label: prop-inertia

For a particle of mass $m$ in a fluid with drag coefficient $\gamma$, Newton's
law $m\ddot x = -\gamma\dot x$ gives $\dot x(t) = \dot x(0)\,e^{-t/\tau_i}$
with the **inertial stopping time** $\tau_i = m/\gamma$. For a crawling cell
$\tau_i \approx 10^{-6}$ s, eight orders of magnitude below the shortest
biological time in the table.
:::

:::{prf:proof}
:class: dropdown

The equation $m\,\dd\dot x/\dd t = -\gamma\dot x$ is separable, with solution
$\dot x(t) = \dot x(0)\,e^{-\gamma t/m}$. For the estimate, take
$m \approx 10^{-12}$ kg (a $10\ \mu$m cell of water density) and estimate
$\gamma$ from Stokes drag on a sphere of radius $a = 5\ \mu$m in a medium of
viscosity $\eta = 10^{-2}$ Pa$\cdot$s:
$$
\gamma = 6\pi\eta a \approx 6\pi \cdot 10^{-2} \cdot 5\times10^{-6}
       \approx 10^{-6}\ \mathrm{kg/s}.
$$
Hence $\tau_i = m/\gamma \approx 10^{-12}/10^{-6} = 10^{-6}$ s. The shortest
biological time scale in the table is the protrusion lifetime,
$T \sim 10^2$ s.
:::

The consequence is the overdamped force balance

$$
\gamma\,\dot{\mathbf x} = \mathbf F(\mathbf x, t) :
$$

force sets velocity, not acceleration. [](../overdamped/index.md) derives it
carefully from Newton's law with inertia included, and computes what dropping
inertia costs. Whenever a cell track looks as though it has inertia, and it
often does, the memory is in an internal state of the cell, not in its mass.

## Forces

A single myosin or actin-polymerization motor stalls at a few piconewtons. A
cell exerts tens of nanonewtons of traction, so the number of motors acting
together is

$$
\frac{F}{F_s} \approx \frac{10^{-8}\ \mathrm{N}}{5\times10^{-12}\ \mathrm{N}}
\approx 10^{4} .
$$

Ten thousand motors act together in a cell, but only about ten protrusions.
The first number is large enough that averaging over motors gives a meaningful
deterministic force-velocity relation for a protrusion. The second is not: ten
is too few to average over, so protrusions must be kept as discrete objects.
