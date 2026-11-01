---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/one-dimensional-exercises.pdf
    id: one-dimensional-exercises-pdf
downloads:
  - id: one-dimensional-exercises-pdf
    title: Download PDF
---

# Exercises: Flows on the Line

The exercises follow the pages of this part. Some are calculations by hand,
some ask you to build a model from a description, and a few ask for a short
computation in Python (or any other language).

## Rates in, rates out

:::{exercise} Messenger RNA
:label: exr-fl-mrna

A gene in *E. coli* is transcribed at a constant rate of $2$ mRNA molecules
per minute, and each mRNA molecule is degraded with a half-life of $5$
minutes.

1. Write the model {eq}`eq-ir-production-decay` for the number $C(t)$ of mRNA
   molecules in the cell. Give $a$ and $b$ with their units.
2. How many mRNA molecules does the cell hold at steady state?
3. The gene is switched on at $t = 0$, with no mRNA present. How long does it
   take to reach $90\%$ of the steady state?
4. Show, from {eq}`eq-ir-solution`, that the time to cover half of the
   remaining distance to the steady state is always the half-life, whatever
   $a$ and the starting value are. Why does the production rate not enter?
:::

:::{exercise} A drug cleared two ways
:label: exr-fl-infusion

A drug is infused into the blood at a constant rate $I$ (mg/L per hour). The
kidneys remove it at a rate proportional to its concentration $C$, with rate
constant $k$, and a liver enzyme breaks it down at the saturating rate
$V C/(K_M + C)$.

1. Write the word equation and the model for $C(t)$. Give the units of $I$,
   $k$, $V$ and $K_M$.
2. Explain each loss term by its mechanism, in one sentence each. Which one
   matters more at low concentrations, and which at high ones?
3. Show, without solving for it, that there is exactly one steady state
   $C^\ast > 0$ and that it is stable. (Hint: what is the sign of the
   derivative of the right-hand side?)
4. The infusion rate is doubled. Does the steady-state concentration more than
   double, double, or less than double? Answer from a sketch of the two loss
   terms against $C$.
:::

:::{exercise} Flagellar length
:label: exr-fl-flagella

The alga *Chlamydomonas* keeps its two flagella at about $10\ \mu$m, although
they are taken apart and rebuilt all the time. Material is carried to the tip
by a fixed number of transport trains that shuttle up and down the flagellum
at constant speed, and the tip falls apart at a constant speed $B$. A longer
flagellum means a longer round trip, so each train visits the tip less often.
This gives the model of Marshall and Rosenbaum {cite:p}`marshall2001`,

$$
\frac{\dd L}{\dd t} = \frac{A}{L} - B, \qquad L > 0 .
$$

1. What are the units of $A$ and $B$? A flagellum that is no longer supplied
   shortens by $10\ \mu$m in about $15$ min. Estimate $B$, and then $A$ from
   the steady length $L_\infty = 10\ \mu$m.
2. Find the steady state, show that it is stable, and find the relaxation time
   $1/|f'(L_\infty)|$ in terms of $L_\infty$ and $B$.
3. With $\ell = L/L_\infty$ and $s = Bt/L_\infty$, show that
   $\dd\ell/\dd s = 1/\ell - 1$. What does it mean that no parameter is left?
4. Separate variables to show that the solution with $\ell(0) = 0$ satisfies
   $s = -\ell - \ln(1 - \ell)$. How many minutes does a flagellum that was cut
   off take to regrow to $90\%$ of its length?
5. The model says a flagellum of length zero grows at infinite speed. Which
   assumption fails for a very short flagellum?
:::

:::{exercise} Checking for artifacts
:label: exr-fl-artifacts

Each model below is proposed for a nonnegative quantity. Check it against
conventions 3 and 5 of [](index.md#modeling-conventions): can the quantity
become negative or blow up, and does it have a stable steady state? If it
fails, say what is wrong with the mechanism and propose a fix.

1. A protein that is made at rate $a$ and removed in pairs (it dimerizes and
   the dimer is degraded): $\dot C = a - bC^2$.
2. A fish population that is harvested at a fixed rate:
   $\dot N = rN(1 - N/K) - H$.
3. A population whose cells each divide at rate $\beta$ and are eaten by
   predators at the fixed total rate $P$: $\dot N = \beta N - P$.
:::

## Populations and fractions

:::{exercise} A strong Allee effect
:label: exr-fl-allee

Some species grow only above a minimum density: below it, individuals fail to
find mates or to defend themselves. A model is

$$
\dot N = f(N) = rN\Big(\frac{N}{A} - 1\Big)\Big(1 - \frac{N}{K}\Big),
\qquad 0 < A < K .
$$

1. Find the steady states and classify them with the slope test.
2. Draw the phase line. What happens to a population introduced at a density
   below $A$, and to one introduced above?
3. Show that the per-capita growth rate $\dot N/N$ is largest at $N =
   (A + K)/2$. Compare with the tolerated dose $h(N)$ of the inoculum model,
   {eq}`eq-bf-tolerance`: what plays the role of the drug, and what is the
   analogue of $N_-$?
:::

:::{exercise} The rare type wins
:label: exr-fl-frequency

In {prf:ref}`prop-pf-fraction` the growth rates may depend on the fraction
$p$. Suppose two types compete for two foods, and each type is better at
using one of them. The rarer type then has more of its food to itself. A
simple model is $r_1 - r_2 = s\,(1 - 2p)$ with $s > 0$.

1. Write the equation for $p$, and find and classify its steady states.
2. Interpret the result. Contrast it with {prf:ref}`cor-pf-takeover`.
3. In {prf:ref}`ex-pf-resistance`, with constant rates and
   $r_1 - r_2 = -0.1$/h, how long does the resistant fraction take to fall
   from $99\%$ to $1\%$?
:::

## Bifurcations

:::{exercise} Critical slowing down in numbers
:label: exr-fl-slowing

In the shared-good model {eq}`eq-bf-shared` take $c = 0.1$/h. Using
{eq}`eq-bf-slopes`, find the doses $a$ at which the time scale $1/|f'|$ of the
stable steady state is one day, on each side of the bifurcation value
$a = c$. How long is the time scale at $a = 0.11$/h? What would an
experimenter who measured the resistant fraction after one day of treatment
conclude at that dose?
:::

:::{exercise} Two saddle-nodes
:label: exr-fl-exam

Consider $\dot x = 2 + \mu x + x^2$ with the parameter $\mu$.

1. Find all steady states and their stability.
2. Find all bifurcation points, using the recipe of
   {prf:ref}`rem-bf-slope-fails`: $f = 0$ and $f_x = 0$ together. What kind of
   bifurcation are they?
3. Draw the phase lines for $|\mu| < 2\sqrt2$, $|\mu| = 2\sqrt2$ and
   $|\mu| > 2\sqrt2$.
4. Draw the bifurcation diagram. (Solve the steady-state equation for the
   parameter, $\mu = -(x^2 + 2)/x$, and plot $\mu$ against $x$ with the axes
   swapped, as in [](saddle-node.md).)
:::

:::{exercise} Inoculum and dose
:label: exr-fl-inoculum

Take the inoculum model {eq}`eq-bf-inoculum` with $r = 1$/day, $K = 10$ and
$N_c = 1$ (in $10^8$ cells/mL), at the dose $a = 2$/day.

1. Find the threshold inoculum $N_-$ and the established state $N_+$ by
   solving $h(N) = a$.
2. Use {eq}`eq-bf-slope-h` to find the rates at $N_-$ and $N_+$, and the
   corresponding time scales.
3. An infection starts with $2 \times 10^8$ cells/mL. Is it cleared at this
   dose? What is the smallest dose that clears it?
4. A patient can tolerate at most $a = 2.8$/day. Which inocula can be cleared?
:::

:::{exercise} The pitchfork
:label: exr-fl-pitchfork

Many systems are symmetric under $x \mapsto -x$: a cell that can polarize to
the left or to the right, or a switch between two equivalent states. The
simplest such equation with a bifurcation is

$$
\dot x = \mu x - x^3 .
$$

1. Find the steady states and their stability for $\mu < 0$ and $\mu > 0$, and
   draw the bifurcation diagram. Why is this bifurcation called a
   **pitchfork**?
2. Near $\mu = 0$, how does the size of the nonzero steady states grow with
   $\mu$? Compare with the transcritical bifurcation
   ({prf:ref}`prop-bf-transcritical`), where the new state grew like
   $\mu - \mu_c$.
3. A small bias breaks the symmetry: $\dot x = h + \mu x - x^3$ with $h > 0$
   small. Fix $\mu > 0$ and use the graphical method: the steady states are
   where the graph of $x^3 - \mu x$ meets the horizontal line at height $h$.
   For which $h$ are there three steady states? (At $\mu = 1$ the answer is
   $|h| < 0.385$.)
4. Turning $h$ slowly up and down at fixed $\mu > 0$ gives hysteresis. Explain
   why, using the bifurcation diagram in $h$.
:::

