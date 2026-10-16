---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/what-is-a-model.pdf
    id: modeling-what-is-a-model-pdf
downloads:
  - id: modeling-what-is-a-model-pdf
    title: Download PDF
---

# What Is Mathematical Modeling?

:::{tip} Big Idea
A model is a simplification built for a question, made precise enough to
settle that question true or false, or precise enough to show what could work.
Change the question and you change the model. Every model has a domain of
validity, and knowing where it breaks is as much a result as anything it
predicts.
:::

## What is a model?

Ask how a crawling cell, with no brain, no eyes, and no map, finds its way
around a corner, and the answers come quickly: it follows a chemical signal, it
senses the wall, it moves randomly and gets lucky, it feels stiffness. Each
answer is a proposed mechanism stated informally. Turning such a hypothesis
into something that can be tested is what modeling is.

:::{prf:definition} Mathematical model
:label: def-model

A mathematical model is a deliberately simplified description of part of the
world, written in mathematical language, built to answer a specific question.
:::

Each part of the definition matters.

- **Deliberately simplified.** A model is not a copy of reality; a map at 1:1
  scale is useless. The simplification is the content of the model, not a
  regrettable compromise.
- **Mathematical language.** Ambiguity that survives in English does not
  survive in an equation. "The cell senses the wall" leaves open how it senses,
  over what time it integrates, and against what reference it compares.
- **A specific question.** There is no such thing as *the* model of a cell,
  only a model *for* a question. Change the question and you change the model.
  Forgetting this is the most common beginner error.

Modeling makes a proposed mechanism precise enough to settle it, true or
false, or, when it cannot be settled outright, precise enough to see which
mechanisms are even capable of producing the observed behavior. "The cell
senses the wall" does neither: it is compatible with every observation, so it
can be neither confirmed nor ruled out, and it says nothing about what would
suffice. Box's aphorism, "all models are wrong, but some are useful," points at
the right axis. Wrong is not the interesting question; *useful for what* is.

## Kinds of model

| Kind | What it does | Example |
|---|---|---|
| Descriptive | Fits and summarizes data | Logistic curve fit to tumor growth |
| Mechanistic | Encodes a hypothesized cause | Hodgkin and Huxley's model of the nerve impulse |
| Predictive | Forecasts, need not explain | A network predicting protein structure |
| Normative | What *should* happen under a criterion | Optimal foraging |

Only a mechanistic model can be falsified by an experiment on the mechanism.
Fitting a curve well is not evidence for a mechanism.

## Successes

:::{note} Why biology is a hard place to model
Physics offers exact, universal laws and parameters known to many digits;
biology offers neither.

| | Classical physics | Biology |
|---|---|---|
| Laws | Known, exact, universal | No comparable laws; rules are approximate and contextual |
| Parameters | Measured to many digits | Often unknown to an order of magnitude |
| Components | Identical particles | Every cell differs |
| Noise | Usually a nuisance | Often functional; the system uses it |
| Scales | Frequently separable | Coupled across seconds to days, nm to cm |
| History | Irrelevant | Everything is as it is because it evolved |

The last three rows are not obstacles to remove. Heterogeneity, noise, and
coupled scales are frequently the phenomenon itself. A crawling cell, for
instance, turns its protrusions over in minutes while taking about twenty
minutes to move its own length, so the two processes cannot be separated.
:::

Despite this, mathematical biology has a long list of successes, built by
identifiable people who were often ignored, often working on something else,
and occasionally far ahead of their evidence. Each carries its equation,
because the equations are the argument.

- **Malthus (1798).** Population grows geometrically, food arithmetically,
  with almost no data behind the claim:

  $$
  \frac{\dd N}{\dd t} = rN
  \qquad\Longrightarrow\qquad
  N(t) = N_0\,e^{rt} .
  $$

  $N$ is the number of individuals and $r$ the per-capita growth rate, births
  minus deaths per individual per unit time. Nothing in the equation can stop
  the growth. Darwin and Wallace both read Malthus independently and both
  arrived at natural selection: a crude demographic model triggered the central
  theory of biology.

- **Verhulst (1838).** He added the saturation term Malthus lacked and named
  the result the logistic:

  $$
  \frac{\dd N}{\dd t} = rN\left(1 - \frac{N}{K}\right) .
  $$

  $K$ is the carrying capacity, the population the environment can support.
  The new factor is a switch: near $N = 0$ it is one and growth is Malthusian,
  and it shuts growth off as $N \to K$. The correction was forgotten and
  rediscovered by Pearl and Reed in 1920.

- **Ronald Ross (1911).** Having shown that mosquitoes transmit malaria, he
  did the mathematics and found something no experiment could hand him:

  $$
  \frac{\dd I}{\dd t} = \beta I(1-I) - \gamma I,
  \qquad
  R_0 = \frac{\beta}{\gamma} .
  $$

  $I$ is the infected fraction of the population, $\beta$ the transmission
  rate, $\gamma$ the recovery rate. The disease grows from rare only if
  $R_0 > 1$, so control is a **threshold** and not a total: there is no need to
  kill every mosquito, only to push the population past the line.

- **Lotka (1925) and Volterra (1926).** When fishing in the Adriatic largely
  stopped during the First World War, the proportion of predatory fish went
  *up*. Why should less fishing favor sharks?

  $$
  \frac{\dd V}{\dd t} = aV - bVP - fV,
  \qquad
  \frac{\dd P}{\dd t} = c\,bVP - dP - fP .
  $$

  $V$ is prey and $P$ predators; $a$ is prey growth, $b$ the rate at which
  predators eat prey, $c$ how efficiently eaten prey becomes new predators,
  $d$ the predator death rate, and $f$ fishing, which kills **both**. The
  steady state is the whole answer:

  $$
  V^\ast = \frac{d+f}{cb}, \qquad P^\ast = \frac{a-f}{b} .
  $$

  Less fishing means smaller $f$, hence *fewer* prey and *more* predators. Each
  species is controlled by the other's parameters, which intuition does not
  offer.

- **Fisher, Haldane, and Wright (1918 to 1932).** They showed that Mendelian
  genetics and Darwinian selection, widely believed incompatible, were the
  same theory:

  $$
  \frac{\dd p}{\dd t} = s\,p(1-p) .
  $$

  $p$ is the frequency of one allele in the population and $s$ its selection
  coefficient, the fitness advantage it confers. Discrete inheritance and
  continuous selection fit in one line.

- **Hodgkin and Huxley (1952).** The squid giant axon is about a millimeter
  across, wide enough to push a wire down the inside:

  $$
  C\frac{\dd V}{\dd t} = I - \bar g_{\rm Na}m^3h\,(V - E_{\rm Na})
     - \bar g_{\rm K}n^4(V - E_{\rm K}) - \bar g_{\rm L}(V - E_{\rm L}) .
  $$

  $V$ is the voltage across the membrane, $C$ its capacitance, $I$ injected
  current. Each remaining term is a current through one type of ion channel,
  with $\bar g$ its maximum conductance and $E$ the voltage at which it
  reverses. The gating variables $m$, $h$, $n$ are fractions of open channels,
  each obeying its own ODE, so the model has four equations. Huxley integrated
  them by hand on a mechanical calculator over roughly three weeks; the work
  won the 1963 Nobel Prize.

- **Turing (1952).** "The Chemical Basis of Morphogenesis":

  $$
  \partial_t u = D_u\nabla^2 u + f(u,v),
  \qquad
  \partial_t v = D_v\nabla^2 v + g(u,v) .
  $$

  $u$ and $v$ are concentrations of two chemicals, $D_u$ and $D_v$ their
  diffusion coefficients, and $f$, $g$ their reactions. Diffusion, the process
  that smooths concentrations out, can, when $D_v \gg D_u$, *destabilize* a
  uniform state and produce stripes and spots: local activation with
  long-range inhibition.

- **Keller and Segel (1970).** Slime mold aggregation, with chemotaxis as an
  instability:

  $$
  \partial_t u = \nabla\cdot\bigl(D\nabla u - \chi\,u\nabla c\bigr),
  \qquad
  \partial_t c = D_c\nabla^2 c + \alpha u - \beta c .
  $$

  $u$ is the density of cells and $c$ the concentration of the chemical they
  follow. The first term spreads cells out; the second, with chemotactic
  sensitivity $\chi$, moves them *up* the gradient of $c$, which by the second
  equation they produce themselves at rate $\alpha$. Cells attract cells.

- **Purcell (1977).** *Life at Low Reynolds Number* {cite}`purcell1977`:

  $$
  \mathrm{Re} = \frac{\rho u L}{\eta} .
  $$

  $\rho$ is the fluid density, $u$ a speed, $L$ a size, $\eta$ the viscosity.
  The number compares inertia with viscous drag. For a crawling cell it is
  about $10^{-8}$ ({prf:ref}`ex-reynolds`), and a bacterium that stops
  swimming coasts about a tenth of an ångström.

Each of these models has few symbols. A model earns its keep by being small
enough to reason about, and every one of them is a deliberate simplification
with a domain of validity.

The next page, [](why-inertia-is-irrelevant.md), sets up the first such
simplification, dropping inertia from the motion of a cell, by computing how
far a cell coasts; [](../overdamped/index.md) makes it and computes exactly
what it costs.
