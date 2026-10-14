---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/one-dimensional.pdf
    id: one-dimensional-index-pdf
downloads:
  - id: one-dimensional-index-pdf
    title: Download PDF
---

# Flows on the Line

:::{tip} Big Idea
Most models of a single changing quantity are a balance: rate of change equals
rate in minus rate out. For one such equation the graph of the right-hand side
decides the long-time behaviour, without solving anything. Turning a parameter,
such as the dose of an antibiotic, moves the steady states, and at special
values the outcome changes abruptly. The same bookkeeping, applied to gene
copies instead of cells, explains why a gene variant that kills children who
inherit two copies of it is carried by a quarter of the people in parts of
Africa.
:::

## Overview

This part builds the habits of one-dimensional modeling on the simplest
examples: a protein that is made and destroyed, and a population of cells that
divide and die. It then turns the dose of an antibiotic in two bacterial models
and meets the basic bifurcations, uses the same tools on the genetics of sickle
cell anemia, and ends with the bifurcation forced by symmetry. Each page uses
only the ones before it. From Part I it uses only the integrating factor
([](../overdamped/integrating-factor.md)) and the idea of scaling a variable
by its natural unit ([](../overdamped/scaling-the-equation.md)).

| Notes | Topics | Strogatz |
|---|---|---|
| [](rates-in-rates-out.md) | balance equations, units, turnover time, steady state, phase line, slope test | 2.1--2.4 |
| [](populations-and-fractions.md) | per-capita rates, exponential and logistic growth, the fraction of one type, a stable mixture of two strains | 2.3 |
| [](transcritical.md) | bifurcation diagrams, the transcritical bifurcation, critical slowing down | 3.1--3.2 |
| [](saddle-node.md) | the saddle-node bifurcation, bistability, threshold inoculum, hysteresis, linearization and normal forms | 3.1--3.2 |
| [](sickle-cell.md) | gene copies as a population, random mating, a stable mixture of alleles | |
| [](sickle-cell-without-malaria.md) | a steady state the slope test cannot classify, slow algebraic decay | |
| [](pitchfork.md) | symmetry and cell polarity, the pitchfork bifurcation, the imperfect pitchfork | 3.4 |

[](populations-and-fractions.md) ends with a question that
[](sickle-cell.md) answers, and [](transcritical.md) and [](saddle-node.md)
turn the dose in the bacterial models of [](populations-and-fractions.md).

## Modeling conventions

Every model in these notes follows the same conventions.

1. **Write the word equation first.** Rate of change $=$ rate in $-$ rate out.
   Every term is one process that you can name.
2. **Every term has a mechanism.** A loss proportional to the amount means
   that every molecule, or every cell, has the same chance per unit time of
   being lost, independently of the others.
3. **Units must match.** Every term has the units of the left-hand side. This
   catches most errors and tells you what each parameter means.
4. **Rates and times are reciprocals.** A rate $b$ in $1/\text{time}$ comes
   with a time scale $1/b$, the first estimate of how fast anything happens.
5. **Check the model for artifacts.** Densities must not become negative or
   blow up. A process that is known to be regulated needs a stable steady
   state.
6. **State the assumptions in words.** Each assumption should be a sentence a
   biologist would accept or reject.

## Biology you need

None in advance. Each page has short **Biology background** notes that explain
everything its model uses: proteins, cell division, antibiotics, genes,
alleles, inheritance, and malaria, at the level of a first-year biology course.
