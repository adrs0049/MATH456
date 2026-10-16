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
values the outcome changes abruptly.
:::

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
everything its model uses: proteins, cell division, and antibiotics, at the
level of a first-year biology course.
