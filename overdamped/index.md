---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/overdamped-index.pdf
    id: overdamped-index-pdf
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/overdamped.pdf
    id: overdamped-part-pdf
    title: "From Newton's Second Law to a Force Balance"
    articles:
      - scaling-the-equation.md
      - singular-limit.md
      - integrating-factor.md
      - error-analysis.md
downloads:
  - id: overdamped-part-pdf
    title: Download the part (PDF)
  - id: overdamped-index-pdf
    title: Download this page (PDF)
---

# From Newton's Second Law to a Force Balance

:::{tip} Big Idea
[](../modeling/why-inertia-is-irrelevant.md) found that a cell's velocity is
$F/\gamma$ except for a few microseconds after each change in force, because
the inertial stopping time $\tau_i = m/\gamma$ is $10^{-9}$ of the time the
cell takes to move its own length. This part turns that observation into the
force balance $\gamma\,\dot x = F$. Measuring time in units of the behavior
we want to predict puts that ratio, $\ve = \tau_i/T$, in front of the
inertial term, and setting $\ve = 0$ gives the force balance. The limit
lowers the order of the equation and discards an initial condition, and the
point of this part is to see what the reduction throws away and to compute its
error in the quantity we actually measure.
:::

## The question

How much error do we make by using the overdamped equation in place of
Newton's law? The sharpest test is a cell whose traction is switched off. A
swimmer who stops stroking glides a body length or two; the overdamped
equation says the cell stops dead. How far does it really glide? This part
answers that in {prf:ref}`ex-coasting-distance`.
