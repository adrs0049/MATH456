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
      - newton-with-drag.md
      - scaling-the-equation.md
      - singular-limit.md
      - force-sets-velocity.md
      - integrating-factor.md
downloads:
  - id: overdamped-part-pdf
    title: Download the part (PDF)
  - id: overdamped-index-pdf
    title: Download this page (PDF)
---

# From Newton's Second Law to a Force Balance

:::{tip} Big Idea
Newton's second law is a statement about acceleration, and for a cell it is
the wrong equation. Solving it with the forcing switched off produces a time,
$\tau_i = m/\gamma$, a microsecond. Measuring time in units of the behavior we
want to predict puts a dimensionless number $\ve = \tau_i/T \approx 10^{-9}$ in
front of the inertial term, and setting $\ve = 0$ gives $\gamma\,\dot x = F$.
The limit lowers the order of the equation and discards an initial condition,
and the point of this part is to see what the reduction throws away and to
compute its error in the quantity we actually measure.
:::

## The question

A swimmer who stops stroking glides a body length or two. How far does a cell
glide after its traction is switched off? The answer
({prf:ref}`ex-coasting-distance`) takes four pages, and each uses a technique
that recurs in every later model.

1. [](newton-with-drag.md): the equation of motion with every term present,
   its units, and the two times hidden in its constants.
2. [](scaling-the-equation.md): reference units, the substitution, and the
   single parameter $\ve$ the four constants collapse into.
3. [](singular-limit.md): what setting $\ve = 0$ loses, the two time scales,
   inner and outer solutions, and matching.
4. [](force-sets-velocity.md): the error of the reduction in the observable,
   the overdamped force balance, and what it says about the physics.

The one solution method used throughout is in [](integrating-factor.md).
