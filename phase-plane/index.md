---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/phase-plane.pdf
    id: phase-plane-index-pdf
downloads:
  - id: phase-plane-index-pdf
    title: Download PDF
---

# The Phase Plane

:::{tip} Big Idea
Two variables, one plane. Near a steady state the plane is linear, and
eigenvectors, the trace and the determinant say what it looks like. Keeping
the drug that the saddle-node model eliminated loses nothing: the plane
reproduces the phase line and shows why the fast/slow reduction was allowed.
Fast/slow reduction is then made into a general method with a recipe and an
error estimate. Flipping one sign, so that an activator makes its own
inhibitor, turns the folds of an S-shaped switch into oscillations and
pulses: a trapping region and the Poincaré–Bendixson theorem give a closed
orbit.
:::

## How the two models compare

| | drug and cells | FitzHugh–Nagumo |
|---|---|---|
| cross effects | same sign: mutual harm | opposite signs: activator and inhibitor |
| fast nullcline | a curve without folds | an S with two folds |
| pull toward the fast nullcline | $-\delta(1 + N/N_c)$, never zero | $(1 - v^2)/\ve$, zero at the folds |
| what can happen | switch: nodes, saddle, separatrix | oscillation, excitable pulses |
| threshold | separatrix of a saddle; permanent outcome | middle branch of the S; transient pulse |
