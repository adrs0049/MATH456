# MATH 456 Notes

## About these notes

These are notes for MATH 456, Mathematical Modeling, given by
[Andreas Buttenschoen](https://www.buttenschoen.ca) at the University of
Massachusetts Amherst. They introduce nonlinear dynamics through models from
cell and population biology. The emphasis is on building a model, reducing it
to something simple enough to analyze, and checking what the reduction cost.

## Contents

0. **[What Is a Model](modeling/index.md).** What a mathematical model is, and
   how estimating the size of each term tells us which terms to keep.

1. **[From Newton's Second Law to a Force Balance](overdamped/index.md).** A
   cell's equation of motion, its scaling, and the singular limit in which
   inertia is dropped: two time scales, inner and outer solutions, and the
   error of the reduction.

2. **[Flows on the Line](one-dimensional/index.md).** One-dimensional
   nonlinear ODEs: balance equations, the phase line, steady states and their
   stability, and the transcritical and saddle-node bifurcations, with
   examples from antibiotic resistance.

## Background

Calculus and an introductory course in ordinary differential equations. No
biology is assumed: the biology each model needs is explained where it is
used.

## Computation

The figures are produced by Python code (`numpy`, `scipy`, `matplotlib`) that
is included on each page and can be expanded.

## License

This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>
