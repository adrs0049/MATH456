---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/modeling-index.pdf
    id: modeling-index-pdf
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/modeling.pdf
    id: modeling-part-pdf
    title: "What Is a Model"
    articles:
      - what-is-a-model.md
      - scales-and-forces.md
downloads:
  - id: modeling-part-pdf
    title: Download the part (PDF)
  - id: modeling-index-pdf
    title: Download this page (PDF)
---

# What Is a Model

:::{tip} Big Idea
A model is a deliberate simplification built to answer a specific question,
written precisely enough that the question can be settled. The simplification
is the content of the model, and every simplification has a domain of validity
that can be computed. The first one we make, dropping inertia from the motion
of a crawling cell, is justified by a single ratio of two numbers read off a
table of sizes, speeds, and forces.
:::

This part has two pages.

1. [](what-is-a-model.md): what a mathematical model is, the kinds of model,
   a short history of successful models in biology, and the moves of a modeling
   argument, worked through on the question of why inertia is irrelevant for a
   cell.
2. [](scales-and-forces.md): the orders of magnitude for a crawling cell, the
   Reynolds number, the inertial stopping time, and the forces a cell exerts.
