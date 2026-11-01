---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/overdamped-exercises.pdf
    id: overdamped-exercises-pdf
downloads:
  - id: overdamped-exercises-pdf
    title: Download PDF
---

# Exercises

The exercises follow the pages of this part. Most are calculations by hand
with the integrating factor ({prf:ref}`prop-integrating-factor`), and a few ask
for a short computation in Python (or any other language). Each one is a
problem with two time scales: find the small parameter, reduce, and ask what
the reduction gets wrong.

## Nondimensionalization

:::{exercise} A sinking cell
:label: exr-od-sinking

A diatom, a single-celled alga, is denser than seawater and sinks. Treat it as
a sphere of radius $a$ and density $\rho$ in water of density
$\rho_f = \rho - \Delta\rho$ and viscosity $\eta$. Its weight minus the
buoyancy is $\Delta m\,g$ with $\Delta m = \tfrac43\pi a^3\Delta\rho$, the drag
is Stokes drag {eq}`eq-stokes-drag` with $\gamma = 6\pi\eta a$, and

$$
m\,\frac{\dd v}{\dd t} = \Delta m\,g - \gamma v, \qquad v(0) = 0 .
$$

1. Find the sinking speed $v_s$ at which drag balances the buoyant weight, the
   inertial time $\tau_i = m/\gamma$, and the time $T = 2a/v_s$ to sink one
   diameter, in terms of $\rho$, $\Delta\rho$, $a$, $\eta$ and $g$.
2. Scale with $v = v_s\tilde v$ and $t = T\tilde t$, as in
   {prf:ref}`def-epsilon-inertia`, and show that
   $\ve\,\dd\tilde v/\dd\tilde t = 1 - \tilde v$ with
   $\ve = 2\rho\,\Delta\rho\,g\,a^3/(81\eta^2)$.
3. Take $a = 10\ \mu$m, $\rho = 1050$ kg/m³, $\Delta\rho = 50$ kg/m³,
   $\eta = 10^{-3}$ Pa s and $g = 9.81$ m/s². Find $v_s$ in metres per day,
   $\tau_i$ and $\ve$. How long does the cell take to sink through a surface
   layer $100$ m deep?
4. $\ve$ grows like $a^3$. How large must the cell be for $\ve = 0.1$? Show
   that $\ve = (\rho/\rho_f)\,\mathrm{Re}/18$, where
   $\mathrm{Re} = \rho_f v_s\,2a/\eta$ is the Reynolds number. Stokes drag
   holds only for $\mathrm{Re} \ll 1$. What does this say about when the
   overdamped equation can be used?
:::

:::{exercise} A bead in an optical trap
:label: exr-od-trap

An optical trap, a focused laser beam, holds a small bead near its centre
with the restoring force $-kx$, where $x$ is the distance from the centre.
In water,

$$
m\,\ddot x = -k\,x - \gamma\,\dot x .
$$

1. Find three combinations of $m$, $\gamma$ and $k$ with units of time, and
   show that one of them is the geometric mean of the other two.
2. The bead is held at rest at $x_0$, and released at $t = 0$. Scale with
   $x = x_0\tilde x$ and $t = (\gamma/k)\,\tilde t$ and show that

   $$
   \ve\,\frac{\dd^2\tilde x}{\dd\tilde t^2} + \frac{\dd\tilde x}{\dd\tilde t} + \tilde x = 0,
   \qquad \tilde x(0) = 1, \quad \frac{\dd\tilde x}{\dd\tilde t}(0) = 0,
   \qquad \ve = \frac{mk}{\gamma^2} .
   $$

   Write $\ve$ as a ratio of two of the times of part 1.
3. A silica bead has radius $0.5\ \mu$m and density $2000$ kg/m³, and a
   typical trap has $k = 10$ pN/$\mu$m. Find $\tau_i = m/\gamma$, $\gamma/k$
   and $\ve$.
4. Scale time with $\sqrt{m/k}$ instead, the period (over $2\pi$) of the bead
   without drag. Show that the equation becomes
   $\tilde x'' + \ve^{-1/2}\tilde x' + \tilde x = 0$. Which term now has the
   largest prefactor, and what does that say about this choice of scale?
:::

## The singular limit

:::{exercise} Regular or singular?
:label: exr-od-regular-singular

For each problem, decide whether the limit $\ve \to 0$ is regular or singular
({prf:ref}`def-singular-limit`), find the reduced solution, and compare it
with the exact one.

1. $\dd x/\dd t = -x + \ve x^2$, $x(0) = 1$. Show that the exact solution is
   $x = 1/\big(\ve + (1 - \ve)e^{t}\big)$ (hint: $y = 1/x$ solves a linear
   equation). Then substitute $x = x_0 + \ve x_1 + O(\ve^2)$, collect powers
   of $\ve$, solve for $x_0$ and $x_1$ with the integrating factor, and check
   the result against the exact solution.
2. $\ve\,\dd x/\dd t = 1 - x$, $x(0) = 0$.
3. $\ve\,\dd x/\dd t = x - 1$, $x(0) = 1 + \delta$. How long does it take until
   $|x - 1| = 1$? Evaluate this for $\ve = 0.01$ and $\delta = 10^{-6}$. Is
   the reduced solution of any use?
:::

:::{exercise} The trap, reduced
:label: exr-od-trap-limit

Continue with the bead of the exercise *A bead in an optical trap*, dropping
the tildes:

$$
\ve\,\ddot x + \dot x + x = 0, \qquad x(0) = 1, \quad \dot x(0) = 0 .
$$

1. Write the equation as the system $\dot x = u$, $\ve\,\dot u = -x - u$.
   Which variable is fast? Set $\ve = 0$ and solve the outer problem. Which of
   the two initial conditions is lost?
2. Write the inner problem for $u$ on the fast time $\sigma = t/\ve$, with
   $x$ frozen at its initial value, and solve it. Match it with the outer
   velocity and form the composite velocity.
3. Integrate the composite velocity from $x(0) = 1$. Once the layer is over,
   how far is the bead from the outer solution? Compare with the coasting
   distance of {prf:ref}`ex-coasting-distance`.
4. Find the roots of $\ve\lambda^2 + \lambda + 1 = 0$. For which $\ve$ are
   they real? Show that $\lambda_s = -1 - \ve + O(\ve^2)$ and
   $\lambda_f = -1/\ve + 1 + O(\ve)$, and that the exact solution is
   $x = Ae^{\lambda_s t} + Be^{\lambda_f t}$ with $A = 1 + \ve + O(\ve^2)$
   and $B = -\ve + O(\ve^2)$. Which term is the layer? Compare with part 3.
5. What does the condition for real roots require of the trap stiffness $k$,
   for the bead of the previous exercise?
:::

:::{exercise} An oscillating force
:label: exr-od-oscillating

A cell that extends and retracts a protrusion pushes forward and back. Take
the scaled equation of motion {eq}`eq-eom-scaled` with

$$
\ve\,\frac{\dd\tilde v}{\dd\tilde t} = \cos\tilde t - \tilde v,
\qquad \tilde v(0) = \tilde v_0 .
$$

1. Find the outer and inner solutions and the composite.
2. Solve the problem exactly with {prf:ref}`prop-integrating-factor`. Show
   that after the layer the velocity oscillates with amplitude
   $1/\sqrt{1 + \ve^2}$ and lags behind the force by the phase
   $\arctan\ve \approx \ve$. What is this lag in dimensional time?
3. Compute the difference between the exact solution and the composite.
   Compare it with the bound of {prf:ref}`thm-composite-uniform`. Is the bound
   nearly attained?
4. Now let the force oscillate on the fast time scale,
   $\tilde F = \cos(\tilde t/\ve)$. Why does
   {prf:ref}`thm-composite-uniform` say nothing useful? Find the periodic
   solution, and compare its amplitude and phase with the outer solution.
   Which observables does the reduced model still get right?
5. Plot the exact solution and the outer solution for both forces with
   $\ve = 0.1$ and $\tilde v_0 = 0$.
:::

:::{exercise} A protein removed in pairs
:label: exr-od-dimer

A protein is made at the rate $A(t)$ and removed by dimerization: two
molecules bind and the dimer is degraded, so the loss is $bC^2$. Production
is switched on at $t = 0$ and then increases slowly, over a time $T$:

$$
\frac{\dd C}{\dd t} = A_0\Big(1 + \frac{t}{T}\Big) - bC^2, \qquad C(0) = 0 .
$$

1. Give the units of $A_0$ and $b$, and show that $1/\sqrt{A_0 b}$ is a time.
   Scale with $C = \sqrt{A_0/b}\;c$ and $t = T\tau$ to get

   $$
   \ve\,\frac{\dd c}{\dd\tau} = 1 + \tau - c^2,
   \qquad \ve = \frac{1}{T\sqrt{A_0 b}} .
   $$

2. First ignore the ramp: $\ve\,\dd c/\dd\tau = 1 - c^2$. The outer problem
   has two solutions; which one is relevant? On the fast time
   $\sigma = \tau/\ve$, show from the sign of $1 - c^2$ that every
   $c(0) \ge 0$ tends to $1$, and that $c = -1$ repels. Solve for $c(0) = 0$.
3. With the ramp, find the outer and inner solutions and the composite.
4. Solve the full problem numerically for $\ve = 0.1$, $0.05$, $0.025$ and
   $0.0125$, and compute the largest difference from the composite on
   $0 \le \tau \le 3$. How does it scale with $\ve$?
5. Improve the outer solution: substitute $c = c_0 + \ve c_1 + O(\ve^2)$ and
   show that $c_1 = -1/\big(4(1 + \tau)\big)$. Why is the correction
   negative?
6. Item 4 of {prf:ref}`rmk-reading-bound` needs the outer solution to be
   attracting. Check that it is.
:::

## Error analysis

:::{exercise} A slow-release tablet
:label: exr-od-tablet

A slow-release tablet delivers a dose $D$ of drug into the gut, from which it
is absorbed into the blood at rate $k_a$. The blood, of volume $V$, clears it
at rate $k_e$. With $G$ the amount in the gut and $C$ the concentration in the
blood,

$$
\frac{\dd G}{\dd t} = -k_aG, \quad G(0) = D,
\qquad
\frac{\dd C}{\dd t} = \frac{k_aG}{V} - k_eC, \quad C(0) = 0 .
$$

The tablet is designed so that absorption is slow: $k_a \ll k_e$.

1. Scale with $\tau = k_at$, $g = G/D$ and $c = C/C_u$, where
   $C_u = k_aD/(Vk_e)$. Show that

   $$
   \frac{\dd g}{\dd\tau} = -g,
   \qquad
   \ve\,\frac{\dd c}{\dd\tau} = g - c,
   \qquad \ve = \frac{k_a}{k_e} .
   $$

   Compare with {eq}`eq-eom-scaled`: what plays the role of the force?
2. Find the outer and inner solutions and the composite.
3. Solve exactly with {prf:ref}`prop-integrating-factor`, and show that the
   exact solution and the composite differ by exactly $\ve c(\tau)$.
4. Compare the exact and the outer solution in three observables:
   (a) the area under the curve $\int_0^\infty c\,\dd\tau$; show that its
   dimensional value $\int_0^\infty C\,\dd t = D/(Vk_e)$ follows from the
   equations without solving them; (b) the time and height of the peak
   concentration; (c) the rate at which $c$ decays at late times.
5. Take an absorption half-life of $6$ h and an elimination half-life of
   $1$ h. Find $\ve$, the time of the peak in hours and the peak height
   relative to the outer solution. A pharmacologist fits the late decay of
   the blood concentration to measure how fast the drug is eliminated. What
   do they measure?
:::

:::{exercise} Position errors do not accumulate
:label: exr-od-no-accumulation

{prf:ref}`prop-reduction-error` computed the position error of the reduced
model for a constant force. Here the force is arbitrary: $F(t)$ piecewise
continuous, with

$$
m\,\frac{\dd v}{\dd t} = F(t) - \gamma v, \quad v(0) = v_0,
\qquad
\gamma\,\frac{\dd x_\mathrm{out}}{\dd t} = F(t),
$$

and $x$ and $x_\mathrm{out}$ starting from the same $x(0)$.

1. Integrate both equations from $0$ to $t$ and show that, for every $t$,

   $$
   x(t) - x_\mathrm{out}(t) = \tau_i\,\big(v_0 - v(t)\big) .
   $$

2. Suppose $|F(t)| \le \gamma V$ for all $t$ and $|v_0| \le V$. Use
   {eq}`eq-if-motion` to show that $|v(t)| \le V$, and conclude that
   $|x(t) - x_\mathrm{out}(t)| \le 2V\tau_i$ for all $t$.
3. A cell switches its traction on and off every minute for a day. How many
   layers does its velocity go through? Use
   {prf:ref}`prop-reduction-error` to find the sign of the position error
   after a switch on and after a switch off, and explain why the errors do not
   add up.
4. Show that the windowed velocity {eq}`eq-windowed-velocity` of the exact
   and the outer solution differ by at most $2V\tau_i/\Delta t$, for every
   force. Compare with {prf:ref}`prop-windowed-error`.
:::
