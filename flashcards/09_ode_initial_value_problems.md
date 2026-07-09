+++
order = 9
subject = "Mathematics"
tags = ["math", "numerical-methods", "ode", "runge-kutta", "euler", "multistep", "adaptive"]
+++

# Numerical Methods — ODE Initial Value Problems

## 9.1 The IVP Problem

Q: What is an [initial value problem] (IVP)?
A: Find $\mathbf{y}: [t_0, t_f] \to \mathbb{R}^n$ satisfying $\mathbf{y}'(t) = \mathbf{f}(t, \mathbf{y}(t))$ with $\mathbf{y}(t_0) = \mathbf{y}_0$. "Given the state now and its derivative law, evolve it forward in time." Pervasive: population dynamics, celestial mechanics, chemical kinetics, epidemiology, neural network training dynamics.

Q: Why do practical ODEs usually lack closed-form solutions?
A: Because closed forms exist only for a handful of special classes (linear constant-coefficient, separable, a few integrable systems). Generic nonlinear ODEs — Lorenz, Van der Pol, $N$-body gravity for $N \geq 3$ — are provably chaotic or have no elementary primitives. Numerical integration is the universal tool.

## 9.2 Euler's Method

C: [Forward Euler] approximates $\mathbf{y}'= \mathbf{f}(t, \mathbf{y})$ by $\mathbf{y}_{n+1} = \mathbf{y}_n + h \mathbf{f}(t_n, \mathbf{y}_n)$ — take the current slope and step forward by $h$.

Q: What is the [local truncation error] of forward Euler?
A: $O(h^2)$ per step — from the Taylor series $\mathbf{y}(t + h) = \mathbf{y}(t) + h\mathbf{y}'(t) + \frac{h^2}{2}\mathbf{y}''(t) + \dots$, truncated at first order.

Q: What is the [global error] of forward Euler over a fixed interval?
A: $O(h)$ — the local $O(h^2)$ errors accumulate over $1/h$ steps. Forward Euler is "first-order accurate": halving $h$ halves the error.

Q: What is [backward Euler] and why is it useful for stiff problems?
A: $\mathbf{y}_{n+1} = \mathbf{y}_n + h \mathbf{f}(t_{n+1}, \mathbf{y}_{n+1})$ — the derivative is evaluated at the new (unknown) time. Implicit: each step solves a nonlinear system for $\mathbf{y}_{n+1}$. Unconditionally stable (A-stable), so enormous time steps are permitted for stiff systems where forward Euler would explode.

## 9.3 Consistency, Stability, Convergence

Q: State the [Dahlquist equivalence theorem] (also called Lax–Richtmyer for ODEs).
A: For linear multistep methods: CONSISTENCY (local truncation error $\to 0$ as $h \to 0$) + [zero-stability] (bounded response to perturbations) $\iff$ CONVERGENCE (numerical solution $\to$ true solution as $h \to 0$). The unifying principle: an order-$p$ consistent zero-stable method has global error $O(h^p)$.

## 9.4 Runge–Kutta Methods

C: [Runge–Kutta] methods approximate $\mathbf{y}(t + h)$ by evaluating $\mathbf{f}$ at multiple "stages" within $\lbrack t, t + h\rbrack $ and combining linearly — no memory of past steps needed.

Q: Why does [classical RK4] achieve fourth-order accuracy with four $\mathbf{f}$-evaluations per step?
A: Because the coefficients $k_i$ and weights are chosen so that the Taylor expansion of $\mathbf{y}_{n+1}$ matches $\mathbf{y}(t_{n+1})$ through $O(h^4)$ — the first five terms. Result: $\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$ where $k_i$ are slopes at sampled intermediate points. Global error $O(h^4)$.

Q: What is the general form of an [explicit Runge–Kutta] method?
A: $k_i = \mathbf{f}(t_n + c_i h, \mathbf{y}_n + h \sum_{j < i} a_{ij} k_j)$ for $i = 1, \dots, s$; $\mathbf{y}_{n+1} = \mathbf{y}_n + h \sum_i b_i k_i$. Coefficients $(c_i, a_{ij}, b_i)$ arranged in a [Butcher tableau]. Explicit means $a_{ij} = 0$ for $j \geq i$ — each stage uses only earlier stages, no equation solving.

Q: What do [implicit Runge–Kutta] methods trade off?
A: Allow $a_{ij} \neq 0$ for $j \geq i$ — stages depend on each other, requiring a nonlinear solve per step. Cost: expensive each step (typically Newton iteration with Jacobian). Reward: A-stability achievable (impossible for explicit RK), letting implicit RK handle stiff problems with large steps. Gauss–Legendre RK methods achieve order $2s$ with $s$ stages.

## 9.5 Multistep Methods

Q: How do [multistep methods] differ from Runge–Kutta?
A: Use HISTORY: $\mathbf{y}_{n+1}$ is computed from several previous values $\mathbf{y}_n, \mathbf{y}_{n-1}, \dots$ and their derivatives — one new $\mathbf{f}$ evaluation per step (plus storage). More memory-efficient per accuracy order than RK. Drawback: need a startup phase (first few steps via RK) and restarts are painful for adaptive step size.

Q: What are [Adams–Bashforth] methods?
A: Explicit multistep: $\mathbf{y}_{n+1} = \mathbf{y}_n + h \sum_{j=0}^{k-1} \beta_j \mathbf{f}_{n-j}$ where $\beta_j$ are chosen so the method exactly integrates polynomials through the past $\mathbf{f}$-values. Adams–Bashforth of order $k$ achieves $O(h^k)$ global error using $k$ past $\mathbf{f}$ evaluations and ONE new one. Cheap per step.

Q: What are [Adams–Moulton] methods?
A: Implicit multistep counterpart: $\mathbf{y}_{n+1} = \mathbf{y}_n + h \sum_{j=-1}^{k-1} \beta_j \mathbf{f}_{n-j}$ — includes $\mathbf{f}_{n+1}$ on the right-hand side. Achieves order $k+1$ with same memory, better stability than Adams–Bashforth. Often paired with AB as a [predictor–corrector]: AB predicts $\mathbf{y}_{n+1}$, AM corrects.

## 9.6 Adaptive Step Size

Q: Why is adaptive step sizing nearly essential for practical ODE integration?
A: Because stiffness, rapid transients, and slow quasi-steady regimes coexist in most real problems. Fixed-step integration chooses either a step size too small for efficiency or too large for accuracy. Adaptive methods estimate local error each step and expand or shrink $h$ to maintain user-specified tolerance — orders of magnitude fewer steps for same accuracy.

Q: How do [embedded Runge–Kutta] methods estimate local error cheaply?
A: Use two RK methods sharing stages but with different weights $b_i, \tilde{b}_i$ giving orders $p$ and $p + 1$. The difference $\|\mathbf{y}_{n+1} - \tilde{\mathbf{y}}_{n+1}\|$ estimates the local error at (essentially) no extra cost beyond one RK step. Classic examples: [RK45 / Dormand–Prince] (Matlab's `ode45`, SciPy's `RK45`), [RK23] for lower orders.

Q: What is the standard step-size adjustment formula?
A: Given current step $h$, local error estimate $\mathbf{e}$, tolerance $\tau$, and method order $p$: $h_{\text{new}} = h \cdot \text{safety} \cdot (\tau/\|\mathbf{e}\|)^{1/(p+1)}$ with a safety factor $\sim 0.9$ and bounded shrink/grow ratios. If $\|\mathbf{e}\| > \tau$: reject the step, retry with smaller $h$. Rooted in the error's known dependence on $h^{p+1}$.

## 9.7 Stiffness

C: A [stiff ODE] has vastly different time scales: fast transients (large negative eigenvalues of the Jacobian) coexist with slow dynamics, forcing explicit methods to take tiny steps for stability even when accuracy allows larger ones.

Q: Why do explicit methods fail on stiff systems?
A: Because explicit stability requires $h \cdot |\lambda_{\max}| < C$ where $|\lambda_{\max}|$ is the largest-magnitude eigenvalue of $\partial \mathbf{f}/\partial \mathbf{y}$. For stiff problems with $|\lambda_{\max}| \sim 10^6$, this forces $h \lesssim 10^{-6}$ — perhaps $10^6$ steps to traverse a unit time interval. Accuracy needs much larger $h$; stability forbids it. Use implicit methods.

Q: What is the [stiffness ratio] and how stiff is "stiff"?
A: $|\lambda_{\max}|/|\lambda_{\min}|$ of the Jacobian. Ratios of $10^3$ or more typically call for implicit methods; ratios of $10^6$+ (combustion chemistry, circuit simulation) demand them absolutely. Diagnostic: if an explicit method takes thousands of tiny steps and the solution "looks smooth," suspect stiffness and switch to `BDF` or `Radau`.

## 9.8 BDF Methods

Q: What are [BDF] methods?
A: [Backward Differentiation Formulas]: implicit multistep methods $\sum_{j=0}^{k} \alpha_j \mathbf{y}_{n+1-j} = h \beta_0 \mathbf{f}_{n+1}$ — fit a polynomial through the past $k$ values PLUS the unknown $\mathbf{y}_{n+1}$, match its derivative at $t_{n+1}$ to $\mathbf{f}$. BDF1 is backward Euler; BDF2, BDF3, ..., BDF6 are the workhorses for stiff ODEs. Orders $\geq 7$ are unstable (second Dahlquist barrier).

Q: Why are BDF methods the go-to for stiff systems?
A: Because they are A-stable (BDF1, BDF2) or $A(\alpha)$-stable for higher orders — large stability regions in the left half-plane. Combined with good newton-solver implementations and sparse Jacobian handling, BDF methods dominate production stiff ODE codes: LSODA, DASSL, Sundials CVODE, SciPy's `BDF`.

## 9.9 Symplectic Integrators

Q: Why use [symplectic integrators] for Hamiltonian systems?
A: Because they conserve a modified Hamiltonian EXACTLY (up to roundoff), giving bounded energy error over arbitrarily long times. Standard non-symplectic methods (RK4, BDF) drift in energy over millions of steps — ruining long-time orbital dynamics simulations. Symplectic leapfrog, Verlet, and Yoshida methods are standard in molecular dynamics and celestial mechanics.

Q: What is the [leapfrog method]?
A: For $\ddot{\mathbf{x}} = \mathbf{F}(\mathbf{x})$: $\mathbf{x}_{n+1} = \mathbf{x}_n + h \mathbf{v}_{n+1/2}$, $\mathbf{v}_{n+3/2} = \mathbf{v}_{n+1/2} + h \mathbf{F}(\mathbf{x}_{n+1})$ — velocities staggered half a step from positions. Second-order accurate, symplectic, time-reversible. Foundational method for molecular dynamics (Verlet form) and $N$-body codes.

## 9.10 Differential-Algebraic Equations

Q: What is a [DAE] (differential-algebraic equation)?
A: A mixed system $\mathbf{F}(t, \mathbf{y}, \mathbf{y}') = \mathbf{0}$ that includes algebraic constraints (equations not involving $\mathbf{y}'$). Arises in constrained mechanics ("particle on a sphere"), circuit simulation (Kirchhoff laws + capacitor dynamics), and index reduction. BDF methods and IDA/DASSL handle DAEs directly via implicit formulation.

Q: What is the [index] of a DAE?
A: The number of times algebraic constraints must be differentiated to get a system of ODEs. Index-1 DAEs (most common) behave essentially like ODEs. Index-$\geq 2$ DAEs are notoriously hard — standard solvers often require preprocessing to reduce index. Constrained multibody dynamics is typically index-3.

## 9.11 Sensitivity and Adjoints

Q: How do you compute the sensitivity of an ODE solution to its parameters?
A: [Forward sensitivity]: solve the augmented ODE $(\mathbf{y}', \mathbf{s}'_i = (\partial \mathbf{f}/\partial \mathbf{y}) \mathbf{s}_i + \partial \mathbf{f}/\partial p_i)$ — one extra IVP per parameter. [Adjoint sensitivity]: solve the ODE forward, then solve a backward adjoint ODE — one extra IVP total regardless of parameter count. Crucial for inverse problems, neural ODEs, optimal control.

## 9.12 Choosing an Integrator

Q: How do you pick an ODE solver for a given problem?
A: Non-stiff, moderate accuracy → [RK45] (`scipy.integrate.solve_ivp(method='RK45')`, Matlab `ode45`). Stiff → [BDF] or [Radau] (`method='BDF'` / `'Radau'`, Matlab `ode15s`). Hamiltonian long-time dynamics → [symplectic] (leapfrog, Yoshida). Near-discontinuous RHS → low-order with tight tolerance. Rule: try `RK45` first; if it takes $10^5+$ steps, switch to `BDF`.
