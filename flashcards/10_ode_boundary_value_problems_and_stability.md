+++
order = 10
subject = "Mathematics"
tags = ["math", "numerical-methods", "ode", "bvp", "shooting", "stability", "a-stability", "stiff"]
+++

# Numerical Methods — ODE Boundary Value Problems and Stability Theory

## 10.1 The Boundary Value Problem

Q: How does a [boundary value problem] (BVP) differ from an initial value problem?
A: Conditions are specified at multiple times (usually both endpoints of the interval) rather than all at one initial time. Example: $y''(x) = f(x, y, y')$ with $y(a) = \alpha$, $y(b) = \beta$. Arises naturally in elliptic PDEs reduced to 1D, in steady-state problems (no time evolution), and in optimal control via Pontryagin's principle.

Q: Why is a BVP fundamentally harder than an IVP?
A: Because no "direction of integration" marches from one fully-specified state. The solution at every interior point depends on BOTH boundary conditions simultaneously — you cannot step forward using only local information. Either guess missing initial data (shooting) or solve globally (finite differences, collocation, finite elements).

## 10.2 The Shooting Method

Q: What is the [shooting method] for solving BVPs?
A: Reformulate as an IVP: guess the missing initial conditions (e.g., $y'(a)$), integrate forward using any IVP solver, check whether the terminal value matches $y(b) = \beta$. Iterate on the guess using a root finder (Newton, secant). "Shoot" repeatedly until you hit the target. Simple and reuses high-quality IVP machinery.

Q: When does the shooting method fail?
A: When the IVP is extremely sensitive to initial conditions — small guess errors explode over the integration interval. Linear BVPs with growing/decaying modes are classic examples: a tiny error in the initial slope becomes enormous at the far boundary. Fix: [multiple shooting] splits $[a, b]$ into subintervals, integrates each separately, and enforces continuity + boundary conditions as a coupled nonlinear system.

## 10.3 Finite Difference Discretization

Q: How do finite differences turn a BVP into a linear system?
A: Discretize $[a, b]$ into $N$ nodes $x_i$; approximate $y''(x_i)$ by $(y_{i+1} - 2y_i + y_{i-1})/h^2$; substitute into the ODE at each interior node. Boundary conditions fix $y_0, y_N$. Result: tridiagonal linear system (or nonlinear for nonlinear BVPs, solved by Newton). Cost: $O(N)$ via Thomas algorithm; error $O(h^2)$.

Q: What's the tradeoff between shooting and finite-difference BVP methods?
A: [Shooting]: simple, reuses adaptive IVP solvers, poor conditioning for sensitive BVPs. [Finite difference]: handles sensitive problems robustly (global constraint), trivially parallelizable, higher upfront cost (solve linear system per iteration). FD is the default for production BVP software (`scipy.integrate.solve_bvp`, `bvp4c`).

## 10.4 Collocation Methods

Q: What is a [collocation] method for BVPs?
A: Represent the solution as a polynomial (or piecewise polynomial / spline) with unknown coefficients. Enforce the ODE exactly at a set of [collocation points]. Boundary conditions fix remaining coefficients. `scipy.integrate.solve_bvp` uses collocation on cubic splines with adaptive mesh refinement. Naturally gives a continuous, differentiable solution (not just at grid points).

## 10.5 Linear Stability — The Test Equation

C: The [test equation] $y' = \lambda y$ with $y(0) = 1$ has true solution $y(t) = e^{\lambda t}$ — decaying iff $\text{Re}(\lambda) < 0$. Numerical stability regions are analyzed by applying each method to this equation.

Q: What is the [stability region] of a numerical ODE method?
A: The set of complex $z = h\lambda$ for which the method's numerical solution to $y' = \lambda y$ does not grow. Applied to $y' = \lambda y$, each method produces $y_{n+1} = R(h\lambda) y_n$ for some [amplification factor] $R$. Stability region: $\{z : |R(z)| \leq 1\}$. Method is stable on this region, unstable outside.

Q: What is the amplification factor for forward Euler?
A: $R(z) = 1 + z$. Stability region: $|1 + z| \leq 1$, a disk of radius $1$ centered at $-1$ in the complex plane. For $\lambda < 0$ real: requires $h |\lambda| \leq 2$ — the origin of the "$h$ must be small enough" constraint.

Q: What is the amplification factor for backward Euler?
A: $R(z) = 1/(1 - z)$. Stability region: $|1 - z| \geq 1$ — the entire complex plane EXCEPT a disk of radius $1$ around $z = 1$. Includes the whole left half-plane: any $h > 0$ is stable for $\text{Re}(\lambda) < 0$. This is why backward Euler can take enormous steps for decaying modes.

## 10.6 A-Stability

C: A method is [A-stable] if its stability region contains the entire left half-plane $\{z : \text{Re}(z) \leq 0\}$ — meaning no step-size restriction for decaying modes $\text{Re}(\lambda) < 0$.

Q: What is the [Dahlquist second barrier]?
A: No explicit linear multistep method is A-stable. Moreover, no A-stable linear multistep method has order greater than $2$. Trapezoidal rule is the optimal A-stable multistep method — achieves order $2$ at the edge of the barrier. To exceed order $2$ with A-stability, must use Runge–Kutta or non-standard methods.

Q: Why is A-stability alone not sufficient for stiff problems?
A: Because A-stability permits slow damping of fast transients (e.g., trapezoidal rule has $|R(z)| = 1$ on the imaginary axis — no damping). A stiff system's fast modes should die FAST numerically, not oscillate. Hence stronger notions: [L-stability] requires $R(z) \to 0$ as $\text{Re}(z) \to -\infty$ — fast modes decay to zero in one step.

Q: Which methods are L-stable?
A: Backward Euler ($R(z) = 1/(1-z) \to 0$ as $z \to -\infty$) and its higher-order generalizations (BDF, Radau IIA). Trapezoidal rule is A-stable but NOT L-stable — famous for producing oscillating tails on problems with large negative eigenvalues. For highly stiff systems, prefer L-stable methods.

## 10.7 Absolute vs. Zero Stability

Q: What's the difference between [zero-stability] and [absolute stability]?
A: [Zero-stability]: behavior as $h \to 0$ — perturbations remain bounded in the limit. Required for convergence (Dahlquist equivalence). [Absolute stability]: behavior at FIXED nonzero $h$ for the test equation — whether $|R(h\lambda)| \leq 1$. Different concerns: zero-stability is about "can the method converge at all?"; absolute stability is about "what $h$ can I actually use?"

## 10.8 Stability for Systems

Q: How does stability analysis extend from scalar to systems $\mathbf{y}' = A\mathbf{y}$?
A: Diagonalize $A = V\Lambda V^{-1}$; in the eigenbasis, the system decouples into scalar test equations $z_i' = \lambda_i z_i$. The overall method is stable iff each eigenvalue's product $h\lambda_i$ lies in the stability region. Stiffness = large spread of $|\lambda_i|$ in the left half-plane; stability dictated by the most extreme $\lambda_i$.

## 10.9 Stiff Solvers in Depth

Q: Why are stiff ODE solvers dramatically more expensive per step than explicit ones?
A: Because each step requires solving a nonlinear system $\mathbf{G}(\mathbf{y}_{n+1}) = \mathbf{0}$ via Newton iteration — each Newton step requires evaluating (or approximating) the Jacobian $\partial \mathbf{f}/\partial \mathbf{y}$ and solving a linear system. Per-step cost: $O(n^3)$ for dense Jacobian vs. $O(n)$ for explicit RK. Trade expensive steps for far fewer of them.

Q: How do production stiff solvers minimize Jacobian cost?
A: (i) Reuse the same Jacobian for many Newton iterations ("modified Newton"). (ii) Factor the Jacobian once, use the factorization for multiple steps. (iii) Approximate the Jacobian by finite differences. (iv) Exploit sparsity: many physical systems have banded or sparse Jacobians — pass sparsity pattern to the solver for orders-of-magnitude speedup.

## 10.10 LSODA and Automatic Stiffness Detection

Q: What does [LSODA] (Livermore Solver for ODEs, Automatic) do?
A: Automatically switches between Adams–Bashforth–Moulton (non-stiff) and BDF (stiff) based on runtime diagnostics of stiffness. Starts with Adams methods; if performance suggests stiffness (small step sizes despite smooth solution), switches to BDF. User need not know the problem type in advance. Standard in many legacy codes; `scipy.integrate.odeint` wraps LSODA.

## 10.11 Convergence Analysis

Q: What's the relation between [local truncation error] and global error for a $p$-th order method?
A: Local truncation error per step is $O(h^{p+1})$; accumulated over $O(1/h)$ steps on a fixed interval, it becomes $O(h^p)$ global. This is WHY "order" is measured by global error. Example: RK4 has local error $O(h^5)$, global error $O(h^4)$. Valid provided the method is zero-stable and the ODE is Lipschitz.

## 10.12 Practical Guidelines

Q: When do you use shooting vs. global methods for a BVP?
A: [Shooting]: the BVP is well-conditioned (initial-value sensitivity is modest), or you already have a great IVP solver. Linear BVPs with non-growing modes. [Global (FD / collocation)]: stiff or sensitive BVPs, higher-dimensional BVPs, when you want guaranteed convergence. `solve_bvp` is the default in SciPy; `bvp4c` in Matlab.

Q: How do you detect stiffness empirically if unsure?
A: Run an explicit solver (RK45) with a loose tolerance. Symptoms of stiffness: (i) step size shrinks to tiny values despite smooth-looking solution, (ii) solver takes $10^5+$ steps for a short interval, (iii) solution barely changes yet the solver crawls. Switch to BDF or Radau; expect dramatic speedup if the problem was indeed stiff.
