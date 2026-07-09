+++
order = 8
subject = "Mathematics"
tags = ["math", "numerical-methods", "optimization", "newton", "bfgs", "gradient-descent", "trust-region"]
+++

# Numerical Methods — Nonlinear Systems and Optimization

## 8.1 The Nonlinear Landscape

Q: How does optimization relate to nonlinear root finding?
A: Unconstrained optimization of $f: \mathbb{R}^n \to \mathbb{R}$ reduces to finding stationary points: $\nabla f(\mathbf{x}^*) = \mathbf{0}$ — an $n$-dimensional root-finding problem. Conversely, root finding can be recast as optimization by minimizing $\|\mathbf{F}(\mathbf{x})\|^2$. The two problems share algorithms and analysis; optimization adds the descent requirement (decreasing $f$ each iteration).

Q: Why is local optimization fundamentally easier than global optimization?
A: Because local methods need only decrease $f$ each iteration to make progress, and convergence analysis is Taylor-local. Global optimization must rule out every other minimum — generally NP-hard without structural assumptions (convexity). Most "optimization software" is local; genuine global methods (branch-and-bound, genetic algorithms, simulated annealing) handle smaller problems.

## 8.2 Descent Directions

C: A [descent direction] at $\mathbf{x}_k$ is any vector $\mathbf{p}_k$ with $\nabla f(\mathbf{x}_k)^T \mathbf{p}_k < 0$ — moving in direction $\mathbf{p}_k$ initially decreases $f$.

Q: Why does the negative gradient point in the steepest descent direction?
A: Because the directional derivative $\nabla f^T \mathbf{p}$ is minimized over unit $\mathbf{p}$ at $\mathbf{p} = -\nabla f / \|\nabla f\|$ (Cauchy–Schwarz). So $-\nabla f$ is the direction of fastest decrease — but "steepest" only locally and only in the Euclidean metric. Different metrics (Newton's method's Hessian metric) give different steepest directions.

## 8.3 Gradient Descent

Q: What is the iteration of [gradient descent]?
A: $\mathbf{x}_{k+1} = \mathbf{x}_k - \alpha_k \nabla f(\mathbf{x}_k)$ for step size $\alpha_k > 0$. Simple, only requires gradients, scales to huge problems (stochastic variants power deep learning). Convergence rate: linear, with rate depending on the ratio $\kappa = L/\mu$ of smoothness to strong-convexity constants — slow for ill-conditioned problems.

Q: Why does gradient descent zigzag on elongated valleys?
A: Because the gradient points perpendicular to level curves, which for an elongated valley means almost perpendicular to the valley axis. Each step traverses the valley's width rather than progressing along its length. Classic example: $f(x, y) = x^2 + 100 y^2$ — gradient descent oscillates across the $y$-axis, making $O(\kappa \log(1/\varepsilon))$ steps vs. Newton's $O(\log \log(1/\varepsilon))$.

Q: What does [momentum] add to gradient descent?
A: Accumulates past gradient directions: $\mathbf{v}_{k+1} = \beta \mathbf{v}_k - \alpha \nabla f(\mathbf{x}_k)$, $\mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{v}_{k+1}$. Smooths zigzagging; accelerates progress along consistent directions. [Nesterov's accelerated gradient] achieves provably optimal $O(1/k^2)$ convergence on smooth convex functions — better than plain gradient descent's $O(1/k)$.

## 8.4 Line Search

Q: Why is [line search] essential for robust optimization?
A: Because even a correct descent direction can overshoot: a too-large step produces $f(\mathbf{x}_{k+1}) > f(\mathbf{x}_k)$, ruining monotone descent. Line search chooses $\alpha_k$ to ensure sufficient decrease. Prevents catastrophic steps in Newton's method when Hessians are poorly conditioned or non-positive-definite.

Q: What are the [Armijo / Wolfe conditions]?
A: [Armijo]: $f(\mathbf{x}_k + \alpha \mathbf{p}_k) \leq f(\mathbf{x}_k) + c_1 \alpha \nabla f_k^T \mathbf{p}_k$ with $c_1 \in (0, 1)$ — "sufficient decrease." [Wolfe] adds curvature: $\nabla f_{k+1}^T \mathbf{p}_k \geq c_2 \nabla f_k^T \mathbf{p}_k$ with $c_2 \in (c_1, 1)$ — "step not too short." Standard criteria ensuring convergence of quasi-Newton methods.

Q: What is [backtracking line search]?
A: Start with a generous initial $\alpha$ (often $\alpha = 1$ for Newton steps). While the Armijo condition fails, shrink $\alpha \leftarrow \rho \alpha$ (typical $\rho = 0.5$). Simple, no derivative-based root finding required, guaranteed to terminate for descent directions and bounded $f$ below. The default line-search algorithm in most libraries.

## 8.5 Newton's Method for Optimization

C: [Newton's method] for minimization iterates $\mathbf{x}_{k+1} = \mathbf{x}_k - H(\mathbf{x}_k)^{-1} \nabla f(\mathbf{x}_k)$ where $H$ is the [Hessian] — giving quadratic convergence near a minimum with positive-definite Hessian.

Q: Why does Newton's method converge quadratically near a strict local minimum?
A: Because Taylor expansion of $\nabla f$ around $\mathbf{x}^*$ gives $\nabla f(\mathbf{x}) \approx H^* (\mathbf{x} - \mathbf{x}^*)$. The Newton step $-H^{-1}\nabla f$ cancels the linear part exactly, leaving error of order $\|\mathbf{x}_k - \mathbf{x}^*\|^2$. Same quadratic rate as scalar Newton, requires invertible $H^*$.

Q: Why can Newton's method fail far from a minimum?
A: Because far from $\mathbf{x}^*$ the Hessian may not be positive-definite: the Newton direction $-H^{-1} \nabla f$ could be an ascent direction (away from a minimum, toward a saddle or maximum). Safeguards: modify $H$ to $H + \lambda I$ until positive-definite (Levenberg–Marquardt), use trust regions, or fall back to gradient descent when the Newton step fails a descent check.

Q: What is the computational drawback of Newton's method in high dimensions?
A: Hessian assembly costs $O(n^2)$ evaluations or $n^2$ memory; solving $H \mathbf{p} = -\mathbf{g}$ costs $O(n^3)$. For $n = 10^6$ (deep learning), explicit Hessians are infeasible. Hence quasi-Newton methods that approximate $H$ with low-rank updates or Hessian-free methods that use $H \mathbf{v}$ products instead of forming $H$.

## 8.6 Quasi-Newton Methods

Q: What do [quasi-Newton methods] approximate?
A: They maintain an approximation $B_k \approx H_k$ (or $H^{-1}_k$) updated each iteration using gradient differences: the [secant condition] $B_{k+1}(\mathbf{x}_{k+1} - \mathbf{x}_k) = \nabla f_{k+1} - \nabla f_k$ enforces correct behavior along the last step. Avoid explicit Hessians while retaining superlinear convergence.

Q: What is the [BFGS] update?
A: [Broyden–Fletcher–Goldfarb–Shanno]: rank-two update to $B_k$ that preserves symmetry and positive-definiteness while satisfying the secant condition. $B_{k+1} = B_k - \frac{B_k \mathbf{s}_k \mathbf{s}_k^T B_k}{\mathbf{s}_k^T B_k \mathbf{s}_k} + \frac{\mathbf{y}_k \mathbf{y}_k^T}{\mathbf{y}_k^T \mathbf{s}_k}$ where $\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k$, $\mathbf{y}_k = \nabla f_{k+1} - \nabla f_k$. Default method for smooth unconstrained optimization — superlinear convergence, robust.

Q: Why is [L-BFGS] preferred over BFGS for very large problems?
A: Because BFGS stores a dense $n \times n$ matrix $B_k$ — infeasible for $n > 10^5$. [Limited-memory BFGS] stores only the last $m$ (typically $m = 5$–$20$) pairs $(\mathbf{s}_k, \mathbf{y}_k)$ and computes $B_k^{-1} \nabla f$ via a two-loop recursion. Memory $O(mn)$, per-step cost $O(mn)$ — scales to $n = 10^9$. Default in SciPy and machine-learning second-order methods.

## 8.7 Trust Region Methods

Q: How does a [trust region] method differ from line search?
A: Instead of choosing a direction then a step size along it, trust-region methods constrain the step to a ball $\|\mathbf{p}\| \leq \Delta_k$ around $\mathbf{x}_k$, then choose the best direction within that region. Minimize a quadratic model of $f$ subject to the trust-region constraint. Better for highly nonlinear or non-convex problems — the ball size adapts to observed agreement between model and function.

Q: How is the trust-region radius $\Delta_k$ adapted?
A: Measure the [agreement ratio] $\rho_k = (\text{actual decrease})/(\text{predicted decrease})$. If $\rho_k$ is large (good agreement), expand $\Delta_{k+1}$. If $\rho_k$ is small or negative (model unreliable), shrink $\Delta_{k+1}$ and reject the step. Self-tuning — no hyperparameters to set.

## 8.8 Nonlinear Least Squares

Q: What is the [nonlinear least squares] problem?
A: Minimize $f(\mathbf{x}) = \frac{1}{2} \sum_i r_i(\mathbf{x})^2 = \frac{1}{2}\|\mathbf{r}(\mathbf{x})\|^2$ for residuals $\mathbf{r}: \mathbb{R}^n \to \mathbb{R}^m$. Pervasive in model fitting, curve fitting, computer vision (bundle adjustment), inverse problems. Special structure: $\nabla f = J^T \mathbf{r}$ and $H = J^T J + \sum r_i \nabla^2 r_i$ where $J$ is the Jacobian of $\mathbf{r}$.

Q: What is the [Gauss–Newton] approximation?
A: Drop the second-derivative terms in $H$, keeping only $J^T J$. Iteration: $\mathbf{x}_{k+1} = \mathbf{x}_k - (J^T J)^{-1} J^T \mathbf{r}$. Equivalent to a linear least-squares problem at each step. Fast (superlinear) when residuals are small at the solution; can fail when $J^T J$ is singular or residuals are large.

Q: What does [Levenberg–Marquardt] add to Gauss–Newton?
A: Replace $J^T J$ with $J^T J + \lambda I$ for damping $\lambda$: interpolates between Gauss–Newton ($\lambda = 0$) and gradient descent ($\lambda \to \infty$). Adaptive $\lambda$ behaves like a trust region — shrinks toward Newton near the solution, reverts to gradient descent far away. Workhorse for nonlinear curve fitting (`scipy.optimize.least_squares`).

## 8.9 Constrained Optimization

Q: What are the [KKT conditions] for constrained optimization?
A: [Karush–Kuhn–Tucker] conditions generalize "$\nabla f = 0$" to problems with constraints $g_i(\mathbf{x}) \leq 0$, $h_j(\mathbf{x}) = 0$: (i) stationarity of the Lagrangian $\nabla f = \sum \mu_i \nabla g_i + \sum \lambda_j \nabla h_j$, (ii) primal feasibility, (iii) dual feasibility $\mu_i \geq 0$, (iv) complementary slackness $\mu_i g_i = 0$. Necessary conditions for a local minimum under regularity.

Q: What is the [Lagrangian] and why use dual variables?
A: $\mathcal{L}(\mathbf{x}, \mathbf{\lambda}, \mathbf{\mu}) = f(\mathbf{x}) + \sum \mu_i g_i(\mathbf{x}) + \sum \lambda_j h_j(\mathbf{x})$. Introduces [multipliers] (dual variables) that encode constraint "prices." The KKT conditions say we seek a saddle point: minimum in $\mathbf{x}$, maximum in $\mathbf{\lambda}, \mathbf{\mu}$. Foundation for duality theory and all modern constrained optimization algorithms.

## 8.10 Penalty and Barrier Methods

Q: How do [penalty methods] handle constraints?
A: Replace the constrained problem by an unconstrained one: minimize $f(\mathbf{x}) + \rho \sum \max(0, g_i(\mathbf{x}))^2$ with $\rho \to \infty$. Violating constraints becomes very expensive. Simple but ill-conditioned for large $\rho$ — the unconstrained problem's Hessian blows up near the feasible boundary.

Q: How do [interior-point methods] differ from penalty methods?
A: Use a [barrier] $-\sum \log(-g_i(\mathbf{x}))$ that is finite inside the feasible region and $+\infty$ at the boundary. Minimize $f - (1/t) \sum \log(-g_i)$ for increasing $t$. Stays strictly feasible — avoids the ill-conditioning of penalty methods. Polynomial-time for LP and convex QP (Karmarkar 1984, Nesterov–Nemirovski 1988). Default in `scipy.optimize.linprog(method='interior-point')`, CVXPY.

## 8.11 Convex Optimization

Q: Why is [convex optimization] tractable in a way general nonlinear optimization is not?
A: Because for convex $f$ and convex feasible set: every local minimum IS a global minimum (convexity eliminates local-vs-global distinction). Interior-point methods solve LP, QP, SOCP, SDP in polynomial time. Duality gives provable optimality certificates. The dividing line between "optimization we can solve" and "NP-hard."

Q: How do you recognize a convex optimization problem?
A: $f$ convex AND feasible set convex. Sufficient conditions: $f$ has positive-semidefinite Hessian everywhere; constraints are affine equalities plus convex inequalities. LP, QP with PSD Hessian, SOCP, SDP, geometric programs after log-transform. Convex optimization theory unifies them all; solvers (CVXPY, Mosek, Gurobi) handle the catalog transparently.

## 8.12 Stochastic Optimization

Q: What is [stochastic gradient descent] (SGD) and why does it dominate machine learning?
A: Replace the true gradient $\nabla f = \frac{1}{N}\sum \nabla f_i$ (over $N$ data points) by an estimate using a random mini-batch. Per-iteration cost independent of $N$. Converges in expectation to stationary points at $O(1/\sqrt{k})$ rate. Combined with momentum (Adam, RMSProp), SGD trains billion-parameter neural networks — where any full-batch second-order method is infeasible.

Q: Why doesn't SGD converge to the minimum even after many iterations?
A: Because stochastic gradient noise doesn't vanish — even at the minimum, random mini-batches produce a nonzero gradient estimate, so SGD iterates orbit within a noise ball. Fix: decrease learning rate over time ($\alpha_k \sim 1/\sqrt{k}$) to shrink the noise ball, or use variance-reduction methods (SVRG, SAGA) that converge to the true minimum exactly.
