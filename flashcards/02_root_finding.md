+++
order = 2
subject = "mathematics"
tags = ["math", "numerical-methods", "root-finding", "bisection", "newton", "secant", "fixed-point"]
+++

# Numerical Methods — Root Finding

## 2.1 The Root-Finding Problem

Q: What is the general [root-finding problem]?
A: Given a continuous function $f: \mathbb{R} \to \mathbb{R}$ (or higher dimensional), find $x^*$ with $f(x^*) = 0$. Also called "solving $f(x) = 0$" or "finding zeros." Many problems reduce to root finding: optimization (find where $f' = 0$), equation solving (rewrite $g(x) = h(x)$ as $g - h = 0$), implicit time stepping.

Q: Why is root finding usually done numerically rather than symbolically?
A: Because closed-form solutions exist only for special cases (quadratics, cubics, some transcendentals). General polynomials of degree $\geq 5$ have no radical solutions (Abel–Ruffini). Transcendentals ($\cos x = x$, $e^x = x^2$) have no elementary form. Numerical iteration is the universal tool.

## 2.2 Bisection Method

C: The [bisection method] solves $f(x) = 0$ on $\lbrack a, b\rbrack $ where $f(a)f(b) < 0$ by repeatedly halving the interval and keeping the half containing a sign change.

Q: Why is the bisection method guaranteed to converge?
A: By the intermediate value theorem: continuous $f$ with $f(a)f(b) < 0$ has a root in $[a, b]$. Each iteration halves the interval containing the root, so the interval length shrinks as $(b - a)/2^n$. Guaranteed convergence to within any tolerance — no choice of starting guess can prevent it.

Q: What is the [convergence rate] of bisection?
A: [Linear] with rate $1/2$: error halves per iteration. After $n$ iterations, error $\leq (b-a)/2^{n+1}$. To get $10^{-k}$ accuracy, need $\sim 3.3 k$ iterations. Slow but rock-solid — a "safe" method often used as a fallback or to generate good starting guesses for faster methods.

P: Use bisection to find $\sqrt{2}$ to $3$ decimal places via $f(x) = x^2 - 2$ on $[1, 2]$.
S:
**IDENTIFY**: $f(1) = -1, f(2) = 2$ — opposite signs, root exists in $[1, 2]$.

**PLAN**: Iterate until interval width $< 10^{-3}$.

**EXECUTE**:
- $[1, 2]$: mid $1.5$, $f = 0.25 > 0$ → $[1, 1.5]$.
- $[1, 1.5]$: mid $1.25$, $f = -0.4375 < 0$ → $[1.25, 1.5]$.
- $[1.25, 1.5]$: mid $1.375$, $f \approx -0.109 < 0$ → $[1.375, 1.5]$.
- Continue until width $< 10^{-3}$. After $\sim 10$ iterations: $x \approx 1.4142$.

**EVALUATE**: $\sqrt{2} \approx 1.41421$. Bisection finds this slowly but reliably; Newton's method would get $15$ digits in the same number of steps. The reliability-vs-speed tradeoff is central to root finding.

## 2.3 Newton's Method

C: [Newton's method] iterates $x_{n+1} = x_n - f(x_n)/f'(x_n)$ — using the tangent line at $x_n$ to predict the root.

Q: Why does Newton's method converge quadratically near a simple root?
A: Because Taylor expansion: $0 = f(x^*) \approx f(x_n) + f'(x_n)(x^* - x_n) + \frac{1}{2}f''(\xi)(x^* - x_n)^2$. Newton's update sets the linear term to zero, leaving error $\sim (x^* - x_n)^2$. Formally $|x_{n+1} - x^*| \leq C |x_n - x^*|^2$ — the number of correct digits approximately doubles each iteration.

Q: Why can Newton's method diverge from a bad initial guess?
A: Because Newton is only locally convergent: the quadratic convergence guarantee requires starting within a "basin" around the root. Far from the root, the tangent extrapolation can overshoot wildly or point away from the root entirely. Starting in a basin of attraction is essential — hence the use of bisection or bracketing to produce a good guess first.

Q: What happens to Newton's iterate when $f'(x_n) \approx 0$?
A: The update $\Delta x = -f(x_n)/f'(x_n)$ blows up to huge magnitudes, throwing the iterate far from the root. Near flat regions (e.g., a local extremum between guess and root) are therefore danger zones. Safeguarded Newton variants detect this and fall back to bisection for that step.

Q: Why does Newton's method converge only linearly at a multiple root where $f'(x^*) = 0$?
A: Because the quadratic convergence argument requires $f'(x^*) \neq 0$ to cancel the linear Taylor term cleanly. At a multiple root both $f(x_n)$ and $f'(x_n)$ vanish, so their ratio approaches a finite nonzero limit times $(x_n - x^*)$ — giving linear, not quadratic, convergence. Covered in §2.9.

Q: Why do practical root-finders combine Newton with a bracketing method like bisection?
A: Because Newton is fast (quadratic) but unreliable (can diverge); bisection is slow but always converges. A [safeguarded Newton] accepts a Newton step only when it lies inside the current bracket and reduces $|f|$; otherwise it falls back to bisection. Best of both worlds — fast when possible, reliable always.

## 2.4 Secant Method

C: The [secant method] approximates the derivative by a finite difference: $x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}$.

Q: What's the convergence rate of the secant method?
A: [Superlinear] with order $\varphi = (1 + \sqrt{5})/2 \approx 1.618$ (the golden ratio). Slower than Newton's quadratic, but doesn't require derivative evaluation — useful when $f'$ is expensive or unavailable. Per-iteration cost: one function evaluation (vs. two for Newton). Often cheaper in practice.

## 2.5 Fixed-Point Iteration

C: [Fixed-point iteration]: rewrite $f(x) = 0$ as $x = g(x)$, then iterate $x_{n+1} = g(x_n)$. A solution of $x = g(x)$ is a root of $f$.

Q: When does fixed-point iteration converge?
A: The [contraction mapping theorem]: if $|g'(x)| \leq L < 1$ on an interval containing the fixed point, iteration converges with linear rate $L$. If $g'(x^*) = 0$, convergence is quadratic (Newton's method is a fixed-point iteration with this property). Diverges if $|g'(x^*)| > 1$.

P: Solve $x = \cos x$ by fixed-point iteration.
S:
**IDENTIFY**: $g(x) = \cos x$, fixed point in $[0, 1]$.

**PLAN**: Iterate from $x_0 = 0.5$ until convergence.

**EXECUTE**:
1. $x_0 = 0.5$; $g'(x) = -\sin x$; $|g'(0.5)| \approx 0.48 < 1$ ✓.
2. $x_1 = \cos(0.5) \approx 0.8776$.
3. $x_2 = \cos(0.8776) \approx 0.6390$.
4. $x_3 \approx 0.8027$; $x_4 \approx 0.6948$; ... oscillates toward $x^* \approx 0.7391$.

**EVALUATE**: Converges linearly with rate $\sim 0.67$ — takes $\sim 50$ iterations for $10$ digits. Newton on $f(x) = x - \cos x$ with $f'(x) = 1 + \sin x$ converges in $\sim 5$ iterations. Fixed-point iteration is simple but slow.

## 2.6 Brent's Method

Q: What does [Brent's method] combine, and why is it the practitioner's default?
A: It combines [bisection] (always-converges), [secant method] (superlinear), and [inverse quadratic interpolation] (even faster). Uses bisection as a safety net when secant/IQI misbehave. Practical implementations (SciPy's `brentq`, Matlab's `fzero`) default to Brent's method — guaranteed convergence AND superlinear speed. Best of both worlds.

## 2.7 Multidimensional Newton

Q: How does Newton's method generalize to systems $\mathbf{F}(\mathbf{x}) = \mathbf{0}$ in $\mathbb{R}^n$?
A: $\mathbf{x}_{k+1} = \mathbf{x}_k - J_F(\mathbf{x}_k)^{-1} \mathbf{F}(\mathbf{x}_k)$, where $J_F$ is the Jacobian matrix. Each iteration requires solving a linear system $J_F \Delta \mathbf{x} = -\mathbf{F}$ — no explicit inversion. Quadratic convergence is preserved under nonsingular Jacobian and sufficient smoothness.

Q: Why do large-scale nonlinear solvers avoid forming and factoring the Jacobian explicitly?
A: Because $J_F$ has $n^2$ entries and $n^3$ factorization cost. For $n \sim 10^6$ (PDE discretizations), this is prohibitive. Instead: [Newton–Krylov methods] solve the linear system iteratively (GMRES) using only Jacobian-vector products, which can be approximated by finite differences without explicitly forming $J_F$.

## 2.8 Stopping Criteria

Q: What are sensible stopping criteria for iterative root finding?
A: Combination: (i) $|f(x_n)| < \tau_f$ (residual small). (ii) $|x_{n+1} - x_n| < \tau_x$ (step small, relative to $|x_n|$). (iii) Maximum iteration count as a safety cap. Using only one criterion is fragile — $|f(x)|$ can be small far from a root for flat functions; $|\Delta x|$ can be small due to stagnation, not convergence.

## 2.9 Multiplicity and Degraded Convergence

Q: Why does Newton's method converge only linearly at a root of multiplicity $m > 1$?
A: Because at a multiple root, $f(x) \approx (x - x^*)^m \cdot g(x^*)$ with $g(x^*) \neq 0$. The Newton step $f/f' \approx (x - x^*)/m$, giving linear convergence with rate $1 - 1/m$. Fix: use the [modified Newton's method] $x_{n+1} = x_n - m \cdot f(x_n)/f'(x_n)$ — restores quadratic convergence if $m$ is known.

## 2.10 Polynomial Root Finding

Q: Why do polynomials get their own root-finding algorithms instead of using Newton on each root?
A: Because polynomials have structure that generic root-finders miss: all roots can be found simultaneously, complex roots come in conjugate pairs (for real coefficients), and deflation (dividing out known roots) reduces the problem. Dedicated polynomial methods exploit this for robustness and efficiency.

Q: How does the [Durand–Kerner] method find all roots of a polynomial simultaneously?
A: Start with $n$ distinct initial guesses $z_1^{(0)}, \dots, z_n^{(0)}$ (spread around a circle in $\mathbb{C}$). Iterate each guess by $z_i^{(k+1)} = z_i^{(k)} - p(z_i^{(k)}) / \prod_{j \neq i}(z_i^{(k)} - z_j^{(k)})$. Converges cubically to all roots simultaneously when well-initialized. Advantage over Newton: no root needs to be deflated.

Q: Why do classic Fortran/C libraries use [Jenkins–Traub] for polynomial root finding?
A: Because it combines three-stage shifting strategies to achieve global convergence (from any starting point), works directly with complex arithmetic, and handles near-multiple roots gracefully. `numpy.roots` uses a companion-matrix eigenvalue approach instead; Fortran/C libraries (SLATEC, RPOLY) use Jenkins–Traub.

Q: What does [Wilkinson's polynomial] demonstrate?
A: $W(x) = \prod_{k=1}^{20}(x - k)$ has simple integer roots $1, 2, \dots, 20$. Perturbing the coefficient of $x^{19}$ by $2^{-23}$ (one bit in single precision) moves some roots by over $10$! Shows that polynomial roots can be wildly ill-conditioned with respect to coefficients. Moral: root-finding accuracy depends on how the polynomial is REPRESENTED, not just its roots.

## 2.11 Pattern Recognition

Q: You need a root and you have a bracket $[a,b]$ with $f(a)f(b) < 0$ but no derivative. What method?
A: Brent's method — guaranteed convergence (bisection safety net) plus superlinear speed (secant + inverse quadratic interpolation). Default in SciPy's `brentq`.

Q: You need a root, you have $f'$ available cheaply, AND you have a good initial guess. What method?
A: Newton's method — quadratic convergence per iteration, fastest in this regime.

Q: You need a root and $f'$ is expensive or unavailable. What method?
A: Secant method — superlinear (order $\varphi \approx 1.618$) without derivative evaluations. One $f$ call per iteration vs. Newton's two.

Q: You need ALL roots of a polynomial of degree $n$. What method?
A: Companion-matrix eigenvalues (`numpy.roots`) or Jenkins–Traub. Avoid Newton-with-deflation — deflation amplifies roundoff and degrades later roots.

Q: You suspect a root has multiplicity $m > 1$ and Newton is converging linearly. What's the fix?
A: Modified Newton: $x_{n+1} = x_n - m \cdot f(x_n)/f'(x_n)$. Restores quadratic convergence when $m$ is known.

Q: Newton diverges or oscillates from your initial guess. What's the systematic fix?
A: Safeguarded Newton — accept Newton step only if it stays inside a bracket and reduces $|f|$, else fall back to bisection for that step.
