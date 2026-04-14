+++
order = 4
subject = "Math"
tags = ["math", "numerical-methods", "differentiation", "integration", "quadrature", "simpson", "gauss"]
+++

# Numerical Methods — Differentiation and Integration

## 4.1 Numerical Differentiation

Q: Why is numerical differentiation inherently ill-conditioned?
A: Because the definition $f'(x) \approx (f(x + h) - f(x))/h$ combines two errors: [truncation error] $O(h)$ (Taylor-series neglect) decreases with smaller $h$; [roundoff error] $O(\varepsilon_{\text{mach}}/h)$ INCREASES with smaller $h$ (catastrophic cancellation in the numerator). Optimal $h$ balances both, giving accuracy $\sim \sqrt{\varepsilon_{\text{mach}}} \approx 10^{-8}$ (not $10^{-16}$).

## 4.2 Finite Difference Formulas

C: The [forward difference] approximation is $f'(x) \approx (f(x+h) - f(x))/h$ with error $O(h)$.

C: The [centered difference] approximation is $f'(x) \approx (f(x+h) - f(x-h))/(2h)$ with error $O(h^2)$.

Q: Why does the centered difference achieve $O(h^2)$ accuracy vs. forward's $O(h)$?
A: Because in the Taylor expansions of $f(x \pm h)$, the odd-order terms cancel when subtracted: $f(x+h) - f(x-h) = 2h f'(x) + \frac{h^3}{3}f'''(x) + \dots$. Dividing by $2h$ leaves leading error of order $h^2$. Symmetric stencils kill odd-order derivatives — a general principle.

Q: What's the second-derivative central-difference formula?
A: $f''(x) \approx (f(x+h) - 2f(x) + f(x-h))/h^2$ with error $O(h^2)$. Derivation: add the Taylor expansions of $f(x \pm h)$ — the first-derivative terms cancel, leaving $2f(x) + h^2 f''(x) + O(h^4)$. Rearrange to isolate $f''$.

## 4.3 Richardson Extrapolation

Q: How does [Richardson extrapolation] improve a finite-difference approximation?
A: If $D(h) = f'(x) + ch^p + O(h^q)$ with $q > p$, then combining $D(h)$ and $D(h/2)$ eliminates the $h^p$ term: $\hat{D} = (2^p D(h/2) - D(h))/(2^p - 1) = f'(x) + O(h^q)$. Doubles the convergence order using just one extra evaluation. Used iteratively in Romberg integration.

## 4.4 The Integration Problem

Q: What is [numerical quadrature]?
A: Approximating a definite integral $\int_a^b f(x) \, dx$ by a weighted sum of function evaluations $\sum_{i=0}^{n} w_i f(x_i)$. The [nodes] $x_i$ and [weights] $w_i$ define the quadrature rule. Different choices give different accuracy and efficiency tradeoffs.

## 4.5 Newton–Cotes Rules

C: [Newton–Cotes] quadrature rules use equispaced nodes and derive weights by integrating the Lagrange interpolating polynomial.

Q: What are the [trapezoidal rule] and its error?
A: $\int_a^b f(x) \, dx \approx \frac{b - a}{2}(f(a) + f(b))$, exact for linear $f$. Error: $-\frac{(b-a)^3}{12} f''(\xi)$ for some $\xi \in (a, b)$. On a uniform grid with $n$ subintervals, composite trapezoidal has error $O(h^2)$ where $h = (b-a)/n$.

Q: What are [Simpson's rule] and its error?
A: $\int_a^b f(x) \, dx \approx \frac{b-a}{6}\left[f(a) + 4f\left(\tfrac{a+b}{2}\right) + f(b)\right]$ — integrates the quadratic through the three points. Error: $-\frac{(b-a)^5}{2880} f^{(4)}(\xi)$. Exact for cubics (bonus order from symmetric cancellation). Composite Simpson's: $O(h^4)$.

Q: Why do Newton–Cotes rules fail at high order?
A: Because they use EQUISPACED nodes, triggering Runge's phenomenon: for $n \geq 8$ or so, some weights become negative, amplifying roundoff and making composite rules stable only at low order. In practice: use composite low-order rules (trapezoidal, Simpson's) or switch to Gauss quadrature for high-accuracy single-interval integration.

## 4.6 Composite Rules

C: A [composite quadrature rule] applies a low-order rule on subintervals of a partition and sums, giving accuracy that scales with the number of subintervals.

Q: How does composite Simpson's rule achieve $O(h^4)$?
A: Partition $[a, b]$ into $n$ even subintervals of width $h = (b-a)/n$. On each pair of adjacent subintervals apply Simpson's rule (cost: $n/2$ quadratic fits). Total: $\int \approx \frac{h}{3}\left[f_0 + 4(f_1 + f_3 + \dots) + 2(f_2 + f_4 + \dots) + f_n\right]$. Local error $O(h^5)$ per pair × $n/2$ pairs → global error $O(h^4)$.

## 4.7 Gauss Quadrature

C: [Gauss quadrature] chooses nodes $x_i$ AND weights $w_i$ to maximize the polynomial degree exactly integrated — a $n$-point Gauss rule is exact for polynomials of degree $\leq 2n - 1$.

Q: Why does Gauss quadrature double the achievable accuracy for the same number of points?
A: Because Newton–Cotes fixes the nodes, giving $n$ degrees of freedom (the weights). Gauss adjusts BOTH nodes and weights, giving $2n$ degrees of freedom — enough to exactly integrate polynomials of degree $2n - 1$. The optimal nodes turn out to be roots of orthogonal polynomials (Legendre on $[-1, 1]$, Hermite on $\mathbb{R}$ with weight $e^{-x^2}$, etc.).

Q: What are [Gauss–Legendre nodes and weights]?
A: For the interval $[-1, 1]$ with weight $w(x) = 1$: nodes are roots of the Legendre polynomial $P_n(x)$; weights $w_i = 2/[(1 - x_i^2)(P_n'(x_i))^2]$. Tabulated (`scipy.special.roots_legendre`) or computed via the Golub–Welsch algorithm (eigenvalues of a tridiagonal Jacobi matrix). Generalize to any $[a, b]$ via affine change of variable.

## 4.8 Adaptive Quadrature

Q: What is [adaptive quadrature] and when is it essential?
A: A method that recursively subdivides intervals where the integrand has sharp features, using error estimates to decide where to refine. Compare the result on $[a, b]$ to the sum on $[a, (a+b)/2] + [(a+b)/2, b]$; if they differ by more than tolerance, split and recurse. Essential for integrands with narrow peaks, logarithmic singularities, or highly nonuniform behavior.

P: Integrate $\int_0^1 \sqrt{x} \, dx$ using Simpson's rule with 2 subintervals.
S:
**IDENTIFY**: $f(x) = \sqrt{x}$, $[a, b] = [0, 1]$, $n = 2$ subintervals → $h = 0.5$, use nodes $x_0 = 0, x_1 = 0.5, x_2 = 1$.

**PLAN**: Apply Simpson's formula $\frac{h}{3}[f_0 + 4f_1 + f_2]$.

**EXECUTE**:
1. $f(0) = 0$; $f(0.5) = \sqrt{0.5} \approx 0.7071$; $f(1) = 1$.
2. Simpson: $\frac{0.5}{3}[0 + 4(0.7071) + 1] = \frac{1}{6}[3.8284] \approx 0.6381$.

**EVALUATE**: Exact answer $\int_0^1 \sqrt{x} \, dx = 2/3 \approx 0.6667$. Simpson's error $\approx 0.029$ — notable because $\sqrt{x}$ has an unbounded derivative at $x = 0$, breaking Simpson's $O(h^4)$ error bound (requires $f^{(4)}$ bounded). Endpoint singularities slash Simpson's order; adaptive methods or variable transformations are needed.

## 4.9 Romberg Integration

Q: How does [Romberg integration] combine trapezoidal rule with Richardson extrapolation?
A: Compute trapezoidal approximations $T_h, T_{h/2}, T_{h/4}, \dots$. Build a Richardson table: $R_{k, j} = R_{k, j-1} + (R_{k, j-1} - R_{k-1, j-1})/(4^j - 1)$. The diagonal $R_{n,n}$ achieves error $O(h^{2n+2})$. Equivalent to integrating the unique polynomial through the trapezoidal approximations. Elegant and high-order without Gauss-rule node tables.

## 4.10 Multidimensional Integration

Q: Why does multidimensional quadrature suffer the [curse of dimensionality]?
A: Because a tensor-product rule with $n$ points per dimension uses $n^d$ total points in $d$ dimensions. For $n = 10$ and $d = 6$: $10^6$ function evaluations; for $d = 20$: $10^{20}$ — infeasible. Even worse, error decays more slowly relative to work. Classical quadrature is practical only for $d \leq 4$ or so.

Q: What's the alternative for high-dimensional integrals?
A: [Monte Carlo integration]: $\int_\Omega f \, dx \approx |\Omega| \cdot \frac{1}{N} \sum_{i=1}^{N} f(X_i)$ with $X_i$ uniformly random in $\Omega$. Error decays as $O(1/\sqrt{N})$ — INDEPENDENT of dimension. Slow per-dimension but unbeatable for $d \gg 5$. [Quasi-Monte Carlo] uses low-discrepancy sequences (Sobol, Halton) for better rates ($O((\log N)^d / N)$) in moderate dimensions.

## 4.11 Improper and Singular Integrals

Q: How does [truncation] handle an infinite integration range?
A: Cut off the tail at a large $R$ and approximate $\int_0^\infty f(x) \, dx \approx \int_0^R f(x) \, dx$. Works only if $f$ decays fast enough that $\int_R^\infty f \, dx$ is negligible. Simple but requires estimating the truncation error — typically by increasing $R$ until the result stabilizes to the desired tolerance.

Q: How can a change of variable turn an improper integral into a proper one?
A: Substitute a transformation that maps $[0, \infty)$ to a bounded interval, e.g., $x = t/(1 - t)$ sends $[0, 1) \to [0, \infty)$ with $dx = dt/(1-t)^2$. After substitution, apply any standard quadrature rule on $[0, 1]$. Effective when the integrand is smooth after the transform; careful with singularities introduced at $t = 1$.

Q: What is [Gauss–Laguerre quadrature]?
A: A specialized rule for $\int_0^\infty f(x) e^{-x} \, dx \approx \sum_i w_i f(x_i)$ where nodes $x_i$ are roots of the Laguerre polynomial $L_n(x)$. Absorbs exponential decay into the weight function, so an $n$-point rule is exact for $f$ polynomial of degree $\leq 2n - 1$. Ideal for integrals of the form $\int_0^\infty g(x) e^{-x} \, dx$ or integrands with similar tail behavior.

Q: What is [Gauss–Hermite quadrature]?
A: A specialized rule for $\int_{-\infty}^{\infty} f(x) e^{-x^2} \, dx \approx \sum_i w_i f(x_i)$ where nodes are roots of the Hermite polynomial $H_n(x)$. Natural for integrals weighted by a Gaussian — statistics, quantum mechanics (harmonic oscillator), and Bayesian inference with Gaussian priors.

## 4.12 Error Estimation in Practice

Q: How do you estimate the error of a computed integral WITHOUT knowing the exact answer?
A: Compute with two rules (or two step sizes) and compare. If rule $A$ has error $O(h^p)$ and rule $B$ has error $O(h^q)$ with $q > p$, then $|A - B|$ is an estimate of $A$'s error. Adaptive quadrature uses this idea recursively; embedded rules (like Gauss–Kronrod, which extends a Gauss rule with extra nodes to give two estimates cheaply) are the practical choice.
