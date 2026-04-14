+++
order = 1
subject = "Math"
tags = ["math", "numerical-methods", "floating-point", "error-analysis", "conditioning", "stability"]
+++

# Numerical Methods — Floating Point and Error Analysis

## 1.1 Why Numerical Methods

Q: Why do scientific problems usually require numerical methods rather than analytical solutions?
A: Because most real models — nonlinear ODEs, high-dimensional PDEs, large linear systems, transcendental equations — have no closed-form solution, or closed forms that are impossible to evaluate by hand. Computers compute approximations with controlled error. Every field of computational science (physics, biology, engineering, finance, ML) runs on numerical algorithms.

Q: What's the difference between "numerical analysis" and "scientific computing"?
A: [Numerical analysis] studies the mathematics of approximate algorithms — convergence rates, stability, error bounds. [Scientific computing] is the practical discipline of using those algorithms at scale — performance, software engineering, hardware utilization. The line is blurry; most practitioners do both. This deck emphasizes the analysis side, with computational consequences noted.

## 1.2 Floating-Point Representation

C: A [floating-point number] is represented as $(-1)^s \cdot m \cdot \beta^e$ with sign bit $s$, mantissa (significand) $m$, exponent $e$, and base $\beta$ — typically $\beta = 2$ for binary hardware.

Q: Why did IEEE 754 standardize binary floating point across CPUs?
A: Because pre-standard, different machines used different bases, rounding rules, and special-value handling, making numerical code non-portable. IEEE 754 (1985, revised 2008) mandates binary representations, deterministic rounding, and defined behavior for $\pm \infty$, NaN, and subnormals — so the same program gives the same bit-exact result on any conforming machine.

Q: What are the sizes of IEEE 754 single and double precision?
A: [Single (`float32`)]: 1 sign bit + 8 exponent bits + 23 mantissa bits = 32 bits total, ~7 decimal digits of precision, range $\sim 10^{\pm 38}$. [Double (`float64`)]: 1 + 11 + 52 = 64 bits, ~16 decimal digits, range $\sim 10^{\pm 308}$. Double is the default for scientific computing; single is used for memory-bound ML workloads.

## 1.3 Machine Epsilon

C: [Machine epsilon] $\varepsilon_{\text{mach}}$ is the smallest positive number such that $1 + \varepsilon_{\text{mach}} \neq 1$ in floating-point arithmetic.

Q: What is $\varepsilon_{\text{mach}}$ for IEEE double precision?
A: $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}$. So double-precision arithmetic carries about $16$ correct decimal digits. Any two doubles differing by less than $\varepsilon_{\text{mach}}$ relative to their magnitude are indistinguishable after rounding.

Q: Why is $\varepsilon_{\text{mach}}$ a fundamental limit on algorithm accuracy?
A: Because rounding any real number to its nearest floating-point representation introduces relative error $\leq \varepsilon_{\text{mach}}/2$. Every arithmetic operation incurs at least this error. No algorithm can do better than $O(\varepsilon_{\text{mach}})$ relative accuracy per operation, so error accumulation sets the floor of achievable precision.

## 1.4 Rounding Errors

Q: What is the [unit roundoff]?
A: $u = \varepsilon_{\text{mach}}/2$ — the maximum relative error in representing a real number as a floating-point number under round-to-nearest. Every floating-point number $\tilde{x}$ satisfies $\tilde{x} = x(1 + \delta)$ with $|\delta| \leq u$. This is the building block of error analysis.

Q: What does IEEE 754's "round to nearest, ties to even" mean?
A: When a real number falls exactly between two floating-point numbers, round to the one with an even last bit. Alternatives (ties to away, ties to up) introduce systematic bias; "ties to even" is unbiased on average. Subtle but crucial for statistical accuracy over many operations.

## 1.5 Special Values

Q: What does IEEE 754's $\pm \infty$ represent, and when is it produced?
A: It represents unrepresentable magnitudes. Produced by overflow (e.g., $10^{308} \cdot 10$), division of a finite nonzero number by $\pm 0$, and $\log 0$. Crucially, $\infty$ obeys sensible rules: $\infty + 1 = \infty$, $1/\infty = 0$, so computations continue rather than crashing.

Q: Why does IEEE 754 distinguish $+0$ from $-0$?
A: Because signed zero preserves the direction of approach through intermediate computations. $1/(+0) = +\infty$ but $1/(-0) = -\infty$; $\text{atan2}(+0, -1) = \pi$ but $\text{atan2}(-0, -1) = -\pi$. The sign bit carries information even when the magnitude is zero — valuable in branch cuts of complex functions.

C: [NaN] (Not a Number) is IEEE 754's result for undefined operations like $0/0$, $\infty - \infty$, and $\sqrt{-1}$.

Q: Why does NaN propagate through arithmetic?
A: Because any operation with NaN returns NaN by design — so a single undefined step taints the whole computation, making errors detectable at the end rather than silently corrupting results. Contrast with returning a dummy value like $0$, which would hide the bug.

C: [Subnormal numbers] fill the gap between zero and the smallest normal float, giving [gradual underflow] that avoids abrupt loss of precision near zero.

Q: Why is $\text{NaN} \neq \text{NaN}$?
A: By IEEE 754 design — NaN represents "invalid/undefined," so no NaN equals anything, not even itself. This property lets you detect NaN via `x != x`. Also: NaNs in comparisons return false for $<, >, =, \leq, \geq$ — a subtlety that breaks sort algorithms if not handled.

## 1.6 Catastrophic Cancellation

C: [Catastrophic cancellation] occurs when subtracting two nearly-equal floating-point numbers, dramatically amplifying their relative errors.

Q: Why is $\sqrt{x + 1} - \sqrt{x}$ numerically unstable for large $x$?
A: Because both operands are nearly equal for large $x$, so their subtraction loses most of the significant digits. Example: $x = 10^{16}$ gives $\sqrt{10^{16} + 1} - 10^8 = 0$ in double precision — all accuracy gone. Fix: rationalize to $1/(\sqrt{x+1} + \sqrt{x})$, avoiding the cancellation.

P: Evaluate $f(x) = (1 - \cos x)/x^2$ at $x = 10^{-8}$ stably.
S:
**IDENTIFY**: Potential catastrophic cancellation in $1 - \cos x$ for small $x$.

**PLAN**: Use trigonometric identity to rewrite without cancellation.

**EXECUTE**:
1. Direct: $\cos(10^{-8}) \approx 1$ to machine precision, so $1 - \cos x \approx 0$, and $f(x)$ returns a garbage value.
2. Identity: $1 - \cos x = 2 \sin^2(x/2)$.
3. So $f(x) = 2 \sin^2(x/2) / x^2 = (\sin(x/2) / (x/2))^2 / 2$.
4. For $x = 10^{-8}$: $\sin(5 \times 10^{-9})/(5 \times 10^{-9}) \to 1$, so $f \to 1/2$ — the true limit value.

**EVALUATE**: Cancellation avoided by rewriting the expression. Real lesson: symbolic rearrangement can dramatically improve numerical stability — always check formulas for near-cancellations before computing.

## 1.7 Forward and Backward Error

C: For a computed solution $\tilde{y}$ to a problem $y = f(x)$: the [forward error] is $|\tilde{y} - y|$ (or relative). The [backward error] is the smallest $|\Delta x|$ such that $\tilde{y} = f(x + \Delta x)$ exactly — the perturbation to inputs needed to make $\tilde{y}$ the true answer.

Q: Why is backward error analysis more informative than forward error?
A: Because backward error separates the ALGORITHM's quality from the PROBLEM's sensitivity. A small backward error says "my algorithm solved a slightly-perturbed problem exactly," which is the best we can hope for. The forward error is then bounded by (backward error) × (condition number), decoupling algorithm correctness from problem conditioning.

## 1.8 Condition Number

C: The [condition number] of a problem measures how sensitive the output is to input perturbations. For $f: \mathbb{R} \to \mathbb{R}$ at $x$: $\kappa(x) = \left| \frac{x f'(x)}{f(x)} \right|$ — relative output change per unit relative input change.

Q: Why is a problem "ill-conditioned" fundamentally harder than one that is "well-conditioned"?
A: Because ill-conditioning is a property of the PROBLEM, not the algorithm. Even with a perfect algorithm, computing $f(x)$ at an ill-conditioned $x$ amplifies input noise by $\kappa(x)$ in the output. For $\kappa = 10^{10}$ and double precision, expect only $\sim 6$ correct digits no matter what algorithm you use. Diagnosing ill-conditioning is the first step.

Q: What is the condition number of matrix inversion?
A: For a linear system $Ax = b$, $\kappa(A) = \|A\| \|A^{-1}\|$. If $\kappa(A) = 10^k$, expect to lose $\sim k$ digits of accuracy in $x$ relative to input precision. Near-singular matrices have huge $\kappa$; Hilbert matrices are classic examples with $\kappa$ growing like $e^{3.5n}$.

## 1.9 Stability

C: An algorithm is [backward stable] if for every input $x$, the computed answer $\tilde{y}$ satisfies $\tilde{y} = f(x + \Delta x)$ for some small $\Delta x$ (of order $\varepsilon_{\text{mach}}$).

Q: What's the canonical example of an unstable algorithm?
A: Computing roots of $ax^2 + bx + c$ via $x = (-b \pm \sqrt{b^2 - 4ac})/(2a)$ when $b^2 \gg 4ac$. The root "$-b + \sqrt{\cdots}$" involves catastrophic cancellation. Stable alternative: compute the larger-magnitude root first, then the smaller via Viète's relation $x_1 x_2 = c/a$.

## 1.10 Absolute vs. Relative Error

Q: When should you report absolute error and when relative?
A: [Absolute error] $|e| = |\tilde{y} - y|$: appropriate when $y$ is near zero or lives on a fixed scale (e.g., probabilities, angles). [Relative error] $|e|/|y|$: appropriate when $y$ can vary over many magnitudes (physical measurements, scientific data). Relative error matches floating-point precision, which is itself relative.

## 1.11 Error Accumulation

Q: Why might summing $n$ numbers in naive order amplify error?
A: Because repeated addition adds $O(n)$ rounding errors, each of order $u \|x\|_\infty$. For large $n$ or poorly scaled inputs (some huge, some tiny), the small numbers get "lost" under the big ones, and final accuracy degrades as $O(n u)$. [Kahan summation] reduces this to $O(u)$ (independent of $n$) by tracking the running compensation term.

## 1.12 Choosing Precision

Q: When is `float32` sufficient vs. when do you need `float64`?
A: [`float32`] (single): 7 digits — OK for rendering, deep learning training (with careful loss scaling), signal processing with bounded dynamic range. [`float64`] (double): 16 digits — default for scientific computing, needed whenever you accumulate many values (Monte Carlo), solve ill-conditioned systems, or simulate long-time dynamics. For extreme precision: arbitrary-precision libraries (`mpmath`, MPFR) trade speed for arbitrary digits.
