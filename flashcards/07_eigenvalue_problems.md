+++
order = 7
subject = "Mathematics"
tags = ["math", "numerical-methods", "linear-algebra", "eigenvalues", "qr-algorithm", "power-iteration", "lanczos", "arnoldi"]
+++

# Numerical Methods — Eigenvalue Problems

## 7.1 The Eigenvalue Problem

C: The eigenvalue problem finds scalars $\lambda$ and nonzero vectors $\mathbf{v}$ with [$A \mathbf{v} = \lambda \mathbf{v}$]. $\lambda$ is an eigenvalue and $\mathbf{v}$ is the corresponding eigenvector.

Q: Why is eigenvalue computation inherently iterative, unlike linear-system solving?
A: Because eigenvalues are roots of the characteristic polynomial $\det(A - \lambda I) = 0$, which for degree $\geq 5$ has no closed-form solution (Abel–Ruffini). Every eigenvalue algorithm is therefore iterative: build a sequence converging to the spectrum. No algorithm can find eigenvalues exactly in finite time.

Q: Why are eigenvalues computed by factoring $A$ rather than by finding roots of $\det(A - \lambda I)$?
A: Because expanding the characteristic polynomial and then root-finding is catastrophically unstable (Wilkinson's polynomial — §2.10): eigenvalues can be well-conditioned while polynomial coefficients are not. Modern algorithms (QR, divide-and-conquer) work on $A$ directly via orthogonal transformations, preserving eigenvalues and numerical stability.

## 7.2 Power Iteration

C: [Power iteration] computes the dominant eigenvalue of $A$ by repeatedly multiplying a starting vector by $A$ and normalizing: $\mathbf{v}_{k+1} = A\mathbf{v}_k / \|A\mathbf{v}_k\|$.

Q: Why does power iteration converge to the eigenvector of the largest-magnitude eigenvalue?
A: Expand the starting $\mathbf{v}_0$ in the eigenbasis: $\mathbf{v}_0 = \sum_i c_i \mathbf{x}_i$. Then $A^k \mathbf{v}_0 = \sum_i c_i \lambda_i^k \mathbf{x}_i$. The $\lambda_1^k$ term (largest $|\lambda|$) dominates exponentially, so $A^k \mathbf{v}_0 / \lambda_1^k \to c_1 \mathbf{x}_1$. Normalization keeps magnitudes bounded.

Q: What is the convergence rate of power iteration?
A: Linear with ratio $|\lambda_2/\lambda_1|$, where $\lambda_1, \lambda_2$ are the two largest-magnitude eigenvalues. Fast convergence requires a well-separated dominant eigenvalue. When $|\lambda_1| \approx |\lambda_2|$, convergence is painfully slow, motivating [shift-and-invert] and [Rayleigh quotient iteration].

Q: Why does power iteration fail when the largest eigenvalues are complex conjugates of real $A$?
A: Because $|\lambda_1| = |\lambda_2|$ exactly — no unique dominant eigenvalue. Iterates oscillate rather than converge. Fix: use power iteration on $A^2$ (real but with squared separation), or switch to methods that handle complex pairs natively (QR algorithm, subspace iteration).

## 7.3 Inverse Iteration

Q: How does [inverse iteration] find the smallest-magnitude eigenvalue?
A: Apply power iteration to $A^{-1}$ — whose eigenvalues are $1/\lambda_i$, so the "dominant" eigenvalue of $A^{-1}$ is the smallest of $A$. Implementation: each step solves $A \mathbf{v}_{k+1} = \mathbf{v}_k$ (reuse an LU factorization). Converges at rate $|\lambda_{n-1}/\lambda_n|$ — the ratio of the two smallest eigenvalues.

Q: What does [shifted inverse iteration] add?
A: Replace $A$ with $A - \sigma I$ for a shift $\sigma$ close to a target eigenvalue. $A^{-1}$-iteration on this shifted matrix converges to the eigenvalue closest to $\sigma$, with rate depending on how much closer the target is than any other eigenvalue. Essential for finding specific eigenvalues — e.g., "find the eigenvalue closest to $5$."

Q: Why is [Rayleigh quotient iteration] cubically convergent for symmetric $A$?
A: Because it adaptively updates the shift each step: $\sigma_k = \mathbf{v}_k^T A \mathbf{v}_k / \mathbf{v}_k^T \mathbf{v}_k$ (the Rayleigh quotient approximating the eigenvalue). As the iterate approaches an eigenvector, the shift approaches the eigenvalue, and the shifted matrix becomes increasingly ill-conditioned — in exactly the direction that accelerates convergence. Near an eigenvector: each iteration roughly CUBES the error.

## 7.4 Deflation

Q: What is [deflation] in eigenvalue computation?
A: After finding an eigenvalue $\lambda_1$ and eigenvector $\mathbf{x}_1$, modify $A$ to remove them: e.g., $A' = A - \lambda_1 \mathbf{x}_1 \mathbf{x}_1^T$ (for symmetric $A$). Remaining eigenvalues of $A'$ match those of $A$ except $\lambda_1$ is replaced by $0$. Repeatedly apply power iteration / deflate to extract multiple eigenvalues. Loses numerical stability as errors accumulate.

## 7.5 The QR Algorithm

C: The QR algorithm iterates $A_k = Q_k R_k$, [$A_{k+1} = R_k Q_k$] — converging under mild assumptions to an upper-triangular (or quasi-triangular) matrix whose diagonal holds the eigenvalues.

Q: Why is the QR algorithm the gold standard for computing all eigenvalues of a dense matrix?
A: Because it is globally convergent (from any starting matrix), works on the matrix directly without forming characteristic polynomials, reveals the full spectrum, and with shifts achieves cubic convergence per eigenvalue. LAPACK's `dgeev` and `dsyev` are all QR-algorithm-based. Cost: $O(n^3)$ — practical up to $n \sim 10^4$.

Q: What preprocessing makes the QR algorithm practical?
A: First reduce $A$ to [upper Hessenberg] form (nonzero entries only on and above the first subdiagonal) via orthogonal similarity transforms — cost $O(n^3)$, but done once. Subsequent QR iterations preserve Hessenberg structure, reducing each iteration from $O(n^3)$ to $O(n^2)$. For symmetric $A$: reduce to tridiagonal — QR step becomes $O(n)$.

Q: Why are [shifts] essential in practical QR iteration?
A: Because unshifted QR converges only linearly with rate $|\lambda_{i+1}/\lambda_i|$ — slow when eigenvalues are close. Wilkinson shifts (chosen from the trailing $2 \times 2$ block) give cubic convergence for the trailing eigenvalue. Typical cost with shifts: $\sim 2{-}3$ QR steps per eigenvalue — the reason it's practical.

## 7.6 Subspace Methods for Large Matrices

Q: Why is QR inadequate for huge sparse eigenvalue problems?
A: Because QR requires $O(n^3)$ operations and $O(n^2)$ memory — infeasible for $n \geq 10^5$. Also destroys sparsity: reducing to Hessenberg fills the matrix densely. For large sparse problems needing only a few eigenvalues (dominant or near a shift), subspace methods like Lanczos and Arnoldi are used instead.

Q: What problem do [Lanczos] and [Arnoldi] iterations solve?
A: Extracting a few extremal eigenvalues (largest, smallest, or near a shift) of a huge sparse $A$ using only matrix-vector products. Build an orthonormal basis for a Krylov subspace; the small Hessenberg/tridiagonal matrix representing $A$ in this basis has eigenvalues (Ritz values) that converge to the extremal eigenvalues of $A$.

## 7.7 Lanczos Iteration

C: [Lanczos iteration] builds an orthonormal basis for $\mathcal{K}_k(A, \mathbf{v})$ for SYMMETRIC $A$, producing a tridiagonal matrix $T_k$ whose eigenvalues approximate $A$'s extremal eigenvalues.

Q: Why does Lanczos use a three-term recurrence?
A: Because symmetry forces the Hessenberg matrix to be tridiagonal — orthogonalizing $A\mathbf{v}_k$ against previous basis vectors requires only subtracting components against $\mathbf{v}_k$ and $\mathbf{v}_{k-1}$. Memory cost: three length-$n$ vectors regardless of iteration count. Crucial for huge problems.

Q: What is Lanczos's main numerical pitfall?
A: [Loss of orthogonality]: in finite precision, basis vectors gradually lose orthogonality to each other, causing spurious "ghost" eigenvalues (apparent multiplicities). Fixes: full reorthogonalization (expensive), selective reorthogonalization (ARPACK), or accepting ghosts and filtering them post-hoc.

## 7.8 Arnoldi Iteration

Q: What does [Arnoldi iteration] do differently from Lanczos?
A: Arnoldi extends Lanczos to NON-symmetric $A$: builds Krylov basis with full Gram–Schmidt orthogonalization (not just three-term). Produces a Hessenberg matrix $H_k$ (not tridiagonal). Memory and per-iteration cost grow linearly with $k$ — hence restarting strategies (IRAM: implicitly restarted Arnoldi) are used for long runs.

Q: Why is ARPACK the standard library for large sparse eigenvalue problems?
A: Because it implements [implicitly restarted Arnoldi] (and Lanczos for symmetric cases) with sophisticated reorthogonalization, shift-invert transformations for interior eigenvalues, and excellent stability. `scipy.sparse.linalg.eigs` wraps ARPACK. De facto standard from computational chemistry to PageRank.

## 7.9 Generalized Eigenvalue Problems

Q: What is the [generalized eigenvalue problem] $A\mathbf{x} = \lambda B \mathbf{x}$?
A: Find $\lambda, \mathbf{x}$ such that $(A - \lambda B)\mathbf{x} = \mathbf{0}$ with $B$ a second given matrix. Arises in: finite-element vibration analysis ($K\mathbf{x} = \lambda M \mathbf{x}$ with stiffness $K$ and mass $M$), control theory, and Fisher discriminant analysis. Generally reduced to a standard eigenproblem $B^{-1}A$ (if $B$ is SPD) or solved directly by QZ algorithm.

## 7.10 Singular Value Decomposition

C: The singular value decomposition (SVD) of $A \in \mathbb{R}^{m \times n}$ writes [$A = U\Sigma V^T$] with orthogonal $U, V$ and diagonal $\Sigma$ of nonnegative singular values.

Q: Why is SVD more fundamental than eigendecomposition for general matrices?
A: Because SVD exists for EVERY matrix (rectangular, singular, etc.) — eigendecomposition exists only for square matrices and may fail for non-diagonalizable ones. SVD reveals numerical rank, condition number ($\sigma_{\max}/\sigma_{\min}$), optimal low-rank approximations (Eckart–Young), and the [pseudo-inverse]. The swiss-army-knife of numerical linear algebra.

Q: How is SVD computed in practice?
A: Not by forming $A^T A$ (which squares condition number). Instead: bidiagonalize $A$ via Householder reflections ($A = U_1 B V_1^T$ with $B$ bidiagonal); apply a specialized QR-like algorithm to $B$. LAPACK's `dgesvd` implements this; one-sided Jacobi variants are used when extremely high accuracy is needed. Cost: $O(mn \min(m, n))$.

## 7.11 Sensitivity and Conditioning

Q: Why are eigenvalues of symmetric matrices extremely well-conditioned?
A: Because [Bauer–Fike theorem] for normal $A$ (symmetric $\Rightarrow$ normal) states that perturbing $A$ by $\Delta A$ shifts each eigenvalue by at most $\|\Delta A\|_2$. No amplification by condition number. For nonsymmetric matrices, eigenvalues can be wildly ill-conditioned (defective eigenvalues, Jordan blocks) — tiny perturbations cause large eigenvalue shifts.

Q: What is a [defective matrix] and why does it matter for eigenvalue computation?
A: A matrix whose eigenvectors do not span $\mathbb{R}^n$ — has a Jordan block of size $\geq 2$. Eigenvalues near a defective one are ultra-sensitive: a perturbation of $\varepsilon$ can shift eigenvalues by $\varepsilon^{1/k}$ (where $k$ is Jordan block size). Numerical algorithms effectively produce approximately diagonalizable matrices — treating defectiveness as a limiting case.

## 7.12 Choosing an Algorithm

Q: Which eigenvalue algorithm computes ALL eigenvalues of a small dense matrix?
A: The [QR algorithm] (`numpy.linalg.eig`).

Q: Which eigenvalue algorithms extract a few extremal eigenvalues of a large sparse matrix?
A: [Lanczos] (symmetric) or [Arnoldi] (general), both via ARPACK / `scipy.sparse.linalg.eigs`.

Q: Which eigenvalue method targets interior eigenvalues near a given shift for a large sparse matrix?
A: [Shift-invert] Arnoldi.

Q: Which eigenvalue method suffices when only the single dominant eigenvalue is needed?
A: [Power iteration].

Q: Which decomposition should be used for singular values or numerical rank of any matrix?
A: The [SVD] (`numpy.linalg.svd`).
