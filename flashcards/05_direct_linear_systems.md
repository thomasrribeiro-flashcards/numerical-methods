+++
order = 5
subject = "mathematics"
tags = ["math", "numerical-methods", "linear-algebra", "gaussian-elimination", "lu", "cholesky", "qr"]
+++

# Numerical Methods — Direct Solution of Linear Systems

## 5.1 The Linear System Problem

Q: What is the core problem of numerical linear algebra?
A: Solve $A \mathbf{x} = \mathbf{b}$ for $\mathbf{x}$ given a matrix $A \in \mathbb{R}^{n \times n}$ and vector $\mathbf{b} \in \mathbb{R}^n$. Fundamental because nearly every numerical algorithm reduces to a linear solve at some stage: finite-difference PDE discretizations, Newton iterations, least squares, eigenvalue problems.

Q: Why are linear systems not simply solved by $\mathbf{x} = A^{-1} \mathbf{b}$?
A: Because computing $A^{-1}$ explicitly costs $O(n^3)$ AND applies it inaccurately — explicit inversion is less stable than factoring $A$ and solving by back-substitution. Rule of thumb: never compute $A^{-1}$ just to apply it. Use LU / Cholesky / QR factorizations and solve with triangular systems.

## 5.2 Triangular Systems

C: A triangular system $T \mathbf{x} = \mathbf{b}$ can be solved in $O(n^2)$ by [back-substitution] when $T$ is upper-triangular, or by [forward-substitution] when $T$ is lower-triangular.

Q: Why is solving a triangular system cheap compared to a general system?
A: Because each equation introduces exactly one new unknown: in upper-triangular, the last row gives $x_n$ directly; substitute into row $n-1$ to find $x_{n-1}$; and so on. No search, no factoring — just $n$ steps of scalar division and multiplication each touching one row. Total cost $\sim n^2$ vs. $n^3$ for elimination.

Q: When does back-substitution on $T \mathbf{x} = \mathbf{b}$ fail numerically?
A: When a diagonal entry $T_{ii}$ is zero or extremely small relative to the other entries in its row, causing division blow-up. Singular or near-singular triangular systems signal that the original factorization encountered a zero (or near-zero) pivot — typically a sign that pivoting was needed upstream.

## 5.3 Gaussian Elimination

C: [Gaussian elimination] transforms $A \mathbf{x} = \mathbf{b}$ into an upper-triangular system $U \mathbf{x} = \mathbf{c}$ by applying row operations that zero out entries below each pivot.

Q: What is the cost of Gaussian elimination on an $n \times n$ system?
A: $O(n^3)$ for the elimination phase (more precisely $\sim \frac{2}{3} n^3$ floating-point operations), plus $O(n^2)$ for back-substitution. Dominant cost is the triple-nested loop over rows, columns, and eliminated entries. Doubling $n$ multiplies work by $8$.

Q: What is a [pivot] in Gaussian elimination?
A: The diagonal entry $a_{kk}^{(k)}$ used to zero out column $k$ in rows below — the divisor in the multipliers $a_{ik}/a_{kk}$. If the pivot is zero, elimination fails; if it's tiny relative to other column entries, the multipliers are huge and amplify roundoff catastrophically.

## 5.4 Pivoting

Q: Why is [partial pivoting] essential for numerical stability of Gaussian elimination?
A: Because without pivoting, a tiny pivot causes huge multipliers and severe roundoff amplification — even for well-conditioned $A$. Partial pivoting swaps rows at each step so the pivot is the largest-magnitude entry in its column below. Keeps multipliers bounded by $1$ in magnitude, making elimination backward-stable in practice.

Q: What does [complete pivoting] add over partial pivoting?
A: Complete pivoting searches the entire remaining submatrix (not just the current column) for the largest entry and swaps both rows AND columns. Better theoretical stability bounds but rarely used — the $O(n^3)$ search overhead and column permutation bookkeeping outweigh marginal stability gains on typical problems.

Q: When does partial pivoting fail to prevent catastrophic growth?
A: On specially constructed matrices (Wilkinson's "pivot-growth" example) where the growth factor reaches $2^{n-1}$. Pathological in theory; essentially never observed in practice for real-world matrices. Partial pivoting is therefore the industry-standard default — LAPACK's `dgesv` uses it.

## 5.5 LU Factorization

C: An LU factorization writes $A = LU$ where $L$ is [unit lower-triangular] and $U$ is [upper-triangular]. With partial pivoting: $PA = LU$ for a permutation matrix $P$.

Q: How does LU factorization reuse work when solving $A \mathbf{x}_i = \mathbf{b}_i$ for many right-hand sides?
A: Factor $A = LU$ once (cost $O(n^3)$). For each new $\mathbf{b}_i$: solve $L \mathbf{y} = \mathbf{b}_i$ by forward-substitution and $U \mathbf{x}_i = \mathbf{y}$ by back-substitution (cost $O(n^2)$ per solve). Crucial for problems where the same matrix appears repeatedly — implicit time stepping, Newton iteration with fixed Jacobian.

Q: Why is LU factorization just a bookkeeping reorganization of Gaussian elimination?
A: Because the multipliers $\ell_{ik} = a_{ik}^{(k)}/a_{kk}^{(k)}$ computed during elimination ARE the below-diagonal entries of $L$, and the eliminated matrix is exactly $U$. LU "stores the elimination history" as a factorization. Same cost, same stability behavior, but reusable.

## 5.6 Cholesky Factorization

C: A [Cholesky factorization] of a symmetric positive-definite (SPD) matrix $A$ writes $A = LL^T$ with $L$ lower-triangular — like LU but exploiting symmetry.

Q: Why is Cholesky factorization twice as fast as LU for SPD matrices?
A: Because symmetry lets us compute only $L$ (half the factors), halving both memory and work. Cost: $\sim \frac{1}{3} n^3$ vs. $\frac{2}{3} n^3$ for LU. Bonus: no pivoting required for SPD matrices — diagonals stay strictly positive throughout, so Cholesky is automatically stable.

Q: Why does Cholesky fail if $A$ is not positive-definite?
A: Because the algorithm takes a square root of a diagonal entry at each step: $\ell_{kk} = \sqrt{a_{kk} - \sum_{j < k} \ell_{kj}^2}$. If $A$ is not SPD, this quantity becomes zero or negative at some step. The failure of Cholesky factorization is itself a certificate of non-positive-definiteness — used in optimization to detect saddle points.

## 5.7 QR Factorization

C: A QR factorization writes $A = QR$ where $Q$ has [orthonormal columns ($Q^T Q = I$)] and $R$ is [upper-triangular].

Q: Why is QR factorization preferred over the normal equations for least squares?
A: Because solving $A^T A \mathbf{x} = A^T \mathbf{b}$ squares the condition number: $\kappa(A^T A) = \kappa(A)^2$. QR computes least squares directly by solving $R \mathbf{x} = Q^T \mathbf{b}$ with condition number only $\kappa(A)$. For ill-conditioned $A$, QR preserves digits that normal equations would lose.

Q: How is QR computed, and what is its cost?
A: Apply successive [Householder reflections] that zero out below-diagonal entries one column at a time. Each reflection is a rank-one update $H = I - 2\mathbf{v}\mathbf{v}^T/\|\mathbf{v}\|^2$. Total cost: $\sim \frac{4}{3} n^3$ for an $n \times n$ matrix (roughly twice LU). Stable without pivoting due to orthogonal transformations preserving norms.

Q: What's the [Gram–Schmidt] process and why isn't it used in production QR?
A: Gram–Schmidt orthogonalizes columns one at a time: subtract projections onto previously-computed $Q$ columns. Intuitive but numerically unstable — [classical Gram–Schmidt] loses orthogonality rapidly. [Modified Gram–Schmidt] is better but still inferior to Householder. Production QR uses Householder reflections or Givens rotations.

## 5.8 Sparse Linear Systems

Q: Why do sparse matrices demand specialized algorithms?
A: Because dense LU on an $n = 10^6$ matrix would need $8$ TB just to store the factors and $\sim 10^{18}$ operations. Sparse matrices (few nonzeros per row) arise everywhere — PDE discretizations, graph Laplacians, circuit analysis. Specialized sparse solvers exploit nonzero structure to reduce both memory and work by orders of magnitude.

Q: What is [fill-in] during sparse LU factorization?
A: New nonzero entries created in $L$ or $U$ at positions where $A$ was zero, produced when elimination adds multiples of a row to another. Fill-in is the enemy of sparse factorization: a matrix with $O(n)$ nonzeros can have $O(n^2)$ nonzeros in its LU factors if the elimination order is bad. Fill-minimizing orderings (approximate minimum degree, nested dissection) are essential.

Q: What does [nested dissection] ordering do for sparse factorization?
A: Recursively partition the graph of the matrix: find a small vertex separator whose removal splits the graph into roughly equal pieces; order separator vertices last; recurse on each piece. For 2D finite-element matrices, nested dissection achieves $O(n^{3/2})$ factorization cost vs. $O(n^3)$ for naïve orderings. State-of-the-art for large structured problems.

## 5.9 Special Structures

Q: How does [banded matrix] structure reduce LU cost?
A: A bandwidth-$p$ matrix has nonzeros only within $p$ diagonals of the main diagonal. LU preserves this bandwidth (no fill outside the band under partial pivoting with band-preserving algorithms). Cost drops from $O(n^3)$ to $O(n p^2)$; storage from $O(n^2)$ to $O(np)$. Finite-difference 1D problems produce tridiagonal ($p = 1$) systems solvable in $O(n)$.

Q: What is the [Thomas algorithm]?
A: Specialized Gaussian elimination for tridiagonal systems — three nonzero entries per row. Forward sweep eliminates subdiagonal; back substitution recovers solution. Total cost: $O(n)$ with tiny constants. Basis for every 1D PDE solver and for ADI (alternating-direction implicit) splitting in 2D/3D.

## 5.10 Backward Error and Residual

Q: Why is checking the [residual] $\|\mathbf{b} - A\tilde{\mathbf{x}}\|$ misleading for ill-conditioned systems?
A: Because a small residual does NOT guarantee a small forward error in $\tilde{\mathbf{x}}$: relative forward error is bounded by $\kappa(A)$ times relative residual. For $\kappa(A) = 10^{10}$, a residual of $10^{-14}$ (near machine precision) can coexist with solution error $10^{-4}$. Condition number, not residual, bounds achievable accuracy.

Q: What does [iterative refinement] do to a computed solution?
A: Given $\tilde{\mathbf{x}}$: compute residual $\mathbf{r} = \mathbf{b} - A\tilde{\mathbf{x}}$ in higher precision; solve $A \mathbf{d} = \mathbf{r}$ using the existing LU factorization; update $\tilde{\mathbf{x}} \leftarrow \tilde{\mathbf{x}} + \mathbf{d}$; repeat. Reduces forward error toward machine precision even for moderately ill-conditioned systems. Cheap because it reuses the factorization.

## 5.11 Condition Number of a Matrix

C: The [condition number] of a matrix is $\kappa(A) = \|A\| \|A^{-1}\|$, measuring how much $A \mathbf{x} = \mathbf{b}$ amplifies input perturbations into solution errors.

Q: How is $\kappa(A)$ estimated without explicitly forming $A^{-1}$?
A: Via [condition number estimators] (LAPACK's `dgecon`): apply $A^{-1}$ to carefully chosen test vectors using the existing LU factorization, and use a power-iteration-like algorithm on $A^{-1}$ to estimate $\|A^{-1}\|$. Cost: $O(n^2)$ — much cheaper than the $O(n^3)$ needed to form $A^{-1}$.

## 5.12 Choosing a Direct Solver

Q: Which direct solver should be used for a symmetric positive-definite linear system?
A: [Cholesky] — cheapest and most stable of the direct solvers.

Q: Which direct solver is the default for a general square linear system?
A: [LU with partial pivoting] — the default in `numpy.linalg.solve` and LAPACK's `dgesv`.

Q: When is QR the right direct solver for a linear system?
A: When $A$ is rectangular (least squares), very ill-conditioned, or when you want to avoid forming normal equations.

Q: When should sparse direct solvers (SuperLU, MUMPS, UMFPACK) be used?
A: Whenever $A$ has $\ll n^2$ nonzeros.
