+++
order = 3
subject = "Math"
tags = ["math", "numerical-methods", "interpolation", "lagrange", "splines", "chebyshev", "runge"]
+++

# Numerical Methods — Interpolation

## 3.1 The Interpolation Problem

Q: What is the [interpolation problem]?
A: Given $n + 1$ data points $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$ with distinct $x_i$, find a function $p$ (usually from a specific family — polynomials, splines, trigonometric) with $p(x_i) = y_i$ for every $i$. Interpolation reconstructs a function from discrete samples.

Q: How does interpolation differ from regression?
A: Interpolation: $p(x_i) = y_i$ EXACTLY for every data point. Regression: $p$ minimizes some measure of error (e.g., sum of squared residuals) across data points but doesn't necessarily match them. Use interpolation when data is exact (physical constants, computer-generated); use regression when data is noisy.

## 3.2 Polynomial Interpolation: Existence and Uniqueness

Q: Why is there exactly one polynomial of degree $\leq n$ interpolating $n + 1$ points?
A: Existence: Lagrange's construction below gives one explicitly. Uniqueness: two such polynomials $p, q$ would have $p - q$ of degree $\leq n$ with $n + 1$ zeros — forcing $p - q = 0$ by the fundamental theorem of algebra. Interpolating polynomial is unique up to representation.

## 3.3 Lagrange Interpolation

C: The [Lagrange form] of the interpolating polynomial is $p(x) = \sum_{i=0}^{n} y_i L_i(x)$, where $L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$.

Q: Why does each Lagrange basis function $L_i$ satisfy $L_i(x_j) = \delta_{ij}$?
A: Because $L_i$'s numerator contains $(x - x_j)$ for every $j \neq i$, so $L_i(x_j) = 0$ for $j \neq i$. At $x = x_i$, the factors all become $(x_i - x_j)/(x_i - x_j) = 1$, giving $L_i(x_i) = 1$. The "indicator polynomials" structure makes the interpolation formula immediate.

## 3.4 Newton's Divided Differences

Q: Why does the [Newton form] of interpolation using [divided differences] generalize Lagrange?
A: Newton's form: $p(x) = \sum_{i=0}^{n} [y_0, \dots, y_i] \prod_{j=0}^{i-1}(x - x_j)$, where $[y_0, \dots, y_i]$ is the $i$-th divided difference. Advantage over Lagrange: adding a new data point $(x_{n+1}, y_{n+1})$ requires only ONE new coefficient, not re-deriving everything. Better for incremental data.

## 3.5 The Runge Phenomenon

C: [Runge's phenomenon]: high-degree polynomial interpolation at equispaced points can produce large oscillations near the interval endpoints, causing $\|f - p_n\|_\infty \to \infty$ as $n \to \infty$.

Q: Why does Runge's phenomenon make equispaced high-degree interpolation dangerous?
A: Because the interpolation error grows dramatically near the ends of the interval even for smooth $f$. Classic example: $f(x) = 1/(1 + 25x^2)$ on $[-1, 1]$ with equispaced nodes — the interpolant's oscillation amplitude at the endpoints grows exponentially with $n$. Adding more points makes the fit WORSE, not better.

Q: How do [Chebyshev nodes] fix Runge's phenomenon?
A: By placing nodes at $x_k = \cos((2k+1)\pi/(2n+2))$ — clustered near the endpoints where equispaced interpolation fails. The clustering matches the inverse density of Chebyshev polynomial zeros, minimizing the max of $|\prod(x - x_k)|$ and giving uniform convergence as $n \to \infty$ for smooth $f$.

Q: How does [spline interpolation] sidestep Runge's phenomenon?
A: By using LOW-DEGREE polynomials on subintervals rather than one high-degree global polynomial. Runge oscillations are a high-degree artifact; cubic splines limit degree $\leq 3$ per subinterval, refining accuracy by adding more subintervals rather than raising degree. Global artifacts can't accumulate.

## 3.6 Chebyshev Interpolation

C: The [Chebyshev nodes of the first kind] on $\lbrack -1, 1\rbrack $ are $x_k = \cos\left(\frac{(2k + 1)\pi}{2(n+1)}\right)$ for $k = 0, 1, \dots, n$ — zeros of the degree-$(n+1)$ Chebyshev polynomial.

Q: Why do Chebyshev nodes give near-optimal polynomial interpolation?
A: Because they MINIMIZE the maximum of the error polynomial $\prod (x - x_k)$ over $[-1, 1]$. By the error formula $f - p = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod(x - x_k)$, minimizing the product minimizes the worst-case error. For smooth $f$, Chebyshev interpolation converges uniformly with rate comparable to the best polynomial approximation.

## 3.7 Piecewise Linear and Cubic Splines

C: A [spline] is a piecewise polynomial: on each subinterval $\lbrack x_i, x_{i+1}\rbrack $, the spline is a polynomial of a fixed degree $k$, joined at the nodes with continuity up to derivative order $k - 1$.

Q: Why are [cubic splines] $C^2$-continuous at the interior nodes?
A: Because the piecewise cubics are joined with matching value, first derivative, AND second derivative at each interior $x_i$. Two continuous derivatives make the curve visually smooth — no visible kinks or curvature jumps — which is why splines look more natural than piecewise linear or piecewise quadratic interpolants.

Q: What minimization property do cubic splines satisfy?
A: Among all $C^2$ interpolants, the natural cubic spline minimizes the bending energy $\int_a^b (p''(x))^2 \, dx$. This is the "thin flexible beam" property: a physical draftsman's spline (a thin wooden strip pinned at data points) takes exactly this shape because it minimizes elastic bending energy.

Q: Why are cubic splines "local" in a way high-degree global polynomials are not?
A: Because changing one data point only perturbs the cubic pieces in a small neighborhood around that point — the influence decays exponentially away from the change. High-degree global interpolation has the opposite property: every coefficient depends on every data point, so a single perturbation reshapes the entire curve.

C: A [natural cubic spline] enforces $p''(x_0) = p''(x_n) = 0$ — zero curvature at the endpoints.

Q: What endpoint condition does a [clamped cubic spline] impose?
A: Specified values of the FIRST derivative at the endpoints: $p'(x_0) = y'_0$ and $p'(x_n) = y'_n$. Use clamped splines when you know the slope at the boundary (e.g., from physical symmetry or a separate derivative measurement); they eliminate the boundary curvature artifacts of natural splines.

Q: What is the [not-a-knot] endpoint condition?
A: Force the third derivative to be continuous at $x_1$ and $x_{n-1}$ — i.e., the first two cubic pieces are actually the SAME cubic, and likewise for the last two. This removes the "knot" between them, giving smoother behavior near boundaries when no derivative data is available. `scipy.interpolate.CubicSpline` uses not-a-knot by default.

## 3.8 Hermite Interpolation

Q: What extra information does [Hermite interpolation] use?
A: Both function values AND derivatives at each node. For $n + 1$ points with values $y_i$ and derivatives $y'_i$, Hermite interpolation produces the unique polynomial of degree $\leq 2n + 1$ matching all of them. Natural for trajectory fitting where position and velocity are both known.

## 3.9 Trigonometric Interpolation

Q: What is [trigonometric interpolation] and when do you use it?
A: Fitting a finite Fourier series $p(x) = \sum_{k=-N}^{N} c_k e^{ikx}$ to $2N + 1$ equispaced samples on $[0, 2\pi]$. Used for periodic data (signals, climate cycles). The [discrete Fourier transform] (DFT) computes the coefficients in $O(N \log N)$ via the FFT. Avoids Runge's phenomenon for periodic functions.

## 3.10 Multivariate Interpolation

Q: How does [tensor-product interpolation] generalize 1D interpolation to a 2D grid?
A: Interpolate along one axis first (e.g., rows), producing a 1D function per column; then interpolate those along the other axis. Using linear 1D pieces gives [bilinear interpolation]; using cubic pieces gives [bicubic interpolation]. Standard in image resampling and texture mapping because grid-aligned data makes tensor products cheap.

Q: What is a [radial basis function] (RBF) interpolant for scattered data?
A: $p(\mathbf{x}) = \sum_i c_i \phi(\|\mathbf{x} - \mathbf{x}_i\|)$ — a weighted sum of radially symmetric basis functions $\phi$ (Gaussian, multiquadric, thin-plate spline) centered at the data points. Coefficients $c_i$ solve a linear system enforcing $p(\mathbf{x}_i) = y_i$. Works for arbitrary scattered points in any dimension.

Q: How does [Delaunay-triangulation-based] interpolation handle scattered 2D data?
A: Build a Delaunay triangulation of the data points; on each triangle, use linear (or higher-order) interpolation among its three vertex values. Fast, piecewise-linear, and produces a continuous $C^0$ interpolant. The go-to method for unstructured 2D/3D data in finite-element preprocessors and `scipy.interpolate.griddata`.

Q: Why does interpolation become expensive in high dimensions?
A: Because tensor-product grids grow as $N^d$ with dimension $d$ — a modest $N=20$ along each axis gives $20^{10} \approx 10^{13}$ points in 10 dimensions. The "[curse of dimensionality]." High-dimensional interpolation requires either sparse grids, low-rank tensor methods, or assumed structure (low intrinsic dimension, additivity).

## 3.11 Error Analysis

Q: What is the error formula for polynomial interpolation?
A: If $f \in C^{n+1}[a, b]$ and $p$ interpolates $f$ at nodes $x_0, \dots, x_n$: for any $x \in [a, b]$, there exists $\xi \in (a, b)$ with $f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^{n}(x - x_i)$. Error depends on $f$'s smoothness AND the node product. Chebyshev nodes minimize the node-product contribution.

## 3.12 Interpolation vs. Approximation

Q: When should you NOT use interpolation?
A: When data is noisy — interpolation forces $p$ through every data point, amplifying noise into wild oscillations. Use regression (least squares) or smoothing splines instead. Also: when you have far more data than needed to determine a simple model — interpolation wastes information. Interpolation is for exact, sparse, clean data.
