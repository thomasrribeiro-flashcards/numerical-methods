+++
order = 12
subject = "Math"
tags = ["math", "numerical-methods", "fem", "spectral", "galerkin", "fft", "weak-form"]
+++

# Numerical Methods — Finite Element and Spectral Methods

## 12.1 Why Go Beyond Finite Differences

Q: What are the main limitations of finite differences that motivate FEM and spectral methods?
A: (i) Finite differences assume structured grids — awkward on complex geometries (airfoils, arteries, gears). (ii) Boundary conditions are fiddly, especially on curved boundaries. (iii) Local polynomial stencils cap accuracy at low order without special effort. FEM handles unstructured meshes; spectral methods achieve exponential convergence on smooth problems.

## 12.2 The Weak Formulation

Q: What is a [weak formulation] of a PDE?
A: Multiply the PDE by a test function $v$, integrate over the domain, and use integration by parts to shift derivatives onto $v$. Example: $-\nabla^2 u = f$ becomes $\int \nabla u \cdot \nabla v \, d\Omega = \int f v \, d\Omega$ for all $v$ with $v = 0$ on the Dirichlet boundary. The weak form demands less smoothness of $u$ and naturally encodes Neumann BCs.

Q: Why does integration by parts in the weak form reduce regularity requirements?
A: Because the strong form requires $u \in C^2$ (two classical derivatives), while integration by parts transfers one derivative onto $v$. The weak form only requires $u \in H^1$ (one weak derivative in $L^2$) — admitting continuous-but-non-smooth solutions like piecewise-linear interpolants. This relaxation is what makes FEM possible.

## 12.3 The Galerkin Method

C: The [Galerkin method] approximates $u \approx u_h = \sum_j c_j \phi_j$ in a finite-dimensional subspace $V_h = \text{span}\{\phi_j\}$, and enforces the weak form $\int \nabla u_h \cdot \nabla v \, d\Omega = \int f v$ for every basis function $v = \phi_i$.

Q: How does Galerkin projection produce a linear system?
A: Substitute $u_h = \sum_j c_j \phi_j$ into the weak form and use each $\phi_i$ as test function: $\sum_j c_j \int \nabla \phi_j \cdot \nabla \phi_i \, d\Omega = \int f \phi_i \, d\Omega$. This is $K \mathbf{c} = \mathbf{b}$ with [stiffness matrix] $K_{ij} = \int \nabla \phi_j \cdot \nabla \phi_i$ and load vector $b_i = \int f \phi_i$. Solve for coefficients $c_j$.

Q: Why does Galerkin give the [best approximation] in a natural energy norm?
A: Because the orthogonality condition $\int \nabla(u - u_h) \cdot \nabla \phi_i = 0$ for all $\phi_i$ means $u_h$ is the orthogonal projection of $u$ onto $V_h$ in the energy inner product $\langle u, v\rangle_E = \int \nabla u \cdot \nabla v$. Among all $v_h \in V_h$, $u_h$ minimizes the energy-norm error. Elegant optimality built into the method.

## 12.4 FEM Basis Functions

Q: What are standard FEM basis functions on a mesh?
A: [Piecewise polynomials]: on each element (triangle, quadrilateral, tetrahedron), $u_h$ is a polynomial of chosen degree $p$; across elements, $u_h$ is continuous ($C^0$). Lowest-order: piecewise-linear "hat" functions $\phi_i$ equal to $1$ at node $i$, $0$ at all other nodes, linear in between. Higher-order ($p = 2, 3, ...$) uses polynomial edge/face/interior nodes.

Q: Why do FEM basis functions have local support?
A: Because each $\phi_i$ is nonzero only on elements touching node $i$ — a handful of elements in any dimension. Consequence: the stiffness matrix is sparse ($O(1)$ nonzeros per row in 2D/3D). Local support is what makes FEM scalable: $O(n)$ memory and matrix-vector products for $n$ unknowns.

## 12.5 Mesh Assembly

Q: How is the global stiffness matrix [assembled] in FEM?
A: Loop over elements. For each element, compute a small [local stiffness matrix] $K^{(e)}$ of size $n_e \times n_e$ (e.g., $3 \times 3$ for a linear triangle) using quadrature of $\int \nabla \phi_j \cdot \nabla \phi_i$ over that element. Scatter-add $K^{(e)}$ entries into the global $K$ via a local-to-global node map. Fully parallelizable at the element level.

Q: What is [isoparametric mapping] in FEM?
A: Every physical element is the image of a fixed [reference element] (the "master triangle" or "master square") under a polynomial map using the same basis functions that represent the solution. Lets FEM handle curved elements, complex geometries, and computations via quadrature rules defined once on the reference element. The "iso" = same map for geometry and solution.

## 12.6 Quadrature in FEM

Q: Why must FEM integrals be computed by quadrature rather than symbolically?
A: Because after mapping curved physical elements to reference coordinates, the Jacobian of the map introduces rational functions — no closed form in general. Gaussian quadrature on the reference element handles this exactly for polynomial integrands of known degree. Choosing quadrature order to integrate the polynomial parts exactly is essential for optimal convergence rates.

## 12.7 A Priori Error Analysis

Q: What convergence rate does FEM with degree-$p$ elements achieve?
A: For smooth solutions on quasi-uniform meshes with element size $h$: error in the energy norm is $O(h^p)$; error in the $L^2$ norm is $O(h^{p+1})$. Higher degree = higher order, but also more unknowns per element and denser stiffness matrix. Tradeoff between [h-refinement] (more elements, fixed $p$) and [p-refinement] (fewer elements, higher $p$).

Q: What is [hp-adaptive] FEM?
A: A strategy that refines the mesh ($h$) in regions of low solution regularity (near singularities, boundary layers) and increases polynomial degree ($p$) in regions of high regularity. Optimally handles solutions with mixed smoothness — exponential convergence in smooth regions, algebraic convergence with local refinement near singularities. State-of-the-art for many engineering applications.

## 12.8 Spectral Methods

Q: What distinguishes [spectral methods] from FEM?
A: Spectral methods use GLOBAL basis functions (Chebyshev polynomials, Fourier modes, Legendre polynomials) spanning the entire domain — no mesh. Each basis function is nonzero everywhere. For smooth solutions on simple geometries: error decays faster than any polynomial in $N$ ([exponential convergence]), vs. FEM's algebraic $O(h^p)$ at fixed $p$.

Q: When do spectral methods give exponential convergence?
A: When the solution $u$ is smooth (analytic, in fact) AND the geometry is compatible with the chosen basis (periodic box for Fourier, interval or box for Chebyshev, sphere for spherical harmonics). Just a handful of modes can capture the solution to machine precision. Spectral methods dominate global climate models, turbulence DNS, and high-accuracy astrophysics simulations.

Q: What limits spectral methods' applicability?
A: (i) Complex geometries — no natural global polynomial basis on an airfoil-shaped region. (ii) Non-smooth solutions — Gibbs phenomenon near discontinuities destroys exponential convergence. (iii) Dense operators — global basis means dense matrices, so operation count per step scales worse than FEM in some cases. [Spectral element method] = FEM with high-order spectral bases on each element, merging the strengths.

## 12.9 FFT-Based Methods

Q: Why is the [FFT] ubiquitous in spectral methods for periodic problems?
A: Because the [fast Fourier transform] computes the discrete Fourier transform in $O(N \log N)$ time vs. $O(N^2)$ for naïve DFT. Spectral differentiation in Fourier space is exact and trivial ($\widehat{\partial_x u} = ik \hat{u}$). Spectral methods for periodic PDEs reduce to: FFT → multiply by $ik$ → inverse FFT. A single time step costs $O(N \log N)$ — nearly linear.

Q: What is a [pseudospectral] method?
A: Compute spatial derivatives in Fourier (or Chebyshev) coefficient space for exactness; compute nonlinear products in physical space to avoid convolution. Switch between spaces via FFT. Example: Burgers' equation $u_t + u u_x = \nu u_{xx}$ — compute $u_x$ and $u_{xx}$ spectrally, compute $u \cdot u_x$ in physical space. Avoids the $O(N^2)$ spectral convolution while retaining spectral accuracy.

## 12.10 Discontinuous Galerkin

Q: What is the [discontinuous Galerkin] (DG) method?
A: FEM-like but with polynomials on each element that are DISCONTINUOUS across element boundaries. Fluxes across element boundaries are handled by Riemann solvers (as in finite volume). Combines FEM's geometric flexibility with finite volume's conservation and shock-capturing. Standard for hyperbolic conservation laws: compressible Navier–Stokes, relativistic hydrodynamics, Maxwell's equations.

Q: Why is DG naturally parallel?
A: Because each element is coupled to neighbors ONLY through boundary fluxes — not through shared basis functions (unlike continuous FEM). Element-local computations dominate; inter-element communication is minimal. Scales to millions of cores; a favorite for exascale CFD codes. Newer methods (hybridizable DG, HDG) reduce the global problem size back toward continuous-FEM efficiency.

## 12.11 Boundary Element and Meshless Methods

Q: What is the [boundary element method] (BEM)?
A: Reformulate a PDE (typically Laplace or Helmholtz) as an integral equation on the boundary using Green's functions. Discretize only the BOUNDARY — reducing a 3D volume problem to a 2D surface problem. Unknowns: boundary values and normal derivatives. Cost: dense $O(N^2)$ matrix (Green's functions are nonlocal) — mitigated by fast methods (FMM, $\mathcal{H}$-matrices). Ideal for exterior problems (acoustic scattering).

Q: What are [meshless methods]?
A: Methods that discretize without connectivity-based meshes — using only scattered nodes + local approximation schemes (moving least squares, RBFs, smoothed-particle hydrodynamics). Avoid mesh-generation headaches, handle large deformations natively, useful for fluid–structure interaction and explosive flows. Less mature than FEM; numerical analysis still active research.

## 12.12 Choosing a Discretization

Q: How do you choose between FEM, spectral, finite difference, and finite volume?
A: [Finite difference]: simple geometries, moderate accuracy, quick prototyping. [Finite volume]: conservation laws (shocks, fluids). [FEM]: complex geometries, elliptic/parabolic PDEs, engineering. [Spectral]: simple geometries + smooth solutions + demand for extreme accuracy. [DG]: hyperbolic on complex geometries, massively parallel. [BEM]: exterior domains, homogeneous PDEs. Reality: production codes often mix methods via domain decomposition.

Q: What is [domain decomposition] in numerical PDE?
A: Split the computational domain into subdomains, solve on each subdomain independently, and iterate to match boundary conditions at subdomain interfaces. Parallelizes huge problems across cores/nodes; accommodates different physics or discretizations per subdomain (e.g., compressible fluid on one side, solid structure on the other). [Schwarz methods] and [FETI] are standard algorithms.
