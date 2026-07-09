+++
order = 11
subject = "Mathematics"
tags = ["math", "numerical-methods", "pde", "finite-differences", "heat-equation", "wave-equation", "cfl"]
+++

# Numerical Methods — Finite Differences for PDEs

## 11.1 Why PDEs Matter

Q: Why do PDEs demand dedicated numerical methods beyond ODE techniques?
A: Because PDEs involve functions of multiple independent variables — $u(x, t)$ or $u(x, y, z, t)$ — producing discretizations with $O(N^d)$ unknowns for $N$ grid points per dimension. Memory, solver cost, and stability analysis all change qualitatively. Modeling heat flow, fluid dynamics, electromagnetics, quantum mechanics — every continuum physical process.

Q: How are PDEs classified?
A: Second-order linear: (i) [elliptic] — steady-state (Poisson $\nabla^2 u = f$); (ii) [parabolic] — diffusion-like (heat equation $u_t = \alpha u_{xx}$); (iii) [hyperbolic] — wave-like with finite propagation (wave equation $u_{tt} = c^2 u_{xx}$). Classification determines the right numerical method: elliptic needs global solves, parabolic combines IVP + BVP, hyperbolic needs CFL-respecting schemes.

## 11.2 Discretizing the Laplacian

Q: What is the standard finite-difference stencil for $\nabla^2 u$ in 2D?
A: The [5-point stencil]: $\nabla^2 u(x_i, y_j) \approx (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4 u_{ij})/h^2$ on a uniform grid with spacing $h$. Error $O(h^2)$. In 3D: 7-point stencil. Higher-order: 9-point (2D) or "wide" stencils for $O(h^4)$ — at the cost of more boundary complications.

Q: Why does the 5-point Laplacian give a sparse linear system?
A: Because each equation couples only 5 neighbors (center + 4 directions), so on an $N \times N$ grid the $N^2 \times N^2$ matrix has exactly $5 N^2$ nonzeros — $O(1)$ per row. Sparsity is essential: dense storage of a $10^6 \times 10^6$ matrix would need 8 TB, but the sparse representation fits in MB. Iterative solvers (CG, multigrid) exploit this.

## 11.3 Poisson's Equation

Q: How is [Poisson's equation] $-\nabla^2 u = f$ solved numerically on a rectangle?
A: Discretize with 5-point stencil + Dirichlet BCs → sparse symmetric positive-definite linear system $A \mathbf{u} = \mathbf{f}$. Solve with: (i) direct sparse factorization (MUMPS, PARDISO) for $n \lesssim 10^5$; (ii) [preconditioned CG] for medium $n$; (iii) [multigrid] for mesh-independent $O(n)$ complexity at any $n$. Multigrid is the gold standard.

## 11.4 Explicit Scheme for Heat Equation

C: The [FTCS (Forward-Time, Centered-Space)] scheme for $u_t = \alpha u_{xx}$: $u_i^{n+1} = u_i^n + \alpha \frac{\Delta t}{h^2}(u_{i+1}^n - 2 u_i^n + u_{i-1}^n)$.

Q: What is the [CFL-like] stability condition for FTCS on the heat equation?
A: $\alpha \Delta t / h^2 \leq 1/2$. Violating it: $u$ oscillates with exploding amplitude. Since $\Delta t$ scales as $h^2$, refining the grid by $2\times$ forces $4\times$ more time steps — $8\times$ total work for 1D, worse in higher dimensions. This scaling makes FTCS impractical for long-time heat simulation on fine grids.

Q: Why does the explicit heat scheme have $\Delta t \sim h^2$ stability whereas explicit wave schemes have $\Delta t \sim h$?
A: Because the heat equation has no finite propagation speed — information diffuses as $x \sim \sqrt{\alpha t}$, so to match spatial scales $h$ the time scale is $h^2/\alpha$. The wave equation has propagation speed $c$, so $\Delta t \sim h/c$. Parabolic PDEs are inherently stiff on fine grids; hyperbolic PDEs are not.

## 11.5 Implicit Schemes for Heat

Q: What is the [backward Euler scheme] for heat?
A: $u_i^{n+1} = u_i^n + \alpha \frac{\Delta t}{h^2}(u_{i+1}^{n+1} - 2 u_i^{n+1} + u_{i-1}^{n+1})$ — evaluate spatial derivative at new time. Implicit: each step requires solving a tridiagonal linear system. Unconditionally stable — any $\Delta t$ works — but only first-order accurate in time ($O(\Delta t)$).

Q: What is the [Crank–Nicolson scheme]?
A: Average of forward and backward Euler: $\frac{u_i^{n+1} - u_i^n}{\Delta t} = \frac{\alpha}{2 h^2}\left[(u_{i+1}^{n+1} - 2u_i^{n+1} + u_{i-1}^{n+1}) + (u_{i+1}^n - 2u_i^n + u_{i-1}^n)\right]$. Implicit, unconditionally stable, AND second-order accurate in both space and time ($O(\Delta t^2 + h^2)$). The default heat-equation scheme.

Q: Why can Crank–Nicolson produce oscillating solutions for rough initial data?
A: Because it is A-stable but NOT L-stable — the amplification factor $(1 - z/2)/(1 + z/2)$ approaches $-1$ (not $0$) as $z \to -\infty$. Sharp initial features excite high-frequency modes that should decay to zero quickly; Crank–Nicolson makes them alternate sign instead. Fix: one step of backward Euler first (the "Rannacher smoothing") or use DIRK schemes.

## 11.6 The Wave Equation

Q: How is the 1D [wave equation] $u_{tt} = c^2 u_{xx}$ discretized with leapfrog?
A: Centered differences in both time and space: $u_i^{n+1} = 2u_i^n - u_i^{n-1} + (c \Delta t / h)^2 (u_{i+1}^n - 2u_i^n + u_{i-1}^n)$. Three-time-level explicit scheme, $O(\Delta t^2 + h^2)$ accurate. Requires knowing two initial time levels; $u^0$ from the initial condition, $u^1$ from a Taylor expansion using $u_t(0)$.

Q: What is the [CFL condition] for the explicit wave scheme?
A: $c \Delta t / h \leq 1$ — the numerical domain of dependence must contain the physical one. If violated: numerical waves travel faster than physical waves, scheme unstable. The ratio $c\Delta t / h$ is the [Courant number]; optimal accuracy at Courant $= 1$ (exact propagation over one grid step per time step).

Q: Why is the CFL condition fundamental, not just technical?
A: Because a hyperbolic PDE has a well-defined [domain of dependence] — the past region whose data determines $u(x, t)$. An explicit numerical scheme uses only $\pm k$ grid points as its numerical domain of dependence. If the numerical domain is smaller than the physical one (i.e., $\Delta t > h/c$), the scheme literally lacks information to approximate $u$ — no scheme can be both explicit AND stable without CFL.

## 11.7 Upwind and Advection

Q: What is the [upwind scheme] for the advection equation $u_t + a u_x = 0$ with $a > 0$?
A: $u_i^{n+1} = u_i^n - a \frac{\Delta t}{h}(u_i^n - u_{i-1}^n)$ — one-sided difference in the direction information flows FROM. First-order accurate, stable under CFL $a \Delta t / h \leq 1$. The downwind scheme (using $u_{i+1}$) is unconditionally unstable — a warning that choice of stencil is as important as choice of order.

Q: Why does upwind introduce [numerical diffusion]?
A: Because the scheme is equivalent (to leading order) to the modified equation $u_t + a u_x = \frac{a h}{2}(1 - \nu) u_{xx}$ where $\nu = a\Delta t/h$ is the Courant number. A diffusion term appears as a truncation-error artifact. Sharp wavefronts smear over time — physically wrong for pure advection, physically natural-looking on under-resolved fluid problems.

Q: What makes higher-order advection schemes (Lax–Wendroff, MUSCL, WENO) valuable?
A: They reduce numerical diffusion by using more grid points per stencil. Lax–Wendroff is second-order accurate with no dissipation on smooth solutions. MUSCL and WENO adapt near discontinuities via slope limiters or essentially-non-oscillatory reconstructions, preventing spurious oscillations near shocks while preserving high order in smooth regions.

## 11.8 Operator Splitting

Q: What is [dimension splitting] for multidimensional PDEs?
A: Reduce a multidimensional implicit solve to a sequence of 1D solves: solve along $x$ (tridiagonal per row), then along $y$ (tridiagonal per column). Each 1D solve is $O(N)$ via Thomas, giving total $O(N^2)$ per time step for a 2D problem — vastly cheaper than solving the 2D implicit system directly. Classic example: [ADI] (Alternating-Direction Implicit) for the heat equation.

Q: What is the [Strang splitting] for $u_t = (A + B) u$?
A: Advance $\Delta t/2$ with operator $A$; $\Delta t$ with $B$; $\Delta t/2$ with $A$ again. Second-order accurate in time (vs. first-order for naïve sequential $A$-then-$B$ splitting). Preserves structural properties (conservation, symplecticity) better. Used when $A$ and $B$ capture physically distinct phenomena (advection + diffusion; kinetic + potential in quantum).

## 11.9 Stability Analysis

Q: What is [von Neumann stability analysis]?
A: For a linear scheme on a uniform grid: plug in a Fourier mode $u_i^n = g^n e^{ik x_i}$ and derive the amplification factor $g(k, \Delta t, h)$. Scheme is stable iff $|g(k)| \leq 1$ for all wavenumbers $k$. The workhorse tool — converts PDE scheme analysis into algebraic manipulation of a single mode. Standard for heat, wave, advection schemes.

Q: What is the [Lax equivalence theorem]?
A: For consistent finite-difference schemes applied to well-posed linear PDEs: STABILITY $\iff$ CONVERGENCE. Consistency alone is not enough — unstable schemes diverge even as $h \to 0$. Same spirit as the Dahlquist theorem for ODEs. The reason stability analysis matters so much: it's precisely what decides whether the scheme works at all.

## 11.10 Boundary Conditions

Q: How are [Neumann boundary conditions] $\partial u/\partial n = g$ discretized?
A: Two standard approaches: (i) [ghost points] — introduce a fictitious point outside the domain, set its value so the centered difference matches the BC, then use the interior stencil uniformly. (ii) One-sided differences at the boundary. Ghost points preserve $O(h^2)$ accuracy; naïve one-sided differences drop to $O(h)$ unless higher-order one-sided stencils are used.

Q: What are [absorbing boundary conditions] and why are they hard?
A: BCs that let waves LEAVE the computational domain without reflection — needed when simulating infinite domains on finite grids (electromagnetics, seismology). Perfect absorption is nonlocal in space and time; common approximations: [PML] (perfectly matched layer — absorbing buffer region), [Engquist–Majda] conditions, [Dirichlet-to-Neumann] maps. Imperfect absorption causes artificial reflections that contaminate solutions.

## 11.11 Finite Volume Methods

Q: How do [finite volume] methods differ from finite differences?
A: Integrate the PDE over control volumes (grid cells) rather than sampling at points. Unknowns are cell averages $\bar{u}_i$; the conservative form tracks fluxes across cell faces. Naturally handles conservation laws — discrete conservation holds exactly regardless of resolution. Standard for fluid dynamics (compressible flow, shocks) where mass/momentum/energy conservation is physically essential.

## 11.12 Practical Choices

Q: How do you pick a time-stepping scheme for a parabolic PDE in practice?
A: Short simulations on coarse grids → [FTCS explicit] (simple, fast if CFL permits). Long simulations, fine grids → [Crank–Nicolson] (unconditionally stable, 2nd-order). Very stiff, rough initial data → [backward Euler] or [BDF2] with smoothing. Multidimensional → [ADI] or [operator-splitting]. For nonlinear problems: use implicit schemes with Newton iteration; exploit sparsity aggressively.
