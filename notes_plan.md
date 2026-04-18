# Optimal Control Notes -- Implementation Plan

## Context

Comprehensive, self-contained LaTeX reference notes on Optimal Control. Covers all topics from the question bank (`/scratch/phd/prep/oc/claude.tex`, 128 questions across 9 parts) plus missing topics: constrained optimization, Pontryagin's Minimum Principle, HJB equation, shooting methods, and MPC. Template from `Note/` directory.

---

## Lecture Plan (9 Lectures)

### Lecture 1: Convex Sets, Convex Functions, and Optimization Foundations

**Covers:** Bank Part A (Q1--Q14) + Part B (Q15--Q26)

**Narrative:** Build the mathematical bedrock. Convexity is the single property that makes optimization tractable -- start here, then immediately connect to optimality conditions.

#### Subsections

**1.1 Convex Sets**
- `definition`: Convex set
- `proposition` + `proof`: Intersection of convex sets is convex (Q1)
- `eg`: Closed ball is convex; unit circle is not; solution set of $Ax = b$; PSD cone (Q4)
- `remark`: Union of convex sets is generally not convex; condition for convexity of union (Q5)

**1.2 Convex Functions**
- `definition`: Convex function, strictly convex function (Q2)
- `eg`: $f(x) = |x|$ is convex but not strictly convex; affine functions
- `theorem` + `proof`: Epigraph characterization -- $f$ convex iff $\epi(f)$ is a convex set (Q3)
- `theorem` + `proof`: First-order condition -- $f(y) \geq f(x) + \nabla f(x)^T(y-x)$ for differentiable convex $f$ (Q7)
- `intuition`: Geometric meaning -- tangent hyperplane is a global underestimator
- `theorem` + `proof`: Second-order condition -- $f$ convex iff $\nabla^2 f(x) \succeq 0$ for $C^2$ functions (Q6)
- `eg`: Checking convexity of specific functions via Hessian (Q11)

**1.3 Operations Preserving Convexity**
- `proposition` + `proof` for each: non-negative weighted sums, pointwise max, composition with affine map, composition rules (Q9)
- `remark`: Why these matter -- building complex convex functions from simple ones

**1.4 Strong Convexity and Conjugate Functions**
- `definition`: $m$-strong convexity (Q10)
- `theorem` + `proof`: Strong convexity lower bound $f(y) \geq f(x) + \nabla f(x)^T(y-x) + \frac{m}{2}\|y-x\|^2$
- `corollary`: Uniqueness of minimizer for strongly convex functions
- `definition`: Legendre-Fenchel conjugate $f^*(y) = \sup_x \{y^T x - f(x)\}$ (Q14)
- `proposition` + `proof`: $f^*$ is always convex
- `eg`: Conjugate of $f(x) = \frac{1}{2}x^T Q x$ with $Q \succ 0$

**1.5 Optimality Conditions for Unconstrained Problems**
- `theorem` + `proof`: FONC -- $\nabla f(x^*) = 0$ at local min (Q15)
- `theorem`: SONC -- $\nabla^2 f(x^*) \succeq 0$; SOSC -- $\nabla^2 f(x^*) \succ 0$ (Q16)
- `proof`: SOSC implies strict local minimum (Q20)
- `remark`: The gap between SONC and SOSC; $f(x)=x^4$ at origin (Q8); $f(x)=x^3$ trap (Q25)
- `theorem` + `proof`: Every local minimum of a convex function is global (Q12)
- `eg`: Classify critical points of $f(x_1,x_2) = 2x_1^2 + x_2^2 - 2x_1x_2 + 4x_1 - 6x_2$ (Q18)
- `eg`: Quadratic $f(x) = \frac{1}{2}x^TQx - b^Tx$: conditions for unique minimizer (Q23)
- `proposition` + `proof`: Saddle point characterization when Hessian has mixed eigenvalue signs (Q22)

---

### Lecture 2: Numerical Methods for Unconstrained Optimization

**Covers:** Bank Part C (Q27--Q42) + Part D (Q43--Q66)

**Narrative:** Now that we know *what* an optimum looks like, how do we *find* it? Progress from the simplest method (steepest descent) through Newton and quasi-Newton, with rigorous convergence analysis.

#### Subsections

**2.1 Descent Directions and Line Search Framework**
- `definition`: Descent direction (Q27)
- `proposition` + `proof`: $d_k$ is descent iff $\nabla f(x_k)^T d_k < 0$
- `notation`: Generic line search iteration $x_{k+1} = x_k + \alpha_k d_k$ (Q28)
- `definition`: Armijo condition (sufficient decrease) (Q29)
- `intuition`: Geometric interpretation of Armijo
- `definition`: Wolfe conditions (Armijo + curvature) (Q30)
- `remark`: Why Armijo alone is insufficient; role of $c_1$ (Q31)
- `definition`: Backtracking line search algorithm (Q32)
- `eg`: Exact line search on $(x-3)^2$ and $x^4 - 4x^2$ (Q33)

**2.2 Steepest Descent**
- `definition`: Steepest descent iteration $d_k = -\nabla f(x_k)$ (Q43)
- `proposition` + `proof`: $-\nabla f$ solves $\min_{\|d\|=1} \nabla f^T d$ (Q43)
- `theorem` + `proof`: Exact line search step size for quadratics: $\alpha_k = g_k^T g_k / (g_k^T Q g_k)$ (Q35)
- `eg`: Two iterations on $f(x_1,x_2) = 5x_1^2 + x_2^2$ showing zig-zagging (Q44)
- `proposition` + `proof`: Consecutive gradients are orthogonal under exact line search (Q47)
- `theorem` + `proof`: Convergence rate for quadratics: $\left(\frac{\kappa-1}{\kappa+1}\right)^2$ (Q46)
- `remark`: What happens as $\kappa \to 1$ vs $\kappa \to \infty$ (Q45)

**2.3 Global Convergence Theory**
- `definition`: Lipschitz continuity of $\nabla f$ (Q41)
- `theorem` + `proof`: Quadratic upper bound from Lipschitz gradient
- `theorem` + `proof`: Gradient descent with $\alpha = 1/L$ gives sufficient decrease (Q42)
- `theorem` (Zoutendijk) + `proof`: $\sum \cos^2\theta_k \|\nabla f(x_k)\|^2 < \infty$ (Q36)
- `corollary` + `proof`: If $\cos\theta_k \geq \delta > 0$, then $\liminf \|\nabla f(x_k)\| = 0$ (Q37)
- `remark`: Limitations of Zoutendijk -- no guarantee of local min, no rate, no non-smooth (Q38, Q39)

**2.4 Convergence Rate Analysis**
- `theorem` + `proof`: $O(1/k)$ sublinear rate for $L$-smooth convex functions (Q63)
- `theorem` + `proof`: Linear convergence for $m$-strongly convex: $(1 - m/L)^k$ rate (Q64)
- `remark`: Connection between condition number $\kappa = L/m$ and convergence speed

**2.5 Newton's Method**
- `intuition`: Approximate $f$ by second-order Taylor, minimize the quadratic model
- `definition`: Newton direction $d_k^N = -[\nabla^2 f(x_k)]^{-1} \nabla f(x_k)$ (Q49, Q50)
- `proposition`: Newton direction is descent when $\nabla^2 f(x_k) \succ 0$
- `theorem` + `proof`: Quadratic convergence near $x^*$ with $\nabla^2 f(x^*) \succ 0$ (Q52)
- `eg`: One Newton step on $f(x_1,x_2) = x_1^2 + 2x_2^2 + x_1x_2 - 4x_1$ (Q51)
- `proposition` + `proof`: Newton converges in one step for quadratics (Q56)
- `remark`: Disadvantages -- cost, far-from-minimum behavior, indefinite Hessian (Q53)
- `remark`: Can converge to saddle point; example with $f(x)=x^3$ (Q54)
- `definition`: Damped Newton (Newton + line search) (Q55)

**2.6 Quasi-Newton Methods (BFGS)**
- `intuition`: Motivation -- avoid Hessian computation (Q57)
- `definition`: Secant condition $B_{k+1}s_k = y_k$ (Q58)
- `definition`: BFGS update formula for $H_{k+1}$ (Q59)
- `remark`: First BFGS step with $H_0 = I$ is steepest descent (Q61)
- `theorem`: BFGS achieves superlinear convergence (Q60)
- `theorem` (Dennis-More condition): characterization of superlinear convergence (Q62)

**2.7 Comparison of Methods**
- Summary table: steepest descent vs Newton vs BFGS -- rate, condition number dependence, cost per iteration (Q65)
- `remark`: Self-concordance and Newton's domain of quadratic convergence (Q66, brief)

---

### Lecture 3: Constrained Optimization

**Covers:** NOT in the question bank -- entirely new content. Essential for MPC and general OC.

**Narrative:** Real problems have constraints. Build the theory of constrained optimization from Lagrangians through KKT, then connect to duality. This is the bridge between pure optimization and control.

#### Subsections

**3.1 Problem Formulation**
- `definition`: General constrained optimization problem (equality and inequality constraints)
- `definition`: Feasible set, active constraints, inactive constraints
- `eg`: Simple 2D constrained problems for geometric intuition

**3.2 Equality-Constrained Optimization and Lagrange Multipliers**
- `intuition`: Why unconstrained FONC fails -- gradient need not be zero, only zero in feasible directions
- `definition`: Lagrangian $\mathcal{L}(x, \lambda) = f(x) + \lambda^T h(x)$
- `theorem` + `proof`: Lagrange multiplier theorem (first-order necessary conditions for equality constraints)
- `intuition`: Geometric interpretation -- gradient of objective parallel to gradient of constraint
- `eg`: Minimize $x_1^2 + x_2^2$ subject to $x_1 + x_2 = 1$
- `definition`: Constraint qualification (LICQ)
- `remark`: What goes wrong without constraint qualification

**3.3 Inequality Constraints and KKT Conditions**
- `definition`: Lagrangian with inequality constraints
- `theorem` (KKT Necessary Conditions) + `proof`: Stationarity, primal feasibility, dual feasibility, complementary slackness
- `intuition`: Complementary slackness -- either a constraint is active or its multiplier is zero
- `eg`: Minimize $f(x) = (x_1 - 2)^2 + (x_2 - 1)^2$ subject to $x_1 + x_2 \leq 1$, $x_1, x_2 \geq 0$
- `theorem` (KKT Sufficient Conditions): KKT + convexity implies global optimality

**3.4 Second-Order Conditions for Constrained Problems**
- `definition`: Critical cone / tangent cone at a KKT point
- `theorem`: SONC -- Hessian of Lagrangian PSD on critical cone
- `theorem`: SOSC -- Hessian of Lagrangian PD on critical cone

**3.5 Lagrangian Duality**
- `definition`: Lagrange dual function
- `proposition` + `proof`: Dual function is concave (even for non-convex primal)
- `definition`: Dual problem
- `theorem` + `proof`: Weak duality
- `theorem`: Strong duality (Slater's condition for convex problems)
- `remark`: Duality gap and its interpretation

**3.6 Quadratic Programming**
- `definition`: QP formulation
- `proposition`: KKT conditions for QP are a linear system (equality-constrained case)
- `remark`: QP as the workhorse of MPC
- `eg`: Solve a small QP by hand using KKT

---

### Lecture 4: Dynamic Programming and the Calculus of Optimal Control

**Covers:** New content (HJB, Pontryagin) + motivation for Bank Parts E, F

**Narrative:** Shift from static optimization to sequential decision-making. Two fundamental approaches: Bellman's DP (backward, sufficient) and Pontryagin's PMP (forward-backward, necessary). These are the twin pillars on which everything else rests.

#### Subsections

**4.1 Optimal Control Problem Formulation**
- `definition`: Continuous-time OCP (Bolza form)
- `definition`: Discrete-time OCP
- `notation`: Bolza, Mayer, and Lagrange forms; equivalence between them
- `remark`: Connection to robotics -- trajectory optimization, energy minimization

**4.2 Dynamic Programming (Discrete-Time)**
- `intuition`: Principle of optimality -- "a tail of an optimal trajectory is itself optimal"
- `theorem` (Bellman's Principle of Optimality) + `proof`
- `definition`: Value function
- `theorem`: Bellman equation (backward recursion)
- `eg`: Simple discrete-time example
- `remark`: Curse of dimensionality

**4.3 Hamilton-Jacobi-Bellman Equation (Continuous-Time)**
- `definition`: Continuous-time value function
- `theorem` + `proof` (sketch): HJB PDE
- `remark`: HJB is sufficient -- if you can solve it, you have the global optimum
- `remark`: Generally intractable for nonlinear systems (curse of dimensionality)
- Connection to LQR: HJB yields the Riccati equation for linear-quadratic problems (preview)

**4.4 Pontryagin's Minimum Principle**
- `definition`: Hamiltonian $H(x, u, p) = L(x, u) + p^T f(x, u)$
- `definition`: Costate (adjoint) variable $p(t)$
- `theorem` (PMP): State equation, costate equation, optimality condition, transversality
- `intuition`: PMP as infinite-dimensional KKT -- costates are Lagrange multipliers for dynamics constraints
- `proof` (outline): Derivation via needle variation argument
- `eg`: Minimum-energy control of $\dot{x} = u$
- `eg`: Bang-bang control: minimize time with bounded control
- `remark`: PMP is necessary, not sufficient (contrast with HJB)
- `remark`: Connection between PMP Hamiltonian and HJB -- when $V$ is smooth, $p = \nabla_x V$

**4.5 Calculus of Variations (Brief)**
- `definition`: Euler-Lagrange equation
- `proof`: Derivation via first variation
- `remark`: PMP generalizes Euler-Lagrange to handle control constraints

---

### Lecture 5: Numerical Methods for Trajectory Optimization

**Covers:** New content -- not in the question bank. Bridges Pontryagin/HJB theory and practical computation.

**Narrative:** PMP gives necessary conditions as a boundary value problem; HJB gives a PDE. Neither is directly computable for general systems. This lecture covers the numerical methods that actually solve trajectory optimization problems.

#### Subsections

**5.1 The Two-Point Boundary Value Problem from PMP**
- `remark`: PMP yields a TPBVP -- state equation forward, costate equation backward, coupled through optimality condition
- `definition`: The TPBVP structure
- `remark`: Why this is hard -- split boundary conditions, nonlinear coupling

**5.2 Indirect Shooting Methods**
- `intuition`: Guess the initial costate, integrate forward, check terminal conditions, iterate
- `definition`: Single shooting for TPBVP
- `remark`: Reduces TPBVP to root-finding, solved via Newton's method
- `remark`: Sensitivity -- small changes in $p(t_0)$ cause large changes at $t_f$
- `eg`: Indirect shooting for minimum-energy control of a double integrator

**5.3 Direct Single Shooting**
- `intuition`: Parameterize controls and optimize directly
- `definition`: Direct single shooting -- discretize controls, simulate forward, optimize
- `remark`: States are implicit functions of controls -- no state decision variables
- `remark`: Advantages and disadvantages

**5.4 Direct Multiple Shooting**
- `intuition`: Break trajectory into segments, shoot independently, stitch with continuity constraints
- `definition`: Multiple shooting -- segment-wise decision variables, continuity constraints
- `remark`: Better conditioning than single shooting
- `remark`: State constraints easy to add at segment boundaries
- `proposition`: Multiple shooting produces a sparse NLP (block structure)
- `remark`: Used by ACADOS, CasADi

**5.5 Direct Collocation**
- `intuition`: Approximate trajectory with polynomials, enforce dynamics at collocation points
- `definition`: Trapezoidal collocation
- `definition`: Hermite-Simpson collocation
- `remark`: States and controls are both decision variables -- dynamics become equality constraints
- `proposition`: Collocation produces a large but very sparse NLP
- `remark`: Handles stiff systems better than explicit methods

**5.6 Transcription to Nonlinear Programming (NLP)**
- `definition`: The general NLP from direct transcription
- `remark`: Connection to Lecture 3 -- KKT conditions, solved by SQP or interior point
- `remark`: Sparsity structure -- banded Jacobians/Hessians
- `remark`: Warm-starting and real-time applicability

**5.7 Comparison of Methods**
- Summary table: indirect shooting vs direct single shooting vs multiple shooting vs collocation
- `remark`: In robotics practice, direct collocation and multiple shooting dominate
- `remark`: Connection to iLQR (Lecture 8) -- iLQR is a specialized solver for unconstrained direct trajectory optimization

---

### Lecture 6: Linear Quadratic Regulator

**Covers:** Bank Part E (Q67--Q82) + Part F (Q83--Q92)

**Narrative:** LQR is *the* solvable optimal control problem. HJB gives Riccati, PMP gives the same Riccati, and we get a closed-form linear feedback law.

#### Subsections

**6.1 Problem Statement and Structure**
- `definition`: Finite-horizon discrete-time LQR (Q67)
- `notation`: Assumptions on $Q \succeq 0$, $R \succ 0$, $Q_f \succeq 0$
- `intuition`: Why $Q \succeq 0$ and $R \succ 0$ (Q82); what goes wrong if $R$ is only PSD (Q77)
- `remark`: Role of $Q$ vs $R$ -- state penalty vs control effort tradeoff (Q71)

**6.2 Finite-Horizon LQR via Dynamic Programming**
- `theorem` + `proof`: Derivation of optimal control $u_k^* = -K_k x_k$ via Bellman equation (Q69)
- `theorem`: Discrete-time Riccati recursion (Q70)
- `definition`: Gain matrix $K_k$
- `eg`: Scalar LQR solved by hand (Q72)

**6.3 Infinite-Horizon LQR and the Algebraic Riccati Equation**
- `definition`: Infinite-horizon LQR problem (Q73)
- `remark`: Existence requires stabilizability/detectability (covered in linear controls notes)
- `theorem`: Conditions for unique stabilizing solution of DARE (Q74, Q84)
- `theorem` + `proof` (outline): Closed-loop stability (Q75)
- `eg`: 2D system setup of ARE (Q81)

**6.4 The Riccati Equation -- Properties and Structure**
- `definition`: DARE vs CARE (Q83)
- `theorem` + `proof` (outline): Monotonicity and convergence (Q85)
- `proposition` + `proof`: Symmetry of $P_k$ by induction (Q87)
- `proposition`: Invertibility of $(R + B^T P B)$ (Q88)
- `theorem`: Hamiltonian matrix connection (Q86)
- `remark`: Schur complement interpretation (Q89)
- `eg`: Scalar DARE (Q90)
- `remark`: Limiting behavior as $R \to 0^+$ and $R \to \infty$ (Q91)

**6.5 Continuous-Time LQR**
- `definition`: CT-LQR problem, DRE, ARE (Q76)
- `eg`: Scalar CT-LQR (Q78)
- `theorem`: Robustness margins -- gain $\geq$ 6 dB, phase $\geq$ 60 deg (Q80)

**6.6 LQR and Pontryagin's Minimum Principle**
- Show PMP yields the same Riccati equation
- `remark`: Two roads to the same destination -- DP vs PMP for LQR

---

### Lecture 7: Estimation -- The Kalman Filter

**Covers:** Bank Part H (Q104--Q118)

**Narrative:** The Kalman filter is the dual of LQR -- same Riccati equation, different roles.

#### Subsections

**7.1 Problem Setup**
- `definition`: Linear stochastic system (Q104)
- `notation`: Noise assumptions
- `definition`: State estimation problem

**7.2 The Kalman Filter Algorithm**
- `definition`: Prediction step (Q105)
- `definition`: Correction step (Q105)
- `intuition`: Predict from dynamics, correct with measurements
- `theorem`: Optimality (Q106)

**7.3 Derivation of the Kalman Gain**
- `theorem` + `proof`: Kalman gain by minimizing $\tr(P_{k|k})$ (Q107)

**7.4 Error Covariance and Riccati Recursion**
- `theorem`: Riccati recursion for error covariance (Q108)
- `remark`: Connection to LQR Riccati
- `proposition`: Limiting behavior of Kalman gain (Q110, Q111)

**7.5 Steady-State Kalman Filter**
- `theorem`: Convergence conditions (Q112, Q114)
- `definition`: Steady-state KF (Q113)
- `remark`: Steady-state equations are a DARE in disguise (Q113)

**7.6 Duality Between LQR and Kalman Filter**
- `theorem` + `proof`: The duality mapping (Q115, Q92)
- `intuition`: Solving one problem immediately solves the other

**7.7 Innovation Sequence**
- `definition`: Innovation (Q116)
- `proposition` + `proof`: Innovations are white noise for correctly tuned filter

**7.8 Extended Kalman Filter**
- `definition`: EKF (Q117)
- `remark`: Limitations
- `remark`: Information filter (Q118)

**7.9 Scalar Kalman Filter Example**
- `eg`: Full worked example (Q109)

---

### Lecture 8: LQG Control, Trajectory Optimization, and iLQR

**Covers:** Bank Part I (Q119--Q128) + Part G (Q93--Q103)

**Narrative:** Combine estimation and control (LQG), then handle nonlinear systems via trajectory optimization. iLQR is the practical algorithm that makes nonlinear optimal control tractable for robotics.

#### Subsections

**8.1 Linear Quadratic Gaussian (LQG) Control**
- `definition`: LQG problem (Q119)
- `theorem` (Separation Principle) + `proof` (outline) (Q120, Q121)
- `remark`: Block diagram (Q122)
- `eg`: Scalar LQG design (Q126)

**8.2 Robustness of LQG**
- `theorem` (Doyle's counterexample): LQG does NOT inherit LQR margins (Q123)
- `intuition`: Why the observer destroys robustness
- `definition`: Loop Transfer Recovery (LTR) (Q124)
- `remark`: LQG as $\mathcal{H}_2$ optimal control (Q125)
- `remark`: Limitations and extensions (Q127, Q128)

**8.3 Nonlinear Trajectory Optimization**
- `definition`: General nonlinear trajectory optimization problem
- `remark`: Why LQR/LQG are insufficient
- `intuition`: Linearize, solve LQR subproblem, re-linearize, iterate

**8.4 Iterative LQR (iLQR)**
- `definition`: What iLQR solves vs standard LQR (Q93)
- Algorithm description: forward pass, backward pass, gains (Q94)
- `definition`: Linearization and quadratization about nominal trajectory (Q96)
- `proposition`: iLQR as approximate Newton on trajectory space (Q95)
- `remark`: No guarantee of global optimality (Q97)

**8.5 iLQR: Practical Considerations**
- `definition`: Regularization (Q98)
- `remark`: Line search (Q101)
- `remark`: Computational complexity (Q102)

**8.6 Differential Dynamic Programming (DDP)**
- `definition`: DDP vs iLQR (Q99)
- `remark`: iLQR is the Gauss-Newton approximation to DDP's full Newton

**8.7 Handling Constraints in Trajectory Optimization**
- `remark`: iLQR/DDP do not natively handle constraints (Q103)
- `definition`: Augmented Lagrangian iLQR (AL-iLQR)
- `remark`: Barrier methods, projected methods

**8.8 Worked iLQR Example**
- `eg`: 1D nonlinear system, one full iteration (Q100)

---

### Lecture 9: Model Predictive Control

**Covers:** NOT in question bank -- entirely new content. Critical for modern robotics.

**Narrative:** MPC bridges optimization and real-time control. It repeatedly solves finite-horizon optimal control problems online, naturally handling constraints. This is where constrained optimization (Lecture 3) meets LQR (Lecture 6).

#### Subsections

**9.1 The MPC Concept**
- `intuition`: Receding horizon
- `definition`: MPC problem formulation
- `remark`: Why MPC vs LQR -- LQR cannot handle constraints (Q79)
- `definition`: Receding horizon strategy
- `remark`: MPC reduces to LQR when unconstrained and $P$ solves DARE

**9.2 Linear MPC and QP Formulation**
- `definition`: Linear MPC
- `theorem`: Condensing MPC into a QP (eliminate states, stack, substitute)
- `eg`: Condensed QP for a 2-step MPC problem
- `remark`: Sparse vs dense QP formulations

**9.3 Stability of MPC**
- `intuition`: Why stability is not automatic
- `definition`: Terminal cost and terminal constraint
- `theorem` + `proof` (sketch): MPC stability via Lyapunov argument
- `remark`: Terminal set, unconstrained tail = LQR

**9.4 Feasibility and Recursive Feasibility**
- `definition`: Recursive feasibility
- `proposition`: Terminal constraint ensures recursive feasibility
- `remark`: Soft constraints as practical relaxation

**9.5 Nonlinear MPC**
- `definition`: NMPC
- `remark`: Computational challenge -- NLP at each step
- `remark`: Connection to iLQR
- `remark`: Real-time iteration scheme

**9.6 MPC in Practice**
- `remark`: QP solvers, timing constraints
- `remark`: Warm-starting
- `remark`: Explicit MPC

---

## Overall Narrative Flow

```
Lecture 1: Convexity + Unconstrained Optimality
    "What does an optimum look like?"
         |
         v
Lecture 2: Numerical Methods (Steepest Descent, Newton, BFGS)
    "How do we find an optimum?"
         |
         v
Lecture 3: Constrained Optimization (Lagrangians, KKT, QP)
    "What if there are constraints?"
         |
         v
Lecture 4: Dynamic Programming, HJB, Pontryagin
    "What if the problem evolves over time?"
         |
         v
Lecture 5: Shooting Methods, Collocation, Direct Transcription
    "How do we numerically solve trajectory optimization?"
         |
         v
Lecture 6: LQR and Riccati Equations
    "The one problem we can solve exactly"
         |
         v
Lecture 7: Kalman Filter
    "What if we can't measure the state?"
         |
         v
Lecture 8: LQG + iLQR/DDP
    "Combining estimation & control; handling nonlinearity"
         |
         v
Lecture 9: Model Predictive Control
    "Real-time constrained optimal control"
```

---

## Execution Order

1. **Setup**: Copy template from `Note/`, update `master.tex` with OC-specific commands
2. **Lecture 1**: Convexity + Unconstrained Optimality
3. **Lecture 2**: Numerical Methods
4. **Lecture 3**: Constrained Optimization
5. **Lecture 4**: DP, HJB, Pontryagin
6. **Lecture 5**: Shooting Methods + Collocation
7. **Lecture 6**: LQR + Riccati
8. **Lecture 7**: Kalman Filter
9. **Lecture 8**: LQG + iLQR/DDP
10. **Lecture 9**: MPC
11. **Final**: Full compilation check
