# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

`0.1.0` has not shipped, so all of this may be folded into it; it is kept separate
because the KdV sign convention below changes what every number in the package means.

### Added — `verify_zeitlin_three_bracket.jl`

Whether the Bialynicki-Birula--Morrison three-bracket `[A, B, S] = {A, B}` reproduces the Zeitlin
bracket on `su(N)`. It does, at `N = 3, 5`, to `2 × 10⁻¹⁵`.

The question is not routine, which is why the script exists rather than a remark. Salmon (2005)
§5 reports that this mechanism fails for the two-dimensional vorticity bracket, the continuum
algebra carrying no non-degenerate invariant form and hence no quadratic Killing Casimir.
**The sine-bracket truncation supplies exactly that**: the Killing form of `sine_algebra(N)` is
non-degenerate, and `κ^{mn} ∝ δ_{m+n,0}` to `2 × 10⁻¹⁵`, so Zeitlin's `I₂ = Σ_k ω_k ω_{-k}` *is*
the Killing Casimir rather than merely something of the same shape. The mechanism fails before
truncation and works after it.

**What semi-simplicity buys is non-degeneracy of `κ`, and nothing else here.** The lowered
structure constants `c_ijk = c_ij^m κ_mk` are totally antisymmetric for *every* Lie algebra, by
ad-invariance of the Killing form; the script measures a `(j,k)` defect of exactly zero on
`se(3)` and on `so(3)` with abelian directions, both of which have a singular `κ`. That check is
kept as a verification of the index conventions, not as evidence for the hypothesis.

Two negative controls, and they are half the point of the script:

- `I₃`, the cubic invariant, in the third slot gives `1.019` relative deviation after optimal
  rescaling — the quadratic hypothesis is load-bearing, and the reason is checked directly as a
  homogeneity mismatch, degree 1 against degree 2 in `ω`.
- **`se(3)`**, where `κ` is singular *and* its degenerate directions are non-central, so even
  `pinv(κ)` — the most generous repair available — recovers nothing: `max|B3 - J| = 2.175`
  against `max|J| = 2.175`. `so3(Float64, n)` for `n > 3` is run alongside it and is a **false**
  control: `κ` is singular there too, but the degenerate directions are central, `J` vanishes on
  them anyway, and `pinv` recovers the bracket exactly. Same trap the `so3` docstring already
  warns about for Jacobi checks, in a second guise — a control that fails for a reason that does
  not generalise is no better than one that cannot fail.

Everything is computed from the structure constants alone — the Killing form as
`κ_mn = Σ_ij c_mi^j c_nj^i` — so nothing depends on the matrix realisation `sine_algebra`
happens to use. No new exports and no new dependencies; the script uses `sine_algebra`,
`lie_poisson_matrix`, `lie_poisson_derivative`, `jacobi_residual`, `so3` and `se3` as they stand.

### Added — an adversarial re-check of the structure-function bracket

`scripts/fable/III1_refutation.jl` attacks the same statement as
`scripts/fable/III1_structure_function_4d.jl` — that `{A,B}_c = ∫ c(f){A_f,B_f} μ` is Poisson for
every smooth `c` — from an independent implementation that shares no code with it or with
`fabletools.jl`, so that a defect in the shared machinery cannot hide in both. It computes the
Jacobiator by carrying each field's first variation as a linear differential operator and
integrating that operator by parts, rather than by directional derivatives and cyclicity, and
compares the result against a closed form derived separately. 99 checks, exact over `ℚ(i)`.

The statement survives, in dimensions 2, 4, 5, 6 and 8, for degenerate bivectors and on an
odd-dimensional torus with no symplectic form at all. Two things it establishes that the earlier
script does not: the identity `Q μ = α∧β∧γ∧δ∧ω^{n−2}/(n−2)!` holds with sign `+1` **only** against
the Liouville form, since `ω^n/n!` and `dx_1∧…∧dv_n` differ by `(−1)^{n(n−1)/2}`; and the Jacobiator
is non-zero — `−127697733/64` — as soon as the volume form is not Liouville, so that hypothesis is
load-bearing rather than decorative. Every vanishing verdict is guarded by an assertion that the
integrand could have been non-zero, after three separate accidental zeros during development.

### Added — an adversarial re-check of the matrix structure-function bracket

`scripts/fable/III2_refutation.jl` attacks the same two statements as
`scripts/fable/III2_matrix_structure_function.jl` — that `tr(c(W)[F_W,G_W])` on `gl(N)*`/`u(N)*` is
Poisson iff `c` is affine, and the spectral-kernel classification behind the pushforward family —
from an implementation that shares no code with it. 118 checks, exact over `ℚ`.

The point of the script is one step that no closed form can cross-check. A spectral-kernel bracket
is defined through the eigendecomposition of `V`, so its Jacobiator needs the derivative of the
**eigenvectors**; the classification theorem rests entirely on that step. Here the gradient is
instead built by brute force from first-order perturbation theory at a diagonal `V` with distinct
rational eigenvalues, which is rational in the entries, and the theorem's closed form is then a
prediction compared against it. The perturbation machinery is itself validated first against a route
that never mentions an eigenvector — the polynomial Fréchet derivative
`Dc_W(E) = Σ_k a_k Σ_{i+j=k−1} W^i E W^j` — which agrees with it exactly.

Both statements survive. The closed form matches brute force on all twenty kernel/dimension pairs
tried, including three kernels that are neither `c`-brackets nor pushforwards, where the two special
families give no cross-check at all; five pushforwards give exactly zero and five non-pushforwards
exactly non-zero. `N = 2` is confirmed to vanish for every `c`, and the script's own derivation shows
why without Cayley–Hamilton: the obstruction needs three *distinct* eigenvalue indices, so `N = 2` is
not a weak control but no control. The log-mean bracket is Poisson to `8e-15` relative through the
same code path that returns an exact `176` for `c = v²`, and the two kernels are shown to separate at
`O(gap²)`, which is the reported `N⁻²` rate.

One thing the script records that the earlier one does not: `c' > 0` restricts `c = v²` to `v > 0`,
the flat coordinate `ψ = ½ln v` does not cross zero, and the log-mean kernel `2(a−b)/log(a/b)` is
real only when `a` and `b` share a sign. The pushforward representative of the enstrophy bracket
therefore exists on definite spectra only — which a Zeitlin truncation of 2-D Euler vorticity is not.

### Added — shared verification machinery for the fable branch

`scripts/fable/fabletools.jl` provides exact-arithmetic machinery for the fable verification suite: trigonometric polynomials over `ℚ(i)` on an arbitrary-dimensional torus T^D, first-order jets for computing Gateaux derivatives, and the canonical Poisson bracket on T^{2n}. Plain definitions in a shared include file, so scripts extending the machinery to their own field types (as III.1 does for floating-point spectral grids) do so without module boundaries. Used by `III1_structure_function_4d.jl` and `III4_multifield_derivative.jl`.

### Added — Vlasov–Poisson on a bounded domain as Lie-Poisson on `u(N)*`

`scripts/fable/I1_vlasov_uN.jl` quantises the Vlasov–Poisson system on a bounded spatial domain over the space of density matrices ρ on u(N)*, using compressed Weyl quantisation. The Lie-Poisson bracket `{F,G}(ρ) = (1/iℏ) tr(ρ [F_ρ, G_ρ])` is the underlying structure. The main theorem: for symbols of degree ≤ 1 in momentum on the inner block (away from boundaries), the homomorphism defect `D(a,b) = (1/iℏ)[T_N a, T_N b] − T_N{a,b}` vanishes exactly, and for degree 2 it equals the Wigner-equation correction with its exact coefficient `(ℏ²/24) φ''' ∂_p³`. 89 checks, including exact arithmetic on gl(N,ℚ) for N=3,4,5, Casimirs and continuity in dimensions and boundary conditions (periodic S¹ and Dirichlet [0,π]), and floating-point refinement studies showing convergence at order two. Not yet independently re-checked: the adversarial pass has covered III.1 and III.2 only.

### Added — semidirect-product brackets via pointwise flat coordinates

`scripts/fable/II1_semidirect_flat.jl` verifies that the push-forward of the constant bracket K₂ = [0 K; K 0] under the pointwise chart Φ(ρ,u) = (ρ, ρu) produces a Poisson bracket on the semidirect product — the 1-D compressible Euler and shallow-water systems. The construction `J = DΦ K₂ DΦ^T` generalises to any antisymmetric K and is exactly Poisson by the Jacobi theorem for pushforwards. 147 checks verify: Jacobi over ℚ exactly for random antisymmetric K and rational data; Casimirs and rank; consistency of the three block structure against SymPy; refinement convergence; and boundary conditions on [0,L] with two antisymmetrisation conventions, both exactly Poisson.

### Added — the structure-function bracket on `T^{2n}`

`scripts/fable/III1_structure_function_4d.jl` proves that the bracket `{A,B}_c = ∫_M c(f) {A_f, B_f} μ` satisfies the Jacobi identity for every smooth `c` and on any closed symplectic manifold (M, ω) of dimension 2n, ruling out the provisional hypothesis that the two-dimensional Plücker relation was essential. The reduction of the Jacobiator to the obstruction term survives dimension unchanged (cyclicity, skew-adjointness, chain rule); the integrand is not zero in 4-D and 6-D but is an exact divergence. 58 checks, exact over ℚ(i), verify: the obstruction on T², T⁴, T⁶ for polynomial `c`; the divergence identity term-by-term; the full Jacobiator by jets for positive cases and negative controls that must fail; and floating-point spectral refinement N = 12…32 for non-polynomial `c`.

### Added — the matrix structure-function bracket on `gl(N)`

`scripts/fable/III2_matrix_structure_function.jl` establishes the finite-dimensional analogue of III.1: the bracket `tr(c(W)[F_W, G_W])` on gl(N) with matrix spectral functions `c`. The Jacobi identity holds exactly iff `c` is affine, or equivalently iff the antisymmetric edge weight `κ_ab (w_b − w_a)` in the eigenbasis Jacobiator forms a cocycle on the complete graph of eigenvalues. The spectral-kernel classification then says which brackets in that wider class *are* Poisson: exactly the pushforwards of Lie-Poisson under a monotone spectral `ψ`, with slope `ϖ = 1/ψ^{[1]}` — an infinite family, one member per `ψ`, of which the log-mean bracket (`ψ = ½ln v`) is the worked example. Restricted to symmetric polynomial slopes of degree at most two the survivors are only `ϖ = λ₁` and `ϖ = λ₂ ab`. 114 checks over ℚ; semiclassical scaling `N^{−2}` is measured on Zeitlin's su(N) for four spectral functions.

### Added — multi-field structure functions, and derivative dependence

`scripts/fable/III4_multifield_derivative.jl` settles the two remaining open questions from the four-bracket paper: (a) the bracket `{A,B} = ∫ c^{ij}(f){A_i, B_j} μ` for several fields f¹…f^m is Poisson iff the product `(ξ⋆η)_k = ∂_k c^{ij} ξ_i η_j` on T*ℝ^m is associative and (in dimensions 2n ≥ 4) the target one-forms close, and (b) structure functions depending on ∇f never produce Poisson brackets — the principal symbol forces all first-order derivatives to vanish pointwise. 77 checks prove both theorems exactly over ℚ(i) by jets, verify the necessity lemma (the three Jacobi functionals obey exactly one relation), test positive and negative cases of multi-field coupling, and rule out derivative dependence in 2-D and 4-D.

### Added — the four-bracket manuscript

`poisson-brackets-from-four-brackets.tex` arrived with a standalone Julia verification suite of
its own: a `BracketChecks.jl` carrying spectral and 8th-order finite-difference differentiation
on the periodic two-torus, a `CheckSet` harness, and four drivers. It is now part of the package,
on the same terms as the two manuscripts converted from Python in `f897fb0` — shared machinery in
`src/`, diagnoses in `scripts/`.

Three new source files, and no new dependencies:

| file | contains |
|:--|:--|
| `src/torus.jl` | `TorusGrid`, `spectral_grid`, `finite_difference_grid`, `∂x`, `∂y`, `sample`, `integrate`, `canonical_bracket` |
| `src/fourbrackets.jl` | the Gardner and symmetric operators, the Gardner/symmetric/weighted two- and four-brackets with their densities, `lie_poisson_2bracket`, `antisymmetry_residuals`, `plucker_residual`, and a `jacobi_residual` method for the continuum obstruction of Theorem 4.5 |
| `src/metriplectic.jl` | `kulkarni_nomizu` and `metriplectic_bracket`, Section 2 |

`TorusGrid` stores a differentiation **matrix** rather than the original's `Dx::Function` /
`Dy::Function` closures: a derivative is `D * f` or `f * transpose(D)`, so the two schemes share
one concrete type, the spectral derivative costs `O(N³)` through BLAS instead of `O(N⁴)` per call,
and the matrix is antisymmetric — exactly, for the centred stencil — which is where the
antisymmetry of `canonical_bracket` now comes from.

The refactoring exposed one structural fact the original's copies obscured: a two-bracket **is**
its four-bracket with the entropy in the second and fourth slots. `gardner_2bracket` is defined
as `gardner_4bracket(g, a, s, b, s)` and `symmetric_2bracket` likewise, and the weighted brackets
are the symmetric densities times the weight. The tests assert both identifications with `==`.

Four new scripts — `verify_fourbracket_metriplectic.jl`, `_identities.jl`, `_convergence.jl` and
`_log_entropy.jl` — plus `scripts/torustools.jl` for the fixed test fields, which the drivers
carried in three copies. 124 new tests in `test/torus_tests.jl`, `test/fourbrackets_tests.jl` and
`test/metriplectic_tests.jl`.

Every finite-difference residual and observed rate reproduced the original suite to all printed
digits. The spectral residuals moved within `1e-14`–`1e-16` and uniformly downward, the matrix
form accumulating less round-off than a per-call DFT. `docs/src/verification.md` records the
comparison, three defects the fold-in turned up, and the one claim that was added.

### Changed — `scripts/check.jl` learned the convergence idiom

`check_exact` and `check_refined` join `check`, for claims about a continuum identity rather than
a matrix: one settles against an absolute threshold, the other against the factor by which the
residual falls when the resolution is doubled. `check_refined` prints
`3.14e-08 -> 1.25e-10   251x`, because against a predicted 256 the factor *is* the evidence.
`relerr` and `normerr` came along with them. The twenty existing scripts are untouched, and the
new functions know nothing about grids — resolutions are the calling manuscript's business.

`CheckSet`, `check!`, `check_pointwise!` and `report` are gone; the LaTeX labels their table
carried in a column (`eq:cyclicity`, `thm:jacobi`) are now in the check labels, so a line of
output still names the equation it settles.

### Fixed — a relative `[sources]` path does not compose

`scripts/Project.toml` now carries its own `SimpleSolvers` entry. A relative source path is
resolved against the *active* project, so the package's own `../SimpleSolvers` became
`PoissonBrackets/SimpleSolvers` whenever `scripts/` was the active environment and resolution
failed outright — which left the scripts manifest stranded without `Sparspak` and every script
unable to load the package. The entry uses `path` rather than `{rev, url}`: the url form resolves
against a cached clone, which had pinned that environment to 0.12.2 while the working tree was
already at the 0.13 the package requires. `docs/Project.toml` will need the same when it next
resolves.

### The mixed two-field formulation

`Integrator(...; formulation = :mixed)` solves an implicit-midpoint step as `2N` equations in
`(y, v)` instead of `N` in `y`:

```
M (y - uⁿ) - Δt K(ū) v = 0        M v - ∂H/∂u(ū) = 0        ū = (uⁿ + y)/2
```

Eliminating `v` recovers the dense residual exactly, so it is the same step — verified to 12–14
digits over single steps and 100-step runs, with the conserved quantities holding to the same
round-off as before. What it avoids is `Minv * K * Minv`: **the dense Newton matrix is dense
only because the inverse mass matrix is**, and every dynamical use of `Minv` was a mass solve
written as a dense multiply. All four blocks of the mixed Jacobian are banded with circular
bandwidth `p`, and the pattern is fixed for the whole run, so one ordering and symbolic
factorisation serve every step.

One step of implicit midpoint on the second KdV flow, `p = 3`, in milliseconds, both
formulations converging in two Newton iterations:

| N | dense | mixed | speedup |
|---:|---:|---:|---:|
| 64 | 0.192 | 0.259 | 0.74 |
| 128 | 0.518 | 0.485 | 1.07 |
| 384 | 3.34 | 1.56 | **2.1** |
| 1024 | 22.9 | 4.22 | **5.4** |
| 1536 | 58.2 | 6.29 | **9.3** |

The crossover is near `N = 125`, so `:dense` remains the default and the small cases are
unaffected. The gain grows because the dense factorisation is `O(N³)` and the assembly's two
dense `N × N` products are too, where the mixed form is `O(N)` in both.

Restricted to `ImplicitMidpoint` and to an `AffineBracket`; every other combination is an
`ArgumentError` rather than a silently different method. `AverageVectorField` would need an
auxiliary per quadrature node, the discrete-gradient methods carry a rank-one term that would
have to be bordered, and the remaining brackets store the sandwiched `Minv*K*Minv` rather than
`K`.

### `K0` is no longer densified

`AffineBracket`'s constructor called `Matrix(_skew(K0))`, throwing away the sparsity of a
matrix that is banded with circular bandwidth `p` — `mixed_matrix(s, 1, 2)` for KdV, zero for
Camassa-Holm. It is now kept as given, exactly as `Ψ` already was, which is what makes
`kernel_operator` banded and hence the mixed formulation possible at all. It also makes
`poisson_apply`'s `K0 * v` an `O(Np)` matvec instead of `O(N²)`, which the explicit runs pay
four times a step for millions of steps.

`Minv` stays dense, and the analysis paths that genuinely need it — `poisson_matrix`,
`poisson_tensor`, `kernel_tensor`, the Jacobi-identity and Casimir diagnostics — are untouched.

### `kernel_operator`

New: `kernel_operator(b::AffineBracket, û)` is `K0 + kernel_matrix(b, û)`, the weak-form block
without the surrounding inverse mass matrices — `poisson_matrix`'s middle factor, and sparse
where `poisson_matrix` is dense. Together with `kernel_directional`, which the bracket already
provided and which is already `∂(K v)/∂û`, it is everything the mixed Jacobian needs.

### The linear solver for the mixed system is a correctness constraint

`SparspakLU`, not the `UmfpackLU` that SimpleSolvers would otherwise select for a sparse
`Float64` Jacobian. UMFPACK mis-handles this block structure: from `N = 768` upward it returns
a solution wrong by a factor of 150 while reporting success, which surfaces as Newton diverging
from a starting point whose residual was already `1e-5`. Sparspak and dense LAPACK both solve
the same matrices to the accuracy their condition number allows. `Sparspak` is therefore a new
dependency — pure Julia, and light.

### Requires SimpleSolvers 0.13

For the sparse-Jacobian plumbing (`jacobian_prototype`) and the sparse linear solvers.

### Changed — the KdV sign convention is now the textbook one

The equation is `u_t + 6uu_x + u_xxx = 0`, so the `sech²` solitons are **elevations**. The
older `u_t = 6uu_x - u_xxx` is reached by `u → -u`.

The one step that does not follow from a global sign flip is the second structure:
`4u∂ₓ + 2u_x - ∂ₓ³` becomes `-4u∂ₓ - 2u_x - ∂ₓ³` — the transport part flips, the third
derivative does not. That was derived and checked symbolically in jet variables
(`verify_kdv_continuous.py`) *before* any code or prose was touched; in the package it is one
character, `AffineBracket`'s `scale` from `2` to `-2`, `K0` being the `-∂ₓ³` block and keeping
its own sign.

- `H₁ = ½∫(u_x² - 2u³)`, `δH₁/δu = -3u² - u_xx`; `H₂` unchanged;
  `H₃ = ∫(5/2 u⁴ - 5uu_x² + ½u_xx²)`.
- The Magri rung `P²g = 2P¹ ∂H₂/∂û` becomes `-2`, for the same reason.
- `soliton` and `two_soliton` return elevations; the Shi-Fu-Liu correspondence maps by
  `u → u/6` rather than `u → -u/6`, and their mass is `6C₀` rather than `-6C₀`.
- Each Hamiltonian carries the matching sign, so `H₁`, `H₂` and `H₃` take the *same numerical
  values* on corresponding data as before — the elevation soliton here has the `H₁` the
  depression soliton had there, `-6.39998180` at `N = 128`, agreeing with
  `plot_kdv_energies.py` to eight digits. Only `C₀` changes, and it changes sign.
- The Miura map becomes `u = -(v² + v_x)` and its Jacobian `-(2v + ∂ₓ)`. The overall sign
  cancels in `L P¹ Lᵀ`, so the bracket is unchanged; the chart is not, and the Hill potential
  becomes `-u_h`. The Hill threshold for `u = c + sin x + 0.4cos2x` at `N = 20`, `p = 3` moves
  from `+0.470894` to `-0.328741`, the admissible side becoming the negative one.

### Added — a spectral parameter for the Miura chart

`miura_map(space, v̂, λ)` implements `u = -(v² + v_x) - λ`, with `miura_invert(…; λ)` and the
new `miura_lambda`; `MiuraHamiltonian` and `MiuraSystem` take `λ`.

At `λ = 0` the image of the map is a half-space and, sharply, the fields whose Hill operator
is positive definite — which excludes every standard example. **The sign convention was never
what caused that**: flipping it moves the half-space and the solitons together. The
obstruction is spectral, a soliton being a reflectionless potential *with* a bound state.

The parameter is not a shift bolted on but the one the Riccati substitution always had:
`u = -(v² + v_x) - λ` turns into `(-∂ₓ² - u)ψ = λψ`, so `λ` is an eigenvalue parameter of the
very Hill operator that decides invertibility, and a preimage exists exactly when `λ < λ₀`.
That can always be arranged. Measured, with `miura_lambda` a tenth below each threshold, `cos
x`, `sin x + 0.4cos2x`, the single soliton and the two-soliton all enter the chart with a
round-trip error of `10⁻¹⁵`, and the pushforward stays exactly antisymmetric with its
Jacobiator at `10⁻¹⁶`.

`L` does not see a constant, so `MiuraBracket` is unchanged. What moves is *which* bracket the
pushforward is: `L P¹ Lᵀ = P² - 4λP¹`, a member of the bi-Hamiltonian pencil. On the `u` side
the extra term is a constant-speed translation, so the trajectory is the KdV solution in a
uniformly moving frame and all three invariants survive. `scripts/kdv.jl` therefore carries
Miura runs for every case now, not only the `miura` one.

Note that `L P¹ Lᵀ` is **not** the Galerkin `P² - 4λP¹`; measured, they differ by 30-57 %.
That difference is the Jacobi defect the Miura construction exists to remove.

### Changed — `LagrangeSpace` stores its tabulation sparsely

Following `SimpleSplines.SplineQuadrature`: `Φ` is a `SparseMatrixCSC` with `(p+1)` structural
entries per quadrature column (0.8 % of a dense table at `p = 2`, `ne = 192`), the mass matrix
is sparse and wrapped in a `FactorizedMass`, `mixed_matrix` is memoised, and the basis
integrals and a scratch buffer are assembled once. `inverse_mass_matrix` stays dense.

`AffineBracket` no longer densifies the tabulation it is handed — its `Ψ` field was declared
`Matrix{T}` and its constructor called `Matrix`, which threw away the *spline* space's
already-sparse table for the second KdV bracket and the first Camassa-Holm one, both of which
contract it on every Newton iteration. One `jacobian` of the second flow at `N = 384`: **7.24
ms before, 1.26 ms after**, identical results; `bracket_directional` alone 25× faster, the
`LagrangeSpace` assemblies about 30×. The two `O(N³)` tensor routines densify locally and
deliberately — they random-access the table rather than contracting it, and are off the time
loop by construction.

### Removed — the local `LapackLU`

SimpleSolvers 0.12.2 ships it, with a superset of the behaviour: a `BlasFloat` guard,
feature-detected allocation-free `getrf!`, `singular_index`, `factorization` and the full
`solve!` surface. `src/linearsolver.jl` is gone and the name is re-exported from SimpleSolvers,
so `linear_solver_method = LapackLU()` keeps working for anyone who has only done
`using PoissonBrackets`. `[compat]` requires `SimpleSolvers = "0.12.2"`.

The local `ldiv!` threw a message-carrying exception; upstream throws
`LinearAlgebra.SingularException(info)`. The explanatory text now lives in `_step!`, which
already owns the "step size is probably too large" warning.

### Fixed

- **`project!` never worked on a `LagrangeSpace`.** `mass_factorization` returned a bare dense
  `Cholesky`, and the generic `project!` goes through `mass_solve!`, which has methods only for
  a `MassOperator`. Every in-place projection onto a Lagrange space was a `MethodError`, and
  nothing exercised it. Now tested.
- **`integrate` dropped the tail of every run.** `nout = nsteps ÷ stride` ran `nout*stride`
  steps while the tables reported `nsteps` — 766400 of a claimed 766781 for the explicit
  cosine run, ending at `T = 99.95`. Now a ceiling with a final short window, as
  `plot_kdv_energies.py` does.
- **The Burgers `√u` convention was inconsistent with its own docstring.**
  `to_sqrt_variables` used `ū = 2√u` while the displayed bracket was `K/4`, which belongs to
  `ū = √u`. Both are now the manuscript's `ū = √u`. The transformation and the factor move
  together, and pairing them wrongly makes the flow four times too slow — which no
  conservation test can see, a constant rescaling of a Poisson field preserving every Casimir,
  every energy bound and every convergence *rate*. There is now a test comparing against the
  pushforward of the `u`-space field, which is the only thing that catches it. The same error
  was found and fixed in `verify_burgers_discretisation.py`.
- `scripts/kdv.jl` used `min(dt_rk4, case.dt)` for the explicit runs where
  `plot_kdv_energies.py` uses `dt_rk4` outright. Invisible on four of six cases; on `soliton5`
  it forced 0.0625 instead of 0.065 and hid the point the benchmark makes. The output stride
  now rounds as the Python does.
- Claims inherited from the manuscript's own list of known defects: the Miura image bound is
  `C₀ ≥ 0` and not strict; the critical set of `Dℳ_h` is quoted as evidence rather than proof,
  with the exact statement that *does* hold added — `∫v_h` is a Casimir of `P¹` in that chart,
  so the flow cannot reach the degeneracy from data off it, now a test; and Ge-Marsden is
  given as the reason to expect the energy/structure trade-off rather than a proof of it for
  these brackets, whose non-degeneracy hypotheses are unchecked.

### Added — documentation

- **Seven new pages, for the material the port added and never explained.** Commit `f897fb0`
  added four source modules, twenty-four scripts and three plotting routines, and touched
  `docs/` only with bare `@autodocs` blocks — everything was reachable and almost nothing was
  explained.

  | page | covers |
  |:--|:--|
  | `diagnostics.md` | the theory of every measurement the package makes, organised by what is being measured: whether the *bracket* is Poisson, whether the *flow* conserves what it should, whether the *method* preserves the structure. Why `jacobi_residual` is normalised the way it is, why a flat residual is the answer and not a failure to converge, why `integrate` reports windowed maxima, why `drift` throws on a vanishing reference, and how to tell a drift from accumulated round-off |
  | `liepoisson.md` | the second manuscript: the Dubrovin-Novikov factorisation, prescribing the Casimir, the structure-constant condition and the closure/consistency split, the algebra zoo, and the two four-bracket ansätze |
  | `dirac.md` | the Dirac companion: the Dirac bracket, the tensorial transport of the Jacobiator that makes reduction faithful but never curative, the four realisations, the broken hierarchical space, and the repair search |
  | `bea.md` | §6 of the KdV notes: the modified Hamiltonian, `H₂` as the momentum of the projected translation, the closed-form defect, the exact increment lemma, and how to read `sweepplot` |
  | `aliasing.md` | the zero-mode theorem and the cocycle obstruction — why Zeitlin-style aliasing cannot repair the second KdV bracket |
  | `nambu.md` | the Vandermonde lemma, the wedge tensor that repackages nothing, the one genuinely new three-bracket, and the decomposability theorem |
  | `scripts.md` | an index of all twenty-two verification scripts against the claims they check and the pages that carry the theory |
  | `fourbrackets.md` | the third manuscript: Kulkarni-Nomizu positivity, the two antisymmetry conditions and why there are two families, the Gardner operators, both reductions to `∫u[A_u,B_u]` and the factor of two between them, the Plücker relation that collapses the Jacobi obstruction, the log-entropy weight's singularity, and Section 6's construction that generates nothing |

- **Three module-overview docstrings were being dropped from the build without a word.** The
  headers of `src/algebras.jl`, `src/dirac.jl` and `src/hierarchical.jl` were `@doc raw"…"`
  blocks attached to no binding — the `C[m,i,j]` index convention, the definition of the Dirac
  bracket and the second-class condition, and the whole rationale for the broken space reached
  no reader of the manual. The prose now lives on the pages above and the file headers are plain
  comments.

- `DiscreteGradient` was exported with no docstring, so it rendered nowhere and
  `[`DiscreteGradient`](@ref)` could not resolve. It has one.

- **`docs/src/discretisation.md` contradicted `docs/src/verification.md`.** Testing trap 1 still
  asserted that every three-dimensional antisymmetric bracket satisfies the Jacobi identity — the
  claim recorded three sections away as an erratum. Corrected to Bianchi class A.

- The `exact.jl`, `algebras.jl`, `dirac.jl`, `hierarchical.jl` and `diagnostics.jl` `@autodocs`
  blocks move from `library.md` to their topic pages, as the equation modules already had.
  `@index` is unchanged and complete.

- `docs/src/verification.md`, recording what the package has been checked against claim by
  claim: the RK4 comparison, the `ρ` cross-check against three independent sources, the
  divergences found and fixed, the errata found in the Python, every item of the manuscript's
  error list dispositioned as inherited-and-fixed, paper-only or still open, and the fact that
  **Camassa-Holm has no reference at all** — neither manuscript covers it and neither
  `Scripts/` directory mentions it, so it is a prototype in the strict sense and must not be
  read as verified.
- A test asserting `stability_limit` against the `ρ = 115, 274, 907, 2170` of the notes, so a
  change in the assembly fails a test rather than drifting away from the manuscript.

### The RK4 question

The RK4 results for the `cos` and single-soliton runs were suspected of disagreeing with the
Python. **They do not.** `Δt` is identical to full precision (766781 and 23833 steps on both
sides) and every quantity above round-off agrees to every digit the tables print. What differs
is confined to round-off-dominated columns, and `|ΔC₀|` — conserved exactly in exact
arithmetic — is the internal control that shows it: 7.7e-12 in Python against 4.7e-12 in Julia
over 766781 steps. Any `|ΔH|` reported at that order is arithmetic, not dynamics.

## [0.1.0]

Initial release.

### Added

- Discrete function spaces: `SplineSpace`, wrapping the periodic B-spline basis of
  SimpleSplines, and `LagrangeSpace`, a periodic nodal Lagrange finite element space built
  on the reference-element bases of CompactBasisFunctions.
- Discrete brackets: `ConstantBracket`, `AffineBracket`, `GaugedBracket` and `MiuraBracket`,
  with a common interface of `poisson_matrix`, `poisson_apply`, `isantisymmetric` and
  `jacobi_residual`.
- Discrete Hamiltonians and Casimirs for the Burgers, KdV and Camassa-Holm equations, each
  with an analytic gradient and Hessian.
- Geometric integrators: implicit midpoint, the average vector field method, the Gonzalez
  discrete gradient in the Euclidean and mass metrics, explicit Euler and Runge-Kutta 4,
  and a projection method that imposes several invariants at once.
- Diagnostics: fused evaluation of the invariants, integration with windowed error
  envelopes, and the Poisson-map defect.
- The Miura chart: `MiuraBracket`, the pushforward of the first bracket along the Miura map
  (`u = v² + v_x` at the time; see Unreleased for the sign convention and for `λ`), antisymmetric *and* Poisson to round-off at every resolution where the
  Galerkin second bracket fails the Jacobi identity at order one; `MiuraSystem`, carrying
  the mKdV degrees of freedom in which the bracket is constant; `miura_map`, `miura_invert`
  and `hill_lambda0` for moving between the two charts and for deciding whether a given
  field has a preimage at all.
- Initial conditions for every example of the manuscript: `cosine`, `soliton`,
  `two_soliton`, `three_solitons`, `five_solitons` and `miura_initial_v`, with
  `scripts/kdv.jl` running all six cases and writing a markdown summary per case.
- `LapackLU`, a LAPACK-backed factorisation plugged into SimpleSolvers' linear-solver
  extension points; submitted upstream as
  [JuliaGNI/SimpleSolvers.jl#183](https://github.com/JuliaGNI/SimpleSolvers.jl/pull/183),
  merged and released in SimpleSolvers 0.12.2, and consequently removed again here — see
  Unreleased.

### The Newton residual tolerance

`f_abstol` is an **absolute** tolerance, and it has to sit above the round-off floor of the
residual — which scales with the number of degrees of freedom *and* with the amplitude of the
field. Both ways of getting it wrong were measured, on implicit midpoint on the second KdV
flow, and neither is acceptable:

| `f_abstol` | iterations, ‖u‖≈4 | \|ΔH₂\|, ‖u‖≈4 | iterations, ‖u‖≈1 | \|ΔH₂\|, ‖u‖≈1 |
|:--|--:|--:|--:|--:|
| `max(8,N)·eps` | 40.0 | 5.8e-15 | 40.0 | 5.7e-16 |
| **`4·max(8,N)·eps·‖û₀‖`** | **3.0** | **6.4e-15** | **2.5** | **7.1e-16** |
| `16·…` | 2.9 | 6.2e-14 | 2.0 | 5.7e-16 |
| `1e-11` | 2.6 | 9.6e-13 | 2.0 | 5.7e-16 |

Below the floor Newton reaches its limit in three iterations and then spins to the iteration
cap making no progress, reporting non-convergence on a step it has in fact solved. Above it
Newton stops early and the conservation laws pay — an order of magnitude of `H₂` per factor
of sixteen. The window is not wide, which is why this is a formula and not a constant, and
why `Integrator` takes `û₀`: an absolute tolerance cannot serve a field of size one and a
field of size four at once.

`min_iterations = 1`, and a step whose first attempt does not converge is redone with a line
search rather than accepted. That fallback is a safeguard against an over-large step size,
not a remedy for a mis-set tolerance, and it does not fire at these settings.

### Performance

An implicit step at `N = 384` went from 45 ms to 4.2 ms. In order of what each change bought:

- SimpleSolvers' hand-written scalar LU accounted for 74 % of a step; `LapackLU` replaces
  the factorisation while keeping the whole nonlinear driver. (It now lives upstream.)
- The basis tabulation is now sparse in SimpleSplines, so the assembly of one Hessian at
  `N = 384` fell from 6.0 ms to 0.4 ms, and the constant assemblies are memoised rather
  than rebuilt on every Newton iteration.
- The default line search is `Static`: at these step sizes Newton converges in two or three
  iterations and backtracking only adds residual evaluations.
- Mass solves go through SimpleSplines' `MassOperator`, which uses planned real FFTs on a
  uniform mesh, where the matrix is circulant. This is the right representation of the
  operator but not where the time went: 0.064 ms to 0.001 ms, against a 4.2 ms step.
