# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

`0.1.0` has not shipped, so all of this may be folded into it; it is kept separate
because the KdV sign convention below changes what every number in the package means.

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
  | `scripts.md` | an index of all eighteen verification scripts against the claims they check and the pages that carry the theory |

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
