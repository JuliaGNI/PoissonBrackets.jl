# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- The Miura chart: `MiuraSystem` carries the mKdV degrees of freedom, where the bracket is
  the constant `P¹`, with `miura_map`, `miura_invert` and `hill_lambda0` for moving between
  the charts and for telling in advance whether a given field has a preimage at all.
- `LapackLU`, a LAPACK-backed linear solver for the Newton iteration. SimpleSolvers' own
  scalar `LU` accounted for 74 % of an implicit step at `N = 384`; submitted upstream as
  [JuliaGNI/SimpleSolvers.jl#183](https://github.com/JuliaGNI/SimpleSolvers.jl/pull/183).

### Performance

An implicit step at `N = 384` went from 45 ms to 4.2 ms, from three changes, in order of
what they were worth:

1. the basis tabulation is stored sparsely, so an assembly costs `O(N p² nq)` rather than
   `O(N² n nq)` — one Hessian went from 6.0 ms to 0.4 ms;
2. the Newton factorisation goes through LAPACK;
3. the constant assemblies are memoised instead of being rebuilt on every Newton iteration.

A fourth, the circulant FFT mass solve on uniform meshes, is correct and cheap but was never
where the time went: at `N = 384` a mass solve is 0.001 ms.
- The Miura chart: `MiuraBracket`, the pushforward of the first bracket along
  `u = v² + v_x`, which is antisymmetric *and* Poisson to round-off at every resolution
  where the Galerkin second bracket fails the Jacobi identity at order one; `MiuraSystem`,
  carrying the mKdV degrees of freedom in which the bracket is constant; `miura_map`,
  `miura_invert` and `hill_lambda0` for moving between the two charts and for deciding
  whether a given field has a preimage at all.
- Initial conditions for every example of the manuscript: `cosine`, `soliton`,
  `two_soliton`, `three_solitons`, `five_solitons` and `miura_initial_v`, with
  `scripts/kdv.jl` running all six cases.
- `LapackLU`, a LAPACK-backed factorisation plugged into SimpleSolvers'
  linear-solver extension points.

### Performance

An implicit step at `N = 384` went from 45 ms to 4.2 ms. In order of what it bought:

- SimpleSolvers' hand-written scalar LU accounted for 74 % of a step; `LapackLU` replaces
  the factorisation while keeping the whole nonlinear driver.
- The basis tabulation is now sparse in SimpleSplines, so the assembly of one Hessian at
  `N = 384` fell from 6.0 ms to 0.4 ms, and the constant assemblies are memoised rather
  than rebuilt on every Newton iteration.
- The default line search is `Static`: at these step sizes Newton converges in two or
  three iterations and backtracking only adds residual evaluations.
- Mass solves go through SimpleSplines' `MassOperator`, which uses planned real FFTs on a
  uniform mesh where the matrix is circulant. This is the right representation but not
  where the time went: 0.064 ms to 0.001 ms, against a 4.2 ms step.
