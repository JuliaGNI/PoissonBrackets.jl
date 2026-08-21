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
