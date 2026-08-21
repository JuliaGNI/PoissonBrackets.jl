# PoissonBrackets

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNI.github.io/PoissonBrackets.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaGNI.github.io/PoissonBrackets.jl/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build Status](https://github.com/JuliaGNI/PoissonBrackets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaGNI/PoissonBrackets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaGNI/PoissonBrackets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/PoissonBrackets.jl)

Structure-preserving discretisations of the Poisson brackets of one-dimensional Hamiltonian
partial differential equations, with the discrete brackets, Hamiltonians, Casimirs and
geometric integrators kept as separate, composable objects.

The package collects the prototype discretisations and numerical experiments behind two
manuscripts on discrete Poisson brackets, and turns the structural claims they make — exact
antisymmetry, the Jacobi identity or its failure, conservation of Hamiltonians and Casimirs —
into a test suite.

Three equations are covered:

| equation | bracket | basis |
|:--|:--|:--|
| Burgers | `J(u) = g(u) K g(u)` with `g = √u`, made constant by `ū = √u` | periodic nodal Lagrange |
| KdV | the bi-Hamiltonian pair `P¹` (constant, Poisson) and `P²(û)` (antisymmetric, not Jacobi) | periodic B-spline |
| Camassa-Holm | the standard bi-Hamiltonian pair — *prototype* | periodic B-spline |

## Development

To run the test suite before every push, enable the repository's git hooks:

```sh
git config core.hooksPath .githooks
```
