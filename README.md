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
into a test suite and eighteen verification scripts.

Three equations are covered:

| equation | bracket | basis |
|:--|:--|:--|
| Burgers | `J(u) = g(u) K g(u)` with `g = √u`, made constant by `ū = √u` | periodic nodal Lagrange |
| KdV | the bi-Hamiltonian pair `P¹` (constant, Poisson) and `P²(û)` (antisymmetric, not Jacobi) | periodic B-spline |
| Camassa-Holm | the standard bi-Hamiltonian pair — *prototype* | periodic B-spline |

and, for the second manuscript and its Dirac companion, the finite-dimensional machinery those
brackets are special cases of:

| topic | what is here |
|:--|:--|
| Lie-Poisson brackets | structure constants `c[m,i,j]` and the closure condition; `se(3)`, `so(N)`, Zeitlin's sine algebra, the truncated Witt, polynomial and torus algebras |
| Dirac reduction | Schur complement, Dirac projector and tensor, and the tensorial transport of the Jacobiator that makes reduction faithful but never curative |
| the broken hierarchical space | the DG realisation in which `[V₁, V₁] ⊆ V₁ ⊕ V₂` holds exactly, assembled over `Rational{BigInt}` with no CAS |
| exact linear algebra | `rref`, `exact_rank` and `kernel` over ℚ, where the SVD-based `rank` and `nullspace` do not run |

The [documentation](https://JuliaGNI.github.io/PoissonBrackets.jl/latest/) carries the theory:
what each diagnostic measures, the backward error analysis of the implicit midpoint rule, and
the no-go results for aliasing and for antisymmetric three-brackets.

## Development

To run the test suite before every push, enable the repository's git hooks:

```sh
git config core.hooksPath .githooks
```
