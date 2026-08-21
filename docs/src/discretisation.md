```@meta
CurrentModule = PoissonBrackets
```

# Discretisation

## Spaces

Two discrete function spaces are provided, and the choice between them is not free — it is
the discretisation each problem actually requires.

[`SplineSpace`](@ref) wraps the periodic B-spline basis of SimpleSplines. Being
``\mathcal{C}^{p-1}``, it can carry the third-derivative term of the KdV brackets: after one
integration by parts, ``-\int \phi_k \phi_l''' = \int \phi_k' \phi_l''``, degree ``p \ge 2``
suffices.

[`LagrangeSpace`](@ref) is a periodic nodal Lagrange finite element space. It is only
``\mathcal{C}^0``, but it is *nodal*, and that is what the Burgers construction needs: the
transformation ``\bar{u} = 2\sqrt{u}`` is applied coefficient by coefficient, which is only
justified when the coefficients are values at points.

Both answer the same small interface — [`basis_values`](@ref), [`quadrature_weights`](@ref),
[`mass_matrix`](@ref) — so every assembly downstream is written once.

## Assembly

Every matrix here is a weighted contraction of tabulated basis derivatives:

```math
\int_\Omega f(x) \, D^a \phi_k \, D^b \phi_l \, dx
  = \Phi_a \, \mathrm{diag}(f \odot w) \, \Phi_b^T .
```

That is [`weighted_matrix`](@ref); [`mass_matrix`](@ref), [`stiffness_matrix`](@ref) and
[`derivative_matrix`](@ref) are the constant-coefficient cases.

## Brackets

```math
\{A, B\}_d = \sum_{i,j} \frac{\partial A}{\partial \hat{u}_i} \,
             \mathbb{P}_{ij} (\hat{u}) \, \frac{\partial B}{\partial \hat{u}_j}
```

| type | antisymmetric | Jacobi |
|:--|:--|:--|
| [`ConstantBracket`](@ref) | yes | yes, automatically |
| [`GaugedBracket`](@ref) | yes | yes, by the theorem it is built on |
| [`MiuraBracket`](@ref) | yes | yes, being a pushforward |
| [`AffineBracket`](@ref) | yes | **no** |

## Why the skew-symmetrised form matters

Every bracket here is assembled in explicitly skew-symmetrised form. This is not cosmetic.

The unsymmetrised assembly of the second KdV bracket is the *same matrix* in exact
arithmetic — the term that differs is symmetric in ``k \leftrightarrow l`` and drops out of
the skew part. But it is antisymmetric only if the quadrature integrates a total derivative
of degree ``3p-1`` exactly. Written skew-symmetrised, antisymmetry holds **identically**, at
any number of quadrature points and on any mesh — and with it the exact conservation of the
generating Hamiltonian, which follows from antisymmetry alone:

```math
\dot{H} = \left( \frac{\partial H}{\partial \hat{u}} \right)^T \mathbb{P}
          \left( \frac{\partial H}{\partial \hat{u}} \right) = 0 .
```

Antisymmetry is thus the one property that has been made free. Accuracy has not: consistency
still needs degree ``3p-1``, and the mass still needs ``2p-1``.

## The Jacobi identity

[`AffineBracket`](@ref) does not satisfy it, and this is a result rather than a defect. The
operator ``4u\partial_x + 2u_x - \partial_x^3 = 2(u\partial_x + \partial_x u) -
\partial_x^3`` is the Virasoro Lie-Poisson structure, so the Jacobi identity would require
the discrete coefficients to close into a Lie algebra — a quadratic condition on the basis
that no choice of quadrature and no regrouping of the integrand can deliver. The normalised
[`jacobi_residual`](@ref) is of order one and *flat* under refinement:

```@example jac
using PoissonBrackets, Random
Random.seed!(1)
[jacobi_residual(kdv_bracket_2(SplineSpace(N, 3)), randn(N)) for N in (12, 16, 20, 24)]
```

[`MiuraBracket`](@ref) is the construction that does give an exactly Poisson second bracket,
by pushing the first forward along ``u = v^2 + v_x``. The price is locality: it is dense
where the Galerkin assembly is banded, and those two interior inverse mass matrices are
precisely the nonlocality that carries the Jacobi identity.

## Testing traps

Four of these are worth knowing before writing a test against any of this:

1. **Never use ``\mathfrak{so}(3)`` as a positive control** for a Jacobi or
   structure-constant routine. Every three-dimensional antisymmetric bracket satisfies the
   Jacobi identity identically, so it passes even after a generator has been rescaled off
   the algebra. Use ``\mathfrak{se}(3)``, with a random antisymmetric `c` in dimension five
   as the negative control.
2. **An antisymmetric matrix has even rank.** For an even number of degrees of freedom the
   corank is two, and the extra sawtooth kernel vector is a *spurious* Casimir. Use an odd
   number.
3. **A uniform mesh confirms some identities for the wrong reason** — its assemblies are
   circulant. Use `RandomMesh` for anything claimed to hold on any mesh.
4. **Testing on a single smooth mode hides the cubic obstruction.** Excite the full spectrum
   with random degrees of freedom.
