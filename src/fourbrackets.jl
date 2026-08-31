#
# The brackets of `poisson-brackets-from-four-brackets.tex`, on the periodic two-torus.
#
# Everything here is a functional of sampled fields on a `TorusGrid`, not a matrix on a
# `DiscreteSpace`: this manuscript works in the continuum and asks which four-brackets reduce
# to a Lie-Poisson bracket, whereas `brackets.jl` and `algebras.jl` ask what survives a
# discretisation. The two lines of argument meet only in the answer, so they share no code.
#
# The narrative — which reduction holds for which weight, and why the log-entropy weight is
# singular — is on the `Poisson brackets from four-brackets` page of the documentation.
#

@doc raw"""
    gardner_x(g, f, k)
    gardner_y(g, f, k)

The Gardner operators of equation (4.3),

```math
G_x(f, k) = \tfrac{1}{2} \left( f k_x - k f_x \right) ,
\qquad
G_y(f, k) = \tfrac{1}{2} \left( f k_y - k f_y \right) ,
```

the antisymmetric bilinear combination of a field and one derivative of another.

Antisymmetric in `(f, k)` pointwise. Two rewritings are used constantly and are what the
auxiliary identities of Section 4.2 record: ``G_x(f,k) = f k_x - \tfrac{1}{2}\partial_x(fk)``
exhibits it as a flux plus a total derivative, so that ``\int F\, G_x`` integrates by parts
cleanly; and ``G_x(f,k) = -\tfrac{1}{2} k^2 \partial_x(f/k)`` is the quotient form of
Lemma 4.1, which is what turns the Gardner two-bracket into a canonical bracket of ratios.

Contrast [`symmetric_x`](@ref), the symmetric combination, which is a total derivative and
therefore generates nothing.
"""
gardner_x(g::TorusGrid, f, k) = (f .* ∂x(g, k) .- k .* ∂x(g, f)) ./ 2

@doc (@doc gardner_x)
gardner_y(g::TorusGrid, f, k) = (f .* ∂y(g, k) .- k .* ∂y(g, f)) ./ 2

@doc raw"""
    symmetric_x(g, f, k)
    symmetric_y(g, f, k)

The symmetric counterparts of the Gardner operators, ``M(f,k) = f k_x + k f_x`` and
``N(f,k) = f k_y + k f_y``, used to build the two-bracket of Section 6.

# Why the Section 6 bracket vanishes

These are total derivatives: ``M(f,k) = \partial_x(fk)`` and ``N(f,k) = \partial_y(fk)``.
The two-bracket assembled from them,
``\int M(A_u,S_u) N(B_u,S_u) - M(B_u,S_u) N(A_u,S_u)``, is therefore
``\int [A_u S_u, B_u S_u]`` — a canonical bracket of two functions, whose integral over a
periodic domain vanishes identically. Section 6's negative result is that construction
producing nothing at all, for any ``S``.

They are written out as ``f k_x + k f_x`` rather than as ``\partial_x(fk)`` deliberately.
Collapsing them to a single derivative would make the vanishing a consequence of how the
code is written rather than a fact about the bracket, and the check in
`scripts/verify_fourbracket_identities.jl` would then be testing nothing.
"""
symmetric_x(g::TorusGrid, f, k) = f .* ∂x(g, k) .+ k .* ∂x(g, f)

@doc (@doc symmetric_x)
symmetric_y(g::TorusGrid, f, k) = f .* ∂y(g, k) .+ k .* ∂y(g, f)

@doc raw"""
    gardner_4bracket_density(g, a, b, c, d)
    gardner_4bracket(g, a, b, c, d)

The Gardner-like four-bracket

```math
(A, B; C, D) = \int_\Omega \left[ G_x(A_u, B_u) G_y(C_u, D_u)
                                - G_x(C_u, B_u) G_y(A_u, D_u) \right] ,
```

as its pointwise integrand and as the integral.

The `_density` form exists because some claims are about the integrand rather than the
integral — whether it stays bounded under refinement, for instance, which is how
`scripts/verify_fourbracket_log_entropy.jl` separates a genuine pole from a merely large
value. See [`weighted_4bracket`](@ref) for the weighted version.

!!! note "Only condition (i)"
    This bracket is antisymmetric under exchanging the first and third slots but **not**
    under exchanging the pairs `(1,2)` with `(3,4)`. The symmetric four-bracket
    [`symmetric_4bracket`](@ref) satisfies exactly the other one. That is Lemma 3.1, and it
    is why the manuscript carries two families rather than one — measure it with
    [`antisymmetry_residuals`](@ref).
"""
gardner_4bracket_density(g::TorusGrid, a, b, c, d) = gardner_x(g, a, b) .*
                                                     gardner_y(g, c, d) .-
                                                     gardner_x(g, c, b) .*
                                                     gardner_y(g, a, d)

@doc (@doc gardner_4bracket_density)
gardner_4bracket(g::TorusGrid, a, b, c, d) = integrate(g, gardner_4bracket_density(
    g, a, b, c, d))

@doc raw"""
    symmetric_4bracket_density(g, a, b, c, d)
    symmetric_4bracket(g, a, b, c, d)

The four-bracket of Section 5 built from the symmetric pairing,

```math
(A, B; C, D) = \int_\Omega \left[ A_u B_u [C_u, D_u] - C_u D_u [A_u, B_u] \right] ,
```

as its pointwise integrand and as the integral.

!!! note "Only condition (ii)"
    The mirror image of [`gardner_4bracket`](@ref): antisymmetric under exchanging the pairs
    `(1,2)` with `(3,4)`, and not under exchanging the first and third slots. Lemma 3.1 again.
"""
symmetric_4bracket_density(g::TorusGrid, a, b, c, d) = a .* b .*
                                                       canonical_bracket(g, c, d) .-
                                                       c .* d .* canonical_bracket(g, a, b)

@doc (@doc symmetric_4bracket_density)
symmetric_4bracket(g::TorusGrid, a, b, c, d) = integrate(g, symmetric_4bracket_density(
    g, a, b, c, d))

@doc raw"""
    gardner_2bracket_density(g, a, b, s)
    gardner_2bracket(g, a, b, s)
    symmetric_2bracket_density(g, a, b, s)
    symmetric_2bracket(g, a, b, s)

The two-brackets induced by fixing the entropy: `s` is the sampled ``S_u``, and each is its
four-bracket with the entropy in both the second and the fourth slot,

```math
\{A, B\}_S = (A, S; B, S) .
```

That substitution is the whole mechanism of Sections 4 and 5, so these are defined as
exactly that call rather than written out again.

# The reductions, and the factor of two

Both collapse to a canonical bracket weighted by the structure function, and the two
normalisations differ:

```math
\{A,B\}_S^{\mathrm{Gardner}} = \tfrac{1}{2} \int S_u^2 [A_u, B_u] ,
\qquad
\{A,B\}_S^{\mathrm{sym}} = \int S_u^2 [A_u, B_u] .
```

So the Lie-Poisson bracket ``\int u [A_u, B_u]`` is reached from the Gardner family at
``S_u = \sqrt{2u}``, that is ``S = \tfrac{2}{3}\sqrt{2}\int u^{3/2}``, and from the symmetric
family already at ``S_u = \sqrt{u}``. Both are exact, not asymptotic — see
[`lie_poisson_2bracket`](@ref).

Equivalently, in quotient form, the Gardner two-bracket is
``\tfrac{1}{4}\int S_u^4 [A_u/S_u, B_u/S_u]``, which is Proposition 4.3. That form involves
the non-band-limited quotient ``A_u/S_u`` and is therefore certified by refinement rather
than at roundoff.
"""
gardner_2bracket_density(g::TorusGrid, a, b, s) = gardner_4bracket_density(g, a, s, b, s)

@doc (@doc gardner_2bracket_density)
gardner_2bracket(g::TorusGrid, a, b, s) = gardner_4bracket(g, a, s, b, s)

@doc (@doc gardner_2bracket_density)
symmetric_2bracket_density(g::TorusGrid, a, b, s) = symmetric_4bracket_density(
    g, a, s, b, s)

@doc (@doc gardner_2bracket_density)
symmetric_2bracket(g::TorusGrid, a, b, s) = symmetric_4bracket(g, a, s, b, s)

@doc raw"""
    weighted_2bracket_density(g, a, b, s, w)
    weighted_2bracket(g, a, b, s, w)
    weighted_4bracket_density(g, a, b, c, d, w)
    weighted_4bracket(g, a, b, c, d, w)

The weighted brackets of Section 5.1: the symmetric bracket of the same arity with the
sampled weight ``\omega(u)`` inside the integral.

```math
\{A, B\}_S^\omega = \int_\Omega \omega(u)
        \left[ A_u S_u [B_u, S_u] - B_u S_u [A_u, S_u] \right] .
```

The weight is a plain factor on the integrand, which is why these are defined as `w .*` the
corresponding [`symmetric_2bracket_density`](@ref) rather than written out again.

# What the weight buys, and what it costs

The weighted two-bracket reduces to ``2 \int \Phi(u) [A_u, B_u]`` with
``\Phi' = \omega s' s''``, so demanding the Lie-Poisson bracket ``\int u [A_u,B_u]`` fixes

```math
\omega = \frac{1}{2 s' s''} .
```

For ``s = u \log u`` that is ``\omega = u / (2(1 + \log u))``, with a simple pole at
``u = e^{-1}`` where ``s'`` vanishes. The pole is not removable and not an artefact: it
cancels against the zero of ``S_u`` in the two-bracket, so ``\{A,B\}_S`` is regular across
``u = e^{-1}``, but the four-bracket with four independent arguments has no such zero to
cancel it and genuinely diverges. And no weight regular at the critical level can be
substituted, because ``c'(u_*) = 2\omega(u_*) s'(u_*) s''(u_*) = 0`` there whereas the
Lie-Poisson bracket needs ``c' \equiv 1``.

What *is* available is a shift: ``\int u`` is a Casimir of every bracket in the family, so
``s_\alpha = u \log u + \alpha u`` generates the same bracket while moving the critical level
to ``u_* = e^{-(1+\alpha)}``. On a state space whose fields take values in a compact
subinterval of ``(0,\infty)``, some ``\alpha`` puts ``u_*`` outside it and the weight is
smooth on the whole range.

This is the same obstruction that the discrete Lie-Poisson brackets meet when the Casimir is
prescribed — see the discussion of the intrinsic singularity on the
[Discrete Lie-Poisson brackets](@ref) page — seen here from the four-bracket side.
"""
weighted_2bracket_density(g::TorusGrid, a, b, s, w) = w .*
                                                      symmetric_2bracket_density(g, a, b, s)

@doc (@doc weighted_2bracket_density)
weighted_2bracket(g::TorusGrid, a, b, s, w) = integrate(g, weighted_2bracket_density(
    g, a, b, s, w))

@doc (@doc weighted_2bracket_density)
weighted_4bracket_density(g::TorusGrid, a, b, c, d, w) = w .* symmetric_4bracket_density(
    g, a, b, c, d)

@doc (@doc weighted_2bracket_density)
weighted_4bracket(g::TorusGrid, a, b, c, d, w) = integrate(g, weighted_4bracket_density(
    g, a, b, c, d, w))

@doc raw"""
    lie_poisson_2bracket(g, a, b, u)

The two-dimensional Lie-Poisson bracket ``\int_\Omega u\, [A_u, B_u]``, the target every
reduction in the manuscript is aimed at.

It is the bracket of the ideal-fluid vorticity equation, and the point of the four-bracket
constructions is that it is *reached exactly* — at ``S_u = \sqrt{2u}`` for the Gardner family
and ``S_u = \sqrt{u}`` for the symmetric one — rather than approximated. See
[`gardner_2bracket`](@ref).
"""
lie_poisson_2bracket(g::TorusGrid, a, b, u) = integrate(g, u .* canonical_bracket(g, a, b))

@doc raw"""
    antisymmetry_residuals(fourbracket, g, a, b, c, d) -> (; slot_exchange, pair_exchange)

The two antisymmetry conditions of Lemma 3.1, measured for the four-bracket `fourbracket`
and normalised by the size of the bracket itself:

```math
\text{(i)}\quad K(A,B;C,D) + K(C,B;A,D) = 0 ,
\qquad
\text{(ii)}\quad K(A,B;C,D) + K(C,D;A,B) = 0 .
```

Condition (i) is antisymmetry under exchanging the first and third slots; condition (ii)
under exchanging the pair `(1,2)` with the pair `(3,4)`.

# Both are sharp, and neither bracket has both

Pass [`gardner_4bracket`](@ref) and `slot_exchange` vanishes while `pair_exchange` is of
order one; pass [`symmetric_4bracket`](@ref) and it is the other way round. A test that only
checked the vanishing residual would be satisfied by a bracket that was accidentally
antisymmetric in every slot, which is why the non-vanishing one is asserted too — and why
the test fields must be generic rather than symmetric under any exchange of their own.
"""
function antisymmetry_residuals(fourbracket, g::TorusGrid, a, b, c, d)
    base = abs(fourbracket(g, a, b, c, d))
    return (
        slot_exchange = abs(fourbracket(g, a, b, c, d) + fourbracket(g, c, b, a, d)) / base,
        pair_exchange = abs(fourbracket(g, a, b, c, d) + fourbracket(g, c, d, a, b)) / base)
end

@doc raw"""
    plucker_residual(g, a, b, c, u)

The pointwise Plücker relation of Lemma 4.4,

```math
[A_u, B_u][C_u, u] + [B_u, C_u][A_u, u] + [C_u, A_u][B_u, u] = 0 ,
```

returned as the array whose largest entry is the violation.

It holds identically, for any four fields on a two-dimensional domain, and it is the reason
the Jacobi obstruction of [`jacobi_residual`](@ref) collapses from twenty-four terms to the
single condition on the structure function. Two dimensions is essential: the relation is the
vanishing of a ``3 \times 3`` determinant of gradients that live in a two-dimensional space,
so it has no analogue in higher dimension.
"""
plucker_residual(g::TorusGrid, a, b, c, u) = canonical_bracket(g, a, b) .*
                                             canonical_bracket(g, c, u) .+
                                             canonical_bracket(g, b, c) .*
                                             canonical_bracket(g, a, u) .+
                                             canonical_bracket(g, c, a) .*
                                             canonical_bracket(g, b, u)

@doc raw"""
    jacobi_residual(g, a, b, c, cu, dcu; normalised = true)

The Jacobi obstruction of Theorem 4.5 for the bracket
``\{A,B\}_c = \int c(u) [A_u, B_u]``, evaluated on the three functional derivatives
`a`, `b`, `c` with the structure function sampled as `cu` ``= c(u)`` and
`dcu` ``= c'(u)`:

```math
\sum_{\text{cyclic}} \int c(u) \left[ c'(u) [A_u, B_u], C_u \right] .
```

Normalised by the integral of the absolute value of the first term's integrand, since for
particular test fields the terms can individually vanish and a relative error against them
would be meaningless. Returns `(residual, scale)` when `normalised = false`.

# Why linear functionals are the sharp test

Taking ``A = \int \alpha u`` linear makes the second variations vanish identically, so the
cyclic sum above *is* the whole obstruction rather than a piece of it. The residual then
vanishes for every ``c`` — the bracket is Poisson for any structure function whatever, which
is the content of the theorem and the reason the Casimir family of
[`plucker_residual`](@ref) is as large as it is. `c(u) = u`, ``u^3``, ``\exp u``,
``\sin u`` and ``\tfrac{1}{2}(1 + \log u)^2`` are all checked in
`scripts/verify_fourbracket_identities.jl`.

Distinct from [`jacobi_residual(::DiscreteBracket, ::AbstractVector)`](@ref), which asks the
same question of a discrete structure matrix; here there is no matrix, only fields.
"""
function jacobi_residual(g::TorusGrid, a, b, c, cu, dcu; normalised::Bool = true)
    term(p, q, r) = integrate(g, cu .*
                                 canonical_bracket(g, dcu .* canonical_bracket(g, p, q), r))
    res = abs(term(a, b, c) + term(b, c, a) + term(c, a, b))
    scale = integrate(g, abs.(cu .*
                              canonical_bracket(g, dcu .* canonical_bracket(g, a, b), c)))
    normalised || return (res, scale)
    return iszero(scale) ? zero(res) : res / scale
end
