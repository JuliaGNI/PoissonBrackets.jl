
#
# The broken (discontinuous Galerkin) hierarchical space, and the direct discretisation of the
# Burgers bracket on it.
#
# This is the fourth realisation of the Dirac companion note. The continuous degree-p Lagrange
# space V1 sits inside the broken degree-(2p-1) space, and V2 is taken as its L2-orthogonal
# complement there. The point of the construction is that [V1, V1] is contained in V1 + V2
# EXACTLY -- a product of two degree-p polynomials with one derivative has degree 2p-1 --
# which is the closure premise that a continuous finite element basis cannot supply, since
# phi' jumps across element interfaces.
#
# Everything here is exact over Q. The shifted Legendre polynomials L~_k(x) = P_k(2x - 1) have
# INTEGER coefficients,
#
#     L~_k(x) = sum_{j=0}^{k} (-1)^{k+j} binom(k,j) binom(k+j,j) x^j ,
#
# and int_0^1 x^n dx = 1/(n+1), so every integral in the assembly is rational arithmetic on
# coefficient vectors. The Python prototype reached for SymPy here, and for the Lagrange shape
# functions went through `float64` and `limit_denominator`; neither is needed, and dropping
# them is what keeps a symbolic dependency out of this package altogether.
#
# Polynomials are coefficient vectors in ASCENDING powers throughout.
#
# This was a file-level `@doc` block, which Julia attached to `_polymul` below and Documenter
# then dropped: none of it reached the manual. The prose now lives in `docs/src/dirac.md`.

_polymul(a::AbstractVector{T}, b::AbstractVector{T}) where {T} = begin
    c = zeros(T, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    c
end

"a - b, padded to the longer of the two; the products below differ in degree."
function _polysub(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
    n = max(length(a), length(b))
    c = zeros(T, n)
    c[eachindex(a)] .= a
    c[eachindex(b)] .-= b
    c
end

_polyder(a::AbstractVector{T}) where {T} =
    length(a) ≤ 1 ? zeros(T, 1) : T[(k - 1) * a[k] for k in 2:length(a)]

"∫₀¹ of a polynomial given by its coefficients in ascending powers."
_polyint01(a::AbstractVector{T}) where {T} = sum(a[k] // k for k in eachindex(a); init = zero(T))

@doc raw"""
    shifted_legendre([T = Rational{BigInt}], kmax) -> Vector{Vector{T}}

The shifted Legendre polynomials ``\tilde{L}_0, \dots, \tilde{L}_{k_{\max}}`` on ``[0,1]``,
as coefficient vectors in ascending powers.

Orthogonal with ``\int_0^1 \tilde{L}_a \tilde{L}_b = \delta_{ab}/(2a+1)``, which is what
makes the local mass matrix of [`dg_local_algebra`](@ref) diagonal.
"""
function shifted_legendre(::Type{T}, kmax::Int) where {T}
    [T[(-1)^(k + j) * binomial(big(k), big(j)) * binomial(big(k + j), big(j)) for j in 0:k]
     for k in 0:kmax]
end

shifted_legendre(kmax::Int) = shifted_legendre(Rational{BigInt}, kmax)

@doc raw"""
    lagrange_coefficients([T = Rational{BigInt}], p) -> Vector{Vector{T}}

The degree-``p`` Lagrange shape functions at the ``p+1`` equispaced nodes ``a/p`` of
``[0,1]``, as exact coefficient vectors.

Built by multiplying out ``\prod_{b \ne a}(\xi - \xi_b)/(\xi_a - \xi_b)`` over
``\mathbb{Q}``, so the nodal property holds identically rather than to fifteen digits.
"""
function lagrange_coefficients(::Type{T}, p::Int) where {T}
    p ≥ 1 || throw(ArgumentError("lagrange_coefficients needs p ≥ 1, got $p"))
    ξ = [T(a) // T(p) for a in 0:p]
    map(1:(p + 1)) do a
        c = T[one(T)]
        for b in 1:(p + 1)
            b == a && continue
            c = _polymul(c, T[-ξ[b], one(T)]) ./ (ξ[a] - ξ[b])
        end
        c
    end
end

lagrange_coefficients(p::Int) = lagrange_coefficients(Rational{BigInt}, p)

@doc raw"""
    dg_local_algebra([T = Rational{BigInt}], p) -> (Mloc, Tloc)

The local mass matrix and commutator tensor on one element, in the shifted Legendre basis of
the degree-``2p-1`` broken space:

```math
M^{\mathrm{loc}}_{ab} = \int_0^1 \tilde{L}_a \tilde{L}_b , \qquad
T^{\mathrm{loc}}_{mab} = \int_0^1 \tilde{L}_m \left( \tilde{L}_a \tilde{L}_b'
                                                   - \tilde{L}_b \tilde{L}_a' \right) .
```

The element size ``h`` drops out of ``T^{\mathrm{loc}}`` — the physical derivative carries
``1/h`` and the physical integral carries ``h`` — so only ``M^{\mathrm{loc}}`` scales with
it. That is why the assembly below can set ``h = 1`` and lose nothing.
"""
function dg_local_algebra(::Type{T}, p::Int) where {T}
    n = 2p
    L = shifted_legendre(T, n - 1)
    dL = _polyder.(L)
    M = [_polyint01(_polymul(L[a], L[b])) for a in 1:n, b in 1:n]
    Tt = Array{T,3}(undef, n, n, n)
    for m in 1:n, a in 1:n, b in 1:n
        Tt[m, a, b] = _polyint01(_polymul(L[m],
            _polysub(_polymul(L[a], dL[b]), _polymul(L[b], dL[a]))))
    end
    return M, Tt
end

dg_local_algebra(p::Int) = dg_local_algebra(Rational{BigInt}, p)

@doc raw"""
    coarse_in_dg([T = Rational{BigInt}], p, ne) -> Q₁

`Q₁[i, A]`: the periodic continuous degree-``p`` Lagrange basis expanded in the global
broken shifted-Legendre basis, exactly — the continuous space *is* contained in the broken
one, so the expansion terminates.

`ne` elements give ``\dim V_1 = p \cdot n_e`` and ``\dim V_{\mathrm{DG}} = 2p \cdot n_e``.
The coefficients come from orthogonality, ``\varphi_a = \sum_k (2k+1)
\left(\int_0^1 \varphi_a \tilde{L}_k\right) \tilde{L}_k``.
"""
function coarse_in_dg(::Type{T}, p::Int, ne::Int) where {T}
    L = shifted_legendre(T, 2p - 1)
    φ = lagrange_coefficients(T, p)
    coef = [[(2k - 1) * _polyint01(_polymul(φ[a], L[k])) for k in 1:2p] for a in 1:(p + 1)]

    N1, Ndg = p * ne, 2p * ne
    Q = zeros(T, N1, Ndg)
    for e in 0:(ne - 1), a in 1:(p + 1)
        i = mod(e * p + (a - 1), N1) + 1          # periodic wrap of the shared end node
        for k in 1:2p
            Q[i, e * 2p + k] += coef[a][k]
        end
    end
    return Q
end

coarse_in_dg(p::Int, ne::Int) = coarse_in_dg(Rational{BigInt}, p, ne)

@doc raw"""
    assemble_dg_hierarchical([T = Rational{BigInt}], p, ne) -> (; C, keep, con, M)

Structure constants of the Burgers bracket on the broken degree-``2p-1`` space, in a basis
adapted to the split ``V_1`` (continuous, degree ``p``) ``\oplus`` ``V_2 = V_1^\perp``.

`C[m,i,j]` are exact rationals, so ``\mathbb{J}_{ij}(u) = \sum_m c_{ij}^m u_m`` is the direct
discretisation of section 4 of the notes for this basis, and `keep`/`con` index the two
blocks ready for the Dirac reduction of `src/dirac.jl`.

``V_2`` is the ``L^2``-orthogonal complement, obtained as `kernel(Q₁ M_DG)` — exactly, where
the floating-point prototype needed an SVD. The transformed commutator tensor is contracted
one index at a time; the direct sextuple sum of the prototype is ``O(N_{DG}^6)`` and
needlessly so.
"""
function assemble_dg_hierarchical(::Type{T}, p::Int, ne::Int) where {T}
    Mloc, Tloc = dg_local_algebra(T, p)
    nl, Ndg = 2p, 2p * ne

    Mdg = zeros(T, Ndg, Ndg)
    Tdg = zeros(T, Ndg, Ndg, Ndg)
    for e in 0:(ne - 1)
        o = e * nl
        Mdg[(o + 1):(o + nl), (o + 1):(o + nl)] = Mloc          # h = 1; h scales M only
        Tdg[(o + 1):(o + nl), (o + 1):(o + nl), (o + 1):(o + nl)] = Tloc
    end

    Q1 = coarse_in_dg(T, p, ne)
    N1 = size(Q1, 1)
    Q2 = transpose(kernel(Q1 * Mdg))                            # V₂ = V₁^⊥ in the DG space
    Q  = vcat(Q1, Q2)
    n  = size(Q, 1)

    M = Q * Mdg * transpose(Q)

    # contract Tdg against Q on each slot in turn: O(n N³) rather than O(n³ N³)
    S1 = zeros(T, n, Ndg, Ndg)
    for i in 1:n, a in 1:Ndg
        iszero(Q[i, a]) && continue
        @views S1[i, :, :] .+= Q[i, a] .* Tdg[a, :, :]
    end
    S2 = zeros(T, n, n, Ndg)
    for j in 1:n, b in 1:Ndg
        iszero(Q[j, b]) && continue
        @views S2[:, j, :] .+= Q[j, b] .* S1[:, b, :]
    end
    T3 = zeros(T, n, n, n)
    for k in 1:n, c in 1:Ndg
        iszero(Q[k, c]) && continue
        @views T3[:, :, k] .+= Q[k, c] .* S2[:, :, c]
    end

    C = galerkin_c(inv(M), T3)
    return (; C, keep = collect(1:N1), con = collect((N1 + 1):n), M)
end

assemble_dg_hierarchical(p::Int, ne::Int) = assemble_dg_hierarchical(Rational{BigInt}, p, ne)
