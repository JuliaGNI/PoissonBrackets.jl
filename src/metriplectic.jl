#
# Section 2 of `poisson-brackets-from-four-brackets.tex`: the metriplectic four-bracket built
# from the Kulkarni-Nomizu product of two symmetric tensors.
#
# Finite-dimensional and purely algebraic, unlike the rest of that manuscript's material —
# no grid and no fields, just two symmetric matrices and two vectors. Kept in its own file
# for that reason, and because the documentation groups reference material by source file.
#

@doc raw"""
    kulkarni_nomizu(σ, μ) -> R

The Kulkarni-Nomizu product of two symmetric matrices, as the four-index array

```math
R^{ijkl} = \sigma^{ik}\mu^{jl} + \sigma^{jl}\mu^{ik}
         - \sigma^{il}\mu^{jk} - \sigma^{jk}\mu^{il} .
```

This is the tensor whose contraction against ``\alpha \otimes h \otimes \alpha \otimes h``
is the metriplectic two-bracket of equation (2.6) — see
[`metriplectic_bracket`](@ref), which evaluates that contraction in ``O(n^2)`` without ever
forming `R`.

Building `R` explicitly costs ``O(n^4)`` and is not how the bracket should be evaluated. It
exists so that the closed form can be checked against the definition it is supposed to
implement, which is what `test/metriplectic_tests.jl` does; and because the symmetries of
`R` — antisymmetry in ``i \leftrightarrow j``, in ``k \leftrightarrow l``, and symmetry
under exchanging the pairs — are properties of the array rather than of the contraction, so
they can only be inspected here.
"""
function kulkarni_nomizu(σ::AbstractMatrix, μ::AbstractMatrix)
    n = LinearAlgebra.checksquare(σ)
    n == LinearAlgebra.checksquare(μ) || throw(DimensionMismatch(
        "the two tensors have sizes $(size(σ)) and $(size(μ))"))
    T = promote_type(eltype(σ), eltype(μ))
    R = zeros(T, n, n, n, n)
    @inbounds for l in 1:n, k in 1:n, j in 1:n, i in 1:n
        R[i, j, k, l] = σ[i, k] * μ[j, l] + σ[j, l] * μ[i, k] -
                        σ[i, l] * μ[j, k] - σ[j, k] * μ[i, l]
    end
    return R
end

@doc raw"""
    metriplectic_bracket(σ, μ, α, h)
    metriplectic_bracket(R, α, h)

The metriplectic two-bracket of equation (2.6), ``(A, A)_H`` for the four-bracket generated
by the Kulkarni-Nomizu product of `σ` and `μ`:

```math
(A, A)_H = (\alpha^\top \sigma \alpha)(h^\top \mu h)
         + (\alpha^\top \mu \alpha)(h^\top \sigma h)
         - 2 (\alpha^\top \sigma h)(\alpha^\top \mu h) ,
```

where `α` is the gradient of the observable and `h` that of the Hamiltonian.

The four-argument form evaluates the closed expression above in ``O(n^2)``; the
three-argument form contracts an explicitly built tensor from
[`kulkarni_nomizu`](@ref) instead, at ``O(n^4)``, and agrees with it identically.

# Proposition 2.1, and why the hypothesis is not decoration

If `σ` and `μ` are symmetric **positive semi-definite** then this quantity is non-negative,
so the bracket is dissipative in the direction it must be. The paper proves it by
Cauchy-Schwarz and AM-GM; `scripts/verify_fourbracket_metriplectic.jl` cross-checks it by
Monte-Carlo, deliberately including singular `σ` and `μ` of every rank pair, since it is the
rank-deficient case that an inequality chain is most likely to lose.

Semi-definiteness is genuinely needed rather than convenient: drop it for a merely symmetric
indefinite `σ` and negative values appear within a few hundred random draws. The same script
demonstrates that, so that the proposition is known not to be vacuous.
"""
function metriplectic_bracket(σ::AbstractMatrix, μ::AbstractMatrix,
                              α::AbstractVector, h::AbstractVector)
    return dot(α, σ, α) * dot(h, μ, h) +
           dot(α, μ, α) * dot(h, σ, h) -
           2 * dot(α, σ, h) * dot(α, μ, h)
end

function metriplectic_bracket(R::AbstractArray{T,4}, α::AbstractVector,
                              h::AbstractVector) where {T}
    n = size(R, 1)
    length(α) == n && length(h) == n || throw(DimensionMismatch(
        "the tensor is $(size(R)) but the vectors have lengths $(length(α)) and $(length(h))"))
    s = zero(promote_type(T, eltype(α), eltype(h)))
    @inbounds for l in 1:n, k in 1:n, j in 1:n, i in 1:n
        s += R[i, j, k, l] * α[i] * h[j] * α[k] * h[l]
    end
    return s
end
