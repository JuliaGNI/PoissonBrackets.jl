
#
# Finite-dimensional Lie algebras, as the structure constants of a Lie-Poisson bracket.
#
# Throughout this file a three-tensor `C` holds c_ij^m at `C[m, i, j]`, so that
#
#     [e_i, e_j] = sum_m c_ij^m e_m ,      J_ij(z) = sum_m c_ij^m z_m .
#
# The leading-index convention is the one `structure_constant_residual` already uses, and is
# the reason it is `C[m,i,j]` rather than the `c[i][j][m]` of the Python prototypes:
# `C[m, :, :]` is then a matrix, and the whole tensor contracts against a coefficient vector
# with a single `lie_poisson_matrix` call.
#
# The mode algebras at the end -- Witt, one-sided Witt, torus -- return a named tuple
# `(; C, keep, con, labels)` splitting the generators into a coarse block V1 = keep and a
# constrained block V2 = con, which is what the Dirac reduction of `src/dirac.jl` consumes.
#
# This was a file-level `@doc` block attached to no binding, which Documenter drops without
# a word: none of it reached the manual. The prose now lives in `docs/src/liepoisson.md`,
# where it is rendered; the comment stays so that a reader of the source is not sent
# elsewhere for the index convention.

@doc raw"""
    lie_poisson_matrix(C, z)

The Lie-Poisson structure matrix ``\mathbb{J}_{ij} = \sum_m c_{ij}^m z_m`` of the structure
constants `C[m,i,j]` at the point `z`.

This bracket is Poisson exactly when the `C` are the structure constants of a Lie algebra,
which [`structure_constant_residual`](@ref) measures; the two are the same statement seen
from either end.
"""
function lie_poisson_matrix(C::AbstractArray{T, 3}, z::AbstractVector{S}) where {T, S}
    N = size(C, 2)
    R = promote_type(T, S)
    J = zeros(R, N, N)
    @inbounds for j in 1:N, i in 1:N

        s = zero(R)
        for m in eachindex(z)
            s += C[m, i, j] * z[m]
        end
        J[i, j] = s
    end
    return J
end

"""
    lie_poisson_derivative(C)

The derivative tensor `dJ[m,i,j]` ``= \\partial \\mathbb{J}_{ij} / \\partial z_m`` of
[`lie_poisson_matrix`](@ref), which for a bracket linear in the field is `C` itself.

Provided under its own name so that a caller assembling a Jacobiator does not have to know
that the two coincide.
"""
lie_poisson_derivative(C::AbstractArray{T, 3}) where {T} = C

const _EPS3 = ((1, 2, 3, 1), (2, 3, 1, 1), (3, 1, 2, 1),
    (3, 2, 1, -1), (2, 1, 3, -1), (1, 3, 2, -1))

@doc raw"""
    so3([T = Rational{BigInt}], n = 3)

``\mathfrak{so}(3)``, optionally padded with `n - 3` abelian directions.

!!! warning "A degenerate control"
    Do not validate a Jacobi or structure-constant routine against this algebra alone. It
    lies inside the six-parameter family ``c_{ij}^k = \epsilon_{ijl} n^{lk}`` with ``n``
    symmetric — Bianchi class A — every member of which satisfies the Jacobi identity, and
    the perturbations one reaches for stay inside it: rescaling a generator, or rescaling a
    single structure constant, both keep ``n`` symmetric and diagonal. So
    ``\mathfrak{so}(3)`` goes on passing after it has apparently been broken.

    Antisymmetry alone does **not** force Jacobi in three dimensions — a general
    antisymmetric `c` has nine parameters against this family's six, and random ones fail
    comfortably. Use [`se3`](@ref) as the positive control and
    [`random_antisymmetric_c`](@ref) in dimension five as the negative one.
"""
function so3(::Type{T} = Rational{BigInt}, n::Int = 3) where {T}
    n ≥ 3 || throw(ArgumentError("so3 needs n ≥ 3, got $n"))
    C = zeros(T, n, n, n)
    for (i, j, k, v) in _EPS3
        C[k, i, j] = T(v)
    end
    return C
end

@doc raw"""
    se3([T = Rational{BigInt}])

``\mathfrak{se}(3)``, the Euclidean algebra, as a six-dimensional structure-constant tensor:

```math
[J_i, J_j] = \epsilon_{ijk} J_k , \qquad
[J_i, P_j] = \epsilon_{ijk} P_k , \qquad
[P_i, P_j] = 0 .
```

Indices `1:3` are the rotations and `4:6` the translations. This is the discriminating
positive control — it is semidirect rather than simple, and it is large enough that the
Jacobi identity is a real constraint rather than an automatic consequence of antisymmetry.
"""
function se3(::Type{T} = Rational{BigInt}) where {T}
    C = zeros(T, 6, 6, 6)
    for (i, j, k, v) in _EPS3
        C[k, i, j] += T(v)
        C[3 + k, i, 3 + j] += T(v)
        C[3 + k, 3 + j, i] -= T(v)
    end
    return C
end

@doc raw"""
    so_n([T = Rational{BigInt}], N) -> (basis, C)

``\mathfrak{so}(N)`` in the basis ``E_{pq} - E_{qp}`` for ``p < q``, returning both the
matrices and their structure constants `C[m,i,j]`.

The constants are read off with the trace form ``\langle X, Y \rangle = -\tfrac12
\operatorname{tr}(XY)``, in which this basis is orthonormal. Used for the spectral Casimirs:
``\operatorname{tr}(W^k)``, and hence ``\operatorname{tr}(\eta(W))`` for any analytic
``\eta``, is a Casimir of the Lie-Poisson bracket, which is the finite-dimensional shadow of
prescribing the Casimir density in the continuum.
"""
function so_n(::Type{T}, N::Int) where {T}
    idx = [(p, q) for p in 1:N for q in (p + 1):N]
    basis = map(idx) do (p, q)
        E = zeros(T, N, N)
        E[p, q], E[q, p] = one(T), -one(T)
        E
    end
    d = length(basis)
    C = zeros(T, d, d, d)
    for j in 1:d, i in 1:d

        comm = basis[i] * basis[j] - basis[j] * basis[i]
        for m in 1:d
            C[m, i, j] = -tr(comm * basis[m]) // 2
        end
    end
    return basis, C
end

so_n(N::Int) = so_n(Rational{BigInt}, N)

"""
    random_antisymmetric_c([rng], [T = Rational{BigInt}], n)

Random `C[m,i,j]` with ``c_{ij}^m = -c_{ji}^m``, which is generically **not** a Lie algebra.

The negative control. Entries are rationals with numerator drawn from `-6:6` and
denominator from `1:4`, matching the Python prototypes closely enough that the residuals
land in the same range — though not the same numbers, since the two random streams differ.
Take `n ≥ 5`: below that, antisymmetry alone is close enough to forcing Jacobi that the
control stops discriminating.
"""
function random_antisymmetric_c(rng::AbstractRNG, ::Type{T}, n::Int) where {T}
    C = zeros(T, n, n, n)
    for i in 1:n, j in (i + 1):n, m in 1:n
        v = T(rand(rng, -6:6)) / T(rand(rng, 1:4))
        C[m, i, j] = v
        C[m, j, i] = -v
    end
    return C
end

function random_antisymmetric_c(rng::AbstractRNG, n::Int)
    random_antisymmetric_c(rng, Rational{BigInt}, n)
end
random_antisymmetric_c(n::Int) = random_antisymmetric_c(Random.default_rng(), n)

@doc raw"""
    sine_algebra(N) -> (modes, C)

The ``\mathfrak{su}(N)`` sine algebra of Zeitlin's truncation, built from clock and shift.

With ``\zeta = e^{4\pi i/N}`` primitive for odd `N`, ``g = \operatorname{diag}(\zeta^a)``,
``h`` the cyclic shift and

```math
T_m = e^{2\pi i m_1 m_2 / N} \, g^{m_1} h^{m_2} , \qquad m \in \mathbb{Z}_N^2 \setminus \{0\} ,
```

the commutator closes on the ``N^2 - 1 = \dim\mathfrak{su}(N)`` generators,

```math
[T_m, T_n] = -2i \sin\!\left( \tfrac{2\pi}{N} \, m \times n \right) T_{m+n} ,
```

and the returned `C` is normalised so that its entries tend to ``m \times n`` as
``N \to \infty`` — that is, to the structure constants of the torus algebra of
[`torus_truncation`](@ref), at second order in ``1/N``.

`N` must be odd: the half-integer power ``\zeta^{1/2} = e^{2\pi i/N}`` implicit in the
prefactor is single-valued on ``\mathbb{Z}_N^2`` only then.

This is the one discretisation in these notes that closes into a Lie algebra **exactly** at
finite `N`, which is what makes it the useful comparison for the finite element brackets
that do not.
"""
function sine_algebra(N::Int)
    isodd(N) ||
        throw(ArgumentError("sine_algebra needs N odd for ζ = exp(4πi/N) to be primitive, got N = $N"))
    ζ = cis(4π / N)
    g = Diagonal([ζ^a for a in 0:(N - 1)])
    h = zeros(ComplexF64, N, N)
    for i in 1:N
        h[i, mod1(i + 1, N)] = 1
    end

    modes = [(m1, m2) for m1 in 0:(N - 1) for m2 in 0:(N - 1) if (m1, m2) != (0, 0)]
    pos = Dict(m => i for (i, m) in enumerate(modes))
    T = Dict(m => cis(2π * m[1] * m[2] / N) * (g^m[1]) * (h^m[2]) for m in modes)

    d = length(modes)
    λ = im * N / (4π)                       # normalisation: coefficients → m × n
    C = zeros(ComplexF64, d, d, d)
    for m in modes, n in modes

        s = (mod(m[1] + n[1], N), mod(m[2] + n[2], N))
        s == (0, 0) && continue             # the commutator is central there
        comm = T[m] * T[n] - T[n] * T[m]
        # expand in the T basis, which is orthogonal: tr(T_s' T_r) = N δ_sr
        C[pos[s], pos[m], pos[n]] = λ * tr(T[s]' * comm) / N
    end
    return modes, C
end

@doc raw"""
    sine_coefficient(N, m, n) -> (s, c)

The single structure constant of [`sine_algebra`](@ref) at the mode pair `m`, `n`: the
target mode ``s = m + n \bmod N`` and

```math
c = \frac{N}{2\pi} \sin\!\left( \frac{2\pi}{N} \, m \times n \right) ,
\qquad m \times n = m_1 n_2 - m_2 n_1 ,
```

in the same normalisation. Returns `(s, 0.0)` when the commutator is central, ``s = 0``.

For each pair of modes the sine algebra has exactly **one** nonzero structure constant, so
the dense ``d \times d \times d`` tensor `sine_algebra` returns holds ``d^2`` numbers in
``d^3`` slots. That is affordable for the closure check at small `N` and not affordable at
all beyond it — ``N = 41`` would ask for 75 GB. Use this when only the coefficients are
wanted, as in the ``N^{-2}`` convergence to the torus algebra of [`torus_truncation`](@ref):
expanding the sine gives

```math
c - m \times n = -\frac{(2\pi)^2}{6 N^2} (m \times n)^3 + O(N^{-4}) ,
```

so the deviation is second order in ``1/N`` *at fixed mode*, which is why the comparison has
to be made on a window of low modes held fixed as `N` grows. Taken over all modes the sine
saturates and nothing converges.
"""
function sine_coefficient(N::Int, m::Tuple{Int, Int}, n::Tuple{Int, Int})
    s = (mod(m[1] + n[1], N), mod(m[2] + n[2], N))
    s == (0, 0) && return s, 0.0
    cross = m[1] * n[2] - m[2] * n[1]
    return s, (N / 2π) * sin(2π * cross / N)
end

@doc raw"""
    witt_truncation([T = Rational{BigInt}], K; q = 2, drop_zero = false)

The Witt algebra ``[L_k, L_l] = (l - k) L_{k+l}`` truncated to ``|k| \le qK``, returned as
`(; C, keep, con, labels)`.

This is the Fourier realisation of the Burgers / ``\mathrm{Vect}(S^1)`` bracket: with
``e_k = e^{ikx}`` and ``[F,G] = FG' - GF'`` one has ``[e_k, e_l] = i(l-k)e_{k+l}``. The
global ``i`` is dropped, which rescales the bracket by a constant and so leaves the Jacobi
identity alone. The ``u_k`` are then independent formal coordinates: Jacobi is a polynomial
identity in them and holds over ``\mathbb{R}`` iff it holds over the complexification, so
the reality constraint ``u_{-k} = \bar{u}_k`` need not be imposed.

`keep` indexes the coarse block ``V_1 = \{|k| \le K\}`` and `con` the constrained block
``V_2 = \{K < |k| \le qK\}``. For ``k, l \in V_1`` we have ``|k+l| \le 2K``, so the coarse
block closes exactly once `q ≥ 2` — the premise the hierarchical-basis idea needs, and one
that a ``C^0`` finite element basis cannot supply.
"""
function witt_truncation(::Type{T}, K::Int; q::Int = 2, drop_zero::Bool = false) where {T}
    M = q * K
    modes = [k for k in (-M):M if !(drop_zero && k == 0)]
    pos = Dict(k => i for (i, k) in enumerate(modes))
    n = length(modes)
    C = zeros(T, n, n, n)
    for (a, k) in enumerate(modes), (b, l) in enumerate(modes)

        haskey(pos, k + l) && (C[pos[k + l], a, b] = T(l - k))
    end
    keep = [pos[k] for k in modes if abs(k) ≤ K]
    con = [pos[k] for k in modes if abs(k) > K]
    return (; C, keep, con, labels = ["L($k)" for k in modes])
end

witt_truncation(K::Int; kwargs...) = witt_truncation(Rational{BigInt}, K; kwargs...)

@doc raw"""
    graded_witt([T = Rational{BigInt}], lo, D) -> (; C, degs)

``\mathrm{span}\{L_a : lo \le a \le D\}`` with ``[L_a, L_b] = (b-a)L_{a+b}``, dropping
whatever leaves the window.

The general truncated-Witt builder, with the cut-offs on both sides as free parameters;
[`witt_truncation`](@ref) is the symmetric case and [`poly_truncation`](@ref) the case
``lo = -1``. Which is worth having separately because the *position* of the lower cut-off is
exactly what decides the tension of the Dirac note:

- `lo ≥ 0` gives a non-negative additive grading, so ``\{\deg > p\}`` is an **ideal**, the
  truncation closes and is a genuine Lie algebra — but then ``a, b \in V_2`` forces
  ``\deg(a+b) > \max(\deg a, \deg b)``, so ``[V_2,V_2]`` never reaches ``V_1``, the
  constraint matrix ``C`` vanishes identically and the constraints are **first** class. Dirac
  reduction is unavailable, and unnecessary: the plain restriction is already Poisson.
- `lo < 0` includes ``L_{-1} = \partial_x`` and destroys the grading. Now ``C`` is invertible
  and the constraints are second class — but the truncation is no longer a Lie algebra.

Closure and second-classness exclude each other, and this family is where that is cleanest.
"""
function graded_witt(::Type{T}, lo::Int, D::Int) where {T}
    degs = collect(lo:D)
    pos = Dict(a => i for (i, a) in enumerate(degs))
    n = length(degs)
    C = zeros(T, n, n, n)
    for (i, a) in enumerate(degs), (j, b) in enumerate(degs)

        haskey(pos, a + b) && (C[pos[a + b], i, j] = T(b - a))
    end
    return (; C, degs)
end

graded_witt(lo::Int, D::Int) = graded_witt(Rational{BigInt}, lo, D)

@doc raw"""
    poly_truncation([T = Rational{BigInt}], p; q = 2)

The one-sided Witt algebra, the global-polynomial realisation, as `(; C, keep, con, labels)`.

For ``[F,G] = FG' - GF'`` on monomials, ``[x^a, x^b] = (b-a)x^{a+b-1}``, so with
``L_a = x^{a+1}`` one recovers ``[L_a, L_b] = (b-a)L_{a+b}`` with ``a \ge -1``. `keep` is the
polynomials of degree ``\le p``, `con` those of degree ``p+1`` to ``qp``.

The one-sidedness matters. For ``a \ge 0`` the high modes *are* an ideal, since
``\deg(a+b) \ge \max(\deg a, \deg b)``, and the naive truncation is already a nilpotent Lie
algebra; including ``L_{-1} = \partial_x`` destroys that, because ``[L_{-1}, L_{D+1}] =
(D+2)L_D`` re-enters from above. The truncation returned here **does** include ``L_{-1}``.
"""
function poly_truncation(::Type{T}, p::Int; q::Int = 2) where {T}
    M = q * p
    degs = collect(-1:(M - 1))              # L_a ↔ x^(a+1), of degree a+1
    pos = Dict(a => i for (i, a) in enumerate(degs))
    n = length(degs)
    C = zeros(T, n, n, n)
    for (i, a) in enumerate(degs), (j, b) in enumerate(degs)

        haskey(pos, a + b) && (C[pos[a + b], i, j] = T(b - a))
    end
    keep = [pos[a] for a in degs if a + 1 ≤ p]
    con = [pos[a] for a in degs if a + 1 > p]
    return (; C, keep, con, labels = ["x^$(a + 1)" for a in degs])
end

poly_truncation(p::Int; kwargs...) = poly_truncation(Rational{BigInt}, p; kwargs...)

@doc raw"""
    torus_truncation([T = Rational{BigInt}], K; q = 2)

The 2D vorticity / Vlasov bracket on ``T^2`` in the Fourier basis, as
`(; C, keep, con, labels)`.

With ``e_m = e^{i m \cdot x}`` and ``[F,G] = F_x G_y - F_y G_x``,
``[e_m, e_n] = -(m \times n) e_{m+n}`` where ``m \times n = m_1 n_2 - m_2 n_1``; the overall
sign is dropped as in [`witt_truncation`](@ref). `keep` is
``\{|m|_\infty \le K\} \setminus \{0\}`` and `con` the shell ``K < |m|_\infty \le qK``. The
coarse block closes exactly for ``q \ge 2``: sums stay within ``|m|_\infty \le 2K``, and the
``m + n = 0`` term carries the coefficient ``m \times (-m) = 0``.

This is the same Lie algebra ``C^\infty(T^2)`` that Zeitlin's truncation quantises, so it is
the case where [`sine_algebra`](@ref) provides a genuinely closing comparison.
"""
function torus_truncation(::Type{T}, K::Int; q::Int = 2) where {T}
    M = q * K
    modes = [(m1, m2) for m1 in (-M):M for m2 in (-M):M if (m1, m2) != (0, 0)]
    pos = Dict(m => i for (i, m) in enumerate(modes))
    n = length(modes)
    C = zeros(T, n, n, n)
    for (i, m) in enumerate(modes), (j, nn) in enumerate(modes)

        s = (m[1] + nn[1], m[2] + nn[2])
        haskey(pos, s) && (C[pos[s], i, j] = T(m[1] * nn[2] - m[2] * nn[1]))
    end
    sup(m) = max(abs(m[1]), abs(m[2]))
    keep = [pos[m] for m in modes if sup(m) ≤ K]
    con = [pos[m] for m in modes if sup(m) > K]
    return (; C, keep, con, labels = ["e$(m)" for m in modes])
end

torus_truncation(K::Int; kwargs...) = torus_truncation(Rational{BigInt}, K; kwargs...)

"""
    closes_on(C, keep, con)

Whether ``[V_1, V_1] \\subseteq V_1 + V_2``, i.e. whether the coarse block closes on the
retained and constrained generators together with nothing left over.

This is the premise of the hierarchical-basis construction, and [`leak`](@ref) measures how
badly it fails when it does.
"""
function closes_on(C::AbstractArray{T, 3}, keep, con) where {T}
    inside = Set{Int}(keep) ∪ Set{Int}(con)
    all(iszero(C[m, i, j]) for i in keep, j in keep for m in axes(C, 1) if m ∉ inside)
end

"""
    leak(C, keep, con)

``\\max |c_{ij}^m|`` over ``i, j \\in V_1`` and ``m`` outside ``V_1 + V_2`` — the closure
defect of [`closes_on`](@ref). Zero exactly when that returns `true`.
"""
function leak(C::AbstractArray{T, 3}, keep, con) where {T}
    inside = Set{Int}(keep) ∪ Set{Int}(con)
    worst = zero(real(T))
    for i in keep, j in keep, m in axes(C, 1)
        m ∈ inside && continue
        worst = max(worst, abs(C[m, i, j]))
    end
    return worst
end

@doc raw"""
    galerkin_c(Minv, T3)

The Galerkin structure constants of a general, non-nodal basis,

```math
c_{ij}^m = \sum_{pq} (M^{-1})_{ip} (M^{-1})_{jq} T_{mpq} , \qquad
T_{mpq} = \int \varphi_m \, [\varphi_p, \varphi_q] ,
```

returned in the `C[m,i,j]` convention of this file.

The direct discretisation of section 4 of the notes states ``c_{ij}^m = b_m \sum_{pq}
(M^{-1})_{ip} [\varphi_p, \varphi_q](x_m) (M^{-1})_{qj}``, which presumes a Lagrange basis
with quadrature collocated at its nodes — something a hierarchical basis does not have.
This formula is the one that survives without that assumption.

If the span is closed under the commutator, ``[\varphi_p, \varphi_q] = \sum_r g_{pq}^r
\varphi_r`` exactly, then ``T_{mpq} = \sum_r g_{pq}^r M_{mr}`` and `c` is the Lie-Poisson
tensor of ``g`` in the coordinates ``u = M^{-1}w`` — a *linear* change of coordinates, under
which Jacobi is invariant. So `c` is a Lie algebra iff `g` is. That change of coordinates is
not block diagonal unless ``M`` is, which is why the orthogonal-basis realisations are the
faithful ones and the finite element realisation has to go through this formula in the
primal coefficients.
"""
function galerkin_c(Minv::AbstractMatrix{S}, T3::AbstractArray{R, 3}) where {S, R}
    n = size(Minv, 1)
    U = promote_type(S, R)
    # Contract the two slots in turn. The direct double sum inside a triple loop is O(n^5),
    # which at the sizes the refinement study of the Dirac note needs (n = 48) is minutes of
    # exact rational arithmetic rather than seconds.
    S1 = zeros(U, n, n, n)                      # S1[m,i,q] = Σ_p Minv[i,p] T3[m,p,q]
    @inbounds for q in 1:n, p in 1:n

        for i in 1:n
            iszero(Minv[i, p]) && continue
            f = Minv[i, p]
            for m in 1:n
                S1[m, i, q] += f * T3[m, p, q]
            end
        end
    end
    C = zeros(U, n, n, n)                       # C[m,i,j] = Σ_q Minv[j,q] S1[m,i,q]
    @inbounds for q in 1:n, j in 1:n

        iszero(Minv[j, q]) && continue
        f = Minv[j, q]
        for i in 1:n, m in 1:n

            C[m, i, j] += f * S1[m, i, q]
        end
    end
    return C
end
