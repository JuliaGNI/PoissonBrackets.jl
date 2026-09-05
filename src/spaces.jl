
@doc raw"""
    DiscreteSpace{T}

A finite-dimensional space of functions on the periodic domain ``\Omega = [0,L)``, together
with the quadrature its assemblies are carried out with.

The whole package is written against one small interface, which every concrete space
answers:

  - `nbasis` — the number of degrees of freedom ``N``;
  - `quadrature_nodes`, `quadrature_weights` — the global rule, `x` and `w`;
  - `basis_values``(space, d)` — the table ``\Phi_d[i,r] = D^d \phi_i (x_r)``;
  - `mass_matrix` and its factorization.

Everything else — the discrete brackets, the Hamiltonians, their gradients and Hessians — is
a weighted contraction of those tables and is therefore written once, generically. A matrix
``\int_\Omega f(x) \, D^a \phi_k \, D^b \phi_l \, dx`` is

```math
\Phi_a \, \mathrm{diag}(f \odot w) \, \Phi_b^T ,
```

whatever the basis underneath happens to be.

Two spaces are provided, and the choice between them is not free: it is the discretisation
each paper actually uses.

  - [`SplineSpace`](@ref), a periodic B-spline basis, ``\mathcal{C}^{p-1}`` and hence able
    to carry the third derivative of the KdV brackets;
  - [`LagrangeSpace`](@ref), a periodic nodal Lagrange finite element space, only
    ``\mathcal{C}^0``, but *nodal*, which is what the pointwise transformation
    ``\bar{u}_i = \sqrt{u_i}`` of the Burgers bracket needs to be justified.
"""
abstract type DiscreteSpace{T} end

Base.eltype(::DiscreteSpace{T}) where {T} = T
Base.length(s::DiscreteSpace) = nbasis(s)

@doc raw"""
    inverse_mass_matrix(space)

The inverse mass matrix ``\mathbb{M}^{-1}``, formed once when the space is built.

Every discrete bracket here is a weak-form operator sandwiched between two of these, so it
appears in every assembly and in every Jacobian. Forming it on each call — even as one
``O(N^3)`` inversion — is a measurable fraction of a run at the step counts these
discretisations reach.
"""
function inverse_mass_matrix end

"""
    mixed_matrix(space, a, b)

The matrix ``\\int_\\Omega D^a \\phi_k \\, D^b \\phi_l \\, dx``.
"""
function mixed_matrix(s::DiscreteSpace, a::Integer, b::Integer)
    basis_values(s, a) * Diagonal(quadrature_weights(s)) * basis_values(s, b)'
end

"""
    weighted_matrix(space, f, a, b)

The matrix ``\\int_\\Omega f(x) \\, D^a \\phi_k \\, D^b \\phi_l \\, dx``, with `f` either a
function of the coordinate or a vector already sampled at `quadrature_nodes`.
"""
function weighted_matrix(s::DiscreteSpace, f, a::Integer, b::Integer)
    weighted_matrix(s, f.(quadrature_nodes(s)), a, b)
end

function weighted_matrix(s::DiscreteSpace, f::AbstractVector, a::Integer, b::Integer)
    length(f) == length(quadrature_weights(s)) || throw(DimensionMismatch(
        "the coefficient was sampled at $(length(f)) points but the quadrature has " *
        "$(length(quadrature_weights(s)))"))
    basis_values(s, a) * Diagonal(f .* quadrature_weights(s)) * basis_values(s, b)'
end

"""
    derivative_matrix(space)

The matrix ``S_{kl} = \\int_\\Omega \\phi_k \\phi_l' \\, dx``, antisymmetric on a periodic
domain.
"""
derivative_matrix(s::DiscreteSpace) = mixed_matrix(s, 0, 1)

"""
    stiffness_matrix(space)

The matrix ``\\mathbb{K}^1_{kl} = \\int_\\Omega \\phi_k' \\phi_l' \\, dx``.
"""
stiffness_matrix(s::DiscreteSpace) = mixed_matrix(s, 1, 1)

"""
    basis_integrals(space)

The vector ``\\int_\\Omega \\phi_i \\, dx``, which equals ``\\mathbb{M} \\mathbf{1}`` for a
basis that is a partition of unity.
"""
basis_integrals(s::DiscreteSpace) = basis_values(s, 0) * quadrature_weights(s)

"""
    field(space, û, d = 0)

The `d`-th derivative of ``u_h = \\sum_j \\hat{u}_j \\phi_j`` sampled at the quadrature
points — the one operation every assembly with a variable coefficient starts from.
"""
field(s::DiscreteSpace, û::AbstractVector, d::Integer = 0) = basis_values(s, d)' * û

"""
    project(space, f)

The coefficients of the ``L^2`` projection of `f` onto the space, with `f` a function or a
vector of values at `quadrature_nodes`.
"""
project(s::DiscreteSpace, f) = project(s, f.(quadrature_nodes(s)))

function project(s::DiscreteSpace, f::AbstractVector)
    mass_factorization(s) \ (basis_values(s, 0) * (quadrature_weights(s) .* f))
end

"""
    project!(û, space, f)

In-place [`project`](@ref). On a uniform mesh the solve goes through the planned transforms
of `CirculantMass` and allocates nothing beyond the right-hand side.
"""
function project!(û::AbstractVector, s::DiscreteSpace, f::AbstractVector)
    û .= basis_values(s, 0) * (quadrature_weights(s) .* f)
    mass_solve!(û, mass_factorization(s), û)
end

## Spline space

@doc raw"""
    SplineSpace(mesh, p; nq = quadrature_order(p), dmax = 3)
    SplineSpace(n, p; L = 2π, kwargs...)
    SplineSpace(quadrature::SplineQuadrature)

The periodic B-spline space of degree `p`, a thin wrapper around the assembly table of
SimpleSplines.

```jldoctest
julia> s = SplineSpace(16, 3);

julia> nbasis(s), degree(s)
(16, 3)
```

The basis is ``\mathcal{C}^{p-1}``, which is what lets the third-derivative term of the
second KdV bracket be carried at all: after one integration by parts it needs only
``p \ge 2``, and the two derivatives of ``\int \phi_k' \phi_l''`` then both exist.

`dmax = 3` tabulates derivatives up to third order, one more than the integrated-by-parts
form needs, so that the identity ``-\int \phi_k \phi_l''' = \int \phi_k' \phi_l''`` can be
checked rather than assumed.
"""
struct SplineSpace{T, QT <: SplineQuadrature{T}} <: DiscreteSpace{T}
    quadrature::QT
    Minv::Matrix{T}

    # an inner constructor, so that the one-argument outer form below is not also generated
    # by default and then overwritten
    # The dense inverse is built by solving against the identity rather than by inverting
    # the assembled matrix: the mass matrix is sparse now, and on a uniform mesh each of
    # those N solves goes through the planned transforms of `CirculantMass`.
    SplineSpace(q::QT) where {T, QT <: SplineQuadrature{T}} = new{T, QT}(
        q, mass_operator(q) \ Matrix{T}(I, nbasis(q), nbasis(q)))
end

function SplineSpace(mesh::Mesh, p::Integer; kwargs...)
    SplineSpace(SplineQuadrature(PeriodicBSplineBasis(mesh, p); kwargs...))
end

function SplineSpace(n::Integer, p::Integer; L = 2π, kwargs...)
    SplineSpace(UniformMesh(n, L), p; kwargs...)
end

quadrature(s::SplineSpace) = s.quadrature
basis(s::SplineSpace) = basis(s.quadrature)
nbasis(s::SplineSpace) = nbasis(s.quadrature)
degree(s::SplineSpace) = degree(s.quadrature)
order(s::SplineSpace) = order(s.quadrature)
ncells(s::SplineSpace) = ncells(s.quadrature)
domainlength(s::SplineSpace) = domainlength(s.quadrature)
nodes(s::SplineSpace) = nodes(basis(s.quadrature))
quadrature_nodes(s::SplineSpace) = quadrature_nodes(s.quadrature)
quadrature_weights(s::SplineSpace) = quadrature_weights(s.quadrature)
basis_values(s::SplineSpace, d::Integer = 0) = basis_values(s.quadrature, d)
mass_matrix(s::SplineSpace) = mass_matrix(s.quadrature)

# Delegated rather than recomputed from `basis_values`: the quadrature memoises these, and
# the stiffness matrix in particular would otherwise be reassembled on every Newton
# iteration of every step, which at N = 384 is the single largest cost in a run.
mixed_matrix(s::SplineSpace, a::Integer, b::Integer) = mixed_matrix(s.quadrature, a, b)
mass_factorization(s::SplineSpace) = mass_factorization(s.quadrature)
inverse_mass_matrix(s::SplineSpace) = s.Minv

"""
    evaluate(space, û, x, d = 0)

The `d`-th derivative of the discrete field with coefficients `û` at the point or points
`x`, which need not be quadrature points.
"""
function evaluate(s::SplineSpace, û::AbstractVector, x, d::Integer = 0)
    evaluate(basis(s.quadrature), û, x, d)
end

## Nodal Lagrange space

@doc raw"""
    LagrangeSpace(p, ne; L = 2π, nq = 2p + 4)

The periodic nodal Lagrange finite element space of degree `p` on a uniform mesh of `ne`
elements, with ``N = p \, n_e`` degrees of freedom.

```jldoctest
julia> s = LagrangeSpace(2, 12);

julia> nbasis(s), degree(s)
(24, 2)
```

The local basis is the equispaced Lagrange basis of `CompactBasisFunctions` on the reference
element ``[0,1]``, scattered with periodic connectivity: local node `p` of element `e` and
local node `0` of element `e+1` are the same degree of freedom, and the last element wraps
onto the first.

# Why a nodal basis

This is the space the Burgers bracket is discretised in, and the reason is the
transformation ``\bar{u} = \sqrt{u}`` that makes that bracket constant. On a nodal basis
the coefficients *are* the values of the field at the nodes, so
``\bar{u}_i = \sqrt{u_i}`` is the transformation applied coefficient by coefficient, and
the chain rule that turns the constant bracket ``\mathbb{K}`` back into
``\mathbb{J}_{ij} = \sqrt{u_i} \, \mathbb{K}_{ij} \sqrt{u_j}`` is exact. In a modal or
B-spline basis the same manipulation has no such justification.

The price is regularity: the basis is only ``\mathcal{C}^0``, so derivatives beyond the
first do not exist across element boundaries and only `d = 0, 1` are tabulated.

# Storage

``\Phi`` is stored **sparse**, as `SimpleSplines.SplineQuadrature` stores the spline
tabulation: a basis function is supported on at most two elements, so only ``p+1`` of the
``N`` rows are structurally nonzero in any quadrature column — 0.8 % of a dense table at
`p = 2`, `ne = 192`. The mass matrix that comes out of it is sparse too, and is wrapped in a
`SimpleSplines.FactorizedMass`, i.e. a sparse Cholesky, so that [`project!`](@ref) and
`mass_solve!` work here as they do on a spline space. `CirculantMass` is *not* usable: it
verifies circulance, and for ``p > 1`` the matrix is block-circulant with ``p \times p``
blocks rather than circulant, the nodes within an element not being translates of one
another. Only [`inverse_mass_matrix`](@ref) is dense, and it is built by solving against
the identity.
[`mixed_matrix`](@ref) is memoised, so the stiffness matrix is assembled once per space
rather than once per Newton iteration.

# A trap

The matrix ``\mathbb{K} = \mathbb{M}^{-1} \mathbb{D} \mathbb{M}^{-1}`` is antisymmetric, and
an antisymmetric matrix has *even* rank. For an even number of degrees of freedom its corank
is therefore two rather than one, and the extra kernel vector — the sawtooth mode
``(-1)^i`` — is a **spurious** discrete Casimir with no continuum counterpart. Choose `p`
and `ne` so that ``N = p \, n_e`` is odd.
"""
struct LagrangeSpace{T, MO <: MassOperator{T}} <: DiscreteSpace{T}
    p::Int
    ne::Int
    L::T
    nq::Int
    x::Vector{T}                        # the nodes, i.e. the coordinates of the dofs
    xq::Vector{T}                       # global quadrature nodes
    w::Vector{T}                        # global quadrature weights
    Φ::Vector{SparseMatrixCSC{T, Int}}  # Φ[d+1][i,r]
    mass::MO
    Minv::Matrix{T}
    integrals::Vector{T}
    scratch::Vector{T}
    cache::Dict{Tuple{Int, Int}, SparseMatrixCSC{T, Int}}

    function LagrangeSpace{T}(p::Integer, ne::Integer; L = 2π,
            nq::Integer = 2p + 4) where {T}
        p ≥ 1 || throw(ArgumentError(
            "the degree of a Lagrange element must be at least one, got p = $(p)"))
        ne ≥ 1 || throw(ArgumentError(
            "at least one element is needed, got ne = $(ne)"))

        N = p * ne
        Lx = convert(T, L)
        h = Lx / ne

        # the local basis on the reference element, equispaced nodes, from the sibling
        # package rather than re-derived here
        ξloc = collect(range(zero(T), one(T), length = p + 1))
        ℓ = Lagrange(ξloc)
        dℓ = ℓ'

        ξ = gauss_legendre_nodes(T, nq)
        ω = gauss_legendre_weights(T, nq)

        x = zeros(T, N)
        xq = Vector{T}(undef, ne * nq)
        w = Vector{T}(undef, ne * nq)

        # Only the p+1 basis functions supported on an element are evaluated there, and only
        # those entries are stored. A dense tabulation is (p+1)/N full -- under one per cent
        # at the resolutions these runs reach -- and makes every contraction
        # Φ diag(f w) Φᵀ cost N² per quadrature point instead of (p+1)².
        nnzΦ = ne * (p + 1) * nq
        Is = Vector{Int}(undef, nnzΦ)
        Js = Vector{Int}(undef, nnzΦ)
        Vs = [zeros(T, nnzΦ) for _ in 0:1]

        t = 0
        for e in 0:(ne - 1)
            gidx = [mod1(e * p + a + 1, N) for a in 0:p]
            # Only the first p local nodes are *owned* by this element; local node p is
            # local node 0 of the next one. Assigning it here too would be harmless inside
            # the domain but wrong at the wrap, where the last element would overwrite the
            # coordinate of node 1 with L instead of 0.
            for a in 0:(p - 1)
                x[gidx[a + 1]] = e * h + ξloc[a + 1] * h
            end
            for r in 1:nq
                q = e * nq + r
                xq[q] = e * h + h * ξ[r]
                w[q] = h * ω[r]
                for a in 0:p
                    t += 1
                    Is[t] = gidx[a + 1]
                    Js[t] = q
                    Vs[1][t] = ℓ[ξ[r], a + 1]
                    # the chain rule of the map from the reference element: d/dx = (1/h) d/dξ
                    Vs[2][t] = dℓ[ξ[r], a + 1] / h
                end
            end
        end

        # `sparse` SUMS duplicate (i,q) pairs, which is the `+=` the dense assembly used to
        # do. The only way to get a duplicate here is ne == 1, where the wrap makes local
        # nodes 0 and p the same degree of freedom on the one element.
        Φ = [sparse(Is, Js, Vs[d + 1], N, ne * nq) for d in 0:1]

        M = Φ[1] * Diagonal(w) * Φ[1]'
        M = (M + M') / 2                # symmetric by construction; enforce it exactly

        # A `FactorizedMass`, i.e. a sparse Cholesky, chosen here rather than through
        # `mass_operator(M, basis)`: that reads the representation off a `SimpleSplines`
        # basis, and a Lagrange space is not built on one. Neither alternative would fit in
        # any case. `CirculantMass` verifies circulance, which this matrix has only at
        # p == 1; for p > 1 it is block-circulant in p×p blocks. `BandedMass` takes its
        # half-bandwidth from the stored entries, and the corner entries of the periodic
        # wrap put that at N - 1, so it would factorise a dense matrix and call it banded.
        # The mass solve is not where the time goes -- the assembly is, and that is what
        # the sparsity above is for.
        mass = try
            FactorizedMass(M)
        catch err
            err isa ArgumentError && throw(ArgumentError(
                "the mass matrix assembled with nq = $(nq) points per element is not " *
                "positive definite for degree-$(p) elements; raise nq"))
            rethrow()
        end

        # ∫ φ_i dx and the f ⊙ w buffer are constants of the discretisation, assembled once
        # rather than rebuilt on every call, as in `SimpleSplines.SplineQuadrature`.
        integrals = Φ[1] * w
        scratch = Vector{T}(undef, ne * nq)

        # The dense inverse is built by solving against the identity rather than by inverting
        # an assembled matrix, matching `SplineSpace`.
        Minv = mass \ Matrix{T}(I, N, N)

        new{T, typeof(mass)}(Int(p), Int(ne), Lx, Int(nq), x, xq, w, Φ, mass, Minv,
            integrals, scratch,
            Dict{Tuple{Int, Int}, SparseMatrixCSC{T, Int}}())
    end
end

function LagrangeSpace(p::Integer, ne::Integer; L = 2π, kwargs...)
    LagrangeSpace{typeof(float(L))}(p, ne; L = L, kwargs...)
end

nbasis(s::LagrangeSpace) = s.p * s.ne
degree(s::LagrangeSpace) = s.p
order(s::LagrangeSpace) = s.p + 1
ncells(s::LagrangeSpace) = s.ne
domainlength(s::LagrangeSpace) = s.L
nodes(s::LagrangeSpace) = s.x
quadrature_nodes(s::LagrangeSpace) = s.xq
quadrature_weights(s::LagrangeSpace) = s.w
mass_matrix(s::LagrangeSpace) = mass_matrix(s.mass)
mass_operator(s::LagrangeSpace) = s.mass
inverse_mass_matrix(s::LagrangeSpace) = s.Minv
basis_integrals(s::LagrangeSpace) = s.integrals

# The `MassOperator` itself, so that the generic `project`, `project!` and `mass_solve!` of
# this file work on a Lagrange space at all: they go through `mass_solve!`, which is defined
# for a `MassOperator` and not for a bare `Cholesky`.
mass_factorization(s::LagrangeSpace) = s.mass

"""
    mixed_matrix(space::LagrangeSpace, a, b)

The matrix ``\\int_\\Omega D^a \\phi_k \\, D^b \\phi_l \\, dx``, memoised.

These are constants of the discretisation, and the stiffness matrix in particular would
otherwise be reassembled on every Newton iteration of every step. `SplineSpace` gets the
same memoisation by delegating to `SimpleSplines.SplineQuadrature`.
"""
function mixed_matrix(s::LagrangeSpace{T}, a::Integer, b::Integer) where {T}
    get!(s.cache, (Int(a), Int(b))) do
        A = basis_values(s, a) * Diagonal(s.w) * basis_values(s, b)'
        SparseMatrixCSC{T, Int}(A)
    end
end

function basis_values(s::LagrangeSpace, d::Integer = 0)
    0 ≤ d ≤ 1 || throw(ArgumentError(
        "a Lagrange finite element basis is only C⁰, so derivatives beyond the first do " *
        "not exist across element boundaries; order $(d) was requested"))
    s.Φ[d + 1]
end

function evaluate(s::LagrangeSpace, û::AbstractVector, x::Number, d::Integer = 0)
    _lagrange_evaluate(s, û, x, d)
end

function evaluate(s::LagrangeSpace, û::AbstractVector, X::AbstractVector, d::Integer = 0)
    [_lagrange_evaluate(s, û, x, d) for x in X]
end

function _lagrange_evaluate(s::LagrangeSpace{T}, û::AbstractVector, x::Number,
        d::Integer) where {T}
    length(û) == nbasis(s) || throw(DimensionMismatch(
        "the coefficient vector has $(length(û)) entries but the space has $(nbasis(s))"))
    0 ≤ d ≤ 1 || throw(ArgumentError(
        "a Lagrange finite element basis is only C⁰; order $(d) was requested"))

    N = nbasis(s)
    h = s.L / s.ne
    x̃ = mod(x, s.L)
    e = min(floor(Int, x̃ / h), s.ne - 1)        # the element containing x̃
    ξ = (x̃ - e * h) / h

    ℓ = Lagrange(collect(range(zero(T), one(T), length = s.p + 1)))
    B = d == 0 ? ℓ : ℓ'

    v = zero(T)
    for a in 0:s.p
        v += û[mod1(e * s.p + a + 1, N)] * B[ξ, a + 1]
    end
    d == 1 ? v / h : v
end

@doc raw"""
    nodal_derivative_matrix(space::LagrangeSpace)

The matrix ``E_{qm} = \phi_q'(x_m)`` of basis-function derivatives evaluated at the nodes,
with the two one-sided limits averaged.

A ``\mathcal{C}^0`` basis has a jump in ``\phi'`` across every element boundary, and the
nodes sit exactly on those jumps. Each node is therefore visited once per element that
contains it, and the value recorded is the mean of what those elements say — a one-sided
limit being zero on a side where ``\phi_q`` vanishes identically.

This is not one convention among several. It is the only choice that leaves the structure
coefficients built from `E` antisymmetric, which is what the direct discretisation of a
Lie-Poisson bracket needs before any question about the Jacobi identity can even be posed.
"""
function nodal_derivative_matrix(s::LagrangeSpace{T}) where {T}
    N = nbasis(s)
    h = s.L / s.ne
    ξloc = collect(range(zero(T), one(T), length = s.p + 1))
    dℓ = Lagrange(ξloc)'

    E = zeros(T, N, N)
    count = zeros(Int, N)
    for e in 0:(s.ne - 1), b in 0:s.p

        count[mod1(e * s.p + b + 1, N)] += 1
    end
    for e in 0:(s.ne - 1)
        gidx = [mod1(e * s.p + a + 1, N) for a in 0:s.p]
        for a in 0:s.p, b in 0:s.p

            E[gidx[a + 1], gidx[b + 1]] += dℓ[ξloc[b + 1], a + 1] / h / count[gidx[b + 1]]
        end
    end
    return E
end
