
@doc raw"""
    DiscreteHamiltonian{T}

A function of the degrees of freedom, obtained by evaluating a continuous functional on the
discrete field ``u_h = \sum_j \hat{u}_j \phi_j``.

The interface is [`hamiltonian`](@ref), [`gradient`](@ref) and [`hessian`](@ref). Casimirs
are Hamiltonians too as far as this interface is concerned — what makes a quantity a Casimir
is a relation between its gradient and a bracket, not anything about the quantity itself, so
[`MassCasimir`](@ref) is an ordinary `DiscreteHamiltonian`.

Every gradient here is analytic. That is not an optimisation: the Gonzalez discrete gradient
and the Poisson-map test both differentiate the gradient again, and a finite-difference
Jacobian has a floor several orders above the signal those tests are trying to measure.
"""
abstract type DiscreteHamiltonian{T} end

Base.eltype(::DiscreteHamiltonian{T}) where {T} = T

"""
    hamiltonian(H, space, û)

The value of `H` at the degrees of freedom `û`.
"""
function hamiltonian end

"""
    gradient(H, space, û)

The gradient ``\\partial H / \\partial \\hat{u}_i``, computed analytically.
"""
function gradient end

"""
    hessian(H, space, û)

The Hessian ``\\partial^2 H / \\partial \\hat{u}_i \\partial \\hat{u}_j``, symmetric.
"""
function hessian end

"""
    gradient!(g, H, space, û)

In-place [`gradient`](@ref). The fallback allocates; the concrete Hamiltonians that are
evaluated in a time loop override it.
"""
gradient!(g::AbstractVector, H::DiscreteHamiltonian, s::DiscreteSpace, û::AbstractVector) =
    (g .= gradient(H, s, û); g)

value(H::DiscreteHamiltonian, s::DiscreteSpace, û::AbstractVector) = hamiltonian(H, s, û)


@doc raw"""
    QuadraticHamiltonian(A)

The quadratic form ``H = \tfrac{1}{2} \hat{u}^T A \hat{u}`` with `A` symmetric.

With `A` the mass matrix this is ``\tfrac{1}{2} \int_\Omega u_h^2 \, dx``, which is the
second KdV Hamiltonian and the Burgers Hamiltonian both.

The gradient ``A \hat{u}`` is what makes the second KdV flow collapse: the mass matrix in
the gradient cancels one of the two inverse mass matrices of the bracket, and what is left
is the plain Galerkin scheme with no auxiliary variable.
"""
struct QuadraticHamiltonian{T, MT <: AbstractMatrix{T}} <: DiscreteHamiltonian{T}
    A::MT

    function QuadraticHamiltonian(A::MT) where {T, MT <: AbstractMatrix{T}}
        isapprox(A, A'; atol = 1e-12 * max(one(T), maximum(abs, A))) || throw(ArgumentError(
            "the matrix of a quadratic Hamiltonian must be symmetric"))
        new{T, MT}(A)
    end
end

hamiltonian(H::QuadraticHamiltonian, s::DiscreteSpace, û::AbstractVector) =
    dot(û, H.A, û) / 2
gradient(H::QuadraticHamiltonian, s::DiscreteSpace, û::AbstractVector) = H.A * û
hessian(H::QuadraticHamiltonian, s::DiscreteSpace, û::AbstractVector) = H.A


@doc raw"""
    MassCasimir(space)

The total mass ``C_0 = \int_\Omega u \, dx = \sum_i \hat{u}_i \int_\Omega \phi_i \, dx``.

Linear, so its gradient ``g_i = \int_\Omega \phi_i \, dx`` is constant and its Hessian
vanishes.

# Why it is conserved for free

``g`` spans the kernel of the first discrete bracket, ``\mathbb{P}^1 g = 0``. Every method
whose increment lies in the range of ``\mathbb{P}^1`` — that is, every method of the form
``\hat{u}^{n+1} = \hat{u}^n + \Delta t \, \mathbb{P}^1 \xi`` for any ``\xi`` whatever —
therefore leaves ``C_0`` untouched exactly, explicit Euler included. This is the only
conservation law in the package that costs nothing.

Under the second bracket the same conclusion holds, but for a weaker reason: there
``g^T \mathbb{P}^2 (\hat{u}) \, \partial H_2 / \partial \hat{u} = 0`` pointwise in
``\hat{u}`` rather than ``g^T \mathbb{P}^2 = 0``, so the mass is Hamiltonian-specific rather
than Casimir-strength, and only methods built from evaluations of *that* field inherit it.
"""
struct MassCasimir{T} <: DiscreteHamiltonian{T}
    g::Vector{T}
end

MassCasimir(s::DiscreteSpace{T}) where {T} = MassCasimir{T}(basis_integrals(s))

hamiltonian(H::MassCasimir, s::DiscreteSpace, û::AbstractVector) = dot(H.g, û)
gradient(H::MassCasimir, s::DiscreteSpace, û::AbstractVector) = H.g
hessian(H::MassCasimir{T}, s::DiscreteSpace, û::AbstractVector) where {T} =
    zeros(T, length(û), length(û))


@doc raw"""
    GaugedCasimir(g, G)

The Casimir ``C = \sum_i n_i \, G(\hat{u}_i)`` of a [`GaugedBracket`](@ref), with `n` in
the kernel of its constant matrix ``\mathbb{K}`` and ``G' = 1/g``.

For the Burgers bracket, ``g(u) = \sqrt{u}`` gives ``G(u) = 2\sqrt{u}`` and
``n_i = \int_\Omega \phi_i \, dx``, so
``C = 2 \sum_i \left( \int \phi_i \right) \sqrt{u_i}``, an *exact* Casimir of the discrete
bracket and a consistent quadrature of the continuous ``2 \int \sqrt{u} \, dx``.
"""
struct GaugedCasimir{T, GT, DT} <: DiscreteHamiltonian{T}
    n::Vector{T}
    G::GT
    dG::DT
end

hamiltonian(H::GaugedCasimir, s::DiscreteSpace, û::AbstractVector) =
    sum(H.n[i] * H.G(û[i]) for i in eachindex(û))
gradient(H::GaugedCasimir, s::DiscreteSpace, û::AbstractVector) = H.n .* H.dG.(û)
hessian(H::GaugedCasimir, s::DiscreteSpace, û::AbstractVector) =
    throw(ArgumentError("the Hessian of a GaugedCasimir is not implemented"))
