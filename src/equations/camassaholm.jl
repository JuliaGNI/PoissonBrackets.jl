
# Camassa-Holm.
#
# PROTOTYPE. Unlike the Burgers and KdV discretisations, this one does not accompany a
# manuscript: it is here because the same construction applies, and to have somewhere for
# the equation to live when it does. The structural properties below are asserted by the
# test suite in the same way as for KdV, but no claim beyond those tests is made for it.
#
# The evolved variable is the MOMENTUM m = u - u_xx, not u. This is not a presentational
# choice: the Camassa-Holm bi-Hamiltonian structure is written with variational derivatives
# with respect to m, and taking them with respect to u instead gives two flows that do not
# agree with each other and neither of which is the Camassa-Holm equation.

@doc raw"""
    HelmholtzMap(space)

The discrete Helmholtz operator relating the Camassa-Holm velocity to its momentum.

In weak form ``m = u - u_{xx}`` reads
``\int \phi_i m_h = \int (\phi_i u_h + \phi_i' u_{h,x})``, i.e.

```math
\mathbb{A} \hat{u} = \mathbb{M} \hat{m} ,
\qquad \mathbb{A} = \mathbb{M} + \mathbb{K}^1 ,
```

with ``\mathbb{A}`` symmetric positive definite — it is the ``H^1`` Gram matrix — so the map
``\hat{m} \mapsto \hat{u}`` is well defined and its factorization can be formed once.

The one integration by parts this uses is exact for ``p \ge 1``, so the velocity is
recovered without ever differentiating the basis twice.
"""
struct HelmholtzMap{T}
    A::Matrix{T}
    Afact::Cholesky{T, Matrix{T}}
    M::Matrix{T}
end

function HelmholtzMap(s::DiscreteSpace{T}) where {T}
    # Kept dense here. The Helmholtz solve appears inside the gradient AND the Hessian of
    # both Camassa-Holm Hamiltonians, and the Hessian needs A⁻¹ against a MATRIX right-hand
    # side, which a sparse Cholesky will not do.
    M = Matrix(mass_matrix(s))
    A = Matrix(M .+ stiffness_matrix(s))
    A = (A + A') / 2
    HelmholtzMap{T}(A, cholesky(A), M)
end

"""
    velocity(H::HelmholtzMap, m̂)

The velocity coefficients ``\\hat{u} = \\mathbb{A}^{-1} \\mathbb{M} \\hat{m}``.
"""
velocity(H::HelmholtzMap, m̂::AbstractVector) = H.Afact \ (H.M * m̂)

"""
    momentum(H::HelmholtzMap, û)

The momentum coefficients ``\\hat{m} = \\mathbb{M}^{-1} \\mathbb{A} \\hat{u}``, the inverse
of [`velocity`](@ref).
"""
momentum(H::HelmholtzMap, û::AbstractVector) = H.M \ (H.A * û)

"""
    pullback(H::HelmholtzMap, g)

Convert a gradient with respect to ``\\hat{u}`` into one with respect to ``\\hat{m}``.

By the chain rule through ``\\hat{u} = \\mathbb{A}^{-1}\\mathbb{M}\\hat{m}`` this is
``\\mathbb{M} \\mathbb{A}^{-1} g``, the transpose of that map, both factors being symmetric.
"""
pullback(H::HelmholtzMap, g::AbstractVector) = H.M * (H.Afact \ g)


@doc raw"""
    camassa_holm_bracket_1(space)

The constant Camassa-Holm bracket, the discretisation of ``B_1 = \partial_x - \partial_x^3``
acting on the momentum,

```math
\mathbb{P}^1 = \mathbb{M}^{-1} \left(
      \tfrac{1}{2} \left( S - S^T \right)
    + \tfrac{1}{2} \left( C - C^T \right)
  \right) \mathbb{M}^{-1} ,
\qquad
C_{kl} = \int_\Omega \phi_k' \phi_l'' \, dx ,
```

with the third-derivative term integrated by parts once, exactly as for KdV, so that
``p \ge 2`` suffices.

Constant and antisymmetric, hence exactly Poisson. Paired with
[`CamassaHolmHamiltonian2`](@ref).

!!! note "Prototype"
    Camassa-Holm is not covered by either manuscript; this discretisation follows the same
    construction and is provided for experiment.
"""
function camassa_holm_bracket_1(s::DiscreteSpace)
    degree(s) ≥ 2 || throw(ArgumentError(
        "the first Camassa-Holm bracket needs the second derivative of the basis after " *
        "one integration by parts, so degree at least two; got p = $(degree(s))"))
    ConstantBracket(s, derivative_matrix(s) .+ mixed_matrix(s, 1, 2))
end

@doc raw"""
    camassa_holm_bracket_2(space)

The momentum-dependent Camassa-Holm bracket, the discretisation of
``B_2 = m \partial_x + \partial_x m`` in explicitly skew-symmetrised form,

```math
\mathbb{P}^2 (\hat{m}) = \mathbb{M}^{-1} \left(
    \int_\Omega m_h
    \left( \phi_k \partial_x \phi_l - \phi_l \partial_x \phi_k \right) dx
  \right) \mathbb{M}^{-1} .
```

This is the same [`AffineBracket`](@ref) as the second KdV bracket, with `scale` one instead
of two and no constant block; the density is the evolved field itself, so
``\Psi = \Phi_0`` as it is there. Paired with [`CamassaHolmHamiltonian1`](@ref).

Antisymmetric by construction on any mesh and at any quadrature. Like its KdV counterpart it
does **not** satisfy the Jacobi identity, and for the same reason: ``m\partial_x +
\partial_x m`` is again a Lie-Poisson structure of Virasoro type. The residual is of order
one and flat under refinement.

!!! note "Prototype"
    See [`camassa_holm_bracket_1`](@ref).
"""
camassa_holm_bracket_2(s::DiscreteSpace) =
    AffineBracket(s, 1, basis_values(s, 0), zeros(eltype(s), nbasis(s), nbasis(s)))


@doc raw"""
    CamassaHolmHamiltonian1(space)

``H_1 = \tfrac{1}{2} \int_\Omega (u^2 + u_x^2) \, dx = \tfrac{1}{2} \int_\Omega u m \, dx``,
as a function of the **momentum** degrees of freedom.

With ``\hat{u} = \mathbb{A}^{-1}\mathbb{M}\hat{m}`` this is
``\tfrac{1}{2} \hat{u}^T \mathbb{A} \hat{u}``, whose gradient with respect to ``\hat{m}`` is
simply ``\mathbb{M} \hat{u}`` — the discrete form of ``\delta H_1 / \delta m = u``.

That gradient is what makes the second flow reproduce Camassa-Holm: contracted with
[`camassa_holm_bracket_2`](@ref) it gives ``m_t = (m\partial_x + \partial_x m) u``.

!!! note "Prototype"
    See [`camassa_holm_bracket_1`](@ref).
"""
struct CamassaHolmHamiltonian1{T} <: DiscreteHamiltonian{T}
    helmholtz::HelmholtzMap{T}
end

CamassaHolmHamiltonian1(s::DiscreteSpace{T}) where {T} =
    CamassaHolmHamiltonian1{T}(HelmholtzMap(s))

function hamiltonian(H::CamassaHolmHamiltonian1, s::DiscreteSpace, m̂::AbstractVector)
    û = velocity(H.helmholtz, m̂)
    dot(û, H.helmholtz.A, û) / 2
end

gradient(H::CamassaHolmHamiltonian1, s::DiscreteSpace, m̂::AbstractVector) =
    H.helmholtz.M * velocity(H.helmholtz, m̂)

hessian(H::CamassaHolmHamiltonian1, s::DiscreteSpace, m̂::AbstractVector) =
    (B = H.helmholtz.M * (H.helmholtz.Afact \ H.helmholtz.M); (B + B') / 2)


@doc raw"""
    CamassaHolmHamiltonian2(space)

``H_2 = \tfrac{1}{2} \int_\Omega (u^3 + u u_x^2) \, dx``, as a function of the momentum.

Cubic in the velocity, with

```math
\frac{\partial H_2}{\partial \hat{u}_i} = \int_\Omega \left(
      \tfrac{3}{2} \phi_i u_h^2
    + \tfrac{1}{2} \phi_i u_{h,x}^2
    + \phi_i' \, u_h u_{h,x} \right) dx ,
```

pulled back to the momentum by ``\partial H/\partial\hat{m} =
\mathbb{M}\mathbb{A}^{-1} \, \partial H/\partial\hat{u}``.

Contracted with [`camassa_holm_bracket_1`](@ref) this gives the same Camassa-Holm equation
that the other pairing does — which is the bi-Hamiltonian property, and what the test suite
checks.

!!! note "Prototype"
    See [`camassa_holm_bracket_1`](@ref).
"""
struct CamassaHolmHamiltonian2{T} <: DiscreteHamiltonian{T}
    helmholtz::HelmholtzMap{T}
end

CamassaHolmHamiltonian2(s::DiscreteSpace{T}) where {T} =
    CamassaHolmHamiltonian2{T}(HelmholtzMap(s))

function hamiltonian(H::CamassaHolmHamiltonian2, s::DiscreteSpace, m̂::AbstractVector)
    û  = velocity(H.helmholtz, m̂)
    w  = quadrature_weights(s)
    uh = field(s, û, 0)
    ux = field(s, û, 1)
    dot(w, uh .^ 3 .+ uh .* ux .^ 2) / 2
end

function gradient(H::CamassaHolmHamiltonian2, s::DiscreteSpace, m̂::AbstractVector)
    û  = velocity(H.helmholtz, m̂)
    w  = quadrature_weights(s)
    uh = field(s, û, 0)
    ux = field(s, û, 1)
    gu = (3 // 2) .* (basis_values(s, 0) * (w .* uh .^ 2)) .+
         (1 // 2) .* (basis_values(s, 0) * (w .* ux .^ 2)) .+
                     (basis_values(s, 1) * (w .* uh .* ux))
    pullback(H.helmholtz, gu)
end

function hessian(H::CamassaHolmHamiltonian2, s::DiscreteSpace, m̂::AbstractVector)
    û  = velocity(H.helmholtz, m̂)
    uh = field(s, û, 0)
    ux = field(s, û, 1)
    Hu = 3 .* weighted_matrix(s, uh, 0, 0) .+
         weighted_matrix(s, ux, 0, 1) .+ weighted_matrix(s, ux, 1, 0) .+
         weighted_matrix(s, uh, 1, 1)
    Hu = (Hu + Hu') / 2
    B  = H.helmholtz.M * (H.helmholtz.Afact \ Hu) * (H.helmholtz.Afact \ H.helmholtz.M)
    (B + B') / 2
end


@doc raw"""
    peakon(c, x₀, L)

The peakon ``u = c \exp(-|x - x_0|)`` of the Camassa-Holm equation, wrapped onto ``[0,L)``,
as a function of `x`.

A travelling wave of speed `c`, and a weak solution only: it has a corner at its crest, so
``u_x`` jumps there and the momentum ``m = u - u_{xx}`` is a delta. A spline discretisation
cannot represent that, and the momentum it computes is a smoothed version of the delta —
which is worth knowing before reading anything into a convergence study run on this initial
condition.

!!! note "Prototype"
    See [`camassa_holm_bracket_1`](@ref).
"""
function peakon(c::Real, x₀::Real, L::Real)
    centre = mod(x₀, L)
    function u(x::Number)
        v = zero(float(c))
        for r in SOLITON_IMAGES
            v += c * exp(-abs(x + r * L - centre))
        end
        return v
    end
end


@doc raw"""
    CamassaHolmSystem(space)
    CamassaHolmSystem(n, p; L = 2π)

The Camassa-Holm system in the **momentum** variable ``\hat{m}``: both brackets, both
Hamiltonians, the mass Casimir, and the two flows.

```jldoctest
julia> sys = CamassaHolmSystem(SplineSpace(48, 3));

julia> m̂ = momentum(sys, project(sys.space, x -> 1 + 0.5sin(x)));

julia> maximum(abs, vectorfield(sys.flow1, m̂) - vectorfield(sys.flow2, m̂)) < 1e-7
true
```

The bi-Hamiltonian pairing is crossed, as it is for KdV:

  - `flow1` = ``\mathbb{P}^1`` with ``H_2``, using the constant bracket;
  - `flow2` = ``\mathbb{P}^2(\hat{m})`` with ``H_1``, using the momentum-dependent one.

Both reproduce ``m_t = (m \partial_x + \partial_x m) u``, and agree with each other to the
order of the discretisation — which is the bi-Hamiltonian property and the main thing the
tests check.

!!! note "The mass is not a Casimir of the second bracket"
    ``\int_\Omega m \, dx`` is Casimir-strength under ``\mathbb{P}^1``, where
    ``\mathbb{P}^1 g = 0`` identically, and is therefore exact there for any integrator.
    Under ``\mathbb{P}^2`` it is conserved only to the order of the discretisation: its
    continuous conservation rests on ``\int u_x m \, dx = 0``, which holds because both
    terms of ``u_x(u - u_{xx})`` are total derivatives — a statement about what the
    discretisation resolves, not an identity in ``\hat{m}``. KdV, where the mass survives
    under both brackets, is the more fortunate case rather than the general one.

Use [`velocity`](@ref) and [`momentum`](@ref) to move between ``\hat{u}`` and ``\hat{m}``.

!!! note "Prototype"
    See [`camassa_holm_bracket_1`](@ref).
"""
struct CamassaHolmSystem{T, ST <: DiscreteSpace{T}, B1, B2, H1, H2, F1, F2}
    space::ST
    helmholtz::HelmholtzMap{T}
    bracket1::B1
    bracket2::B2
    H1::H1
    H2::H2
    C0::MassCasimir{T}
    flow1::F1
    flow2::F2
end

function CamassaHolmSystem(s::DiscreteSpace{T}) where {T}
    hm = HelmholtzMap(s)
    b1 = camassa_holm_bracket_1(s)
    b2 = camassa_holm_bracket_2(s)
    h1 = CamassaHolmHamiltonian1{T}(hm)
    h2 = CamassaHolmHamiltonian2{T}(hm)
    f1 = HamiltonianFlow(s, b1, h2)      # constant bracket with the cubic Hamiltonian
    f2 = HamiltonianFlow(s, b2, h1)      # momentum-dependent bracket with the quadratic one
    CamassaHolmSystem{T, typeof(s), typeof(b1), typeof(b2), typeof(h1), typeof(h2),
                      typeof(f1), typeof(f2)}(s, hm, b1, b2, h1, h2, MassCasimir(s), f1, f2)
end

CamassaHolmSystem(n::Integer, p::Integer; L = 2π, kwargs...) =
    CamassaHolmSystem(SplineSpace(n, p; L = L, kwargs...))

"""
    velocity(sys::CamassaHolmSystem, m̂)
    momentum(sys::CamassaHolmSystem, û)

Move between the velocity and momentum coefficients of a Camassa-Holm system.
"""
velocity(sys::CamassaHolmSystem, m̂::AbstractVector) = velocity(sys.helmholtz, m̂)
momentum(sys::CamassaHolmSystem, û::AbstractVector) = momentum(sys.helmholtz, û)

Base.eltype(::CamassaHolmSystem{T}) where {T} = T
nbasis(sys::CamassaHolmSystem) = nbasis(sys.space)
