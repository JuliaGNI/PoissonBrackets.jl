
@doc raw"""
    kdv_bracket_1(space)

The first discrete KdV bracket,

```math
\mathbb{P}^1 = \mathbb{M}^{-1} \, \tfrac{1}{2} \left( S - S^T \right) \, \mathbb{M}^{-1} ,
\qquad S_{kl} = \int_\Omega \phi_k \, \partial_x \phi_l \, dx ,
```

the discretisation of ``\{A,B\}_1 = \int_\Omega \frac{\delta A}{\delta u} \partial_x
\frac{\delta B}{\delta u} \, dx``.

Constant and antisymmetric, hence **exactly Poisson**: the Jacobi identity for a constant
bracket follows from antisymmetry alone.

Note that ``S`` is already exactly antisymmetric on a periodic mesh, so ``S - S^T = 2S`` and
the factor of one half does *not* cancel against a doubling. Writing the assembly in the
skew-symmetrised form anyway is what makes the antisymmetry independent of the quadrature.

Its kernel contains the gradient of the mass, so [`MassCasimir`](@ref) is a genuine Casimir.
On a uniform mesh the rank is ``N-2`` for even ``N`` and ``N-1`` for odd ``N`` — an
antisymmetric matrix has even rank — the zero and Nyquist modes lying in the kernel.
"""
kdv_bracket_1(s::DiscreteSpace) = ConstantBracket(s, derivative_matrix(s))

@doc raw"""
    kdv_bracket_2(space)

The second discrete KdV bracket, the discretisation of
``\{A,B\}_2 = \int_\Omega \frac{\delta A}{\delta u} \left( -4u \partial_x - 2u_x -
\partial_x^3 \right) \frac{\delta B}{\delta u} \, dx``, in explicitly skew-symmetrised form,

```math
\mathbb{P}^2_{ij} (\hat{u}) = \mathbb{M}^{-1} \left(
      \tfrac{1}{2} \left( C - C^T \right)
    - 2 \int_\Omega u_h \left( \phi_k \partial_x \phi_l
                             - \phi_l \partial_x \phi_k \right) dx
  \right) \mathbb{M}^{-1} ,
\qquad C_{kl} = \int_\Omega \phi_k' \phi_l'' \, dx .
```

The third-derivative term has been integrated by parts once, ``-\int \phi_k \phi_l''' =
\int \phi_k' \phi_l''``, which lowers the regularity required of the basis from ``p \ge 3``
to ``p \ge 2``. At ``p = 2`` this is not a convenience but the only correct reading: the
un-integrated form vanishes identically there, every degree-2 spline having zero third
derivative inside each cell, while the integrated one does not.

# Not Poisson

This bracket is exactly antisymmetric on any mesh and at any quadrature, and it does **not**
satisfy the Jacobi identity. The normalised residual is about `0.42` and is flat under
refinement.

That is a result, not a defect. The operator ``-4u\partial_x - 2u_x - \partial_x^3 =
-2(u\partial_x + \partial_x u) - \partial_x^3`` is the Virasoro Lie-Poisson structure, so the
Jacobi identity would require the discrete coefficients to close into a Lie algebra — a
quadratic condition on the basis, not something a choice of quadrature or a regrouping of
the integrand can deliver. See [`kdv_miura_bracket`](@ref) for the construction that does
give an exactly Poisson second bracket, and what it costs.
"""
kdv_bracket_2(s::DiscreteSpace) =
    AffineBracket(s, -2, basis_values(s, 0), mixed_matrix(s, 1, 2))

@doc raw"""
    kdv_miura_bracket(space, v̂)

The second KdV bracket obtained by pushing the first forward along the Miura map
``u = -(v^2 + v_x)``, at the mKdV field `v̂`; see [`MiuraBracket`](@ref).

Antisymmetric *and* Poisson to round-off at every `N`, at the cost of being dense. It is a
bracket on ``v``-space; to run it, use [`MiuraSystem`](@ref) rather than forming this matrix.
"""
kdv_miura_bracket(s::DiscreteSpace, v̂::AbstractVector) =
    MiuraBracket(s, v̂, kdv_bracket_1(s))


@doc raw"""
    KdVHamiltonian1()

The first KdV Hamiltonian, ``H_1 = \tfrac{1}{2} \int_\Omega (u_x^2 - 2u^3) \, dx``, whose
variational derivative is ``\delta H_1/\delta u = -3u^2 - u_{xx}``.

Cubic in the degrees of freedom, with

```math
\frac{\partial H_1}{\partial \hat{u}_i}
    = \int_\Omega \left( \partial_x \phi_i \, u_{h,x} - 3 \phi_i \, u_h^2 \right) dx ,
\qquad
\frac{\partial^2 H_1}{\partial \hat{u}_i \partial \hat{u}_j}
    = \mathbb{K}^1_{ij} - 6 \int_\Omega \phi_i \phi_j \, u_h \, dx .
```

No integration by parts is needed for the gradient, so it is exact for ``p \ge 1``.

Being *cubic* is what puts it outside the reach of the theorem that a symplectic
Runge-Kutta method conserves quadratic invariants — which is why implicit midpoint on the
first flow does not conserve it, and why an energy-preserving method has to be used instead
if exact conservation is wanted.
"""
struct KdVHamiltonian1{T} <: DiscreteHamiltonian{T} end

KdVHamiltonian1() = KdVHamiltonian1{Float64}()

function hamiltonian(::KdVHamiltonian1, s::DiscreteSpace, û::AbstractVector)
    w  = quadrature_weights(s)
    uh = field(s, û, 0)
    ux = field(s, û, 1)
    (dot(w, ux .^ 2) - 2 * dot(w, uh .^ 3)) / 2
end

function gradient(::KdVHamiltonian1, s::DiscreteSpace, û::AbstractVector)
    w  = quadrature_weights(s)
    uh = field(s, û, 0)
    ux = field(s, û, 1)
    basis_values(s, 1) * (w .* ux) .- 3 .* (basis_values(s, 0) * (w .* uh .^ 2))
end

function hessian(::KdVHamiltonian1, s::DiscreteSpace, û::AbstractVector)
    uh = field(s, û, 0)
    H  = stiffness_matrix(s) .- 6 .* weighted_matrix(s, uh, 0, 0)
    (H + H') / 2
end


@doc raw"""
    KdVHamiltonian2(space)

The second KdV Hamiltonian, ``H_2 = \tfrac{1}{2} \int_\Omega u^2 \, dx``, whose variational
derivative is ``\delta H_2 / \delta u = u``.

Quadratic, with gradient ``\mathbb{M} \hat{u}``. That mass matrix is what makes the second
flow collapse: it cancels one of the two inverse mass matrices of
[`kdv_bracket_2`](@ref), leaving the plain Galerkin scheme

```math
\int_\Omega \phi_i \, \partial_t u_h \, dx
    = -\int_\Omega \left( 6 u_h u_{h,x} \phi_i - \partial_x \phi_i \, u_{h,xx} \right) dx
```

with no auxiliary variable and a single sparse solve, where the first flow is a *mixed*
two-field method with two.
"""
KdVHamiltonian2(s::DiscreteSpace) = QuadraticHamiltonian(mass_matrix(s))


@doc raw"""
    KdVHamiltonian3()

The Hamiltonian of the fifth-order member of the KdV hierarchy,

```math
H_3 = \int_\Omega \left( \tfrac{5}{2} u^4 - 5 u u_x^2 + \tfrac{1}{2} u_{xx}^2 \right) dx ,
```

whose variational derivative is ``10u^3 + 10 u u_{xx} + 5 u_x^2 + u_{xxxx}``.

Used for the Magri hierarchy tests: ``\mathbb{P}^2 \, \partial H_1/\partial \hat{u}`` should
reproduce this flow, and does so at ``O(h^{2p})`` — the first rung of the hierarchy is exact
at finite ``N``, the second only to the order of the discretisation.
"""
struct KdVHamiltonian3{T} <: DiscreteHamiltonian{T} end

KdVHamiltonian3() = KdVHamiltonian3{Float64}()

function hamiltonian(::KdVHamiltonian3, s::DiscreteSpace, û::AbstractVector)
    w   = quadrature_weights(s)
    uh  = field(s, û, 0)
    ux  = field(s, û, 1)
    uxx = field(s, û, 2)
    dot(w, (5 // 2) .* uh .^ 4 .- 5 .* uh .* ux .^ 2 .+ uxx .^ 2 ./ 2)
end

function gradient(::KdVHamiltonian3, s::DiscreteSpace, û::AbstractVector)
    w   = quadrature_weights(s)
    uh  = field(s, û, 0)
    ux  = field(s, û, 1)
    uxx = field(s, û, 2)
    10 .* (basis_values(s, 0) * (w .* uh .^ 3)) .-
     5 .* (basis_values(s, 0) * (w .* ux .^ 2)) .-
    10 .* (basis_values(s, 1) * (w .* uh .* ux)) .+
          (basis_values(s, 2) * (w .* uxx))
end


## Initial data
#
# The convention is the textbook one, u_t + 6 u u_x + u_xxx = 0, so the sech^2 solitons are
# ELEVATIONS. Getting the sign wrong gives a profile that steepens and blows up instead of
# translating, which is the quickest way to tell the two conventions apart numerically.
#
# On a torus a soliton is not an exact solution, but for the boxes used here its tails are
# below round-off, so the periodic images are summed and the result is periodic to machine
# precision.

const SOLITON_IMAGES = -2:2

# cosh overflows for arguments beyond about 710 in Float64, and 1/cosh(x)^2 is already zero
# to round-off far below that, so the argument is clipped rather than guarded downstream.
const COSH_CLIP = 350.0

@doc raw"""
    soliton(κ, x₀, L; t = 0)

The single soliton of ``u_t + 6uu_x + u_{xxx} = 0``, wrapped onto ``[0,L)``, as a function of
`x`:

```math
u(t,x) = 2 \kappa^2 \operatorname{sech}^2
         \left( \kappa \left( x - 4\kappa^2 t - x_0 \right) \right) .
```

An **elevation** of height ``2\kappa^2``, travelling right at speed ``4\kappa^2``.

The centre is reduced modulo `L` before the periodic images are summed, so the formula stays
valid at large `t`, where the soliton has crossed the box many times and an unreduced
argument would fall outside the reach of any finite image sum.
"""
function soliton(κ::Real, x₀::Real, L::Real; t::Real = 0)
    centre = mod(x₀ + 4κ^2 * t, L)
    function u(x::Number)
        v = zero(float(κ))
        for r in SOLITON_IMAGES
            ξ = clamp(κ * (x + r * L - centre), -COSH_CLIP, COSH_CLIP)
            v += 2κ^2 / cosh(ξ)^2
        end
        return v
    end
end

@doc raw"""
    two_soliton(k₁, k₂, x₁, x₂, L; t = 0)

Hirota's exact two-soliton of ``u_t + 6uu_x + u_{xxx} = 0``, wrapped onto ``[0,L)``.

From the tau function,

```math
f = 1 + e^{\eta_1} + e^{\eta_2} + A_{12} e^{\eta_1 + \eta_2} ,
\qquad \eta_i = 2 k_i \left( x - x_i - 4 k_i^2 t \right) ,
\qquad A_{12} = \left( \frac{k_1 - k_2}{k_1 + k_2} \right)^2 ,
```

with ``u = 2 (\log f)_{xx}``. The second derivative of the logarithm is taken in closed form
rather than numerically, and the exponentials are written through a shifted maximum so that
the tails do not overflow.

!!! warning "Exact at t = 0 only"
    This is the exact solution on the line. On the torus it is exact at `t = 0`, but it
    carries the phase shift of *one* overtaking, whereas on a periodic domain the two
    solitons meet again and again. Use it as initial data, not as a reference solution at
    large `t`.
"""
function two_soliton(k₁::Real, k₂::Real, x₁::Real, x₂::Real, L::Real; t::Real = 0)
    A₁₂ = ((k₁ - k₂) / (k₁ + k₂))^2
    c₁ = mod(x₁ + 4k₁^2 * t, L)
    c₂ = mod(x₂ + 4k₂^2 * t, L)
    logA = log(A₁₂)

    function u(x::Number)
        v = 0.0
        for r in SOLITON_IMAGES
            xs = x + r * L
            e₁ = 2k₁ * (xs - c₁)
            e₂ = 2k₂ * (xs - c₂)
            # f = Σ_a exp(p_a) with p_a linear in x; the shift by the maximum keeps the
            # exponentials in range without changing the ratios f'/f and f''/f
            ps = (0.0, e₁, e₂, e₁ + e₂ + logA)
            ds = (0.0, 2k₁, 2k₂, 2 * (k₁ + k₂))
            m  = maximum(ps)
            ws = ntuple(a -> exp(clamp(ps[a] - m, -700.0, 700.0)), 4)
            f  = sum(ws)
            f1 = sum(ds[a] * ws[a] for a in 1:4)
            f2 = sum(ds[a]^2 * ws[a] for a in 1:4)
            v += 2 * (f2 / f - (f1 / f)^2)
        end
        return v
    end
end

@doc raw"""
    solitons(κs, centres, L; t = 0)

A superposition of single [`soliton`](@ref) profiles — the standard multi-soliton benchmark.

A sum of ``\operatorname{sech}^2`` profiles is *not* an exact solution of the nonlinear
equation, only an approximation good to the overlap of the tails; the waves sort themselves
out into a genuine multi-soliton within the first few time units.

This is the initial condition of the three- and five-soliton problems of Shi, Fu and Liu.
That paper writes KdV as ``u_t = \alpha u u_x + \nu u_{xxx}`` with ``\alpha = \nu = -1`` and
gives elevations ``12\kappa_i^2 \operatorname{sech}^2``; our convention is reached by
``u \to u/6``, which turns those into exactly the elevations ``2\kappa^2
\operatorname{sech}^2`` above, at the same speeds ``4\kappa^2``. Their momentum and energy
are then ``36 H_{2,d}`` and ``36 H_{1,d}``, and their mass is ``6 C_{0,d}``.

The three Hamiltonians take the *same numerical values* on corresponding data as under the
older ``u_t = 6uu_x - u_{xxx}`` convention. Not because they are even — ``H_1`` is not — but
because the two conventions are related by ``u \to -u`` and each Hamiltonian carries the
matching sign: ``H_1^{\mathrm{new}}(-w) = H_1^{\mathrm{old}}(w)``, and likewise for ``H_3``.
So the elevation soliton here has exactly the ``H_1`` the depression soliton had there, which
is checked against `plot_kdv_energies.py` to eight digits. Only ``C_{0,d}`` differs, and it
differs by a sign.
"""
function solitons(κs, centres, L::Real; t::Real = 0)
    length(κs) == length(centres) || throw(DimensionMismatch(
        "got $(length(κs)) wave numbers but $(length(centres)) centres"))
    parts = [soliton(κ, c, L; t = t) for (κ, c) in zip(κs, centres)]
    x -> sum(f(x) for f in parts)
end

@doc raw"""
    three_solitons(; L = 200.0, xleft = -100.0)
    five_solitons(; L = 300.0, xleft = -150.0)

The three- and five-soliton benchmark initial conditions of Shi, Fu and Liu,
*Appl. Math. Comput.* **508** (2026) 129620, Problems 4.2 and 4.3, as functions of `x`.

| | wave numbers ``\kappa`` | centres | box | ``T`` |
|:--|:--|:--|:--|:--|
| `three_solitons` | `0.3, 0.25, 0.2` | `-60, -44, -26` | ``[-100, 100)`` | `250` |
| `five_solitons` | `0.3, 0.25, 0.2, 0.15, 0.1` | `-120, -90, -60, -30, 0` | ``[-150, 150)`` | `500` |

The solver works on ``[0, L)``, so the centres are shifted by the left endpoint here and
only the plots carry it back.

That paper writes KdV as ``u_t = \alpha u u_x + \nu u_{xxx}`` with ``\alpha = \nu = -1`` and
gives elevations ``12\kappa_i^2 \operatorname{sech}^2``; our convention is reached by
``u \to u/6``, which turns those into exactly the elevations ``2\kappa^2
\operatorname{sech}^2`` of [`soliton`](@ref), at the same speeds ``4\kappa^2``. Their
momentum and energy are then ``36 H_{2,d}`` and ``36 H_{1,d}``, and their mass is
``6 C_{0,d}``, so the invariants correspond one to one.

In both, the tallest wave overtakes all the others over the course of the run.
"""
three_solitons(; L = 200.0, xleft = -100.0) =
    solitons((0.3, 0.25, 0.2), (-60.0, -44.0, -26.0) .- xleft, L)

@doc (@doc three_solitons)
five_solitons(; L = 300.0, xleft = -150.0) =
    solitons((0.3, 0.25, 0.2, 0.15, 0.1), (-120.0, -90.0, -60.0, -30.0, 0.0) .- xleft, L)

@doc raw"""
    miura_initial_v(L = 2π)

The initial condition ``v_0 = 1 + \sin(2\pi x / L)`` of the Miura example, posed in the
**mKdV** variable.

It has to be posed there. The image of the discrete Miura map lies in ``C_{0,d} \le 0``,
since ``\int u_h = -\int v_h^2``, with equality only for the trivial field, and none of the
other examples here is in it — ``\cos x`` has zero mass without being zero, and the solitons
are elevations, hence of positive mass. Flipping the sign convention does not help: it moves
the half-space along with the solitons. The corresponding KdV field is
``u_0 = \mathcal{M}_h(v_0)``; see [`MiuraSystem`](@ref).
"""
miura_initial_v(L = 2π) = x -> 1 + sin(2π * x / L)

"""
    cosine(L)

``u_0 = \\cos(2\\pi x / L)``, the single-mode example; `cos(x)` on ``[0, 2\\pi)``.

The one initial condition here whose mass is exactly zero, which is why the mass drift is
reported in absolute rather than relative terms.
"""
cosine(L::Real) = x -> cos(2π * x / L)


@doc raw"""
    KdVSystem(space; kwargs...)

The bi-Hamiltonian KdV system on `space`: both brackets, both Hamiltonians, the mass
Casimir, and the two flows they generate.

```jldoctest
julia> sys = KdVSystem(SplineSpace(16, 3));

julia> û = project(sys.space, sin);

julia> length(vectorfield(sys.flow1, û))
16
```

The two flows are

  - `flow1` — ``\mathbb{P}^1 \, \partial H_1/\partial \hat{u}``, a **mixed** two-field
    Galerkin method, and a genuine finite-dimensional Poisson system;
  - `flow2` — ``\mathbb{P}^2 \, \partial H_2/\partial \hat{u}``, which collapses to the
    **plain** Galerkin discretisation of ``u_t + 6uu_x + u_{xxx} = 0``.

They are not the same. The difference is exactly the projection error that the first inserts
between its two derivatives, ``\int \phi_i' (\mathrm{id} - \Pi) (3u_h^2 - u_{h,xx})``, which
converges at ``O(h^{2p})`` — the same order at which either scheme reproduces the equation
in the first place. The choice between them is therefore structural rather than a matter of
accuracy: only the first is a Poisson system.
"""
struct KdVSystem{T, ST <: DiscreteSpace{T}, B1, B2, H1, H2, F1, F2}
    space::ST
    bracket1::B1
    bracket2::B2
    H1::H1
    H2::H2
    C0::MassCasimir{T}
    flow1::F1
    flow2::F2
end

function KdVSystem(s::DiscreteSpace{T}) where {T}
    b1 = kdv_bracket_1(s)
    b2 = kdv_bracket_2(s)
    h1 = KdVHamiltonian1()
    h2 = KdVHamiltonian2(s)
    f1 = HamiltonianFlow(s, b1, h1)
    f2 = HamiltonianFlow(s, b2, h2)
    KdVSystem{T, typeof(s), typeof(b1), typeof(b2), typeof(h1), typeof(h2),
              typeof(f1), typeof(f2)}(s, b1, b2, h1, h2, MassCasimir(s), f1, f2)
end

KdVSystem(n::Integer, p::Integer; L = 2π, kwargs...) =
    KdVSystem(SplineSpace(n, p; L = L, kwargs...))

Base.eltype(::KdVSystem{T}) where {T} = T
nbasis(sys::KdVSystem) = nbasis(sys.space)


## The Miura chart

@doc raw"""
    MiuraHamiltonian(space)

The second KdV Hamiltonian pulled back along the Miura map,
``\tilde{H}(\hat{v}) = H_{2,d}(\mathcal{M}_h(\hat{v}))
 = \tfrac{1}{2} \| \Pi ( v_h^2 + v_{h,x} ) \|^2``.

Its gradient and Hessian follow from the chain rule through
``\mathbb{L} = D\mathcal{M}_h``:

```math
\frac{\partial \tilde{H}}{\partial \hat{v}} = \mathbb{L}^T \mathbb{M} \hat{u}
    = \big( 2\mathbb{B} - S \big) \hat{u} ,
\qquad
\frac{\partial^2 \tilde{H}}{\partial \hat{v}^2}
    = \mathbb{L}^T \mathbb{M} \mathbb{L} + 2 \int_\Omega u_h \phi_m \phi_n \, dx ,
```

the second term of the Hessian being the ``\hat{v}``-dependence of ``\mathbb{L}``, which is
affine, contracted against ``\mathbb{M}\hat{u}``.

# Quartic, where `H₂` is quadratic

This is the whole point of the chart. ``H_{2,d} = \tfrac{1}{2}\hat{u}^T \mathbb{M} \hat{u}``
is *quadratic* in ``\hat{u}``, so a symplectic Runge-Kutta method conserves it exactly in
that chart. Pulled back, ``\tilde{H}`` is *quartic* in ``\hat{v}`` and the theorem no longer
applies — the same equation, the same method, ten orders of magnitude apart:

| implicit midpoint on the Miura flow | ``\max\lvert\Delta H_{2,d}\rvert / \lvert H_{2,d}\rvert`` |
|:--|:--|
| in the ``\hat{v}`` chart | ``\sim 10^{-5}`` |
| in the ``\hat{u}`` chart | ``\sim 10^{-15}`` |

What is chart-free is the *value* of a function, not its algebraic form. So a discrete
gradient method — [`Gonzalez`](@ref) — does conserve ``\tilde{H}``, and hence ``H_{2,d}``,
exactly in the ``\hat{v}`` chart, because the discrete-gradient property
``\bar{g}\cdot\Delta\hat{v} = \Delta\tilde{H}`` holds by construction whatever the degree.
"""
struct MiuraHamiltonian{T} <: DiscreteHamiltonian{T}
    λ::T
end

MiuraHamiltonian() = MiuraHamiltonian{Float64}(0.0)
MiuraHamiltonian(λ::Real) = MiuraHamiltonian{typeof(float(λ))}(float(λ))
MiuraHamiltonian(::DiscreteSpace{T}; λ::Real = 0) where {T} =
    MiuraHamiltonian{T}(convert(T, λ))

hamiltonian(H::MiuraHamiltonian, s::DiscreteSpace, v̂::AbstractVector) =
    (û = miura_map(s, v̂, H.λ); dot(û, mass_matrix(s), û) / 2)

# ∂H̃/∂v̂ = Lᵀ M û, and BOTH L and û changed sign with the convention, so the product is
# numerically what it always was -- but it has to be written with the minus, because
# `miura_map` and `miura_derivative` now carry theirs.
gradient(H::MiuraHamiltonian, s::DiscreteSpace, v̂::AbstractVector) =
    .-(2 .* miura_moment_matrix(s, v̂) .- derivative_matrix(s)) * miura_map(s, v̂, H.λ)

function hessian(H::MiuraHamiltonian, s::DiscreteSpace, v̂::AbstractVector)
    L  = miura_derivative(s, v̂)
    û  = miura_map(s, v̂, H.λ)
    # LᵀML is even in L, so it is unchanged; the second term carries the sign of
    # ∂²û/∂v̂∂v̂, which is -2 M⁻¹T for û = -(v² + v_x)
    H  = L' * mass_matrix(s) * L .- 2 .* weighted_matrix(s, field(s, û), 0, 0)
    (H + H') / 2
end


@doc raw"""
    MiuraSystem(space; v̂ = nothing)
    MiuraSystem(n, p; L = 2π)

The KdV system carried in the **mKdV degrees of freedom** ``\hat{v}``, where the Poisson
matrix is the constant ``\mathbb{P}^1``:

```math
\dot{\hat{v}} = \mathbb{P}^1 \, \frac{\partial \tilde{H}}{\partial \hat{v}} ,
\qquad \tilde{H}(\hat{v}) = H_{2,d}\big( \mathcal{M}_h(\hat{v}) \big) ,
```

the semi-discrete mKdV equation, whose pushforward along ``\mathcal{M}_h`` is the Miura
bracket flow ``\dot{\hat{u}} = \mathbb{P}^2_{\mathrm{M}} \, \partial H_{2,d}/\partial\hat{u}``.

**The state is `v̂`, not `û`.** [`invariants`](@ref) and evaluation map through
[`miura_map`](@ref) first, so a `MiuraSystem` can be handed to [`integrate`](@ref) and to the
plotting routines unchanged and everything is reported in ``u``.

```jldoctest
julia> sys = MiuraSystem(SplineSpace(32, 3));

julia> v̂ = project(sys.space, x -> 1 + sin(x));

julia> inv = invariants(sys, v̂); inv.C0 < 0        # nonzero v̂, so strictly negative
true
```

# Initial data must be posed in v

They have to be: ``\int u_h = -\int v_h^2``, so the image of ``\mathcal{M}_h`` lies in
``C_{0,d} \le 0``, with equality only at ``v_h \equiv 0``, and none of the standard KdV
examples — ``\cos x`` has zero mass without being zero, the solitons are elevations of
positive mass — has a preimage at all. Use [`miura_invert`](@ref) if a particular
``\hat{u}`` is wanted, and [`hill_lambda0`](@ref) to tell in advance whether it has one.

# The two charts are not the same numerical method

The two systems are ``\mathcal{M}_h``-related exactly, so their *flows* are conjugate. For a
numerical method the corresponding statement is equivariance, and Runge-Kutta methods —
indeed all B-series methods — are equivariant under **affine** changes of variables and no
others. ``\mathcal{M}_h`` is quadratic, so equivariance fails and the two computations are
different methods of the same order, differing at ``O(\Delta t^3)`` per step.

What survives a change of chart is what can be stated without one. The Poisson-map property
is such a statement and it transports: a midpoint step in ``\hat{v}`` preserves
``\mathbb{P}^1``, and its pushforward preserves ``\mathbb{P}^2_{\mathrm{M}}``. This is the
first genuine Poisson integrator for a discrete *second* KdV structure — the Galerkin
[`kdv_bracket_2`](@ref) admits none. What does not survive is the algebraic *form* of an
invariant; see [`MiuraHamiltonian`](@ref).
"""
struct MiuraSystem{T, ST <: DiscreteSpace{T}, BT, HT, FT}
    space::ST
    bracket::BT
    H::HT
    C0::MassCasimir{T}
    flow::FT
end

function MiuraSystem(s::DiscreteSpace{T}; λ::Real = 0) where {T}
    b = kdv_bracket_1(s)
    h = MiuraHamiltonian{T}(convert(T, λ))
    f = HamiltonianFlow(s, b, h)
    MiuraSystem{T, typeof(s), typeof(b), typeof(h), typeof(f)}(s, b, h, MassCasimir(s), f)
end

"""
    miura_lambda(sys::MiuraSystem)

The spectral parameter the system was built with.
"""
miura_lambda(sys::MiuraSystem) = sys.H.λ

MiuraSystem(n::Integer, p::Integer; L = 2π, kwargs...) =
    MiuraSystem(SplineSpace(n, p; L = L, kwargs...))

Base.eltype(::MiuraSystem{T}) where {T} = T
nbasis(sys::MiuraSystem) = nbasis(sys.space)

"""
    miura_map(sys::MiuraSystem, v̂)

The KdV field ``\\hat{u} = \\mathcal{M}_h(\\hat{v})`` of a Miura system state.
"""
miura_map(sys::MiuraSystem, v̂::AbstractVector) = miura_map(sys.space, v̂, sys.H.λ)

"""
    miura_bracket(sys::MiuraSystem, v̂)

The pushforward bracket ``\\mathbb{P}^2_{\\mathrm{M}}(\\hat{v})`` at the current state.
"""
miura_bracket(sys::MiuraSystem, v̂::AbstractVector) =
    MiuraBracket(sys.space, v̂, sys.bracket)

@doc raw"""
    miura_casimir(sys::MiuraSystem, v̂)

``\int_\Omega v_h \, dx``, the Casimir of ``\mathbb{P}^1`` **in this chart**.

It is not the mass of ``u``. The quantity that maps over is
``C_{0,d} = \int_\Omega v_h^2 \, dx``, and that one is *not* a Casimir — it is conserved
because it equals ``2 \tilde{H}``'s companion under the hierarchy, not by any kernel
property.
"""
miura_casimir(sys::MiuraSystem, v̂::AbstractVector) = dot(basis_integrals(sys.space), v̂)
