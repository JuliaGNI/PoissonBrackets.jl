#!/usr/bin/env julia
#
# III.1 -- The bracket  {A,B}_c = <c(f) {A_f, B_f}>  on the torus T^{2n} with the canonical
# bracket {a,b} = sum_k (a_{x_k} b_{v_k} - a_{v_k} b_{x_k}) satisfies the Jacobi identity for
# EVERY smooth structure function c and EVERY n -- not only for n = 1 (Thm 4.5 of the
# four-bracket paper).
#
#     julia --project=scripts scripts/fable/III1_structure_function_4d.jl
#
# The paper's proof reduces the Jacobiator to the obstruction
#
#     O(A,B,C) = < (c'(f))^2 ( {A_f,B_f}{C_f,f} + {B_f,C_f}{A_f,f} + {C_f,A_f}{B_f,f} ) >
#
# and kills the integrand pointwise by the planar Pluecker relation.  In 2n >= 4 dimensions
# the integrand does not vanish, but the cyclic sum is the totally antisymmetric tensor
# T^{ijkl} = Pi^{ij}Pi^{kl} - Pi^{ik}Pi^{jl} + Pi^{il}Pi^{jk} = (1/2)(Pi wedge Pi)^{ijkl}
# contracted with dA_f, dB_f, dC_f, and (c')^2 df = dG(f), G' = (c')^2, so that
#
#     integrand = T^{ijkl} d_i A_f d_j B_f d_k C_f d_l G(f) = d_i ( A_f T^{ijkl} d_j B_f d_k C_f d_l G(f) )
#
# is a divergence and integrates to zero on a closed domain.  Coordinate-free:
# integrand * dmu = dA_f ^ dB_f ^ dC_f ^ dG(f) ^ omega^{n-2}/(n-2)!, an exact top form.
#
#   1. exact trigonometric polynomials over Q(i): the obstruction integrand is identically zero
#      in 2-D (Pluecker), NOT identically zero in 4-D and 6-D, but its mean is exactly zero;
#      the divergence identity holds term by term; controls that must fail: a weight that is
#      not a function of f, and a single non-cyclic term;
#   2. the top-form identity S(a,b,c,d) mu = a^b^c^d^omega^{n-2}/(n-2)! for n = 2, 3;
#   3. the FULL Jacobiator, exactly, in 4-D, for c = f (control), c = f^2, c = f^3 - f, and
#      c = (f + phi(z))^2 (the class Phi(f + phi(z)) -- the only z-dependence Jacobi allows in
#      2n >= 4), against the controls c = f^2 + phi(z), c = phi(z) f and c = f d_{x1} f that
#      must fail; the first two of those pass in 2-D, where the Pluecker relation is pointwise;
#   4. the full Jacobiator in floating point on a 4-D spectral grid for c = f^{3/2} (not a
#      polynomial), refined N = 12..32, with c = f as the round-off control and c = phi(z) f
#      as the negative whose residual must reproduce the exact value of 3 at every N;
#   5. the finite-dimensional counterpart: J_ij = sum_m c_ij^m w(z_m) with w = z^2 fails for
#      se(3) (as Prop. 32-forced says) and for the su(3) sine algebra -- the "any c" property
#      is one of the continuum Poisson algebra, not of its truncations.
#
# All fields have the mean <.> = (2pi)^{-2n} integral, so no pi appears in the exact part.

using PoissonBrackets
using LinearAlgebra
using Random
using Printf

include(joinpath(@__DIR__, "..", "check.jl"));
using .Checks: header, check, check_exact, check_refined, summary, fmt

# ---------------------------------------------------------------------------
# exact trigonometric polynomials over Q(i)
# ---------------------------------------------------------------------------
# A field on T^D is a Dict wavevector => coefficient of exp(i k.z).  Rational{Int128} throws
# on overflow, so a wrong answer cannot be silent; the coefficients here stay far below that.

const Q = Rational{Int128}
const CQ = Complex{Q}

struct TP{D}
    c::Dict{NTuple{D, Int}, CQ}
end
TP{D}() where {D} = TP{D}(Dict{NTuple{D, Int}, CQ}())

clean!(p::TP) = (filter!(kv -> !iszero(kv.second), p.c); p)

function Base.:+(a::TP{D}, b::TP{D}) where {D}
    c = copy(a.c)
    for (k, v) in b.c
        c[k] = get(c, k, zero(CQ)) + v
    end
    clean!(TP{D}(c))
end
Base.:-(a::TP{D}) where {D} = TP{D}(Dict(k => -v for (k, v) in a.c))
Base.:-(a::TP, b::TP) = a + (-b)
function Base.:*(s::Number, a::TP{D}) where {D}
    clean!(TP{D}(Dict(k => CQ(s) * v for (k, v) in a.c)))
end
Base.:*(a::TP, s::Number) = s * a
Base.:+(s::Number, a::TP{D}) where {D} = a + constant(D, s)
Base.:+(a::TP, s::Number) = s + a
Base.:-(s::Number, a::TP) = s + (-a)
Base.:-(a::TP, s::Number) = a + (-s)

function Base.:*(a::TP{D}, b::TP{D}) where {D}
    c = Dict{NTuple{D, Int}, CQ}()
    for (ka, va) in a.c, (kb, vb) in b.c

        k = ka .+ kb
        c[k] = get(c, k, zero(CQ)) + va * vb
    end
    clean!(TP{D}(c))
end

"d/dz_j of a trigonometric polynomial: multiply the k-th coefficient by i k_j."
∂(a::TP{D}, j::Int) where {D} = clean!(TP{D}(Dict(k => (im * k[j]) * v for (k, v) in a.c)))

constant(D, s) = TP{D}(Dict(ntuple(_ -> 0, D) => CQ(s)))
"amp * cos(k.z)"
function cosk(k::NTuple{D, Int}, amp = 1) where {D}
    clean!(TP{D}(Dict(k => CQ(Q(amp) // 2), (.-k) => CQ(Q(amp) // 2))))
end
"amp * sin(k.z)"
function sink(k::NTuple{D, Int}, amp = 1) where {D}
    clean!(TP{D}(Dict(k => CQ(0, -Q(amp) // 2), (.-k) => CQ(0, Q(amp) // 2))))
end

"The mean <p> over the torus: the zero-mode coefficient.  Exact."
function mean(p::TP{D}) where {D}
    v = get(p.c, ntuple(_ -> 0, D), zero(CQ))
    iszero(imag(v)) || error("mean of a real field has an imaginary part: $v")
    real(v)
end
"The mean of p^2 (Parseval), exact; zero iff p is identically zero."
meansq(p::TP) = sum(abs2, values(p.c); init = zero(Q))
nterms(p::TP) = length(p.c)

"Horner evaluation of the polynomial with coefficients `coef` (coef[1] + coef[2] f + ...) at a field."
function polyval(coef, f)
    r = coef[end] + 0 * f
    for i in (length(coef) - 1):-1:1
        r = coef[i] + f * r
    end
    r
end
"Antiderivative coefficients, constant of integration zero."
antiderivative(coef) = vcat(zero(eltype(coef)), [coef[i] // i for i in eachindex(coef)])
polyderivative(coef) = [i * coef[i + 1] for i in 1:(length(coef) - 1)]

# ---------------------------------------------------------------------------
# first-order jets in epsilon, for the Gateaux derivative of a functional
# ---------------------------------------------------------------------------
# F = v + eps d, truncated at first order.  Everything below is generic in the field type, so
# the same code computes the exact Jacobiator on TP and the floating-point one on grid fields.

struct Jet{T}
    v::T
    d::T
end
Base.:+(a::Jet, b::Jet) = Jet(a.v + b.v, a.d + b.d)
Base.:-(a::Jet, b::Jet) = Jet(a.v - b.v, a.d - b.d)
Base.:-(a::Jet) = Jet(-a.v, -a.d)
Base.:*(a::Jet, b::Jet) = Jet(a.v * b.v, a.v * b.d + a.d * b.v)
Base.:*(s::Number, a::Jet) = Jet(s * a.v, s * a.d)
Base.:*(a::Jet, s::Number) = s * a
Base.:+(s::Number, a::Jet) = Jet(s + a.v, a.d)
Base.:+(a::Jet, s::Number) = s + a
∂(a::Jet, j::Int) = Jet(∂(a.v, j), ∂(a.d, j))
mean(a::Jet) = Jet(mean(a.v), mean(a.d))

# ---------------------------------------------------------------------------
# the canonical bracket on T^{2n}, coordinates (x_1..x_n, v_1..v_n)
# ---------------------------------------------------------------------------

phasedim(::TP{D}) where {D} = D
phasedim(a::Jet) = phasedim(a.v)

"{a,b} = sum_k d_{x_k} a d_{v_k} b - d_{v_k} a d_{x_k} b, generic in the field type."
function pb(a, b)
    n = phasedim(a) ÷ 2
    r = ∂(a, 1) * ∂(b, n + 1) - ∂(a, n + 1) * ∂(b, 1)
    for k in 2:n
        r = r + ∂(a, k) * ∂(b, n + k) - ∂(a, n + k) * ∂(b, k)
    end
    r
end

"The constant Poisson bivector Pi^{ij} of the canonical bracket in 2n dimensions."
function bivector(n)
    Π = zeros(Int, 2n, 2n)
    for k in 1:n
        Π[k, n + k] = 1
        Π[n + k, k] = -1
    end
    Π
end

"The nonzero entries of T^{ijkl} = Pi^{ij}Pi^{kl} - Pi^{ik}Pi^{jl} + Pi^{il}Pi^{jk}, as (i,j,k,l,value)."
function pfaffian_tensor(n)
    Π = bivector(n)
    D = 2n
    out = NTuple{5, Int}[]
    for i in 1:D, j in 1:D, k in 1:D, l in 1:D
        t = Π[i, j] * Π[k, l] - Π[i, k] * Π[j, l] + Π[i, l] * Π[j, k]
        iszero(t) || push!(out, (i, j, k, l, t))
    end
    out
end

# ---------------------------------------------------------------------------
# floating-point fields on a 4-D spectral grid
# ---------------------------------------------------------------------------

struct GF
    D::Matrix{Float64}          # one-dimensional spectral differentiation matrix
    a::Array{Float64, 4}
end
Base.:+(a::GF, b::GF) = GF(a.D, a.a .+ b.a)
Base.:-(a::GF, b::GF) = GF(a.D, a.a .- b.a)
Base.:-(a::GF) = GF(a.D, -a.a)
Base.:*(a::GF, b::GF) = GF(a.D, a.a .* b.a)
Base.:*(s::Number, a::GF) = GF(a.D, s .* a.a)
Base.:*(a::GF, s::Number) = s * a
Base.:+(s::Number, a::GF) = GF(a.D, s .+ a.a)
Base.:+(a::GF, s::Number) = s + a
Base.:^(a::GF, p::Real) = GF(a.D, a.a .^ p)
mean(a::GF) = sum(a.a) / length(a.a)
meansq(a::GF) = sum(abs2, a.a) / length(a.a)
phasedim(::GF) = 4

function ∂(F::GF, j::Int)
    A = F.a
    N = size(A, 1)
    j == 1 && return GF(F.D, reshape(F.D * reshape(A, N, :), size(A)))
    perm = collect(1:4)
    perm[1], perm[j] = j, 1                    # a transposition, its own inverse
    B = permutedims(A, perm)
    B = reshape(F.D * reshape(B, N, :), size(B))
    GF(F.D, permutedims(B, perm))
end

"(f + eps d)^p to first order."
jpow(F::Jet{GF}, p) = Jet(F.v^p, p * (F.v^(p - 1)) * F.d)
jpow(F::GF, p) = F^p

"Sample fun(x1, x2, v1, v2) on the N^4 spectral grid."
function sample4(N, fun)
    g = spectral_grid(N)
    x = g.nodes
    GF(g.D, [fun(x[i1], x[i2], x[i3], x[i4]) for i1 in 1:N, i2 in 1:N, i3 in 1:N, i4 in 1:N])
end

# 2-D grid fields, for the floating-point positive control (Thm 4.5).  Same code path: the
# array has two singleton trailing dimensions so that GF's 4-D machinery applies unchanged,
# and pb must see phasedim 2 -- so a thin wrapper carries that.
struct GF2
    g::GF
end
Base.:+(a::GF2, b::GF2) = GF2(a.g + b.g)
Base.:-(a::GF2, b::GF2) = GF2(a.g - b.g)
Base.:-(a::GF2) = GF2(-a.g)
Base.:*(a::GF2, b::GF2) = GF2(a.g * b.g)
Base.:*(s::Number, a::GF2) = GF2(s * a.g)
Base.:*(a::GF2, s::Number) = s * a
Base.:+(s::Number, a::GF2) = GF2(s + a.g)
Base.:+(a::GF2, s::Number) = s + a
Base.:^(a::GF2, p::Real) = GF2(a.g^p)
mean(a::GF2) = mean(a.g)
meansq(a::GF2) = meansq(a.g)
phasedim(::GF2) = 2
∂(F::GF2, j::Int) = GF2(∂(F.g, j))
jpow(F::Jet{GF2}, p) = Jet(F.v^p, p * (F.v^(p - 1)) * F.d)
jpow(F::GF2, p) = F^p
function sample2(N, fun)
    g = spectral_grid(N)
    x = g.nodes
    GF2(GF(g.D, reshape([fun(x[i1], x[i2]) for i1 in 1:N, i2 in 1:N], N, N, 1, 1)))
end

# ---------------------------------------------------------------------------
# the bracket, and its Jacobiator by Gateaux derivative
# ---------------------------------------------------------------------------
# {A,B}_c[f] = < c[f] {A_f, B_f} >.  With <c {X_f, C_f}> = <X_f {C_f, c}> (cyclicity, checked
# separately below) one has {X, C}_c = DX[f].(J C_f), J C_f := {C_f, c[f]}, for ANY functional
# X -- in particular for X = {A,B}_c.  So the Jacobiator is a sum of three directional
# derivatives of the functional f -> {A,B}_c[f], which the jets compute without ever forming
# delta{A,B}_c/delta f.  The check is therefore independent of the paper's variation formula.

bracket(c, Af, Bf, F) = mean(c(F) * pb(Af(F), Bf(F)))

"""
    jacobiator(c, Af, Bf, Cf, f) -> (sum, scale)

The cyclic sum {{A,B}_c,C}_c + {{B,C}_c,A}_c + {{C,A}_c,B}_c at the field `f`, and the
largest of the three terms, by which the sum is to be normalised.
"""
function jacobiator(c, Af, Bf, Cf, f)
    η(X) = pb(X(f), c(f))
    term(X, Y, Z) = bracket(c, X, Y, Jet(f, η(Z))).d
    t = (term(Af, Bf, Cf), term(Bf, Cf, Af), term(Cf, Af, Bf))
    return sum(t), maximum(abs, t)
end

normalised(sum, scale) = iszero(scale) ? abs(sum) : abs(sum) / scale

# ---------------------------------------------------------------------------
# the fields
# ---------------------------------------------------------------------------
# One set per dimension, all sparse (so the exact products stay small), all involving every
# coordinate, none symmetric under any exchange.  The functionals are
#     A = (1/2) <phi_A f^2>,    A_f = phi_A f,          second variation: multiplication
#     B = (1/3) <phi_B f^3>,    B_f = phi_B f^2,        second variation: multiplication
#     C = (1/2) <phi_C f>^2 + <psi_C f>,   C_f = phi_C <phi_C f> + psi_C,   nonlocal
# so that the second-variation cancellation in the Jacobiator is exercised in three ways.

# 2-D, coordinates (x, v)
f2 = 3 + cosk((1, 0)) + sink((1, 1), 1 // 2) + cosk((0, 1), 1 // 3)
φA2 = 1 + cosk((0, 1), 1 // 2) + sink((1, -1), 1 // 3)
φB2 = sink((1, 1), 1 // 2) + cosk((1, 0))
φC2 = cosk((1, -1)) + sink((0, 1), 1 // 2)
ψC2 = sink((1, 0)) + cosk((2, 1), 1 // 2)
φz2 = 1 + cosk((1, 0), 1 // 2) + sink((0, 1), 1 // 3)      # a weight that is NOT a function of f

# 4-D, coordinates (x1, x2, v1, v2)
f4 = 3 + cosk((1, 0, 0, 0)) + sink((0, 1, 1, 0), 1 // 2) + cosk((-1, 0, 0, 1), 1 // 3)
φA4 = 1 + cosk((0, 1, 0, 0), 1 // 2) + sink((0, 0, 1, 0), 1 // 3)
φB4 = sink((1, 0, 0, 1), 1 // 2) + cosk((0, 0, 1, 0))
φC4 = cosk((1, -1, 0, 0)) + sink((0, 0, 0, 1), 1 // 2)
ψC4 = sink((1, 0, 0, 0)) + cosk((0, 0, 1, 1), 1 // 2)
φz4 = 1 + cosk((1, 0, 0, -1), 1 // 2) + sink((0, 1, 0, 0), 1 // 3)
# For the additive control c = f^2 + phi(z) the obstruction is 2 <f dA^dB^dC^dphi>, whose zero
# mode needs five wavevectors (one from each factor) to sum to zero: with phi = phi_z4 and the
# sparse fields above none do, and the Jacobiator vanishes by sparsity, not by structure.  One
# extra mode in phi closes a pentagon of wavevectors and makes the control bite.
φadd4 = φz4 + sink((0, 1, 1, 1), 1 // 4)

# 6-D, coordinates (x1, x2, x3, v1, v2, v3); two terms each, the products are already large
f6 = 3 + cosk((1, 0, 0, 0, 0, 1)) + sink((0, 1, 0, 1, 0, 0), 1 // 2) +
     cosk((0, 0, 1, 0, 1, 0), 1 // 3)
φA6 = 1 + cosk((0, 0, 1, 0, 0, 0), 1 // 2) + sink((0, 0, 0, 1, 0, 0), 1 // 3)
φB6 = sink((1, 0, 0, 0, 1, 0), 1 // 2) + cosk((0, 1, 0, 0, 0, 0))
φC6 = cosk((0, 0, 0, 0, 0, 1)) + sink((0, 0, 0, 1, 0, 1), 1 // 2)
ψC6 = sink((0, 0, 1, 0, 0, 0)) + cosk((1, 0, 0, 0, 0, -1), 1 // 2)

functionals(φA, φB, φC, ψC) = (
    F -> φA * F,
    F -> φB * F * F,
    F -> φC * mean(φC * F) + ψC
)
# `phi_C * mean(phi_C * F)` for a Jet: mean returns a Jet of scalars, so define that product.
Base.:*(a::Union{TP, GF, GF2}, s::Jet) = Jet(a * s.v, a * s.d)
Base.:*(s::Jet, a::Union{TP, GF, GF2}) = a * s
Base.:+(a::Jet, b::Union{TP, GF, GF2}) = Jet(a.v + b, a.d)
Base.:+(b::Union{TP, GF, GF2}, a::Jet) = a + b

Af2, Bf2, Cf2 = functionals(φA2, φB2, φC2, ψC2)
Af4, Bf4, Cf4 = functionals(φA4, φB4, φC4, ψC4)
Af6, Bf6, Cf6 = functionals(φA6, φB6, φC6, ψC6)

# ---------------------------------------------------------------------------
# 1. the obstruction integrand: identically zero in 2-D, a divergence in 2n >= 4
# ---------------------------------------------------------------------------

header("1. The obstruction < w(f) ({A_f,B_f}{C_f,f} + cyc) >, exactly, w polynomial in f")

"The cyclic combination P = w(f) ( {A,B}{C,f} + {B,C}{A,f} + {C,A}{B,f} )."
function obstruction(w, A, B, C, f)
    w * (pb(A, B) * pb(C, f) + pb(B, C) * pb(A, f) + pb(C, A) * pb(B, f))
end

"""
    pfaffian_vector(B, C, G) -> V

V^i = T^{ijkl} d_j B d_k C d_l G.  The (k,l) products are formed once each and the (i,j) sums
accumulated before the multiplication by d_j B, which keeps the exact arithmetic affordable
in six dimensions.
"""
function pfaffian_vector(B, C, G)
    n = phasedim(B) ÷ 2
    D = 2n
    dB = [∂(B, i) for i in 1:D]
    dC = [∂(C, i) for i in 1:D]
    dG = [∂(G, i) for i in 1:D]
    W = Dict{Tuple{Int, Int}, typeof(B)}()
    M = [0 * B for _ in 1:D, _ in 1:D]
    for (i, j, k, l, t) in pfaffian_tensor(n)
        w = get!(() -> dC[k] * dG[l], W, (k, l))
        M[i, j] = M[i, j] + t * w
    end
    [sum(M[i, j] * dB[j] for j in 1:D) for i in 1:D]
end

"T^{ijkl} d_i A d_j B d_k C d_l G, the integrand written through the Pfaffian tensor."
function tensor_form(A, B, C, G)
    V = pfaffian_vector(B, C, G)
    sum(∂(A, i) * V[i] for i in eachindex(V))
end

"sum_i d_i ( A V^i ),  V^i = T^{ijkl} d_j B d_k C d_l G: the divergence the integrand equals."
function divergence_form(A, B, C, G)
    V = pfaffian_vector(B, C, G)
    sum(∂(A * V[i], i) for i in eachindex(V))
end

# weights (c')^2 as polynomials in f:  c = f^2 -> 4 f^2;  c = f^{3/2} -> (9/4) f;  and a generic one
weights = (
    ("c = f^2,      (c')^2 = 4 f^2", Q[0, 0, 4]),
    ("c = f^{3/2},  (c')^2 = 9f/4", Q[0, 9 // 4]),
    ("generic       w = 1 + f + f^2/2 + f^3/3", Q[1, 1, 1 // 2, 1 // 3])
)

for (D, f, Af, Bf, Cf) in ((2, f2, Af2, Bf2, Cf2), (4, f4, Af4, Bf4, Cf4), (
    6, f6, Af6, Bf6, Cf6))
    A, B, C = Af(f), Bf(f), Cf(f)
    for (label, wc) in weights
        w = polyval(wc, f)
        G = polyval(antiderivative(wc), f)
        P = obstruction(w, A, B, C, f)
        ms = meansq(P)
        if D == 2
            check("$(D)-D  $label: integrand identically zero (Pluecker)", iszero(ms),
                "mean square " * fmt(ms) * ", $(nterms(P)) terms")
        else
            check("$(D)-D  $label: integrand NOT identically zero", !iszero(ms),
                "mean square " * fmt(ms) * ", $(nterms(P)) terms")
            check("$(D)-D  $label: mean of the integrand exactly zero", iszero(mean(P)),
                "mean " * fmt(mean(P)))
            Pt = tensor_form(A, B, C, G)
            check("$(D)-D  $label: integrand == T^{ijkl} d_iA d_jB d_kC d_lG",
                iszero(meansq(P - Pt)),
                "mean square of the difference " * fmt(meansq(P - Pt)))
            Pd = divergence_form(A, B, C, G)
            check("$(D)-D  $label: integrand == d_i(A_f V^i), a divergence",
                iszero(meansq(P - Pd)),
                "mean square of the difference " * fmt(meansq(P - Pd)))
        end
    end
end

header("1c. Controls that must fail (4-D, c = f^2)")

A4, B4, C4 = Af4(f4), Bf4(f4), Cf4(f4)
w4 = 4 * f4 * f4
Pz = obstruction(φz4 * w4, A4, B4, C4, f4)
check("weight (c')^2 phi(z), phi not a function of f: mean NONZERO", !iszero(mean(Pz)),
    "mean " * fmt(mean(Pz)))
P1 = w4 * pb(A4, B4) * pb(C4, f4)
check(
    "a single term {A,B}{C,f} instead of the cyclic sum: mean NONZERO", !iszero(mean(P1)),
    "mean " * fmt(mean(P1)))
# in 2-D even the z-dependent weight is killed pointwise: 2-D admits c(z, f)
Pz2 = obstruction(φz2 * 4 * f2 * f2, Af2(f2), Bf2(f2), Cf2(f2), f2)
check("2-D: weight phi(z) (c')^2 still identically zero (Pluecker is pointwise)",
    iszero(meansq(Pz2)),
    "mean square " * fmt(meansq(Pz2)))

header("1d. 4-D: the cyclic combination is det(grad A, grad B, grad C, grad f) in (x1,v1,x2,v2)")

"Leibniz determinant of the 4 x 4 matrix whose rows are the gradients in the order (x1, v1, x2, v2)."
function det4(vs...)
    rows = [[∂(v, i) for i in (1, 3, 2, 4)] for v in vs]      # (x1, v1, x2, v2) of (x1, x2, v1, v2)
    r = 0 * vs[1]
    for p in permutations4()
        s = permsign(p)
        r = r + s * (rows[1][p[1]] * rows[2][p[2]] * rows[3][p[3]] * rows[4][p[4]])
    end
    r
end
function permutations4()
    out = NTuple{4, Int}[]
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4
        length(unique((a, b, c, d))) == 4 && push!(out, (a, b, c, d))
    end
    out
end
function permsign(p)
    s = 1
    for i in 1:length(p), j in (i + 1):length(p)

        p[i] > p[j] && (s = -s)
    end
    s
end
S4 = pb(A4, B4) * pb(C4, f4) + pb(B4, C4) * pb(A4, f4) + pb(C4, A4) * pb(B4, f4)
Δ4 = det4(A4, B4, C4, f4)
check("S(dA,dB,dC,df) == det(grad A, grad B, grad C, grad f)", iszero(meansq(S4 - Δ4)),
    "mean square of the difference " * fmt(meansq(S4 - Δ4)) * "; mean square of S " *
    fmt(meansq(S4)))

# ---------------------------------------------------------------------------
# 2. the top-form identity  S(a,b,c,d) mu = a^b^c^d^omega^{n-2}/(n-2)!,  n = 2, 3
# ---------------------------------------------------------------------------

header("2. Top-form identity S(a,b,c,d) omega^n/n! = a^b^c^d^omega^{n-2}/(n-2)!, random rational covectors")

# exterior algebra on R^D over Q: a form is a Dict sorted-index-tuple => coefficient
const Form = Dict{Vector{Int}, Q}
oneform(v) = Form(Dict([i] => Q(v[i]) for i in eachindex(v) if !iszero(v[i])))
function wedge(a::Form, b::Form)
    r = Form()
    for (ia, va) in a, (ib, vb) in b

        isempty(intersect(ia, ib)) || continue
        idx = vcat(ia, ib)
        s = permsign(sortperm(idx))
        key = sort(idx)
        r[key] = get(r, key, zero(Q)) + s * va * vb
    end
    filter!(kv -> !iszero(kv.second), r)
end
function symplectic_form(n)
    ω = Form()
    for k in 1:n
        ω[[k, n + k]] = one(Q)
    end
    ω
end
Pi_pair(n, α, β) = sum(α[k] * β[n + k] - α[n + k] * β[k] for k in 1:n)
function S_pair(n, α, β, γ, δ)
    Pi_pair(n, α, β) * Pi_pair(n, γ, δ) - Pi_pair(n, α, γ) * Pi_pair(n, β, δ) +
    Pi_pair(n, α, δ) * Pi_pair(n, β, γ)
end

rng = MersenneTwister(2026)
for n in (2, 3)
    D = 2n
    α, β, γ, δ = (Q.(rand(rng, -5:5, D)) .// Q.(rand(rng, 1:3, D)) for _ in 1:4)
    ω = symplectic_form(n)
    top = wedge(wedge(wedge(oneform(α), oneform(β)), oneform(γ)), oneform(δ))
    μ = ω
    for _ in 2:n
        μ = wedge(μ, ω)
    end
    μc = get(μ, collect(1:D), zero(Q)) // factorial(n)          # omega^n / n!
    for _ in 1:(n - 2)
        top = wedge(top, ω)
    end
    topc = get(top, collect(1:D), zero(Q)) // factorial(n - 2)   # a^b^c^d^omega^{n-2}/(n-2)!
    S = S_pair(n, α, β, γ, δ)
    check(
        "n = $n: S(a,b,c,d) * (omega^n/n!) == a^b^c^d^omega^{n-2}/(n-2)!", S * μc == topc,
        "S = " * fmt(S) * ", top/mu = " * fmt(topc // μc))
end

# ---------------------------------------------------------------------------
# 3. the full Jacobiator, exactly, in 4-D
# ---------------------------------------------------------------------------

header("3a. Cyclicity <c {X,Y}> = <X {Y,c}> in 4-D, exactly (what turns {{A,B},C} into a Gateaux derivative)")

cyc_lhs = mean(f4 * f4 * pb(A4, B4))
cyc_rhs = mean(A4 * pb(B4, f4 * f4))
check("<f^2 {A_f,B_f}> == <A_f {B_f, f^2}>", cyc_lhs == cyc_rhs, fmt(cyc_lhs) * " == " *
                                                                 fmt(cyc_rhs))
check("<{A_f,B_f}> == 0 (Hamiltonian vector fields are divergence-free)",
    iszero(mean(pb(A4, B4))),
    fmt(mean(pb(A4, B4))))

header("3b. Full Jacobiator in 4-D, exact, polynomial functionals with local and nonlocal second variations")

cases4 = (
    (:lp, "c = f            (Lie-Poisson, control: must vanish)", F -> F, true),
    (:sq, "c = f^2", F -> F * F, true),
    (:cubic, "c = f^3 - f", F -> F * F * F - F, true),
    (:shifted, "c = (f + phi(z))^2   (the class Phi(f + phi(z)))",
        F -> (F + φz4) * (F + φz4), true),
    (:additive, "c = f^2 + phi(z)  (additive z-dependence, not Phi(f + phi): must FAIL)",
        F -> F * F + φadd4, false),
    (:phiz, "c = phi(z) f      (position-dependent weight: must FAIL)",
        F -> φz4 * F, false),
    (:deriv, "c = f d_{x1} f    (derivative-dependent: must FAIL)",
        F -> F * ∂(F, 1), false)
)
jac4 = Dict{Symbol, Tuple{Q, Q}}()
for (key, label, c, expect_zero) in cases4
    t = @elapsed s, scale = jacobiator(c, Af4, Bf4, Cf4, f4)
    jac4[key] = (s, scale)
    ok = expect_zero ? iszero(s) : !iszero(s)
    check(label, ok,
        "Jacobiator " * fmt(s) * ", largest term " * fmt(scale) *
        @sprintf(", normalised %.3e  (%.1fs)", Float64(normalised(s, scale)), t))
end

header("3b'. z-dependent c: the Jacobiator equals the obstruction <c_f S(dA_f, dB_f, dC_f, d c~)>, c~ = c(z, f(z))")

# The reduction of the Jacobiator to the obstruction goes through verbatim for c = c(z, f):
# J psi = {psi, c~} with c~(z) = c(z, f(z)), and the surviving term is c_f S(dA, dB, dC, dc~).
Sfun(A, B, C, G) = pb(A, B) * pb(C, G) + pb(B, C) * pb(A, G) + pb(C, A) * pb(B, G)
O_phiz = mean(φz4 * Sfun(A4, B4, C4, φz4 * f4))
check("c = phi(z) f: Jacobiator == <phi S(dA,dB,dC,d(phi f))>", jac4[:phiz][1] == O_phiz,
    fmt(jac4[:phiz][1]) * " == " * fmt(O_phiz))
O_add = mean(2 * f4 * Sfun(A4, B4, C4, f4 * f4 + φadd4))
check("c = f^2 + phi(z): Jacobiator == <2f S(dA,dB,dC,d(f^2 + phi))>",
    jac4[:additive][1] == O_add,
    fmt(jac4[:additive][1]) * " == " * fmt(O_add))
# and the smallest example in which the additive obstruction 2 <f dA^dB^dC^dphi> is a single
# number: linear functionals with A_f = sin x1, B_f = sin v1, C_f = sin x2, phi = sin v2, and
# f = 3 + cos x1 cos v1 cos x2 cos v2, so that 2 <f dA^dB^dC^dphi> = 2 <cos^2 ...> = 2/16.
let a = sink((1, 0, 0, 0)), b = sink((0, 0, 1, 0)), cc = sink((0, 1, 0, 0)),
    φ = sink((0, 0, 0, 1)),
    u = 3 +
        cosk((1, 0, 0, 0)) * cosk((0, 0, 1, 0)) * cosk((0, 1, 0, 0)) * cosk((0, 0, 0, 1))

    s, scale = jacobiator(
        F -> F * F + φ, F -> a + 0 * F, F -> b + 0 * F, F -> cc + 0 * F, u)
    check("c = f^2 + sin v2, linear functionals: Jacobiator == 2 <f dA^dB^dC^dphi> == 1/8",
        s == 1 // 8,
        "Jacobiator " * fmt(s) * ", largest term " * fmt(scale))
end

header("3c. The same in 2-D: phi(z) f passes there (Thm 4.5 extends to c(z, f) in 2-D)")

for (label, c, expect_zero) in (
    ("2-D  c = f^2", F -> F * F, true),
    ("2-D  c = phi(z) f  (must vanish in 2-D)", F -> φz2 * F, true),
    ("2-D  c = f^2 + phi(z)  (must vanish in 2-D)", F -> F * F + φz2, true),
    ("2-D  c = f d_x f   (derivative-dependent: must FAIL)", F -> F * ∂(F, 1), false)
)
    s, scale = jacobiator(c, Af2, Bf2, Cf2, f2)
    ok = expect_zero ? iszero(s) : !iszero(s)
    check(label, ok, "Jacobiator " * fmt(s) * ", largest term " * fmt(scale))
end

# ---------------------------------------------------------------------------
# 4. floating point, 4-D spectral grid: c = f^{3/2}, refined
# ---------------------------------------------------------------------------

header("4. Full Jacobiator on the 4-D spectral grid, c = f^{3/2} (not band-limited), N = 12, 16, 24, 32")

fld4 = (
    f = (x1, x2, v1, v2) -> 3 + cos(x1) + 0.5sin(x2 + v1) + cos(v2 - x1) / 3,
    φA = (x1, x2, v1, v2) -> 1 + 0.5cos(x2) + sin(v1) / 3,
    φB = (x1, x2, v1, v2) -> 0.5sin(x1 + v2) + cos(v1),
    φC = (x1, x2, v1, v2) -> cos(x1 - x2) + 0.5sin(v2),
    ψC = (x1, x2, v1, v2) -> sin(x1) + 0.5cos(v1 + v2),
    φz = (x1, x2, v1, v2) -> 1 + 0.5cos(x1 - v2) + sin(x2) / 3
)

function float_jacobiators(N)
    f = sample4(N, fld4.f)
    φA, φB, φC, ψC, φz = (sample4(N, fld4[k]) for k in (:φA, :φB, :φC, :ψC, :φz))
    Af, Bf, Cf = functionals(φA, φB, φC, ψC)
    r(c) = normalised(jacobiator(c, Af, Bf, Cf, f)...)
    (lp = r(F -> F), pow32 = r(F -> jpow(F, 1.5)),
        sq = r(F -> F * F), phiz = r(F -> φz * F))
end

resolutions = (12, 16, 24, 32)
results = Dict{Int, Any}()
for N in resolutions
    t = @elapsed results[N] = float_jacobiators(N)
    r = results[N]
    @printf("  N = %2d:  c=f %.2e   c=f^2 %.2e   c=f^{3/2} %.2e   c=phi(z)f %.2e   (%.1fs)\n",
        N, r.lp, r.sq, r.pow32, r.phiz, t)
end
check_exact("c = f, N = 32: Lie-Poisson at round-off (control)", results[32].lp; atol = 1e-11)
check_exact("c = f^2, N = 32: at round-off (band-limited, so exact quadrature)", results[32].sq; atol = 1e-11)
check_refined("c = f^{3/2}: residual falls under refinement 16 -> 32 (or is at round-off)",
    results[16].pow32, results[32].pow32; atol = 1e-11, minrate = 40.0)
# The grid fields are the trigonometric polynomials of section 3, so the floating-point residual
# for c = phi(z) f must reproduce the exact rational one -- and stay there under refinement.
phiz_exact = Float64(normalised(jac4[:phiz]...))
check(
    "c = phi(z) f: residual equals the exact 4-D value of 3b and is flat under refinement (negative control)",
    abs(results[32].phiz - phiz_exact) < 1e-10 &&
        abs(results[16].phiz - phiz_exact) < 1e-10,
    @sprintf("%.6e -> %.6e, exact %.6e", results[16].phiz, results[32].phiz, phiz_exact))

header("4b. 2-D spectral grid, c = f^{3/2}: the positive control for the floating-point code (Thm 4.5)")

function float_jacobiator_2d(N)
    f = sample2(N, (x, v) -> 3 + cos(x) + 0.5sin(x + v) + cos(v) / 3)
    φA = sample2(N, (x, v) -> 1 + 0.5cos(v) + sin(x - v) / 3)
    φB = sample2(N, (x, v) -> 0.5sin(x + v) + cos(x))
    φC = sample2(N, (x, v) -> cos(x - v) + 0.5sin(v))
    ψC = sample2(N, (x, v) -> sin(x) + 0.5cos(2x + v))
    Af, Bf, Cf = functionals(φA, φB, φC, ψC)
    normalised(jacobiator(F -> jpow(F, 1.5), Af, Bf, Cf, f)...)
end
r2 = (float_jacobiator_2d(16), float_jacobiator_2d(32))
check_refined("2-D, c = f^{3/2}: residual 16 -> 32", r2...; atol = 1e-11, minrate = 40.0)

# ---------------------------------------------------------------------------
# 5. the finite-dimensional counterpart: J_ij = sum_m c_ij^m w(z_m) with nonlinear w
# ---------------------------------------------------------------------------

header("5. Finite dimensions: J_ij = sum_m c_ij^m w(z_m), w = z^2, for se(3) and the su(3) sine algebra")

"J and dJ[l,i,j] = d J_ij / d z_l for J_ij = sum_m C[m,i,j] w(z_m)."
function weighted_lie_poisson(C, z, w, w′)
    d = size(C, 1)
    P = [sum(C[m, i, j] * w(z[m]) for m in 1:d) for i in 1:d, j in 1:d]
    dP = [C[l, i, j] * w′(z[l]) for l in 1:d, i in 1:d, j in 1:d]
    P, dP
end
zq = Q.(rand(rng, 1:9, 6)) .// 2
P, dP = weighted_lie_poisson(se3(), zq, z -> z, z -> one(z))
check("se(3), w = z: Jacobi exactly (control)", iszero(jacobi_residual(P, dP)), fmt(jacobi_residual(P, dP)))
P, dP = weighted_lie_poisson(se3(), zq, z -> z^2, z -> 2z)
check("se(3), w = z^2: Jacobi FAILS (Prop. 32-forced)", !iszero(jacobi_residual(P, dP)),
    "normalised residual " * fmt(jacobi_residual(P, dP)))
modes, Csine = sine_algebra(3)
zf = rand(rng, size(Csine, 1)) .+ 0.5
P, dP = weighted_lie_poisson(Csine, zf, z -> z, z -> one(z))
check("su(3) sine algebra, w = z: Jacobi to round-off (control)",
    jacobi_residual(P, dP) < 1e-12,
    @sprintf("%.2e", jacobi_residual(P, dP)))
P, dP = weighted_lie_poisson(Csine, zf, z -> z^2, z -> 2z)
check(
    "su(3) sine algebra, w = z^2: Jacobi FAILS (the truncation does not inherit 'any c')",
    jacobi_residual(P, dP) > 1e-3, @sprintf("normalised residual %.3e",
        jacobi_residual(P, dP)))

summary("III1_structure_function_4d.jl")
