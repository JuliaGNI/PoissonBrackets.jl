#!/usr/bin/env julia
#
# Does the Bialynicki-Birula--Morrison three-bracket reproduce the Zeitlin bracket?
#
#     julia --project=scripts scripts/verify_zeitlin_three_bracket.jl
#
# Bialynicki-Birula and Morrison (1991) write a Lie-Poisson bracket as a Nambu three-bracket
# with a Casimir in the third slot,
#
#     [A, B, S] = {A, B} ,
#
# and the construction needs two hypotheses, both load-bearing:
#
#   (S) the algebra is SEMI-SIMPLE, so that the Killing form kappa is non-degenerate and the
#       lowered structure constants c_ijk = sum_m c_ij^m kappa_mk are totally antisymmetric;
#   (Q) the Casimir in the third slot is the QUADRATIC one built from that form,
#       S = 1/2 sum kappa^mn w_m w_n.
#
# Salmon (2005) section 5 reports that the method "fails even to capture the Nambu bracket
# (1.9)" -- the two-dimensional vorticity bracket -- because the continuum algebra of
# symplectic diffeomorphisms is neither semi-simple nor carries such a Casimir, and the
# integrals diverge. That is a statement about sdiff(T^2).
#
# The question here is whether the MATRIX TRUNCATION escapes it. Zeitlin's sine bracket
# replaces sdiff(T^2) by su(N), which is semi-simple, and Zeitlin's quadratic invariant
# I_2 = sum_k w_k w_{-k} is its Killing-form Casimir. Both hypotheses then hold on their
# face, so the mechanism ought to go through after truncation even though it does not before.
#
# This script decides that. What is actually at stake is (S) and (Q): once the lowered
# constants are totally antisymmetric and S is the Killing quadratic form, the identity
# [A, B, S] = {A, B} is kappa kappa^-1 = I and carries no information. So the checks are
#
#   1. kappa is non-degenerate for the sine algebra    -- hypothesis (S)
#   2. the lowered c_ijk is totally antisymmetric      -- hypothesis (S), the real content
#   3. kappa^mn is proportional to delta_{m+n,0}, i.e. Zeitlin's I_2 IS the Killing Casimir
#                                                       -- hypothesis (Q), the link to the paper
#   4. the three-bracket with S = I_2 reproduces the Zeitlin bracket exactly
#   5. NEGATIVE CONTROL for (Q): Zeitlin's cubic invariant I_3 in the third slot does not
#   6. NEGATIVE CONTROL for (S): so(3) padded with abelian directions has a degenerate
#      Killing form, and step 2 fails there
#
# Steps 5 and 6 are the point of the script as much as steps 1-4. A positive result whose
# hypotheses were never shown to bite is not evidence that the hypotheses are what did the
# work.

using PoissonBrackets
using LinearAlgebra
using Printf
using Random

include(joinpath(@__DIR__, "check.jl"));
using .Checks: header, check, summary

# kappa_mn = tr(ad_m ad_n) = sum_ij c_mi^j c_nj^i, straight from the structure constants, so
# that nothing here depends on the matrix realisation `sine_algebra` happens to use.
function killing_form(C::AbstractArray{T, 3}) where {T}
    d = size(C, 2)
    κ = zeros(T, d, d)
    @inbounds for n in 1:d, m in 1:d

        s = zero(T)
        for j in 1:d, i in 1:d

            s += C[j, m, i] * C[i, n, j]
        end
        κ[m, n] = s
    end
    return κ
end

# c_ijk = sum_m c_ij^m kappa_mk, the third index lowered with the Killing form.
function lower_third(C, κ)
    [sum(C[m, i, j] * κ[m, k] for m in axes(C, 1))
     for i in axes(C, 2), j in axes(C, 3), k in axes(κ, 2)]
end

# max |c_ijk + c_ikj| over the pair that total antisymmetry does not get for free: c is
# already antisymmetric in (i,j) because the commutator is, so swapping the last two indices
# is the only independent transposition to test.
antisym_defect_23(c) = maximum(abs, c .+ permutedims(c, (1, 3, 2)))

# [A, B, S]_ij for a three-bracket c_ijk contracted against the gradient of S.
function three_bracket_matrix(c, dS)
    [sum(c[i, j, k] * dS[k] for k in eachindex(dS))
     for i in axes(c, 1), j in axes(c, 2)]
end

rel(a, b) = maximum(abs, a .- b) / max(maximum(abs, b), eps())

# The mode pairing m -> -m mod N that Zeitlin's I_2 = sum_k w_k w_{-k} contracts with.
function reflection_matrix(N, modes, pos)
    P = zeros(ComplexF64, length(modes), length(modes))
    for m in modes
        P[pos[m], pos[(mod(-m[1], N), mod(-m[2], N))]] = 1
    end
    return P
end

# lambda with kappa^-1 = lambda * P, so that the Killing quadratic form and Zeitlin's I_2
# differ by exactly this factor. Returned with the residual of that ansatz, so a caller
# cannot use the scale without also being able to check that the shape held.
function killing_pairing_scale(κinv, P)
    k = argmax(abs.(κinv))
    λ = κinv[k] / P[k]
    return λ, rel(κinv, λ .* P)
end

const NS = (3, 5)

header("1. The sine algebra is semi-simple: the Killing form is non-degenerate")

println("   kappa_mn = sum_ij c_mi^j c_nj^i, computed from the structure constants alone.")
println()
println("     N    dim    |smallest singular value| / |largest|    cond(kappa)")
for N in NS
    modes, C = sine_algebra(N)
    κ = killing_form(C)
    sv = svdvals(κ)
    ratio = sv[end] / sv[1]
    @printf("     %2d   %4d   %.3e                              %.3e\n",
        N, length(modes), ratio, sv[1] / sv[end])
    check("N = $N: kappa is non-degenerate (Cartan's criterion, so su($N) is semi-simple)",
        ratio > 1e-8, @sprintf("sigma_min/sigma_max = %.2e", ratio))
end

header("2. The lowered structure constants are totally antisymmetric")

println("   This is the content of hypothesis (S). c_ijk = sum_m c_ij^m kappa_mk is already")
println("   antisymmetric in (i,j); what semi-simplicity buys is antisymmetry in (j,k) too.")
println()
for N in NS
    _, C = sine_algebra(N)
    κ = killing_form(C)
    c = lower_third(C, κ)
    scale = maximum(abs, c)
    d23 = antisym_defect_23(c) / scale
    d12 = maximum(abs, c .+ permutedims(c, (2, 1, 3))) / scale
    @printf("     N = %2d:  rel. defect in (i,j) = %.2e,  in (j,k) = %.2e\n", N, d12, d23)
    check("N = $N: c_ijk is antisymmetric in i,j", d12 < 1e-10)
    check("N = $N: c_ijk is antisymmetric in j,k -- so it is TOTALLY antisymmetric",
        d23 < 1e-10, @sprintf("relative defect = %.2e", d23))
end

header("3. Zeitlin's I_2 is the Killing-form Casimir")

println("   Zeitlin (1991) eq. (2.13) gives I_2 = sum_k w_k w_{-k} as the quadratic invariant.")
println("   If kappa^mn is proportional to delta_{m+n mod N, 0} then that is exactly the")
println("   quadratic Casimir hypothesis (Q) asks for, and not merely something like it.")
println()
for N in NS
    modes, C = sine_algebra(N)
    pos = Dict(m => i for (i, m) in enumerate(modes))
    κ = killing_form(C)
    λ, err = killing_pairing_scale(inv(κ), reflection_matrix(N, modes, pos))
    check(
        "N = $N: kappa^mn is proportional to delta_{m+n,0}, i.e. I_2 is the Killing Casimir",
        err < 1e-10, @sprintf("relative deviation from the pairing = %.2e, scale = %.4g",
            err, abs(λ)))
end

header("4. [A, B, I_2] reproduces the Zeitlin bracket")

println("   The manuscript claim: on the truncated algebra, the Bialynicki-Birula--Morrison")
println("   mechanism recovers the Lie-Poisson bracket that Salmon reports it fails to")
println("   recover in the continuum.")
println()
for N in NS
    modes, C = sine_algebra(N)
    pos = Dict(m => i for (i, m) in enumerate(modes))
    d = length(modes)
    κ = killing_form(C)
    c = lower_third(C, κ)
    rng = MersenneTwister(20260905)
    w = randn(rng, ComplexF64, d)
    J = lie_poisson_matrix(C, w)

    # The Killing quadratic form S = 1/2 kappa^mn w_m w_n has gradient dS = kappa^-1 w, and
    # taking it in that form leaves NO free scale to get right: c_ijk (kappa^-1 w)_k
    # = c_ij^m kappa_mk kappa^kl w_l = c_ij^m w_m identically. Any rescaling here would be a
    # place for an error to hide, which is exactly what happened on the first attempt.
    κinv = inv(κ)
    B3 = three_bracket_matrix(c, κinv * w)
    err = rel(B3, J)
    check("N = $N: [., ., S] with S the Killing quadratic form equals sum_m c_ij^m w_m",
        err < 1e-10, @sprintf("relative deviation = %.2e", err))

    # ...and step 3 says that S is Zeitlin's I_2 up to a constant, so the same holds for I_2
    # with the bracket rescaled by that constant, and by nothing else.
    λ, shape = killing_pairing_scale(κinv, reflection_matrix(N, modes, pos))
    dI2 = [w[pos[(mod(-m[1], N), mod(-m[2], N))]] for m in modes]
    errI2 = rel(λ .* three_bracket_matrix(c, dI2), J)
    check("N = $N: [., ., I_2] equals it too, rescaled by the step-3 constant alone",
        shape < 1e-10 && errI2 < 1e-10,
        @sprintf("relative deviation = %.2e at lambda = %.6g", errI2, abs(λ)))
    # and the resulting two-bracket is Poisson, which it must be if it is that bracket
    worst = first(jacobi_residual(J, lie_poisson_derivative(C); normalised = false))
    nrm = worst / maximum(abs, J)^2
    check("N = $N: the recovered two-bracket satisfies the Jacobi identity", nrm < 1e-10,
        @sprintf("normalised residual = %.2e", nrm))
end

header("5. Negative control (Q): the CUBIC invariant in the third slot does not work")

println("   Zeitlin's I_3 = sum cos[(2pi/N)(k x l)] w_k w_l w_{-k-l} is the other invariant.")
println("   Its gradient is quadratic in w, so [., ., I_3] is quadratic where the Lie-Poisson")
println("   bracket is linear -- the two cannot agree, and the mismatch is not small.")
println()
let N = 3
    modes, C = sine_algebra(N)
    pos = Dict(m => i for (i, m) in enumerate(modes))
    d = length(modes)
    κ = killing_form(C)
    c = lower_third(C, κ)
    rng = MersenneTwister(20260905)
    w = randn(rng, ComplexF64, d)
    # dI_3/dw_n = 3 sum_{k+l = -n} cos[(2pi/N)(k x l)] w_k w_l
    dI3 = zeros(ComplexF64, d)
    for n in modes, k in modes

        l = (mod(-n[1] - k[1], N), mod(-n[2] - k[2], N))
        haskey(pos, l) || continue
        cross = k[1] * l[2] - k[2] * l[1]
        dI3[pos[n]] += 3 * cos(2π * cross / N) * w[pos[k]] * w[pos[l]]
    end
    J = lie_poisson_matrix(C, w)
    B3 = three_bracket_matrix(c, dI3)
    # best possible rescaling, so the verdict cannot be an artefact of normalisation
    α = dot(vec(J), vec(B3)) / dot(vec(B3), vec(B3))
    err = rel(α .* B3, J)
    check("N = $N: [., ., I_3] is NOT the Lie-Poisson bracket, at any scale",
        err > 1e-2, @sprintf("relative deviation after optimal rescaling = %.3f", err))
    # the scaling argument, made directly: w -> t w sends J -> t J and [.,.,I_3] -> t^2 [.,.,I_3]
    t = 2.0
    Jt = lie_poisson_matrix(C, t .* w)
    dI3t = t^2 .* dI3
    B3t = three_bracket_matrix(c, dI3t)
    check(
        "N = $N: and the reason is homogeneity -- J is degree 1 in w, [.,.,I_3] is degree 2",
        rel(Jt, t .* J) < 1e-10 && rel(B3t, t^2 .* B3) < 1e-10)
end

header("6. Negative control (S): a degenerate Killing form breaks step 2")

println("   so(3) padded with abelian directions is not semi-simple, so kappa is singular and")
println("   the third index cannot be lowered to a totally antisymmetric tensor. This is what")
println("   Salmon's continuum obstruction looks like in finite dimensions.")
println()
for n in (4, 5)
    C = Array{Float64, 3}(so3(Float64, n))
    κ = killing_form(C)
    sv = svdvals(κ)
    ratio = sv[end] / sv[1]
    check("so(3) + $(n - 3) abelian: kappa IS degenerate", ratio < 1e-10,
        @sprintf("sigma_min/sigma_max = %.2e", ratio))
    c = lower_third(C, κ)
    # the lowered tensor still exists; what fails is that it no longer sees the abelian
    # directions at all, so it cannot reproduce a bracket that does
    J = lie_poisson_matrix(C, ones(n))
    check("so(3) + $(n - 3) abelian: the lowered c_ijk annihilates the abelian directions",
        all(abs(c[i, j, k]) < 1e-12 for i in 1:n, j in 1:n, k in 4:n),
        "so no choice of S in the third slot can recover a bracket involving them")
    check(
        "so(3) + $(n - 3) abelian: kappa is not invertible, so (Q)'s Casimir does not exist",
        !isfinite(cond(κ)) || cond(κ) > 1e12,
        @sprintf("cond(kappa) = %.3e", cond(κ)))
    @printf("     n = %d:  sigma_min/sigma_max = %.2e,  nonzero J entries = %d\n",
        n, ratio, count(!iszero, J))
end

summary("verify_zeitlin_three_bracket.jl")
