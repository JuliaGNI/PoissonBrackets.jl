#!/usr/bin/env julia
#
# Section 5.2: the collision-like metric bracket's coefficients collapse to global moments.
#
#     D_s(u;x) = \int_Omega Q2(grad phi(x) - grad phi(x')) M(x',u(x')) dmu(x')
#     F_s(u;x) = \int_Omega Q2(grad phi(x) - grad phi(x')) w(x')          dmu(x')
#
# with Q2(z) = |z|^2 I_2 - z (x) z and w = grad u + M d_x d_y s.  Written that way each outer
# quadrature node carries a full inner quadrature, so a residual costs O(N_q^2).  The claim
# under test is that it does not have to: Q2 is *quadratic* in z and the kernel factorises, so
# the inner integral reduces to a fixed set of global moments -- each one a single integral
# over Omega, evaluated once -- after which D_s and F_s are pointwise-local in x, i.e. O(N_q).
#
#     julia --project=scripts scripts/verify_metric_collapse.jl
#
# Only LinearAlgebra, Random and Printf are used, so it also runs with no project at all.
#
# The setup is finite-dimensional throughout.  A quadrature rule is just a list of nodes with
# whatever the rule and the measure together put on them, so nothing here needs a grid, a
# basis or the package:
#
#     node k   ->   g_k = grad phi(x_k)        the gradient, *not* perped
#                   c_k = M(x_k) mu_k          the D_s weight  (measure folded in)
#                   v_k = w(x_k) mu_k          the F_s weight  (measure folded in)
#
# Folding mu into the weights is not a simplification of the problem, it is the reason the
# Grad-Shafranov weight dmu = dr dz / r is harmless -- section 7 makes that explicit rather
# than assuming it.
#
# The collapse rests on one identity: in 2D, Q2(z) = z^perp (x) z^perp exactly, where
# z^perp = (-z_2, z_1).  Perping is linear, so
#
#     Q2(grad phi(x) - grad phi(x')) = (alpha - beta) (x) (alpha - beta)
#
# with alpha = grad phi(x)^perp and beta = grad phi(x')^perp, and the x'-integral of a
# quadratic in beta is a moment expansion.  Recentring on the M-weighted mean betabar and
# writing delta = alpha - betabar, gamma = beta - betabar (so that the gamma moment vanishes):
#
#     D_s = m0 delta (x) delta + Sigma          Sigma = int gamma (x) gamma M dmu'
#     F_s = delta (delta . n0) - delta tr(Bc) - Bc delta + Tc
#
#     m0 = int M dmu'   p1 = int beta M dmu'   betabar = p1 / m0
#     n0 = int w dmu'   Bc = int gamma (x) w dmu'   Tc = int gamma (gamma . w) dmu'
#
# Index convention, load-bearing and asserted in section 3 rather than assumed:
# (b (x) w)_ij = b_i w_j and (A v)_i = A_ij v_j.
#
# Sections:
#
#   1. Q2(z) = z^perp (x) z^perp in 2D, exactly, over random z.  Everything else rests on it.
#   2. The collapsed forms against a brute-force O(N_q^2) evaluation, at many outer points.
#      Brute force builds Q2 from its definition |z|^2 I - z (x) z, so this tests the perp
#      identity and the moment expansion together rather than one rearrangement against another.
#   3. The control.  Swap Bc for its transpose in the one place where that matters and the
#      error must become O(1).  A check that cannot fail proves nothing, so this check asserts
#      the sabotaged error is LARGE.
#   4. The moment count.  Sufficiency is settled by perturbing the node data along the null
#      space of the moment map: 14 numbers are sufficient statistics iff a perturbation that
#      moves no moment moves neither D_s nor F_s.  Necessity is settled by the rank of the
#      (moments -> outputs) Jacobian.  The two claimed redundancies -- m2 = tr M2, c1 = tr B --
#      and the claim that n2 and T occur only as n2 - T are each checked on their own.
#   5. Recentring is a correctness requirement, not a refinement.  Sweep grad phi towards a
#      large mean with small spread; the raw moment form loses positive semi-definiteness,
#      which flips the sign of the entropy production, while the centred form cannot.  Over an
#      *ensemble* at each spread, because the raw form's error is the rounding residue of terms
#      that cancel and its sign is a property of the draw, not of the spread.
#   6. Where the collapse breaks.  The same null-space perturbation that leaves the polynomial
#      kernel invariant must move the true Landau kernel |z|^-3 Q2(z) by O(1) -- otherwise the
#      polynomial hypothesis would be decoration.
#   7. The Grad-Shafranov measure dmu = dr dz / r is x-independent, hence absorbed into the
#      weights and harmless; a genuinely x-dependent measure is the control that breaks it.
#
# Every seed is a local MersenneTwister, so a failure is reproducible and no section can
# disturb another's stream.

using LinearAlgebra
using Printf
using Random

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

# ---------------------------------------------------------------------------------------
# the two forms of Q2, and the node data
# ---------------------------------------------------------------------------------------

"Q2 from the definition, `|z|^2 I_2 - z (x) z`.  The brute-force side uses this one."
Q2def(z) = (z[1]^2 + z[2]^2) * I(2) - z * z'

"The perp of a 2-vector, `z^perp = (-z_2, z_1)`."
perp(z) = [-z[2], z[1]]

"""
A quadrature rule reduced to what the bracket sees: gradients, and the two weights with the
measure already folded in.  `g` and `v` are `2 x N`, `c` is length `N`.
"""
struct Rule
    g::Matrix{Float64}
    c::Vector{Float64}
    v::Matrix{Float64}
end

Base.length(r::Rule) = length(r.c)

"A random rule: `N` nodes, positive `M`, arbitrary gradients and arbitrary `w`."
function random_rule(rng, N)
    g = randn(rng, 2, N)
    mu = 0.2 .+ rand(rng, N)               # quadrature weights, positive
    M = 0.5 .+ rand(rng, N)                # M > 0, as eq:M-condition requires
    w = randn(rng, 2, N)
    Rule(g, M .* mu, w .* mu')
end

# ---------------------------------------------------------------------------------------
# brute force: the O(N_q^2) evaluation, one outer point at a time
# ---------------------------------------------------------------------------------------

function D_brute(r::Rule, a; kernel = Q2def)
    D = zeros(2, 2)
    for k in 1:length(r)
        D .+= r.c[k] .* kernel(a .- r.g[:, k])
    end
    return D
end

function F_brute(r::Rule, a; kernel = Q2def)
    F = zeros(2)
    for k in 1:length(r)
        F .+= kernel(a .- r.g[:, k]) * r.v[:, k]
    end
    return F
end

# ---------------------------------------------------------------------------------------
# the collapsed form: 14 global moments, then O(1) work per outer point
# ---------------------------------------------------------------------------------------

"The six centred D_s moments and the eight centred F_s moments, one pass over the nodes."
function centred_moments(r::Rule)
    N = length(r)
    beta = reduce(hcat, [perp(r.g[:, k]) for k in 1:N])

    m0 = sum(r.c)
    p1 = beta * r.c
    bb = p1 ./ m0                                        # betabar

    Sig = zeros(2, 2)
    n0 = zeros(2)
    Bc = zeros(2, 2)
    Tc = zeros(2)
    for k in 1:N
        gam = beta[:, k] .- bb
        vk = r.v[:, k]
        Sig .+= r.c[k] .* (gam * gam')                   # (b (x) w)_ij = b_i w_j
        n0 .+= vk
        Bc .+= gam * vk'
        Tc .+= gam .* dot(gam, vk)
    end
    return (; m0, bb, Sig, n0, Bc, Tc)
end

"`D_s` at the outer point whose gradient is `a`, from the centred moments."
function D_collapsed(mom, a)
    d = perp(a) .- mom.bb
    return mom.m0 .* (d * d') .+ mom.Sig
end

"""
`F_s` at the outer point whose gradient is `a`, from the centred moments.

`transposeB = true` is the section-3 sabotage: it swaps `Bc` for `Bc'` in the one term where
the two differ.  The trace term does not distinguish them, so there is exactly one such place.
"""
function F_collapsed(mom, a; transposeB = false)
    d = perp(a) .- mom.bb
    B = transposeB ? mom.Bc' : mom.Bc
    return d .* dot(d, mom.n0) .- d .* tr(mom.Bc) .- B * d .+ mom.Tc
end

# ---------------------------------------------------------------------------------------
# the raw, uncentred moments -- the form section 5 shows must not be used
# ---------------------------------------------------------------------------------------

"The 14 raw moments `[m0, p1, M2, n0, B, T]`, in the perped variables."
function raw_moments(r::Rule)
    N = length(r)
    beta = reduce(hcat, [perp(r.g[:, k]) for k in 1:N])

    m0 = sum(r.c)
    p1 = beta * r.c
    M2 = zeros(2, 2)
    n0 = zeros(2)
    B = zeros(2, 2)
    T = zeros(2)
    for k in 1:N
        bk, vk = beta[:, k], r.v[:, k]
        M2 .+= r.c[k] .* (bk * bk')
        n0 .+= vk
        B .+= bk * vk'
        T .+= bk .* dot(bk, vk)
    end
    return (; m0, p1, M2, n0, B, T)
end

D_raw(m, a) = (al = perp(a); m.m0 .* (al * al') .- al * m.p1' .- m.p1 * al' .+ m.M2)

F_raw(m, a) = (al = perp(a);
    al .* dot(al, m.n0) .- al .* tr(m.B) .- m.B * al .+ m.T)

# ---------------------------------------------------------------------------------------
# the moment map as a matrix, for the sufficiency and necessity arguments of section 4
# ---------------------------------------------------------------------------------------

"""
The 14-by-3N matrix `A` with `A * [c; v_1; v_2]` the 14 raw moments.  Every moment is a
*linear* functional of the node weights at fixed nodes, which is what makes the null-space
argument in section 4 available at all.
"""
function moment_map(r::Rule)
    N = length(r)
    beta = reduce(hcat, [perp(r.g[:, k]) for k in 1:N])
    A = zeros(14, 3N)
    for k in 1:N
        b1, b2 = beta[1, k], beta[2, k]
        A[1, k] = 1.0                                    # m0
        A[2, k], A[3, k] = b1, b2                        # p1
        A[4, k], A[5, k], A[6, k] = b1^2, b1 * b2, b2^2  # M2 (11, 12, 22)
        kv1, kv2 = N + k, 2N + k
        A[7, kv1], A[8, kv2] = 1.0, 1.0                  # n0
        A[9, kv1], A[10, kv2] = b1, b1                   # B_11, B_12
        A[11, kv1], A[12, kv2] = b2, b2                  # B_21, B_22
        A[13, kv1], A[13, kv2] = b1 * b1, b1 * b2        # T_1
        A[14, kv1], A[14, kv2] = b2 * b1, b2 * b2        # T_2
    end
    return A
end

"Node weights as the single vector the moment map acts on, and back again."
pack(r::Rule) = vcat(r.c, r.v[1, :], r.v[2, :])
function unpack(r::Rule, x)
    N = length(r)
    Rule(r.g, x[1:N], permutedims(hcat(x[(N + 1):(2N)], x[(2N + 1):(3N)])))
end

# =======================================================================================

header("1. Q2(z) = z^perp (x) z^perp in 2D, exactly")

let rng = MersenneTwister(20260905), worst = 0.0
    for _ in 1:20_000
        z = randn(rng, 2) .* exp(4 * randn(rng))         # sweep magnitude over many decades
        zp = perp(z)
        worst = max(worst, opnorm(Q2def(z) .- zp * zp') / max(dot(z, z), 1e-300))
    end
    check("|z|^2 I - z (x) z = z^perp (x) z^perp, 20000 draws over 1e-8..1e8",
        worst < 1e-14, @sprintf("worst relative difference %.2e", worst))
end

header("2. the collapsed form against brute force, N_q = 512")

const NQ = 512
const NOUT = 64

let rng = MersenneTwister(4711)
    r = random_rule(rng, NQ)
    mom = centred_moments(r)
    outer = randn(rng, 2, NOUT)

    wD, wF = 0.0, 0.0
    for j in 1:NOUT
        a = outer[:, j]
        Db, Fb = D_brute(r, a), F_brute(r, a)
        wD = max(wD, opnorm(Db .- D_collapsed(mom, a)) / opnorm(Db))
        wF = max(wF, norm(Fb .- F_collapsed(mom, a)) / norm(Fb))
    end
    check("D_s: collapsed = O(N_q^2) double sum, $NOUT outer points", wD < 1e-12,
        @sprintf("worst relative error %.2e", wD))
    check("F_s: collapsed = O(N_q^2) double sum, $NOUT outer points", wF < 1e-12,
        @sprintf("worst relative error %.2e", wF))

    # The raw (uncentred) form is the same mathematics, so it must agree too -- at these
    # well-scaled gradients.  Section 5 is where the two part company.
    rm = raw_moments(r)
    wDr, wFr = 0.0, 0.0
    for j in 1:NOUT
        a = outer[:, j]
        Db, Fb = D_brute(r, a), F_brute(r, a)
        wDr = max(wDr, opnorm(Db .- D_raw(rm, a)) / opnorm(Db))
        wFr = max(wFr, norm(Fb .- F_raw(rm, a)) / norm(Fb))
    end
    check("the raw moment form agrees too, at O(1) gradients", max(wDr, wFr) < 1e-12,
        @sprintf("worst relative error %.2e", max(wDr, wFr)))
end

header("3. the index convention is load-bearing, not cosmetic")

# The swap is not uniformly O(1), and saying so is the point of reporting the spread rather
# than one number.  Correct minus swapped is exactly (Bc - Bc') delta, which vanishes with
# delta -- so at an outer point that happens to sit near betabar the wrong formula is nearly
# right.  A control built on the *minimum* over sampled points would therefore fail for a
# reason that has nothing to do with the convention.  The claim that survives is the one
# below: typical error O(1), and an exact algebraic difference that cannot be zero unless Bc
# is symmetric.

let rng = MersenneTwister(97)
    r = random_rule(rng, NQ)
    mom = centred_moments(r)

    worst_ok, errs, ident = 0.0, Float64[], 0.0
    for _ in 1:2000
        a = randn(rng, 2)
        Fb = F_brute(r, a)
        Fok = F_collapsed(mom, a)
        Fbad = F_collapsed(mom, a; transposeB = true)
        worst_ok = max(worst_ok, norm(Fb .- Fok) / norm(Fb))
        push!(errs, norm(Fb .- Fbad) / norm(Fb))
        d = perp(a) .- mom.bb
        ident = max(ident,
            norm((Fbad .- Fok) .- (mom.Bc .- mom.Bc') * d) / max(norm(Fbad .- Fok), 1e-300))
    end
    sort!(errs)
    med = errs[length(errs) ÷ 2]

    check("with (b (x) w)_ij = b_i w_j the collapse is exact", worst_ok < 1e-12,
        @sprintf("worst relative error %.2e over 2000 outer points", worst_ok))
    check("swapping Bc for its transpose gives an O(1) error typically",
        med > 0.1 && errs[end] > 0.5,
        @sprintf("relative error min %.2e  median %.2e  max %.2e; %.1f%% exceed 0.1",
            errs[1], med, errs[end], 100 * count(>(0.1), errs) / length(errs)))
    check("and it is never a silent no-op: swapped - correct = (Bc - Bc') delta exactly",
        ident < 1e-12 && opnorm(mom.Bc .- mom.Bc') > 0.1 * opnorm(mom.Bc),
        @sprintf("identity holds to %.2e, ||Bc - Bc'|| / ||Bc|| = %.3f", ident,
            opnorm(mom.Bc .- mom.Bc') / opnorm(mom.Bc)))

    # Where the two do *not* differ, so that "one place" is a measurement and not a guess.
    check("tr(Bc) = tr(Bc'), so the trace term is not a second such place",
        abs(tr(mom.Bc) - tr(mom.Bc')) <= 1e-14 * abs(tr(mom.Bc)),
        @sprintf("|difference| %.2e", abs(tr(mom.Bc) - tr(mom.Bc'))))
end

header("4. the moment count is 14")

let rng = MersenneTwister(31337)
    r = random_rule(rng, NQ)
    rm = raw_moments(r)

    # (a) and (b): two quantities that look like separate accumulators and are not.
    m2 = sum(r.c[k] * dot(perp(r.g[:, k]), perp(r.g[:, k])) for k in 1:NQ)
    c1 = sum(dot(perp(r.g[:, k]), r.v[:, k]) for k in 1:NQ)
    check("m2 = int |beta|^2 M dmu' is tr(M2), not a 15th moment",
        abs(m2 - tr(rm.M2)) <= 1e-12 * abs(m2), @sprintf("|m2 - tr M2| %.2e",
            abs(m2 - tr(rm.M2))))
    check("c1 = int (beta . w) dmu' is tr(B), not a 15th moment",
        abs(c1 - tr(rm.B)) <= 1e-12 * max(abs(c1), 1.0),
        @sprintf("|c1 - tr B| %.2e", abs(c1 - tr(rm.B))))

    # (c) n2 and T occur only as n2 - T.  Expanding Q2 from its *definition* rather than
    # through the perp identity produces ten F_s accumulators: n0, Braw (4), n2 (2), Traw (2),
    # in the unperped gradients.  Two of them are redundant because they enter only as a
    # difference, which is precisely why the perp route reaches eight directly.
    n0 = zeros(2)
    Braw = zeros(2, 2)
    n2 = zeros(2)
    Traw = zeros(2)
    for k in 1:NQ
        bk, vk = r.g[:, k], r.v[:, k]
        n0 .+= vk
        Braw .+= bk * vk'
        n2 .+= dot(bk, bk) .* vk
        Traw .+= bk .* dot(bk, vk)
    end
    F_ten(a, n2_, T_) = dot(a, a) .* n0 .- 2 .* (Braw' * a) .+ n2_ .-
                        a .* dot(a, n0) .+ a .* tr(Braw) .+ Braw * a .- T_

    worst, worst_shift = 0.0, 0.0
    shift = randn(rng, 2)
    for _ in 1:NOUT
        a = randn(rng, 2)
        Fb = F_brute(r, a)
        worst = max(worst, norm(Fb .- F_ten(a, n2, Traw)) / norm(Fb))
        worst_shift = max(worst_shift,
            norm(F_ten(a, n2, Traw) .- F_ten(a, n2 .+ shift, Traw .+ shift)) / norm(Fb))
    end
    check("the ten-accumulator expansion of Q2 from its definition reproduces F_s",
        worst < 1e-12, @sprintf("worst relative error %.2e", worst))
    check("n2 and T enter only as n2 - T: a common shift changes nothing",
        worst_shift < 1e-14, @sprintf("worst relative change %.2e", worst_shift))

    # Sufficiency: 14 numbers are sufficient statistics for the node data iff moving the data
    # along the null space of the moment map moves neither D_s nor F_s -- measured on the
    # brute-force side, so the moment formulas are not being used to prove themselves.
    A = moment_map(r)
    x = pack(r)
    y = randn(rng, 3NQ)
    dx = y .- A' * ((A * A') \ (A * y))                  # project out the row space
    dx .*= norm(x) / norm(dx)                            # an O(1) perturbation of the data
    rp = unpack(r, x .+ dx)

    worstD, worstF = 0.0, 0.0
    for _ in 1:NOUT
        a = randn(rng, 2)
        Db, Fb = D_brute(r, a), F_brute(r, a)
        worstD = max(worstD, opnorm(Db .- D_brute(rp, a)) / opnorm(Db))
        worstF = max(worstF, norm(Fb .- F_brute(rp, a)) / norm(Fb))
    end
    check("the perturbation really is in the null space of the 14 moments",
        norm(A * dx) <= 1e-9 * norm(A) * norm(dx),
        @sprintf("||A dx|| / (||A|| ||dx||) %.2e", norm(A * dx) / (norm(A) * norm(dx))))
    check(
        "SUFFICIENT: a null-space perturbation of the node data moves neither D_s nor F_s",
        max(worstD, worstF) < 1e-10,
        @sprintf("worst relative change %.2e over %d outer points", max(worstD, worstF),
            NOUT))

    # Necessity: D_s and F_s are linear in the 14 moments, so the (moments -> outputs) map is
    # a constant matrix.  Rank 14 says no moment can be dropped or expressed through the
    # others; a smaller rank would name the redundancy.
    outer = randn(rng, 2, 40)
    basis = [(m0 = e[1], p1 = e[2:3], M2 = [e[4] e[5]; e[5] e[6]],
                 n0 = e[7:8], B = [e[9] e[10]; e[11] e[12]], T = e[13:14])
             for e in eachcol(Matrix{Float64}(I, 14, 14))]
    J = reduce(hcat,
        [vcat([begin
                   D = D_raw(m, outer[:, j])
                   vcat(D[1, 1], D[1, 2], D[2, 2], F_raw(m, outer[:, j]))
               end for j in 1:40]...) for m in basis])
    s = svdvals(J)
    check("NECESSARY: the (14 moments -> D_s, F_s) Jacobian has full rank 14",
        length(s) == 14 && s[14] > 1e-8 * s[1],
        @sprintf("sigma_14 / sigma_1 = %.2e over 40 outer points", s[14] / s[1]))
    check("the established count is 14 = 6 (D_s) + 8 (F_s)", size(A, 1) == 14,
        @sprintf("%d moments: m0(1) p1(2) M2(3) | n0(2) B(4) T(2)", size(A, 1)))
end

header("5. recentring is a correctness requirement: PSD under a large mean")

# Whether a *particular* draw loses definiteness is luck: the raw form's error is a rounding
# residue of terms that cancel, so its sign varies with the draw.  A single trajectory of grad
# phi therefore measures nothing: one fixed draw can survive down to a spread of 1e-8 and
# report "never negative".  The sweep is over an ensemble, and reports the *rate* at which
# definiteness is lost alongside the accuracy loss, which is deterministic.

const NDRAW = 200

let spreads = [10.0^(-e) for e in 2:10]
    raw_err, raw_bad = Float64[], Int[]
    cen_err, cen_bad = Float64[], Int[]

    for eps in spreads
        rng = MersenneTwister(2718)                      # same ensemble at every spread
        re, ce, rb, cb = 0.0, 0.0, 0, 0
        for _ in 1:NDRAW
            gbar = 3.0 .* randn(rng, 2)                  # a large mean gradient
            xi = randn(rng, 2, NQ)
            mu = 0.2 .+ rand(rng, NQ)
            M = 0.5 .+ rand(rng, NQ)
            w = randn(rng, 2, NQ)
            aq = randn(rng, 2)                           # the outer point's own fluctuation

            r = Rule(gbar .+ eps .* xi, M .* mu, w .* mu')
            a = gbar .+ eps .* aq

            Db = D_brute(r, a)                           # accurate: differences taken first
            Dr = D_raw(raw_moments(r), a)
            Dc = D_collapsed(centred_moments(r), a)
            re = max(re, opnorm(Dr .- Db) / opnorm(Db))
            ce = max(ce, opnorm(Dc .- Db) / opnorm(Db))
            lr = eigvals(Symmetric((Dr .+ Dr') ./ 2))
            lc = eigvals(Symmetric((Dc .+ Dc') ./ 2))
            lr[1] < -1e-12 * abs(lr[2]) && (rb += 1)
            lc[1] < -1e-12 * abs(lc[2]) && (cb += 1)
        end
        push!(raw_err, re)
        push!(raw_bad, rb)
        push!(cen_err, ce)
        push!(cen_bad, cb)
    end

    println("     spread   raw rel.err   raw not PSD   centred rel.err   centred not PSD")
    for (i, eps) in enumerate(spreads)
        @printf("     %.0e   %11.2e   %5d/%d      %11.2e   %8d/%d\n",
            eps, raw_err[i], raw_bad[i], NDRAW, cen_err[i], cen_bad[i], NDRAW)
    end

    first_neg = findfirst(>(0), raw_bad)
    check("the RAW moment form loses positive semi-definiteness", first_neg !== nothing,
        first_neg === nothing ? "never, over the whole sweep -- UNEXPECTED" :
        @sprintf("first at spread %.0e (%d/%d draws), rising to %d/%d at %.0e",
            spreads[first_neg], raw_bad[first_neg], NDRAW, maximum(raw_bad), NDRAW,
            spreads[argmax(raw_bad)]))
    check("the RAW form loses all relative accuracy first", maximum(raw_err) > 1.0,
        @sprintf("relative error reaches 1 at spread %.0e, %.2e at %.0e",
            spreads[findfirst(>(1.0), raw_err)], maximum(raw_err), spreads[end]))
    check("the CENTRED form is PSD in every draw at every spread", all(iszero, cen_bad),
        @sprintf("%d/%d failures over %d draws", sum(cen_bad), NDRAW * length(spreads),
            NDRAW * length(spreads)))
    check("and its error degrades gracefully rather than catastrophically",
        maximum(cen_err) < 1e-3,
        @sprintf("worst relative error %.2e over the whole sweep", maximum(cen_err)))
end

header("6. the collapse needs a polynomial kernel: the Landau kernel breaks it")

let rng = MersenneTwister(60606)
    r = random_rule(rng, NQ)
    landau(z) = Q2def(z) ./ (sqrt(dot(z, z))^3)

    A = moment_map(r)
    x = pack(r)
    y = randn(rng, 3NQ)
    dx = y .- A' * ((A * A') \ (A * y))
    dx .*= norm(x) / norm(dx)
    rp = unpack(r, x .+ dx)

    poly_change, landau_change, closest = 0.0, 0.0, Inf
    for _ in 1:NOUT
        a = randn(rng, 2)
        closest = min(closest, minimum(norm(a .- r.g[:, k]) for k in 1:NQ))
        Db = D_brute(r, a)
        poly_change = max(poly_change, opnorm(Db .- D_brute(rp, a)) / opnorm(Db))
        Dl = D_brute(r, a; kernel = landau)
        landau_change = max(landau_change,
            opnorm(Dl .- D_brute(rp, a; kernel = landau)) / opnorm(Dl))
    end
    check(
        "the same perturbation leaves the quadratic kernel invariant", poly_change < 1e-10,
        @sprintf("worst relative change %.2e", poly_change))
    check("it moves the Landau kernel |z|^-3 Q2(z) by O(1)", landau_change > 0.1,
        @sprintf("worst relative change %.2e (closest node %.2e)", landau_change, closest))
end

header("7. the Grad-Shafranov measure dmu = dr dz / r is absorbed into the weights")

let rng = MersenneTwister(881)
    # C1's rectangle, [1,7] x [-9.5,9.5], where 1/r is bounded and the weight is genuinely
    # non-Lebesgue. Only r enters dmu = dr dz / r, so the z coordinate is never needed: the
    # rectangle fixes the range of r and nothing else.
    rr = 1.0 .+ 6.0 .* rand(rng, NQ)
    q = 0.2 .+ rand(rng, NQ)                             # a tensor Gauss rule's own weights
    mu = q ./ rr                                         # dmu = dr dz / r
    g = randn(rng, 2, NQ)
    M = 0.5 .+ rand(rng, NQ)
    w = randn(rng, 2, NQ)
    r = Rule(g, M .* mu, w .* mu')

    mom = centred_moments(r)
    worstD, worstF = 0.0, 0.0
    for _ in 1:NOUT
        a = randn(rng, 2)
        Db, Fb = D_brute(r, a), F_brute(r, a)
        worstD = max(worstD, opnorm(Db .- D_collapsed(mom, a)) / opnorm(Db))
        worstF = max(worstF, norm(Fb .- F_collapsed(mom, a)) / norm(Fb))
    end
    check("the collapse is still exact under dmu = dr dz / r", max(worstD, worstF) < 1e-12,
        @sprintf("worst relative error %.2e, weights spanning %.2f..%.2f",
            max(worstD, worstF), minimum(mu), maximum(mu)))

    # The control: it is x-independence that does the work, not the weight being harmless.
    # Give the measure a dependence on the outer point and one fixed set of moments cannot
    # serve every x -- which is what "dmu(x') must not depend on x" actually forbids.
    worst_dep = 0.0
    for _ in 1:NOUT
        a = randn(rng, 2)
        s = 1.0 ./ (rr .+ abs(a[1]))                     # an x-dependent measure
        rx = Rule(g, M .* q .* s, w .* (q .* s)')
        Db = D_brute(rx, a)
        worst_dep = max(worst_dep, opnorm(Db .- D_collapsed(mom, a)) / opnorm(Db))
    end
    check("an x-DEPENDENT measure does break the collapse", worst_dep > 0.1,
        @sprintf("worst relative error %.2e", worst_dep))
end

summary("verify_metric_collapse.jl")
