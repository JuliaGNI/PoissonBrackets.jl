#!/usr/bin/env julia
#
# III.4 -- Two extensions of the structure-function bracket  {A,B}_c = <c(f) {A_f, B_f}>  on
# T^{2n}: several field variables, and structure functions depending on derivatives of f.
#
#     julia --project=scripts scripts/fable/III4_multifield_derivative.jl
#
# (a) Several fields f^1..f^m,  {A,B} = < c^{ij}(f) {A_i, B_j} >,  c^{ij} = c^{ji}.  The Jacobiator
#     reduces (cyclicity, skew-adjointness, chain rule -- as in III.1) to
#
#         O = sum_cyc < Theta^{ijl}_m {A_i,B_j} {C_l, f^m} >,   Theta^{ijl}_m = d_k c^{ij} d_m c^{kl},
#
#     and the three cyclic terms carry three DIFFERENT components of Theta.  Jacobi for all f iff
#       (A) Theta is totally symmetric in (i,j,l)  -- the product (xi * eta)_k = d_k c^{ij} xi_i eta_j
#           on the cotangent space of the target is associative --      [every dimension]
#       (B) the target one-forms Theta^{ijl}_m du^m are closed          [2n >= 4 only].
#     Checked exactly: the necessity lemma (the three functionals <{a,b}{c,g}>, <{b,c}{a,g}>,
#     <{c,a}{b,g}> obey only the one relation "sum = 0"), and the full Jacobiator by jets for
#     positive cases (decoupled species, the coupled g = f^1 + k, the slaved family
#     c^{ij} = int p^i p^j du^1, the frozen family a(f^1, f^2) with a = Phi(u1 + phi(u2))) against
#     controls that must fail (g = f^1 f^2, g = (f^1)^2 / 2, c = Hess(u1^2 u2 / 2), and the
#     frozen a = u1 u2 which passes in 2-D and fails in 4-D exactly as (B) predicts).
# (b) c = c(f, grad f):  O = < c_f S(dA,dB,dC,dc) > + sum_cyc < {A_f,B_f} c_{,a} d_a {C_f, c} >,
#     whose second term has a principal symbol forcing c_{,a} = 0: no first-order dependence
#     survives, in any dimension.  What is checked here is the CONCLUSION, instance by instance,
#     not that split -- `obstruction_general` forms the unsplit sum, and nothing isolates the two
#     terms.  The instances are c = f_x f_v (the paper's eq:sym-sym-alternative), |grad f|^2,
#     {f, g}, f + {f, g}, and the second-order c = Laplace f, f Laplace f, f_{x1 x1}; the formula
#     for O is checked against the jet Jacobiator; the Casimirs int kappa(f) are shown lost.
# (c) The multi-field 4-bracket < T^{ijkl} (A_i B_j {C_k, D_l} - C_k D_l {A_i, B_j}) > generates
#     c^{ik} = T^{ijkl} s_j s_l from S = int s(f), checked exactly.
#
# All means <.> are (2pi)^{-2n} times the integral, so no pi appears.  Exact over Q(i) throughout.

using Printf

include(joinpath(@__DIR__, "..", "check.jl"));
using .Checks: header, check, summary, fmt
include(joinpath(@__DIR__, "fabletools.jl"))

# ---------------------------------------------------------------------------
# several fields: the bracket and its Jacobiator by Gateaux derivative
# ---------------------------------------------------------------------------
# A field configuration is a Vector of m fields (TP or Jet).  A structure function c(F) returns
# the symmetric m x m Matrix c^{ij}(F); a functional derivative A(F) returns the Vector A_i(F).
# {X, C} = DX[f].(J C) with (J C)^i = sum_j {C_j, c^{ij}} for ANY functional X (cyclicity), so
# the Jacobiator is a sum of three directional derivatives of f -> {A,B}[f].

function mbracket(c, A, B, F)
    C = c(F)
    a, b = A(F), B(F)
    m = length(F)
    mean(sum(C[i, j] * pb(a[i], b[j]) for i in 1:m, j in 1:m))
end

function mjacobiator(c, A, B, C, f)
    m = length(f)
    Cf = c(f)
    η(X) = (x = X(f); [sum(pb(x[j], Cf[i, j]) for j in 1:m) for i in 1:m])
    term(X, Y, Z) = (d = η(Z); mbracket(c, X, Y, [Jet(f[i], d[i]) for i in 1:m]).d)
    t = (term(A, B, C), term(B, C, A), term(C, A, B))
    return sum(t), maximum(abs, t)
end

# Convert before dividing: this is the reporting path, and a Rational{Int128} quotient of two
# perfectly good residuals can overflow while being printed.
normalised(s, scale) = iszero(scale) ? Float64(abs(s)) : Float64(abs(s)) / Float64(scale)

# ---------------------------------------------------------------------------
# fields: sparse, involving every coordinate, no exchange symmetry
# ---------------------------------------------------------------------------
# Functionals with local and nonlocal second variations, exercising the cancellation:
#   A = 1/2 <phiA (f1)^2> + <psiA f1 f2>          A_1 = phiA f1 + psiA f2,  A_2 = psiA f1
#   B = 1/3 <phiB (f2)^3> + <chiB f1 (f2)^2>      B_1 = chiB (f2)^2,  B_2 = phiB (f2)^2 + 2 chiB f1 f2
#   C = 1/2 <phiC f1>^2 + <psiC f2><chiC f1>      C_1 = phiC<phiC f1> + chiC<psiC f2>,  C_2 = psiC<chiC f1>

function mfunctionals(φA, ψA, φB, χB, φC, ψC, χC)
    A(F) = [φA * F[1] + ψA * F[2], ψA * F[1]]
    B(F) = [χB * F[2] * F[2], φB * F[2] * F[2] + 2 * χB * F[1] * F[2]]
    C(F) = [φC * mean(φC * F[1]) + χC * mean(ψC * F[2]), ψC * mean(χC * F[1])]
    A, B, C
end
# The 4-D fields below are dense enough that the cubic B is unaffordable exactly; there B is
# quadratic, B = <chiB f1 f2> + 1/2 <phiB (f2)^2>, with the constant Hessian [0 chiB; chiB phiB].
function mfunctionals_quadratic(φA, ψA, φB, χB, φC, ψC, χC)
    A(F) = [φA * F[1] + ψA * F[2], ψA * F[1]]
    B(F) = [χB * F[2], χB * F[1] + φB * F[2]]
    C(F) = [φC * mean(φC * F[1]) + χC * mean(ψC * F[2]), ψC * mean(χC * F[1])]
    A, B, C
end

# 2-D, coordinates (x, v)
f2 = [3 + cosk((1, 0)) + sink((1, 1), 1 // 2) + cosk((0, 1), 1 // 3),
    2 + sink((1, 0), 1 // 2) + cosk((1, -1), 1 // 3) + sink((0, 2), 1 // 5)]
mf2 = mfunctionals(1 + cosk((0, 1), 1 // 2) + sink((1, -1), 1 // 3),   # phiA
    sink((1, 0)) + cosk((0, 1), 1 // 4),                                 # psiA
    sink((1, 1), 1 // 2) + cosk((1, 0)),                                 # phiB
    cosk((0, 1)) + sink((2, 1), 1 // 3),                                 # chiB
    cosk((1, -1)) + sink((0, 1), 1 // 2),                                # phiC
    sink((1, 0)) + cosk((2, 1), 1 // 2),                                 # psiC
    cosk((1, 1), 1 // 3) + sink((0, 1)))                                 # chiC

# 4-D, coordinates (x1, x2, v1, v2).  III.1's trap: with fields too sparse the zero mode of a
# product of five factors needs five wavevectors summing to zero and none may -- an exact zero
# for combinatorial reasons.  A second trap sits on top: with pure sines and cosines the closed
# pentagon's coefficient is imaginary for an odd number of sines and cancels against its
# conjugate.  Every 4-D field therefore carries a mode with all four components nonzero AND a
# generic phase (sin and cos), and every negative control below is shown to be nonzero.
gen(k, amp) = cosk(k, amp) + sink(k, amp // 2)
# products of three one-coordinate factors: all 8 sign patterns of (1,1,1,0), resp. (1,0,1,1),
# at once (four factors would double every mode count and the cubic cases take minutes)
e1, e2, e3, e4 = (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)
prodc = cosk(e1) * cosk(e2) * cosk(e3)
prods = sink(e1) * cosk(e3) * sink(e4)
f4 = [
    3 + cosk((1, 0, 0, 0)) + sink((0, 1, 1, 0), 1 // 2) + cosk((-1, 0, 0, 1), 1 // 3) +
    gen((1, 1, 1, 1), 1 // 4) + (1 // 4) * prodc,
    2 + sink((0, 0, 1, 0), 1 // 2) + cosk((1, 0, 0, -1), 1 // 3) +
    sink((0, 1, 0, 1), 1 // 5) +
    gen((1, -1, 1, 1), 1 // 4) + (1 // 4) * prods]
# The functional weights enter linearly, so they can afford a three-factor product each (the
# quadratic 4-D functionals below have fewer f-factors per pentagon, and the weights must close it).
mf4 = mfunctionals_quadratic(
    1 + cosk((0, 1, 0, 0), 1 // 2) + sink((0, 0, 1, 0), 1 // 3) +
    gen((1, 1, -1, 1), 1 // 4) +
    (1 // 3) * cosk(e1) * sink(e2) * cosk(e4),
    sink((1, 0, 0, 0)) + cosk((0, 0, 0, 1), 1 // 4) + gen((1, 1, 1, -1), 1 // 5) +
    (1 // 3) * sink(e2) * cosk(e3) * cosk(e4),
    sink((1, 0, 0, 1), 1 // 2) + cosk((0, 0, 1, 0)) + gen((-1, 1, 1, 1), 1 // 4) +
    (1 // 3) * cosk(e1) * cosk(e2) * sink(e4),
    cosk((0, 1, 0, 0)) + sink((0, 0, 1, 1), 1 // 3) + gen((1, 1, 1, 0), 1 // 4) +
    (1 // 3) * sink(e1) * sink(e3) * cosk(e4),
    # phiC, psiC, chiC share a mode with f1, f2, f1 respectively: the nonlocal C_l(f) must not
    # vanish at the base point, see the check below
    cosk((1, -1, 0, 0)) + sink((0, 0, 0, 1), 1 // 2) + gen((1, 0, 1, 1), 1 // 3) +
    (1 // 3) * cosk(e1) * cosk(e3) * cosk(e4) + cosk((1, 0, 0, 0), 1 // 2),
    sink((1, 0, 0, 0)) + cosk((0, 0, 1, 1), 1 // 2) + gen((1, 1, 1, 1), 1 // 3) +
    (1 // 3) * sink(e1) * cosk(e2) * sink(e3) + sink((0, 0, 1, 0), 1 // 3),
    cosk((0, 1, 1, 0), 1 // 3) + sink((0, 0, 0, 1)) + gen((1, 1, 0, -1), 1 // 4) +
    (1 // 3) * cosk(e2) * sink(e3) * cosk(e4) + cosk((-1, 0, 0, 1), 1 // 3))

# ---------------------------------------------------------------------------
# 1. the necessity lemma: J1 + J2 + J3 = 0 is the ONLY relation among the three functionals
# ---------------------------------------------------------------------------

header("1. J1 = <{a,b}{c,g}>, J2 = <{b,c}{a,g}>, J3 = <{c,a}{b,g}>: each nonzero, sum exactly zero")

# 4-D: a, b, c, g carry the wavevectors (1,0,0,1), (0,0,1,0), (0,1,0,0), -(1,1,1,1), which sum
# to zero with {a,b}, {c,g} both nonzero; generic phases, plus generic modes (also used in 2c)
la4 = sink((1, 0, 0, 1)) + cosk((0, 1, 1, 0), 1 // 2)
lb4 = sink((0, 0, 1, 0)) + sink((1, 0, 0, -1), 1 // 3)
lc4 = sink((0, 1, 0, 0)) + cosk((0, 1, 1, 1), 1 // 2)
lg4 = gen((1, 1, 1, 1), 1) + gen((1, 2, 2, 0), 1 // 3) + sink((1, 0, 1, 0), 1 // 4)
for (D, a, b, c, g) in (
    (2, sink((1, 0)) + cosk((1, 1), 1 // 2), cosk((0, 1)) + sink((1, -1), 1 // 3),
    cosk((1, 0), 1 // 2) + sink((0, 1)), sink((1, 1)) + cosk((2, -1), 1 // 4)),
    (4, la4, lb4, lc4, lg4))
    J1, J2, J3 = mean(pb(a, b) * pb(c, g)), mean(pb(b, c) * pb(a, g)),
    mean(pb(c, a) * pb(b, g))
    check("$(D)-D: J1, J2, J3 individually nonzero",
        !iszero(J1) && !iszero(J2) && !iszero(J3),
        fmt(J1) * ", " * fmt(J2) * ", " * fmt(J3))
    check("$(D)-D: J1 + J2 + J3 == 0 (Pluecker in 2-D, exact top form in 4-D)",
        iszero(J1 + J2 + J3), fmt(J1 + J2 + J3))
    # The two specialisations of the proof, each evaluated in all three components rather than
    # asserted in a label.  c = a, g = b sends (J1,J2,J3) to (S, -S, 0) with S = <{a,b}^2>;
    # b = a, g = c sends it to (0, T, -T) with T = <{a,c}^2>.
    S = (mean(pb(a, b) * pb(a, b)), mean(pb(b, a) * pb(a, b)), mean(pb(a, a) * pb(b, b)))
    T = (mean(pb(a, a) * pb(c, c)), mean(pb(a, c) * pb(a, c)), mean(pb(c, a) * pb(a, c)))
    check("$(D)-D: c = a, g = b gives (J1, J2, J3) = (S, -S, 0), S = <{a,b}^2> > 0",
        S[1] > 0 && S[2] == -S[1] && iszero(S[3]), "(" * join(fmt.(S), ", ") * ")")
    check("$(D)-D: b = a, g = c gives (J1, J2, J3) = (0, T, -T), T = <{a,c}^2> > 0",
        T[2] > 0 && T[3] == -T[2] && iszero(T[1]), "(" * join(fmt.(T), ", ") * ")")
    # This is what "exactly one relation" needs: two independent points of the plane J1+J2+J3 = 0,
    # so the plane is the whole of it and no second relation constrains the three functionals.
    check(
        "$(D)-D: the two specialisations are linearly independent -- one relation, not two",
        !iszero(S[1] * T[2] - S[2] * T[1]), "2x2 minor " * fmt(S[1] * T[2] - S[2] * T[1]))
end

# ---------------------------------------------------------------------------
# 2. several fields: the full Jacobiator, exactly, in 2-D and 4-D
# ---------------------------------------------------------------------------

header("2. Several fields: full Jacobiator by jets, m = 2, exact  (A) associativity, (B) closedness")

# A trap this section has to avoid: if every mean <phiC f1>, <psiC f2>, <chiC f1> vanishes, C_l(f) = 0 at the
# base point, and then the obstruction is zero for EVERY bracket -- each cyclic term carries a
# factor C_l -- while the individual Jacobiator terms need not be.  Unless phiC, psiC and chiC each
# carry a mode of f, all four 4-D negatives below give an exact zero for this reason.
let Cf = mf4[3](f4)
    check(
        "4-D: the nonlocal functional derivative C_l(f) is not identically zero at the base point",
        !iszero(meansq(Cf[1])) && !iszero(meansq(Cf[2])),
        "mean squares " * fmt(meansq(Cf[1])) * ", " * fmt(meansq(Cf[2])))
end

# structure functions c^{ij}(F) as 2 x 2 matrices of fields
sym2(a, b, d) = [a b; b d]
mcases = (
    (:lp, "c = diag(f1, f2)   decoupled Vlasov species (control: must vanish)",
        F -> sym2(F[1], 0 * F[1], F[2]), true, true),
    (:dec, "c = diag(f1^2, f2^2 + f2)   decoupled, nonlinear",
        F -> sym2(F[1] * F[1], 0 * F[1], F[2] * F[2] + F[2]), true, true),
    (:lin,
        "c = [f1, f1 + 1/2; f1 + 1/2, f2]   cross term g = f1 + k: (A) holds, algebra R+R",
        F -> sym2(F[1], F[1] + 1 // 2, F[2]), true, true),
    (:nil, "c = [f1, f2; f2, 0]   Lie-Poisson of R[e]/e^2 (associative, non-semisimple)",
        F -> sym2(F[1], F[2], 0 * F[1]), true, true),
    (:slaved,
        "c = [f1, f1^2/2; f1^2/2, f1^3/3]   slaved family int p^i p^j du1, p = (1, u1)",
        F -> sym2(F[1], F[1] * F[1] * (1 // 2), F[1] * F[1] * F[1] * (1 // 3)), true, true),
    (:frozen_ok,
        "c = diag((f1 + f2)^2, 1)   frozen f2, a = Phi(u1 + phi(u2)): (A) and (B)",
        F -> sym2((F[1] + F[2]) * (F[1] + F[2]), 0 * F[1], 1 + 0 * F[1]), true, true),
    (:frozen_bad,
        "c = diag(f1 f2, 1)   frozen f2, a = u1 u2: (A) holds, (B) FAILS -> 2-D yes, 4-D NO",
        F -> sym2(F[1] * F[2], 0 * F[1], 1 + 0 * F[1]), true, false),
    (:gprod,
        "c = [f1, f1 f2; f1 f2, f2]   cross term g = f1 f2: (A) FAILS (must fail in 2-D too)",
        F -> sym2(F[1], F[1] * F[2], F[2]), false, false),
    (:gsq, "c = [f1, f1^2/2; f1^2/2, f2]   cross term g = f1^2/2: (A) FAILS",
        F -> sym2(F[1], F[1] * F[1] * (1 // 2), F[2]), false, false),
    (:hess,
        "c = Hess(u1^2 u2 / 2) = [f2, f1; f1, 0]   d_k c^{ij} totally symmetric, (A) FAILS",
        F -> sym2(F[2], F[1], 0 * F[1]), false, false)
)
for (_, label, c, zero2, zero4) in mcases
    for (D, f, mf, expect_zero) in ((2, f2, mf2, zero2), (4, f4, mf4, zero4))
        t = @elapsed s, scale = mjacobiator(c, mf..., f)
        # `scale` is the largest of the three cyclic terms.  A vanishing verdict is only evidence
        # if it came from cancellation, so require it: without this, anything that zeroes the whole
        # pipeline reports every expect_zero case as a pass.
        ok = expect_zero ? (iszero(s) && !iszero(scale)) : !iszero(s)
        check("$(D)-D  " * label, ok,
            "Jacobiator " * fmt(s) * ", largest term " * fmt(scale) *
            @sprintf(", normalised %.3e  (%.1fs)", Float64(normalised(s, scale)), t))
    end
end

header("2c. Necessity of (A) made exact: c = Hess(u1^2 u2/2), linear functionals, Jacobiator == -J1  (4-D)")

# For A = <a f1>, B = <b f1>, C = <c f2> and f = (2, 3 + g) the three cyclic terms carry
# Theta^{112}_2 = 0, Theta^{121}_2 = 1, Theta^{211}_2 = 1, so the Jacobiator is J2 + J3 = -J1
# of section 1 -- the necessity computation of the theorem, term by term.
let z = 0 * la4, A = F -> [la4 + 0 * F[1], z + 0 * F[1]],
    B = F -> [lb4 + 0 * F[1], z + 0 * F[1]], C = F -> [z + 0 * F[1], lc4 + 0 * F[1]],
    f = [2 + 0 * lg4, 3 + lg4]

    s, _ = mjacobiator(F -> sym2(F[2], F[1], 0 * F[1]), A, B, C, f)
    J1 = mean(pb(la4, lb4) * pb(lc4, lg4))
    check("Hess, linear functionals: Jacobiator == -J1 != 0", s == -J1 && !iszero(J1),
        "Jacobiator " * fmt(s) * ", -J1 = " * fmt(-J1))
    s0, scale0 = mjacobiator(F -> sym2(F[1], F[2], 0 * F[1]), A, B, C, f)
    check("R[e]/e^2 (associative) on the same data: exactly zero (control)",
        iszero(s0) && !iszero(scale0), fmt(s0) * ", largest term " * fmt(scale0))
end

header("2b. Casimirs int kappa(f1,f2): iff sum_j dc^{ij} ^ d(d_j kappa) = 0 for every i  (4-D, exact)")

# {K, A} = < c^{ij} {K_i, A_j} > for K = int kappa(f); A the local functional above
function casimir_bracket(c, Kf, A, f)
    C, k, a = c(f), Kf(f), A(f)
    mean(sum(C[i, j] * pb(k[i], a[j]) for i in 1:2, j in 1:2))
end
A4 = mf4[1]
dec = F -> sym2(F[1], 0 * F[1], F[2])
for (label, c, Kf, expect_zero) in (
    ("decoupled species, kappa = f1^3 + f2^2 (separable): Casimir", dec,
    F -> [3 * F[1] * F[1], 2 * F[2]], true),
    ("decoupled species, kappa = f1 f2 (not separable): NOT a Casimir", dec,
    F -> [F[2], F[1]], false),
    ("slaved family p = (1, u1), kappa = f2 - f1^2/2 (the frozen combination): Casimir",
    F -> sym2(F[1], F[1] * F[1] * (1 // 2), F[1] * F[1] * F[1] * (1 // 3)),
    F -> [-F[1], 1 + 0 * F[1]], true),
    ("slaved family p = (1, u1), kappa = f1 f2: NOT a Casimir",
    F -> sym2(F[1], F[1] * F[1] * (1 // 2), F[1] * F[1] * F[1] * (1 // 3)),
    F -> [F[2], F[1]], false))
    s = casimir_bracket(c, Kf, A4, f4)
    check(label, expect_zero ? iszero(s) : !iszero(s), "{K, A} = " * fmt(s))
end

# ---------------------------------------------------------------------------
# 3. derivative-dependent structure functions, single field
# ---------------------------------------------------------------------------
# {A,B} = < c[f] {A_f, B_f} >, c[f] a local function of the jet of f.  Same jet machinery:
# {X, C} = DX[f].{C_f, c[f]} for any X, because cyclicity does not care what c is.

bracket(c, Af, Bf, F) = mean(c(F) * pb(Af(F), Bf(F)))
function jacobiator(c, Af, Bf, Cf, f)
    η(X) = pb(X(f), c(f))
    term(X, Y, Z) = bracket(c, X, Y, Jet(f, η(Z))).d
    t = (term(Af, Bf, Cf), term(Bf, Cf, Af), term(Cf, Af, Bf))
    return sum(t), maximum(abs, t)
end
"The general obstruction (*): sum_cyc < {A_f,B_f} Dc[f].{C_f, c[f]} >, Dc by a jet."
function obstruction_general(c, Af, Bf, Cf, f)
    A, B, C, q = Af(f), Bf(f), Cf(f), c(f)
    Dc(ψ) = c(Jet(f, ψ)).d
    mean(pb(A, B) * Dc(pb(C, q))) + mean(pb(B, C) * Dc(pb(A, q))) +
    mean(pb(C, A) * Dc(pb(B, q)))
end
functionals(φA, φB, φC, ψC) = (F -> φA * F, F -> φB * F * F, F -> φC * mean(φC * F) + ψC)
laplace(F) = sum(∂(∂(F, i), i) for i in 1:phasedim(F))
gradsq(F) = sum(∂(F, i) * ∂(F, i) for i in 1:phasedim(F))

# III.1's single-field data (2-D: (x, v); 4-D: (x1, x2, v1, v2)), 4-D with a closing mode
sf2 = 3 + cosk((1, 0)) + sink((1, 1), 1 // 2) + cosk((0, 1), 1 // 3)
sfun2 = functionals(1 + cosk((0, 1), 1 // 2) + sink((1, -1), 1 // 3),
    sink((1, 1), 1 // 2) + cosk((1, 0)), cosk((1, -1)) + sink((0, 1), 1 // 2),
    sink((1, 0)) + cosk((2, 1), 1 // 2))
g2 = 1 + cosk((1, 0), 1 // 2) + sink((0, 1), 1 // 3)                     # a fixed function g(z)
sf4 = 3 + cosk((1, 0, 0, 0)) + sink((0, 1, 1, 0), 1 // 2) + cosk((-1, 0, 0, 1), 1 // 3) +
      gen((1, 1, 1, 1), 1 // 4) + (1 // 4) * prodc
sfun4 = functionals(1 + cosk((0, 1, 0, 0), 1 // 2) + sink((0, 0, 1, 0), 1 // 3),
    sink((1, 0, 0, 1), 1 // 2) + cosk((0, 0, 1, 0)), cosk((1, -1, 0, 0)) +
                                                     sink((0, 0, 0, 1), 1 // 2),
    sink((1, 0, 0, 0)) + cosk((0, 0, 1, 1), 1 // 2))
g4 = 1 + cosk((1, 0, 0, -1), 1 // 2) + sink((0, 1, 0, 0), 1 // 3) +
     gen((0, 1, 1, 1), 1 // 4)

header("3. Derivative-dependent c: full Jacobiator, exact, 2-D and 4-D")

dcases = (
    ("c = f^2                (control: must vanish)", (F, g) -> F * F, true),
    ("c = f_x f_v            (eq:sym-sym-alternative, s = u^2/2; first order: must FAIL)",
        (F, g) -> ∂(F, 1) * ∂(F, phasedim(F) ÷ 2 + 1), false),
    ("c = |grad f|^2         (first order: must FAIL)", (F, g) -> gradsq(F), false),
    ("c = {f, g}, g fixed    (first order, symplectic invariant: must FAIL)",
        (F, g) -> pb(F, g), false),
    ("c = f + {f, g}         (Vlasov plus a first-order term: must FAIL)",
        (F, g) -> F + pb(F, g), false),
    ("c = Laplace f          (second order, linear in f: must FAIL)",
        (F, g) -> laplace(F), false),
    ("c = f Laplace f        (second order: must FAIL)", (F, g) -> F * laplace(F), false),
    ("c = f_{x1 x1}          (second order: must FAIL)", (F, g) -> ∂(∂(F, 1), 1), false)
)
djac = Dict{Tuple{Int, Int}, Tuple{Q, Q}}()
for (idx, (label, c, expect_zero)) in enumerate(dcases)
    for (D, f, fun, g) in ((2, sf2, sfun2, g2), (4, sf4, sfun4, g4))
        cg = F -> c(F, g)
        t = @elapsed s, scale = jacobiator(cg, fun..., f)
        djac[(idx, D)] = (s, scale)
        # As in section 2: a vanishing verdict counts only if the terms that cancelled were there.
        check("$(D)-D  " * label, expect_zero ? (iszero(s) && !iszero(scale)) : !iszero(s),
            "Jacobiator " * fmt(s) * ", largest term " * fmt(scale) *
            @sprintf(", normalised %.3e  (%.1fs)", Float64(normalised(s, scale)), t))
    end
end

header("3a. The general obstruction (*) sum_cyc <{A_f,B_f} Dc.{C_f,c}> equals the Jacobiator (nonlinear functionals)")

for (idx, (label, c, _)) in enumerate(dcases)
    for (D, f, fun, g) in ((2, sf2, sfun2, g2), (4, sf4, sfun4, g4))
        O = obstruction_general(F -> c(F, g), fun..., f)
        check("$(D)-D  (*) == Jacobiator for " * strip(split(label, "  ")[1]),
            O == djac[(idx, D)][1],
            fmt(O) * " == " * fmt(djac[(idx, D)][1]))
    end
end

header("3b. Minimal exact example, 2-D, linear functionals A_f = sin(2x + v), B_f = sin v, C_f = cos(x - v)")

# f = 3 + cos x cos v, c = f_x f_v: exactly one of the three cyclic terms survives and equals 1/8.
# A_f must carry a generic mode rather than a single low one: at A_f = sin x the orders of the
# operator R(a, b) cancel against each other and the Jacobiator is exactly zero for a reason that
# has nothing to do with the theorem.
let a = sink((2, 1)), b = sink((0, 1)), γ = cosk((1, -1)),
    u = 3 + cosk((1, 0)) * cosk((0, 1)),
    lin = (F -> a + 0 * F, F -> b + 0 * F, F -> γ + 0 * F)

    s, scale = jacobiator(F -> ∂(F, 1) * ∂(F, 2), lin..., u)
    check("c = f_x f_v: Jacobiator of the linear functionals == 1/8", s == 1 // 8,
        "Jacobiator " * fmt(s) * ", largest term " * fmt(scale))
    O = obstruction_general(F -> ∂(F, 1) * ∂(F, 2), lin..., u)
    check(
        "c = f_x f_v: equals sum_cyc <{a,b} (d_x{c,q} f_v + f_x d_v{c,q})>, q = f_x f_v", O ==
                                                                                          s,
        fmt(O))
    s2, scale2 = jacobiator(F -> F * F, lin..., u)
    check(
        "c = f^2 on the same data: exactly zero (control)", iszero(s2) && !iszero(scale2),
        fmt(s2) * ", largest term " * fmt(scale2))
end

header("3c. Casimirs int kappa(f): lost under any derivative dependence  ({K, A}, K = int f^3/3, exact, 2-D)")

for (label, c, expect_zero) in (
    ("c = f^2 (control): {K, A} = 0", F -> F * F, true),
    ("c = f_x f_v: {K, A} != 0", F -> ∂(F, 1) * ∂(F, 2), false),
    ("c = |grad f|^2: {K, A} != 0", F -> gradsq(F), false),
    ("c = Laplace f: {K, A} != 0", F -> laplace(F), false),
    ("c = {f, g}: {K, A} != 0", F -> pb(F, g2), false))
    s = bracket(c, F -> F * F, sfun2[1], sf2)
    check(label, expect_zero ? iszero(s) : !iszero(s), "{K, A} = " * fmt(s))
end
# Mass survives every c in this family, and structurally: K_f = 1 and {1, A_f} = 0 whatever c is.
# On its own that cannot fail, so it is paired here with int f^3/3 on the SAME c, which must not
# vanish -- the pair is what carries the section's point, that derivative dependence costs the
# f-dependent Casimirs and leaves the mass.
let c = F -> ∂(F, 1) * ∂(F, 2)
    mass = bracket(c, F -> 1 + 0 * F, sfun2[1], sf2)
    cubic = bracket(c, F -> F * F, sfun2[1], sf2)
    check("c = f_x f_v: mass int f survives where int f^3/3 does not",
        iszero(mass) && !iszero(cubic),
        "{mass, A} = " * fmt(mass) * ", {int f^3/3, A} = " * fmt(cubic))
end

# ---------------------------------------------------------------------------
# 4. the four-bracket side: < T^{ijkl} (A_i B_j {C_k, D_l} - C_k D_l {A_i, B_j}) >
# ---------------------------------------------------------------------------

header("4. Multi-field 4-bracket with constant T: {A,S;B,S} == < T^{ijkl} s_j s_l {A_i, B_k} >  (2-D, exact)")

function fourbracket(T, A, B, C, D)
    m = length(A)
    mean(sum(T[i, j, k, l] * (A[i] * B[j] * pb(C[k], D[l]) - C[k] * D[l] * pb(A[i], B[j]))
    for i in 1:m, j in 1:m, k in 1:m, l in 1:m))
end
Tdiag = zeros(Int, 2, 2, 2, 2)
Tdiag[1, 1, 1, 1] = 1;
Tdiag[2, 2, 2, 2] = 1;                    # delta^{ij} delta^{kl}: one entropy per species
Tnil = zeros(Int, 2, 2, 2, 2)
Tnil[1, 1, 1, 1] = 1;
Tnil[1, 2, 2, 2] = 1;
Tnil[2, 2, 1, 2] = 1;                     # c^{11} = s_1^2, c^{12} = s_2^2, c^{22} = 0: Lie-Poisson of R[e]/e^2 for s = 2/3 sum u^{3/2}
for (label, T, S, cexp, poisson) in (
    ("T diagonal, s = u1^2 + u2^3 (separable): c = diag(4 f1^2, 9 f2^4)", Tdiag,
    F -> [2 * F[1], 3 * F[2] * F[2]], F -> sym2(4 * F[1] * F[1], 0 * F[1], 9 * F[2] * F[2] *
                                                                           F[2] * F[2]), true),
    ("T diagonal, s = u1^2 + u1 u2 (not separable): c = diag((2f1 + f2)^2, f1^2)", Tdiag,
    F -> [2 * F[1] + F[2], F[1]], F -> sym2(
        (2 * F[1] + F[2]) * (2 * F[1] + F[2]), 0 *
                                               F[1], F[1] *
                                                     F[1]), false),
    ("T of R[e]/e^2, s = u1^2 + u2^2: c = [4 f1^2, 4 f2^2; 4 f2^2, 0]", Tnil,
    F -> [2 * F[1], 2 * F[2]], F -> sym2(4 * F[1] * F[1], 4 * F[2] * F[2], 0 * F[1]), false))
    A, B = mf2[1](f2), mf2[2](f2)
    Sf = S(f2)
    lhs = fourbracket(T, A, Sf, B, Sf)
    rhs = mbracket(cexp, mf2[1], mf2[2], f2)
    check(label * ": {A,S;B,S} == {A,B}_c", lhs == rhs, fmt(lhs) * " == " * fmt(rhs))
    # T^{ijkl} = T^{klij} for both arrays here, and that symmetry is what forces the antisymmetry;
    # assert it alongside, so the check is about the four-bracket and not only about two literals.
    check(label * ": antisymmetric in A, B (forced by T^{ijkl} = T^{klij})",
        lhs == -fourbracket(T, B, Sf, A, Sf) &&
            all(T[i, j, k, l] == T[k, l, i, j] for i in 1:2, j in 1:2, k in 1:2, l in 1:2),
        "sum with the swap " * fmt(lhs + fourbracket(T, B, Sf, A, Sf)))
    # Checked here rather than recomputed after the loop: the same three structure functions were
    # entered twice and every Jacobiator evaluated twice, with the printed one gating nothing.
    s, scale = mjacobiator(cexp, mf2..., f2)
    check(
        label * (poisson ? ": the generated bracket is Poisson" :
         ": the generated bracket FAILS Jacobi"),
        poisson ? (iszero(s) && !iszero(scale)) : !iszero(s),
        "Jacobiator " * fmt(s) * ", largest term " * fmt(scale))
end

summary("III4_multifield_derivative.jl")
