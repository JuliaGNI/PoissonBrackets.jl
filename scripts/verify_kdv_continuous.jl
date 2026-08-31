#!/usr/bin/env julia
#
# The continuous bi-Hamiltonian structure of Section 1 of the KdV notes.
#
#     julia --project=scripts scripts/verify_kdv_continuous.jl
#
# Differential polynomials in u are represented in the jet variables u0 = u, u1 = u_x,
# u2 = u_xx, ...  The total derivative and the Euler operator
#
#     D f  = sum_k u_{k+1} df/du_k ,
#     E(f) = sum_k (-D)^k df/du_k ,
#
# give an EXACT test for the two things the section asserts: E is the variational derivative
# delta/delta u, and f integrates to zero over a periodic domain iff E(f) = 0 in every
# dependent variable, i.e. iff f is a total x-derivative.
#
# The convention is the textbook one, u_t + 6uu_x + u_xxx = 0, so that the sech^2 solitons are
# elevations. It is reached from the older u_t = 6uu_x - u_xxx by u -> -u, under which
# H2 = 1/2 int u^2 is unchanged, the cubic term of H1 changes sign, and the SECOND structure
# picks up a sign on its u-dependent part but NOT on its third derivative:
# 4u d_x + 2u_x - d_x^3 becomes -4u d_x - 2u_x - d_x^3. That asymmetry is the one thing here
# that cannot be guessed, and it is checked below rather than asserted.
#
#   1. delta H1/delta u = -3u^2 - u_xx  and  delta H2/delta u = u;
#   2. both structures generate the same equation, u_t = -6uu_x - u_xxx;
#   3. D2 = -4u d_x - 2u_x - d_x^3 is skew-adjoint;
#   4. {H2,H1}_1 = 0 and {H1,H2}_2 = 0, with the closed forms the text displays;
#   5. the intermediate step of {H2,H1}_1 is one integration by parts from the exact
#      pointwise integrand, and the difference is exactly -d_x ( u u_xx ). An early draft
#      displayed a form that was a spurious d_x ( 3u^3 ) away from it; the check is kept
#      because that slip was invisible to the conclusion, every one of these expressions
#      being a total derivative.
#
# Needs SymPy. The jet indices below are the mathematical ones -- u[k] means u_k, the k-th
# x-derivative -- so every array access carries a +1 for Julia's 1-based indexing.

using SymPyPythonCall

include(joinpath(@__DIR__, "check.jl"));
using .Checks: header, check, summary

const NJ = 8
const Uj = [Sym("u$i") for i in 0:(NJ - 1)]
const Aj = [Sym("a$i") for i in 0:(NJ - 1)]
const Bj = [Sym("b$i") for i in 0:(NJ - 1)]
const FAMILIES = (Uj, Aj, Bj)

"Total derivative with respect to x, applied `n` times."
function D(f, n = 1)
    for _ in 1:n
        f = expand(sum(fam[k + 2] * diff(f, fam[k + 1])
        for fam in FAMILIES for k in 0:(NJ - 2)))
    end
    f
end

"Euler operator of `f` with respect to the family `fam`."
euler(f, fam) = expand(sum((-1)^k * D(diff(f, fam[k + 1]), k) for k in 0:(NJ - 1)))

"`f` is a total x-derivative iff its Euler operator vanishes identically."
is_total_derivative(f) = all(simplify(euler(f, fam)) == 0 for fam in FAMILIES)

u(k) = Uj[k + 1]
a(k) = Aj[k + 1]
b(k) = Bj[k + 1]

header("0. the Euler-operator machinery")
check("D and E are consistent: E(D f) = 0 for a sample f",
    is_total_derivative(D(u(0)^2 * u(1) + u(2)^3)))
check("a non-total-derivative is detected: E(u^2) != 0", !is_total_derivative(u(0)^2))

header("1. functional derivatives")
h1 = (u(1)^2 - 2 * u(0)^3) / 2
h2 = u(0)^2 / 2
dH1, dH2 = euler(h1, Uj), euler(h2, Uj)
check("delta H1 / delta u = -3u^2 - u_xx", simplify(dH1 - (-3 * u(0)^2 - u(2))) == 0,
    "got $dH1")
check("delta H2 / delta u = u", simplify(dH2 - u(0)) == 0, "got $dH2")

header("2. both structures generate the same KdV equation")
kdv = -6 * u(0) * u(1) - u(3)
D2(f) = expand(-4 * u(0) * D(f) - 2 * u(1) * f - D(f, 3))
check("bracket 1: d_x ( delta H1 / delta u ) = -6uu_x - u_xxx", simplify(D(dH1) - kdv) == 0)
check("bracket 2: D2 ( delta H2 / delta u ) = -6uu_x - u_xxx", simplify(D2(dH2) - kdv) == 0)

header("3. skew-adjointness of D2 = -4u d_x - 2u_x - d_x^3")
D2b = expand(-4 * u(0) * D(b(0)) - 2 * u(1) * b(0) - D(b(0), 3))
D2a = expand(-4 * u(0) * D(a(0)) - 2 * u(1) * a(0) - D(a(0), 3))
check("int ( a D2 b + b D2 a ) dx = 0 for arbitrary a, b",
    is_total_derivative(expand(a(0) * D2b + b(0) * D2a)))

header("4. conservation of the respective other Hamiltonian")
I1 = expand(dH2 * D(dH1))
I2 = expand(dH1 * D2(dH2))
check("{H2,H1}_1 = 0", is_total_derivative(I1))
check("{H1,H2}_2 = 0", is_total_derivative(I2))
check("the {H1,H2}_2 computation is exact pointwise, integrand = 1/2 d_x (3u^2 + u_xx)^2",
    simplify(I2 - D((3 * u(0)^2 + u(2))^2) / 2) == 0)

header("5. the intermediate step of {H2,H1}_1")
stale = D(u(1)^2 + 2 * u(0)^3) / 2
corrected = D(u(1)^2 - 4 * u(0)^3) / 2
exact = D(u(1)^2 / 2 - u(0) * u(2) - 2 * u(0)^3)
check(
    "the integrand is exactly d_x ( u_x^2/2 - u u_xx - 2u^3 ), with no integration by parts",
    simplify(I1 - exact) == 0)
check(
    "the displayed 1/2 d_x ( u_x^2 - 4u^3 ) is exactly one integration by parts away: " *
    "the difference is -d_x ( u u_xx )",
    simplify((I1 - corrected) + D(u(0) * u(2))) == 0)
check(
    "a form off by a spurious d_x ( 3u^3 ) would NOT be: this is the shape of the slip " *
    "an early draft made, in the other sign convention",
    simplify((corrected - stale) + D(3 * u(0)^3)) == 0)
check("and it is not equal to the preceding line", simplify(I1 - stale) != 0,
    "difference = $(simplify(I1 - stale))")
# the conclusion of the section survives the slip: every one of these expressions is a total
# derivative, so all of them integrate to zero.
check(
    "both forms still integrate to zero, which is why the slip was invisible to the " *
    "conclusion {H2,H1}_1 = 0",
    is_total_derivative(stale) && is_total_derivative(corrected))

summary("verify_kdv_continuous.jl")
