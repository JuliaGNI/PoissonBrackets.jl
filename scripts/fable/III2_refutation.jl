#
# Adversarial re-check of III.2, on gl(N)*/u(N)* with the trace pairing.
#
#   negative half:  {F,G}_c(W) = tr(c(W) [F_W, G_W])  is Poisson iff c is affine, for N >= 3
#   positive half:  thm:kernel -- an Ad-invariant spectral-kernel bracket
#                   {F,G}_theta(V) = sum_{pq} theta(v_p,v_q) Xt_pq Yt_qp
#                   is Poisson iff rho_pq = (v_p-v_q)/varpi_pq is a coboundary, iff it is a
#                   pushforward of Lie-Poisson under a spectral psi with varpi = 1/psi^[1].
#
# An adversarial re-check, independent of III2_matrix_structure_function.jl: it includes nothing
# and imports nothing beyond LinearAlgebra and Test.  What is genuinely independent is
# `bracket_gradient`, the perturbation-theory route below -- that is where the weight sits.  Three
# closed forms are re-derived here rather than borrowed and are deliberately the same mathematics
# as their counterparts there: `psi_prediction` against `jacobiator_kernel`, `dpoly_mat` against
# `fderiv(::Poly)`, `closed_form_jacobiator` against `jacobiator_closed`.  They are predictions to
# be falsified against `bracket_gradient`, not corroborations of it.
#
# Scope: the matrix c-bracket is evaluated at a DIAGONAL W only, and sections 5, 6b and 7d of
# III2_matrix_structure_function.jl -- the semiclassical N -> infinity claims -- are not re-checked.
#
# Everything is exact over Q, so a claimed vanishing is a vanishing.
#
# THE POINT OF THIS SCRIPT.  A spectral-kernel bracket is defined through the eigendecomposition of
# V, so its Jacobiator needs the derivative of the EIGENVECTORS, not only of the eigenvalues.  That
# is the one step thm:kernel's proof rests on and the one step a closed-form cross-check cannot
# reach.  Here the gradient of V |-> {F,G}_theta(V) is instead computed by brute force from
# first-order perturbation theory at a diagonal V with distinct rational eigenvalues,
#
#     v_p -> v_p + eps E_pp ,      U -> I + eps Gamma ,   Gamma_pq = E_pq/(v_p - v_q),  Gamma_pp = 0,
#
# which is rational in the entries of E.  No closed form is assumed anywhere; Psi is a prediction
# that is compared against the result, on kernels that are neither c-brackets nor pushforwards.
#
# Run:  julia --project=scripts scripts/fable/III2_refutation.jl
#

using LinearAlgebra
using Test

const Q = Rational{BigInt}

# ---------------------------------------------------------------------------
# the spectral-kernel bracket and its Jacobiator, brute force
# ---------------------------------------------------------------------------
# A kernel is given by theta(x,y) together with its two partial derivatives, all as exact rational
# functions of the eigenvalues.  V is diagonal with entries `v`; the "matrices" are plain
# Matrix{Q} and the pairing is <A,B> = tr(A*B).

"""
Gradient of the scalar function V |-> {F,G}_theta(V) at the diagonal V = diag(v), for linear F, G
with constant gradients X, Y.  Returned as the matrix M with <M,E> = D{F,G}(V).E for every E.

Only the product rule and first-order eigenvalue/eigenvector perturbation theory are used.
"""
function bracket_gradient(v, θ, ∂1θ, ∂2θ, X, Y)
    N = length(v)
    T = promote_type(eltype(X), eltype(Y), typeof(θ(v[1], v[2])))
    M = zeros(T, N, N)
    for i in 1:N, j in 1:N

        E = zeros(T, N, N)
        E[i, j] = one(T)
        # eigenvalue motion: v_p -> v_p + eps E_pp
        # p == q contributes nothing for an antisymmetric kernel: theta(x,x) = 0 and
        # d1theta + d2theta = 0 there.  It is skipped rather than evaluated because a kernel like
        # 2(a-b)^2/log(a/b) is 0/0 on the diagonal and would poison the sum with a NaN.
        d = zero(T)
        for p in 1:N, q in 1:N

            p == q && continue
            d += (∂1θ(v[p], v[q]) * E[p, p] + ∂2θ(v[p], v[q]) * E[q, q]) * X[p, q] * Y[q, p]
        end
        # eigenbasis rotation: Xt -> X + eps [X, Gamma]
        # U^{-1}(V + eps E)U diagonal to first order needs E_pq + (v_p - v_q) Gamma_pq = 0,
        # i.e. Gamma_pq = E_pq/(v_q - v_p).  The sign is easy to get the wrong way round and the
        # answer depends on it.  The formula and its sign are checked on an explicit 2x2 in the
        # first testset; that this function uses them is pinned downstream, by route A agreeing
        # with route B, which never mentions an eigenvector at all.
        Γ = zeros(T, N, N)
        for p in 1:N, q in 1:N

            p == q && continue
            Γ[p, q] = E[p, q] / (v[q] - v[p])
        end
        XΓ = X * Γ - Γ * X
        YΓ = Y * Γ - Γ * Y
        for p in 1:N, q in 1:N

            p == q && continue
            d += θ(v[p], v[q]) * (XΓ[p, q] * Y[q, p] + X[p, q] * YΓ[q, p])
        end
        # <M,E> = tr(M E) = M[j,i] for E = e_i e_j^T
        M[j, i] = d
    end
    M
end

"""
The Jacobiator of the spectral-kernel bracket at the diagonal V = diag(v), brute force.

{{F,G},H} = sum_pq theta_pq M_pq Z_qp with M = grad{F,G}: at a diagonal V the eigenbasis is the
identity, so no further rotation enters this outer evaluation.
"""
function kernel_jacobiator(v, θ, ∂1θ, ∂2θ, X, Y, Z)
    N = length(v)
    total = zero(promote_type(eltype(X), eltype(Y), eltype(Z)))
    for (A, B, C) in ((X, Y, Z), (Y, Z, X), (Z, X, Y))
        M = bracket_gradient(v, θ, ∂1θ, ∂2θ, A, B)
        for p in 1:N, q in 1:N

            p == q && continue
            total += θ(v[p], v[q]) * M[p, q] * C[q, p]
        end
    end
    total
end

"thm:kernel's predicted Jacobiator, as a prediction to be tested -- not used to compute anything."
function psi_prediction(v, θ, X, Y, Z)
    N = length(v)
    ϖ(p, q) = θ(v[p], v[q]) / (v[p] - v[q])
    total = zero(promote_type(eltype(X), eltype(Y), eltype(Z)))
    for p in 1:N, q in 1:N, r in 1:N
        (p == q || q == r || r == p) && continue
        Ψ = ϖ(p, q) * ϖ(q, r) * (v[p] - v[r]) +
            ϖ(q, r) * ϖ(r, p) * (v[q] - v[p]) +
            ϖ(r, p) * ϖ(p, q) * (v[r] - v[q])
        total += Ψ * X[p, q] * (Y[q, r] * Z[r, p] - Z[q, r] * Y[r, p])
    end
    total
end

# ---------------------------------------------------------------------------
# the literal matrix c-bracket, by a route that never enters the eigenbasis
# ---------------------------------------------------------------------------
# For polynomial c the Frechet derivative is Dc_W(E) = sum_k a_k sum_{i+j=k-1} W^i E W^j, which
# needs neither divided differences nor eigenvectors.  This is the independent route against which
# the perturbation-theory machinery above is checked.

polyval_mat(a::Vector{Q}, W::Matrix{Q}) = sum(a[k + 1] * W^k for k in 0:(length(a) - 1))

function dpoly_mat(a::Vector{Q}, W::Matrix{Q}, E::Matrix{Q})
    N = size(W, 1)
    out = zeros(Q, N, N)
    for k in 1:(length(a) - 1)
        iszero(a[k + 1]) && continue
        for i in 0:(k - 1)
            out += a[k + 1] * (W^i) * E * (W^(k - 1 - i))
        end
    end
    out
end

comm(A, B) = A * B - B * A

"Gradient of V |-> tr(c(V)[X,Y]) at W, from the polynomial Frechet derivative alone."
function cbracket_gradient(a::Vector{Q}, W::Matrix{Q}, X::Matrix{Q}, Y::Matrix{Q})
    N = size(W, 1)
    M = zeros(Q, N, N)
    XY = comm(X, Y)
    for i in 1:N, j in 1:N

        E = zeros(Q, N, N)
        E[i, j] = one(Q)
        M[j, i] = tr(dpoly_mat(a, W, E) * XY)
    end
    M
end

"Jacobiator of the literal matrix c-bracket, brute force, no eigenbasis anywhere."
function cbracket_jacobiator(a::Vector{Q}, W::Matrix{Q}, X, Y, Z)
    cW = polyval_mat(a, W)
    total = zero(Q)
    for (A, B, C) in ((X, Y, Z), (Y, Z, X), (Z, X, Y))
        M = cbracket_gradient(a, W, A, B)
        total += tr(cW * comm(M, C))
    end
    total
end

"The refutation note's closed form (star): sum_cyc tr([X,Y] Dc_W^2([Z,W])), at diagonal W."
function closed_form_jacobiator(a::Vector{Q}, v::Vector{Q}, X, Y, Z)
    N = length(v)
    cv = [sum(a[k + 1] * v[p]^k for k in 0:(length(a) - 1)) for p in 1:N]
    dd(p, q) = p == q ?
               sum(k * a[k + 1] * v[p]^(k - 1) for k in 1:(length(a) - 1); init = zero(Q)) :
               (cv[p] - cv[q]) // (v[p] - v[q])
    W = Matrix{Q}(Diagonal(v))
    total = zero(Q)
    for (A, B, C) in ((X, Y, Z), (Y, Z, X), (Z, X, Y))
        ZW = comm(C, W)
        D2 = [dd(p, q)^2 * ZW[p, q] for p in 1:N, q in 1:N]
        total += tr(comm(A, B) * D2)
    end
    total
end

# ---------------------------------------------------------------------------
# test data
# ---------------------------------------------------------------------------

"A reproducible rational matrix with no symmetry, entries small enough to read."
function testmat(N, seed)
    pool = Q[
        2, -3, 5, 1, -1, 4, -2, 7, 3, -5, 6, -4, 1, 2, -7, 5, 3, -1, 8, -3, 2, 4, -6, 1, 5]
    [pool[mod1(seed * 7 + 3 * i + 5 * j, length(pool))] // (1 + mod(i + j, 3))
     for i in 1:N, j in 1:N]
end

# ---------------------------------------------------------------------------
# exact partial derivatives of a rational kernel, by dual numbers over Q
# ---------------------------------------------------------------------------
# Hand-differentiating each kernel is the obvious place to introduce an error that then agrees
# with itself, so the partials are taken from the kernel expression itself.

struct Dual{T <: Number} <: Number
    a::T
    b::T
end
Dual(a::T) where {T} = Dual{T}(a, zero(a))
Base.convert(::Type{Dual{T}}, x::Number) where {T} = Dual{T}(convert(T, x), zero(T))
# without this, the rule above matches a Dual and tries to convert it to T
Base.convert(::Type{Dual{T}}, x::Dual) where {T} = Dual{T}(convert(T, x.a), convert(T, x.b))
function Base.promote_rule(::Type{Dual{T}}, ::Type{S}) where {T, S <: Number}
    Dual{promote_type(T, S)}
end
# without this, the rule above matches Dual against itself and builds Dual{Dual{...}}
Base.promote_rule(::Type{Dual{T}}, ::Type{Dual{S}}) where {T, S} = Dual{promote_type(T, S)}
Base.:+(x::Dual, y::Dual) = Dual(x.a + y.a, x.b + y.b)
Base.:-(x::Dual, y::Dual) = Dual(x.a - y.a, x.b - y.b)
Base.:-(x::Dual) = Dual(-x.a, -x.b)
Base.:*(x::Dual, y::Dual) = Dual(x.a * y.a, x.a * y.b + x.b * y.a)
Base.:/(x::Dual, y::Dual) = Dual(x.a / y.a, (x.b * y.a - x.a * y.b) / (y.a * y.a))
Base.zero(::Type{Dual{T}}) where {T} = Dual(zero(T), zero(T))
Base.one(::Type{Dual{T}}) where {T} = Dual(one(T), zero(T))
Base.log(x::Dual) = Dual(log(x.a), x.b / x.a)

∂1(θ) = (x, y) -> θ(Dual(x, one(x)), Dual(y, zero(y))).b
∂2(θ) = (x, y) -> θ(Dual(x, zero(x)), Dual(y, one(y))).b

"rho_pq + rho_qr + rho_rp, the cocycle defect thm:kernel's criterion is stated in."
function cocycle_defect(v, θ, p, q, r)
    ρ(i, j) = (v[i] - v[j])^2 / θ(v[i], v[j])
    ρ(p, q) + ρ(q, r) + ρ(r, p)
end

println("=" ^ 78)
println("III.2 adversarial re-check -- exact arithmetic over Q")
println("=" ^ 78)

@testset "III.2 refutation attempt" begin

    # -----------------------------------------------------------------------
    @testset "machinery: the eigenvector-perturbation formula and its sign" begin
        # The formula Gamma_pq = E_pq/(v_q - v_p) and its sign, against an explicitly diagonalised
        # 2x2.  That `bracket_gradient` uses this formula is not tested here -- it is pinned by the
        # route A / route B agreement in the next testset, route B having no eigenvectors in it.
        v2 = Q[1, 2]
        E = Q[0 1; 0 0]
        Γ = Q[0 (E[1, 2]//(v2[2] - v2[1])); (E[2, 1]//(v2[1] - v2[2])) 0]
        @test Γ == Q[0 1; 0 0]           # matches U = [[1,eps],[0,1]] computed by hand
        # and the same by construction: (I - eps G)(V + eps E)(I + eps G) is diagonal to O(eps)
        Vd = Matrix{Q}(Diagonal(v2))
        @test all(iszero, (E + comm(Vd, Γ)) - Matrix{Q}(Diagonal(diag(E + comm(Vd, Γ)))))
        println("  eigenvector perturbation: Gamma_pq = E_pq/(v_q - v_p), verified on a 2x2")
    end

    # -----------------------------------------------------------------------
    @testset "machinery: the two routes to the c-bracket Jacobiator agree" begin
        # Route A: spectral kernel theta = c(x) - c(y), gradient by eigenvalue AND eigenvector
        #          perturbation.
        # Route B: the literal matrix bracket, gradient from the polynomial Frechet derivative.
        # Route C: the closed form (star) of the refutation note.
        # A shares no line with B; B never mentions an eigenvector.
        for (N, a) in [
            (3, Q[0, 0, 1]), (3, Q[1, -2, 0, 3]), (4, Q[0, 0, 1]), (4, Q[2, 1, -1, 1])]
            v = Q[1, 2, 4, 7][1:N]
            W = Matrix{Q}(Diagonal(v))
            X, Y, Z = testmat(N, 1), testmat(N, 2), testmat(N, 3)
            cv(x) = sum(a[k + 1] * x^k for k in 0:(length(a) - 1))
            dcv(x) = sum(k * a[k + 1] * x^(k - 1) for k in 1:(length(a) - 1); init = zero(Q))
            θ(x, y) = cv(x) - cv(y)
            jA = kernel_jacobiator(v, θ, (x, y) -> dcv(x), (x, y) -> -dcv(y), X, Y, Z)
            jB = cbracket_jacobiator(a, W, X, Y, Z)
            jC = closed_form_jacobiator(a, v, X, Y, Z)
            @test jA == jB
            @test jB == jC
            # Three routes agreeing on zero would agree for nothing.  Every `a` here is non-affine
            # and every N is at least 3, so the Jacobiator must not vanish.
            @test !iszero(jA)
            println("  N = $N, c = $(a):  kernel route = matrix route = closed form = ", jA)
        end
    end

    # -----------------------------------------------------------------------
    @testset "negative half: affine only, and N = 2 is not a control" begin
        cases = [("c = v (Lie-Poisson)", Q[0, 1]), ("c = 3v - 2 (affine)", Q[-2, 3]),
            ("c = v^2", Q[0, 0, 1]), ("c = v^3", Q[0, 0, 0, 1]),
            ("c = 1 + v - v^2 + v^3/2", Q[1, 1, -1, 1 // 2])]
        for (name, a) in cases
            affine = length(a) <= 2
            row = String[]
            for N in 2:5
                v = Q[1, 2, 4, 7, 11][1:N]
                W = Matrix{Q}(Diagonal(v))
                X, Y, Z = testmat(N, 1), testmat(N, 2), testmat(N, 3)
                j = cbracket_jacobiator(a, W, X, Y, Z)
                push!(row, "N=$N: $j")
                if affine || N == 2
                    @test j == 0
                else
                    @test j != 0        # the control must be able to fail, and does
                end
            end
            println("  ", rpad(name, 26), join(row, "   "))
        end
        println("  -> N = 2 vanishes for every c: no triple of distinct eigenvalues exists, so")
        println("     Phi_pqr never appears.  N = 2 is a degeneracy, not evidence.")
    end

    # -----------------------------------------------------------------------
    @testset "negative half: the conclusion survives on u(N) itself" begin
        # The tests above run over gl(N,Q).  u(N) is a real form of gl(N,C) and the Jacobiator is
        # complex-trilinear, so non-vanishing on gl transfers; here it is confirmed directly on
        # anti-Hermitian data over Q(i), with V Hermitian (the standard u(N)* identification).
        N = 3
        v = Q[1, 2, 4]
        function antiherm(seed)
            A = Complex{Q}.(testmat(N, seed)) + im * Complex{Q}.(testmat(N, seed + 10))
            A - adjoint(A)
        end
        X, Y, Z = antiherm(1), antiherm(2), antiherm(3)
        for (name, a) in [
            ("c = v", Q[0, 1]), ("c = v^2", Q[0, 0, 1]), ("c = v^3", Q[0, 0, 0, 1])]
            cv(x) = sum(a[k + 1] * x^k for k in 0:(length(a) - 1))
            dcv(x) = sum(k * a[k + 1] * x^(k - 1) for k in 1:(length(a) - 1); init = zero(Q))
            j = kernel_jacobiator(Complex{Q}.(v), (x, y) -> cv(x) - cv(y),
                (x, y) -> dcv(x), (x, y) -> -dcv(y), X, Y, Z)
            length(a) <= 2 ? (@test j == 0) : (@test j != 0)
            println("  u(3), anti-Hermitian X,Y,Z, ", rpad(name, 8), " Jacobiator = ", j)
        end
    end

    # -----------------------------------------------------------------------
    @testset "thm:kernel: Psi against brute force, off the two special families" begin
        # Psi is verified above only where varpi is a divided difference, because there Psi and the
        # note's Phi are algebraically equal.  The theorem's content is the general kernel, whose
        # Jacobiator needs the eigenvector rotation.  These kernels are chosen so that some are
        # pushforwards (varpi = 1/psi^[1], predict zero) and some are neither pushforwards nor
        # c-brackets (predict nonzero, and equal to Psi).
        kernels = [
            ("Lie-Poisson         psi = v", (x, y) -> x - y, :zero),
            ("pushforward         psi = v^2", (x, y) -> (x - y) / (x + y), :zero),
            ("pushforward         psi = v^3",
                (x, y) -> (x - y) / (x^2 + x * y + y^2), :zero),
            ("pushforward         psi = -1/v", (x, y) -> x * y * (x - y), :zero),
            ("pushforward         psi = -1/(1+v)",
                (x, y) -> (1 + x) * (1 + y) * (x - y), :zero),
            ("c-bracket           c = v^2", (x, y) -> x^2 - y^2, :nonzero),
            ("c-bracket           c = v^3", (x, y) -> x^3 - y^3, :nonzero),
            ("neither             varpi = x^2+y^2",
                (x, y) -> (x^2 + y^2) * (x - y), :nonzero),
            ("neither             varpi = 1+x^2y^2",
                (x, y) -> (1 + x^2 * y^2) * (x - y), :nonzero),
            ("neither             varpi = 1/(1+x^2+y^2)",
                (x, y) -> (x - y) / (1 + x^2 + y^2), :nonzero)
        ]
        for N in (3, 4)
            v = Q[1, 2, 4, 7][1:N]
            X, Y, Z = testmat(N, 4), testmat(N, 5), testmat(N, 6)
            for (name, θ, expect) in kernels
                j = kernel_jacobiator(v, θ, ∂1(θ), ∂2(θ), X, Y, Z)
                ψpred = psi_prediction(v, θ, X, Y, Z)
                @test j == ψpred                       # the theorem's closed form, tested
                expect === :zero ? (@test j == 0) : (@test j != 0)
                # and the criterion it is turned into
                defect = cocycle_defect(v, θ, 1, 2, 3)
                expect === :zero ? (@test defect == 0) : (@test defect != 0)
                N == 3 && println("  ", rpad(name, 36), "Jac = ", rpad(string(j), 22),
                    "Psi = ", rpad(string(ψpred), 22), "cocycle defect = ", defect)
            end
        end
        println("  -> brute force equals Psi on every kernel, including the three that are")
        println("     neither c-brackets nor pushforwards; the cocycle criterion tracks it.")
    end

    # -----------------------------------------------------------------------
    @testset "the log-mean bracket: Poisson, and its continuum limit" begin
        # (a) The log-mean kernel is not rational, so it is tested in floating point, on the same
        #     data and through the same code path as the exact c = v^2 case above, which gives 176.
        N = 3
        vf = [1.0, 2.0, 4.0]
        Xf, Yf, Zf = Float64.(testmat(N, 4)), Float64.(testmat(N, 5)),
        Float64.(testmat(N, 6))
        θlog(x, y) = 2 * (x - y)^2 / log(x / y)          # varpi = 2(a-b)/log(a/b), psi = log(v)/2
        θsq(x, y) = x^2 - y^2                             # the literal c-bracket, c = v^2
        jlog = kernel_jacobiator(vf, θlog, ∂1(θlog), ∂2(θlog), Xf, Yf, Zf)
        jsq = kernel_jacobiator(vf, θsq, ∂1(θsq), ∂2(θsq), Xf, Yf, Zf)
        @test abs(jsq - 176) < 1e-9                       # the float path reproduces the exact 176
        # Absolute, not relative to jsq: |jsq| is 176, so the old 1e-9*|jsq| admitted 1.8e-7 on a
        # quantity whose round-off is ~1e-13.  "Poisson" means this is zero.
        @test abs(jlog) < 1e-9                            # the log-mean bracket is Poisson
        println(
            "  N = 3:  c = v^2 Jacobiator = ", jsq, " (exact 176);  log-mean Jacobiator = ",
            jlog)

        # (b) Why the two share a continuum limit: their kernels agree on the diagonal and differ
        #     at second order in the eigenvalue gap.  With a = v + d, b = v - d,
        #         2 L(a,b) - (a+b) = -2 d^2 / (3v) + O(d^4).
        #     Under the Zeitlin correspondence d = O(hbar) = O(1/N), which is the N^{-2} rate.
        setprecision(BigFloat, 200) do
            v = BigFloat(3)
            for k in 3:6
                d = BigFloat(10)^(-k)
                a, b = v + d, v - d
                diff = 2 * (a - b) / log(a / b) - (a + b)
                predicted = -2 * d^2 / (3 * v)
                @test abs(diff - predicted) < 1e-3 * abs(predicted)
                k == 4 && println("  kernel gap at d = 1e-$k:  measured ", Float64(diff),
                    "   predicted -2d^2/3v = ", Float64(predicted))
            end
        end
        # and the diagonal values agree exactly, which is what fixes the limit to be c' = 2v
        setprecision(BigFloat, 200) do
            v = BigFloat(3)
            d = BigFloat(10)^(-20)
            @test abs(2 * ((v + d) - (v - d)) / log((v + d) / (v - d)) - 2v) < 1e-30
        end
        println("  -> both kernels have diagonal value c'(v) = 2v and differ at O(gap^2),")
        println("     so the two brackets share the continuum limit and separate at N^{-2}.")
    end

    # -----------------------------------------------------------------------
    @testset "the continuum step behind the flat coordinate" begin
        # The positive half rests on: the c-bracket is the pushforward of the continuum Lie-Poisson
        # bracket under the pointwise field change g = Psi(f) with Psi' = 1/c'.  Expanding
        # {A_f/Psi', B_f/Psi'} produces, besides u^2{A_f,B_f}, the terms
        #     u u' ( A_f {f,B_f} - B_f {f,A_f} ) ,
        # and the whole claim turns on those reducing to a multiple of {A_f,B_f} under the integral:
        #     < h(f) A_f {f,B_f} > = - < H(f) {A_f,B_f} > ,   H' = h.
        # Without that reduction c' = 1/Psi' does not come out.  Checked exactly on T^2 below.
        D = 2
        R = Rational{BigInt}
        CQ = Complex{R}
        mul(p, q) = (r = Dict{NTuple{2, Int}, CQ}();
            for (ka, va) in p, (kb, vb) in q
                k = ka .+ kb
                r[k] = get(r, k, zero(CQ)) + va * vb
            end;
            Dict(k => v for (k, v) in r if !iszero(v)))
        der(p, j) = Dict(k => CQ(0, k[j]) * v for (k, v) in p if !iszero(k[j] * v))
        add(p, q) = (r = copy(p);
            for (k, v) in q
                r[k] = get(r, k, zero(CQ)) + v
            end;
            Dict(k => v for (k, v) in r if !iszero(v)))
        scal(s, p) = Dict(k => CQ(s) * v for (k, v) in p if !iszero(CQ(s) * v))
        mean0(p) = real(get(p, (0, 0), zero(CQ)))
        pb(a, b) = add(mul(der(a, 1), der(b, 2)), scal(-1, mul(der(a, 2), der(b, 1))))
        cosk(k, amp) = Dict(k => CQ(R(amp) // 2), (.-k) => CQ(R(amp) // 2))
        sink(k, amp) = Dict(k => CQ(0, -R(amp) // 2), (.-k) => CQ(0, R(amp) // 2))
        one0 = Dict((0, 0) => one(CQ))
        polyv(coef, p) = (r = scal(coef[end], one0);
            for i in (length(coef) - 1):-1:1
                r = add(scal(coef[i], one0), mul(p, r))
            end;
            r)

        f = add(add(scal(3, one0), cosk((1, 0), 1)), add(sink((0, 1), 2), cosk((1, 1), 1)))
        A = add(cosk((0, 1), 1), sink((1, 0), 3))
        B = add(sink((1, 1), 2), cosk((1, -1), 1))
        for h in ([R(1)], [R(0), R(1)], [R(2), R(-1), R(3)], [R(0), R(0), R(0), R(1)])
            H = vcat(R(0), [h[i] // i for i in eachindex(h)])          # antiderivative
            lhs = mean0(mul(mul(polyv(h, f), A), pb(f, B)))
            rhs = -mean0(mul(polyv(H, f), pb(A, B)))
            @test lhs == rhs
            @test !iszero(lhs)          # each h in turn, not only the linear one below
        end
        println("  < h(f) A_f {f,B_f} > = - < H(f) {A_f,B_f} > exactly, h up to cubic, on T^2")
        println("  -> the reduction that yields c' = 1/Psi' holds; the flat-coordinate claim is not")
        println("     an approximation.")
    end

    # -----------------------------------------------------------------------
    @testset "the pushforward family needs a definite spectrum" begin
        # thm:kernel asks for c' > 0, and for c = v^2 that is v > 0.  So the flat coordinate
        # psi' = 1/c' = 1/(2v) blows up at v = 0 and the log-mean kernel is not real across it.
        # The theorem is consistent; what this bounds is the *application* claim, since a Zeitlin
        # truncation of 2-D Euler vorticity has eigenvalues of both signs.
        ϖlog(a, b) = 2 * (a - b) / log(a / b)
        # Across zero the kernel is not a real number at all -- log of a negative ratio.
        for (a, b) in ((1.0, -2.0), (3.0, -0.5), (-1.0, 4.0))
            @test_throws DomainError ϖlog(a, b)
        end
        # Within one sign it is finite, which is what makes this a restriction rather than an
        # outright failure.
        for (a, b) in ((1.0, 2.0), (0.5, 4.0), (-2.0, -3.0))
            @test isfinite(ϖlog(a, b))
        end
        # And the agreement with the c-bracket kernel a + b, which is what "same continuum limit"
        # rests on, degrades as the spectrum approaches zero: 2L -> 0 while a + b -> a.
        gap(a, b) = abs(ϖlog(a, b) - (a + b)) / (a + b)
        @test gap(1.0, 0.9) < 0.01             # comparable eigenvalues: the two kernels agree
        @test gap(1.0, 1.0e-6) > 0.8           # one eigenvalue near zero: they do not
        @test gap(1.0, 1.0e-12) > gap(1.0, 1.0e-6)
        println("  c = v^2 has c' > 0 only on v > 0: the log-mean kernel 2(a-b)/log(a/b) is real")
        println("  only when a and b share a sign, so the pushforward representative of the")
        println("  enstrophy bracket exists on definite spectra only.")
    end
end
