#
# Adversarial re-check of III.1: the structure-function bracket
#
#     {A,B}_c[f] = < w(z) c(z,f) Pi(dA_f, dB_f) >          < . > the mean over T^D
#
# and its Jacobi identity.  Written by a session whose brief was to overturn III.1, independently
# of III1_structure_function_4d.jl: it shares no code with that script and no code with
# fabletools.jl, so that a defect in the shared machinery cannot hide in both.
#
# Everything is exact over Q(i) with BigInt numerators, so a claimed vanishing is a vanishing.
#
# The Jacobiator is computed *brute force*.  The functional derivative of a functional is obtained
# by carrying, alongside every field, the linear differential operator that is its first variation,
# and then integrating by parts on the torus (< m d^alpha v > = (-1)^|alpha| < d^alpha m  v >).
# Only the product rule, the chain rule and integration by parts are used.  In particular the
# result is NOT computed from the closed-form Jacobiator of III.1 or of the refutation note; that
# closed form is a separate prediction, and the two are compared.
#
# Predicted master formula, derived independently in fable/III.1-refutation.md:
#
#     Jac = - < w c S(dA_f, dB_f, dC_f, D(w c_u)) >,
#     S(a,b,g,d) = Pi(a,b)Pi(g,d) + Pi(b,g)Pi(a,d) + Pi(g,a)Pi(b,d),
#
# with D the total differential on M of z |-> (w c_u)(z, f(z)).  The tests below check this against
# brute force in cases where it is nonzero, and only then use it where it vanishes.
#
# Run:  julia --project=. scripts/fable/III1_refutation.jl
#

using LinearAlgebra
using Test

# ---------------------------------------------------------------------------
# exact trigonometric polynomials on T^D over Q(i), BigInt numerators
# ---------------------------------------------------------------------------

const R = Rational{BigInt}
const C = Complex{R}

struct Trig{D}
    c::Dict{NTuple{D, Int}, C}
end

Trig{D}() where {D} = Trig{D}(Dict{NTuple{D, Int}, C}())

function prune(p::Trig{D}) where {D}
    Trig{D}(Dict(k => v for (k, v) in p.c if !iszero(v)))
end

function Base.:+(a::Trig{D}, b::Trig{D}) where {D}
    c = copy(a.c)
    for (k, v) in b.c
        c[k] = get(c, k, zero(C)) + v
    end
    prune(Trig{D}(c))
end
Base.:-(a::Trig{D}) where {D} = Trig{D}(Dict(k => -v for (k, v) in a.c))
Base.:-(a::Trig{D}, b::Trig{D}) where {D} = a + (-b)

function Base.:*(a::Trig{D}, b::Trig{D}) where {D}
    c = Dict{NTuple{D, Int}, C}()
    for (ka, va) in a.c, (kb, vb) in b.c

        k = ka .+ kb
        c[k] = get(c, k, zero(C)) + va * vb
    end
    prune(Trig{D}(c))
end

function scale(s, a::Trig{D}) where {D}
    prune(Trig{D}(Dict(k => C(s) * v for (k, v) in a.c)))
end

"Coordinate derivative d/dz_j: the k-th coefficient picks up i k_j."
function der(a::Trig{D}, j::Int) where {D}
    prune(Trig{D}(Dict(k => C(0, k[j]) * v for (k, v) in a.c)))
end

konst(D, s) = prune(Trig{D}(Dict(ntuple(_ -> 0, D) => C(s))))
zerotrig(D) = Trig{D}()

function harm(k::NTuple{D, Int}, amp, phase::Symbol) where {D}
    # amp*cos(k.z) or amp*sin(k.z), real by construction
    if phase === :cos
        prune(Trig{D}(Dict(k => C(R(amp) // 2), (.-k) => C(R(amp) // 2))))
    else
        prune(Trig{D}(Dict(k => C(0, -R(amp) // 2), (.-k) => C(0, R(amp) // 2))))
    end
end

"The mean over the torus.  Real by construction for the real fields used here."
function meanval(p::Trig{D}) where {D}
    v = get(p.c, ntuple(_ -> 0, D), zero(C))
    @assert iszero(imag(v)) "mean acquired an imaginary part: $v"
    real(v)
end

"Exact zero test, and its negation: < p^2 > by Parseval, zero iff p vanishes identically."
sqnorm(p::Trig) = sum(abs2, values(p.c); init = zero(R))
isz(p::Trig) = iszero(sqnorm(p))

# ---------------------------------------------------------------------------
# fields carrying their first variation as a linear differential operator
# ---------------------------------------------------------------------------
# Var(v, L) is the field v together with L : delta f |-> sum_alpha m_alpha d^alpha (delta f),
# stored as multi-index => coefficient field.

struct Var{D}
    v::Trig{D}
    L::Dict{NTuple{D, Int}, Trig{D}}
end

lift(p::Trig{D}) where {D} = Var{D}(p, Dict{NTuple{D, Int}, Trig{D}}())
"The field f itself: value f, first variation the identity operator."
seed(f::Trig{D}) where {D} = Var{D}(f, Dict(ntuple(_ -> 0, D) => konst(D, 1)))

function addop(a::Dict{NTuple{D, Int}, Trig{D}}, b::Dict{NTuple{D, Int}, Trig{D}}) where {D}
    c = copy(a)
    for (k, v) in b
        c[k] = haskey(c, k) ? c[k] + v : v
    end
    Dict(k => v for (k, v) in c if !isz(v))
end

function mulop(p::Trig{D}, L::Dict{NTuple{D, Int}, Trig{D}}) where {D}
    Dict(k => p * m for (k, m) in L if !isz(p * m))
end

Base.:+(a::Var{D}, b::Var{D}) where {D} = Var{D}(a.v + b.v, addop(a.L, b.L))
Base.:-(a::Var{D}) where {D} = Var{D}(-a.v, Dict(k => -m for (k, m) in a.L))
Base.:-(a::Var{D}, b::Var{D}) where {D} = a + (-b)
function Base.:*(a::Var{D}, b::Var{D}) where {D}
    Var{D}(a.v * b.v, addop(mulop(a.v, b.L), mulop(b.v, a.L)))
end
function scale(s, a::Var{D}) where {D}
    Var{D}(scale(s, a.v), Dict(k => scale(s, m) for (k, m) in a.L))
end

"d/dz_j applied to a field and, by the product rule, to its variation operator."
function der(a::Var{D}, j::Int) where {D}
    L = Dict{NTuple{D, Int}, Trig{D}}()
    for (α, m) in a.L
        L = addop(L, Dict(α => der(m, j)))
        β = ntuple(i -> α[i] + (i == j ? 1 : 0), D)
        L = addop(L, Dict(β => m))
    end
    Var{D}(der(a.v, j), L)
end

"Horner evaluation of a polynomial in the field; coefficients are fields (so c may depend on z)."
function horner(coef::Vector{Trig{D}}, f::Var{D}) where {D}
    r = lift(coef[end])
    for i in (length(coef) - 1):-1:1
        r = lift(coef[i]) + f * r
    end
    r
end
function horner(coef::Vector{Trig{D}}, f::Trig{D}) where {D}
    isempty(coef) && return zerotrig(D)      # c_u of an f-independent c
    r = coef[end]
    for i in (length(coef) - 1):-1:1
        r = coef[i] + f * r
    end
    r
end

"""
The functional derivative of the functional f |-> < rho(f) > from the density's variation operator.

< sum_alpha m_alpha d^alpha v > = < (sum_alpha (-1)^|alpha| d^alpha m_alpha) v >, exactly on the
torus, because the mean of a coordinate derivative vanishes.  Integration by parts and nothing else.
"""
function fderiv(rho::Var{D}) where {D}
    out = zerotrig(D)
    for (α, m) in rho.L
        t = m
        for j in 1:D, _ in 1:α[j]

            t = der(t, j)
        end
        out = out + scale((-1)^sum(α), t)
    end
    out
end

# ---------------------------------------------------------------------------
# the constant Poisson bivector and the four-vector combination S
# ---------------------------------------------------------------------------

"Pi(dg,dh) = Pi^{ij} d_i g d_j h for a constant antisymmetric Pi.  Generic in the field type."
function pib(Π, g, h)
    D = size(Π, 1)
    acc = nothing
    for i in 1:D, j in 1:D

        iszero(Π[i, j]) && continue
        t = scale(Π[i, j], der(g, i) * der(h, j))
        acc = acc === nothing ? t : acc + t
    end
    acc === nothing ? (g isa Var ? lift(zerotrig(D)) : zerotrig(D)) : acc
end

"S(a,b,g,d) = Pi(a,b)Pi(g,d) + Pi(b,g)Pi(a,d) + Pi(g,a)Pi(b,d), on differentials of fields."
function Sfour(Π, A, B, Cc, G)
    pib(Π, A, B) * pib(Π, Cc, G) + pib(Π, B, Cc) * pib(Π, A, G) +
    pib(Π, Cc, A) * pib(Π, B, G)
end

"Canonical Pi on T^{2n} in coordinates (x_1..x_n, v_1..v_n)."
function canonical_bivector(n)
    Π = zeros(Int, 2n, 2n)
    for k in 1:n
        Π[k, n + k] = 1
        Π[n + k, k] = -1
    end
    Π
end

"Constant Pi on T^D of rank 2r, supported on the first 2r coordinates: degenerate for 2r < D."
function degenerate_bivector(D, r)
    Π = zeros(Int, D, D)
    for k in 1:r
        Π[2k - 1, 2k] = 1
        Π[2k, 2k - 1] = -1
    end
    Π
end

# ---------------------------------------------------------------------------
# a minimal exterior algebra over Q, for the pointwise identity and Lefschetz
# ---------------------------------------------------------------------------
# A p-form is a Dict from a strictly increasing index tuple to its coefficient.

const Form = Dict{Vector{Int}, R}

function wedge(a::Form, b::Form)
    out = Form()
    for (ia, va) in a, (ib, vb) in b

        idx = vcat(ia, ib)
        length(unique(idx)) == length(idx) || continue
        s = 1
        # insertion sort, counting transpositions
        idx = copy(idx)
        for i in 2:length(idx)
            j = i
            while j > 1 && idx[j] < idx[j - 1]
                idx[j], idx[j - 1] = idx[j - 1], idx[j]
                s = -s
                j -= 1
            end
        end
        out[idx] = get(out, idx, zero(R)) + s * va * vb
    end
    Form(k => v for (k, v) in out if !iszero(v))
end

covector(v::Vector) = Form([i] => R(v[i]) for i in eachindex(v) if !iszero(v[i]))
formpow(a::Form, k::Int) = k == 0 ? Form(Int[] => one(R)) : wedge(a, formpow(a, k - 1))

"Strictly increasing p-tuples drawn from 1:D."
function subsets(D::Int, p::Int)
    p == 0 && return [Int[]]
    out = Vector{Int}[]
    function rec(start, acc)
        length(acc) == p && (push!(out, copy(acc)); return)
        for i in start:D
            push!(acc, i)
            rec(i + 1, acc)
            pop!(acc)
        end
    end
    rec(1, Int[])
    out
end

"Rank by exact Gaussian elimination over Q; `rank` from LinearAlgebra needs BLAS floats."
function exactrank(Mx::Matrix{R})
    A = copy(Mx)
    m, n = size(A)
    r = 0
    for c in 1:n
        piv = findfirst(i -> !iszero(A[i, c]), (r + 1):m)
        piv === nothing && continue
        r += 1
        piv += r - 1
        A[[r, piv], :] = A[[piv, r], :]
        A[r, :] = A[r, :] .// A[r, c]
        for i in 1:m
            i == r && continue
            iszero(A[i, c]) || (A[i, :] = A[i, :] .- A[i, c] .* A[r, :])
        end
    end
    r
end

"omega = sum_k dx_k ^ dy_k in coordinates (x_1..x_n, v_1..v_n)."
symplectic_form(n) = Form([k, n + k] => one(R) for k in 1:n)

"S(alpha,beta,gamma,delta) for constant Pi and covectors given as coefficient vectors."
function Snum(Π, α, β, γ, δ)
    p(u, v) = sum(Π[i, j] * u[i] * v[j] for i in axes(Π, 1), j in axes(Π, 2); init = zero(R))
    p(α, β) * p(γ, δ) + p(β, γ) * p(α, δ) + p(γ, α) * p(β, δ)
end

"""
The (D-4)-form Theta forced by  S(a,b,g,d) mu = a ^ b ^ g ^ d ^ Theta  on basis covectors,
with `vol` the coefficient of the chosen volume form on e^1 ^ ... ^ e^D.

Reading Theta off the basis is a definition, not a check.  The check is that the resulting Theta
then reproduces S for *arbitrary* covectors, and that it equals the closed form III.1 claims.
Getting `vol` right matters: e^1^...^e^{2n} is NOT the Liouville form omega^n/n! in the coordinate
order (x_1..x_n, v_1..v_n) -- they differ by (-1)^{n(n-1)/2}, and reading the identity against the
wrong one manufactures a spurious n-dependent sign.
"""
function theta_from_S(Π, vol::R = one(R))
    D = size(Π, 1)
    Θ = Form()
    e(i) = R[j == i ? 1 : 0 for j in 1:D]
    for i in 1:D, j in (i + 1):D, k in (j + 1):D, l in (k + 1):D
        s = Snum(Π, e(i), e(j), e(k), e(l))
        iszero(s) && continue
        comp = setdiff(1:D, [i, j, k, l])
        sgn = wedge(Form([i, j, k, l] => one(R)), Form(comp => one(R)))
        Θ[comp] = get(Θ, comp, zero(R)) + s * vol * sgn[collect(1:D)]
    end
    Form(k => v for (k, v) in Θ if !iszero(v))
end

# -----------------------------------------------------------------------
# attack 1: the Jacobiator itself, brute force
# -----------------------------------------------------------------------
# A model is (Pi, w, c) with c(z,u) = sum_j ccoef[j] u^{j-1}; ccoef[j] are fields, so c may
# depend on z as well as on f.  w = 1 and z-independent c is the case III.1 claims.

struct Model{D}
    Π::Matrix{Int}
    w::Trig{D}
    ccoef::Vector{Trig{D}}
end

cu_coef(M::Model{D}) where {D} = [scale(j, M.ccoef[j + 1]) for j in 1:(length(M.ccoef) - 1)]

"The bracket integrand as a Var, so that its variation operator comes with it."
function integrand(M::Model{D}, f::Var{D}, X::Var{D}, Y::Var{D}) where {D}
    lift(M.w) * horner(M.ccoef, f) * pib(M.Π, X, Y)
end

"""
The Jacobiator, brute force.

delta{A,B}/delta f is obtained from the variation operator of the bracket integrand by
integration by parts (`fderiv`) -- no closed form is assumed anywhere.  `Afd` etc. are the
functional derivatives A_f as expressions in f, so their own variation operators are A_ff.
"""
function jacobiator(M::Model{D}, f0::Trig{D}, Afd, Bfd, Cfd) where {D}
    f = seed(f0)
    A, B, Cc = Afd(f), Bfd(f), Cfd(f)
    cval = horner(M.ccoef, f0)
    total = zero(R)
    for (X, Y, Z) in ((A, B, Cc), (B, Cc, A), (Cc, A, B))
        DXY = fderiv(integrand(M, f, X, Y))            # delta{X,Y}/delta f, exactly
        total += meanval(M.w * cval * pib(M.Π, DXY, Z.v))
    end
    total
end

"The master formula of the refutation note: - < w c S(dA_f,dB_f,dC_f, D(w c_u)) >."
function jac_predicted(M::Model{D}, f0::Trig{D}, Afd, Bfd, Cfd) where {D}
    f = seed(f0)
    A, B, Cc = Afd(f).v, Bfd(f).v, Cfd(f).v
    wcu = M.w * horner(cu_coef(M), f0)
    -meanval(M.w * horner(M.ccoef, f0) * Sfour(M.Π, A, B, Cc, wcu))
end

"The functionals used throughout: A_f, B_f, C_f as expressions in f, with C_f nonlocal."
function testfunctionals(D, a1, a2, b1, b2, g1, g2)
    Afd = f -> lift(a1) * f * f + lift(a2) * f
    Bfd = f -> lift(b1) * f + lift(b2)
    Cfd = f -> lift(g1) * f * f + der(lift(g2) * der(f, 1), 1)
    # the densities they come from, for the consistency check below
    ρA = f -> scale(1 // 3, lift(a1) * f * f * f) + scale(1 // 2, lift(a2) * f * f)
    ρB = f -> scale(1 // 2, lift(b1) * f * f) + lift(b2) * f
    ρC = f -> scale(1 // 3, lift(g1) * f * f * f) -
              scale(1 // 2, lift(g2) * der(f, 1) * der(f, 1))
    (Afd, Bfd, Cfd), (ρA, ρB, ρC)
end

println("=" ^ 78)
println("III.1 adversarial re-check -- exact arithmetic over Q(i)")
println("=" ^ 78)

@testset "III.1 refutation attempt" begin

    # -----------------------------------------------------------------------
    @testset "attack 3: the pointwise identity and its combinatorial factor" begin
        # S(a,b,g,d) mu = a^b^g^d^Theta with Theta constant; and in the symplectic case
        # Theta = omega^{n-2}/(n-2)!.  Both halves checked exactly.
        rng = [3, -2, 5, 1, -4, 2, 7, -1, 3, -5, 2, 4, -3, 6, 1, -2]
        for n in 2:4
            D = 2n
            Π = canonical_bivector(n)
            ω = symplectic_form(n)
            # mu = Liouville = omega^n/n!, expressed against the coordinate top form
            liou = get(formpow(ω, n), collect(1:D), zero(R)) // factorial(n)
            @test liou == (-1)^(n * (n - 1) ÷ 2)
            Θ = theta_from_S(Π, liou)
            # (i) Theta reproduces S on arbitrary (not merely basis) covectors
            for trial in 1:6
                vs = [R[rng[mod1(trial * 4 + 3i + q, length(rng))] for i in 1:D]
                      for q in 1:4]
                lhs = Snum(Π, vs...) * liou
                w4 = wedge(wedge(covector(vs[1]), covector(vs[2])),
                    wedge(covector(vs[3]), covector(vs[4])))
                rhs = get(wedge(w4, Θ), collect(1:D), zero(R))
                @test lhs == rhs
            end
            # (ii) the claimed closed form of Theta, factorial and sign included
            claimed = Form(k => v // factorial(n - 2) for (k, v) in formpow(ω, n - 2))
            @test Θ == claimed
            println(
                "  n = $n (dim $D):  S (omega^n/n!) = a^b^g^d ^ omega^$(n-2)/$(n - 2)!",
                "   exact, sign +1")
        end
        # The same statement for a *degenerate* constant Pi on a torus: Theta = iota_{Pi^Pi/2} mu is
        # still a constant, hence closed, (D-4)-form, which is all the Stokes argument needs.
        for (D, r) in [(5, 2), (6, 2), (7, 3), (6, 3)]
            Π = degenerate_bivector(D, r)
            Θ = theta_from_S(Π)
            for trial in 1:5
                vs = [R[rng[mod1(trial * 5 + 2i + q, length(rng))] for i in 1:D]
                      for q in 1:4]
                w4 = wedge(wedge(covector(vs[1]), covector(vs[2])),
                    wedge(covector(vs[3]), covector(vs[4])))
                @test Snum(Π, vs...) == get(wedge(w4, Θ), collect(1:D), zero(R))
            end
            println(
                "  D = $D, rank Pi = $(2r):  S mu = a^b^g^d ^ Theta with Theta constant",
                r == 1 ? "" : "", "  (Theta ", isempty(Θ) ? "= 0" : "!= 0", ")")
        end
    end

    # -----------------------------------------------------------------------
    @testset "attack 2: Lefschetz injectivity and the graph lemma" begin
        for n in 2:4
            D = 2n
            ω = symplectic_form(n)
            ωk = formpow(ω, n - 2)
            basis2 = subsets(D, 2)
            basis_hi = subsets(D, 2n - 2)
            Mx = zeros(R, length(basis_hi), length(basis2))
            for (c, b) in enumerate(basis2)
                img = wedge(Form(b => one(R)), ωk)
                for (r, h) in enumerate(basis_hi)
                    Mx[r, c] = get(img, h, zero(R))
                end
            end
            @test exactrank(Mx) == length(basis2)     # wedging with omega^{n-2} is injective
            println("  n = $n:  ^omega^$(n-2) : Lambda^2 -> Lambda^$(2n-2) has rank ",
                exactrank(Mx), " = dim Lambda^2 = ", length(basis2))
        end

        # A 2-form on M x R vanishing on the tangent space of every graph {(v, p(v))} vanishes.
        # Basis of Lambda^2(R^{D+1}): pairs from 1..D, and (i, D+1).
        for D in 2:5
            basis2 = subsets(D + 1, 2)
            ps = Vector{R}[]
            push!(ps, R[0 for _ in 1:D])
            for i in 1:D
                q = R[0 for _ in 1:D]
                q[i] = 1
                push!(ps, q)
            end
            rows = Vector{R}[]
            for p in ps
                # pullback under v -> (v, p.v): basis e_a of R^D maps to e_a + p_a e_{D+1}
                for a in 1:D, b in (a + 1):D

                    row = R[]
                    for bb in basis2
                        # evaluate the basis 2-form dz^i ^ dz^j on (e_a + p_a e_{D+1}, e_b + p_b e_{D+1})
                        i, j = bb
                        ci(k) = k == a ? one(R) : (k == D + 1 ? p[a] : zero(R))
                        cj(k) = k == b ? one(R) : (k == D + 1 ? p[b] : zero(R))
                        push!(row, ci(i) * cj(j) - ci(j) * cj(i))
                    end
                    push!(rows, row)
                end
            end
            Mx = permutedims(reduce(hcat, rows))
            @test exactrank(Mx) == length(basis2)
            println(
                "  D = $D:  2-forms on M x R vanishing on all graphs: kernel is trivial ",
                "(rank ", exactrank(Mx), " = ", length(basis2), ")")
        end
    end

    @testset "machinery: the checker computes what it says it does" begin
        D = 4
        f0 = harm((1, 0, 0, 0), 2, :cos) + harm((0, 1, 1, 0), 3, :sin) +
             harm((0, 0, 1, 1), 1, :cos)
        a1 = harm((0, 1, 0, 0), 1, :cos) + konst(D, 2)
        a2 = harm((1, 0, 1, 0), 3, :sin)
        b1 = harm((0, 0, 0, 1), 2, :sin) + konst(D, 1)
        b2 = harm((1, 1, 0, 0), 1, :cos)
        g1 = harm((0, 1, 0, 1), 2, :cos) + konst(D, 3)
        g2 = harm((1, 0, 0, 1), 1, :sin) + konst(D, 2)
        (Afd, Bfd, Cfd), (ρA, ρB, ρC) = testfunctionals(D, a1, a2, b1, b2, g1, g2)
        f = seed(f0)

        # (i) each A_f really is the functional derivative of the density it is claimed to come from
        for (fd, ρ) in zip((Afd, Bfd, Cfd), (ρA, ρB, ρC))
            @test isz(fderiv(ρ(f)) - fd(f).v)
        end
        # (ii) each second variation is symmetric -- the property Step 1 of the note leans on
        u = harm((1, 1, 1, 0), 5, :cos) + harm((0, 1, 0, 0), 2, :sin)
        v = harm((0, 0, 1, 0), 3, :sin) + harm((1, 0, 0, 1), 4, :cos)
        applyop(L, x) = begin
            out = zerotrig(D)
            for (α, m) in L
                t = x
                for j in 1:D, _ in 1:α[j]

                    t = der(t, j)
                end
                out = out + m * t
            end
            out
        end
        for fd in (Afd, Bfd, Cfd)
            L = fd(f).L
            @test meanval(u * applyop(L, v)) == meanval(v * applyop(L, u))
        end
        # (iii) nothing here is accidentally zero
        @test !isz(f0) && !isz(Afd(f).v) && !isz(Bfd(f).v) && !isz(Cfd(f).v)
        @test !isz(Sfour(canonical_bivector(2), Afd(f).v, Bfd(f).v, Cfd(f).v, f0))
        println("  machinery: A_f matches its density, A_ff symmetric, ",
            "S(dA_f,dB_f,dC_f,df) not identically zero")
    end

    # -----------------------------------------------------------------------
    @testset "attack 1: the Jacobiator, brute force" begin
        """
        Fields for T^D.  The wavevectors must touch EVERY coordinate: with Pi canonical on T^6 and
        fields supported on z_1..z_4 only, the (2,5) and (3,6) planes drop out, Pi acts at rank 2,
        and S vanishes identically -- an accidental zero that makes the test vacuous.  This script
        hit exactly that on the first run, so `run` now asserts non-vacuity.
        """
        function setup(D)
            # Nine wavevectors: the D coordinate directions first, then skew pairs.  Two things
            # must hold and neither is automatic.  (a) Every coordinate is excited -- otherwise the
            # unexcited half of a canonical pair drops out of Pi and its rank collapses.  (b) The
            # wavevectors span at least four independent directions -- the pool e_i + e_{i+1} looks
            # fine under (a) and fails (b), since those four vectors satisfy
            # (e1+e2) - (e2+e3) + (e3+e4) - (e4+e1) = 0, so every 4-fold antisymmetrisation dies.
            # Both failures produce S == 0 and a Jacobiator of exactly zero that means nothing.
            k(i) = i <= D ? ntuple(l -> l == i ? 1 : 0, D) :
                   ntuple(l -> (l == mod1(i, D)) + (l == mod1(i + 2, D)), D)
            f0 = harm(k(1), 2, :cos) + harm(k(2), 3, :sin) + harm(k(3), 1, :cos)
            a1 = harm(k(4), 1, :cos) + konst(D, 2)
            a2 = harm(k(5), 3, :sin)
            b1 = harm(k(6), 2, :sin) + konst(D, 1)
            b2 = harm(k(7), 1, :cos)
            g1 = harm(k(8), 2, :cos) + konst(D, 3)
            g2 = harm(k(9), 1, :sin) + konst(D, 2)
            f0, testfunctionals(D, a1, a2, b1, b2, g1, g2)[1]
        end

        "Run one case and report; `expect` is :zero, :nonzero, or :zero_vacuous (S is expected to
        vanish identically, so the case carries no information beyond consistency)."
        function run(name, M, f0, fds, expect)
            jac = jacobiator(M, f0, fds...)
            pred = jac_predicted(M, f0, fds...)
            A, B, Cc = (fd(seed(f0)).v for fd in fds)
            wcu = M.w * horner(cu_coef(M), f0)
            snz = !isz(Sfour(M.Π, A, B, Cc, f0))
            @test jac == pred                       # brute force vs the master formula
            if expect === :nonzero
                @test jac != 0                      # the checker is capable of a nonzero verdict
                @test snz
            else
                @test jac == 0
                # a vanishing verdict is worth nothing unless the integrand could have been nonzero
                expect === :zero ? (@test snz) : (@test !snz)
            end
            println("  ", rpad(name, 52), "Jac = ", jac,
                "   (master formula: ", pred, ")",
                snz ? "" : "   [S identically zero -- vacuous]")
            jac
        end

        # --- the claim itself -------------------------------------------------------------
        for (D, Π, tag, ex) in [(4, canonical_bivector(2), "T^4 canonical", :zero),
            (6, canonical_bivector(3), "T^6 canonical", :zero),
            (8, canonical_bivector(4), "T^8 canonical", :zero),
            (6, degenerate_bivector(6, 2), "T^6 rank-4 degenerate", :zero),
            (5, degenerate_bivector(5, 2),
                "T^5 rank-4 degenerate (no symplectic form)", :zero),
            (4, degenerate_bivector(4, 1), "T^4 rank-2 degenerate", :zero_vacuous),
            (2, canonical_bivector(1), "T^2 canonical (Plucker)", :zero_vacuous)]
            f0, fds = setup(D)
            # c cubic in f, with an awkward shape on purpose: c(0) = 0, c' with a real zero
            cc = [konst(D, 0), konst(D, 3), konst(D, -2), konst(D, 5)]
            run("$tag, w = 1, c cubic", Model{D}(Π, konst(D, 1), cc), f0, fds, ex)
        end

        # c of high degree, large coefficients: "c with a big derivative" changes nothing,
        # because f has compact range on a closed M and c c'' df is still exact.
        let D = 4
            f0, fds = setup(D)
            cc = [konst(D, 7), konst(D, -11), konst(D, 13), konst(D, -17), konst(D, 19),
                konst(D, 23), konst(D, -29)]
            run("T^4 canonical, w = 1, c of degree 6",
                Model{D}(canonical_bivector(2),
                    konst(D, 1), cc), f0, fds, :zero)
        end

        # --- the controls: cases where the mechanism is broken on purpose --------------------
        # (a) a non-Liouville volume form.  Theta = w omega^{n-2}/(n-2)! is then not closed and the
        #     Jacobiator must NOT vanish.  This is what proves the checker is not a zero machine.
        let D = 4
            f0, fds = setup(D)
            cc = [konst(D, 0), konst(D, 3), konst(D, -2), konst(D, 5)]
            w = konst(D, 5) + harm((1, 0, 0, 0), 2, :cos) + harm((0, 0, 1, 0), 3, :sin)
            j = run("CONTROL T^4, non-constant weight w(z)",
                Model{D}(canonical_bivector(2), w, cc),
                f0, fds, :nonzero)
            # and the two candidate densities are genuinely different functions here
            A, B, Cc = (fd(seed(f0)).v for fd in fds)
            cpr = horner(cu_coef(Model{D}(canonical_bivector(2), w, cc)), f0)
            cval = horner(cc, f0)
            cpp = horner(
                [scale(j2, cu_coef(Model{D}(canonical_bivector(2), w, cc))[j2 + 1])
                 for j2 in 1:(length(cc) - 2)],
                f0)
            S = Sfour(canonical_bivector(2), A, B, Cc, f0)
            naive_ccpp = -meanval(w * w * cval * cpp * S)
            naive_cpr2 = meanval(w * w * cpr * cpr * S)
            println("      (c')^2 form: ", naive_cpr2, "   c c'' form: ", naive_ccpp,
                "   -- equal only when the weight is constant")
        end

        # (b) c depending on z as well as on f, violating dc ^ dc_u = 0: must NOT vanish.
        let D = 4
            f0, fds = setup(D)
            # The z-dependence must be rich enough for a wavevector pentagon to close: with one
            # harmonic per coefficient this control returned an exact zero, i.e. it failed to be a
            # control at all -- the trap III.1 documents, hit here in the other negative.
            γ0 = konst(D, 1) + harm((1, 0, 0, 0), 2, :cos) + harm((1, 1, 0, 0), 1, :sin)
            γ1 = konst(D, 2) + harm((0, 1, 0, 0), 3, :sin) + harm((0, 1, 1, 0), 1, :cos) +
                 harm((1, 0, 1, 0), 1, :sin)
            γ2 = harm((0, 0, 1, 0), 1, :cos) + harm((0, 0, 1, 1), 1, :sin) +
                 harm((1, 0, 0, 1), 1, :cos)
            run("CONTROL T^4, c(z,f) with dc ^ dc_u != 0",
                Model{D}(canonical_bivector(2), konst(D, 1), [γ0, γ1, γ2]), f0, fds, :nonzero)
        end

        # --- the companion proposition's admissible families ---------------------------------
        # c = Phi(f + phi(z)) with Phi(s) = s^3 : ccoef from the binomial expansion.
        let D = 4
            f0, fds = setup(D)
            φ = harm((1, 0, 0, 0), 2, :sin) + harm((0, 0, 1, 0), 1, :cos) + konst(D, 3)
            cc = [φ * φ * φ, scale(3, φ * φ), scale(3, φ), konst(D, 1)]
            run("T^4, c = (f + phi(z))^3   [dc ^ dc_u = 0]",
                Model{D}(canonical_bivector(2), konst(D, 1), cc), f0, fds, :zero)
        end
        # c = c(z) alone: c_u = 0, so the Jacobiator vanishes identically.
        let D = 4
            f0, fds = setup(D)
            cc = [harm((1, 1, 0, 0), 2, :cos) + konst(D, 4)]
            run("T^4, c = c(z) only", Model{D}(canonical_bivector(2), konst(D, 1), cc),
                f0, fds, :zero)
        end
    end

    # -----------------------------------------------------------------------
    @testset "attack 2b: the classification clause of the c(z,f) proposition" begin
        # The proposition's condition dc_u ^ dc = 0 is fine.  Its reading of that condition --
        # "where c_u != 0 this means c(z,u) = h(u + phi(z))" -- is a LOCAL statement asserted
        # globally.  On a manifold with H^1 != 0 it is false, and T^{2n} is such a manifold.
        #
        # Counterexample: c(z,u) = sin(u + z_1).  Both u and z_1 enter 2pi-periodically, so c is a
        # genuine smooth function on T^D x R (indeed on T^D x S^1); here u is carried as the extra
        # coordinate z_{D+1} so that everything stays a trigonometric polynomial and exact.
        D = 4
        E = D + 1
        u = E                                              # the u slot
        c = harm(ntuple(j -> (j == u) + (j == 1), E), 1, :sin)     # sin(u + z_1)
        cu = der(c, u)                                     # cos(u + z_1)
        @test !isz(cu)
        # dc_u ^ dc = 0 on M x R, componentwise and exactly
        allzero = true
        for i in 1:E, j in (i + 1):E

            allzero &= isz(der(cu, i) * der(c, j) - der(cu, j) * der(c, i))
        end
        @test allzero
        # It is NOT of the form h(u + phi(z)) for any global smooth phi on T^D: that would force
        # phi = z_1 mod 2pi, which has winding number one around the first circle and so admits no
        # continuous lift.  What fails concretely is exactness of eta = c_u d(c o graph): take the
        # constant field f = 0, so eta = cos^2(z_1) dz_1, whose z_1-component has mean 1/2.  Every
        # component of an exact one-form d(psi) has mean zero on the torus, so eta is closed but
        # not exact -- and the sufficiency proof only ever needs it closed.
        ctilde = harm(ntuple(j -> j == 1 ? 1 : 0, D), 1, :sin)         # sin(z_1), i.e. f = 0
        cu_on = harm(ntuple(j -> j == 1 ? 1 : 0, D), 1, :cos)          # cos(z_1)
        η1 = cu_on * der(ctilde, 1)                                     # cos^2(z_1)
        @test meanval(η1) == 1 // 2
        @test meanval(der(ctilde, 1)) == 0          # the same functional kills any exact one-form
        println("  c(z,u) = sin(u + z_1) on T^4 x R:  dc_u ^ dc = 0 exactly, c_u != 0,")
        println("      but eta = c_u dc~ has non-zero period (mean of its z_1 component = 1/2),")
        println("      so c is admissible and is NOT h(u + phi(z)) for any global phi.")
    end
end
