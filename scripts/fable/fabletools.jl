#
# Exact-arithmetic machinery shared by the fable scripts: trigonometric polynomials over Q(i)
# on the torus T^D, first-order jets for Gateaux derivatives, and the canonical bracket on
# T^{2n}.  Plain definitions, not a module, so that a script can extend `∂`, `mean`, `phasedim`
# to its own field types (III1 does so for its floating-point grid fields).
#
#     include(joinpath(@__DIR__, "fabletools.jl"))
#
# Written for III1_structure_function_4d.jl and factored out for III4_multifield_derivative.jl.

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
# the same code computes the exact Jacobiator on TP and a floating-point one on grid fields.

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

# A field times a Jet of scalars (`phi * mean(phi * F)` for a nonlocal functional derivative),
# and a Jet plus a field.
Base.:*(a::TP, s::Jet) = Jet(a * s.v, a * s.d)
Base.:*(s::Jet, a::TP) = a * s
Base.:+(a::Jet, b::TP) = Jet(a.v + b, a.d)
Base.:+(b::TP, a::Jet) = a + b

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

"The cyclic combination S(dA,dB,dC,dG) = {A,B}{C,G} + {B,C}{A,G} + {C,A}{B,G}."
Sfun(A, B, C, G) = pb(A, B) * pb(C, G) + pb(B, C) * pb(A, G) + pb(C, A) * pb(B, G)
