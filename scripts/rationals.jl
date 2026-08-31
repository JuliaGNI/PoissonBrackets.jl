#
# Random exact rationals, the test data the algebraic verification scripts run on.
#
#     include(joinpath(@__DIR__, "rationals.jl"))
#     using .Rationals
#
# These mirror `rnd`/`rnd_vec`/`rnd_mat` from the `common.py` they replace, including the
# default numerator range -6:6 and denominator range 1:4, so residuals land in the same
# ballpark as the Python's. They will NOT be the same numbers: Julia's random stream is not
# Python's, and no attempt is made to make it so. What is compared between the two is the
# verdict of each check and the order of magnitude of each residual, never its digits.
#
# Where a claim can be settled on deterministic data instead -- se(3), so(N), a fixed mesh --
# prefer that; those DO agree with the Python exactly, and are the better witness anyway.

module Rationals

using Random

export rnd, rnd_vec, rnd_mat, positive_vec

const Q = Rational{BigInt}

"""
    rnd(rng; lo = -6, hi = 6, dmax = 4)

A random rational with numerator in `lo:hi` and denominator in `1:dmax`. May be zero.
"""
function rnd(rng::AbstractRNG; lo::Int = -6, hi::Int = 6, dmax::Int = 4)
    Q(rand(rng, lo:hi)) // Q(rand(rng, 1:dmax))
end

rnd_vec(rng::AbstractRNG, n::Int; kwargs...) = Q[rnd(rng; kwargs...) for _ in 1:n]

"""
    rnd_mat(rng, n, kind = :general)

A random `n × n` rational matrix. `kind` is one of `:general`, `:symmetric`,
`:antisymmetric` or `:diagonal`.
"""
function rnd_mat(rng::AbstractRNG, n::Int, kind::Symbol = :general; kwargs...)
    M = zeros(Q, n, n)
    if kind === :general
        for i in 1:n, j in 1:n

            M[i, j] = rnd(rng; kwargs...)
        end
    elseif kind === :diagonal
        for i in 1:n
            M[i, i] = rnd(rng; kwargs...)
        end
    elseif kind === :symmetric
        for i in 1:n, j in i:n

            M[i, j] = M[j, i] = rnd(rng; kwargs...)
        end
    elseif kind === :antisymmetric
        for i in 1:n, j in (i + 1):n

            v = rnd(rng; kwargs...)
            M[i, j], M[j, i] = v, -v
        end
    else
        throw(ArgumentError("unknown kind $kind; expected :general, :symmetric, " *
                            ":antisymmetric or :diagonal"))
    end
    return M
end

"""
    positive_vec(rng, n; kwargs...)

A random rational vector with every entry strictly positive — the admissible region for the
brackets gauged by `√u`, where `u_i ≤ 0` is not merely inconvenient but outside the domain.
"""
function positive_vec(rng::AbstractRNG, n::Int; kwargs...)
    Q[(abs(numerator(x)) + 1) // denominator(x) for x in rnd_vec(rng, n; kwargs...)]
end

end # module
