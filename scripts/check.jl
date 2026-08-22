#
# The PASS/FAIL harness the verification scripts report through.
#
# The output format is copied deliberately from `common.py` in the manuscripts' `Scripts/`
# folders, down to the two leading spaces and the three spaces before the detail, so that
# a converted script can be checked against the Python it replaces with a plain `diff`.
# Once the Python is retired that constraint is free to go, but until then it is what makes
# the comparison mechanical rather than a matter of reading two tables side by side.
#
#     include(joinpath(@__DIR__, "check.jl"))
#     using .Checks: header, check, summary
#
# Import the three names explicitly rather than with a bare `using`: `summary` also exists in
# `Base`, and a bare `using` leaves the two ambiguous rather than preferring ours.
#
#     header("1. the bracket is antisymmetric")
#     check("J + J' = 0", norm(J + J') < 1e-14, @sprintf("%.2e", norm(J + J')))
#     summary("verify_something.jl")
#
# `summary` exits 1 if any check failed, so `run_all.jl` and CI see a nonzero status.

module Checks

export header, check, summary, fmt, failures, reset_failures!

const _failures = String[]

"""
    header(text)

A blank line, `text`, and a rule of matching length.  Sections of a verification script are
numbered in the Python; keep the numbering so the two outputs line up.
"""
function header(text)
    println()
    println(text)
    println("-"^length(text))
    return nothing
end

"""
    check(label, condition, detail = "")

Record and print one claim.  Returns the boolean, so a check can gate what follows it.

`detail` is the measured quantity behind the verdict -- a residual, a rate, a rank.  A bare
PASS says only that someone chose the tolerance well; the number says what was actually
seen, and it is the number that gets compared against the Python.
"""
function check(label, condition, detail = "")
    ok = Bool(condition)
    ok || push!(_failures, label)
    mark = ok ? "PASS" : "FAIL"
    println("  [$mark] $label" * (isempty(detail) ? "" : "   $detail"))
    return ok
end

"""
    summary(title)

Print the tally and exit 1 if anything failed.  Call it once, at the end of a script.
"""
function summary(title)
    println()
    if !isempty(_failures)
        println("$title: $(length(_failures)) FAILURE(S): " * join(_failures, ", "))
        exit(1)
    end
    println("$title: all checks passed.")
    return nothing
end

"""
    fmt(x)

Render a value for the `detail` field of [`check`](@ref).

Rationals print as `p/q`, and as `p` when the denominator is one, rather than in Julia's
`p//q` form. That is not a cosmetic whim: it is the form the Python prototypes print, and
keeping it means a `diff` of the two outputs shows only the numbers that genuinely differ.
"""
fmt(x::Rational) = isone(denominator(x)) ? string(numerator(x)) :
                   string(numerator(x), "/", denominator(x))
fmt(x::AbstractFloat) = string(x)
fmt(x) = string(x)

"The labels of the checks that have failed so far."
failures() = copy(_failures)

"Forget the recorded failures.  Only `run_all.jl` needs this, between scripts."
reset_failures!() = (empty!(_failures); nothing)

end # module
