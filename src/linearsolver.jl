
@doc raw"""
    LapackLU()

A LAPACK-backed LU factorisation, plugged into SimpleSolvers through its linear-solver
extension points.

SimpleSolvers ships one linear solver, a hand-written scalar LU. It is portable, works for
static matrices and arbitrary number types, and is the right default for the small dense
systems that package is usually pointed at. At the sizes these discretisations reach it is
the wrong tool: at ``N = 384`` the factorisation of the Newton matrix is ``N^3/3 \approx 2
\times 10^7`` bounds-checked scalar operations, which measured **74 % of the cost of an
implicit step** — around 17 ms against LAPACK's 0.6 ms for the same factorisation.

This method keeps SimpleSolvers as the nonlinear driver, with its quasi-Newton
refactorisation logic, convergence tests and status reporting, and replaces only the
factorisation:

```julia
Integrator(flow, ImplicitMidpoint(), Δt; linear_solver_method = LapackLU())
```

which is the default here. Pass `SimpleSolvers.LU()` to get the original behaviour back.

The factorisation is `lu!` on a cached copy of the matrix, so it allocates only the pivot
vector that LAPACK requires, and reuses the same storage on every refactorisation.
"""
struct LapackLU <: SimpleSolvers.LinearSolverMethod end

"""
    LapackLUCache{T}

Storage for [`LapackLU`](@ref): the working copy of the matrix that `lu!` overwrites, and
the factorisation object built on it.

Mutable, because the factorisation is replaced on every refactorisation while the matrix
storage it points into is reused.
"""
mutable struct LapackLUCache{T} <: SimpleSolvers.LinearSolverCache{T}
    A::Matrix{T}
    fact::Union{Nothing, LinearAlgebra.LU{T, Matrix{T}, Vector{Int}}}
end

function SimpleSolvers.LinearSolverCache(::LapackLU, A::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    LapackLUCache{T}(Matrix{T}(undef, n, n), nothing)
end

function SimpleSolvers.factorize!(ls::SimpleSolvers.LinearSolver{T, LapackLU},
                                  A::AbstractMatrix) where {T}
    c = SimpleSolvers.cache(ls)
    copyto!(c.A, A)
    c.fact = LinearAlgebra.lu!(c.A; check = false)
    return ls
end

SimpleSolvers.factorize!(ls::SimpleSolvers.LinearSolver{T, LapackLU},
                         lp::SimpleSolvers.LinearProblem) where {T} =
    SimpleSolvers.factorize!(ls, lp.A)

function LinearAlgebra.ldiv!(x::AbstractVector, ls::SimpleSolvers.LinearSolver{T, LapackLU},
                             b::AbstractVector) where {T}
    c = SimpleSolvers.cache(ls)
    c.fact === nothing && throw(ArgumentError(
        "the linear solver has not been factorized yet; call factorize! first"))
    # A singular Newton matrix is a statement about the step size, not a bug, so it is
    # reported as such rather than left to LAPACK's index-of-the-zero-pivot message.
    LinearAlgebra.issuccess(c.fact) || throw(SingularException(
        "the Newton matrix is singular; the step size is probably too large for the " *
        "stiffness of this discretisation"))
    copyto!(x, b)
    LinearAlgebra.ldiv!(c.fact, x)
    return x
end

struct SingularException <: Exception
    msg::String
end

Base.showerror(io::IO, e::SingularException) = print(io, "SingularException: ", e.msg)
