using PoissonBrackets
using Random
using Test

# Several tests excite the full spectrum with random degrees of freedom, which is
# deliberate: a cross-conservation identity tested on a single smooth mode reports round-off
# and hides the cubic obstruction entirely. The seed is fixed so that a failure is
# reproducible.
Random.seed!(0x5c1e9a3b)

include("exact_tests.jl")
include("spaces_tests.jl")
include("brackets_tests.jl")
include("algebras_tests.jl")
include("dirac_tests.jl")
include("hamiltonians_tests.jl")
include("flows_tests.jl")
include("integrators_tests.jl")
include("kdv_tests.jl")
include("miura_tests.jl")
include("burgers_tests.jl")
include("camassaholm_tests.jl")
include("diagnostics_tests.jl")
