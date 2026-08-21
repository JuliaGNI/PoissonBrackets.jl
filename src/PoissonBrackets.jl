module PoissonBrackets

    using CompactBasisFunctions
    using LinearAlgebra
    using QuadratureRules
    using SimpleSolvers
    using SimpleSplines

    # extended rather than defined here, so that the packages of the ecosystem share one
    # generic function per accessor
    import GeometricBase: integrate, value
    import SimpleSplines: basis, basis_integrals, basis_values, degree, domainlength,
                          evaluate, l2_projection, mass_factorization, mass_matrix, nbasis,
                          ncells, nodes, order, quadrature_nodes, quadrature_weights

    export DiscreteSpace, SplineSpace, LagrangeSpace
    export basis, nbasis, degree, order, nodes, ncells, domainlength,
           mass_matrix, mass_factorization, basis_values, basis_integrals,
           quadrature_nodes, quadrature_weights, project, evaluate,
           derivative_matrix, stiffness_matrix, mixed_matrix, weighted_matrix

    include("spaces.jl")

    export DiscreteBracket, ConstantBracket, AffineBracket, GaugedBracket, MiuraBracket
    export poisson_matrix, poisson_apply, poisson_tensor, isantisymmetric,
           jacobi_residual, structure_constant_residual, casimir_gradient

    include("brackets.jl")

    export DiscreteHamiltonian, hamiltonian, gradient, hessian
    export MassCasimir, QuadraticHamiltonian

    include("hamiltonians.jl")

    export HamiltonianFlow, vectorfield, vectorfield!, jacobian

    include("flows.jl")

    export IntegratorMethod, ExplicitEuler, RungeKutta4, ImplicitMidpoint,
           AverageVectorField, DiscreteGradient, Gonzalez, GonzalezMass,
           ProjectionMethod
    export Integrator, integrate_step!, tangent_map

    include("integrators.jl")

    export KdVSystem, kdv_bracket_1, kdv_bracket_2, kdv_miura_bracket,
           KdVHamiltonian1, KdVHamiltonian2, KdVHamiltonian3,
           MiuraSystem, MiuraHamiltonian, miura_map, miura_invert, miura_bracket,
           miura_casimir, miura_derivative, miura_moment_matrix, hill_lambda0,
           soliton, two_soliton, solitons, cosine,
           three_solitons, five_solitons, miura_initial_v

    include("equations/kdv.jl")

    export BurgersSystem, burgers_bracket, BurgersHamiltonian, burgers_casimir

    include("equations/burgers.jl")

    export CamassaHolmSystem, camassa_holm_bracket_1, camassa_holm_bracket_2,
           CamassaHolmHamiltonian1, CamassaHolmHamiltonian2, HelmholtzMap,
           velocity, momentum, peakon

    include("equations/camassaholm.jl")

    export invariants, invariant_names, integrate, Trajectory, deviation,
           drift, absolute_drift, growth, poisson_defect
    export energyplot, stateplot, INTEGRATOR_COLORS, FLOW_STYLES, ERROR_FLOOR

    include("diagnostics.jl")

end
