using PoissonBrackets
using LinearAlgebra
using SimpleSplines: UniformMesh, RandomMesh
using Test

@testset "$(rpad("Miura Tests",80))" begin

    @testset "$(rpad("the mass identity int u_h = int v_h^2 is EXACT",76))" begin
        # Pairing û = M_h(v̂) with the partition of unity gives this to the last bit, and it
        # is what confines the image of the map to positive mass.
        for (nm, mk) in SPLINE_MESHES, N in (16, 24, 32)
            s = SplineSpace(mk(N), 3)
            v̂ = project(s, miura_initial_v(2π))
            û = miura_map(s, v̂)
            lhs = dot(basis_integrals(s), û)
            rhs = dot(quadrature_weights(s), PoissonBrackets.field(s, v̂) .^ 2)
            @test lhs ≈ rhs atol = 1e-12
            @test lhs > 0
        end
    end

    @testset "$(rpad("Miura bracket is antisymmetric and Poisson at every N",76))" begin
        # The Jacobiator is COMPUTED here, through the chart -- it is not zero by fiat.
        # Against ≈ 0.42, flat, for the Galerkin bracket.
        for N in (12, 16, 24)
            s = SplineSpace(UniformMesh(N, 2π), 3)
            v̂ = project(s, miura_initial_v(2π))
            b = kdv_miura_bracket(s, v̂)
            û = miura_map(s, v̂)
            P = poisson_matrix(b, û)
            @test maximum(abs, P + P') < 1e-10 * maximum(abs, P)
            @test jacobi_residual(b, û) < 1e-12
            @test jacobi_residual(kdv_bracket_2(s), û) > 0.3
        end
    end

    @testset "$(rpad("the derivative tensor against finite differences",76))" begin
        # the chain rule through the chart, ∂P/∂û = (∂P/∂v̂) L⁻¹
        N = 12
        s = SplineSpace(UniformMesh(N, 2π), 3)
        v̂ = project(s, miura_initial_v(2π))
        b = kdv_miura_bracket(s, v̂)
        û = miura_map(s, v̂)
        dP = PoissonBrackets.poisson_derivative(b, û)
        L = miura_derivative(s, v̂)
        h = 1e-6
        for m in (1, 5, 11)
            e = zeros(N); e[m] = h
            vp = v̂ .+ L \ e
            vm = v̂ .- L \ e
            fd = (poisson_matrix(kdv_miura_bracket(s, vp), û) .-
                  poisson_matrix(kdv_miura_bracket(s, vm), û)) ./ 2h
            @test dP[m, :, :] ≈ fd atol = 1e-4 * max(1, maximum(abs, fd))
        end
    end

    @testset "$(rpad("K_Miura is DENSE where the Galerkin assembly is banded",76))" begin
        # the two interior inverse mass matrices are the nonlocality that carries the
        # Jacobi identity; the comparison must be on the weak-form blocks K = M P M
        for N in (16, 24)
            s = SplineSpace(UniformMesh(N, 2π), 3)
            v̂ = project(s, miura_initial_v(2π))
            M = mass_matrix(s)
            û = miura_map(s, v̂)
            Km = M * poisson_matrix(kdv_miura_bracket(s, v̂), û) * M
            Kg = M * poisson_matrix(kdv_bracket_2(s), û) * M
            frac(A) = count(x -> abs(x) > 1e-8 * maximum(abs, A), A) / length(A)
            @test frac(Km) > 0.8
            @test frac(Kg) < 0.6
        end
    end

    @testset "$(rpad("the Hill criterion is sharp",76))" begin
        # A preimage exists exactly when the lowest eigenvalue of -∂ₓ² + u_h is positive.
        # The threshold for u = c + sin x + 0.4 cos 2x at N = 20, p = 3 is c = 0.470894.
        s = SplineSpace(UniformMesh(20, 2π), 3)
        u(c) = project(s, x -> c + sin(x) + 0.4cos(2x))
        for c in (0.3, 0.45, 0.47)
            @test hill_lambda0(s, u(c)) < 0
            @test miura_invert(s, u(c)) === nothing        # no preimage at all
        end
        for c in (0.48, 0.6, 1.0)
            @test hill_lambda0(s, u(c)) > 0
            @test miura_invert(s, u(c)) !== nothing
        end
        # the threshold itself, to four digits
        @test abs(hill_lambda0(s, u(0.470894))) < 1e-5
    end

    @testset "$(rpad("miura_invert round-trips, and the branches differ",76))" begin
        s = SplineSpace(UniformMesh(20, 2π), 3)
        û = project(s, x -> 1.0 + sin(x) + 0.4cos(2x))
        v̂ = miura_invert(s, û)
        @test v̂ !== nothing
        @test miura_map(s, v̂) ≈ û atol = 1e-12
        # where there is a preimage there are two, with opposite ∫v
        w = miura_invert(s, û; seed = -1.0)
        @test w !== nothing
        @test miura_map(s, w) ≈ û atol = 1e-12
        @test sign(dot(basis_integrals(s), v̂)) == -sign(dot(basis_integrals(s), w))
    end

    @testset "$(rpad("none of the standard examples has a preimage",76))" begin
        # cos x has zero mass and the solitons are depressions, so all have C₀ ≤ 0
        s = SplineSpace(UniformMesh(32, 2π), 3)
        for f in (cosine(2π), x -> sin(x) + 0.4cos(2x))
            û = project(s, f)
            @test dot(basis_integrals(s), û) < 1e-10
            @test miura_invert(s, û) === nothing
        end
        sl = SplineSpace(UniformMesh(64, 40.0), 3)
        û = project(sl, soliton(1.0, 10.0, 40.0))
        @test dot(basis_integrals(sl), û) < 0
        @test miura_invert(sl, û) === nothing
    end

    @testset "$(rpad("MiuraHamiltonian gradient and Hessian",76))" begin
        s = SplineSpace(UniformMesh(12, 2π), 3)
        H = MiuraHamiltonian(s)
        v̂ = project(s, miura_initial_v(2π))
        h = 1e-6
        g = gradient(H, s, v̂)
        fd = [(hamiltonian(H, s, v̂ .+ h .* e) - hamiltonian(H, s, v̂ .- h .* e)) / 2h
              for e in eachcol(Matrix(I, 12, 12))]
        @test g ≈ fd atol = 1e-5 * max(1, maximum(abs, g))
        A = hessian(H, s, v̂)
        @test A ≈ A'
        fdH = zeros(12, 12)
        for m in 1:12
            e = zeros(12); e[m] = h
            fdH[:, m] = (gradient(H, s, v̂ .+ e) .- gradient(H, s, v̂ .- e)) ./ 2h
        end
        @test A ≈ fdH atol = 1e-4 * maximum(abs, A)
        # H̃ is the same number as H₂ at the corresponding û -- the value of a function is
        # chart-free, its algebraic form is not
        @test hamiltonian(H, s, v̂) ≈ hamiltonian(KdVHamiltonian2(s), s, miura_map(s, v̂))
    end

    @testset "$(rpad("MiuraSystem reports in u",76))" begin
        s = SplineSpace(UniformMesh(24, 2π), 3)
        sys = MiuraSystem(s)
        v̂ = project(s, miura_initial_v(2π))
        û = miura_map(s, v̂)
        kdv = KdVSystem(s)
        @test invariants(sys, v̂).H1 ≈ invariants(kdv, û).H1
        @test invariants(sys, v̂).H2 ≈ invariants(kdv, û).H2
        @test invariants(sys, v̂).C0 ≈ invariants(kdv, û).C0
        @test invariants(sys, v̂).C0 > 0
        @test invariant_names(sys) == (:H1, :H2, :C0)
        @test nbasis(sys) == 24
    end

    @testset "$(rpad("the Miura flow reproduces the KdV flow of the second bracket",76))" begin
        for N in (32, 64)
            s = SplineSpace(UniformMesh(N, 2π), 3)
            sys = MiuraSystem(s)
            v̂ = project(s, miura_initial_v(2π))
            û = miura_map(s, v̂)
            # v̇ = P¹ ∂H̃/∂v̂ pushes forward to u̇ = P²_M ∂H₂/∂û
            v̇ = vectorfield(sys.flow, v̂)
            pushed = miura_derivative(s, v̂) * v̇
            direct = poisson_apply(kdv_miura_bracket(s, v̂), û, gradient(KdVHamiltonian2(s), s, û))
            @test pushed ≈ direct atol = 1e-9 * maximum(abs, direct)
        end
    end

    @testset "$(rpad("CHART: midpoint loses H2 in v, a discrete gradient does not",76))" begin
        # H₂ is quadratic in û and QUARTIC in v̂, so the theorem that a symplectic
        # Runge-Kutta method conserves quadratic invariants does not apply in this chart.
        # What is chart-free is the VALUE of a function, so a discrete gradient method,
        # whose property is ḡ ⋅ Δv̂ = ΔH̃ by construction, does conserve it.
        s   = SplineSpace(UniformMesh(20, 2π), 3)
        sys = MiuraSystem(s)
        v0  = project(s, miura_initial_v(2π))
        dt  = 1.6e-3
        NT  = round(Int, 0.2 / dt)

        mid = integrate(sys, Integrator(sys.flow, ImplicitMidpoint(), dt), v0, NT; stride = 5)
        @test drift(mid, :H2) > 1e-6            # ~7e-5
        @test drift(mid, :H2) < 1e-3

        dg = integrate(sys, Integrator(sys.flow, Gonzalez(), dt), v0, NT; stride = 5)
        @test drift(dg, :H2) < 1e-11            # exact
    end

    @testset "$(rpad("the Poisson property transports through the chart",76))" begin
        # This is the first genuine Poisson integrator for a discrete SECOND KdV structure;
        # the Galerkin bracket admits none.
        s   = SplineSpace(UniformMesh(20, 2π), 3)
        sys = MiuraSystem(s)
        v0  = project(s, miura_initial_v(2π))
        dt  = 1.6e-3

        integ = Integrator(sys.flow, ImplicitMidpoint(), dt)
        v1 = copy(v0); integrate_step!(v1, integ)
        DΦv = tangent_map(integ, copy(v0))

        P1 = poisson_matrix(sys.bracket)
        @test maximum(abs, DΦv * P1 * DΦv' - P1) / maximum(abs, P1) < 1e-13

        # pushed forward, DΦᵤ = L(v¹) DΦᵥ L(v⁰)⁻¹ must map P²_M(v⁰) to P²_M(v¹)
        L0 = miura_derivative(s, v0); L1 = miura_derivative(s, v1)
        DΦu = L1 * DΦv * inv(L0)
        PM0 = poisson_matrix(kdv_miura_bracket(s, v0), v0)
        PM1 = poisson_matrix(kdv_miura_bracket(s, v1), v1)
        @test maximum(abs, DΦu * PM0 * DΦu' - PM1) / maximum(abs, PM1) < 1e-12
    end

    @testset "$(rpad("miura_casimir is the Casimir of P1 in the v chart",76))" begin
        s = SplineSpace(UniformMesh(20, 2π), 3)
        sys = MiuraSystem(s)
        v̂ = project(s, miura_initial_v(2π))
        @test miura_casimir(sys, v̂) ≈ dot(basis_integrals(s), v̂)
        # ∫v is conserved by the flow; ∫u = ∫v² is not a Casimir here
        @test abs(dot(basis_integrals(s), vectorfield(sys.flow, v̂))) < 1e-10
    end

    @testset "$(rpad("new initial conditions",76))" begin
        u3 = three_solitons()
        @test u3(-60.0 - (-100.0)) ≈ -2 * 0.3^2 atol = 1e-3
        u5 = five_solitons()
        @test u5(-120.0 - (-150.0)) ≈ -2 * 0.3^2 atol = 1e-3
        v0 = miura_initial_v(2π)
        @test v0(0.0) ≈ 1
        @test v0(π/2) ≈ 2
    end

end
