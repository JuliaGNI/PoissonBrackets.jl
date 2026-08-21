using PoissonBrackets
using LinearAlgebra
using SimpleSplines: UniformMesh
using Test

@testset "$(rpad("Integrator Tests",80))" begin

    s   = SplineSpace(UniformMesh(20, 2π), 3)
    sys = KdVSystem(s)
    u0  = project(s, cosine(2π))

    @testset "$(rpad("one step of every method runs and stays finite",76))" begin
        for meth in (ExplicitEuler(), RungeKutta4(), ImplicitMidpoint(),
                     AverageVectorField(), Gonzalez(), GonzalezMass())
            for flow in (sys.flow1, sys.flow2)
                integ = Integrator(flow, meth, 1e-4)
                û = copy(u0)
                integrate_step!(û, integ)
                @test all(isfinite, û)
                @test û != u0
            end
        end
    end

    @testset "$(rpad("implicit midpoint reduces to the fixed point it defines",76))" begin
        integ = Integrator(sys.flow1, ImplicitMidpoint(), 1e-3)
        û = copy(u0)
        integrate_step!(û, integ)
        r = û .- u0 .- 1e-3 .* vectorfield(sys.flow1, (u0 .+ û) ./ 2)
        @test maximum(abs, r) < 1e-11
    end

    @testset "$(rpad("the discrete gradient really is one",76))" begin
        # the defining property: ḡ ⋅ Δû = ΔH exactly, which is what makes energy
        # conservation follow from antisymmetry alone
        for m in (Gonzalez(), GonzalezMass())
            x = copy(u0)
            y = u0 .+ 0.05 .* randn(length(u0))
            ḡ = PoissonBrackets.discrete_gradient(m, sys.flow1, x, y)
            ΔH = hamiltonian(sys.flow1, y) - hamiltonian(sys.flow1, x)
            @test dot(ḡ, y .- x) ≈ ΔH atol = 1e-12 * max(1, abs(ΔH))
        end
    end

    @testset "$(rpad("the rank-one correction vanishes for a quadratic Hamiltonian",76))" begin
        # which is where the name "midpoint discrete gradient" comes from: for a quadratic
        # Hamiltonian the method IS the plain midpoint rule
        x = copy(u0)
        y = u0 .+ 0.05 .* randn(length(u0))
        for m in (Gonzalez(), GonzalezMass())
            ḡ = PoissonBrackets.discrete_gradient(m, sys.flow2, x, y)
            @test ḡ ≈ gradient(sys.flow2, (x .+ y) ./ 2) atol = 1e-10
        end
    end

    @testset "$(rpad("discrete gradient Jacobian against finite differences",76))" begin
        x = copy(u0)
        y = u0 .+ 0.05 .* randn(length(u0))
        for m in (Gonzalez(), GonzalezMass())
            J = PoissonBrackets.discrete_gradient_jacobian(m, sys.flow1, x, y)
            h = 1e-6
            for k in (1, 7, 15)
                e = zeros(length(y)); e[k] = h
                fd = (PoissonBrackets.discrete_gradient(m, sys.flow1, x, y .+ e) .-
                      PoissonBrackets.discrete_gradient(m, sys.flow1, x, y .- e)) ./ 2h
                @test J[:, k] ≈ fd atol = 1e-4 * max(1, maximum(abs, fd))
            end
        end
    end

    @testset "$(rpad("tangent map against finite differences",76))" begin
        for meth in (ImplicitMidpoint(), AverageVectorField(), Gonzalez(), ExplicitEuler())
            integ = Integrator(sys.flow1, meth, 1e-3)
            DΦ = tangent_map(integ, copy(u0))
            h = 1e-6
            for k in (1, 9, 17)
                e = zeros(length(u0)); e[k] = h
                yp = u0 .+ e; ym = u0 .- e
                integrate_step!(yp, integ); integrate_step!(ym, integ)
                @test DΦ[:, k] ≈ (yp .- ym) ./ 2h atol = 1e-5
            end
        end
    end

    @testset "$(rpad("implicit midpoint IS a Poisson map, AVF is not",76))" begin
        # P¹ is constant, so implicit midpoint is the Cayley transform of P¹ times a
        # symmetric matrix, which preserves P¹ exactly. Ge-Marsden then says no method can
        # be both this and exactly energy-preserving.
        mid = Integrator(sys.flow1, ImplicitMidpoint(), 2e-3)
        avf = Integrator(sys.flow1, AverageVectorField(), 2e-3)
        dgr = Integrator(sys.flow1, Gonzalez(), 2e-3)
        eul = Integrator(sys.flow1, ExplicitEuler(), 2e-3)
        @test poisson_defect(mid, copy(u0)) < 1e-13
        @test poisson_defect(avf, copy(u0)) > 1e-7
        @test poisson_defect(dgr, copy(u0)) > 1e-9
        @test poisson_defect(eul, copy(u0)) > 1e-2
    end

    @testset "$(rpad("stability limit and explicit stability",76))" begin
        dtmax = PoissonBrackets.stability_limit(sys.flow1, copy(u0))
        @test dtmax > 0
        # stable at 40% of the limit, unstable at 160%
        for (factor, stable) in ((0.4, true), (1.6, false))
            integ = Integrator(sys.flow1, RungeKutta4(), factor * dtmax)
            û = copy(u0)
            for _ in 1:400
                integrate_step!(û, integ)
                all(isfinite, û) || break
            end
            grew = !all(isfinite, û) || maximum(abs, û) > 100 * maximum(abs, u0)
            @test grew == !stable
        end
    end

    @testset "$(rpad("projection method holds several invariants at once",76))" begin
        # projecting onto both Hamiltonians alone LOSES the Casimir, because the correction
        # direction is not mass-neutral; a third multiplier along g restores all three
        both  = ProjectionMethod(ImplicitMidpoint(), (sys.H1, sys.H2))
        three = ProjectionMethod(ImplicitMidpoint(), (sys.H1, sys.H2, sys.C0))

        t2 = integrate(sys, Integrator(sys.flow1, both, 2e-3), u0, 100; stride = 10)
        @test drift(t2, :H1) < 1e-12
        @test drift(t2, :H2) < 1e-12

        t3 = integrate(sys, Integrator(sys.flow1, three, 2e-3), u0, 100; stride = 10)
        @test drift(t3, :H1) < 1e-12
        @test drift(t3, :H2) < 1e-12
        @test absolute_drift(t3, :C0) < 1e-12
    end

    @testset "$(rpad("the LAPACK linear solver agrees with the SimpleSolvers one",76))" begin
        # LapackLU replaces only the factorisation inside SimpleSolvers' Newton loop, so
        # the two must produce the same step -- it is a performance substitution, not a
        # different method.
        import SimpleSolvers
        for meth in (ImplicitMidpoint(), AverageVectorField(), Gonzalez())
            a = copy(u0); b = copy(u0)
            integrate_step!(a, Integrator(sys.flow1, meth, 1e-3;
                                          linear_solver_method = LapackLU()))
            integrate_step!(b, Integrator(sys.flow1, meth, 1e-3;
                                          linear_solver_method = SimpleSolvers.LU()))
            @test a ≈ b atol = 1e-11
        end
        # and over a run
        ta = integrate(sys, Integrator(sys.flow1, ImplicitMidpoint(), 1e-3;
                                       linear_solver_method = LapackLU()), u0, 100; stride = 10)
        tb = integrate(sys, Integrator(sys.flow1, ImplicitMidpoint(), 1e-3;
                                       linear_solver_method = SimpleSolvers.LU()), u0, 100; stride = 10)
        @test ta.final ≈ tb.final atol = 1e-10
        @test drift(ta, :H1) ≈ drift(tb, :H1) rtol = 1e-6
    end

    @testset "$(rpad("solver is built once and reused across steps",76))" begin
        integ = Integrator(sys.flow1, ImplicitMidpoint(), 1e-3)
        solver = integ.solver
        û = copy(u0)
        for _ in 1:5
            integrate_step!(û, integ)
        end
        @test integ.solver === solver
    end

end
