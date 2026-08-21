#!/usr/bin/env julia
#
# The numerical experiments of the KdV manuscript: the relative error of BOTH Hamiltonians
# against time, and the initial and final states, for each of the standard cases.
#
#     julia --project=scripts scripts/kdv_energies.jl [case ...]
#
# with cases `cos`, `soliton1`, `soliton2`, `soliton3`, `soliton5` and `miura`. Writes two
# PDFs per case into scripts/figures/. With no arguments it runs all six.
#
# In every case each geometric method holds the Hamiltonian generating its flow to round-off
# and lets the other drift by orders of magnitude more. RK4's error is the only one that
# genuinely grows, and on the narrow examples it also pays for the structure it lacks in
# step count.
#
# The `miura` case is different in kind from the other five and adds a THIRD vector field.
# It is posed in v, not in u, and it has to be: int u_h = int v_h^2, so the image of the
# discrete Miura map lies in C_0,d > 0, and none of the other five -- cos x has zero mass,
# the solitons are depressions -- has a preimage at all.

using CairoMakie
using PoissonBrackets

const FIGDIR = joinpath(@__DIR__, "figures")

const T_END = 100.0          # the four small examples run to the same final time
const RK4_SAFETY = 0.8       # fraction of the RK4 stability limit actually used
const NOUT = 400             # output windows per run

Base.@kwdef struct Case
    key::String
    title::String
    n::Int
    p::Int = 3
    L::Float64
    T::Float64
    dt::Float64
    u0 = nothing             # posed in u ...
    v0 = nothing             # ... or, for the Miura case, in v
    xleft::Float64 = 0.0     # only the plots carry the box back to its own coordinates
    note::String = ""
end

const LSOL = 40.0

const CASES = [
    Case(key = "cos", title = "u₀ = cos x",
         n = 64, L = 2π, T = T_END, dt = 1e-3,
         u0 = cosine(2π)),

    Case(key = "soliton1", title = "single soliton",
         n = 128, L = LSOL, T = T_END, dt = 5e-3,
         u0 = soliton(1.0, 10.0, LSOL), note = "ten full traversals"),

    Case(key = "soliton2", title = "two-soliton",
         n = 128, L = LSOL, T = T_END, dt = 5e-3,
         u0 = two_soliton(1.2, 0.8, 12.0, 22.0, LSOL), note = "eight overtakings"),

    # Shi, Fu and Liu, Appl. Math. Comput. 508 (2026) 129620, Problem 4.2
    Case(key = "soliton3", title = "three solitons",
         n = 256, L = 200.0, T = 250.0, dt = 0.1,
         u0 = three_solitons(), xleft = -100.0,
         note = "the tall wave overtakes both"),

    # ... and Problem 4.3
    Case(key = "soliton5", title = "five solitons",
         n = 384, L = 300.0, T = 500.0, dt = 1/16,
         u0 = five_solitons(), xleft = -150.0,
         note = "the tall wave overtakes all four"),

    # posed in v; u₀ = M_h(v₀)
    Case(key = "miura", title = "Miura, v₀ = 1 + sin x",
         n = 64, L = 2π, T = T_END, dt = 1e-3,
         v0 = miura_initial_v(2π), note = "u₀ = M_h(v₀)"),
]

find_case(k) = (i = findfirst(c -> c.key == k, CASES);
                i === nothing ? error("unknown case $(k); known: " *
                                      join((c.key for c in CASES), ", ")) : CASES[i])

"""The six runs posed in u: colour carries the integrator, line style the vector field."""
runs_u(sys) = (("midpoint, flow 1",          sys.flow1, ImplicitMidpoint(),   false),
               ("midpoint, flow 2",          sys.flow2, ImplicitMidpoint(),   false),
               ("discrete gradient, flow 1", sys.flow1, Gonzalez(),           false),
               ("AVF, flow 1",               sys.flow1, AverageVectorField(), false),
               ("RK4, flow 1",               sys.flow1, RungeKutta4(),        true),
               ("RK4, flow 2",               sys.flow2, RungeKutta4(),        true))

"""The two Miura runs, carried in v. They cannot be added to `runs_u`: a MiuraSystem holds
v̂ as its state, so the two families never appear on the same system object — only on the
same figure, through the same initial u_h."""
runs_v(sysM) = (("midpoint, Miura",          sysM.flow, ImplicitMidpoint(), false),
                ("discrete gradient, Miura", sysM.flow, Gonzalez(),         false))

function run_one(sys, flow, meth, x0, dt, T)
    NT = max(round(Int, T / dt), NOUT)
    stride = max(NT ÷ NOUT, 1)
    traj = integrate(sys, Integrator(flow, meth, dt), x0, NT; stride = stride)
    (traj, NT)
end

report(label, traj) = @info(
    "  " * rpad(label, 26) *
    " |ΔH₁| = $(round(drift(traj, :H1), sigdigits = 2))" *
    "  |ΔH₂| = $(round(drift(traj, :H2), sigdigits = 2))" *
    "  |ΔC₀| = $(round(absolute_drift(traj, :C0), sigdigits = 2))" *
    "  g(H₁) = $(round(growth(traj, :H1), sigdigits = 2))")

function run_case(case::Case)
    @info "running $(case.key): $(case.title)" * (isempty(case.note) ? "" : "  ($(case.note))")
    s   = SplineSpace(case.n, case.p; L = case.L)
    sys = KdVSystem(s)

    # the Miura case is posed in v, and the u system is then seeded from u₀ = M_h(v₀) so
    # that all eight runs start from the same field
    sysM = case.v0 === nothing ? nothing : MiuraSystem(s)
    v0   = case.v0 === nothing ? nothing : project(s, case.v0)
    u0   = case.v0 === nothing ? project(s, case.u0) : miura_map(s, v0)

    dt_rk4 = RK4_SAFETY * min(PoissonBrackets.stability_limit(sys.flow1, u0),
                              PoissonBrackets.stability_limit(sys.flow2, u0))

    trajs, labels, finals = Trajectory[], String[], Vector{Float64}[]

    for (label, flow, meth, explicit) in runs_u(sys)
        dt = explicit ? min(dt_rk4, case.dt) : case.dt
        traj, _ = run_one(sys, flow, meth, u0, dt, case.T)
        push!(trajs, traj); push!(labels, label); push!(finals, traj.final)
        report(label, traj)
    end

    if sysM !== nothing
        for (label, flow, meth, _) in runs_v(sysM)
            traj, _ = run_one(sysM, flow, meth, v0, case.dt, case.T)
            push!(trajs, traj); push!(labels, label)
            # the state of a Miura run is v̂; map it before plotting alongside the others
            push!(finals, miura_map(s, traj.final))
            report(label, traj)
        end
    end

    mkpath(FIGDIR)
    save(joinpath(FIGDIR, "kdv-energy-$(case.key).pdf"),
         energyplot(trajs, labels; names = (:H1, :H2), title = case.title))
    save(joinpath(FIGDIR, "kdv-state-$(case.key).pdf"),
         stateplot(sys, vcat([u0], finals), vcat(["initial"], labels);
                   xleft = case.xleft,
                   title = "$(case.title): initial and final states"))
    @info "  wrote kdv-energy-$(case.key).pdf and kdv-state-$(case.key).pdf"
end

for k in (isempty(ARGS) ? [c.key for c in CASES] : ARGS)
    run_case(find_case(k))
end
