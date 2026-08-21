module PoissonBracketsMakieExt

using CairoMakie
using PoissonBrackets
using PoissonBrackets: ERROR_FLOOR, FLOW_STYLES, INTEGRATOR_COLORS

"""
    _series_style(label)

Resolve a run label of the form `"<integrator>, flow <n>"` into a colour and a line style.

Colour carries the integrator and line style the vector field, so that two runs of the same
method on different flows are told apart by shape rather than by hue, and no series in the
figure is identified by colour alone.
"""
function _series_style(label::AbstractString)
    parts = split(label, ',')
    name  = strip(parts[1])
    color = get(INTEGRATOR_COLORS, name, "#666666")
    flow = 1
    for p in parts[2:end]
        q = strip(p)
        occursin("Miura", q) && (flow = :M; continue)
        m = match(r"flow\s*(\d+)", q)
        m === nothing || (flow = parse(Int, m.captures[1]))
    end
    (color = color, linestyle = get(FLOW_STYLES, flow, :solid))
end

_clamped(v) = max.(v, ERROR_FLOOR)

function PoissonBrackets.energyplot(trajs, labels;
                                    name::Symbol = :H1,
                                    title = "",
                                    ylabel = nothing,
                                    size = (620, 380),
                                    resolution = nothing)
    isempty(trajs) && throw(ArgumentError("no trajectories to plot"))
    length(trajs) == length(labels) || throw(DimensionMismatch(
        "got $(length(trajs)) trajectories but $(length(labels)) labels"))

    # A reference value that vanishes -- the mass of the cosine initial condition is zero to
    # round-off -- makes a relative error meaningless, so the whole panel switches to the
    # absolute one rather than reporting a 100-fold "drift" of a conserved quantity.
    refs = [abs(getproperty(t.reference, name)) for t in trajs]
    relative = all(r -> r > sqrt(eps(Float64)), refs)

    fig = Figure(size = resolution === nothing ? size : resolution)
    ax = Axis(fig[1, 1];
              xlabel = "t",
              ylabel = ylabel === nothing ?
                       (relative ? "relative error" : "absolute error") : ylabel,
              title = isempty(title) ? string(name) : "$(title):  $(name)",
              yscale = log10)

    for (traj, label) in zip(trajs, labels)
        st = _series_style(label)
        d = deviation(traj, name)
        y = relative ? d ./ abs(getproperty(traj.reference, name)) : d
        lines!(ax, traj.t, _clamped(y);
               color = st.color, linestyle = st.linestyle, label = label)
    end

    Legend(fig[2, 1], ax; orientation = :horizontal, nbanks = 2, framevisible = false)
    fig
end

function PoissonBrackets.stateplot(system, states, labels;
                                   npoints = 601,
                                   title = "",
                                   xleft = 0.0,
                                   size = (900, 400),
                                   resolution = nothing)
    length(states) == length(labels) || throw(DimensionMismatch(
        "got $(length(states)) states but $(length(labels)) labels"))

    s = system.space
    L = domainlength(s)
    xs = collect(range(0, L, length = npoints))

    fig = Figure(size = resolution === nothing ? size : resolution)
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "u")

    # Nested strokes of decreasing width, so that curves lying on top of one another all
    # stay visible instead of the last one drawn hiding the rest. The widths are spread
    # over the whole range for however many runs there are -- zipping against a fixed
    # tuple would silently truncate, and the extra runs of the Miura case would vanish
    # from the figure without a word.
    widths = length(states) == 1 ? [2.0] :
             collect(range(2.8, 0.6, length = length(states)))
    for (i, (û, label)) in enumerate(zip(states, labels))
        st = _series_style(label)
        lines!(ax, xs .+ xleft, evaluate(s, û, xs);
               color = st.color, linestyle = st.linestyle,
               linewidth = widths[i], label = label)
    end

    # The legend goes BESIDE the axes, not inside them: the final states of these runs lie
    # almost on top of one another and fill the frame, so an inset legend covers the very
    # part of the curve the figure is about.
    Legend(fig[1, 2], ax; framevisible = false)
    isempty(title) || Label(fig[0, 1:2], title; fontsize = 16)
    fig
end

end
