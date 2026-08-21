```@meta
CurrentModule = PoissonBrackets
```

# Integrators

Three properties are in play, and no method has all three.

| method | Poisson map | generating `H` | mass |
|:--|:--|:--|:--|
| [`ImplicitMidpoint`](@ref) | yes, for a constant bracket | only if `H` is quadratic | yes |
| [`AverageVectorField`](@ref) | no | yes, for a constant bracket | yes |
| [`Gonzalez`](@ref), [`GonzalezMass`](@ref) | no | yes, for any `H` | yes |
| [`ExplicitEuler`](@ref), [`RungeKutta4`](@ref) | no | no | yes |
| [`ProjectionMethod`](@ref) | no | yes, by construction | only if asked for |

## The mass is free

Its gradient spans the kernel of the first bracket, ``\mathbb{P}^1 g = 0``, so *any*
increment lying in that bracket's range leaves it untouched — explicit Euler included. This
is the only conservation law here that costs nothing.

## Energy or structure, but not both

That no method here is both a Poisson map and exactly energy-preserving is the
**Ge-Marsden theorem**: such a method would reproduce the exact flow up to a reparametrisation
of time. The theorem's non-degeneracy hypotheses are not checked for these systems, so it is
the reason to expect the trade-off rather than a proof of it here. The trade-off itself is
measured directly, through [`poisson_defect`](@ref).

## Why the tangent map is analytic

[`tangent_map`](@ref) differentiates the step *implicitly* rather than by finite
differences. A central difference has a floor around ``4\times10^{-10}`` on these problems,
while the defect of a genuine Poisson integrator sits at round-off — so a differenced
tangent map would measure its own truncation error and report it as a violation of the
Poisson property, three orders too large.

For implicit midpoint the implicit differentiation gives the Cayley transform

```math
D\Phi_{\Delta t} = \left( \mathbb{I} - \tfrac{1}{2}\Delta t \, \mathbb{P} H'' \right)^{-1}
                   \left( \mathbb{I} + \tfrac{1}{2}\Delta t \, \mathbb{P} H'' \right) ,
```

and a Cayley transform of ``\mathbb{P}`` times a symmetric matrix preserves ``\mathbb{P}``
exactly.

## The nonlinear solve

[`Integrator`](@ref) builds one `SimpleSolvers.NewtonSolver` and reuses it across every
step. Its `refactorize` option is the quasi-Newton scheme these problems want: the Jacobian
is factorised once and reused for several iterations.

The Jacobian supplied is analytic, including for the Gonzalez discrete gradient, whose
rank-one correction is differentiated in closed form.

Two departures from SimpleSolvers' defaults are worth knowing about, both measured rather
than assumed:

  - **No line search.** The residual is a small perturbation of the identity at these step
    sizes and Newton converges in two or three iterations from the previous state, so the
    default `Backtracking` spends several extra residual evaluations per iteration probing a
    step length that is always accepted at one. The default here is `Static`. Pass
    `linesearch = Backtracking(Float64)` if a run is pushed to a step size that needs it.
  - **`LapackLU` for the factorisation.** SimpleSolvers' default is a hand-written scalar
    LU, portable and the right choice for the small dense systems it is usually pointed at.
    At ``N = 384`` it accounted for **74 % of the cost of an implicit step** — about 17 ms
    against LAPACK's 0.6 ms. `SimpleSolvers.LapackLU` keeps the whole nonlinear driver and
    replaces only the factorisation; the test suite checks that the two produce the same
    step. It is restricted to the element types LAPACK provides, so `SimpleSolvers.LU()`
    remains the only option for e.g. `BigFloat`.

Together with the assembly changes of the [Discretisation](@ref) page this took an implicit
step at ``N = 384`` from 45 ms to 4.2 ms.

## Diagnostics

[`integrate`](@ref) evaluates the invariants at **every** step but reports **windowed
maxima**. The error envelope of a geometric integrator oscillates with a period of tens of
steps; sampling it every `stride` steps over a run of millions aliases that oscillation into
noise an order of magnitude below the true envelope, and the resulting figure shows a method
conserving an invariant far better than it does.

Use [`drift`](@ref) for the Hamiltonians and [`absolute_drift`](@ref) for a quantity whose
reference value is zero — the mass of the cosine initial condition, for one.
[`growth`](@ref) separates a bounded oscillation from a secular drift, which a single
end-of-run number cannot.
