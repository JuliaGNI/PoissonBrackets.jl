```@meta
CurrentModule = PoissonBrackets
```

# PoissonBrackets.jl

Structure-preserving discretisations of the Poisson brackets of one-dimensional Hamiltonian
partial differential equations.

The package keeps the discrete bracket, the discrete Hamiltonian and the time integrator as
three separate, composable objects, so that which structure survives a discretisation can be
asked and answered one property at a time.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGNI/PoissonBrackets.jl")
```

## The pieces

| object | what it is |
|:--|:--|
| [`DiscreteSpace`](@ref) | [`SplineSpace`](@ref) or [`LagrangeSpace`](@ref), with its quadrature |
| [`DiscreteBracket`](@ref) | ``\mathbb{P}(\hat{u})``, the discrete structure matrix |
| [`DiscreteHamiltonian`](@ref) | a functional evaluated on the discrete field, with analytic gradient |
| [`HamiltonianFlow`](@ref) | the two contracted: ``\dot{\hat{u}} = \mathbb{P}(\hat{u}) \, \partial H/\partial \hat{u}`` |
| [`Integrator`](@ref) | a one-step method with its nonlinear solver built once |

## A first run

```@example intro
using PoissonBrackets

sys = KdVSystem(SplineSpace(32, 3))
û₀  = project(sys.space, cosine(2π))

integ = Integrator(sys.flow1, ImplicitMidpoint(), 1e-3)
traj  = integrate(sys, integ, û₀, 2000; stride = 20)

(H1 = drift(traj, :H1), H2 = drift(traj, :H2), C0 = absolute_drift(traj, :C0))
```

Implicit midpoint on the first flow is a **Poisson map** but does not conserve the cubic
``H_1``. Swapping the method for an energy-preserving one trades one for the other:

```@example intro
avf = Integrator(sys.flow1, AverageVectorField(), 1e-3)
t2  = integrate(sys, avf, û₀, 2000; stride = 20)

(H1 = drift(t2, :H1), poisson_defect_midpoint = poisson_defect(integ, copy(û₀)),
 poisson_defect_avf = poisson_defect(avf, copy(û₀)))
```

That no method achieves both is the Ge-Marsden theorem, not a gap in the list: a Poisson
integrator that also conserved the Hamiltonian exactly would reproduce the exact flow up to
a reparametrisation of time.

## What is preserved by what

| | Poisson map | generating `H` | mass |
|:--|:--|:--|:--|
| [`ImplicitMidpoint`](@ref) | yes, for a constant bracket | only if `H` is quadratic | yes |
| [`AverageVectorField`](@ref) | no | yes, for a constant bracket | yes |
| [`Gonzalez`](@ref) | no | yes, for any `H` | yes |
| [`ExplicitEuler`](@ref), [`RungeKutta4`](@ref) | no | no | yes |

The mass comes free to *every* method, explicit Euler included, because its gradient spans
the kernel of the first bracket. See [`MassCasimir`](@ref).

## Plotting and the experiments

Load `CairoMakie` to get [`energyplot`](@ref) and [`stateplot`](@ref).

`scripts/kdv.jl` runs the full set of experiments:

```sh
julia --project=scripts scripts/kdv.jl          # all six cases
julia --project=scripts scripts/kdv.jl cos      # or just one
```

Per case it writes into `scripts/figures/`:

| file | contents |
|:--|:--|
| `kdv-<case>-{H1,H2,C0}-flow1.pdf` | the four flow-1 runs |
| `kdv-<case>-{H1,H2,C0}-other.pdf` | the two flow-2 runs, and the two Miura runs where the case has them |
| `kdv-<case>-state-{flow1,other}.pdf` | initial and final states, split the same way |
| `<case>.md` | a table of the maximum error in each invariant, and what it says |

One invariant per figure and one family of vector fields per figure: the three invariants'
errors sit orders of magnitude apart, so a shared axis flattens them, and the comparison that
matters is *within* a family rather than across.

Only the `miura` case has Miura runs, and it is the only one that can: the image of the
discrete Miura map lies in ``C_{0,d} > 0``, and the other five all have ``C_{0,d} \le 0`` —
``\cos x`` has zero mass and the solitons of this sign convention are depressions. Its
initial data are therefore posed in ``v``. See [`MiuraSystem`](@ref).
