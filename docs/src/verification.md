# Verification

This page records what the package has been checked against, claim by claim, and what came
out of it. The two references are the manuscripts —
`discrete-kdv-brackets{,-notes}.tex` and `discrete-lie-poisson-brackets{,-notes}.tex` — and
the Python prototypes in their `Scripts/` directories, which the Julia is a port of.

It is written so that the *absence* of an entry is informative: everything examined is listed,
including the things that turned out to be correct and the one system that has nothing to be
checked against.

## Method

  - Python: `./run_all.sh` in both `Scripts/` directories (nine verification scripts for KdV,
    ten for Lie-Poisson; each exits nonzero on a failed check), plus
    `plot_kdv_energies.py cos soliton1` for the drift tables, which exist only on stdout.
  - Julia: the test suite, then `julia --project=scripts scripts/kdv.jl cos soliton1`.
  - Both manuscripts recompile with XeTeX — not pdfTeX or LuaTeX; `fontspec` rules out the
    first and a `\boundary` clash with a LuaTeX primitive in `definitions_math.tex` the
    second.

## The Runge-Kutta 4 question

The RK4 results for the `cos` and single-soliton runs were suspected of disagreeing with the
Python. **They do not.** Every quantity above round-off agrees; the columns that differ are
round-off-dominated, and the mass column is the built-in control that shows it.

The step size is identical to full precision, which rules out the obvious explanation first:

| case | ``\Delta t_{\mathrm{RK4}}`` (Python) | steps | Julia |
|:--|:--|--:|:--|
| `cos` | `1.3041539392303622e-4` | 766781 | identical |
| `soliton1` | `4.195777470400295e-3` | 23833 | identical |

The drift tables, Python against Julia:

| `cos` | ``\lvert\Delta H_1\rvert`` | ``\lvert\Delta H_2\rvert`` | ``\lvert\Delta C_0\rvert`` |
|:--|:--|:--|:--|
| RK4, flow 1 | 1.1e-11 / 6.2e-12 | 1.1e-10 / **1.1e-10** | 7.7e-12 / 4.7e-12 |
| RK4, flow 2 | 3.4e-07 / **3.4e-07** | 1.3e-13 / 2.3e-13 | 2.2e-13 / 6.3e-13 |

| `soliton1` | ``\lvert\Delta H_1\rvert`` | ``\lvert\Delta H_2\rvert`` | ``\lvert\Delta C_0\rvert`` |
|:--|:--|:--|:--|
| RK4, flow 1 | 5.9e-08 / **5.9e-08** | 3.7e-08 / **3.6e-08** | 3.9e-12 / 6.5e-12 |
| RK4, flow 2 | 4.4e-07 / **4.4e-07** | 3.8e-08 / **3.8e-08** | 1.3e-14 / 1.5e-14 |

The bold entries are the genuine dynamical drifts, and they agree to every digit the tables
print. What differs — `cos` flow 1's ``\lvert\Delta H_1\rvert``, and the mass column
throughout — sits at ``10^{-11}`` or below.

That those are round-off is not an assumption. ``C_0`` is conserved *exactly* in exact
arithmetic, by [`MassCasimir`](@ref), so its measured drift is a pure meter of accumulated
round-off: 7.7e-12 in Python and 4.7e-12 in Julia over 766781 steps of the cosine run. Any
``\lvert\Delta H\rvert`` reported at that same order is therefore arithmetic, not dynamics —
which is exactly what `cos` RK4 flow 1's 1.1e-11 against 6.2e-12 is. The Python's own closing
note makes the same point: a large growth ratio on an error still at round-off is "only the
random walk of the arithmetic".

## Three independent sources agree on the spectral radius

The explicit step size of every run in the manuscripts comes from
``\Delta t \le 2\sqrt{2}/\rho``, so ``\rho`` is the single number the whole RK4 comparison
rests on. It is now pinned from three directions:

| ``N`` | notes §5 | `kdvsim.py` | `stability_limit` |
|--:|--:|--:|--:|
| 12 | 115 | 115.4 | 115.4 |
| 16 | 274 | 274.2 | 274.2 |
| 24 | 907 | 907.2 | 907.2 |
| 32 | 2170 | 2170.0 | 2170.0 |

at ``p = 3`` on the first flow, from an ``O(1)`` single-mode initial condition. The value is
sensitive to the amplitude — ``u_0 = 0.1\sin x`` gives 112.5, 270.0, 900.2, 2160.1, off by two
per cent — which is why the notes now state their initial datum. Both the Julia test suite and
`kdvsim.py` assert these four figures, so a change in either assembly fails a test rather than
drifting silently away from the manuscript.

## Divergences found, and fixed

  - **`integrate` dropped the tail of every run.** `nout = nsteps ÷ stride` ran
    `nout * stride` steps while the tables went on reporting `nsteps`: 766400 of a claimed
    766781 for the cosine run, ending at ``T = 99.95`` rather than 100. Now `cld` with a final
    short window, as `plot_kdv_energies.py` does. The reported numbers do not move at the
    printed precision, but they were not measuring what they said.
  - **The explicit step size was clamped.** `scripts/kdv.jl` used
    `min(dt_rk4, case.dt)` where the Python uses `dt_rk4` outright. Invisible on four of the
    six cases, because there the stability limit is the smaller; on `soliton5` it forced
    0.0625 instead of 0.065 and hid the point the benchmark makes, which is that a wide,
    small-amplitude field costs an explicit method no step-count penalty at all.
  - **Window rounding.** `NT ÷ NOUT` against the Python's `round(NT / NOUT)`.

## Errata found in the Python

  - **`verify_burgers_discretisation.py` integrated four times too slowly.** Section 4
    transformed to ``\bar{u} = 2\sqrt{u}`` but used the bracket ``\mathbb{K}/4`` that belongs
    to ``\bar{u} = \sqrt{u}``. The transformation and the factor move together:
    ``\sqrt{u}`` with ``\mathbb{K}/4``, or ``2\sqrt{u}`` with ``\mathbb{K}`` and no factor.

    Nothing in the script could see it, and this is the interesting part: a constant rescaling
    of a Poisson vector field preserves every Casimir, every energy bound and every
    convergence *rate*, so all three of its checks passed. Only a comparison against the
    pushforward of the ``u``-space field catches it. That comparison is now a test here — see
    the `LagrangeSpace` Burgers testset — and a check there.

    The package did not inherit the error, its dynamics being carried in ``u``; but
    `to_sqrt_variables` was in the ``2\sqrt{u}`` convention while the docstring beside it
    displayed ``\mathbb{K}/4``. Both are now the manuscript's ``\bar{u} = \sqrt{u}``.
  - A stale comment in `kdvsim.py` put the explicit cosine run at "two million steps"; it is
    766781.

## The sign convention was flipped, and both sides agree afterwards

The package and both manuscripts now use the textbook ``u_t + 6uu_x + u_{xxx} = 0``, so the
``\operatorname{sech}^2`` solitons are elevations. The older ``u_t = 6uu_x - u_{xxx}`` is
reached by ``u \to -u``.

The one step that cannot be guessed is the second structure: it does **not** simply change
sign. ``4u\partial_x + 2u_x - \partial_x^3`` becomes ``-4u\partial_x - 2u_x -
\partial_x^3`` — the transport part flips, the third derivative does not. That was derived
and checked in `verify_kdv_continuous.py`, which works in jet variables, *before* any prose
or code was touched; in the package it is the single character of `AffineBracket`'s `scale`,
from `2` to `-2`, `K0` being the ``-\partial_x^3`` block and keeping its own sign.

Consequences, all verified:

  - ``H_1 = \tfrac12\int(u_x^2 - 2u^3)``, ``\delta H_1/\delta u = -3u^2 - u_{xx}``;
    ``H_2`` unchanged; ``H_3 = \int(\tfrac52 u^4 - 5uu_x^2 + \tfrac12 u_{xx}^2)``.
  - The Magri rung ``\mathbb{P}^2 g = 2\mathbb{P}^1 \partial H_2/\partial\hat{u}`` becomes
    ``-2``, for the same reason.
  - The Miura map becomes ``u = -(v^2 + v_x)``, its Jacobian ``-(2v + \partial_x)``. The
    overall sign cancels in ``\mathbb{L}\mathbb{P}^1\mathbb{L}^T``, so the bracket is
    unchanged; the chart is not, and the Hill potential becomes ``-u_h``.
  - The image of the Miura map moves to ``C_{0,d} \le 0`` — and the solitons move with it, so
    they are excluded exactly as before. **The sign convention was never what restricted the
    chart.** The obstruction is spectral: a soliton is a reflectionless potential *with* a
    bound state, and the chart covers the fields whose Hill operator has none.
  - The Hill threshold for ``u = c + \sin x + 0.4\cos 2x`` at ``N = 20``, ``p = 3`` moves
    from ``+0.470894`` to ``-0.328741``, the admissible side becoming the negative one. The
    magnitude differs because ``\pm(\sin x + 0.4\cos 2x)`` are not translates of one another.
    Measured by bisection it is ``-0.32893``, four digits rather than six: at the threshold
    ``\psi`` acquires a zero and ``v = \psi_x/\psi`` is unbounded, so the Newton iteration
    and not the criterion sets the resolution.

**Julia and Python agree after the flip**, on `soliton1` at ``N = 128``:
``H_1 = -6.39998180``, ``H_2 = +2.66666660``, ``C_0 = +4.00000000`` in both, and
``\rho = 115.4, 274.2, 907.2, 2170.0`` in both. ``H_1`` and ``H_2`` are the values the *old*
convention reported on the *old* data, since each Hamiltonian carries the matching sign; only
``C_0`` flipped.

## The Miura chart now covers every example

The restriction was never the sign convention, so flipping it did not help — and the fix is
not a shift bolted on but the spectral parameter the Riccati substitution always had:
``u = -(v^2 + v_x) - \lambda`` turns into ``(-\partial_x^2 - u)\psi = \lambda\psi``, so
``\lambda`` is an eigenvalue parameter of the Hill operator that decides invertibility and a
preimage exists exactly when ``\lambda < \lambda_0``.

Measured, at ``p = 3``, with `miura_lambda` a tenth below each threshold:

| example | ``C_{0,d}`` | ``\lambda_0`` | at ``\lambda = 0`` | at ``\lambda`` | round trip |
|:--|--:|--:|:--|:--|--:|
| ``\cos x``, ``N = 20`` | 0.000 | -0.3785 | none | reached | 1.2e-15 |
| ``\sin x + 0.4\cos 2x``, ``N = 20`` | 0.000 | -0.3287 | none | reached | 1.3e-15 |
| single soliton, ``N = 64`` | +4.000 | -0.9999 | none | reached | 1.1e-15 |
| two-soliton, ``N = 64`` | +8.000 | -1.4394 | none | reached | 1.6e-15 |

and the pushforward stays *exactly* antisymmetric with its Jacobiator at ``10^{-16}`` for
``\lambda = -0.5, -2, -5``. Both the Julia test suite and `verify_kdv_miura.py` §6b assert
this.

Two caveats worth stating rather than discovering:

  - ``\mathbb{L}\mathbb{P}^1\mathbb{L}^T`` is **not** the Galerkin ``\mathbb{P}^2 -
    4\lambda\mathbb{P}^1``; measured, they differ by 30-57 %. That difference *is* the Jacobi
    defect the Miura construction exists to remove, so the identity must not be asserted.
  - In the ``\hat{v}`` chart ``C_{0,d}`` is the mKdV momentum, not a Casimir, so a Miura run
    loses it — 2.1e-6 on the single soliton. The two methods exchange roles exactly as at
    ``\lambda = 0``: ``\tilde{H}`` is quartic in ``\hat{v}`` and ``C_{0,d}`` quadratic, so the
    discrete gradient holds the energy (1.3e-13) and the symplectic midpoint rule the mass.

### What the flip cost: four duplicated formulas

Each of these reimplements something for speed and had to be changed in step. Every one
produced a plausible-looking wrong answer rather than an error, and each was caught by a
test rather than by reading:

  - `invariants(::KdVSystem, …)` and `invariants(::MiuraSystem, …)` in `src/diagnostics.jl`,
    which assemble ``H_1`` from one field evaluation instead of calling
    `hamiltonian(KdVHamiltonian1(), …)`. Left stale, they made implicit midpoint report a
    373 % drift in a quantity the flow conserves exactly.
  - `kdvsim.Sys.invariants` and `kdvsim.Sys.H1`, the same duplication on the Python side.
  - `kdvsim.MiuraSys.Mmap`, which reimplements `splinetools.miura_map`.
  - Local ``H_1`` lambdas in `verify_kdv_semidiscrete.py` and `verify_kdv_timedisc.py`.

The lesson is the same as the Burgers factor above: a sign that moves in one copy and not
the other is invisible to every structural test and shows up only as a drift.

## The notes' own error list, dispositioned

`discrete-kdv-brackets-notes.tex` carries a standing list of known defects. Every item is
accounted for here.

| item | in the package? | disposition |
|:--|:--|:--|
| Errors 1 — mass-Casimir index and factor | no | paper-only. `MassCasimir` sidesteps it, arguing from ``g^T \mathbb{P}^2 \partial H_2/\partial\hat{u} = 0`` pointwise rather than ``g^T\mathbb{P}^2 = 0``, which is the correct weaker claim. |
| Errors 2 — dangling RK4 cross-reference | n/a | **closed.** The measurement was missing, not the reference; §5 now carries it, with the ``\rho`` table above. |
| Errors 3 — strict ``C_{0,d} > 0`` | **yes**, six sites | **fixed** in both, to ``C_{0,d} \ge 0`` with the equality case named. |
| Errors 4 — zero-mode theorem skips ``N = 4`` | no | paper-only; no counterpart here. |
| Numbers 1-3 | no | paper-only; each needs a measurement in the notes' own sections. |
| Numbers 4 — unstated initial data | partly | the time-discretisation datum is now stated. The disagreement between the *Miura* runs' data is still open. |
| Overstatements 1, 5 | no | paper-only. |
| Overstatements 2 — critical set of ``D\mathcal{M}_h`` | **yes** | **fixed.** Condition numbers are now given as evidence, and the exact statement that does hold is added: ``\int v_h`` is a Casimir of ``\mathbb{P}^1`` in this chart, so the flow cannot reach the degeneracy from data off it. Now a test. |
| Overstatements 3 — Ge-Marsden hypotheses | **yes**, four sites | **fixed** in both; the theorem is now the reason to expect the trade-off, not a proof of it here. |
| Overstatements 4 — modified Hamiltonian | no | paper-only. |
| Wording 1-5 | no | paper-only. |

## Confirmed correct

Checked and reproduced, and worth recording so that their absence above is not read as
neglect:

  - the Burgers bracket ``\mathbb{J} = \sqrt{u}\,\mathbb{K}\sqrt{u}`` reproducing
    ``u_t = 3uu_x`` at second order, and its exact Casimir ``2\int\sqrt{u}``;
  - the even-rank argument and the spurious sawtooth Casimir at even ``N``, with the
    odd-``N`` advice;
  - the whole Miura chain: ``C_{0,d} = \int v_h^2`` exactly, the sharp criterion
    ``\lambda_0(-\partial_x^2 + u_h) > 0``, two Floquet preimages with opposite ``\int v``,
    the threshold ``-\lambda_0 = 0.470894``, and that none of the standard examples is in the
    image;
  - both KdV flows reproducing the right-hand side of the KdV equation, and the Magri
    hierarchy step;
  - the Jacobiator of the Galerkin second bracket flat at ``\approx 0.42`` under refinement,
    against round-off for the Miura bracket.

## Two latent faults found by making the tabulation sparse

Neither was reachable before, which is why neither had been noticed.

  - **`project!` never worked on a `LagrangeSpace`.** `mass_factorization` returned a bare
    dense `Cholesky`, and `mass_solve!` — which the generic `project!` calls — has methods only
    for a `MassOperator`. Every in-place projection onto a Lagrange space was a `MethodError`.
    Nothing in the suite exercised it. It works now, and is tested.
  - **`AffineBracket` densified the tabulation it was handed.** Its `Ψ` field was declared
    `Matrix{T}` and its constructor called `Matrix(Ψ)`, so the *spline* space's already-sparse
    table was thrown away for the second KdV bracket and the first Camassa-Holm one — both of
    which contract it on every Newton iteration. Measured at ``N = 384``, one `jacobian` of
    the second flow: 7.24 ms before, 1.26 ms after, identical results. `bracket_directional`
    alone is 25 times faster; the `LagrangeSpace` assemblies about 30 times.

    The two ``O(N^3)`` tensor routines densify locally and on purpose — they random-access the
    table rather than contracting it — and are off the time loop by construction.

## What cannot be verified

**Camassa-Holm has no reference.** Grep for `camassa`, `peakon`, `holm` or `helmholtz` over
both `Scripts/` directories returns nothing, and neither manuscript covers it. The
implementation here is a *prototype* in the strict sense: it is internally consistent and its
structural tests pass, but there is no independent statement of what it should reproduce. It
should not be read as verified.
