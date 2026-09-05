#!/usr/bin/env sh
#
# Run the fable verification scripts.  Exits nonzero if any check fails.
#
#     sh scripts/run_fable.sh [script ...]
#
# These sit apart from `run_all.jl` for one reason: `fable/III4_multifield_derivative.jl` takes
# about eleven minutes, longer than the whole of `run_all.jl` put together, so folding them in
# would turn a two-minute driver into a fifteen-minute one.  They are verification scripts on the
# same terms as the rest -- not exploratory, unlike `search_dirac_variants.jl` -- and are run
# before anything in `Papers/Poisson Brackets from Four Brackets/` is changed.
#
# The order below is the order the results were established: I.1, II.1, then the III series, each
# refutation immediately after the script it attacks.  `fable/fabletools.jl` is an include shared
# by III1 and III4, not a driver, and so is not listed.
#
# Each script runs in its own process, so a failure is contained and reported rather than taking
# the runner down.  Six are pure Julia; `fable/II1_semidirect_flat.jl` also needs SymPy, which
# CondaPkg provisions through SymPyPythonCall on first use, so a cold first run is much longer
# than the thirteen minutes the seven take once that environment exists.

set -u

# CDPATH= so a user's CDPATH cannot redirect the cd.  No `--`: POSIX gives `dirname` no options at
# all, and $0 is a path to this script, so there is nothing for an option parser to mistake.
root=$(CDPATH= cd "$(dirname "$0")/.." && pwd)

if [ "$#" -eq 0 ]; then
    set -- fable/I1_vlasov_uN.jl \
        fable/II1_semidirect_flat.jl \
        fable/III1_structure_function_4d.jl \
        fable/III1_refutation.jl \
        fable/III2_matrix_structure_function.jl \
        fable/III2_refutation.jl \
        fable/III4_multifield_derivative.jl
fi

rule="==================================================================="
failed=""

for s in "$@"; do
    path="$root/scripts/$s"
    if [ ! -f "$path" ]; then
        printf '%s\n%s\n' "$rule" "$s -- no such script"
        failed="$failed $s(missing)"
        continue
    fi
    printf '%s\n%s\n' "$rule" "$s"
    if julia --project="$root/scripts" "$path"; then
        :
    else
        failed="$failed $s"
    fi
done

printf '%s\n' "$rule"
if [ -n "$failed" ]; then
    printf 'FAILED:%s\n' "$failed"
    exit 1
fi
printf 'all fable scripts passed.\n'
