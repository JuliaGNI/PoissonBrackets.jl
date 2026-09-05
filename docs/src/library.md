```@meta
CurrentModule = PoissonBrackets
```

# Library

The docstrings that belong to a topic live with it and are not repeated here: the three
equation modules on their own pages — [Korteweg-de Vries](@ref), [Camassa-Holm](@ref) and
[Burgers](@ref) — the Lie algebras and the exact linear algebra under
[Discrete Lie-Poisson brackets](@ref), the reduction and the broken hierarchical space under
[Dirac reduction](@ref), and the measurement layer under [Diagnostics](@ref).

The index below is complete regardless.

```@index
```

## Spaces — `spaces.jl`

`Pages` is matched against the *end* of each source path, so this filter has to be qualified
too: a bare `"spaces.jl"` also matches `tensorspaces.jl`, whose docstrings belong to the
section below and would otherwise be emitted twice.

```@autodocs
Modules = [PoissonBrackets]
Pages = ["src/spaces.jl"]
```

## Tensor-product spaces — `tensorspaces.jl`

```@autodocs
Modules = [PoissonBrackets]
Pages = ["tensorspaces.jl"]
```

## Brackets — `brackets.jl`

`Pages` is matched against the *end* of each source path, so this filter has to be qualified:
a bare `"brackets.jl"` also matches `fourbrackets.jl`, whose docstrings belong to
[Poisson brackets from four-brackets](@ref) and would otherwise be emitted twice.

```@autodocs
Modules = [PoissonBrackets]
Pages = ["src/brackets.jl"]
```

## Hamiltonians — `hamiltonians.jl`

```@autodocs
Modules = [PoissonBrackets]
Pages = ["hamiltonians.jl"]
```

## Flows — `flows.jl`

```@autodocs
Modules = [PoissonBrackets]
Pages = ["flows.jl"]
```

## Integrators — `integrators.jl`

```@autodocs
Modules = [PoissonBrackets]
Pages = ["integrators.jl"]
```

## The mixed formulation — `mixed.jl`

```@autodocs
Modules = [PoissonBrackets]
Pages = ["mixed.jl"]
```
