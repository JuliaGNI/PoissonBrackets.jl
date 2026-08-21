using PoissonBrackets
using Documenter

DocMeta.setdocmeta!(PoissonBrackets, :DocTestSetup,
    :(using PoissonBrackets); recursive = true)

makedocs(;
    modules = [PoissonBrackets],
    authors = "Michael Kraus",
    repo = Remotes.GitHub("JuliaGNI", "PoissonBrackets.jl"),
    sitename = "PoissonBrackets.jl",
    checkdocs = :exports,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaGNI.github.io/PoissonBrackets.jl",
        # set explicitly: without a configured remote, Documenter cannot read the branch
        # from `git remote` and would default to "master"
        edit_link = "main",
        # the library page collects every docstring in the package and is legitimately
        # large; the default warning threshold is aimed at pages that are large by accident
        size_threshold_warn = 400 * 1024,
        size_threshold = 800 * 1024,
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Discretisation" => "discretisation.md",
        "Equations" => ["Korteweg-de Vries" => "kdv.md",
                        "Camassa-Holm" => "camassaholm.md",
                        "Burgers" => "burgers.md"],
        "Integrators" => "integrators.md",
        "Verification" => "verification.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaGNI/PoissonBrackets.jl",
    devurl = "latest",
    devbranch = "main",
)
