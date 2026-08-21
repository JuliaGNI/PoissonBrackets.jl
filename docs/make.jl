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
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Discretisation" => "discretisation.md",
        "Equations" => ["Korteweg-de Vries" => "kdv.md",
                        "Camassa-Holm" => "camassaholm.md",
                        "Burgers" => "burgers.md"],
        "Integrators" => "integrators.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaGNI/PoissonBrackets.jl",
    devurl = "latest",
    devbranch = "main",
)
