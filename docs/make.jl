using Documenter, PMP

CI = get(ENV, "CI", nothing) == "true"

makedocs(sitename = "PMP.jl",
    format = Documenter.HTML(
    prettyurls = CI,
    ),
    pages = [
       "GettingStarted.md",
       "MoistureMaximization.md",
       "OtherMethods.md",
       "index.md"]
)

if CI
    deploydocs(
    repo   = "github.com/JuliaExtremes/PMP.jl.git",
    devbranch = "dev",
    versions = ["stable" => "v^", "v#.#", "main"],
    push_preview = false,
    target = "build"
    )
end