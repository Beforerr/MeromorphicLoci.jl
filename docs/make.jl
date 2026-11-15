using Documenter, MeromorphicLoci
using CairoMakie

readme = read(joinpath(@__DIR__, "..", "README.md"), String)

write(joinpath(@__DIR__, "src", "index.md"),
    replace(readme, "](docs/src/" => "](", "](./" => "](https://github.com/Beforerr/MeromorphicLoci.jl/blob/main/",
        r"```julia\n(?!findzeros\(f, zdomain)" => "```@example readme\n"))

makedocs(;
    sitename="MeromorphicLoci.jl",
    modules=[MeromorphicLoci],
    checkdocs=:exports,
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", size_threshold=nothing),
    pages=[
        "Home" => "index.md",
        "AAA" => "aaa.md",
        "Continuation" => "continuation.md",
        "ArgumentPrinciple" => "argument_principle.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo="github.com/Beforerr/MeromorphicLoci.jl.git", push_preview=true)
