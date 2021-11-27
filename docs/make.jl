push!(LOAD_PATH,"../src/")
using Documenter, WildBootTests

makedocs(sitename="WildBootTests.jl",
         authors="David Roodman",
  pages = [
        "Home" => "index.md",
        "Examples" => "README.md"
        ]
    )
