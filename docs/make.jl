push!(LOAD_PATH,"../src/")
using Documenter, WildBootTests

makedocs(sitename="WildBootTests.jl",
         authors="David Roodman",
         pages = [
                  "Overview" => "index.md",
                  "OLS examples" => "OLS examples.md",
                  "IV/2SLS examples" => "IV-2SLS examples.md",
                  "Public functions and types" => "exported.md"
                 ]
        )

deploydocs(
  deploy_config = Documenter.GitHubActions(),
  repo = "github.com/droodman/WildBootTests.jl.git",
)