push!(LOAD_PATH,"../src/")
using Documenter, WildBootTests

makedocs(sitename="WildBootTests.jl",
         authors="David Roodman",
         pages = [
                  "Overview" => "index.md",
                  "OLS examples" => "OLSexamples.md",
                  "IV/2SLS examples" => "IVexamples.md",
                  "Public functions and types" => "exported.md"
                 ]
        )

deploydocs(
  deploy_config = Documenter.GitHubActions(),
  repo = "github.com/droodman/WildBootTests.jl.git",
)

