# WildBootTests.jl

## Documentation
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://droodman.github.io/WildBootTests.jl/dev)

## Example

```
julia> using WildBootTests, CSV, DataFrames, GLM, Plots

julia> d = download("https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv");

julia> df = CSV.read(d, DataFrame);

julia> f = @formula(y ~ 1 + x);  # state OLS model

julia> f = apply_schema(f, schema(f, df));  # link model to data

julia> lm(f, df)  # run OLS for illustration; not needed for following lines
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Vector{Float64}}, GLM.DensePredChol{Float64, LinearAlgebra.CholeskyPivoted{Float64, Matrix{Float64}}}}, Matrix{Float64}}

y ~ 1 + x

Coefficients:
─────────────────────────────────────────────────────────────────────────
                 Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────────
(Intercept)  0.0296797   0.0283593   1.05    0.2954  -0.025917  0.0852764
x            1.03483     0.0285833  36.20    <1e-99   0.978798  1.09087
─────────────────────────────────────────────────────────────────────────

julia> resp, predexog = modelcols(f, df);  # extract response & (exogenous) predictor variables

julia> clustid = df.firm;  # extract clustering variable

julia> R = [0 1]; r = [1];  # put null that coefficient on x = 1 in Rβ = r form, where β is parameter vector

julia> test = wildboottest(R, r; resp=resp, predexog=predexog, clustid=clustid)
WildBootTests.BoottestResult{Float32}

p  = 0.492
CI = Float32[0.93461335 1.1347668]

julia> test = wildboottest(R, r; resp, predexog, clustid);  # same, using Julia syntactic sugar

julia> p(test)  # programmatically extract p value
0.49459493f0

julia> CI(test)  # programmatically extract confidence interval
1×2 Matrix{Float32}:
 0.934961  1.13469

julia> plot(plotpoints(test)...)  # plot confidence curve
```
