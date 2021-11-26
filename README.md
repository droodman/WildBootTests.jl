# WildBootTests.jl

This package performs wild bootstrap-based hypothesis tests in Julia at extreme speed. It is intended mainly for linear models: ordinary least squares (OLS) and instrumental variables/two-stage least squares (IV/2SLS). For an introduction to the wild bootstrap and the algorithms deployed here, see [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf).

The package offers/supports:
* Confidence intervals formed by inverting the test and iteratively searching for bounds.
* Multiway clustering.
* Arbitrary and multiple linear hypotheses in the parameters.
* Maintained linear constraints on the model (restricted OLS, IV/2SLS/LIML).
* One-way fixed effects.
* The wild bootstrap for OLS ([Wu 1986](https://doi.org/10.1214/aos/1176350142)).
* The Wild Restricted Efficient bootstrap (WRE) for IV/2SLS/LIML ([Davidson and MacKinnon 2010](https://doi.org/10.1198/jbes.2009.07221)).
* The subcluster bootstrap ([MacKinnon and Webb 2018]( https://doi.org/10.1111/ectj.12107)).
* Non-bootstrapped Wald, Rao, and Anderson-Rubin tests, optionally with multiway clustering.
* Generation of data for plotting of confidence curves or surfaces after one- or two-dimensional hypothesis tests.

WildBootTests.jl incorporates order-of-magnitude algorithmic speed-ups developed since [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf) for [OLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1586107-boottest-just-as-wild-10x-faster) and [IV/2SLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1597888-boottest-~100x-faster-after-iv-gmm). And it exploits the efficiency of Julia, such as by offering single-precision (Float32) computation.

The interface is low-level: the exported function `wildboottest()` accepts scalars, vectors, and matrices, not [DataFrame](https://github.com/JuliaData/DataFrames.jl)s or results from estimation functions such as [lm()](https://juliastats.org/GLM.jl/v1.5/). This design minimizes the package's dependency footprint while making the core functionality available to multiple programming environments, including Julia, R (through [JuliaCall](https://cran.r-project.org/web/packages/JuliaCall/index.html) and [wildboottestjlr](https://github.com/s3alfisc/wildboottestjlr)), and Python (through [PyJulia](https://github.com/JuliaPy/pyjulia)). A separate package will provide a higher-level Julia interface.

`wildboottest()` accepts many optional arguments. Most correspond to options of the Stata package `boottest`, which are documented in [Roodman et al. (2019), §7](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf#page=28). Julia-specific additions include an optional first argument `T`, which can be `Float32` or `Float64` to specify the precision of computation; and `rng`, which takes a random number generator such as `MersenneTwister(2302394)`.

# Example with output

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
# Examples omitting output
```
# use Webb instead of Rademacher weights, 99,999 replications instead of 999
test = wildboottest(R, r; resp, predexog, clustid, reps=99999, auxwttype=WildBootTests.webb)

# bootstrap in double-precision (Float64) instead of single (Float32)
# slow on first use because of recompile
test = wildboottest(Float64, R, r; resp, predexog, clustid)

# use guaranteed-stable random number generator for exact replicability
using StableRNGs
test = wildboottest(R, r; resp, predexog, clustid, rng=StableRNG(23948572))

# test that coefficient on intercept = 0 and coefficient on x = 1; plot confidence surface
test = wildboottest([1 0; 0 1], [0;1]; resp, predexog, clustid, reps=9999)
plot(plotpoints(test).X..., plotpoints(test).p, st=:contourf)

# multiway-cluster error variance by firm and year; bootstrap by firm
test = wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:firm, :year]]), nerrclustvar=2, nbootclustvar=1)

# same but bootstrap by year
test = wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:year, :firm]]), nerrclustvar=2, nbootclustvar=1)

# same but bootstrap by year-firm pair
test = wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:year, :firm]]), nerrclustvar=2, nbootclustvar=2)
```
