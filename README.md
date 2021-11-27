# WildBootTests.jl

WildBootTests.jl performs wild bootstrap-based hypothesis tests at extreme speed. It is intended mainly for linear models: ordinary least squares (OLS) and instrumental variables/two-stage least squares (IV/2SLS). For an introduction to the wild bootstrap and the algorithms deployed here, see [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf).

The package offers and/orsupports:
* The wild bootstrap for OLS ([Wu 1986](https://doi.org/10.1214/aos/1176350142)).
* The Wild Restricted Efficient bootstrap (WRE) for IV/2SLS/LIML ([Davidson and MacKinnon 2010](https://doi.org/10.1198/jbes.2009.07221)).
* The subcluster bootstrap ([MacKinnon and Webb 2018](https://doi.org/10.1111/ectj.12107)).
* Non-bootstrapped Wald, Rao, and Anderson-Rubin tests, optionally with multiway clustering.
* Confidence intervals formed by inverting the test and iteratively searching for bounds.
* Multiway clustering.
* Arbitrary and multiple linear hypotheses in the parameters.
* Maintained linear constraints on the model (restricted OLS, IV/2SLS/LIML).
* One-way fixed effects.
* Generation of data for plotting of confidence curves or surfaces after one- or two-dimensional hypothesis tests.

WildBootTests.jl incorporates order-of-magnitude algorithmic speed-ups developed since [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf) for [OLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1586107-boottest-just-as-wild-10x-faster) and [IV/2SLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1597888-boottest-~100x-faster-after-iv-gmm). And it exploits the efficiency of Julia, for example by offering single-precision (Float32) computation.

The interface is low-level: the exported function `wildboottest()` accepts scalars, vectors, and matrices, not [DataFrame](https://github.com/JuliaData/DataFrames.jl)s or results from estimation functions such as [lm()](https://juliastats.org/GLM.jl/v1.5/). This design minimizes the package's dependency footprint while making the core functionality available to multiple programming environments, including Julia, R (through [JuliaConnectoR](https://cran.r-project.org/web/packages/JuliaConnectoR/index.html)), and Python (through [PyJulia](https://github.com/JuliaPy/pyjulia)). A separate package will provide a higher-level Julia interface.

`wildboottest()` accepts many optional arguments. Most correspond to options of the Stata package `boottest`, which are documented in [Roodman et al. (2019), §7](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf#page=28). Julia-specific additions include an optional first argument `T`, which can be `Float32` or `Float64` to specify the precision of computation; and `rng`, which takes a random number generator such as `MersenneTwister(2302394)`.

# OLS example with output

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
# OLS examples omitting output
```
# use Webb instead of Rademacher weights, 99,999 bootstrap replications instead of 999
wildboottest(R, r; resp, predexog, clustid, reps=99999, auxwttype=WildBootTests.webb)

# bootstrap in double-precision (Float64) instead of single (Float32)
# slow on first use because of recompile
wildboottest(Float64, R, r; resp, predexog, clustid)

# use guaranteed-stable random number generator for exact replicability
using StableRNGs
wildboottest(R, r; resp, predexog, clustid, rng=StableRNG(23948572))

# test that coefficient on intercept = 0 and coefficient on x = 1; plot confidence surface
test = wildboottest([1 0; 0 1], [0;1]; resp, predexog, clustid, reps=9999)
plot(plotpoints(test).X..., plotpoints(test).p, st=:contourf)

# multiway-cluster errors by firm and year; bootstrap by firm
wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:firm, :year]]), nerrclustvar=2, nbootclustvar=1)

# same but bootstrap by year
wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:year, :firm]]), nerrclustvar=2, nbootclustvar=1)

# same but bootstrap by year-firm pair
wildboottest(R, r; resp, predexog, clustid=Matrix(df[:,[:year, :firm]]), nerrclustvar=2, nbootclustvar=2)

# Rao/score test with multiway clustering of errors but no bootstrap
wildboottest(R, r; resp, predexog, predendog, inst, Matrix(df[:,[:year, :firm]]), reps=0)

# Same but Wald test: i.e., conventional, multiway clustered errors
wildboottest(R, r; resp, predexog, predendog, inst, clustid=Matrix(df[:,[:year, :firm]]), reps=0, imposenull=false)

# add year fixed effects to model; cluster by firm
wildboottest(R, r; resp, predexog, feid=df.year, clustid=df.firm)

# test hypotheses, while imposing model constraint that constant term = 0.2
R1 = [1 0]; r1 = [.2]
wildboottest(R, r; R1, r1, resp, predexog, clustid=df.firm)
```
# IV/2SLS examples omitting output
```
# specify exactly identified model: regress wage on on tenure, instrumented by union,
# controlling for ttl_exp and collgrad
d = download("http://www.stata-press.com/data/r8/nlsw88.dta", tempname() * ".dta")
df = DataFrame(load(d))[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
predendog, inst = modelcols(ivf, df)

# test that coefficient on tenure = 0, clustering errors by industry
R = [0 0 0 1]; r = [0]
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry)

# use equal-tailed instead of symmetric p value
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry, ptype=WildBootTests.equaltail)

# perform bootstrap-c instead of bootstrap-t, as advocated by Young (2019)
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry, bootstrapc=true)

# Rao/score test without bootstrap
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry, reps=0)

# Wald test without bootstrap
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry, reps=0, imposenull=false)

# Anderson-Rubin test that hypothesis holds and instrument is valid
wildboottest(R, r; resp, predexog, predendog, inst, clustid=df.industry, ARubin=true)

# modify model to drop controls and make ttl_exp an instrument
f = @formula(wage ~ 1)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
ivf = @formula(tenure ~ collgrad + ttl_exp)
ivf = apply_schema(ivf, schema(ivf, df))
predendog, inst = modelcols(ivf, df)

# test same hypothesis in context of LIML regression
R = [0 1]; r = [0]
wildboottest(R, r; resp, predexog, predendog, inst, LIML=true, clustid=df.industry)
```
