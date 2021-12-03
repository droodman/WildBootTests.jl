# OLS examples
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

julia> R = [0 1]; r = [1];  # put null that coefficient on x = 1 in Rβ̂ = r form, where β̂ is parameter vector

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
## Further examples
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