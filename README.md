# WildBootTests.jl
WildBootTests.jl performs wild bootstrap-based hypothesis tests at extreme speed. It is intended mainly for linear models: ordinary least squares (OLS) and instrumental variables/two-stage least squares (IV/2SLS). For an introduction to the wild bootstrap and the algorithms deployed here, see [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf).

## Documentation
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://droodman.github.io/WildBootTests.jl/dev)

## Example

```
using WildBootTests, CSV, DataFrames, StatsModels, Plots
d = download("https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv");
df = CSV.read(d, DataFrame);
f = @formula(y ~ 1 + x);                             # state OLS model
f = apply_schema(f, schema(f, df));                  # link model to data
resp, predexog = modelcols(f, df);                   # extract response & (exogenous) predictor variables
clustid = df.firm;                                   # extract clustering variable
R = [0 1]; r = [1];                                  # put null in Rβ = r form, where β is parameter vector

test = wildboottest(R, r; resp, predexog, clustid);  # run test
plot(plotpoints(test)...)                            # plot confidence curve
```
