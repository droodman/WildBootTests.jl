# WildBootTests.jl
WildBootTests.jl performs wild bootstrap-based hypothesis tests at extreme speed. It is intended mainly for linear models: ordinary least squares (OLS) and instrumental variables/two-stage least squares (IV/2SLS). For an introduction to the wild bootstrap and the algorithms deployed here, see [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf). It is a Julia program, but can be accessed from other environments as demonstrated below.

## Documentation
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://droodman.github.io/WildBootTests.jl/dev/)

## Examples

In all cases, to install, start Julia and type `using Pkg; Pkg.add("WildBootTests")`.

### From Julia

```julia
using WildBootTests, CSV, DataFrames, StatsModels, Plots
d = download("https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv");
df = CSV.read(d, DataFrame);
f = @formula(y ~ 1 + x);                            # state OLS model
f = apply_schema(f, schema(f, df));                 # link model to data
resp, predexog = modelcols(f, df);                  # extract response & (exogenous) predictor variables
clustid = df.firm;                                  # extract clustering variable
R = [0 1]; r = [1];                                 # put null in Rβ = r form, where β is parameter vector

test = wildboottest(R, r; resp, predexog, clustid)  # run test & display results
plot(plotpoints(test)...)                           # plot confidence curve
```

### From R via fwildclusterboot
```R
library(fwildclusterboot)
df <- read.csv("https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv")
lm_fit <- lm(y ~ x, data = df)
boot_lm <- boottest(lm_fit, clustid = "firm", param = "x", r = 1, B = 999,  boot_algo = "WildBootTests.jl")
summary(boot_lm)
plot(boot_lm)
```

### From R via JuliaConnectoR
```julia
library(JuliaConnectoR)
startJuliaServer()
WildBootTests <- juliaImport("WildBootTests")
df <- read.csv(file = 'https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv')
R <- matrix(c(0,1), nrow=1); r <- c(1)
test <- WildBootTests$wildboottest(R, r, resp=df$y, predexog=cbind(1, df$x), clustid=df$firm)
test
WildBootTests$teststat(test)
WildBootTests$p(test)
WildBootTests$ci(test)
plotpoints <- WildBootTests$plotpoints(test)
plot(plotpoints$X[[1]], plotpoints$p, type="l")
```

### From Python via PyJulia
```python
from julia import WildBootTests as wbt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(r'https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv')
R = np.array([[0, 1]]); r = np.array([1])
resp = df.y.values
predexog = np.c_[np.ones(df.firm.size), df.x]
clustid = df.firm.values
test = wbt.wildboottest(R, r, resp=resp, predexog=predexog, clustid=clustid)
wbt.teststat(test)
wbt.p(test)
wbt.ci(test)
plotpoints = wbt.plotpoints(test)
plt.plot(plotpoints.X[0], plotpoints.p)
```

### From Stata via boottest
```
import delimited https://raw.github.com/vincentarelbundock/Rdatasets/master/csv/sandwich/PetersenCL.csv
regress y x, cluster(firm)
boottest x, reps(99999) julia
```
