```
using WildBootTests, CSV, DataFrames, StatsModels, GLM, Plots

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
