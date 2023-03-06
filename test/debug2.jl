# pushfirst!(LOAD_PATH, ".")
using WildBootTests, StatFiles, DataFrames, CategoricalArrays, StatsModels
try cd("test") catch end
df = DataFrame(load("collapsed.dta"))
dropmissing!(df)
f = @formula(hasinsurance ~ 1 + selfemployed + post + post_self)
f = apply_schema(f, schema(f, df, Dict(:hasinsurance => ContinuousTerm)))
resp, predexog = modelcols(f, df)
test = wildboottest([0 0 0 1], [.04]; resp, predexog, clustid=Int32.(df.year), auxwttype=:webb)
