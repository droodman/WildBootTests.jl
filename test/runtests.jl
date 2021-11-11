using WildBootTests
using StatFiles, StatsModels, DataFrames, DataFramesMeta, BenchmarkTools, Plots, CategoricalArrays, Random, StableRNGs

open("unittests.log", "w") do log  # use Github (Desktop) to detect changes in output

df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\collapsed.dta"))
dropmissing!(df)
f = @formula(hasinsurance ~ 1 + selfemployed + post + post_self)
f = apply_schema(f, schema(f, df, Dict(:hasinsurance => ContinuousTerm)))
resp, predexog = modelcols(f, df)

println(log, "\nboottest post_self=.04, weight(webb)")
test = wildboottest(([0 0 0 1.], [.04]); resp, predexog, clustid=df.year, auxwttype=WildBootTest.webb, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nboottest post_self=.04, weight(webb) reps(9999999) noci")
test = wildboottest(([0 0 0 1.], [.04]); resp, predexog, clustid=df.year, reps=9999999, auxwttype=WildBootTest.webb, getCI=false, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test))")


println(log, "\nregress hasinsurance selfemployed post post_self, cluster(year)")
println(log, "boottest (post_self=.05) (post=-.02), reps(9999) weight(webb)")
test = wildboottest(([0 0 0 1.; 0 0 1. 0], [0.05; -0.02]); resp, predexog, clustid=df.year, reps=9999, auxwttype=WildBootTest.webb, rng=StableRNG(1231))
println(log, "F=$(teststat(test)) p=$(p(test))")

df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))[:,[:wage; :tenure; :ttl_exp; :collgrad; :industry]]
dropmissing!(df)
desc = describe(df, :eltype)
for i in axes(desc, 1)  # needed only for lm()
  if desc[i,:eltype] == Float32
    sym = desc[i,:variable]
    @transform!(df, @byrow $sym=Float64($sym))
  end
end
f = @formula(wage ~ 1 + tenure + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)

println(log, "\nconstraint 1 ttl_exp = .2")
println(log, "cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(industry)")
println(log, "boottest tenure")
test = wildboottest(([0 1. 0 0], [.0]); H₁=([0 0 1. 0], [.2]), resp, predexog, clustid=df.industry, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)")
println(log, "boottest tenure, ptype(equaltail)")
df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
predendog, inst = modelcols(ivf, df)
test = wildboottest(([0 0 0 1.], [.0]); resp, predexog, predendog, inst, clustid=df.industry, small=false, reps=9999, ptype=WildBootTest.equaltail, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nboottest tenure, ptype(equaltail) reps(99999) weight(webb) stat(c)")
test = wildboottest(([0 0 0 1.], [.0]); resp, predexog, predendog, inst, clustid=df.industry, small=false, reps=9999, auxwttype=WildBootTest.webb, bootstrapc=true, ptype=WildBootTest.equaltail, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
plot(plotpoints(test)...)

println(log, "\nboottest, ar")
test = wildboottest(([0 0 0 1.], [.0]); resp, predexog, predendog, inst, clustid=df.industry, small=false, ARubin=true, reps=9999, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
plot(plotpoints(test)...)

println(log, "\nscoretest tenure")
test = wildboottest(([0 0 0 1.], [.0]); resp, predexog, predendog, inst, clustid=df.industry, small=false, reps=0, scorebs=true, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nwaldtest tenure")
test = wildboottest(([0 0 0 1.], [.0]); resp, predexog, predendog, inst, clustid=df.industry, small=false, reps=0, imposenull=false, scorebs=true, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nivregress liml wage (tenure = collgrad ttl_exp), cluster(industry)")
println(log, "boottest tenure")
df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :tenure, :ttl_exp, :collgrad, :industry]]
dropmissing!(df)
f = @formula(wage ~ 1)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
ivf = @formula(tenure ~ collgrad + ttl_exp)
ivf = apply_schema(ivf, schema(ivf, df))
predendog, inst = modelcols(ivf, df)
test = wildboottest(([0 1.], [.0]); resp, predexog, predendog, inst, LIML=true, clustid=df.industry, small=false, reps=999, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nivreg2 wage collgrad smsa race age (tenure = union married), cluster(industry) fuller(1)")
println(log, "boottest tenure, nograph weight(webb) reps(9999)")
df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry]]
dropmissing!(df)
f = @formula(wage ~ 1 + collgrad + smsa + race + age)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
ivf = @formula(tenure ~ union + married)
ivf = apply_schema(ivf, schema(ivf, df))
predendog, inst = modelcols(ivf, df)
test = wildboottest(([0 0 0 0 0 1.], [0.]); resp, predexog, predendog, inst, Fuller=1, clustid=df.industry, small=false, reps=9999, auxwttype=WildBootTest.webb, rng=StableRNG(1231))
println(log, "z=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nareg wage ttl_exp collgrad tenure [aw=hours] if occupation<., cluster(age) absorb(industry)")
println(log, "boottest tenure, cluster(age occupation) bootcluster(occupation)")
df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :ttl_exp, :collgrad, :tenure, :age, :industry, :occupation, :hours]]
dropmissing!(df)
f = @formula(wage ~ ttl_exp + collgrad + tenure)  # constant unneeded in FE model
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)
test = wildboottest(([0 0 1.], [0.]); resp, predexog, clustid=Matrix(df[:, [:occupation, :age]]), nbootclustvar=1, nerrclustvar=2, obswt=df.hours, feid=df.industry, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

println(log, "\nglobal pix lnkm pixpetro pixdia pixwaterd pixcapdist pixmal pixsead pixsuit pixelev pixbdist")
println(log, "global geo lnwaterkm lnkm2split mean_elev mean_suit malariasuit petroleum diamondd")
println(log, "global poly capdistance1 seadist1 borderdist1")
println(log, "encode pixwbcode, gen(ccode)  // make numerical country identifier")
println(log, "qui areg lnl0708s centr_tribe lnpd0 \$pix \$geo \$poly, absorb(ccode)")
println(log, "boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(ccode)")
println(log, "boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(pixcluster)")
println(log, "boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(ccode pixcluster)")
df = DataFrame(load(raw"d:\OneDrive\Documents\Work\Econometrics\Wild cluster\pixel-level-baseline-final.dta"))
pix  = [:lnkm, :pixpetro, :pixdia, :pixwaterd, :pixcapdist, :pixmal, :pixsead, :pixsuit, :pixelev, :pixbdist]
geo  = [:lnwaterkm, :lnkm2split, :mean_elev, :mean_suit, :malariasuit, :petroleum, :diamondd]
poly = [:capdistance1, :seadist1, :borderdist1]
df = df[:,[pix; geo; poly; :lnl0708s; :centr_tribe; :lnpd0; :pixwbcode; :pixcluster]]
dropmissing!(df)
df.ccode = levelcode.(categorical(df.pixwbcode, compress=true))
df.pixcode = levelcode.(categorical(df.pixcluster, compress=true))
f = Term(:lnl0708s) ~ sum(term.([:centr_tribe; :lnpd0; pix; geo; poly]))
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)

test = wildboottest(([1 zeros(1,size(predexog,2)-1)], [0.]); resp, predexog, clustid=Matrix(df[:, [:ccode, :pixcode]]), nbootclustvar=1, nerrclustvar=2, feid=df.ccode, reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([1 zeros(1,size(predexog,2)-1)], [0.]); resp, predexog, clustid=Matrix(df[:, [:pixcode, :ccode]]), nbootclustvar=1, nerrclustvar=2, feid=df.ccode, reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([1 zeros(1,size(predexog,2)-1)], [0.]); resp, predexog, clustid=Matrix(df[:, [:pixcode, :ccode]]), nbootclustvar=2, nerrclustvar=2, feid=df.ccode, reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")


println(log, "\ninfile coll merit male black asian year state chst using regm.raw, clear")
println(log, "qui regress coll merit male black asian i.year i.state if !inlist(state,34,57,59,61,64,71,72,85,88), cluster(state)	")
println(log, "generate individual = _n  // unique ID for each observation")
println(log, "boottest merit, nogr reps(9999) gridpoints(10)  // defaults to bootcluster(state)")
println(log, "boottest merit, nogr reps(9999) gridpoints(10) nonull")
println(log, "boottest merit, nogr reps(9999) gridpoints(10) bootcluster(state year)")
println(log, "boottest merit, nogr reps(9999) gridpoints(10) nonull bootcluster(state year)")
println(log, "boottest merit, nogr reps(9999) gridpoints(10) bootcluster(individual)")
println(log, "boottest merit, nogr reps(9999) gridpoints(10) nonull bootcluster(individual)")
df = DataFrame(load(raw"d:\OneDrive\Documents\Work\Econometrics\Wild cluster\regm.dta"))
df = DataFrame(coll=Bool.(df.coll), merit=Bool.(df.merit), male=Bool.(df.male), black=Bool.(df.black), asian=Bool.(df.asian), state=categorical(Int8.(df.state)), year=categorical(Int16.(df.year)))
dropmissing!(df)
df = df[df.state .∉ Ref([34,57,59,61,64,71,72,85,88]),:]
f = @formula(coll ~ 1 + merit + male + black + asian + year + state)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)

test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid=levelcode.(df.state), gridpoints=[10], reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid=levelcode.(df.state), reps=9999, imposenull=false, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid=[levelcode.(df.year) levelcode.(df.state)], nbootclustvar=2, nerrclustvar=1, reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid=[levelcode.(df.year) levelcode.(df.state)], nbootclustvar=2, nerrclustvar=1, reps=9999, imposenull=false, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid= [collect(1:nrow(df)) levelcode.(df.state)], nbootclustvar=1, nerrclustvar=1, reps=9999, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")
test = wildboottest(([0. 1 zeros(1,size(predexog,2)-2)], [0.]); resp, predexog, clustid= [collect(1:nrow(df)) levelcode.(df.state)], nbootclustvar=1, nerrclustvar=1, reps=9999, imposenull=false, rng=StableRNG(1231))
println(log, "t=$(teststat(test)) p=$(p(test)) CI=$(CI(test))")

# df = DataFrame(load(raw"D:\OneDrive\Documents\Macros\WildBootTest.jl"))
# df = df[:, [:proposition_vote, :treatment, :ideology1, :log_income, :Q1_immigration, :group_id1, :group_id2]]
# dropmissing!(df)
# f = @formula(proposition_vote ~ treatment + ideology1 + log_income + Q1_immigration + 1)
# f = apply_schema(f, schema(f, df, Dict(:Q1_immigration  => CategoricalTerm)))
# resp, predexog = modelcols(f, df)
# test = WildBootTest.wildboottest(([0 0 0 0 0 0 0 0 0 0 0 0 1.], [0.]); resp, predexog, clustid=Matrix(df[:, [:Q1_immigration, :group_id1, :group_id2]]), rng=StableRNG(1231))

end