module WildBootTests
export BoottestResult, wildboottest, AuxWtType, PType, MAdjType, DistStatType, 
       teststat, stattype, p, padj, reps, repsfeas, nbootclust, dof, dof_r, plotpoints, peak, CI, dist, statnumer, statvar, auxweights

using LinearAlgebra, Random, Distributions, SortingAlgorithms, LoopVectorization

include("utilities.jl")
include("StrBoottest.jl")
include("estimators.jl")
include("init.jl")
include("WRE.jl")
include("nonWRE.jl")
include("plot-CI.jl")
include("interface.jl")

# top-level computation routine for OLS/ARubin (and score BS on IV/2SLS); split off to reduce latency when just doing WRE
function boottestOLSARubin!(o::StrBootTest{T}) where T
  if !o.initialized
		Init!(o)
	elseif !o.null
		NoNullUpdate!(o)
		return
  end

  if o.Nw > 1  # if more than one weight group to save memory, make on every call to boottest(), not just once in Init!()
		Random.seed!(o.rng,o.seed)
		MakeWildWeights!(o, last(o.WeightGrp[1])-1, first=true)
	end

  MakeInterpolables!(o)  # make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for w ∈ 1:o.Nw  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		w > 1 && MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)

		MakeNonWREStats!(o, w)

		!o.bootstrapt && UpdateBootstrapcDenom!(o, w)
  end

  o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
  o.distCDR = zeros(T,0,0)
  nothing
end

# top-level computation routine for non-ARubin WRE; split off to reduce latency when just doing other tests
function boottestWRE!(o::StrBootTest{T}) where T
  if !o.initialized
		Init!(o)
	elseif !o.null
		NoNullUpdate!(o)
		return
  end

  if o.Nw > 1  # if more than one weight group to save memory, make on every call to boottest(), not just once in Init!()
		Random.seed!(o.rng,o.seed)
		MakeWildWeights!(o, last(o.WeightGrp[1])-1, first=true)
	end

  PrepWRE!(o)

  for w ∈ 1:o.Nw  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		w > 1 && MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)

		MakeWREStats!(o, w)

		!o.bootstrapt && UpdateBootstrapcDenom!(o, w)
  end

  o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
  o.distCDR = zeros(T,0,0)
  nothing
end

# if not imposing null and we have returned to boottest!(), then dof=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
function NoNullUpdate!(o::StrBootTest{T} where T)
  if o.WREnonARubin
		o.numer[:,1] = o.R * o.DGP.Rpar * o.β̂s[1] - o.r
  elseif o.ARubin
		EstimateARubin!(o.DGP, o, o.r)
		o.numer[:,1] = o.v_sd * @view o.DGP.β̂[o.kX₁+1:end,:]  # coefficients on excluded instruments in ARubin OLS
  else
		o.numer[:,1] = o.v_sd * (o.R * (o.ML ? o.β̂ : o.M.β̂) - o.r)  # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this
  end
  o.dist[1] = isone(o.dof) ? o.numer[1] / sqrt(o.statDenom[1]) : o.numer[:,1]'invsym(o.statDenom)*o.numer[:,1]
	nothing
end

# compute bootstrap-c denominator from all bootstrap numerators
function UpdateBootstrapcDenom!(o::StrBootTest{T} where T, w::Integer)
  if isone(w)
		tmp = o.numer[:,1]
		o.statDenom = o.numer * o.numer' - tmp * tmp'
		o.numersum = rowsum(o.numer) - tmp
  else
		statDenom .+= o.numer * o.numer'
		numersum .+= rowsum(o.numer)
  end
  if w == o.Nw  # last weight group?
		o.statDenom .= (o.statDenom - o.numersum * o.numersum' / o.B) / o.B
		if o.sqrt
			o.dist .= o.numer ./ sqrtNaN.(o.statDenom)
		else
			colquadform!(o.dist, invsym(o.statDenom), o.numer)
		end
  end
	nothing
end

if Base.VERSION >= v"1.4.2"  # source: https://timholy.github.io/SnoopCompile.jl/stable/snoopi_deep_parcel/#SnoopCompile.write
	include("../src/precompile_WildBootTests.jl")
	_precompile_()
end
 
end

# using StatFiles, StatsModels, DataFrames, DataFramesMeta, CategoricalArrays
# df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))
# df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
# dropmissing!(df)
# f = @formula(wage ~ 1 + ttl_exp + collgrad)
# f = apply_schema(f, schema(f, df))
# resp, predexog = modelcols(f, df)
# ivf = @formula(tenure ~ union)
# ivf = apply_schema(ivf, schema(ivf, df))
# predendog, inst = modelcols(ivf, df)
# test = WildBootTests.wildboottest([0 0 0 1], [0]; resp, predexog, predendog, inst, clustid=df.industry, small=false, reps=0)
