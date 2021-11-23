# WildBootTest.jl 0.1 2 November 2021

# MIT License
#
# Copyright (C) 2015-21: David Roodman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#   * The above copyright notice and this permission notice shall be included in all copies or substantial
#     portions of the Software.
#   * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
#     NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


module WildBootTests
export BoottestResult, wildboottest, AuxWtType, PType, MAdjType, DistStatType, teststat, stattype, p, padj, reps, repsfeas, NBootClust, dof, dof_r, plotpoints, peak, CI, dist, statnumer, statvar, auxweights

using LinearAlgebra, Random, Distributions, LoopVectorization, SortingAlgorithms

include("utilities.jl")
include("estimators.jl")
include("StrBoottest.jl")
include("init.jl")
include("WRE.jl")
include("nonWRE.jl")
include("plot-CI.jl")
include("interface.jl")

# main routine
function boottest!(o::StrBootTest{T}) where T
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

  o.WREnonARubin ? PrepWRE!(o) :
				           MakeInterpolables!(o)  # make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for w ∈ 1:o.Nw  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		w > 1 && MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)

		o.WREnonARubin ? MakeWREStats!(o, w) :
										 MakeNonWREStats!(o, w)

		!o.bootstrapt && UpdateBootstrapcDenom!(o, w)
  end

  o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
  o.distCDR = zeros(T,0,0)
  o.dirty = false
  o.initialized = true
end

# if not imposing null and we have returned to boottest!(), then dof=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
function NoNullUpdate!(o::StrBootTest{T} where T)
  if o.WREnonARubin
		o.numer[:,1] = o.R * o.DGP.Rpar * o.βs[1] - o.r
  elseif o.ARubin
		o.DGP.Estimate(o.r)
		o.numer[:,1] = o.v_sd * o.DGP.Rpar * @view o.DGP.β[o.kX₁+1:end,:] # coefficients on excluded instruments in ARubin OLS
  else
		o.numer[:,1] = o.v_sd * (o.R * (o.ML ? o.β : o.M.Rpar * o.M.β) - o.r) # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this
  end
  o.dist[1] = isone(o.dof) ? o.numer[1] / sqrt(o.statDenom[1]) : o.numer[:,1]'invsym(o.statDenom)*o.numer[:,1]
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
			o.dist[:,:] .= o.numer' ./ sqrtNaN.(o.statDenom)
		else
			o.dist = colquadform(invsym(o.statDenom), o.numer)
		end
  end
end

end

# using StatFiles, StatsModels, DataFrames, DataFramesMeta, BenchmarkTools, Plots, CategoricalArrays, Random, StableRNGs
# df = DataFrame(load(raw"d:\OneDrive\Documents\Work\Econometrics\Wild cluster\regm.dta"))
# df = DataFrame(coll=Bool.(df.coll), merit=Bool.(df.merit), male=Bool.(df.male), black=Bool.(df.black), asian=Bool.(df.asian), state=categorical(Int8.(df.state)), year=categorical(Int16.(df.year)))
# dropmissing!(df)
# df = df[df.state .∉ Ref([34,57,59,61,64,71,72,85,88]),:]
# f = @formula(coll ~ 1 + merit + male + black + asian + year + state)
# f = apply_schema(f, schema(f, df))
# resp, predexog = modelcols(f, df)
# test = WildBootTests.wildboottest(([0 1 zeros(1,size(predexog,2)-2)], [0]); resp, predexog, clustid= [collect(1:nrow(df)) levelcode.(df.state)], nbootclustvar=1, nerrclustvar=1, reps=9999, rng=StableRNG(1231))
