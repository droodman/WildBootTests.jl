module WildBootTests
export BootTestResult, wildboottest, wildboottest!, waldtest, scoretest, teststat, stattype, p, padj, reps, repsfeas, nbootclust, dof, dof_r, plotpoints, peak, ci, dist, statnumer, statvar, auxweights

using LinearAlgebra, Random, Distributions, SortingAlgorithms, Printf, LoopVectorization, ThreadsX, SharedArrays, SparseArrays

include("structs.jl")
include("utilities.jl")
include("estimators.jl")
include("init.jl")
include("WRE.jl")
include("nonWRE.jl")
include("plot-CI.jl")
include("interface.jl")

# top-level computation routine for OLS/arubin (and score BS on IV/2SLS); split off to reduce latency when just doing WRE
function boottestOLSARubin!(o::StrBootTest{T}) where T
	if !o.initialized
		Init!(o)
	elseif !o.null
		NoNullUpdate!(o)
		return
	end

	if !isone(o.Nw)  # if more than one weight group to save memory, make on every call to boottest(), not just once in Init!()
		Random.seed!(o.rng,o.seed)
		MakeWildWeights!(o, last(o.WeightGrp[1])-1, first=true)
	end

	MakeInterpolables!(o)  # make stuff that depends linearly on r, possibly by interpolating, for first weight group

	MakeNonWREStats!(o, 1)  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
	for w ∈ 2:o.Nw
		MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)
		MakeNonWREStats!(o, w)
	end
	!o.bootstrapt && UpdateBootstrapcDenom!(o)

	o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
	o.distCDR = zeros(T,0,0)
	nothing
end

# top-level computation routine for non-arubin WRE; split off to reduce latency when just doing other tests
function boottestWRE!(o::StrBootTest{T}) where T
	if !o.initialized
		Init!(o)
	elseif !o.null
		NoNullUpdate!(o)
		return
	end

	if !isone(o.Nw)  # if more than one weight group to save memory, make on every call to boottest(), not just once in Init!()
		Random.seed!(o.rng,o.seed)
		MakeWildWeights!(o, last(o.WeightGrp[1])-1, first=true)
	end

	PrepWRE!(o)

	MakeWREStats!(o, 1)
	for w ∈ 2:o.Nw  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)
		MakeWREStats!(o, w)
	end
	!o.bootstrapt && UpdateBootstrapcDenom!(o)

	o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
	o.distCDR = zeros(T,0,0)
	nothing
end

# if not imposing null and we have returned to boottest!(), then dof=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
function NoNullUpdate!(o::StrBootTest{T} where T)
	if o.WREnonARubin
		o.numer[:,1] = o.R * o.DGP.Rpar * (isone(o.Repl.kZ) ? o.β̈sAs[1:1,1] : o.β̈s[:,1]) - o.r
	elseif o.arubin
		EstimateARubin!(o.DGP, o, false, o.r)
		o.numer[:,1] = @view o.DGP.β̈[o.kX₁+1:end,1]  # coefficients on excluded instruments in arubin OLS
	else
		o.numer[:,1] = o.R * (o.ml ? o.β̈  : iszero(o.κ) ? view(o.M.β̈  ,:,1) : o.M.Rpar * view(o.M.β̈  ,:,1)) - o.r  # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this
	end
	o.dist[1] = isone(o.dof) ? o.numer[1] / sqrtNaN(o.statDenom[1]) : o.numer[:,1]'invsym(o.statDenom)*o.numer[:,1]
	nothing
end

# compute bootstrap-c denominator from all bootstrap numerators
function UpdateBootstrapcDenom!(o::StrBootTest{T} where T)
	numer1 = o.numer[:,1]
	numersum = rowsum(o.numer) - numer1
	o.statDenom = (o.numer * o.numer' - numer1 * numer1' - numersum * numersum' / o.B) / o.B
	if o.sqrt
		o.dist = o.numer ./ sqrtNaN.(o.statDenom)
	else
		o.dist = colquadform(Matrix(invsym(o.statDenom)), o.numer)  # to reduce latency by minimizing @tturbo instances, work with negative of colquadform in order to fuse code with colquadformminus!
	end
	nothing
end
 
include("precompiler.jl")

end