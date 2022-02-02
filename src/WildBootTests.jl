module WildBootTests
export BootTestResult, wildboottest, AuxWtType, PType, MAdjType, DistStatType, 
       teststat, stattype, p, padj, reps, repsfeas, nbootclust, dof, dof_r, plotpoints, peak, CI, dist, statnumer, statvar, auxweights

using LinearAlgebra, Random, Distributions, SortingAlgorithms, LoopVectorization

include("StrBootTest.jl")
include("utilities.jl")
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
			negcolquadform!(o.dist, -invsym(o.statDenom), o.numer)  # to reduce latency by minimizing @tturbo instances, work with negative of colquadform in order to fuse code with colquadformminus!
		end
  end
	nothing
end

include("precompile_WildBootTests.jl")  # source: https://timholy.github.io/SnoopCompile.jl/stable/snoopi_deep_parcel/#SnoopCompile.write
_precompile_()
 
end


using StableRNGs, Random, BenchmarkTools
N=1_000_00; G=40; k=2; l=4
Random.seed!(1231)
β=rand(); γ=rand(k); Π=rand(l)
W = rand(N,l)
u₂ = randn(N)
y₂ = W * Π + u₂
u₁ = u₂ + randn(N)
Z = rand(N,k)
y₁ = y₂ * β + Z * γ + u₁
ID = floor.(Int8, collect(0:N-1) / (N/G))
R = [zeros(1,k) 1]; r = [0]
# FEID = rand([1,2,3,4,5],N)
WildBootTests.wildboottest(Float32, R,#=.36=#.46; R1 = [1 zeros(1,k)], r1=[0], resp=y₁, predexog=Z, predendog=y₂, inst=W, clustid=ID, reps=999, issorted=true, rng=StableRNG(1231), LIML=true)
WildBootTests.wildboottest(Float64, R,.46; R1 = [1 zeros(1,k)], r1=[0], resp=y₁, predexog=Z, predendog=y₂, inst=W, clustid=ID, reps=999, issorted=true, rng=StableRNG(1231), LIML=true)

N=1_000_000; G=40; k=12; l=40
Random.seed!(1231)
β=rand(); γ=rand(k); Π=rand(l)
W = rand(N,l)
u₂ = randn(N)
y₂ = W * Π + u₂
u₁ = u₂ + randn(N)
Z = rand(N,k)
y₁ = y₂ * β + Z * γ + u₁
ID = floor.(Int8, collect(0:N-1) / (N/G))
R = [zeros(1,k) 1]; r = [0]
@btime WildBootTests.wildboottest(Float32, R,.4; resp=y₁, predexog=Z, predendog=y₂, inst=W, clustid=ID, reps=9999, issorted=true, getCI=false)

# @btime WildBootTests.wildboottest(Float32, R,.4; resp=y₁, predexog=Z, predendog=y₂, inst=W, clustid=ID, reps=9999, issorted=true, getCI=true);