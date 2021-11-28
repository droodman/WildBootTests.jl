# low-level interface, meaning works with vectors and matrices, not data frames and estimation objects

"Structure to store wild bootstrap test results"
struct BoottestResult{T}
  stat::T; stattype::String
  p::T; padj::T
  reps::Int64; repsfeas::Int64
  NBootClust::Int64
  dof::Int64; dof_r::Int64
  plot::Union{Nothing, NamedTuple{(:X, :p), Tuple{Tuple{Vararg{Vector{T}, N} where N},Vector{T}}}}
  peak::Union{Nothing, NamedTuple{(:X, :p), Tuple{Vector{T}, T}}}
  CI::Union{Nothing, Matrix{T}}
  dist::Matrix{T}
  b::Vector{T}
  V::Matrix{T}
  auxweights::Union{Nothing,Matrix{T}}
  # M::StrBootTest
end

"Return test statistic subject to wild bootstrap test"
teststat(o::BoottestResult) = o.stat

"Return numerator of test statistic in wild bootstrap test"
statnumer(o::BoottestResult) = o.b

"Return denominator of test statistic in wild bootstrap test"
statvar(o::BoottestResult) = o.V

"""Return type of test statistic subject to wild bootstrap test: "t", "z", "F", or "χ²" """
stattype(o::BoottestResult) = o.stattype

"Return p value from wild bootstrap test"
p(o::BoottestResult) = o.p

"Returnp p value from wild bootstrap test after multiple-hypothesis adjustment, if any"
padj(o::BoottestResult) = o.padj

"Return requested number of replications in wild bootstrap test"
reps(o::BoottestResult) = o.reps

"Return actual number of replications in wild bootstrap test, subject to enumeration of Rademacher draws"
repsfeas(o::BoottestResult) = o.repsfeas

"Return number of bootstrapping clusters in wild bootstrap test"
NBootClust(o::BoottestResult) = o.NBootClust

"Return degrees of freedom wild bootstrap test"
dof(o::BoottestResult) = o.dof

"Return residual degrees of freedom wild bootstrap test"
dof_r(o::BoottestResult) = o.dof_r

"""
Return data for confidence plot of wild bootstrap test.
Return value is a 2-tuple with named entries `X` and `p` holding
the confidence sampling locations and p values respectively. `X` is in turn
a 1- or 2-tuple of vectors of sampling coordinates for each 
dimension of the tested hypothesis.
"""
plotpoints(o::BoottestResult) = o.plot

"Return parameter value with peak p value in wild bootstrap test"
peak(o::BoottestResult) = o.peak

"Return confidence interval matrix from wild bootstrap test, one row per disjoint piece"
CI(o::BoottestResult) = o.CI

"Return bootstrap distribution of statistic or statistic numerator in wild bootstrap test"
dist(o::BoottestResult) = o.dist

"Return auxilliary weight matrix for wild bootstrap"
auxweights(o::BoottestResult) = o.auxweights

using Printf
function Base.show(io::IO, o::BoottestResult{T}) where T
	print(io, "WildBootTests.BoottestResult{$T}\n\n")
	Printf.@printf(io, "p  = %5.3f\n", p(o))
	isdefined(o, :CI) && !isnothing(o.CI) && length(o.CI)>0 && print(io, "CI = $(CI(o))\n")
end

"""
wildboottest([T::DataType=Float32,] R::AbstractMatrix, r::AbstractVector; 
             resp, <optional keyword arguments>) -> WildBootTest.BoottestResult

Function to perform wild-bootstrap-based hypothesis test

# Positional arguments
* `T::DataType`: data type for inputs, results, and computations: Float32 (default) or Float64
* `R::AbstractMatrix` and `r::AbstractVector`: required matrix and vector expressesing the null Rβ=r; see Notes

# Required keyword argument
* `resp::AbstractVector`: response/dependent variable (y/y₁)

# Optional keyword arguments
* `predexog::AbstractVecOrMat`: exogenous predictors, including constant term, if any (X/X₁)
* `predendog::AbstractVecOrMat`: endogenous predictors (Y₂)
* `inst::AbstractVecOrMat`: instruments (X₂)
* `R1::AbstractMatrix` and `r1::AbstractVector`: model constraints; same format as for `R` and `r`
* `clustid::AbstractVecOrMat{<:Integer}`: data vector/matrix of error and bootstrapping cluster identifiers; see Notes 
* `nbootclustvar::Integer=1`: number of bootstrap-clustering variables
* `nerrclustvar::Integer=nbootclustvar`: number of error-clustering variables
* `hetrobust::Bool=true`: true unless errors are treated as iid
* `feid::AbstractVector{<:Integer}`: data vector for fixed effect group identifier
* `fedfadj::Bool=true`: true if small-sample adjustment should reflect number of fixed effects
* `obswt::AbstractVector`: observation weight vector; default is equal weighting
* `fweights::Bool=false`: true for frequency weights
* `maxmatsize::Number`: maximum size of auxilliary weight matrix (v), in gigabytes
* `ptype::PType=symmetric`: p value type (`symmetric`, `equaltail`, `lower`, `upper`)
* `bootstrapc::Bool=false`: true for bootstrap-c
* `LIML::Bool=false`: true for LIML or Fuller LIML
* `Fuller::Number`: Fuller factor
* `κ::Number`: fixed κ for _k_-class estimation
* `ARubin::Bool=false`: true for Anderson-Rubin test
* `small::Bool=true`: true for small-sample corrections
* `scorebs::Bool=false`: true for score bootstrap instead of wild bootstrap
* `reps::Integer=999`: number of bootstrap replications; `reps` = 0 requests classical Rao (or Wald) test if `imposenull` = `true` (or `false`)
* `imposenull::Bool=true`: true to impose null
* `auxwttype::AuxWtType=rademacher`: auxilliary weight type (`rademacher`, `mammen`, `webb`, `normal`, `gamma`)
* `rng::AbstractRNG=MersenneTwister()`: randon number generator
* `level::Number=.95`: significance level (0-1)
* `rtol::Number=1e-6`: tolerance for CI bound determination
* `madjtype::MAdjType=nomadj`: multiple hypothesis adjustment (`nomadj`, `bonferroni`, `sidak`)
* `NH₀::Integer=1`: number of hypotheses tested, including one being tested now
* `ML::Bool=false`: true for (nonlinear) ML estimation
* `scores::AbstractVecOrMat`: for ML, pre-computed scores
* `β::AbstractVector`: for ML, parameter estimates
* `A::AbstractMatrix`: for ML, covariance estimates
* `gridmin`: vector of graph lower bounds, max length 2, `missing` entries ask wildboottest() to choose
* `gridmax`: vector of graph upper bounds
* `gridpoints`: vector of number of sampling points
* `diststat::DistStatType=nodiststat`: t to save bootstrap distribution of Wald/χ²/F/t statistics; numer to save numerators thereof
* `getCI::Bool=true`: whether to return CI
* `getplot::Bool=getCI`: whether to generate plot data
* `getauxweights::Bool=false`: whether to save auxilliary weight matrix (v)

# Notes
The columns of `R` in the statement of the null should correspond to those of the matrix [`predexog` `predendog`],
where `predendog` is non-empty only in instrumental variables regression. 

Order the columns of `clustid` this way:
1. Variables only used to define bootstrapping clusters, as in the subcluster bootstrap.
2. Variables used to define both bootstrapping and error clusters.
3. Variables only used to define error clusters.
In the most common case, `clustid` is a single column of type 2.

The code does not handle missing data values: all data and identifier matrices must 
be restricted to the estimation sample.

"""
function wildboottest(T::DataType,
					  R::AbstractMatrix,
						r::AbstractVecOrMat;
					  resp::AbstractVector{<:Real},
					  predexog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  predendog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  inst::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  R1::AbstractMatrix=zeros(T,0,0),
						r1::AbstractVector=zeros(T,0),
					  clustid::AbstractVecOrMat{<:Integer}=zeros(Int,0,0),  # bootstrap-only clust vars, then boot&err clust vars, then err-only clust vars
					  nbootclustvar::Integer=1,
					  nerrclustvar::Integer=nbootclustvar,
						issorted::Bool=false,
					  hetrobust::Bool=true,
					  feid::AbstractVector{<:Integer}=Int8[],
					  fedfadj::Bool=true,
					  obswt::Union{AbstractVector{<:Real},UniformScaling{Bool}}=I,
					  fweights::Bool=false,
					  maxmatsize::Number=0,
					  ptype::PType=symmetric,
					  bootstrapc::Bool=false,
					  LIML::Bool=false,
					  Fuller::Number=0,
					  κ::Number=NaN,
					  ARubin::Bool=false,
					  small::Bool=true,
					  scorebs::Bool=false,
					  reps::Integer=999,
					  imposenull::Bool=true,
					  auxwttype::AuxWtType=rademacher,
					  rng::AbstractRNG=MersenneTwister(),
					  level::Number=.95,
					  rtol::Number=1e-6,
					  madjtype::MAdjType=nomadj,
					  NH₀::Integer=1,
					  ML::Bool=false,
					  scores::AbstractVecOrMat=Matrix{Float32}(undef,0,0),
					  β::AbstractVector=T[],
					  A::AbstractMatrix=zeros(T,0,0),
					  gridmin::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = T[],
					  gridmax::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = T[],
					  gridpoints::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Integer = Int64[],
					  diststat::DistStatType=nodist,
					  getCI::Bool=true,
					  getplot::Bool=getCI,
					  getauxweights::Bool=false)

  @assert length(predexog)==0 || nrows(predexog)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(predendog)==0 || nrows(predendog)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(inst)==0 || nrows(inst)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(feid)==0 || nrows(feid)==nrows(resp) "feid vector must have same height as data matrices"
  @assert length(clustid)==0 || nrows(clustid)==nrows(resp) "clustid must have same height as data matrices"
  @assert obswt==I || nrows(obswt)==nrows(resp) "obswt must have same height as data matrices"
  @assert nrows(R)==nrows(r) "Entries of H₀ tuple must have same height"
  @assert ncols(R)==ncols(predexog)+ncols(predendog) && isone(ncols(r)) "Wrong number of columns in null specification"
  @assert nrows(R1)==nrows(r1) "Entries of H₁ tuple must have same height"
  @assert length(R1)==0 || ncols(R1)==ncols(predexog)+ncols(predendog) "Wrong number of columns in model constraint specification"
  @assert nbootclustvar ≤ ncols(clustid) "nbootclustvar > width of clustid"
  @assert nerrclustvar ≤ ncols(clustid) "nerrclustvar > width of clustid"
  @assert reps ≥ 0 "reps < 0"
  @assert level ≥ 0. && level≤1. "level must be in the range [0,1]"
  @assert rtol > 0. "rtol ≤ 0"
  @assert NH₀ > 0 "NH₀ ≤ 0"
	if getplot || getCI
		@assert iszero(length(gridmin   )) || length(gridmin   )==nrows(R) "Length of gridmin doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridmax   )) || length(gridmax   )==nrows(R) "Length of gridmax doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridpoints)) || length(gridpoints)==nrows(R) "Length of gridpoints doesn't match number of hypotheses being jointly tested"
	end

  M = StrBootTest{T}(R, r, R1, r1, resp, predexog, predendog, inst, obswt, fweights, LIML, Fuller, κ, ARubin,
	                   reps, auxwttype, rng, maxmatsize, ptype, imposenull, scorebs, !bootstrapc, clustid, nbootclustvar, nerrclustvar, issorted, hetrobust, small,
	                   feid, fedfadj, level, rtol, madjtype, NH₀, ML, β, A, scores, getplot,
	                   map(x->ismissing(x) ? missing : T(x), gridmin),
	                   map(x->ismissing(x) ? missing : T(x), gridmax),
	                   map(x->ismissing(x) ? missing : Int32(x), gridpoints))

	if getplot || (level<1 && getCI)
		plot = getplot ? _getplot(M) : nothing
		peak = getpeak(M)
		CI = level<1 & getCI ? _getCI(M) : nothing
  else
		CI = plot = peak = nothing
  end

  BoottestResult{T}(getstat(M),
	                  isone(nrows(R)) ? (small ? "t" : "z") : (small ? "F" : "χ²"),
	                  getp(M), getpadj(M), getreps(M), getrepsfeas(M), getNBootClust(M), getdf(M), getdf_r(M), plot, peak, CI,
	                  getdist(M,diststat),
	                  getb(M), getV(M),
	                  getauxweights && reps>0 ? getauxweights(M) : nothing #=, M=#)
end

wildboottest(T::DataType, R, r::Number; kwargs...) = wildboottest(T, R, [r]; kwargs...)
wildboottest(R::AbstractMatrix, r::Union{Number,AbstractVecOrMat}; kwargs...) = wildboottest(Float32, R, r; kwargs...)

