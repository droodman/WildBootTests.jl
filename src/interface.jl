# low-level interface, meaning works with vectors and matrices, not data frames and estimation objects

"Structure to store test results"
struct BootTestResult{T}
  stat::T; stattype::String
  p::T; padj::T
  reps::Int64; repsfeas::Int64
  nbootclust::Int64
  dof::Int64; dof_r::T
  plot::Union{Nothing, NamedTuple{(:X, :p), Tuple{Tuple{Vararg{Vector{T}, N} where N},Vector{T}}}}
  peak::Union{Nothing, NamedTuple{(:X, :p), Tuple{Vector{T}, T}}}
  CI::Union{Nothing, Matrix{T}}
  dist::Matrix{T}
  b::Vector{T}
  V::Matrix{T}
  auxweights::Union{Nothing,Matrix{T}}
  # M::StrBootTest
end

"Return test statistic"
teststat(o::BootTestResult) = o.stat

"Return numerator of test statistic"
statnumer(o::BootTestResult) = o.b

"Return denominator of test statistic"
statvar(o::BootTestResult) = o.V

"""Return type of test statistic subject: "t", "z", "F", or "χ²" """
stattype(o::BootTestResult) = o.stattype

"Return p value"
p(o::BootTestResult) = o.p

"Returnp p value after multiple-hypothesis adjustment, if any"
padj(o::BootTestResult) = o.padj

"Return requested number of replications"
reps(o::BootTestResult) = o.reps

"Return actual number of replications, subject to enumeration of Rademacher draws"
repsfeas(o::BootTestResult) = o.repsfeas

"Return number of bootstrapping clusters in test"
nbootclust(o::BootTestResult) = o.nbootclust

"Return degrees of freedom test"
dof(o::BootTestResult) = o.dof

"Return residual degrees of freedom test"
dof_r(o::BootTestResult) = o.dof_r

"""
Return data for confidence plot of test.
Return value is a 2-tuple with named entries `X` and `p` holding
the confidence sampling locations and p values respectively. `X` is in turn
a 1- or 2-tuple of vectors of sampling coordinates for each 
dimension of the tested hypothesis.
"""
plotpoints(o::BootTestResult) = o.plot

"Return parameter value with peak p value in test"
peak(o::BootTestResult) = o.peak

"Return confidence interval matrix from test, one row per disjoint piece"
CI(o::BootTestResult) = o.CI

"Return bootstrap distribution of statistic or statistic numerator in bootstrap test"
dist(o::BootTestResult) = o.dist

"Return auxilliary weight matrix for wild bootstrap"
auxweights(o::BootTestResult) = o.auxweights

using Printf
function Base.show(io::IO, o::BootTestResult{T}) where T
	print(io, "WildBootTests.BootTestResult{$T}\n\n")
	Printf.@printf(io, "%s = %5.3f\n", stattype(o)*repeat(' ',2-length(stattype(o))), teststat(o))
	Printf.@printf(io, "p  = %5.3f\n", p(o))
	isdefined(o, :CI) && !isnothing(o.CI) && length(o.CI)>0 && print(io, "CI = $(CI(o))\n")
end

# single entry point with arguments already converted to standardized types, to allow a small set of precompile() calls
function __wildboottest(
	R::Matrix{T},
	r::Vector{T};
	resp::Vector{T},
	predexog::Matrix{T},
	predendog::Matrix{T},
	inst::Matrix{T},
	R1::Matrix{T},
	r1::Vector{T},
	clustid::Matrix{Int64},
	nbootclustvar::Int8,
	nerrclustvar::Int8,
	issorted::Bool,
	hetrobust::Bool,
	feid::Vector{Int64},
	fedfadj::Bool,
	obswt::Vector{T},
	fweights::Bool,
	maxmatsize::Float16,
	ptype::PType,
	bootstrapc::Bool,
	LIML::Bool,
	Fuller::T,
	kappa::T,
	ARubin::Bool,
	small::Bool,
	scorebs::Bool,
	reps::Int64,
	imposenull::Bool,
	auxwttype::AuxWtType,
	rng::AbstractRNG,
	level::T,
	rtol::T,
	madjtype::MAdjType,
	NH0::Int16,
	ML::Bool,
	scores::Matrix{T},
	beta::Vector{T},
	A::Symmetric{T,Matrix{T}},
	gridmin::Vector{T},
	gridmax::Vector{T},
	gridpoints::Vector{Float32},
	diststat::DistStatType,
	getCI::Bool,
	getplot::Bool,
	getauxweights::Bool,
	turbo::Bool) where T

	M = StrBootTest{T}(R, r, R1, r1, resp, predexog, predendog, inst, obswt, fweights, LIML, Fuller, kappa, ARubin,
	                   reps, auxwttype, rng, maxmatsize, ptype, imposenull, scorebs, !bootstrapc, clustid, nbootclustvar, nerrclustvar, issorted, hetrobust, small,
	                   feid, fedfadj, level, rtol, madjtype, NH0, ML, beta, A, scores, getplot,
	                   gridmin, gridmax, gridpoints, turbo)

	if getplot || (level<1 && getCI)
		plot!(M)
		plot = getplot ? (X=Tuple(M.plotX), p=M.plotY) : nothing
		peak = M.peak
		CI = level<1 & getCI ? M.CI : nothing
	else
		CI = plot = peak = nothing
	end
	
	padj = getp(M)  # trigger central (re)computation

	BootTestResult{T}(getstat(M),
	                  isone(nrows(R)) ? (small ? "t" : "z") : (small ? "F" : "χ²"),
	                  M.p, padj, M.B, M.BFeas, M.N✻, M.dof, M.dof_r, plot, peak, CI,
	                  getdist(M, diststat),
	                  getb(M), getV(M),
	                  getauxweights && reps>0 ? getv(M) : nothing #=, M=#)
end

vecconvert(T::DataType, X) = Vector(isa(X, AbstractArray) ? vec(    eltype(X)==T ? X : T.(X)) : X)
matconvert(T::DataType, X) = Matrix(isa(X, AbstractArray) ? reshape(eltype(X)==T ? X : T.(X), size(X,1), size(X,2)) : X)

function _wildboottest(T::DataType,
					  R::AbstractVecOrMat,
						r::AbstractVecOrMat;
					  resp::AbstractVector{<:Real},
					  predexog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  predendog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  inst::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  R1::AbstractVecOrMat=zeros(T,0,0),
						r1::AbstractVector=zeros(T,0),
					  clustid::AbstractVecOrMat{<:Integer}=zeros(Int,0,0),  # bootstrap-only clust vars, then boot&err clust vars, then err-only clust vars
					  nbootclustvar::Integer=1,
					  nerrclustvar::Integer=nbootclustvar,
						issorted::Bool=false,
					  hetrobust::Bool=true,
					  feid::AbstractVector{<:Integer}=Int8[],
					  fedfadj::Bool=true,
					  obswt::AbstractVector{<:Real}=T[],
					  fweights::Bool=false,
					  maxmatsize::Number=0,
					  ptype::PType=symmetric,
					  bootstrapc::Bool=false,
					  LIML::Bool=false,
					  Fuller::Number=0,
					  kappa::Number=NaN,
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
					  NH0::Integer=1,
					  ML::Bool=false,
					  scores::AbstractVecOrMat=Matrix{Float32}(undef,0,0),
					  beta::AbstractVector=T[],
					  A::AbstractMatrix=zeros(T,0,0),
					  gridmin::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = T[],
					  gridmax::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = T[],
					  gridpoints::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Integer = Int32[],
					  diststat::DistStatType=nodist,
					  getCI::Bool=true,
					  getplot::Bool=getCI,
					  getauxweights::Bool=false,
						turbo::Bool=false)

	nrows(R)>2 && (getplot = getCI = false)

  @assert length(predexog)==0 || nrows(predexog)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(predendog)==0 || nrows(predendog)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(inst)==0 || nrows(inst)==nrows(resp) "All data vectors/matrices must have same height"
  @assert length(feid)==0 || nrows(feid)==nrows(resp) "feid vector must have same height as data matrices"
  @assert length(clustid)==0 || nrows(clustid)==nrows(resp) "clustid must have same height as data matrices"
  @assert nrows(obswt)==0 || nrows(obswt)==nrows(resp) "obswt must have same height as data matrices"
  @assert nrows(R)==nrows(r) "R and r must have same height"
  @assert ncols(R)==ncols(predexog)+ncols(predendog) && isone(ncols(r)) "Wrong number of columns in null specification"
  @assert nrows(R1)==nrows(r1) "R₁ and r₁ must have same height"
  @assert length(R1)==0 || ncols(R1)==ncols(predexog)+ncols(predendog) "Wrong number of columns in model constraint specification"
	@assert ncols(r)==1 "r should have one column"
	@assert length(R1)==0 || ncols(r1)==1 "r1 should have one column"
  @assert nbootclustvar ≤ ncols(clustid) "nbootclustvar > width of clustid"
  @assert nerrclustvar ≤ ncols(clustid) "nerrclustvar > width of clustid"
  @assert reps ≥ 0 "reps < 0"
  @assert level ≥ 0. && level≤1. "level must be in the range [0,1]"
  @assert rtol > 0. "rtol ≤ 0"
  @assert NH0 > 0 "NH0 ≤ 0"
	if getplot || getCI
		@assert iszero(length(gridmin   )) || length(gridmin   )==nrows(R) "Length of gridmin doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridmax   )) || length(gridmax   )==nrows(R) "Length of gridmax doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridpoints)) || length(gridpoints)==nrows(R) "Length of gridpoints doesn't match number of hypotheses being jointly tested"
	end

	_gridmin = Vector{T}(undef, length(gridmin))
	_gridmax = Vector{T}(undef, length(gridmax))
	_gridpoints = Vector{Float32}(undef, length(gridpoints))
	for i ∈ 1:length(gridmin)  # cumbersome loops because map() and list comprehensions mess up type inference(?!)
		_gridmin[i] = T(ismissing(gridmin[i]) ? NaN : gridmin[i])
	end
	for i ∈ 1:length(gridmax)
		_gridmax[i] = T(ismissing(gridmax[i]) ? NaN : gridmax[i])
	end
	for i ∈ 1:length(gridpoints)
		_gridpoints[i] = T(ismissing(gridpoints[i]) ? NaN : gridpoints[i])
	end

	__wildboottest(
		matconvert(T,R),
		vecconvert(T,r);
		resp=vecconvert(T,resp),
		predexog=matconvert(T,predexog),
		predendog=matconvert(T,predendog),
		inst=matconvert(T,inst),
		R1=matconvert(T,R1),
		r1=vecconvert(T,r1),
		clustid=matconvert(Int64,clustid),
		nbootclustvar=Int8(nbootclustvar),
		nerrclustvar=Int8(nerrclustvar),
		issorted,
		hetrobust,
		feid=vecconvert(Int64,feid),
		fedfadj,
		obswt=vecconvert(T,obswt),
		fweights,
		maxmatsize=Float16(maxmatsize),
		ptype,
		bootstrapc,
		LIML,
		Fuller=T(Fuller),
		kappa=T(kappa),
		ARubin,
		small,
		scorebs,
		reps=Int64(reps),
		imposenull,
		auxwttype,
		rng,
		level=T(level),
		rtol=T(rtol),
		madjtype,
		NH0=Int16(NH0),
		ML,
		scores=matconvert(T,scores),
		beta=vecconvert(T,beta),
		A=Symmetric(matconvert(T,A)),
		gridmin=_gridmin,
		gridmax=_gridmax,
		gridpoints=_gridpoints,
		diststat,
		getCI,
		getplot,
		getauxweights,
		turbo)
end

_wildboottest(T::DataType, R, r::Number; kwargs...) = _wildboottest(T, R, [r]; kwargs...)
_wildboottest(T::DataType, R::UniformScaling{Bool}, r; kwargs...) = _wildboottest(T, diagm(fill(T(R.λ),nrows(r))), r; kwargs...)
_wildboottest(R::Union{UniformScaling{Bool},AbstractVecOrMat}, r::Union{Number,AbstractVecOrMat}; kwargs...) = _wildboottest(Float32, R, r; kwargs...)


"""
wildboottest([T::DataType=Float32,] R::AbstractMatrix, r::AbstractVector; 
             resp, <optional keyword arguments>) -> WildBootTests.BootTestResult

Function to perform wild-bootstrap-based hypothesis test

# Positional arguments
* `T::DataType`: data type for inputs, results, and computations: Float32 (default) or Float64
* `R::AbstractMatrix` and `r::AbstractVector`: required matrix and vector expressing the null Rβ=r; see notes below

# Required keyword argument
* `resp::AbstractVector`: response/dependent variable (y or y₁ in Roodman et al. (2019))

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
* `fedfadj::Bool=true`: true if small-sample adjustment should reflect number of fixed effects (if any)
* `obswt::AbstractVector=[]`: observation weight vector; default is equal weighting
* `fweights::Bool=false`: true for frequency weights
* `maxmatsize::Number`: maximum size of auxilliary weight matrix (v), in gigabytes
* `ptype::PType=symmetric`: p value type (`symmetric`, `equaltail`, `lower`, `upper`)
* `bootstrapc::Bool=false`: true to request bootstrap-c instead of bootstrap-t
* `LIML::Bool=false`: true for LIML or Fuller LIML
* `Fuller::Number`: Fuller LIML factor
* `kappa::Number`: fixed κ for _k_-class estimation
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
* `NH0::Integer=1`: number of hypotheses tested, including one being tested now
* `ML::Bool=false`: true for (nonlinear) ML estimation
* `scores::AbstractVecOrMat`: for ML, pre-computed scores
* `beta::AbstractVector`: for ML, parameter estimates
* `A::AbstractMatrix`: for ML, covariance estimates
* `gridmin`: vector of graph lower bounds; max length 2, `missing`/`NaN` entries ask wildboottest() to choose
* `gridmax`: vector of graph upper bounds; `missing`/`NaN` entries ask wildboottest() to choose
* `gridpoints`: vector of number of sampling points; `missing`/`NaN` entries ask wildboottest() to choose
* `diststat::DistStatType=nodiststat`: t to save bootstrap distribution of t/z/F/χ² statistics; numer to save numerators thereof
* `getCI::Bool=true`: whether to return CI
* `getplot::Bool=getCI`: whether to generate plot data
* `getauxweights::Bool=false`: whether to save auxilliary weight matrix (v)
* `turbo::Bool=false`: whether to exploit acceleration of the LoopVectorization package: slower on first use in a session, faster after

# Notes
`T`, `ptype`, `auxwttype`, `madjtype`, and `diststat` may also be strings. Examples: `"Float32"` and `"webb"`.

The columns of `R` in the statement of the null should correspond to those of the matrix [`predexog` `predendog`],
where `predendog` is non-empty only in regressions with instruments. 

Order the columns of `clustid` this way:
1. Variables only used to define bootstrapping clusters, as in the subcluster bootstrap.
2. Variables used to define both bootstrapping and error clusters.
3. Variables only used to define error clusters.
`nbootclustvar` is then the number of columns of type 1 or 2; `nerrclustvar` is the number of columns of type 2 or 3. Typically `clustid` is a single column of type 2. 

`wildboottest()` does not handle missing data values: all data and identifier matrices must 
be restricted to the estimation sample.

"""
wildboottest(   R, r; args...) = _wildboottest(                                                  R, r; Dict(a.first => (isa(a.second, AbstractString) ? eval(Meta.parse(a.second)) : a.second) for a ∈ args)...)
wildboottest(T, R, r; args...) = _wildboottest(isa(T, AbstractString) ? eval(Meta.parse(T)) : T, R, r; Dict(a.first => (isa(a.second, AbstractString) ? eval(Meta.parse(a.second)) : a.second) for a ∈ args)...)