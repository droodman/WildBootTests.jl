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
  ci::Union{Nothing, Matrix{T}}
  numerdist::Matrix{T}
  dist::Matrix{T}
  b::Vector{T}
  V::Matrix{T}
  auxweights::Union{Nothing,Matrix{T}}
  # o::StrBootTest
end

"""
`teststat(::WildBootTests.BootTestResult{T}) -> T`

Given a wildboottest() return object, extract test statistic
"""
teststat(o::BootTestResult) = o.stat

"""
`stattype(::WildBootTests.BootTestResult{T}) -> String`

Given a wildboottest() return object, extract type of test statistic: "t", "z", "F", or "χ²"
"""
stattype(o::BootTestResult) = o.stattype

"""
`statnumer(::WildBootTests.BootTestResult{T}) -> T`

Given a wildboottest() return object, extract numerator of test statistic
"""
statnumer(o::BootTestResult) = o.b

"""
`statvar(::WildBootTests.BootTestResult{T}) -> T`

Given a wildboottest() return object, extract squared denominator of test statistic
"""
statvar(o::BootTestResult) = o.V

"""
`numerdist(::WildBootTests.BootTestResult{T}) -> Matrix{T}`

Given a wildboottest() return object, extract bootstrap distribution of numerator of statistic
"""
numerdist(o::BootTestResult) = o.numerdist

"""
`dist(::WildBootTests.BootTestResult{T}) -> Matrix{T}`

Given a wildboottest() return object, extract bootstrap distribution of statistic
"""
dist(o::BootTestResult) = o.dist

"""
`p(::WildBootTests.BootTestResult{T}) -> T`

Given a wildboottest() return object, extract p value
"""
p(o::BootTestResult) = o.p

"""
`padj(::WildBootTests.BootTestResult{T}) -> T`

Given a wildboottest() return object, extract p value after multiple-hypothesis adjustment, if any
"""
padj(o::BootTestResult) = o.padj

"""
`reps(::WildBootTests.BootTestResult{T}) -> Int64`

Given a wildboottest() return object, extract number of replications
"""
reps(o::BootTestResult) = o.reps

"""
`repsfeas(::WildBootTests.BootTestResult{T}) -> Int64`

Given a wildboottest() return object, extract actual number of replications, subject to enumeration of Rademacher draws
"""
repsfeas(o::BootTestResult) = o.repsfeas

"""
`nbootclust(::WildBootTests.BootTestResult{T}) -> Int64`

Given a wildboottest() return object, extract number of bootstrapping clusters in test
"""
nbootclust(o::BootTestResult) = o.nbootclust

"""
`dof(::WildBootTests.BootTestResult{T}) -> Int64`

Given a wildboottest() return object, extract degrees of freedom of test
"""
dof(o::BootTestResult) = o.dof

"""
`dof_r(::WildBootTests.BootTestResult{T}) -> Int64`

Given a wildboottest() return object, extract residual degrees of freedom of test
"""
dof_r(o::BootTestResult) = o.dof_r

"""
`plotpoints(::WildBootTests.BootTestResult{T}) -> NamedTuple{(:X, :p), Tuple{Tuple{Vararg{Vector{T}, N} where N},Vector{T}}}`

Return data for confidence plot of test.
Return value is a 2-tuple with named entries `X` and `p` holding
the confidence sampling locations and p values respectively. `X` is in turn
a 1- or 2-tuple of vectors of sampling coordinates for each 
dimension of the tested hypothesis.
"""
plotpoints(o::BootTestResult) = o.plot

"""
`peak(::WildBootTests.BootTestResult{T}) -> NamedTuple{(:X, :p), Tuple{Vector{T}, T}}`

Given a wildboottest() return object for a one-dimensional test, return the parameter value with peak p value in test
Return value is a 2-tuple with named entries `X` and `p` holding the parameter value and p value.
"""
peak(o::BootTestResult) = o.peak

"""
`ci(::WildBootTests.BootTestResult{T}) -> Matrix{T}`

Given a wildboottest() return object for a one-dimensional test, extract the confidence interval(s)
for test, one row per disjoint piece
"""
ci(o::BootTestResult) = o.ci

"""
`auxweights(::WildBootTests.BootTestResult{T}) -> Matrix{T}`

Given a wildboottest() return object for a one-dimensional test, extract auxilliary weight matrix
"""
auxweights(o::BootTestResult) = o.auxweights

strint(x) = iszero(mod(x,1)) ? "$(Int64(x))" : "$x"
function Base.show(io::IO, o::BootTestResult{T}) where T
	s = stattype(o) * ( iszero(dof_r(o)) ? isone(dof(o)) ? " " : "(" * strint(dof(o)) * ")" :                                                     # z, χ²
	                                       isone(dof(o)) ? "(" * strint(dof_r(o)) * ")" : "(" * strint(dof(o)) * ", " * strint(dof_r(o)) * ")" )  # t, F
	Printf.@printf(io, "%s = %6.4f\n", s, teststat(o))
	Printf.@printf(io, "p%s = %6.4f\n", repeat(" ", length(s)-1), p(o))
	isdefined(o, :ci) && !isnothing(o.ci) && length(o.ci)>0 && print(io, "CI" * repeat(" ", length(s)-2) * " = $(round.(ci(o); sigdigits=4))\n")
end

# single entry point with arguments already converted to standardized types, to allow a smaller set of precompile() calls(?)
function __wildboottest(
	R::Matrix{T},
	r::Vector{T};
	resp::VecOrMat{T},
	predexog::Matrix{T},
	predendog::Matrix{T},
	inst::Matrix{T},
	R1::Matrix{T},
	r1::Vector{T},
	clustid::Matrix{Int64},
	nbootclustvar::Int64,
	nerrclustvar::Int64,
	issorted::Bool,
	hetrobust::Bool,
	feid::VecOrMat{Int64},
	fedfadj::Int64,
	obswt::VecOrMat{T},
	fweights::Bool,
	maxmatsize::Float16,
	ptype::Symbol,
	bootstrapc::Bool,
	liml::Bool,
	fuller::T,
	kappa::T,
	arubin::Bool,
	small::Bool,
	clusteradj::Bool,
	clustermin::Bool,
	jk::Bool,
	scorebs::Bool,
	reps::Int64,
	imposenull::Bool,
	auxwttype::Symbol,
	rng::AbstractRNG,
	level::T,
	rtol::T,
	madjtype::Symbol,
	nH0::Int16,
	ml::Bool,
	scores::Matrix{T},
	beta::Vector{T},
	A::Matrix{T},
	gridmin::VecOrMat{T},
	gridmax::VecOrMat{T},
	gridpoints::VecOrMat{T},
	getdist::Bool,
	getci::Bool,
	getplot::Bool,
	getauxweights::Bool,
	overwrite::Bool,
	v::Matrix{T}) where T

	o = StrBootTest{T}(R, r, R1, r1, resp, iszero(length(predexog)) ? Matrix{T}(undef,nrows(resp),0) : predexog, 
	                   predendog, inst, obswt, fweights, liml, fuller, kappa, arubin,
	                   reps, auxwttype, rng, maxmatsize, ptype, imposenull, jk, scorebs, !bootstrapc, clustid, nbootclustvar, nerrclustvar, issorted, hetrobust, small, clusteradj, clustermin,
	                   feid, fedfadj, level, rtol, madjtype, nH0, ml, beta, A, scores, getplot,
	                   gridmin, gridmax, gridpoints, overwrite, v)

	o.getci = getplot || (level<1 && getci)
	if o.getci
		plot!(o)
		plot = getplot & isdefined(o, :plotX) ? (X=Tuple(o.plotX), p=o.plotY) : nothing
		peak = o.peak
		ci = level<1 & getci ? o.ci : nothing
	else
		ci = plot = peak = nothing
		o.boottest!(o)
	end
	
	padj = getp(o)

	BootTestResult{T}(getstat(o),
	                  isone(nrows(R)) ? (small ? "t" : "z") : (small ? "F" : "χ²"),
	                  o.p, padj, o.B, o.BFeas, o.N✻, o.dof, o.dof_r, plot, peak, ci,
	                  getdists(o, getdist)...,
	                  getb(o), getV(o),
	                  getauxweights && reps>0 ? getv(o) : nothing #=, o=#)
end

function vecconvert(T::DataType, X)
	t = X isa SharedArray ? sdata(X) : X
	if iszero(length(t))
		T[]
	elseif t isa Vector{T}
		t
	elseif t isa AbstractVector{T}
		Vector{T}(t)
	elseif eltype(t)==T
		reshape(t,size(t,1))
	else
		Vector{T}(reshape(t, size(t,1)))
	end
end

function matconvert(T::DataType, X)
	t = X isa SharedArray ? sdata(X) : X
	if t isa Matrix{T}
		t
	elseif t isa AbstractMatrix{T}
		Matrix{T}(t)
	elseif eltype(t)==T
		reshape(t, size(t,1), size(t,2))
	else
		Matrix{T}(reshape(t, size(t,1), size(t,2)))
	end
end

@nospecialize

function _wildboottest(T::DataType,
					  R::AbstractVecOrMat,
						r::AbstractVecOrMat;
					  resp::AbstractVecOrMat{<:Real},
					  predexog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  predendog::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  inst::AbstractVecOrMat{<:Real}=zeros(T,0,0),
					  R1::AbstractVecOrMat=zeros(T,0,0),
						r1::AbstractVecOrMat=zeros(T,0),
					  hetrobust::Bool=true,
					  clustid::AbstractVecOrMat{<:Integer}=zeros(Int,0,0),  # bootstrap-only clust vars, then boot&err clust vars, then err-only clust vars
					  nbootclustvar::Integer=ncols(clustid),
					  nerrclustvar::Integer=nbootclustvar,
						issorted::Bool=false,
					  feid::AbstractVecOrMat{<:Integer}=Int8[],
					  fedfadj::Integer=length(feid)>0 ? -1 : 0,
					  obswt::AbstractVecOrMat{<:Real}=T[],
					  fweights::Bool=false,
					  maxmatsize::Number=0,
					  ptype::Symbol=:symmetric,
					  bootstrapc::Bool=false,
					  liml::Bool=false,
					  fuller::Number=0,
					  kappa::Number=NaN,
					  arubin::Bool=false,
					  small::Bool=true,
						clusteradj::Bool=small,
						clustermin::Bool=false,
					  jk::Bool=false,
					  scorebs::Bool=false,
					  reps::Integer=999,
					  imposenull::Bool=true,
					  auxwttype::Symbol=:rademacher,
					  rng::AbstractRNG=Xoshiro(),
					  level::Number=.95,
					  rtol::Number=1e-3,
					  madjtype::Symbol=:none,
					  nH0::Integer=1,
					  ml::Bool=false,
					  scores::AbstractVecOrMat=Matrix{Float64}(undef,0,0),
					  beta::AbstractVecOrMat=T[],
					  A::AbstractMatrix=zeros(T,0,0),
					  gridmin::Union{VecOrMat{S},VecOrMat{Union{S,Missing}}} where S<:Number = T[],
					  gridmax::Union{VecOrMat{S},VecOrMat{Union{S,Missing}}} where S<:Number = T[],
					  gridpoints::Union{VecOrMat{S},VecOrMat{Union{S,Missing}}} where S<:Number = Int64[],
					  getdist::Bool=true,
					  getci::Bool=true,
					  getplot::Bool=getci,
					  getauxweights::Bool=false,
						overwrite::Bool=false,
						v::AbstractVecOrMat=T[;;])

	nrows(R)>2 && (getplot = getci = false)

	@assert any(auxwttype .== (:rademacher, :mammen, :webb, :gamma, :normal)) "auxwttype shoud be :rademacher, :mammen, :webb, :gamma, or :normal"
	@assert any(ptype .==(:symmetric, :equaltail, :lower, :upper)) "ptype should be :symmetric, :equaltail, :lower, or :upper"
	@assert any(madjtype .== (:none, :bonferroni, :sidak)) "madjtype should be :none, :bonferroni, or :sidak"
	@assert ncols(obswt)≤1 "obswt must have one column"
  @assert nrows(R)==nrows(r) "R and r must have same height"
  @assert (ncols(R) == (ml ? nrows(beta) : ncols(predexog)+ncols(predendog)) && isone(ncols(r))) "Wrong number of columns in null specification"
  @assert length(R1)==0 || ncols(R1)==ncols(predexog)+ncols(predendog) "Wrong number of columns in model constraint specification"
	@assert ncols(r)==1 "r should have one column"
	@assert length(R1)==0 || ncols(r1)==1 "r1 should have one column"
  @assert nbootclustvar ≤ ncols(clustid) "nbootclustvar > width of clustid"
  @assert nerrclustvar ≤ ncols(clustid) "nerrclustvar > width of clustid"
  @assert reps ≥ 0 "reps < 0"
  @assert level ≥ 0. && level≤1. "level must be in the range [0,1]"
  @assert rtol > 0. "rtol ≤ 0"
  @assert nH0 > 0 "nH0 ≤ 0"
	@assert !liml || (ncols(predendog)>0 && ncols(inst)>0) "For liml, non-empty predendog and inst arguments are needed"
	@assert fuller==0 || (ncols(predendog)>0 && ncols(inst)>0) "For Fuller liml, non-empty predendog and inst arguments are needed"
	
	if !ml 
		@assert ncols(resp)==1 "resp should have one column"
	  @assert (length(predexog)==0 || nrows(predexog)==nrows(resp)) && 
		        (length(predendog)==0 || nrows(predendog)==nrows(resp)) &&
						(length(inst)==0 || nrows(inst)==nrows(resp)) "All data vectors/matrices must have same height"
						@assert ncols(inst) >= ncols(predendog) "Model has fewer instruments than instrumented variables"
						@assert length(feid)==0 || nrows(feid)==nrows(resp) "feid vector must have same height as data matrices"
						@assert ncols(feid)≤1 "feid should have one column"
						@assert length(clustid)==0 || nrows(clustid)==nrows(resp) "clustid must have same height as data matrices"
						@assert nrows(obswt)==0 || nrows(obswt)==nrows(resp) "obswt must have same height as data matrices"
						@assert nrows(R1)==nrows(r1) "R₁ and r₁ must have same height"
						@assert iszero(ncols(predendog)) || ncols(inst)>0 "predendog provided without inst"
						@assert !arubin || ncols(predendog)>0 "Anderson-Rubin test requested but predendog not provided"
	end

	if getplot || getci
		@assert iszero(length(gridmin   )) || length(gridmin   )==nrows(R) "Length of gridmin doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridmax   )) || length(gridmax   )==nrows(R) "Length of gridmax doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridpoints)) || length(gridpoints)==nrows(R) "Length of gridpoints doesn't match number of hypotheses being jointly tested"
		@assert iszero(length(gridmin   )) || ncols(gridmin)   ==1 "gridmin should have one column"
		@assert iszero(length(gridmax   )) || ncols(gridmax)   ==1 "gridmax should have one column"
		@assert iszero(length(gridpoints)) || ncols(gridpoints)==1 "gridpoints should have one column"
	end

	_gridmin = Vector{T}(undef, length(gridmin))
	_gridmax = Vector{T}(undef, length(gridmax))
	_gridpoints = Vector{T}(undef, length(gridpoints))
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
		nbootclustvar=Int64(nbootclustvar),
		nerrclustvar=Int64(nerrclustvar),
		issorted,
		hetrobust,
		feid=vecconvert(Int64,feid),
		fedfadj,
		obswt=vecconvert(T,obswt),
		fweights,
		maxmatsize=Float16(maxmatsize),
		ptype,
		bootstrapc,
		liml,
		fuller=T(fuller),
		kappa=T(kappa),
		arubin,
		small,
		clusteradj,
		clustermin,
		jk,
		scorebs,
		reps=Int64(reps),
		imposenull,
		auxwttype,
		rng,
		level=T(level),
		rtol=T(rtol),
		madjtype,
		nH0=Int16(nH0),
		ml,
		scores=matconvert(T,scores),
		beta=vecconvert(T,beta),
		A=matconvert(T,A),
		gridmin=_gridmin,
		gridmax=_gridmax,
		gridpoints=_gridpoints,
		getdist,
		getci,
		getplot,
		getauxweights,
		overwrite,
		v=matconvert(T,v))
end

_wildboottest(T::DataType, R, r::Number; kwargs...) = _wildboottest(T, R, [r]; kwargs...)
_wildboottest(T::DataType, R::UniformScaling{Bool}, r; kwargs...) = _wildboottest(T, diagm(fill(T(R.λ),nrows(r))), r; kwargs...)
_wildboottest(R::Union{UniformScaling{Bool},AbstractVecOrMat}, r::Union{Number,AbstractVecOrMat}; kwargs...) = _wildboottest(Float64, R, r; kwargs...)


"""
`wildboottest([T::DataType=Float64,] R::AbstractMatrix, r::AbstractVector; 
             resp, <optional keyword arguments>) -> WildBootTests.BootTestResult`

Function to perform wild-bootstrap-based hypothesis test

# Positional arguments
* `T::DataType`: data type for inputs, results, and computations: Float32 or Float64 (default)
* `R::AbstractMatrix` and `r::AbstractVector`: required matrix and vector expressing the null Rβ=r; see notes below

# Required keyword argument
* `resp::AbstractVector`: response/dependent variable (y or y₁ in Roodman et al. (2019))

# Optional keyword arguments
* `predexog::AbstractVecOrMat`: exogenous predictors, including constant term, if any (X/X₁)
* `predendog::AbstractVecOrMat`: endogenous predictors (Y₂)
* `inst::AbstractVecOrMat`: instruments (X₂)
* `R1::AbstractMatrix` and `r1::AbstractVector`: model constraints; same format as for `R` and `r`
* `clustid::AbstractVecOrMat{<:Integer}`: data vector/matrix of error and bootstrapping cluster identifiers; see notes 
* `nbootclustvar::Integer=size(clustid,2)`: number of bootstrap-clustering variables
* `nerrclustvar::Integer=nbootclustvar`: number of error-clustering variables
* `issorted:Bool=false`: time-saving flag: data matrices	 are already sorted by column types 2, then 3, then 1 (see notes)
* `hetrobust::Bool=true`: true unless errors are treated as iid
* `feid::AbstractVector{<:Integer}`: data vector for one-way fixed effect group identifier
* `fedfadj::Integer`: degrees of freedom that fixed effects (if any) consume; defaults to number of FEs
* `obswt::AbstractVector=[]`: observation weight vector; default is equal weighting
* `fweights::Bool=false`: true for frequency weights
* `maxmatsize::Number`: maximum size of auxilliary weight matrix (v), in gigabytes
* `ptype::Symbol=:symmetric`: p value type (`:symmetric`, `:equaltail`, `:lower`, `:upper`)
* `bootstrapc::Bool=false`: true to request bootstrap-c instead of bootstrap-t
* `liml::Bool=false`: true for LIML or Fuller LIML
* `fuller::Number`: Fuller LIML factor
* `kappa::Number`: fixed κ for _k_-class estimation
* `arubin::Bool=false`: true for Anderson-Rubin test
* `small::Bool=true`: true to multiply test statistics by G/(G-1) × N/(N-k), where G, N, k are number of clusters, observations, and predictors
* `clusteradj::Bool=true`: false to drop G/(G-1) factor
* `clustermin::Bool=false``: for multiway clustering, true to base G/(G-1) factor for all clusterings ]on the smallest G across clusterings
* `jk::Bool=false`: true to base the bootstrap data-generating process on residuals jackknifed by bootstrap cluster
* `scorebs::Bool=false`: true for score bootstrap instead of wild bootstrap
* `reps::Integer=999`: number of bootstrap replications; `reps` = 0 requests classical Rao (or Wald) test if `imposenull` = `true` (or `false`)
* `imposenull::Bool=true`: true to impose null
* `auxwttype::Symbol=:rademacher`: auxilliary weight type (`:rademacher`, `:mammen`, `:webb`, `:normal`, `:gamma`)
* `rng::AbstractRNG=MersenneTwister()`: randon number generator
* `level::Number=.95`: significance level (0-1)
* `rtol::Number=1e-3`: tolerance for confidence set bound determination
* `madjtype::Symbol=:none`: multiple hypothesis adjustment (`:none`, `:bonferroni`, `:sidak`)
* `nH0::Integer=1`: number of hypotheses tested, including one being tested now
* `ml::Bool=false`: true for (nonlinear) ML estimation
* `scores::AbstractVecOrMat`: for ML, pre-computed scores
* `beta::AbstractVector`: for ML, parameter estimates
* `A::AbstractMatrix`: for ML, covariance estimates
* `gridmin`: vector of graph lower bounds; max length 2, `missing`/`NaN` entries ask wildboottest() to choose
* `gridmax`: vector of graph upper bounds; `missing`/`NaN` entries ask wildboottest() to choose
* `gridpoints`: vector of number of sampling points; `missing`/`NaN` entries ask wildboottest() to choose
* `getdist::Bool=:false`: whether to return bootstrapped distribution for t/z/F/χ² statistics; and their numerators
* `getci::Bool=true`: whether to return confidence interval
* `getplot::Bool=getci`: whether to generate plot data
* `getauxweights::Bool=false`: whether to save auxilliary weight matrix (v)

# Notes
`T`, `ptype`, `auxwttype`, and `madjtype` may also be strings. Examples: `"Float32"` and `"webb"`.

The columns of `R` in the statement of the null should correspond to those of the matrix [`predexog` `predendog`],
where `predendog` is non-empty only in regressions with instruments. 

Order the columns of `clustid` this way:
1. Variables only used to define bootstrapping clusters, as in the subcluster bootstrap.
2. Variables used to define both bootstrapping and error clusters.
3. Variables only used to define error clusters.
`nbootclustvar` is then the number of columns of type 1 or 2; `nerrclustvar` is the number of columns of type 2 or 3. Typically `clustid` is a single column of type 2 and `nbootclustvar` and `nerrclustvar` default to 1.

`wildboottest()` does not handle missing data values: all data and identifier matrices must 
be restricted to the estimation sample.

"""
wildboottest(   R, r; kwargs...) = _wildboottest(                                              R, r; Dict(a.first => a.second isa AbstractString ? Symbol(a.second) : a.second for a ∈ kwargs)...)  # convert any string parameter to symbols
wildboottest(T, R, r; kwargs...) = _wildboottest(isa(T, AbstractString) ? eval(Symbol(T)) : T, R, r; Dict(a.first => a.second isa AbstractString ? Symbol(a.second) : a.second for a ∈ kwargs)...)

wildboottest!(   R, r; kwargs...) = _wildboottest(                                              R, r; overwrite=true, Dict(a.first => a.second isa AbstractString ? Symbol(a.second) : a.second for a ∈ kwargs if a.first ≠ :overwrite)...)  # convert any string parameter to symbols
wildboottest!(T, R, r; kwargs...) = _wildboottest(isa(T, AbstractString) ? eval(Symbol(T)) : T, R, r; overwrite=true, Dict(a.first => a.second isa AbstractString ? Symbol(a.second) : a.second for a ∈ kwargs if a.first ≠ :overwrite)...)

const skipargs = (:imposenull, :reps, :scorebs)
waldtest(args...; kwargs...) = wildboottest(args...; reps=0, imposenull=false, scorebs=true, Dict(a.first => a.second for a ∈ kwargs if a.first ∉ skipargs)...)
scoretest(args...; kwargs...) = wildboottest(args...; reps=0, imposenull=true, scorebs=true, Dict(a.first => a.second for a ∈ kwargs if a.first ∉ skipargs)...)
