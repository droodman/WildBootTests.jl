# Definition of StrBoottest "class" for holding intermediate results, with associated utilities and get functions

"Auxilliary weight types: `rademacher`, `mammen`, `webb`, `normal`, `gamma`"
@enum AuxWtType rademacher mammen webb normal gamma

"p value types: `symmetric`, `equaltail`, `lower`, `upper`"
@enum PType symmetric equaltail lower upper

"Multiple hypothesis adjustment types: `nomadj`, `bonferroni`, `sidak`"
@enum MAdjType nomadj bonferroni sidak

"Bootstrap distribution statistics optionally returned"
@enum DistStatType nodist t numer

struct StrClust{T<:Real}
	N::Int; multiplier::T; even::Bool
	order::Vector{Int64}
	info::Vector{UnitRange{Int64}}
end

struct StrFE{T<:Real}
	is::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}
  wt::Vector{T}
end


@inline matconvert(T::DataType, X) = !isa(X, AbstractArray) || eltype(X)==T ? X : T.(X)

mutable struct StrBootTest{T<:AbstractFloat}
  R::Matrix{T}; r::Vector{T}; R₁::Matrix{T};r₁::Vector{T}
  y₁::Vector{T}; X₁::VecOrMat{T}; Y₂::VecOrMat{T}; X₂::VecOrMat{T}
  wt::Union{Vector{T}, UniformScaling}; fweights::Bool
  LIML::Bool; Fuller::T; κ::T; ARubin::Bool
  B::Int64; auxtwtype::AuxWtType; rng::AbstractRNG; maxmatsize::Float16
  ptype::PType; null::Bool; bootstrapt::Bool
	ID::Union{VecOrMat{Int8},VecOrMat{Int16},VecOrMat{Int32},VecOrMat{Int64}}; nbootclustvar::Int8; nerrclustvar::Int64; issorted::Bool; small::Bool
  FEID::Union{VecOrMat{Int8},VecOrMat{Int16},VecOrMat{Int32},VecOrMat{Int64}}; FEdfadj::Bool
  level::T; rtol::T
  madjtype::MAdjType; NH₀::Int16
  ML::Bool; β̂::Vector{T}; A::Matrix{T}; sc::VecOrMat{T}
  willplot::Bool; gridmin::Vector{Union{T,Missing}}; gridmax::Vector{Union{T,Missing}}; gridpoints::Vector{Union{Int32,Missing}}

  q::Int16; twotailed::Bool; scorebs::Bool; robust::Bool

  WRE::Bool; initialized::Bool; NFE::Int64; FEboot::Bool; granular::Bool; NErrClustCombs::Int16; subcluster::Int8; BFeas::Int64; interpolating::Bool
  dirty::Bool; v_sd::T; notplotted::Bool
  confpeak::Vector{T}
  IDBootData::Vector{Int64}; IDBootAll::Vector{Int64}
  anchor::Vector{T}; poles::Vector{T}; numer::Matrix{T}; dist::Vector{T}
  CI::Matrix{T}
  peak::NamedTuple{(:X, :p), Tuple{Vector{T}, T}}

  sqrt::Bool; Nobs::Int64; _Nobs::Int64; kZ::Int64; kY₂::Int64; kX₁::Int64; sumwt::T; NClustVar::Int8; haswt::Bool; REst::Bool; multiplier::T; smallsample::T
		WREnonARubin::Bool; dof::Int64; dof_r::Int64; p::T; BootClust::Int8
		purerobust::Bool; N✻::Int64; Nw::Int64; enumerate::Bool; interpolable::Bool; interpolate_u::Bool; kX₂::Int64; kX::Int64
  _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; u✻::Matrix{T}; CT_WE::Matrix{T}
  infoBootData::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; infoErrAll::Vector{UnitRange{Int64}}
  JN⋂N✻::Matrix{T}; statDenom::Matrix{T}; uXAR::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; β̂dev::Matrix{T}; δdenom_b::Matrix{T}
	_J⋂::Matrix{T}; YY✻_b::Matrix{T}; YPXY✻_b::Matrix{T}; numerw::Matrix{T}; Zyg::Vector{Matrix{T}}; numer_b::Vector{T}
		
	distCDR::Matrix{T}; plotX::Tuple{Vararg{Vector{T}, N} where N}; plotY::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  numersum::Vector{T}; ü₀::Vector{T}; invFEwt::Vector{T}
	β̂s::Matrix{T}; As::Matrix{T}
	infoAllData::Vector{UnitRange{Int64}}; info⋂Data::Vector{UnitRange{Int64}}; IDAll::Matrix{T}; Ü₂par::Matrix{T}
	ü::Vector{T}
	DGP::StrEstimator{T,E} where E; Repl::StrEstimator{T,E} where E; M::StrEstimator{T,E} where E
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; SCT⋂uXinvXX::Matrix{Matrix{T}}; S✻UU::Matrix{Vector{T}}; CTUX::Matrix{Matrix{T}}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; IDCT⋂✻::Vector{Vector{Int64}}; infoCT⋂✻::Vector{Vector{UnitRange{Int64}}}; S✻UX::Vector{Matrix{T}}; S✻UXinvXX::Vector{Matrix{T}}; S✻UZperpinvZperpZperp::Vector{Matrix{T}}; S✻uY::Vector{Matrix{T}}; S✻UMZperp::Vector{Matrix{T}}; S✻UPX::Vector{Matrix{T}}; S✻UZperp::Vector{Matrix{T}}; CTFEU::Vector{Matrix{T}}
  ∂denom∂r::Vector{Matrix{Matrix{T}}}; ∂Jcd∂r::Vector{Matrix{Matrix{T}}}
  ∂²denom∂r²::Matrix{Matrix{Matrix{T}}}
	FEs::Vector{StrFE{T}}
  T1L::Vector{Matrix{T}}; T1R::Vector{Matrix{T}}
	crosstab⋂✻ind::Vector{Int64}; crosstabBootind::Vector{Int64}
  seed::UInt64

  StrBootTest{T}(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, LIML, 
	               Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, scorebs, bootstrapt, ID, nbootclustvar, nerrclustvar, issorted, robust, small, FEID, FEdfadj, level, rtol, madjtype, NH₀, ML,
								 β̂, A, sc, willplot, gridmin, gridmax, gridpoints) where T<:Real =
	  new(matconvert(T,R), matconvert(T,r), matconvert(T,R₁), matconvert(T,r₁), matconvert(T,y₁), matconvert(T,X₁), matconvert(T,Y₂), matconvert(T,X₂), matconvert(T,wt), fweights, LIML || !iszero(Fuller), 
		    Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, bootstrapt, matconvert(Int64,ID), nbootclustvar, nerrclustvar, issorted, small, FEID, FEdfadj, level, rtol, madjtype, NH₀, ML, 
				matconvert(T,β̂), matconvert(T,A), matconvert(T,sc), willplot, gridmin, gridmax, gridpoints,
		  nrows(R),
		  ptype==symmetric || ptype==equaltail,
		  scorebs || iszero(B) || ML,
		  robust || nerrclustvar>0,
		  false, false, 0, false, false, 0, 0, 0, false,
		  true, one(T), true,
		  [T(0)],
		  Vector{Int64}(undef,0), Vector{Int64}(undef,0),
		  Vector{T}(undef,0), Vector{T}(undef,0), Matrix{T}(undef,0,0), Vector{T}(undef,0),
		  Matrix{T}(undef,0,0),
		  (X = Vector{T}(undef,0), p = T(NaN)))
end


# cross-tab sum of a column vector w.r.t. given panel info and fixed-effect var
# one row per FE, one col per other grouping
function crosstabFE(o::StrBootTest{T}, v::AbstractVector{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, o.NFE, nrows(info))
  if length(info)>0
		@inbounds for i ∈ 1:nrows(info)
	    FEIDi = @view o._FEID[info[i]]
	    vi    = @view       v[info[i]]
	    @inbounds @simd for j in eachindex(vi, FEIDi)
	  	  dest[FEIDi[j],i] += vi[j]
	    end
	  end
  else  # "robust" case, no clustering
	  @inbounds @simd for i in eachindex(v,o._FEID)
	    dest[o._FEID[i],i] = v[i]
	  end
  end
  dest
end
# same, transposed
function crosstabFEt(o::StrBootTest{T}, v::AbstractVector{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, nrows(info), o.NFE)
  if length(info)>0
		@inbounds for i ∈ 1:nrows(info)
	    FEIDi = @view o._FEID[info[i]]
	    vi    = @view       v[info[i]]
	    @simd for j ∈ eachindex(vi, FEIDi)
	  	  dest[i,FEIDi[j]] += vi[j]
	    end
	  end
  else  # "robust" case, no clustering
	  @inbounds @simd for i ∈ eachindex(v,o._FEID)
	    dest[i,o._FEID[i]] = v[i]
	  end
  end
  dest
end
# same, but handles multiple columns in v
# *2nd* dimension of resulting 3-D array corresponds to cols of v;
# this facilitates reshape() to 2-D array in which results for each col of v are stacked vertically
function crosstabFEt(o::StrBootTest{T}, v::AbstractMatrix{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, nrows(info), ncols(v), o.NFE)
  if length(info)>0
	  @inbounds for i ∈ 1:nrows(info)
	    FEIDi = @view o._FEID[info[i]]
	    vi    = @view       v[info[i],:]
	    @simd for j ∈ 1:length(FEIDi)
	  	  dest[i,:,FEIDi[j]] += @view vi[j,:]
	    end
	  end
  else  # "robust" case, no clustering
	  @inbounds @simd for i ∈ 1:length(FEIDi)
	    dest[i,:,o._FEID[i]] = @view v[i,:]
	  end
  end
  dest
end


# partial any fixed effects out of a data matrix
function partialFE!(o::StrBootTest, In::AbstractArray)
  if length(In)>0
	  for f ∈ o.FEs
	    tmp = @view In[f.is,:]
	    tmp .-= f.wt'tmp
    end
  end
end
function partialFE(o::StrBootTest, In::AbstractArray)
  Out = similar(In)
  if length(In)>0
	  for f ∈ o.FEs
	    tmp = @view In[f.is,:]
	    Out[f.is,:] .= tmp .- f.wt'tmp
	  end
  end
  Out
end


macro storeWtGrpResults!(dest, content)  # poor hygiene in referencing caller's o and w
  if dest == :(o.dist)
	return quote
	  if isone($(esc(:o)).Nw)
		$(esc(dest)) = $(esc(content))
	  else
		$(esc(dest))[$(esc(:o)).WeightGrp[$(esc(:w))]] = $(esc(content))
	  end
	end
  else
	  return quote
	    if isone($(esc(:o)).Nw)
	  	  $(esc(dest)) = $(esc(content))
	    else
	  	  $(esc(dest))[:,$(esc(:o)).WeightGrp[$(esc(:w))]] = $(esc(content))
	    end
	  end
  end
end

macro clustAccum!(X, c, Y)  # efficiently add a cluster combination-specific term, factoring in the needed multiplier and sign
  return quote
	  if isone($(esc(c)))
	    if isone($(esc(:o)).clust[1].multiplier)
	  	  $(esc(X)) = $(esc(:o)).clust[1].even ? $(esc(Y)) : -$(esc(Y))
	    else
	  	  $(esc(X)) = $(esc(Y)) * ($(esc(:o)).clust[1].even ? $(esc(:o)).clust[1].multiplier : -$(esc(:o)).clust[1].multiplier)
	    end
	  elseif $(esc(:o)).clust[$(esc(c))].even
	    if isone($(esc(:o)).clust[$(esc(c))].multiplier)
	  	  $(esc(X)) .+= $(esc(Y))
	    else
	  	  $(esc(X)) .+= $(esc(Y)) .* $(esc(:o)).clust[$(esc(c))].multiplier
	    end
	  elseif isone($(esc(:o)).clust[$(esc(c))].multiplier)
	    $(esc(X)) .-= $(esc(Y))
	  else
	    $(esc(X)) .-= $(esc(Y)) .* $(esc(:o)).clust[$(esc(c))].multiplier
	  end
  end
end

function getdist(o::StrBootTest, diststat::DistStatType=nodist)
  o.dirty && boottest!(o)
  if diststat == numer
	  _numer = isone(o.v_sd) ? o.numer : o.numer / o.v_sd
	  o.distCDR = (@view _numer[:,2:end]) .+ o.r
	  sort!(o.distCDR)
  elseif nrows(o.distCDR)==0  # return test stats
    if length(o.dist) > 1
	    o.distCDR = reshape((@view o.dist[2:end]), :, 1) * o.multiplier
	    sort!(o.distCDR, dims=1)
	  else
	    o.distCDR = zeros(0,1)
	  end
  end
  o.distCDR
end

function sumgreater(x, v)
  dest = zero(Int64)
  @inbounds @simd for i in v
	  x > i && (dest += 1)
  end
  dest
end
function sumless(x, v)
  dest = zero(Int64)
  @inbounds @simd for i in v
	  x < i && (dest += 1)
  end
  dest
end
function sumlessabs(x, v)
  dest = zero(Int64)
  @inbounds @simd for i in v
	  x < abs(i) && (dest += 1)
  end
  dest
end

# get p valuo. Robust to missing bootstrapped values interpreted as +infinity.
function getp(o::StrBootTest{T}; classical::Bool=false) where T
  o.dirty && boottest!(o)
  tmp = o.dist[1]
  isnan(tmp) && return tmp
  if o.B>0 && !classical
  	if o.sqrt && o.ptype ≠ upper
  	  if o.ptype==symmetric
  	    n = sumlessabs(abs(tmp), o.dist)   # symmetric p value; do so as not to count missing entries in *dist
  	  elseif o.ptype==equaltail
  	    n = 2min(sumgreater(tmp, o.dist) , sumless(tmp, o.dist))
  	  else
  		  n = sumgreater(tmp,  o.dist)  # lower-tailed p value
      end
  	else
  	  n = sumless(tmp, o.dist)  # upper-tailed p value or p value based on squared stats
    end
    o.p = n / o.BFeas |> T
  else
  	tmp *= o.multiplier
    _p = ccdf(o.small ? FDist(o.dof, o.dof_r) : Chisq(o.dof), Float64(o.sqrt ? tmp^2 : tmp))  |> T
  	if o.sqrt && !o.twotailed
  	  _p /= 2
  	  (ptype==upper) == (tmp<0) && (_p = 1 - _p)
  	end
    classical && return _p
    o.p = _p
  end
  o.p
end

# numerator for full-sample test stat
function getb(o::StrBootTest)
  o.dirty && boottest!(o)
  @views isone(o.v_sd) ? o.numer[:,1] : o.numer[:,1] / o.v_sd
end

# denominator for full-sample test stat
function getV(o::StrBootTest)
  o.dirty && boottest!(o)
  o.statDenom / ((isone(o.v_sd) ? o.smallsample : o.v_sd^2 * o.smallsample) * (o.sqrt ? o.multiplier^2 : o.multiplier) * o.dof)
end

# wild weights
getv(o::StrBootTest) = @views isone(o.v_sd) ? o.v[:,2:end] : o.v[:,2:end] / o.v_sd

# Return number of bootstrap replications with feasible results
# Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
getrepsfeas(o::StrBootTest) = o.BFeas
getNBootClust(o::StrBootTest) = o.N✻
getreps(o::StrBootTest) = o.B  # return number of replications, possibly reduced to 2^G

function getpadj(o::StrBootTest{T}; classical::Bool=false) where T
  _p = o.dirty || classical ? getp(o, classical=classical) : o.p
  if o.madjtype==bonferroni min(one(T), o.NH₀ * _p)
  elseif o.madjtype==sidak  one(T) - (one(T) - _p) ^ o.NH₀
  else _p
  end
end

function getstat(o::StrBootTest)
  o.dirty && boottest!(o)
  o.multiplier * o.dist[1]
end
function getdf(o::StrBootTest)
  o.dirty && boottest!(o)
  o.dof
end
function getdf_r(o::StrBootTest)
  o.dirty && boottest!(o)
  o.dof_r
end
function _getplot(o::StrBootTest)
  o.notplotted && plot(o)
  (X=o.plotX, p=o.plotY)
end
function getpeak(o::StrBootTest)  # x and y values of confidence curve peak (at least in OLS && ARubin)
  o.notplotted && plot(o)
  o.peak
end
function _getCI(o::StrBootTest)
  o.notplotted && plot(o)
  o.CI
end
