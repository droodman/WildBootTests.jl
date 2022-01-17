# Definition of StrBootTest "class" for holding intermediate results, with associated utilities and get functions

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

mutable struct StrEstimator{T<:AbstractFloat}
  isDGP::Bool; LIML::Bool; Fuller::T; κ::T
  R₁perp::Matrix{T}; Rpar::Matrix{T}

  kZ::Int64
  y₁::Vector{T}; ü₁::Vector{T}; u⃛₁::Vector{T}; β̈::Vector{T}; γ̈::Vector{T}; β̈₀::Vector{T}; invXXXy₁par::Vector{T}
  Yendog::Vector{Bool}
  invZperpZperp::Symmetric{T,Matrix{T}}; XZ::Matrix{T}; PXZ::Matrix{T}; YPXY::Symmetric{T,Matrix{T}}; R₁invR₁R₁::Matrix{T}
	restricted::Bool; RperpX::#=Uber=#Matrix{T}; RperpXperp::#=Uber=#Matrix{T}; RRpar::Matrix{T}; RparY::#=Uber=#Matrix{T}; RR₁invR₁R₁::Matrix{T}
	∂β̈∂r::Matrix{T}; YY::Symmetric{T,Matrix{T}}; AR::Matrix{T}; XAR::Matrix{T}; R₁invR₁R₁Y::Matrix{T}; invXXXZ::Matrix{T}; Ü₂::Matrix{T}; XinvXX::Matrix{T}; Rt₁::Vector{T}
	invXX::Symmetric{T,Matrix{T}}; Y₂::Matrix{T}; X₂::Matrix{T}; invH::Symmetric{T,Matrix{T}}
	y₁par::Vector{T}; Xy₁par::Vector{T}
	A::Symmetric{T,Matrix{T}}; Z::Matrix{T}; Zperp::Matrix{T}; X₁::Matrix{T}
	FillingT₀::Matrix{Matrix{T}}
	WXAR::Matrix{T}; S⋂PXYZperp::Vector{Matrix{T}}; S⋂YX::Vector{Matrix{T}}; CT_XAR::Vector{Matrix{T}}; CT_FE⋂PY::Vector{Matrix{T}}

  # IV/GMM only
  ZZ::Symmetric{T,Matrix{T}}; XY₂::Matrix{T}; XX::Symmetric{T,Matrix{T}}; H_2SLS::Symmetric{T,Matrix{T}}; V::Matrix{T}; ZY₂::Matrix{T}; X₂Y₂::Matrix{T}; X₁Y₂::Matrix{T}; ZR₁ZR₁::Symmetric{T,Matrix{T}}; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZZR₁::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; ZXinvXXXZ::Matrix{T}; H_2SLSmZZ::Symmetric{T,Matrix{T}}
  ZXinvXXXy₁par::Vector{T}; t₁Y::Vector{T}
  Y₂y₁::Vector{T}; twoR₁Zy₁::Vector{T}
  y₁y₁::T; y₁pary₁par::T
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}
  Y₂y₁par::Vector{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}
  kX::Int64
	Π̂::Matrix{T}

  StrEstimator{T}(isDGP, LIML, Fuller, κ) where T<:AbstractFloat = new(isDGP, LIML, Fuller, κ, Matrix{T}(undef,0,0))
end

mutable struct StrBootTest{T<:AbstractFloat}
  R::Matrix{T}; r::Vector{T}; R₁::Matrix{T}; r₁::Vector{T}
  y₁::Vector{T}; X₁::Matrix{T}; Y₂::Matrix{T}; X₂::Matrix{T}
  wt::Vector{T}; fweights::Bool
  LIML::Bool; Fuller::T; κ::T; ARubin::Bool
  B::Int64; auxtwtype::AuxWtType; rng::AbstractRNG; maxmatsize::Float16
  ptype::PType; null::Bool; bootstrapt::Bool
	ID::Matrix{Int64}; nbootclustvar::Int8; nerrclustvar::Int8; issorted::Bool; small::Bool
  FEID::Vector{Int64}; FEdfadj::Bool
  level::T; rtol::T
  madjtype::MAdjType; NH₀::Int16
  ML::Bool; β̈::Vector{T}; A::Symmetric{T,Matrix{T}}; sc::Matrix{T}
  willplot::Bool; gridmin::Vector{T}; gridmax::Vector{T}; gridpoints::Vector{Float32}

  q::Int16; twotailed::Bool; scorebs::Bool; robust::Bool

  WRE::Bool; initialized::Bool; NFE::Int64; FEboot::Bool; granular::Bool; NErrClustCombs::Int16; subcluster::Int8; BFeas::Int64; interpolating::Bool
  v_sd::T; notplotted::Bool
  confpeak::Vector{T}
  ID✻::Vector{Int64}; ID✻_✻⋂::Vector{Int64}
  anchor::Vector{T}; poles::Vector{T}; numer::Matrix{T}
  CI::Matrix{T}
  peak::NamedTuple{(:X, :p), Tuple{Vector{T}, T}}

	Nobs::Int64; NClustVar::Int8; kX₁::Int64; kX₂::Int64; kY₂::Int64; WREnonARubin::Bool; boottest!::Function
	coldotplus!::Function; colquadformminus!::Function; matmulplus!::Function; panelsum!::Function

  sqrt::Bool; _Nobs::T; kZ::Int64; sumwt::T; haswt::Bool; REst::Bool; multiplier::T; smallsample::T
		dof::Int64; dof_r::T; p::T; BootClust::Int8
		purerobust::Bool; N✻::Int64; Nw::Int64; enumerate::Bool; interpolable::Bool; interpolate_u::Bool; kX::Int64
  _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; u✻::Matrix{T}; CT_WE::Matrix{T}
  info✻::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; info⋂_✻⋂::Vector{UnitRange{Int64}}
  JN⋂N✻::Matrix{T}; statDenom::Matrix{T}; uXAR::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; β̈dev::Matrix{T}; δdenom_b::Matrix{T}
	_J⋂::Matrix{T}; YY✻_b::Matrix{T}; YPXY✻_b::Matrix{T}; numerw::Matrix{T}; Zyg::Vector{Matrix{T}}; numer_b::Vector{T}; dist::Matrix{T}
		
	distCDR::Matrix{T}; plotX::Vector{Vector{T}}; plotY::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  numersum::Vector{T}; ü₀::Vector{T}; invFEwt::Vector{T}
	β̈s::Matrix{T}; As::Matrix{T}
	info✻⋂::Vector{UnitRange{Int64}}; info⋂::Vector{UnitRange{Int64}}; ID✻⋂::Matrix{T}
	ü::Vector{T}
	DGP::StrEstimator{T}; Repl::StrEstimator{T}; M::StrEstimator{T}
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; S✻⋂u₁XinvXX::Matrix{Matrix{T}}; S✻UU::Matrix{Vector{T}}; CTUX::Matrix{Matrix{T}}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; IDCT⋂✻::Vector{Vector{Int64}}; infoCT⋂✻::Vector{Vector{UnitRange{Int64}}}; S✻XU::Vector{Matrix{T}}; invXXS✻XU::Vector{Matrix{T}}; invZperpZperpS✻ZperpU::Vector{Matrix{T}}; S✻uY::Vector{Matrix{T}}; S✻UMZperp::Vector{Matrix{T}}; S✻UPX::Vector{Matrix{T}}; S✻ZperpU::Vector{Matrix{T}}; CTFEU::Vector{Matrix{T}}
  ∂denom∂r::Vector{Matrix{Matrix{T}}}; ∂Jcd∂r::Vector{Matrix{Matrix{T}}}
  ∂²denom∂r²::Matrix{Matrix{Matrix{T}}}
	FEs::Vector{StrFE{T}}
  T1L::Vector{Matrix{T}}; T1R::Vector{Matrix{T}}
	crosstab⋂✻ind::Vector{Int64}; crosstabBootind::Vector{Int64}
  seed::UInt64

	S✻XY₂::Array{T,3}; S✻XX::Array{T,3}; S✻XZ::Array{T,3}; S✻Xy₁::Array{T,3}; S✻XZR₁::Array{T,3}
	invXXS✻XY₂::Array{T,3}; invXXS✻XX::Array{T,3}; invXXS✻XZ::Array{T,3}; invXXS✻Xy₁::Array{T,3}; invXXS✻XZR₁::Array{T,3}
	S✻⋂XY₂::Array{T,3}; S✻⋂XX::Array{T,3}; S✻⋂XZ::Array{T,3}; S✻⋂Xy₁::Array{T,3}; S✻⋂XZR₁::Array{T,3}
	invXXS✻⋂XY₂::Array{T,3}; invXXS✻⋂XX::Array{T,3}; invXXS✻⋂XZ::Array{T,3}; invXXS✻⋂Xy₁::Array{T,3}; invXXS✻⋂XZR₁::Array{T,3}
	S✻ZperpY₂::Array{T,3}; S✻ZperpX::Array{T,3}; S✻ZperpZ::Array{T,3}; S✻Zperpy₁::Array{T,3}; S✻ZperpZR₁::Array{T,3}
	invZperpZperpS✻ZperpY₂::Array{T,3}; invZperpZperpS✻ZperpX::Array{T,3}; invZperpZperpS✻ZperpZ::Array{T,3}; invZperpZperpS✻Zperpy₁::Array{T,3}; invZperpZperpS✻ZperpZR₁::Array{T,3}
	_ID✻⋂::Vector{Int}

  StrBootTest{T}(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, LIML, 
	               Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, scorebs, bootstrapt, ID, nbootclustvar, nerrclustvar, issorted, robust, small, FEID, FEdfadj, level, rtol, madjtype, NH₀, ML,
								 β̈, A, sc, willplot, gridmin, gridmax, gridpoints, turbo) where T<:Real =

		begin
			kX₂ = ncols(X₂)
			scorebs = scorebs || iszero(B) || ML
			WREnonARubin = !(iszero(kX₂) || scorebs) && !ARubin

			new(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, LIML || !iszero(Fuller), 
					Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, bootstrapt, ID, nbootclustvar, nerrclustvar, issorted, small, FEID, FEdfadj, level, rtol, madjtype, NH₀, ML, 
					β̈, A, sc, willplot, gridmin, gridmax, gridpoints,
				nrows(R), ptype==symmetric || ptype==equaltail, scorebs, robust || nerrclustvar>0,
				false, false, 0, false, false, 0, 0, 0, false,
				one(T), true,
				[zero(T)],
				Vector{Int64}(undef,0), Vector{Int64}(undef,0),
				Vector{T}(undef,0), Vector{T}(undef,0), Matrix{T}(undef,0,0),
				Matrix{T}(undef,0,0),
				(X = Vector{T}(undef,0), p = T(NaN)),
				nrows(X₁), ncols(ID), ncols(X₁), kX₂, ncols(Y₂), WREnonARubin, WREnonARubin ? boottestWRE! : boottestOLSARubin!, 
				turbo ? coldotplus_turbo! : coldotplus_nonturbo!, turbo ? colquadformminus_turbo! : colquadformminus_nonturbo!, turbo ? matmulplus_turbo! : matmulplus_nonturbo!, turbo ? panelsum_turbo! : panelsum_nonturbo!)
		end
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
	nothing
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
				$(esc(dest)) .= $(esc(content))
			else
				$(esc(dest))[$(esc(:o)).WeightGrp[$(esc(:w))]] .= $(esc(content))
			end
		end
  else
	  return quote
			local _content = $(esc(content))
	    if isone($(esc(:o)).Nw)
	  	  $(esc(dest)) = _content
	    else
	  	  $(esc(dest))[:,$(esc(:o)).WeightGrp[$(esc(:w))]] = _content
	    end
	  end
  end
end

macro clustAccum!(X, c, Y)  # efficiently add a cluster combination-specific term, factoring in the needed multiplier and sign
  return quote
		local _Y = $(esc(Y))
	  if isone($(esc(c)))
	    if isone($(esc(:o)).clust[1].multiplier)
	  	  $(esc(X)) = $(esc(:o)).clust[1].even ? _Y : -_Y
	    else
	  	  $(esc(X)) = _Y * ($(esc(:o)).clust[1].even ? $(esc(:o)).clust[1].multiplier : -$(esc(:o)).clust[1].multiplier)
	    end
	  elseif $(esc(:o)).clust[$(esc(c))].even
	    if isone($(esc(:o)).clust[$(esc(c))].multiplier)
	  	  $(esc(X)) .+= _Y
	    else
	  	  $(esc(X)) .+= _Y .* $(esc(:o)).clust[$(esc(c))].multiplier
	    end
	  elseif isone($(esc(:o)).clust[$(esc(c))].multiplier)
	    $(esc(X)) .-= _Y
	  else
	    $(esc(X)) .-= _Y .* $(esc(:o)).clust[$(esc(c))].multiplier
	  end
  end
end

function getdist(o::StrBootTest, diststat::DistStatType=nodist)
  if diststat == numer
	  _numer = isone(o.v_sd) ? o.numer : o.numer / o.v_sd
	  o.distCDR = (@view _numer[:,2:end]) .+ o.r
	  sort!(o.distCDR)
  elseif nrows(o.distCDR)==0  # return test stats
    if length(o.dist) > 1
	    o.distCDR = (@view o.dist[2:end])' * o.multiplier
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

# store p value in o.p. Return optionally-multiple-hypothesis-adjusted p value. Robust to missing bootstrapped values interpreted as +infinity.
function getp(o::StrBootTest{T}) where T
  o.boottest!(o)
  tmp = o.dist[1]
  isnan(tmp) && return tmp
  if o.B>0
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
		o.p = ccdf(o.small ? FDist{T}(T(o.dof), o.dof_r) : Chisq{T}(T(o.dof)), o.sqrt ? tmp^2 : tmp)
		if o.sqrt && !o.twotailed
			o.p /= 2
			(o.ptype==upper) == (tmp<0) && (o.p = 1 - o.p)
		end
  end
	
	if o.madjtype==bonferroni min(one(T), o.NH₀ * o.p)
  elseif o.madjtype==sidak  one(T) - (one(T) - o.p) ^ o.NH₀
  else o.p
  end
end

getb(o::StrBootTest) = isone(o.v_sd) ? o.numer[:,1] : o.numer[:,1] / o.v_sd  # numerator for full-sample test stat
getV(o::StrBootTest) = o.statDenom / ((isone(o.v_sd) ? o.smallsample : o.v_sd^2 * o.smallsample) * (o.sqrt ? o.multiplier^2 : o.multiplier) * o.dof)  # denominator for full-sample test stat
getv(o::StrBootTest) = @views isone(o.v_sd) ? o.v[:,2:end] : o.v[:,2:end] / o.v_sd  # wild weights
@inline getstat(o::StrBootTest) = o.multiplier * o.dist[1]
