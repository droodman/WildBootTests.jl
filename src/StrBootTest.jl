# Definition of StrBootTest "class" for holding intermediate results, with associated utilities and get functions

struct StrClust{T<:Real}
	N::Int; multiplier::T; even::Bool
	order::Vector{Int64}
	info::Vector{UnitRange{Int64}}
end

struct StrFE{T<:Real}
	is::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}
  wt::Vector{T}
  sqrtwt::Vector{T}
end

mutable struct StrEstimator{T<:AbstractFloat}
  const isDGP::Bool; const liml::Bool; const fuller::T; κ::T
  R₁perp::Matrix{T}; Rpar::Matrix{T}

  kZ::Int64
  y₁::Vector{T}; ü₁::Vector{Vector{T}}; u⃛₁::Vector{Vector{T}}; β̈::Vector{T}; γ̈::Vector{T}; β̈₀::Vector{T}; invXXXy₁par::Vector{T}
  Yendog::Matrix{Bool}
  invZperpZperp::#=Symmetric{T,=#Matrix{T}#=}=#; invZperpZperpZperpX::Matrix{T}; XZ::Matrix{T}; YPXY::#=Symmetric{T,=#Matrix{T}#=}=#; R₁invR₁R₁::Matrix{T}
	restricted::Bool; RperpX::Matrix{T}; RperpXperp::Matrix{T}; RRpar::Matrix{T}; RparX::Matrix{T}; RparY::Matrix{T}; RR₁invR₁R₁::Matrix{T}
	∂β̈∂r::Matrix{T}; YY::#=Symmetric{T,=#Matrix{T}#=}=#; AR::Matrix{T}; XAR::Matrix{T}; Ü₂::Vector{Matrix{T}}; Rt₁::Vector{T}
	invXX::#=Symmetric{T,=#Matrix{T}#=}=#; Y₂::Matrix{T}; X₂::Matrix{T}; invH::#=Symmetric{T,=#Matrix{T}#=}=#
	y₁par::Vector{T}; Xy₁par::Vector{T}
	A::#=Symmetric{T,=#Matrix{T}#=}=#; Z::Matrix{T}; Zperp::Matrix{T}; X₁::Matrix{T}
	WXAR::Matrix{T}; CT_XAR::Vector{Matrix{T}}

	S✻XX::Array{T,3}; XinvHjk::Vector{Matrix{T}}; invMjk::Vector{Matrix{T}}; invMjkv::Vector{T}; XXt1jk::Matrix{T}; t₁::Vector{T}

  # IV/GMM only
  ZZ::#=Symmetric{T,=#Matrix{T}#=}=#; XY₂::Matrix{T}; XX::#=Symmetric{T,=#Matrix{T}#=}=#; H_2SLS::#=Symmetric{T,=#Matrix{T}#=}=#; V::Matrix{T}; ZY₂::Matrix{T}; ZR₁ZR₁::#=Symmetric{T,=#Matrix{T}#=}=#; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZR₁Z::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; H_2SLSmZZ::#=Symmetric{T,=#Matrix{T}#=}=#
  ZXinvXXXy₁par::Vector{T}
  Y₂y₁::Vector{T}; twoZR₁y₁::Vector{T}; y₁y₁::T; y₁pary₁par::T; Y₂Y₂::Matrix{T}
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}; Y₂y₁par::Vector{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}
  kX::Int64; kX₁::Int64; kZperp::Int64
	Π̈ ::Matrix{T}
	S✻⋂ZperpZpar::Array{T,3}; S✻⋂ZperpY₂::Array{T,3}; S⋂y₁X₁::Array{T,3}; S⋂y₁X₂::Array{T,3}; S✻⋂ZperpZperp::Array{T,3}; ZperpX::Matrix{T}; S✻⋂Zperpy₁::Array{T,3}; S✻⋂XZR₁::Array{T,3}; S✻⋂XY₂::Array{T,3}; S✻⋂XZperp::Array{T,3}; S✻⋂XX::Array{T,3}; S✻⋂XZpar::Array{T,3}
	S✻⋂X₁Y₂::Array{T,3}; S✻⋂X₂Y₂::Array{T,3}; S✻ZparY₂::Array{T,3}; S✻y₁y₁::Array{T,3}; S✻Zpary₁::Array{T,3}; S✻ZR₁y₁::Array{T,3}; S✻ZR₁Y₂::Array{T,3}; S✻ZR₁ZR₁::Array{T,3}; S✻ZR₁Z::Array{T,3}
	ZperpZR₁::Matrix{T}; ZperpZperp::Matrix{T}; S✻⋂ZperpZR₁::Array{T,3}
	S✻Y₂Y₂::Array{T,3}; S✻ZparZpar::Array{T,3}; S✻Y₂y₁::Array{T,3}; S✻⋂Xy₁::Array{T,3}
	Xpar₁::Matrix{T}; X₁par::Matrix{T}
	invZperpZperpZperpX₁::Matrix{T}; invZperpZperpZperpX₂::Matrix{T}; invZperpZperpZperpy₁::Vector{T}; invZperpZperpZperpY₂::Matrix{T}; S✻UY₂::Matrix{T}; invZperpZperpZperpZpar::Matrix{T}; invZperpZperpZperpZR₁::Matrix{T}
	Ü₂Ü₂::Matrix{T}; γ̈X::Vector{T}; γ̈Y::Vector{T}; γ⃛::Vector{T}; Xȳ₁::Vector{T}; ȳ₁ȳ₁::T; XÜ₂::Matrix{T}; ȳ₁Ü₂::Matrix{T}; Ȳ₂::Matrix{T}; ȳ₁::Vector{T}

	X₁ⱼₖ::Matrix{T}; X₂ⱼₖ::Matrix{T}; y₁ⱼₖ::Vector{T}; Y₂ⱼₖ::Matrix{T}; Zⱼₖ::Matrix{T}; ZR₁ⱼₖ::Matrix{T}; Y₂y₁ⱼₖ::Array{T,3}; X₂y₁ⱼₖ::Array{T,3}; X₁y₁ⱼₖ::Array{T,3}; Zy₁ⱼₖ::Array{T,3}; XZⱼₖ::Array{T,3}; ZZⱼₖ::Array{T,3}; ZY₂ⱼₖ::Array{T,3}; y₁y₁ⱼₖ::Array{T,3}; XY₂ⱼₖ::Array{T,3}; invXXⱼₖ::Array{T,3}; X₁ZR₁ⱼₖ::Array{T,3}; X₂ZR₁ⱼₖ::Array{T,3}; ZZR₁ⱼₖ::Array{T,3}; twoZR₁y₁ⱼₖ::Array{T,3}; ZR₁ZR₁ⱼₖ::Array{T,3}; ZR₁Y₂ⱼₖ::Array{T,3} 
	Y₂y₁parⱼₖ::Array{T,3}; Zy₁parⱼₖ::Array{T,3}; y₁pary₁parⱼₖ::Array{T,3};	Xy₁parⱼₖ::Array{T,3}; y₁parⱼₖ::Vector{T}; XXⱼₖ::Array{T,3}; H_2SLSⱼₖ::Array{T,3}; H_2SLSmZZⱼₖ::Array{T,3}; invHⱼₖ::Array{T,3}
	β̈ⱼₖ::Array{T,3}; κⱼₖ::Array{T,3}; YPXYⱼₖ::Array{T,3}; YYⱼₖ::Array{T,3}; invXXXy₁parⱼₖ::Array{T,3}; ZXinvXXXy₁parⱼₖ::Array{T,3}

  StrEstimator{T}(isDGP, liml, fuller, κ) where T<:AbstractFloat = new(isDGP, liml, fuller, κ, Matrix{T}(undef,0,0))
end

mutable struct StrBootTest{T<:AbstractFloat}
  R::Matrix{T}; r::Vector{T}; R₁::Matrix{T}; r₁::Vector{T}
  y₁::Vector{T}; X₁::Matrix{T}; Y₂::Matrix{T}; X₂::Matrix{T}
  wt::Vector{T}; const fweights::Bool
  liml::Bool; const fuller::T; κ::T; const arubin::Bool
  B::Int64; const auxtwtype::Symbol; const rng::AbstractRNG; maxmatsize::Float16
  const ptype::Symbol; const null::Bool; const bootstrapt::Bool
	ID::Matrix{Int64}; const NBootClustVar::Int8; const NErrClustVar::Int8; const issorted::Bool; const small::Bool; const clusteradj::Bool; const clustermin::Bool
  FEID::Vector{Int64}; FEdfadj::Int64
  const level::T; const rtol::T
  const madjtype::Symbol; const NH₀::Int16
  const ml::Bool; β̈::Vector{T}; A::#=Symmetric{T,=#Matrix{T}#=}=#; sc::Matrix{T}
  const willplot::Bool; gridmin::Vector{T}; gridmax::Vector{T}; gridpoints::Vector{Float32}

  const q::Int16; const twotailed::Bool; const jk::Bool; scorebs::Bool; const robust::Bool

  WRE::Bool; initialized::Bool; NFE::Int64; FEboot::Bool; granular::Bool; granularjk::Bool; NErrClustCombs::Int16; subcluster::Int8; BFeas::Int64; interpolating::Bool
  confpeak::Vector{T}
  ID✻::Vector{Int64}; ID✻_✻⋂::Vector{Int64}
  anchor::Vector{T}; poles::Vector{T}; numer::Matrix{T}
  ci::Matrix{T}
  peak::NamedTuple{(:X, :p), Tuple{Vector{T}, T}}

	const Nobs::Int64; const NClustVar::Int8; const kX₁::Int64; const kX₂::Int64; const kY₂::Int64; const WREnonARubin::Bool; const boottest!::Function
	# end of fields initialized by initializer

  sqrt::Bool; _Nobs::T; kZ::Int64; sumwt::T; haswt::Bool; sqrtwt::Vector{T}; multiplier::T; smallsample::T; getci::Bool
		dof::Int64; dof_r::T; p::T; BootClust::Int8; ncolsv::Int64
		purerobust::Bool; N✻::Int64; N⋂::Int64; N✻⋂::Int64; Nw::Int64; enumerate::Bool; interpolable::Bool; interpolate_u::Bool; kX::Int64
  _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; u✻::Matrix{T}
  info✻::Vector{UnitRange{Int64}}; info✻_✻⋂::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; info⋂_✻⋂::Vector{UnitRange{Int64}}
  JN⋂N✻::Matrix{T}; statDenom::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; β̈dev::Matrix{T}
	numerw::Matrix{T}; numer_b::Vector{T}; dist::Matrix{T}

	distCDR::Matrix{T}; plotX::Vector{Vector{T}}; plotY::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  numersum::Vector{T}; u✻₀::Matrix{T}; invFEwt::Vector{T}
	β̈s::Matrix{T}; As::Array{T,3}
	info✻⋂::Vector{UnitRange{Int64}}; info⋂::Vector{UnitRange{Int64}}; ID✻⋂::Matrix{T}
	DGP::StrEstimator{T}; Repl::StrEstimator{T}; M::StrEstimator{T}
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; 
	S✻UU::Matrix{SubArray{T, 1, Array{T,3}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}; S✻u₁u₁::Array{T,3}; S✻U₂paru₁::Array{T,3}; S✻U₂parU₂par::Array{T,3}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; S✻Xu₁::Array{T,3}; S✻XU₂par::Array{T,3}; S✻XU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}; invXXS✻Xu₁::Array{T,3}; invXXS✻XU₂par::Array{T,3}; invXXS✻XU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	invZperpZperpS✻ZperpU₂par::Array{T,3}; invZperpZperpS✻Zperpu₁::Array{T,3}; invZperpZperpS✻ZperpU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}} 
	negS✻UMZperp::Array{T,3}; S✻UPX::Array{T,3}; S✻Zperpu₁::Array{T,3}; S✻ZperpU₂par::Array{T,3}; S✻ZperpU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	CT✻FEu₁::Array{T,3}; CT✻FEU₂par::Array{T,3}; CT✻FEU::Vector{SubArray{T, 2, Array{T,3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	invFEwtCT✻FEu₁::Array{T,3}; invFEwtCT✻FEU₂par::Array{T,3}; invFEwtCT✻FEU::Vector{SubArray{T, 2, Array{T,3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}; S✻ȲUfold::Array{T,3}
  ∂denom∂r::Array{Matrix{T},3}; ∂Jcd∂r::Array{Matrix{T},3}; ∂²denom∂r²::Array{Matrix{T},4}
	FEs::Vector{StrFE{T}}
  T1L::Matrix{T}; T1R::Matrix{T}; J⋂s::Array{T,3}; Q::Array{T,3 #=2=#}; β̈v::Matrix{T}
	crosstab⋂✻ind::Vector{Int64}; crosstab✻ind::Vector{Int64}
  seed::UInt64

	S✻XY₂::Array{T,3}; S✻XX::Array{T,3}; S✻XDGPZ::Array{T,3}; S✻Xy₁::Array{T,3}; S✻XZR₁::Array{T,3}
	invXXS✻XY₂::Array{T,3}; invXXS✻XX::Array{T,3}; invXXS✻XDGPZ::Array{T,3}; invXXS✻Xy₁::Array{T,3}; invXXS✻XDGPZR₁::Array{T,3}
	S✻⋂XY₂::Array{T,3}; S✻⋂XX::Array{T,3}; S✻⋂XDGPZ::Array{T,3}; S✻⋂Xy₁::Array{T,3}; S✻⋂X_DGPZR₁::Array{T,3}
	invXXS✻⋂XY₂::Array{T,3}; invXXS✻⋂XX::Array{T,3}; invXXS✻⋂XDGPZ::Array{T,3}; invXXS✻⋂Xy₁::Array{T,3}; invXXS✻⋂XDGPZR₁::Array{T,3}
	S✻ZperpY₂::Array{T,3}; S✻ZperpX::Array{T,3}; S✻ZperpDGPZ::Array{T,3}; S✻Zperpy₁::Array{T,3}; S✻ZperpDGPZR₁::Array{T,3}
	invZperpZperpS✻ZperpY₂::Array{T,3}; invZperpZperpS✻ZperpX::Array{T,3}; invZperpZperpS✻ZperpDGPZ::Array{T,3}; invZperpZperpS✻Zperpy₁::Array{T,3}; invZperpZperpS✻ZperpDGPZR₁::Array{T,3}
	_ID✻⋂::Vector{Int}
	S✻Y₂y₁::Array{T,3}; S✻DGPZy₁::Array{T,3}; S✻y₁y₁::Array{T,3}; S✻DGPZR₁y₁::Array{T,3}
	negS✻UMZperpX::Vector{Array{T,3}}; S⋂XZperpinvZperpZperp::Array{T,3}; CT✻FEX::Array{T,3}; CT⋂FEX::Array{T,3}; CT✻FEY₂::Array{T,3}; CT✻FEU₂::Array{T,3}; CT✻FEZ::Array{T,3}; CT✻FEy₁::Array{T,3}; CT✻FEZR₁::Array{T,3}
	S✻Y₂Y₂::Array{T,3}; S✻ZparY₂Z_DGPZ::Array{T,3}; S✻DGPZY₂::Array{T,3}; S✻DGPZR₁Y₂::Array{T,3}; S✻DGPZDGPZ::Array{T,3};  S✻DGPZR₁DGPZR₁::Array{T,3}; S✻DGPZR₁DGPZ::Array{T,3}; S✻DGPZR₁X::Array{T,3}
	XinvXX::Matrix{T}; FillingT₀::Matrix{Matrix{T}}
	S⋂ReplZX::Array{T,3}
	S✻⋂XU₂::Array{T,3}; S✻⋂XU₂par::Array{T,3}; S✻XU₂::Array{T,3}; S✻ZperpU₂::Array{T,3}; invZperpZperpS✻ZperpU₂::Array{T,3}; invXXS✻XU₂::Array{T,3}

	YY₁₁::Matrix{T}; YY₁₂::Matrix{T}; YY₂₂::Matrix{T}; YPXY₁₁::Matrix{T}; YPXY₁₂::Matrix{T}; YPXY₂₂::Matrix{T}
	YY₁₂YPXY₁₂::Matrix{T}; x₁₁::Matrix{T}; x₁₂::Matrix{T}; x₂₁::Matrix{T}; x₂₂::Matrix{T}; κs::Matrix{T}; numerWRE::Matrix{T}
	δnumer::Matrix{T}; YY✻::Array{T,3}; YPXY✻::Array{T,3}; κWRE::Array{T,3}; denomWRE::Array{T,3}; ARpars::Array{T,3}; J⋂ARpars::Array{T,3}; Jc::Vector{Array{T,3}}
	invZperpZperpS✻ZperpUv::Matrix{T}; S✻ZperpUv::Matrix{T}; CT✻FEUv::Matrix{T}; invFEwtCT✻FEUv::Matrix{T}; PXY✻::Matrix{T}; S✻UMZperpv::Matrix{T}
	T₀::Vector{T}; T₁::Matrix{T}; Qv::Matrix{T}; willfill::Bool; S✻diagUX::Array{T,3}
	invXXXZ̄::Matrix{T}; XȲ::Matrix{T}; ZÜ₂par::Matrix{T}; ȲȲ::Matrix{T}
	S✻ȳ₁u₁::Array{T,3}; S✻Z̄u₁::Array{T,3}; S✻ȳ₁U₂par::Array{T,3}; S✻Z̄U₂par::Array{T,3};	PXZ̄::Matrix{T}; S✻ȲU::Matrix{SubArray{T, 1, Array{T, 3}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	Π̈RparY::Matrix{T}; Ü₂par::Matrix{T}; Z̄::Matrix{T}; S⋂ȳ₁X::Array{T,3}; S⋂ReplZ̄X::Array{T,3}

	StrBootTest{T}(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, liml, 
	               fuller, κ, arubin, B, auxtwtype, rng, maxmatsize, ptype, null, jk, scorebs, bootstrapt, ID, NBootClustVar, NErrClustVar, issorted, robust, small, clusteradj, clustermin,
								 NFE, FEID, FEdfadj, level, rtol, madjtype, NH₀, ml,
								 β̈, A, sc, willplot, gridmin, gridmax, gridpoints) where T<:Real =
		begin
			kX₂ = ncols(X₂)
			scorebs = scorebs || iszero(B) || ml
			WREnonARubin = !(iszero(kX₂) || scorebs) && !arubin

			new(R, r, R₁, r₁, 
			    y₁, X₁, Y₂, X₂, 
					wt, fweights, 
					liml || !iszero(fuller), fuller, κ, arubin, 
					B, auxtwtype, rng, maxmatsize, 
					ptype, null, bootstrapt, 
					ID, NBootClustVar, NErrClustVar, issorted, small, clusteradj, clustermin, 
					FEID, FEdfadj,
					level, rtol,
					madjtype, NH₀,
					ml, β̈, A, sc,
					willplot, gridmin, gridmax, gridpoints,
					nrows(R), ptype == :symmetric || ptype == :equaltail, jk, scorebs, robust || NErrClustVar>0,
					false, false, NFE, false, false, false, 0, 0, 0, false,
					[zero(T)],
					Vector{Int64}(undef,0), Vector{Int64}(undef,0),
					Vector{T}(undef,0), Vector{T}(undef,0), Matrix{T}(undef,0,0),
					Matrix{T}(undef,0,0),
					(X = Vector{T}(undef,0), p = T(NaN)),
					nrows(X₁), ncols(ID), ncols(X₁), kX₂, ncols(Y₂), WREnonARubin, WREnonARubin ? boottestWRE! : boottestOLSARubin!)
		end
end

function getdist(o::StrBootTest, diststat::Symbol=:none)
  if diststat == :numer
	  _numer = o.numer
	  o.distCDR = (@view _numer[:,2:end]) .+ o.r
	  # sort!(o.distCDR, dims=1)
  elseif nrows(o.distCDR)==0  # return test stats
    if length(o.dist) > 1
	    o.distCDR = (@view o.dist[1,2:end])' * o.multiplier
	    # sort!(o.distCDR, dims=1)
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
  tmp = o.dist[1]
  isnan(tmp) && return tmp
  if o.B>0
  	if o.sqrt && o.ptype ≠ :upper
  	  if o.ptype == :symmetric
  	    n = sumlessabs(abs(tmp), o.dist)   # symmetric p value; do so as not to count missing entries in *dist
  	  elseif o.ptype == :equaltail
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
			(o.ptype == :upper) == (tmp<0) && (o.p = 1 - o.p)
		end
  end
	
	if o.madjtype == :bonferroni min(one(T), o.NH₀ * o.p)
  elseif o.madjtype == :sidak  one(T) - (one(T) - o.p) ^ o.NH₀
  else o.p
  end
end

getb(o::StrBootTest) = o.numer[:,1]  # numerator for full-sample test stat
getV(o::StrBootTest) = o.statDenom / (o.smallsample * (o.sqrt ? o.multiplier^2 : o.multiplier) * o.dof)  # denominator for full-sample test stat
getv(o::StrBootTest) = @views o.v[:,2:end]  # wild weights
@inline getstat(o::StrBootTest) = o.multiplier * o.dist[1]
