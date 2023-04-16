# Definition of StrBootTest "class" for holding intermediate results, with associated utilities and get functions

struct StrClust{T<:Real}
	N::Int; multiplier::T; even::Bool
	order::Vector{Int64}
	info::Vector{UnitRange{Int64}}
end

# DesignerMatrix type to efficiently represent selection matrices
@enum MatType identity selection regular
struct DesignerMatrix{T} <: AbstractMatrix{T}
	type::MatType
	size::Tuple{Int64,Int64}
	p::Vector{Int64}  # rows with 1's/columns of left-multipled matrix to retain
	q::Vector{Int64}  # corresponding cols with the 1's/destination columns in result (remainder will be all 0)
	M::Matrix{T}
end

AdjointDesignerMatrix{T} = Adjoint{T, DesignerMatrix{T}} where T
# DesignerProduct{T}  = SubArray{T, 2, Matrix{T}, Tuple{Base.Slice{Base.OneTo{Int64}}, Vector{Int64}}, false} where T

function DesignerMatrix(X::AbstractMatrix{T}) where T
	if isa(X,AbstractMatrix) && length(X)>0 && all(sum(isone.(X) .| iszero.(X), dims=1) .== size(X,1)) && all(sum(isone.(X), dims=1) .<= one(T))
		d = diag(X)
		if length(d)==size(X,1)==size(X,2)==sum(isone.(d))
			DesignerMatrix{T}(identity, (length(d),length(d)), Int64[], Int64[], zeros(T,0,0))
		else
			p = getindex.(findall(>(0), sum(X;dims=2)),1)
			DesignerMatrix{T}(selection, size(X), p, [findall(>(0), X[i,:])[1] for i ∈ p], Matrix{T}(undef,0,0))
		end
	else
		DesignerMatrix{T}(regular, size(X), Int64[], Int64[], X)
	end
end

import Base.size, Base.getindex, Base.*
size(X::DesignerMatrix) = X.size
function getindex(X::DesignerMatrix{T}, i, j) where T
	if X.type==identity
		T(i==j)
	elseif X.type==selection
		a = findfirst(==(i),X.p)
		T(!isnothing(a) && j==X.q[a])
	else
		X.M[i,j]
	end
end

function *(X::AbstractMatrix{T}, Y::DesignerMatrix{T}) where T
	if Y.type==regular
		Base.:*(X, Y.M)
	elseif Y.type==identity
		X
	elseif length(Y.p)==Y.size[2]
		 X[:,Y.p[Y.q]]
	else
		dest = zeros(T, size(X,1), Y.size[2])  # Y is selection matrix with some cols all 0
		if length(dest)>0
			@inbounds for j ∈ eachindex(axes(Y.p,1))
				pⱼ, qⱼ = Y.p[j], Y.q[j] 
				@tturbo warn_check_args=false for i ∈ indices((dest,X),1)
					dest[i,qⱼ] = X[i,pⱼ]
				end
			end
		end
		dest
	end
end
function *(X::AbstractArray{T,3}, Y::DesignerMatrix{T}) where T
	if Y.type==regular
		view(Base.:*(X, Y.M),:,:,collect(1:size(Y.M,2)))  # collect() makes the column slicer Vector{Int64}, for type consistency
	elseif Y.type==identity
		view(X,:,:,collect(1:size(X,3)))
	elseif length(Y.q)==Y.size[2]
		 view(X,:,:,Y.p[Y.q])
	else
		dest = zeros(T, size(X,1), size(X,2), Y.size[2])  # Y is selection matrix with some cols all 0
		if length(dest)>0
			@inbounds for j ∈ axes(Y.p,1)
				pⱼ, qⱼ = Y.p[j], Y.q[j] 
				@tturbo warn_check_args=false for i ∈ indices((dest,X),1), g ∈ indices((dest,X),2)
					dest[i,g,qⱼ] = X[i,g,pⱼ]
				end
			end
		end
		view(dest,:,:,collect(1:size(dest,3)))
	end
end
function *(Y::AdjointDesignerMatrix{T}, X::AbstractMatrix{T}) where T
	if Y.parent.type==regular
		view(Base.:*(Y.parent.M', X),collect(1:size(Y.parent.M,2)),:)  # collect() makes the column slicer Vector{Int64}, for type consistency
	elseif Y.parent.type==identity
		view(X,collect(1:size(X,2)),:)
	elseif length(Y.parent.p)==Y.parent.size[2]
		 view(X,Y.parent.p[Y.parent.q],:)
	else
		dest = zeros(T, Y.parent.size[2], size(X,2))  # Y is selection matrix with some cols all 0
		if length(dest)>0
			@inbounds for j ∈ axes(Y.parent.p,1)
				pⱼ, qⱼ = Y.parent.p[j], Y.parent.q[j] 
				@tturbo warn_check_args=false for i ∈ indices((X,dest),2)
				 dest[qⱼ,i] = X[pⱼ,i]
				end
			end
		end
		view(dest,collect(1:size(dest,1)),:)
	end
end
function *(Y::AdjointDesignerMatrix{T}, X::AbstractArray{T,3}) where T
	if Y.parent.type==regular
		view(Base.:*(Y.parent.M', X),collect(1:size(Y.parent.M,2)),:,:)  # collect() makes the column slicer Vector{Int64}, for type consistency
	elseif Y.parent.type==identity
		view(X,collect(1:size(X,3)),:,:)
	elseif length(Y.parent.q)==Y.parent.size[2]
		 view(X,Y.parent.p[Y.parent.q],:,:)
	else
		dest = zeros(T, Y.parent.size[2], size(X,2), size(X,3))  # Y is selection matrix with some cols all 0
		if length(dest)>0
			@inbounds for j ∈ eachindex(axes(Y.parent.p,1))
				pⱼ, qⱼ = Y.parent.p[j], Y.parent.q[j] 
				@tturbo warn_check_args=false for i ∈ indices((dest,X),3), g ∈ indices((dest,X),2)
					dest[qⱼ,g,i] = X[pⱼ,g,i]
				end
			end
		end
		view(dest,collect(1:size(dest,1)),:,:)
	end
end

function *(X::AdjointDesignerMatrix{T}, Y::DesignerMatrix{T})::DesignerMatrix{T} where T 
	return DesignerMatrix(Matrix(X.parent)'Matrix(Y))
end

mutable struct StrEstimator{T<:AbstractFloat}
  const isDGP::Bool; const liml::Bool; const fuller::T; κ::T
  R₁perp::Matrix{T}; Rpar::Matrix{T}

  kZ::Int64
  y₁::Vector{T}; ü₁::Vector{Vector{T}}; u⃛₁::Vector{T}; β̈::Vector{T}; γ̈::Vector{T}; β̈₀::Vector{T}; invXXXy₁par::Vector{T}
  Yendog::Matrix{Bool}
  invZperpZperp::Matrix{T}; invZperpZperpZperpX::Matrix{T}; XZ::Matrix{T}; YPXY::Matrix{T}; R₁invR₁R₁::Matrix{T}; R₁invR₁R₁X::DesignerMatrix{T}; R₁invR₁R₁Y::DesignerMatrix{T}
	restricted::Bool; RperpX::DesignerMatrix{T}; RperpXperp::DesignerMatrix{T}; RRpar::Matrix{T}; RparX::Matrix{T}; RparY::DesignerMatrix{T}; RR₁invR₁R₁::Matrix{T}
	∂β̈∂r::Matrix{T}; YY::Matrix{T}; AR::Matrix{T}; XAR::Matrix{T}; Ü₂::Matrix{T}; Rt₁::Vector{T}
	invXX::Matrix{T}; Y₂::Matrix{T}; X₂::Matrix{T}; invH::Matrix{T}
	y₁par::Vector{T}; Xy₁par::Vector{T}
	A::Matrix{T}; Zpar::Matrix{T}; Zperp::Matrix{T}; X₁::Matrix{T}
	WXAR::Matrix{T}; CT_XAR::Vector{SparseMatrixCSC{T}}
	copyfromDGP::Bool
	S✻XX::Array{T,3}; XinvHjk::Vector{Matrix{T}}; invMjk::Vector{Matrix{T}}; invMjkv::Vector{T}; XXt1jk::Matrix{T}; t₁::Vector{T}

  # IV/GMM only
  ZZ::Matrix{T}; XY₂::Matrix{T}; XX::Matrix{T}; H_2SLS::Matrix{T}; V::Matrix{T}; ZY₂::Matrix{T}; ZR₁ZR₁::Matrix{T}; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZR₁Z::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; H_2SLSmZZ::Matrix{T}
  ZXinvXXXy₁par::Vector{T}
  Y₂y₁::Vector{T}; twoZR₁y₁::Vector{T}; y₁y₁::T; y₁pary₁par::T; Y₂Y₂::Matrix{T}; X₁X₁::Matrix{T}; X₂X₁::Matrix{T}; X₂X₂::Matrix{T}; X₂Y₂::Matrix{T}; X₁Z::Matrix{T}; X₂Z::Matrix{T}
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}; Y₂y₁par::Vector{T}; X₁Y₂::Matrix{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}
  kX::Int64; kX₁::Int64; kZperp::Int64
	Π̈ ::Matrix{T}; t₁Y::Vector{T}; PXZ::Matrix{T}
	S✻⋂ZperpZpar::Array{T,3}; S✻⋂ZperpY₂::Array{T,3}; S⋂y₁X₁::Array{T,3}; S⋂y₁X₂::Array{T,3}; S✻⋂ZperpZperp::Array{T,3}; S✻⋂Zperpy₁::Array{T,3}; S✻⋂XZR₁::Array{T,3}; S✻⋂XY₂::Array{T,3}; S✻⋂XZperp::Array{T,3}; S✻⋂XX::Array{T,3}; S✻⋂XZpar::Array{T,3}
	S✻⋂X₁Y₂::Array{T,3}; S✻⋂X₂Y₂::Array{T,3}; S✻ZparY₂::Array{T,3}; S✻y₁y₁::Array{T,3}; S✻Zpary₁::Array{T,3}; S✻ZR₁y₁::Array{T,3}; S✻ZR₁Y₂::Array{T,3}; S✻ZR₁ZR₁::Array{T,3}; S✻ZR₁Z::Array{T,3}
	S✻⋂ZperpZR₁::Array{T,3}
	S✻Y₂Y₂::Array{T,3}; S✻ZparZpar::Array{T,3}; S✻Y₂y₁::Array{T,3}; S✻⋂Xy₁::Array{T,3}
	ZparX::Matrix{T}; ZperpX₁::Matrix{T}; ZperpX₂::Matrix{T}; Zperpy₁::Vector{T}; ZperpY₂::Matrix{T}; ZperpZpar::Matrix{T}; ZperpZR₁::Matrix{T};
	invZperpZperpZperpX₁::Matrix{T}; invZperpZperpZperpX₂::Matrix{T}; invZperpZperpZperpy₁::Vector{T}; invZperpZperpZperpY₂::Matrix{T}; S✻UY₂::Matrix{T}; invZperpZperpZperpZ::Matrix{T}; invZperpZperpZperpZR₁::Matrix{T}
	Ü₂Ü₂::Matrix{T}; γ̈X::Vector{T}; γ̈Y::Vector{T}; γ⃛::Vector{T}; Xȳ₁::Vector{T}; ȳ₁ȳ₁::T; XÜ₂::Matrix{T}; ȳ₁Ü₂::Matrix{T}; Ȳ₂::Matrix{T}; ȳ₁::Vector{T}
	Xpar₁toZparX::DesignerMatrix{T}

	X₁ⱼₖ::Matrix{T}; X₂ⱼₖ::Matrix{T}; y₁ⱼₖ::Vector{T}; Y₂ⱼₖ::Matrix{T}; Zⱼₖ::Matrix{T}; ZR₁ⱼₖ::Matrix{T}; Y₂y₁ⱼₖ::Array{T,3}; X₂y₁ⱼₖ::Array{T,3}; X₁y₁ⱼₖ::Array{T,3}; Zy₁ⱼₖ::Array{T,3}; XZⱼₖ::Array{T,3}; ZZⱼₖ::Array{T,3}; ZY₂ⱼₖ::Array{T,3}; y₁y₁ⱼₖ::Array{T,3}; XY₂ⱼₖ::Array{T,3}; invXXⱼₖ::Array{T,3}; XXⱼₖ::Array{T,3}; X₁ZR₁ⱼₖ::Array{T,3}; X₂ZR₁ⱼₖ::Array{T,3}; ZR₁Zⱼₖ::Array{T,3}; twoZR₁y₁ⱼₖ::Array{T,3}; ZR₁ZR₁ⱼₖ::Array{T,3}; ZR₁Y₂ⱼₖ::Array{T,3} 
	Y₂y₁parⱼₖ::Array{T,3}; Zy₁parⱼₖ::Array{T,3}; y₁pary₁parⱼₖ::Array{T,3};	Xy₁parⱼₖ::Array{T,3}; y₁parⱼₖ::Vector{T}; H_2SLSⱼₖ::Array{T,3}; H_2SLSmZZⱼₖ::Array{T,3}; invHⱼₖ::Array{T,3}
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
	clustid::Matrix{Int64}; const NBootClustVar::Int8; const NErrClustVar::Int8; const issorted::Bool; const small::Bool; const clusteradj::Bool; const clustermin::Bool
  FEID::Vector{Int64}; FEdfadj::Int64
  const level::T; const rtol::T
  const madjtype::Symbol; const NH₀::Int16
  const ml::Bool; β̈::Vector{T}; A::Matrix{T}; sc::Matrix{T}
  const willplot::Bool; gridmin::Vector{T}; gridmax::Vector{T}; gridpoints::Vector{Float32}; overwrite::Bool

  const q::Int16; const twotailed::Bool; const jk::Bool; scorebs::Bool; const robust::Bool

  WRE::Bool; initialized::Bool; FEboot::Bool; granular::Bool; granularjk::Bool; NErrClustCombs::Int16; subcluster::Int8; BFeas::Int64; interpolating::Bool
  confpeak::Vector{T}
  ID✻::Vector{Int64}; ID⋂::Vector{Int64}; ID✻_✻⋂::Vector{Int64}
  anchor::Vector{T}; poles::Vector{T}; numer::Matrix{T}
  ci::Matrix{T}
  peak::NamedTuple{(:X, :p), Tuple{Vector{T}, T}}

	NFE::Int64; const Nobs::Int64; const NClustVar::Int8; const kX₁::Int64; const kX₂::Int64; const kY₂::Int64; const WREnonARubin::Bool; const boottest!::Function
	# end of fields initialized by initializer

  sqrt::Bool; _Nobs::T; kZ::Int64; sumwt::T; haswt::Bool; sqrtwt::Vector{T}; multiplier::T; smallsample::T; getci::Bool
		dof::Int64; dof_r::T; p::T; BootClust::Int8; ncolsv::Int64
		purerobust::Bool; N✻::Int64; N⋂::Int64; N✻⋂::Int64; Nw::Int64; enumerate::Bool; interpolable::Bool; interpolate_u::Bool; kX::Int64
  _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; u✻::Matrix{T}
  info✻::Vector{UnitRange{Int64}}; info✻_✻⋂::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; info⋂_✻⋂::Vector{UnitRange{Int64}}
  statDenom::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; β̈dev::Matrix{T}
	numerw::Matrix{T}; numer_b::Vector{T}; dist::Matrix{T}

	distCDR::Matrix{T}; plotX::Vector{Vector{T}}; plotY::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  numersum::Vector{T}; u✻₀::Matrix{T}; invsumFEwt::Vector{T}; FEwt::Vector{T}
	β̈s::Matrix{T}; As::Array{T,3}; β̈sAs::Matrix{T}
	info✻⋂::Vector{UnitRange{Int64}}; info⋂::Vector{UnitRange{Int64}}
	DGP::StrEstimator{T}; Repl::StrEstimator{T}; M::StrEstimator{T}
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; 
	S✻UU::Matrix{SubArray{T, 1, Array{T,3}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}; S✻u₁u₁::Array{T,3}; S✻U₂paru₁::Array{T,3}; S✻U₂parU₂par::Array{T,3}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; S✻Xu₁::Array{T,3}; S✻XU₂par::Array{T,3}; S✻XU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}; invXXS✻Xu₁::Array{T,3}; invXXS✻XU₂par::Array{T,3}; invXXS✻XU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	invZperpZperpS✻ZperpU₂par::Array{T,3}; invZperpZperpS✻Zperpu₁::Array{T,3}; invZperpZperpS✻ZperpU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}} 
	S✻Zperpu₁::Array{T,3}; S✻ZperpU₂par::Array{T,3}; S✻ZperpU::Vector{SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	Mjw::Vector{T}; CT✻FEU₂::Vector{SparseMatrixCSC{T}}; CT✻FEU::Vector{SparseMatrixCSC{T}}; invFEwtCT✻FEU::Vector{SparseMatrixCSC{T}}
	S✻ȲUfold::Array{T,3}
  ∂denom∂r::Array{Matrix{T},3}; ∂Jcd∂r::Array{Matrix{T},3}; ∂²denom∂r²::Array{Matrix{T},4}
	T1L::Matrix{T}; T1R::Matrix{T}; J⋂s::Array{T,3}; J⋂s1::Vector{Matrix{T}}; β̈v::Matrix{T}
	crosstab⋂✻ind::Vector{Int64}
  seed::UInt64

	S✻XY₂::Array{T,3}; S✻XX::Array{T,3}; S✻XDGPZ::Array{T,3}; S✻Xy₁::Array{T,3}; S✻XZR₁::Array{T,3}
	invXXS✻XDGPZ::Array{T,3}; invXXS✻Xy₁::Array{T,3}; invXXS✻XDGPZR₁::Array{T,3}
	S✻⋂XY₂::Array{T,3}; S✻⋂XX::Array{T,3}; S✻⋂XDGPZ::Array{T,3}; S✻⋂Xy₁::Array{T,3}; S✻⋂X_DGPZR₁::Array{T,3}
	S✻ZperpY₂::Array{T,3}; S✻ZperpX::Array{T,3}; S✻ZperpDGPZ::Array{T,3}; S✻Zperpy₁::Array{T,3}; S✻ZperpDGPZR₁::Array{T,3}
	invZperpZperpS✻ZperpY₂::Array{T,3}; invZperpZperpS✻ZperpX::Array{T,3}; invZperpZperpS✻ZperpDGPZ::Array{T,3}; invZperpZperpS✻Zperpy₁::Array{T,3}; invZperpZperpS✻ZperpDGPZR₁::Array{T,3}
	S✻Y₂y₁::Array{T,3}; S✻DGPZy₁::Array{T,3}; S✻y₁y₁::Array{T,3}; S✻DGPZR₁y₁::Array{T,3}
	negS✻UMZperpX::Vector{Array{T,3}}; S⋂XZperpinvZperpZperp::Array{T,3}; CT✻FEX::Vector{SparseMatrixCSC{T}}; CT⋂FEX::Vector{SparseMatrixCSC{T}}; CT✻FEY₂::Vector{SparseMatrixCSC{T}}; CT✻FEZ::Vector{SparseMatrixCSC{T}}; CT✻FEy₁::Vector{SparseMatrixCSC{T}}; CT✻FEZR₁::Vector{SparseMatrixCSC{T}}
	S✻Y₂Y₂::Array{T,3}; S✻ZparY₂Z_DGPZ::Array{T,3}; S✻DGPZY₂::Array{T,3}; S✻DGPZR₁Y₂::Array{T,3}; S✻DGPZDGPZ::Array{T,3};  S✻DGPZR₁DGPZR₁::Array{T,3}; S✻DGPZR₁DGPZ::Array{T,3}; S✻DGPZR₁X::Array{T,3}
	XinvXX::Matrix{T}; S⋂XX::Array{T,3}
	S✻⋂XU₂::Array{T,3}; S✻⋂XÜ₂par::Array{T,3}; S✻⋂Xu₁::Array{T,3}; S✻XU₂::Array{T,3}; S✻ZperpU₂::Array{T,3}; invZperpZperpS✻ZperpU₂::Array{T,3}; invXXS✻XU₂::Array{T,3}

	YY₁₁::Matrix{T}; YY₁₂::Matrix{T}; YY₂₂::Matrix{T}; YPXY₁₁::Matrix{T}; YPXY₁₂::Matrix{T}; YPXY₂₂::Matrix{T}
	YY₁₂YPXY₁₂::Matrix{T}; x₁₁::Matrix{T}; x₁₂::Matrix{T}; x₂₁::Matrix{T}; x₂₂::Matrix{T}; κs::Matrix{T}; numerWRE::Matrix{T}
	δnumer::Matrix{T}; YY✻::Array{T,3}; YPXY✻::Array{T,3}; κWRE::Array{T,3}; denomWRE::Array{T,3}; ARpars::Array{T,3}; J⋂ARpars::Vector{Array{T,3}}; Jc::Vector{Array{T,3}}
	invZperpZperpS✻ZperpUv::Matrix{T}; S✻ZperpUv::Matrix{T}; CT✻FEUv::Matrix{T}; S✻UMZperp::Matrix{T}; PXY✻::Matrix{T}; S✻UZperpinvZperpZperpv::Matrix{T}; invFEwtCT✻FEUv::Matrix{T}
	F₁::Matrix{T}; F₁β::Matrix{T}; F₂::Matrix{T}; willfill::Bool; not2SLS::Bool
	invXXXZ̄::Matrix{T}; XȲ::Matrix{T}; ZÜ₂par::Matrix{T}; ȲȲ::Matrix{T}
	S✻ȳ₁u₁::Array{T,3}; S✻Z̄u₁::Array{T,3}; S✻ȳ₁Ü₂par::Array{T,3}; S✻Z̄Ü₂par::Array{T,3};	PXZ̄::Matrix{T}; S✻XUv::Matrix{T}; S✻ȲU::Matrix{SubArray{T, 1, Array{T, 3}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}
	Π̈Rpar::Matrix{T}; Ü₂par::Matrix{T}; Z̄::Matrix{T}; S⋂ȳ₁X::Array{T,3}; S⋂ReplZ̄X::Array{T,3}

	StrBootTest{T}(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, liml, 
	               fuller, κ, arubin, B, auxtwtype, rng, maxmatsize, ptype, null, jk, scorebs, bootstrapt, clustid, NBootClustVar, NErrClustVar, issorted, robust, small, clusteradj, clustermin,
								 FEID, FEdfadj, level, rtol, madjtype, NH₀, ml,
								 β̈, A, sc, willplot, gridmin, gridmax, gridpoints, overwrite) where T<:Real =
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
					clustid, NBootClustVar, NErrClustVar, issorted, small, clusteradj, clustermin, 
					FEID, FEdfadj,
					level, rtol,
					madjtype, NH₀,
					ml, β̈, A, sc,
					willplot, gridmin, gridmax, gridpoints, overwrite,
					nrows(R), ptype == :symmetric || ptype == :equaltail, jk, scorebs, robust || NErrClustVar>0,
					false, false, false, false, false, 0, 0, 0, false,
					[zero(T)],
					Vector{Int64}(undef,0), Vector{Int64}(undef,0), Vector{Int64}(undef,0),
					Vector{T}(undef,0), Vector{T}(undef,0), Matrix{T}(undef,0,0),
					Matrix{T}(undef,0,0),
					(X = Vector{T}(undef,0), p = T(NaN)),
					0, nrows(y₁), ncols(clustid), ncols(X₁), kX₂, ncols(Y₂), WREnonARubin, WREnonARubin ? boottestWRE! : boottestOLSARubin!)
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
