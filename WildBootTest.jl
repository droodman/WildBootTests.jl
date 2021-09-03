# WildBootTest.jl 0.1 3 September 2021

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


module Boottest
export StrBoottest, setsqrt!, setptype!, setstattype!, setX₁!, setX₂!, sety₁!, setY₂!, setobswt!, setsc!, setML!, setLIML!, setARubin!, setFuller!, setκ!,
       setquietly!, setβ!, setA!, setsmall!, setscoreBS!, setB!, setnull!, setWald!, setRao!, setID!, setFEID!, setlevel!, setptol!,
       setwillplot!, setrobust!, setR₁!, setR!, setgrid!, setmadjust!, setauxwttype!, setMaxMatSize!, setrng!, getdist, getp, getb, getV, getv, getrepsFeas,
       getNBootClust, getreps, getpadj, getstat, getdf, getdf_r, getplot, getpeak, getCI, AuxWeightType, PType

using LinearAlgebra, Random, Distributions, LoopVectorization, LazyArrays

@enum AuxWeightType rademacher mammen webb normal gamma
@enum PType symmetric equaltail lower upper
@enum MAdjType none bonferroni sidak
@enum StatType t c

struct StrClust{T<:Real}
	N::Int; multiplier::T; even::Bool
	order::Vector{Int64}
	info::Vector{UnitRange{Int64}}
end

struct StrFE{T<:Real}
	is::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}
  wt::Vector{T}
end

# right-multiply a data matrix by a matrix, with efficient handling of special cases like latter is identity
# pointer (real matrix) scalar XB(real matrix X, real matrix M) {
# 	scalar r
#   r = rows(M)
# 	return (all(colsum(M) .== 1) && all(colsum(!M) .== r-1)?  # is M 0's but for one 1 in each col, so it's just copying/reordering cols?
#              (all(diagonal(M)) && rows(M)==cols(M)?  # is M just the identity matrix?
#                 &X            :
#                 &X[:,colsum(M.*(1::r))]) :  # reorder X
#              &(X * M))
# }

# (wt of type UniformScaling <=> wt=I)
# @inline isone(X::AbstractArray) = X==ones(eltype(X),ones(Int64, ndims(X))...)
@inline sqrtNaN(x) = x<0 ? typeof(x)(NaN) : sqrt(x)
@inline invsym(X) = iszero(length(X)) ? X : inv(Symmetric(X))
@inline eigensym(X) = eigen(Symmetric(X))  # does eigen recognize symmetric matrices?
@inline symcross(X::AbstractArray, wt::AbstractVector) = Symmetric(X'*(wt.*X))  # maybe bad name since it means cross product in Julia
@inline symcross(X::AbstractArray, wt::UniformScaling) = Symmetric(X'X)  # maybe bad name since it means cross product in Julia
# @inline cross(X::AbstractArray, H::AbstractArray, Y::AbstractArray) = Symmetric(X'H*X)  unused?
@inline cross(X::AbstractArray, wt::AbstractVector, Y::AbstractArray) = X'*(wt.*Y)
@inline cross(X::AbstractArray, wt::UniformScaling, Y::AbstractArray) = X'Y
@inline vHadw(v::AbstractArray, w::AbstractVector) = v .* w
@inline vHadw(v::AbstractArray, w::UniformScaling) = v
@inline rows(X::AbstractArray) = size(X,1)
@inline cols(X::AbstractArray) = size(X,2)
@inline colsum(X::AbstractArray) = iszero(cols(X)) ? Matrix{eltype(X)}(undef,1,0) : isone(cols(X)) ? hcat(sum(X)) : sum(X, dims=1)
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))
@inline wtsum(wt::AbstractArray, X::AbstractArray) = wt'X
@inline wtsum(wt::UniformScaling, X::AbstractArray) = sum(X,dims=1)
checkI!(X::AbstractArray) = all(abs.(X - I) .< 10*eps(eltype(X))) ? I : X
@inline X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractMatrix) = @views X₁*B[1:size(X₁,2),:] + X₂*B[size(X₁,2)+1:end,:]
@inline X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractVector) = @views X₁*B[1:size(X₁,2)  ] + X₂*B[size(X₁,2)+1:end  ]

function coldot(A::AbstractMatrix, B::AbstractMatrix)  # colsum(A .* B)
  retval = zeros(promote_type(eltype(A), eltype(B)), 1, size(A,2))
  @turbo for i ∈ axes(A,2), j ∈ axes(A,1)
    retval[1,i] += A[j,i] * B[j,i]
  end
  retval
end
function coldot(A::AbstractMatrix)  # colsum(A .* A)
  retval = zeros(eltype(A), 1, size(A,2))
  @turbo for i ∈ axes(A,2), j ∈ axes(A,1)
    retval[1,i] += A[j,i] ^ 2
  end
  retval
end
function coldot(A::AbstractVector, B::AbstractMatrix)  # colsum(A .* B)
  retval = zeros(promote_type(eltype(A), eltype(B)), 1, size(B,2))
  @turbo for i ∈ axes(B,2), j ∈ axes(B,1)
    retval[1,i] += A[j] * B[j,i]
  end
  retval
end
coldot(A::AbstractMatrix, B::AbstractVector) = coldot(B,A)
coldot(A::AbstractVector, B::AbstractVector) = vec(dot(A,B))

# compute the norm of each col of A using quadratic form Q
function colquadform(Q::AbstractMatrix, A::AbstractMatrix) :: AbstractVector
  retval = zeros(promote_type(eltype(Q), eltype(A)), size(A,2))
  @turbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
    retval[i] += A[j,i] * Q[k,j] * A[k,i]
  end
  retval
end

 # From given row of given matrix, substract inner products of corresponding cols of A & B with quadratic form Q
function colquadformminus!(X::AbstractMatrix, row::Integer, Q::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
  	@turbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
  	  X[row,i] -= A[j,i] * Q[k,j] * B[k,i]
     end
 	X
end
colquadformminus!(X::AbstractMatrix, Q::AbstractMatrix, A::AbstractMatrix) = colquadformminus!(X, 1, Q, A, A)

function matmulplus!(A::Matrix, B::Matrix, C::Matrix)  # add B*C to A in place
  @turbo for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2)), j ∈ eachindex(axes(B,2),axes(C,1))
    A[i,k] += B[i,j] * C[j,k]
  end
end

abstract type Estimator end
struct OLS<:Estimator end
struct ARubin<:Estimator end
struct IVGMM<:Estimator end

mutable struct StrEstimator{T<:Real, E<:Estimator}
	parent
	isDGP::Bool; LIML::Bool; Fuller::T; κ::T
  R₁perp::Matrix{T}; Rpar::Union{UniformScaling{Bool}, Matrix{T}}

  kZ::Int64
	y₁::Vector{T}; ü₁::Vector{T}; u⃛₁::Vector{T}; β::Vector{T}; β₀::Vector{T}; PXy₁::Vector{T}; invXXXy₁par::Vector{T}
	Yendog::Vector{Bool}
  invZperpZperp::Matrix{T}; ZperpinvZperpZperp::Matrix{T}; XZ::Matrix{T}; PXZ::Matrix{T}; YPXY::Matrix{T}; R₁invR₁R₁::Matrix{T}; RperpX::Matrix{T}; RRpar::Matrix{T}; RparY::Matrix{T}; RR₁invR₁R₁::Matrix{T}; ∂β∂r::Matrix{T}; YY::Matrix{T}; AR::Matrix{T}; XAR::Matrix{T}; R₁invR₁R₁Y::Matrix{T}; invXXXZ::Matrix{T}; Ü₂::Matrix{T}; XinvXX::Matrix{T}; Rt₁::Vector{T}; invXX::Matrix{T}; Y₂::Matrix{T}; X₂::Matrix{T}; invH
	y₁par::Vector{T}; Xy₁par::Vector{T}
	A::Matrix{T}; Z::Matrix{T}; Zperp::Matrix{T}; X₁::Matrix{T}
	FillingT₀::Matrix{Vector{T}}
	WXAR::Matrix{T}; ScapPXYZperp::Vector{Matrix{T}}; ScapYX::Vector{Matrix{T}}; CT_XAR::Vector{Matrix{T}}; CT_FEcapPY::Vector{Matrix{T}}

	# IV/GMM only
  ZZ::Matrix{T}; XY₂::Matrix{T}; XX::Matrix{T}; H_2SLS::Matrix{T}; V::Matrix{T}; ZY₂::Matrix{T}; X₂Y₂::Matrix{T}; X₁Y₂::Matrix{T}; ZR₁ZR₁::Matrix{T}; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZZR₁::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; ZXinvXXXZ::Matrix{T}; H_2SLSmZZ::Matrix{T}
  ZXinvXXXy₁par::Vector{T}; t₁Y::Vector{T}
  Y₂y₁::Vector{T}; twoy₁ZR₁::Matrix{T}
  y₁y₁::T; y₁pary₁par::T
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}
	Y₂y₁par::Vector{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}

  StrEstimator{T,E}(parent) where T<:Real where E<:Estimator = new(parent, true, E==IVGMM, false, T(E==IVGMM ? 1 : 0), Matrix{T}(undef,0,0), I)
end

mutable struct StrBoottest{T<:Real}
	ARubin::Bool; LIML::Bool; Fuller::T; WRE::Bool; small::Bool; scoreBS::Bool; auxtwtype::AuxWeightType; ML::Bool; initialized::Bool; quietly::Bool; sqrt::Bool; ptype::PType; robust::Bool; NFE::Int64; FEboot::Bool; granular::Bool; NErrClustCombs::Int16; subcluster::Int8; B::Int64; BFeas::Int64; interpolating::Bool
  twotailed::Bool; null::Bool; dirty::Bool; willplot::Bool; u_sd::T; bootstrapt::Bool; notplotted::Bool; FEdfadj::Bool
  level::T
  rtol::T
  confpeak::Vector{T}; MaxMatSize::Float16
  Y₂::Matrix{T}; X₂::Matrix{T}; X₁::Matrix{T}; y₁::Vector{T}; Sc::Matrix{T}; ID::Matrix{T}; R₁::Matrix{T}; R::Matrix{T}; wt::Union{Vector{T}, UniformScaling}
  r₁::Vector{T}; r::Vector{T}
  IDBootData::Vector{Int64}; IDBootAll::Vector{Int64}
  κ::T
  anchor::Vector{T}; poles::Vector{T}; numer::Matrix{T}; Dist::Vector{T}
  CI::Matrix{T}
  madjtype::MAdjType
  rng::AbstractRNG

  Nobs::Int64; _Nobs::Int64; kZ::Int64; kY₂::Int64; kX₁::Int64; sumwt::T; NClustVar::Int8; haswt::Bool; REst::Bool; multiplier::T; smallsample::T
		WREnonARubin::Bool; df::Int64; df_r::Int64; NumH0s::Int32; p::T; NBootClustVar::Int8; NErrClust::Int64; BootClust::Int8
		purerobust::Bool; Nstar::Int64; Nw::Int64; enumerate::Bool; q::Int16; interpolable::Bool; interpolate_u::Bool; kX₂::Int64; kX::Int64
  FEID::Vector{T}; _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; ustar::Matrix{T}; CT_WE::Matrix{T}
  infoBootData::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; infoErrAll::Vector{UnitRange{Int64}}
  JNcapNstar::Matrix{T}; statDenom::Matrix{T}; uXAR::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; βdev::Matrix{T}; δdenom_b::Matrix{T}; _Jcap::Matrix{T}; YYstar_b::Matrix{T}; YPXYstar_b::Matrix{T}; numerw::Matrix{T}
	DistCDR::Matrix{T}; plotX::Matrix{T}; plotY::Vector{T}; β::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  gridmin::Vector{Union{T,Missing}}; gridmax::Vector{Union{T,Missing}}; gridpoints::Vector{Union{Int32,Missing}}; numersum::Vector{T}; ü₀::Vector{T}; invFEwt::Vector{T}
	peak::Matrix{T}; βs::Matrix{T}; As::Matrix{T}
	fweights::Bool
	infoAllData::Vector{UnitRange{Int64}}; infoCapData::Vector{UnitRange{Int64}}; IDAll::Matrix{T}; Ü₂par::Matrix{T}
	A::Matrix{T}; ü::Vector{T}
	DGP::StrEstimator{T,E} where E; Repl::StrEstimator{T,E} where E; M::StrEstimator{T,E} where E
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; SCTcapuXinvXX::Matrix{Matrix{T}}; SstarUU::Matrix{Vector{T}}; CTUX::Matrix{Matrix{T}}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; IDCTCapstar::Vector{Vector{Int64}}; infoCTCapstar::Vector{Vector{UnitRange{Int64}}}; SstarUX::Vector{Matrix{T}}; SstarUXinvXX::Vector{Matrix{T}}; SstarUZperpinvZperpZperp::Vector{Matrix{T}}; δdenom::Vector{Matrix{T}}; SstaruY::Vector{Matrix{T}}; SstarUMZperp::Vector{Matrix{T}}; SstarUPX::Vector{Matrix{T}}; SstarUZperp::Vector{Matrix{T}}; CTFEU::Vector{Matrix{T}}
  ∂denom∂r::Vector{Matrix{Matrix{T}}}; ∂Jcd∂r::Vector{Matrix{Matrix{T}}}
  ∂²denom∂r²::Matrix{Matrix{Matrix{T}}}
	FEs::Vector{StrFE}

  StrBoottest{T}() where T<:Real = new(false, false, 0, false, false, false, rademacher, false, false, false, false, symmetric, false, 0, false, false, 0, false, 0, 0, false,
                                       true, true, true, true, 1., true, true, true,
                                       .95,
                                       1e-6,
                                       [0], 0,
                                       Matrix{T}(undef,0,0), Matrix{T}(undef,0,0), Matrix{T}(undef,0,0), Vector{T}(undef,0), Matrix{T}(undef,0,0), Matrix{T}(undef,0,0), Matrix{T}(undef,0,0), Matrix{T}(undef,0,0), I,
                                       Vector{T}(undef,0), Vector{T}(undef,0),
                                       Vector{T}(undef,0), Vector{T}(undef,0),
                                       NaN,
                                       Vector{T}(undef,0), Vector{T}(undef,0), Matrix{T}(undef,0,0), Vector{T}(undef,0),
                                       Matrix{T}(undef,0,0),
                                       none,
                                       MersenneTwister(0))
end

function perp(A::AbstractMatrix)
	F = eigensym(A*invsym(A'A)*A')
	F.vectors[:, abs.(F.values) .< 1000*eps(eltype(A))]
  # checkI!(retval)
end

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 rows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T,E}, R₁::AbstractMatrix, R::Union{UniformScaling{Bool},AbstractMatrix}=Matrix{T}(undef,0,0)) where {T,E}
	if length(R₁) > 0
		o.R₁invR₁R₁ = invsym(R₁ * R₁')
		all(iszero.(diag(o.R₁invR₁R₁))) && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
    o.R₁invR₁R₁ = R₁'o.R₁invR₁R₁
		F = eigensym(o.R₁invR₁R₁ * R₁)
		o.R₁perp = F.vectors[:, abs.(F.values) .< 1000*eps(T)]  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
	else
		o.R₁invR₁R₁ = Matrix{T}(undef, o.parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
		RR₁perp = Matrix([R ; zeros(T, o.parent.kY₂, o.parent.kX₁) I])  # rows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

		length(R₁)>0 && (RR₁perp *= o.R₁perp)

		F = eigensym(RR₁perp'pinv(RR₁perp * RR₁perp')*RR₁perp); val = abs.(F.values) .> 1000*eps(T)
		o.Rpar   = F.vectors[:,   val]  # perp and par of RR₁perp
  	o.RperpX = F.vectors[:, .!val]

		if length(R₁) > 0  # fold model constraint factors into Rpar, RperpX
			o.Rpar   = o.R₁perp * o.Rpar
			o.RperpX = o.R₁perp * o.RperpX
		end

		o.RRpar = R * o.Rpar
		o.RperpX = o.RperpX[1:o.parent.kX₁,:]  # Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X₁, since all endogenous ones will be retained
		o.RparY = o.Rpar[o.parent.kX₁+1:end,:]  # part of Rpar that refers to Y₂
		o.R₁invR₁R₁Y = o.R₁invR₁R₁[o.parent.kX₁+1:end,:]
    o.RR₁invR₁R₁ = R * o.R₁invR₁R₁
	end
end

# stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
function InitVars!(o::StrEstimator{T,OLS} where T, Rperp::AbstractMatrix) # Rperp is for replication regression--no null imposed
  o.y₁par = o.parent.y₁
  H = symcross(o.parent.X₁, o.parent.wt)
  o.invH = inv(H)

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.invH : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
  o.β₀  = R₁AR₁ * cross(o.parent.X₁, o.parent.wt, o.y₁par)
  o.∂β∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

	o.A = iszero(length(Rperp)) ? o.invH : Rperp * invsym(Rperp'H*Rperp) * Rperp'   # for replication regression
  o.AR = o.A * o.parent.R'
	(o.parent.scoreBS || o.parent.robust) && (o.XAR = o.parent.X₁ * o.AR)
end

function InitVars!(o::StrEstimator{T,ARubin} where T, Rperp::AbstractMatrix = Matrix{Float64}(undef,0,0))
  X₂X₁ = cross(o.parent.X₂, o.parent.wt, o.parent.X₁)
  H = Symmetric([symcross(o.parent.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.parent.X₂, o.parent.wt)])  # XXX use LazyArrays?
  o.A = inv(H)
  o.AR = o.A * o.parent.R'
	(o.parent.scoreBS || o.parent.robust) && (o.XAR = X₁₂B(o.parent.X₁, o.parent.X₂, o.AR))

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.A : o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp'
  o.β₀   = R₁AR₁ * [cross(o.parent.X₁, o.parent.wt, o.parent.y₁) ; cross(o.parent.X₂, o.parent.wt, o.parent.y₁)]  # XXX use LazyArrays?
  o.∂β∂r = R₁AR₁ * [cross(o.parent.X₁, o.parent.wt, o.parent.Y₂) ; cross(o.parent.X₂, o.parent.wt, o.parent.Y₂)]
end

function InitVars!(o::StrEstimator{T,IVGMM}, Rperp::AbstractMatrix...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

  o.Zperp = o.parent.X₁ * o.RperpX  # XXX XB(o.parent.X₁, o.RperpX)
  o.invZperpZperp = iszero(length(o.Zperp)) ? Matrix{T}(undef,0,0) : inv(symcross(o.Zperp, o.parent.wt))
  o.ZperpinvZperpZperp = o.Zperp * o.invZperpZperp

  o.X₁ = o.parent.X₁ * perp(o.RperpX); o.X₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.X₁)  # FWL-process X₁ XXX XB(o.parent.X₁, perp(RperpX))
  o.X₂ = o.parent.X₂ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.X₂)                # FWL-process X₂
  X₂X₁ = cross(o.X₂, o.parent.wt, o.X₁)
  o.XX = Symmetric([symcross(o.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.X₂, o.parent.wt)])
  o.invXX = invsym(o.XX)

  o.Z   = X₁₂B(o.parent.X₁, o.parent.Y₂, o.Rpar     )  # Zpar
  o.ZR₁ = X₁₂B(o.parent.X₁, o.parent.Y₂, o.R₁invR₁R₁)

  o.Z   .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.Z  )  # partialling out
  o.ZR₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.ZR₁)
  o.Y₂ = o.parent.Y₂ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.Y₂)
  o.y₁ = o.parent.y₁ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.y₁)

  o.X₁Y₂ = cross(o.X₁, o.parent.wt, o.Y₂)
  o.X₂Y₂ = cross(o.X₂, o.parent.wt, o.Y₂)
  o.XY₂ = [o.X₁Y₂ ; o.X₂Y₂]
  o.Y₂y₁ = cross(o.Y₂ , o.parent.wt, o.y₁)
  o.X₂y₁ = cross(o.X₂ , o.parent.wt, o.y₁)
  o.X₁y₁ = cross(o.X₁ , o.parent.wt, o.y₁)
  o.y₁y₁ = cross(o.y₁ , o.parent.wt, o.y₁)
  o.Zy₁  = cross(o.Z  , o.parent.wt, o.y₁)
  o.XZ   = [cross(o.X₁, o.parent.wt, o.Z) ;
            cross(o.X₂, o.parent.wt, o.Z)  ]
  o.ZY₂ =  cross(o.Z  , o.parent.wt, o.Y₂)
  o.ZZ  =  symcross(o.Z  , o.parent.wt      )

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ

  if length(o.R₁invR₁R₁)>0
    o.X₂ZR₁    = cross(o.X₂ , o.parent.wt, o.ZR₁)
    o.X₁ZR₁    = cross(o.X₁ , o.parent.wt, o.ZR₁)
    o.ZZR₁     = cross(o.Z  , o.parent.wt, o.ZR₁)
    o.twoy₁ZR₁ = cross(o.y₁ , o.parent.wt, o.ZR₁) * 2
    o.ZR₁ZR₁   = symcross(o.ZR₁, o.parent.wt)
    o.ZR₁Y₂    = cross(o.ZR₁, o.parent.wt, o.Y₂ )
  else
    o.Y₂y₁par    = o.Y₂y₁
    o.X₂y₁par    = o.X₂y₁
    o.X₁y₁par    = o.X₁y₁
    o.Zy₁par     = o.Zy₁
    o.y₁pary₁par = o.y₁y₁
    o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
    o.y₁par      = o.y₁
  end

  o.V =  o.invXX * o.XZ # in 2SLS case, StrEstimator is (V' XZ)^-1 * (V'Xy₁). Also used in kZ-class and LIML robust VCV by Stata convention
  o.H_2SLS = Symmetric(o.V'o.XZ)  # Hessian
  (o.LIML || o.κ ≠ 1) && (o.H_2SLSmZZ = o.H_2SLS - o.ZZ)

  if o.isDGP
    !o.LIML && MakeH!(o, !isempty(Rperp)) # DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
  else
    o.kZ = cols(o.Rpar)
    o.ScapYX       = Vector{Matrix{T}}(undef, o.kZ+1)
    o.ScapPXYZperp = Vector{Matrix{T}}(undef, o.kZ+1)
    o.XinvXX = X₁₂B(o.X₁, o.X₂, o.invXX)

    o.Yendog = [true; vec(colsum(o.RparY.≠0)).>0]  # columns of Y = [y₁par Zpar] that are endogenous (normally all)

    if o.parent.robust  # for WRE replication regression, prepare for CRVE
      if o.parent.bootstrapt
        o.PXZ = X₁₂B(o.X₁, o.X₂, o.invXXXZ)

        o.FillingT₀ = Matrix{Matrix{T}}(undef, o.kZ+1, o.kZ+1)  # fixed component of groupwise term in sandwich filling
        o.parent.NFE>0 &&
          (o.CT_FEcapPY = Vector{Matrix{T}}(undef, o.kZ+1))
        for i ∈ 1:o.kZ
          uwt = vHadw(o.PXZ[:,i], o.parent.wt)
          o.parent.NFE>0 &&
            (o.CT_FEcapPY[i+1] = crosstabFE(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)
          tmp = panelsum(o.Z, uwt, o.parent.infoCapData)
          for j ∈ 1:o.kZ
            o.FillingT₀[i+1,j+1] = tmp[:,j]  # XXX make this a view or make FillingT₀ 4D array or vector{3D array}
          end
        end

        for i ∈ 1:o.kZ  # precompute various clusterwise sums
          o.ScapPXYZperp[i+1] = panelsum(o.Zperp, vHadw(o.PXZ[:,i], o.parent.wt), o.parent.infoCapData)  # Scap(P_(MZperpX) * Z .* Zperp)
          !o.parent.granular &&
            (o.ScapYX[i+1] = panelsum(o.X₁, o.X₂, vHadw(o.Z[:,i], o.parent.wt), o.parent.infoCapData))  # Scap(M_Zperp[Z or y₁] .* P_(MZperpX)])
        end
      end
    end
  end
end


# do most of estimation; for LIML r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β₀ (β when r=0) and ∂β∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function Estimate!(o::StrEstimator{T,OLS} where T, r₁::AbstractVector)
  o.β = o.β₀ - o.∂β∂r * r₁
end

function Estimate!(o::StrEstimator{T,ARubin} where T, r₁::AbstractVector)
  o.β = o.β₀ - o.∂β∂r * r₁
  o.y₁par = o.parent.y₁ - o.parent.Y₂ * r₁
end

function MakeH!(o::StrEstimator{T,IVGMM} where T, makeXAR::Bool=false)
  H = isone(o.κ) ? o.H_2SLS : o.ZZ + o.κ * o.H_2SLSmZZ
  o.invH = invsym(H)
  if makeXAR  # for replication regression in score bootstrap of IV/GMM
    o.A = length(o.Rperp)>0 ? o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp' : o.invH
    o.AR = o.A * (o.Rpar'o.parent.R')
    o.XAR = X₁₂B(o.X₁, o.X₂, o.V * o.AR)
  end
end

function Estimate!(o::StrEstimator{T,IVGMM} where T, r₁::AbstractVector)
	if length(o.R₁invR₁R₁)>0
		o.y₁pary₁par = o.y₁y₁ - (o.twoy₁ZR₁ * r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
		o.y₁par   = o.y₁ - o.ZR₁ * r₁
		o.Y₂y₁par = o.Y₂y₁  - o.ZR₁Y₂'r₁
		o.X₂y₁par = o.X₂y₁ - o.X₂ZR₁ * r₁
		o.X₁y₁par = o.X₁y₁ - o.X₁ZR₁ * r₁
		o.Zy₁par  = o.Zy₁ -  o.ZZR₁ * r₁
		o.Xy₁par  = [o.X₁y₁par ; o.X₂y₁par]
	end

  o.invXXXy₁par = o.invXX * o.Xy₁par
  o.ZXinvXXXy₁par = o.XZ'o.invXXXy₁par
  o.YY = [o.y₁pary₁par o.Zy₁par' ; o.Zy₁par o.ZZ]
  o.YPXY = [o.invXXXy₁par'o.Xy₁par  o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par o.ZXinvXXXZ]

  if o.isDGP
    if o.LIML
      o.κ = 1/(1 - eigvals(invsym(o.YY) * o.YPXY)[1])  # like Fast & Wild (81), but more stable, at least in Mata
      !iszero(o.Fuller) && (o.κ -= o.Fuller / (o.parent._Nobs - o.parent.kX))
      MakeH!(o)
    end

    o.β = o.invH * (isone(o.κ) ? o.ZXinvXXXy₁par : o.κ * (o.ZXinvXXXy₁par - o.Zy₁par) + o.Zy₁par)
    o.t₁Y = o.R₁invR₁R₁Y * r₁
  elseif o.parent.WREnonARubin  # if not score bootstrap of IV/GMM
    o.Rt₁ = o.RR₁invR₁R₁ * r₁

    if o.parent.robust && o.parent.bootstrapt  # prepare WRE replication regressions
      o.PXy₁ = X₁₂B(o.X₁, o.X₂, o.invXXXy₁par)

      uwt = vHadw(o.PXy₁, o.parent.wt)
      tmp = panelsum(o.Z, uwt, o.parent.infoCapData)
      for i ∈ 1:o.kZ
        o.FillingT₀[1,i+1] = tmp[:,i]
      end
      o.ScapPXYZperp[1] = panelsum(o.Zperp, uwt, o.parent.infoCapData)  # Scap(P_(MZperpX) * y₁ .* Zperp)

      o.parent.NFE>0 &&
        (o.CT_FEcapPY[1] = crosstabFE(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)

      uwt = vHadw(o.y₁par, o.parent.wt)
      tmp = panelsum(o.PXZ, uwt, o.parent.infoCapData)
      for i ∈ 1:o.kZ
        o.FillingT₀[i+1,1] = tmp[:,i]
      end
      o.FillingT₀[1] = panelsum(      o.PXy₁, uwt, o.parent.infoCapData)
      !o.parent.granular &&
        (o.ScapYX[1] = panelsum(o.X₁, o.X₂  , uwt, o.parent.infoCapData))  # Scap(M_Zperp*y₁ .* P_(MZperpX)])
    end
  end
end


@inline function MakeResiduals!(o::StrEstimator{T,<:Union{OLS,ARubin}} where T)
  o.ü₁ = o.y₁par - X₁₂B(o.parent.X₁, o.parent.X₂, o.β)
end

function MakeResiduals!(o::StrEstimator{T,IVGMM} where T)
  o.ü₁ = o.y₁par - o.Z * o.β

  if !o.parent.scoreBS
    _β = [1 ; -o.β]
    uu = _β'o.YY * _β

    Xu = o.Xy₁par - o.XZ * o.β  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
    negXuinvuu = Xu / -uu
    o.Ü₂ = o.Y₂ - X₁₂B(o.X₁, o.X₂, invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β)' + o.XY₂))  # large expression is Pihat

    o.u⃛₁ = o.ü₁ + o.Ü₂ * (o.t₁Y + o.RparY * o.β)
  end
end


# non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
# since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
# depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then κ from DGP LIML affects Hessian, H.
function InitTestDenoms!(o::StrEstimator)
  if o.parent.bootstrapt && (o.parent.scoreBS || o.parent.robust)
    (o.parent.granular || o.parent.purerobust) && (o.WXAR = vHadw(o.XAR, o.parent.wt))

    if o.parent.robust && o.parent.NFE>0 && !(o.parent.FEboot || o.parent.scoreBS) && Int(o.parent.granular) < o.parent.NErrClustCombs  # make first factor of second term of (64) for c=cap (c=1)
      !isdefined(o, :WXAR) && (o.WXAR = vHadw(o.XAR, o.parent.wt))
      o.CT_XAR = [crosstabFE(o.parent, view(o.WXAR,:,d), o.parent.infoCapData) for d ∈ axes(o.WXAR,2)]
    end
	end
end


# partial fixed effects out of a data matrix
function partialFE!(o::StrBoottest, In::AbstractArray)
	if length(In)>0
    for f ∈ o.FEs
			tmp = @view In[f.is,:]
			tmp .-= f.wt'tmp
		end
	end
	return In
end

function setdirty!(o::StrBoottest, _dirty::Bool; noinitialize::Bool=false)
	o.dirty = _dirty
	!noinitialize && (o.initialized = false)
end

function setsqrt!(o::StrBoottest, _sqrt::Bool)
	if _sqrt < o.sqrt
		if !o.dirty
    	Dist .^= 2
      multiplier ^= 2
    end
	else
		setdirty!(o, true)
  end
	o.sqrt = _sqrt
end

function setptype!(o::StrBoottest, ptype::PType)
  o.ptype = ptype
	o.twotailed = o.ptype==symmetric || o.ptype==equaltail
end

function setstattype!(o::StrBoottest, stattype::StatType)
  o.bootstrapt = stattype == t
	setdirty!(o, true)
end

function setX₁!(o::StrBoottest{T}, X₁::AbstractMatrix) where T
	o.X₁ = eltype(X₁)==T ? X₁ : T.(X₁); setdirty!(o, true)
end
function setX₂!(o::StrBoottest{T}, X₂::AbstractMatrix) where T
	o.X₂ = eltype(X₂)==T ? X₂ : T.(X₂); setdirty!(o, true)
end
function sety₁!(o::StrBoottest{T}, y₁::AbstractVector) where T
	o.y₁ = eltype(y₁)==T ? y₁ : T.(y₁); setdirty!(o, true)
end
function setY₂!(o::StrBoottest{T}, Y₂::AbstractMatrix) where T
	o.Y₂ = eltype(Y₂)==T ? Y₂ : T.(Y₂); setdirty!(o, true)
end
setY₂!(o::StrBoottest, Y₂::AbstractVector) = setY₂!(o,Y₂[:,:])
function setobswt!(o::StrBoottest{T}, wt::AbstractVector, fweights::Bool) where T
	o.wt = eltype(wt)==T ? wt : T.(wt); o.fweights = fweights; setdirty!(o, true)
end
function setsc!(o::StrBoottest{T}, Sc::AbstractMatrix) where T
	o.Sc = eltype(Sc)==T ? Sc : T.(Sc); setdirty!(o, true)
end
function setML!(o::StrBoottest; ML::Bool)
	o.ML  = ML; setdirty!(o, true)
	ML && setscoreBS(true)
end
function setLIML!(o::StrBoottest, LIML::Bool)
	o.LIML = LIML; setdirty!(o, true)
end
function setARubin!(o::StrBoottest, ARubin::Bool)
	o.ARubin = ARubin; setdirty!(o, true)
end
function setFuller!(o::StrBoottest, Fuller::Number)
	o.Fuller = Fuller; setdirty!(o, true)
end
function setκ!(o::StrBoottest, κ::Number)  # κ as in k-class
	o.κ = κ; setdirty!(o, true)
end
function setquietly!(o::StrBoottest, quietly::Bool)
	o.quietly = quietly
end
function setβ!(o::StrBoottest, β::AbstractVector)
	o.β = β; setdirty!(o, true)
end
function setA!(o::StrBoottest, A::AbstractMatrix)
	o.A = A; setdirty!(o, true)
end
function setsmall!(o::StrBoottest, small::Bool)
	o.small = small; setdirty!(o, true)
end
function setscoreBS!(o::StrBoottest, scoreBS::Bool)
	o.scoreBS = scoreBS; setdirty!(o, true)
end
function setB!(o::StrBoottest, B::Integer)
	o.B = B
	B==0 && setscoreBS!(o, true)
	setdirty!(o, true)
end
function setnull!(o::StrBoottest, null::Bool)
	o.null = null; setdirty!(o, true)
end
function setWald!(o::StrBoottest) # set-up for classical Wald test
	o.scoreBS = true; o.B = 0; o.null = false; setdirty!(o, true)
end
function setRao!(o::StrBoottest) # set-up for classical Rao test
	o.scoreBS = true; o.B = 0; o.null = true; setdirty!(o, true)
end
function setID!(o::StrBoottest{T}, ID::AbstractArray; NBootClustVar::Integer=1, NErrClust::Integer=NBootClustVar) where T
	o.ID = (eltype(ID)==T ? ID : T.(ID))[:,:]; o.NBootClustVar = NBootClustVar; o.NErrClust = NErrClust; setdirty!(o, true)
	length(ID)>0 && (o.robust = true)
end
function setFEID!(o::StrBoottest{T}, ID::AbstractVector; NFE::Integer=-1, FEdfadj::Bool=true) where T
	o.FEID = eltype(ID)==T ? ID : T.(ID); o.NFE = NFE; o.FEdfadj = FEdfadj; setdirty!(o, true)
end

setlevel!(o::StrBoottest{T}, level::T) where T = (o.level = level)
setptol!(o::StrBoottest{T}, rtol::T) where T = (o.rtol = rtol)  # call "ptol" in Mata
setwillplot!(o::StrBoottest, willplot::Bool) = (o.willplot = willplot)

function setrobust!(o::StrBoottest{T} where T, robust::Bool)
	o.robust = robust
	!o.robust && setID!(o, zeros(T,0,0), 1, 1)
	setdirty!(o, true)
end
function setR₁!(o::StrBoottest, R₁::AbstractMatrix, r₁::AbstractVector)
	o.R₁ = R₁; o.r₁  = r₁; setdirty!(o, true)
end
function setR!(o::StrBoottest, R::AbstractMatrix, r::AbstractVector)
	o.R = R; o.r = r; o.q = rows(R); setdirty!(o, true)  # q can differ from df in ARubin test
end
function setgrid!(o::StrBoottest, gridmin::Union{Vector{T1},Vector{Union{T1,Missing}}}, gridmax::Union{Vector{T2},Vector{Union{T2,Missing}}}, gridpoints::Union{Vector{T3},Vector{Union{T3,Missing}}}) where {T1<:Real, T2<:Real, T3<:Real}
	o.gridmin = gridmin; o.gridmax = gridmax; o.gridpoints = Vector{Union{Int32,Missing}}(gridpoints)
end
function setmadjust!(o::StrBoottest, madjtype::MAdjType, NumH0s::Integer)
	o.NumH0s = NumH0s
  o.madjtype = madjtype
end
function setauxwttype!(o::StrBoottest, auxtwtype::AuxWeightType)
	o.auxtwtype = auxtwtype; setdirty!(o, true)
end
function setrng!(o::StrBoottest, rng::AbstractRNG)
	o.rng = rng
end
function setMaxMatSize!(o::StrBoottest, MaxMatSize::Number)
	o.MaxMatSize = MaxMatSize; setdirty!(o, true)
end
setauxwttype!
function getdist(o::StrBoottest; diststat::String="")
	o.dirty && boottest!(o)
	if diststat == "numer"
		_numer = isone(o.u_sd) ? o.numer : o.numer / o.u_sd
		o.DistCDR = (@view _numer[:,2:end]) .+ o.r
    sort!(o.DistCDR)  # need to specify horizontal sort??
	elseif length(o.DistCDR)==0
		if length(o.Dist) > 1
      o.DistCDR = (@views o.Dist[2:end])[:,:] * o.multiplier
      sort!(o.DistCDR, dims=1)  # need to specify horizontal sort??
		else
			o.DistCDR = zeros(0,1)
    end
  end
	o.DistCDR
end

# get p valuo. Robust to missing bootstrapped values interpreted as +infinity.
function getp(o::StrBoottest{T}; classical::Bool=false) where T
	o.dirty && boottest!(o)
	tmp = o.Dist[1]
	isnan(tmp) && return tmp
	if o.B>0 && !classical
		if o.sqrt && o.ptype ≠ upper
			if o.ptype==symmetric
				o.p = sum(-abs(tmp) .> -abs.(o.Dist)) / o.BFeas  # symmetric p value; do so as not to count missing entries in *Dist
			elseif o.ptype==equaltail
				o.p = 2 * min(sum(tmp .> o.Dist) , sum(-tmp .> -o.Dist)) / o.BFeas
			else
				o.p = sum(tmp .>  o.Dist) / o.BFeas  # lower-tailed p value
      end
		else
			o.p = sum(-tmp .> -o.Dist) / o.BFeas  # upper-tailed p value or p value based on squared stats
    end
	else
		tmp *= o.multiplier
    _p = ccdf(o.small ? FDist(o.df, o.df_r) : Chisq(o.df), Float64(o.sqrt ? tmp^2 : tmp))
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
function getb(o::StrBoottest)
	o.dirty && boottest!(o)
	@views isone(o.u_sd) ? o.numer[:,1] : o.numer[:,1] / o.u_sd
end

# denominator for full-sample test stat
function getV(o::StrBoottest)
	o.dirty && boottest!(o)
	o.statDenom / ((isone(o.u_sd) ? o.smallsample : o.u_sd^2 * o.smallsample)  * (o.sqrt ? o.multiplier^2 : o.multiplier) * o.df)
end

# wild weights
getv(o::StrBoottest) = @views isone(o.u_sd) ? o.v[:,2:end] : o.v[:,2:end] / o.u_sd

# Return number of bootstrap replications with feasible results
# Returns 0 if getp() not yet accessed, or doing non-bootstrapping tests
getrepsFeas(o::StrBoottest) = o.BFeas
getNBootClust(o::StrBoottest) = o.Nstar
getreps(o::StrBoottest) = o.B  # return number of replications, possibly reduced to 2^G

function getpadj(o::StrBoottest{T} where T; classical::Bool=false)
	_p = o.dirty | classical ? getp(o, classical=classical) : o.p
	if o.madjtype==bonferroni min(one(T), o.NumH0s * _p)
  elseif o.madjtype==sidak  one(T) - (one(T) - _p) ^ o.NumH0s
  else _p
	end
end

function getstat(o::StrBoottest)
	o.dirty && boottest!(o)
	o.multiplier * o.Dist[1]
end
function getdf(o::StrBoottest)
	o.dirty && boottest!(o)
	o.df
end
function getdf_r(o::StrBoottest)
	o.dirty && boottest!(o)
	o.df_r
end
function getplot(o::StrBoottest)
	o.notplotted && plot(o)
	(o.plotX, o.plotY)
end
function getpeak(o::StrBoottest)  # x and y values of confidence curve peak (at least in OLS && ARubin)
  o.notplotted && plot(o)
	o.peak
end
function getCI(o::StrBoottest)
	o.notplotted && plot(o)
	o.CI
end

macro storeWtGrpResults!(dest, content)  # poor hygiene in referencing caller's o and w
  if dest == :(o.Dist)
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

function Init!(o::StrBoottest{T}) where T  # for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  o.Nobs = rows(o.X₁)
  o.NClustVar = cols(o.ID)
  o.kX = (o.kX₁ = cols(o.X₁)) + (o.kX₂ = cols(o.X₂))
  o.kX₂==0 && (o.X₂ = zeros(T,o.Nobs,0))
  (o.kY₂ = cols(o.Y₂))==0 && (o.Y₂ = zeros(T,o.Nobs,0))
  o.kZ  = o.kX₁ + o.kY₂
	if o.LIML && o.kX₂==o.kY₂  # exactly identified LIML = 2SLS
		o.κ = one(T)
		o.LIML = false
	end
	if !(o.REst = length(o.R₁)>0)  # base model contains no restrictions?
    o.R₁ = zeros(T,0,o.kZ)
    o.r₁ = zeros(T,0)
  end
  isnan(o.κ) && (o.κ = o.kX₂>0 ? one(T) : zero(T))  # if κ in κ-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  o.WRE = !(iszero(o.κ) || o.scoreBS) || o.ARubin
  o.WREnonARubin = o.WRE && !o.ARubin

  o.haswt = typeof(o.wt) <: AbstractVector
  if o.haswt
    o.sumwt = sum(o.wt)
  else
    o.wt = I
    o.sumwt = 1
  end
  o._Nobs = o.haswt && o.fweights ? o.sumwt : o.Nobs

  if o.WREnonARubin
    if o.NClustVar>0
      o.infoBootData, o.IDBootData = panelsetupID(o.ID, collect(1:o.NBootClustVar))
    else
      o.infoCapData = o.infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # no clustering, so no collapsing by cluster
    end
  elseif o.NClustVar>0
    o.infoBootData = panelsetup(o.ID, collect(1:min(o.NClustVar,o.NBootClustVar)))  # bootstrap cluster grouping defs rel to original data
  else
    infoCapData = infoAllData = infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # causes no collapsing of data in panelsum() calls, only multiplying by weights if any
  end
  o.Nstar = rows(o.infoBootData)

  if o.bootstrapt
    if !iszero(o.NClustVar)
      minN = Inf; sumN = 0

      Combs = combs(o.NErrClust)  # represent all error clustering combinations. First is intersection of all error clustering vars
      o.clust = Vector{StrClust{T}}(undef, rows(Combs)-1)  # leave out no-cluster combination
      o.NErrClustCombs = length(o.clust)
      o.subcluster = o.NClustVar - o.NErrClust

      if o.NClustVar > o.NBootClustVar  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
        if o.WREnonARubin && !o.granular
          o.infoAllData, IDAllData = panelsetupID(o.ID, collect(1:o.NClustVar))
        else
          o.infoAllData            = panelsetup(o.ID, collect(1:o.NClustVar))
        end
      else
        o.infoAllData = o.infoBootData  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
        o.WREnonARubin && !o.granular && (IDAllData = o.IDBootData)
      end

      if o.NClustVar > o.NErrClust  # info for intersections of error clustering wrt data
        if o.WREnonARubin && !o.granular
          o.infoCapData, IDCapData = panelsetupID(o.ID, collect(o.subcluster+1:o.NClustVar))
        else
          o.infoCapData            = panelsetup(o.ID, collect(o.subcluster+1:o.NClustVar))
        end
        IDCap = rows(o.infoCapData)==o.Nobs ? o.ID : @views o.ID[first.(o.infoCapData),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
        o.IDAll = rows(o.infoAllData)==o.Nobs ? o.ID : @views o.ID[first.(o.infoAllData),:]  # version of ID matrix with one row for each all-bootstrap && error cluster-var intersection instead of 1 row for each obs
      else
        o.infoCapData = o.infoAllData  # info for intersections of error clustering wrt data
        o.WREnonARubin && !o.granular && (IDCapData = IDAllData)
        o.IDAll = IDCap = rows(o.infoCapData)==o.Nobs ? o.ID : @views o.ID[first.(o.infoCapData),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
      end

      o.BootClust = 2^(o.NClustVar - o.NBootClustVar)  # location of bootstrap clustering within list of cluster combinations

			for c ∈ 1:o.NErrClustCombs  # for each error clustering combination
        ClustCols = o.subcluster .+ findall(@view Combs[c,:])
        even = isodd(length(ClustCols))  # not a typo

        if isone(c)
          if iszero(o.subcluster)
            order = Vector{Int64}(undef,0)
            info  = Vector{UnitRange{Int64}}(undef, rows(o.infoAllData))  # causes no collapsing of data in panelsum() calls
          else
            order = sortperm(collect(eachrow(@view IDCap[:,ClustCols])))  # XXX slow?
            IDCap = IDCap[order, :]
            info  = panelsetup(IDCap, ClustCols)
          end
        else
          if any(Combs[c, min(findall(Combs[c,:] .≠ Combs[c-1,:])...):end])  # if this sort ordering same as last to some point and missing thereafter, no need to re-sort
            order = sortperm(collect(eachrow(@view IDCap[:,ClustCols])))  # XXX slow?
            IDCap = IDCap[order,:]
          else
            order = Vector{Int64}(undef,0)
          end
          info = panelsetup(IDCap, ClustCols)
        end

        N = rows(info)
        sumN += N

        if o.small
          multiplier = T(N / (N-1))
          N < minN && (minN = N)
        else
          multiplier = one(T)
        end
        o.clust[c] = StrClust(N, multiplier, even, order, info)
      end

      (o.scoreBS || !o.WREnonARubin) &&
        (o.ClustShare = o.haswt ? panelsum(o.wt, o.infoCapData)/o.sumwt : length.(o.infoCapData)./o.Nobs) # share of observations by group

    else  # if no clustering, cast "robust" as clustering by observation
      clust = StrClust{T}(Nobs, small ? _Nobs / (_Nobs - 1) : 1, true, Vector{Int64}(undef,0), Vector{UnitRange{Int64}}(undef,0))
      sumN = o.Nobs
      o.NErrClustCombs = 1
      (o.scoreBS || !o.WREnonARubin) &&
        (o.ClustShare = o.haswt ? o.wt/o.sumwt : 1/o._Nobs)
    end

    o.purerobust = o.robust && !o.scoreBS && iszero(o.subcluster) && o.Nstar==o.Nobs  # do we ever error-cluster *and* bootstrap-cluster by individual?
    o.granular   = o.WREnonARubin ? 2*o.Nobs*o.B*(2*o.Nstar+1) < o.Nstar*(o.Nstar*o.Nobs+o.clust[1].N*o.B*(o.Nstar+1)) :
		                                o.NClustVar>0 && !o.scoreBS && (o.purerobust || (o.clust[1].N+o.Nstar)*o.kZ*o.B + (o.clust[1].N-o.Nstar)*o.B + o.kZ*o.B < o.clust[1].N*o.kZ^2 + o.Nobs*o.kZ + o.clust[1].N * o.Nstar * o.kZ + o.clust[1].N * o.Nstar)

    if o.robust && !o.purerobust
      (o.subcluster>0 || o.granular) &&
        (o.infoErrAll = panelsetup(o.IDAll, collect(Int(o.subcluster)+1:o.NClustVar)))  # info for error clusters wrt data collapsed to intersections of all bootstrapping && error clusters; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusterings
      ((o.scoreBS && o.B>0) || (o.WREnonARubin && !o.granular && o.bootstrapt)) &&
        (o.JNcapNstar = zeros(T, o.clust[1].N, o.Nstar))
    end

    if o.WREnonARubin && o.robust && o.bootstrapt && !o.granular
      if iszero(length(IDAllData))
        _, IDAllData = panelsetupID(o.ID, collect(                  1:o.NClustVar))
        _, IDCapData = panelsetupID(o.ID, collect(Int(o.subcluster)+1:o.NClustVar))
      end
      o.IDCTCapstar   = Vector{Vector{Int64}}(undef, o.Nstar)
      o.infoCTCapstar = Vector{Vector{UnitRange{Int64}}}(undef, o.Nstar)
      for i ∈ 1:o.Nstar
        tmp = IDAllData[o.infoBootData[i]]                        # ID numbers w.r.t. intersection of all bootstrap/error clusterings contained in bootstrap cluster i
        o.infoCTCapstar[i] = o.infoAllData[tmp[1]:tmp[end]]       # for each of those ID's, panel info for the all-bootstrap/error-clusterings data row groupings
        o.IDCTCapstar[i] = IDCapData[first.(o.infoCTCapstar[i])]  # ID numbers of those groupings w.r.t. the all-error-clusterings grouping
      end
    end
  else
    minN = rows(o.infoBootData)
  end

  if isdefined(o, :FEID)
    p = sortperm(o.FEID)
    sortID = o.FEID[p]
    i_FE = 1; o.FEboot = o.B>0 && !o.WREnonARubin && o.NClustVar>0; j = o.Nobs; o._FEID = ones(Int64, o.Nobs)
    o.invFEwt = zeros(T, o.NFE>0 ? o.NFE : o.Nobs)
    o.FEs = Vector{StrFE}(undef, o.NFE>0 ? o.NFE : o.Nobs)
    @inbounds for i ∈ o.Nobs-1:-1:1
      if sortID[i] ≠ sortID[i+1]
        is = @view p[i+1:j]
        if o.haswt
          tmp  = o.wt[is]
          wt = tmp / (sumFEwt = sum(tmp))
        else
          sumFEwt = j - i
          wt = fill(1/sumFEwt, j-i)
        end
        o.FEs[i_FE] = StrFE(is, wt)
        if (o.B>0 && o.robust && Int(o.granular) < o.NErrClust) || (o.WREnonARubin && o.robust && o.granular && o.bootstrapt)
          o.invFEwt[i_FE] = 1 / sumFEwt
        end

        j = i

        if o.FEboot  # are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
          tmp = o.ID[is, 1:o.NBootClustVar]
          o.FEboot = all(tmp .== view(tmp, 1,:)')
        end
        i_FE += 1
      end
      o._FEID[p[i]] = i_FE
    end
    is = @view p[1:j]
    if o.haswt
      tmp = o.wt[is]
      wt = tmp / (sumFEwt = sum(tmp))
    else
      sumFEwt = j
      wt = fill(1/sumFEwt, j)
    end
    o.FEs[i_FE] = StrFE(is, wt)
    o.robust && ((o.B>0 && Int(o.granular) < o.NErrClust) || (o.WREnonARubin && o.granular && o.bootstrapt)) &&
      (o.invFEwt[i_FE] = 1 / sumFEwt)
    o.NFE = i_FE
    resize!(o.invFEwt, o.NFE)
    resize!(o.FEs, o.NFE)
    if o.FEboot  # are all of this FE's obs in same bootstrapping cluster?
      tmp = o.ID[is, 1:o.NBootClustVar]
      o.FEboot = all(tmp .== @view tmp[1,:])
    end

    if o.robust && o.B>0 && o.bootstrapt && !o.FEboot && Int(o.granular) < o.NErrClust
      o.infoBootAll = panelsetup(o.IDAll, collect(1:o.NBootClustVar))  # info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping && error clusters
    end

		partialFE!(o, o.X₁)
		partialFE!(o, o.X₂)
		partialFE!(o, o.y₁)
		partialFE!(o, o.Y₂)
	end

  if o.B>0 && o.robust && o.granular && !o.purerobust && !o.bootstrapt && !o.WREnonARubin
    if o.NFE>0 && !o.FEboot
      _, o.IDBootData = panelsetupID(o.ID   , 1:o.NBootClustVar)
    else
      _, o.IDBootAll  = panelsetupID(o.IDAll, 1:o.NBootClustVar)
    end
  end

  o.enumerate = o.B>0 && o.auxtwtype==rademacher && o.Nstar*log(2) < log(o.B)+1e-6  # generate full Rademacher set?
  o.enumerate && (o.MaxMatSize = 0)

  o.Nw = iszero(o.MaxMatSize) ? 1 : ceil((o.B+1) * Float64(max(rows(o.IDBootData), length(o.IDBootAll), o.Nstar) * sizeof(T)) / o.MaxMatSize / 1073741824) # 1073741824 = giga(byte)
  if isone(o.Nw)
    MakeWildWeights!(o, o.B, first=true)  # make all wild weights, once
    o.enumerate && (o.B = cols(o.v) - 1)  # replications reduced to 2^G
    o.WeightGrp = [1:cols(o.v)]
  else
    _B = ceil(Int64, (o.B+1) / o.Nw)
    o.Nw = ceil(Int64, (o.B+1) / _B)
    o.WeightGrp = [(i-1)*_B+1:i*_B for i ∈ 1:o.Nw]
    o.WeightGrp[end] = first(o.WeightGrp[end]):o.B+1
  end

  if o.ML
    o.df = rows(o.R)
  else
    if o.ARubin
      o.R = ApplyArray(hcat, zeros(o.kX₂,o.kX₁), I(o.kX₂))  # attack surface is all endog vars
      o.R₁ = o.kX₁>0 && rows(o.R₁)>0 ? ApplyArray(hcat, o.R₁[:,1:kX₁], zeros(rows(o.R₁),o.kX₂)) : zeros(0, o.kX)  # and convert model constraints from referring to X₁, Y₂ to X₁, X₂
    end
    o.df = rows(o.R)

    if !o.WRE && iszero(o.κ)  # regular OLS
      o.DGP = StrEstimator{T,OLS}(o)
      o.Repl = StrEstimator{T,OLS}(o)
      o.DGP.LIML = o.LIML; o.DGP.Fuller = o.Fuller; o.DGP.κ = o.κ
      setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁)  # DGP constraints: model constraints + null if imposed
      setR!(o.Repl, o.R₁)  # model constraints only
      InitVars!(o.DGP, o.Repl.R₁perp)
      InitTestDenoms!(o.DGP)
      o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
    elseif o.ARubin
      if o.willplot  # for plotting/CI purposes get original point estimate since not normally generated
        o.DGP = StrEstimator{T,IVGMM}(o)
        o.DGP.LIML = o.LIML; o.DGP.Fuller = o.Fuller; o.DGP.κ = o.κ
        setR!(o.DGP, o.R₁, zeros(T,0,o.kZ))  # no-null model
        InitVars!(o.DGP)
        Estimate!(o.DGP, o.r₁)
        o.confpeak = o.DGP.β  # estimated coordinate of confidence peak
      end

      o.DGP = StrEstimator{T,ARubin}(o)
      setR!(o.DGP, o.R₁)
      InitVars!(o.DGP)
      InitTestDenoms!(o.DGP)
      o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
      o.kZ = o.kX

    elseif o.WREnonARubin

      o.DGP = StrEstimator{T,IVGMM}(o)
      o.DGP.LIML = o.kX₂ ≠ o.kY₂
      setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁, zeros(T,0,o.kZ))  # DGP constraints: model constraints + null if imposed
      InitVars!(o.DGP)
      if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
      	Estimate!(o.DGP, o.r₁)
        MakeResiduals!(o.DGP)
      end

      o.Repl = StrEstimator{T,IVGMM}(o)
      o.Repl.isDGP = false
      o.Repl.LIML, o.Repl.Fuller, o.Repl.κ = o.LIML, o.Fuller, o.κ
      setR!(o.Repl, o.R₁, o.R)
      InitVars!(o.Repl)
      Estimate!(o.Repl, o.r₁)

      o.LIML && o.Repl.kZ==1 && o.Nw==1 && (o.As = o.βs = zeros(1, o.B+1))
      o.SstarUZperpinvZperpZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
      o.SstarUZperp              = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
      o.SstaruY                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
      o.SstarUXinvXX             = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
      o.SstarUX                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
      o.SstarUU                  = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
      if o.bootstrapt
        o.δdenom_b = zeros(o.Repl.kZ, o.Repl.kZ)
        o.SstarUMZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
        o.SstarUPX     = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
        o.δdenom   = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
        o._Jcap = zeros(o.clust[1].N, o.Repl.kZ)
        !o.granular && (o.SCTcapuXinvXX = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Nstar))
        if o.LIML || !o.robust
          YYstar_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
          YPXYstar_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
        end
        o.NFE>0 && (o.bootstrapt || !isone(o.κ) || o.LIML) && (o.CTFEU = Vector{Matrix{T}}(undef, o.Repl.kZ+1))
      end

    else  # the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R₁, one with FWL, one not

      o.DGP = StrEstimator{T,IVGMM}(o)
      o.DGP.LIML = o.LIML; o.DGP.Fuller = o.Fuller; o.DGP.κ = o.κ
      setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁, zeros(T,0,o.kZ))  # DGP constraints: model constraints + null if imposed
      InitVars!(o.DGP)
      o.Repl = StrEstimator{T,IVGMM}(o)
      o.Repl.LIML = o.LIML; o.Repl.Fuller = o.Fuller; o.Repl.κ = o.κ
      setR!(o.Repl, o.R₁, I)  # process replication restraints = model constraints only
      InitVars!(o.Repl, o.Repl.R₁perp)
      Estimate!(o.Repl, o.r₁)  # bit inefficient to estimate in both objects, but maintains the conformity
      InitTestDenoms!(o.Repl)
      o.M = o.Repl  # StrEstimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scoreBS for IV/GMM mixes the two
      if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
      	Estimate!(o.DGP, o.r₁)
        MakeResiduals!(o.DGP)
      end
    end
  end

  if !o.WREnonARubin && o.bootstrapt
    o.denom = [Matrix{T}(undef,0,0) for _ in 1:o.df, _ in 1:o.df]
    if o.robust
      o.Kcd =                       Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.df)
      o.Jcd = iszero(o.B) ? o.Kcd : Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.df)  # if B = 0, Kcd will be multiplied by v, which is all 1's, and will constitute Jcd
    end
  end

  o.small && (o.df_r = o.NClustVar>0 ? minN - 1 : o._Nobs - o.kZ - o.NFE)

  isone(o.df) && setsqrt!(o,true)  # work with t/z stats instead of F/chi2

  if o.small
    o.multiplier = (o.smallsample = (o._Nobs - o.kZ - Int(o.FEdfadj) * o.NFE) / (o._Nobs - Int(o.robust))) / o.df  # divide by # of constraints because F stat is so defined
  else
    o.multiplier = o.smallsample = 1
  end

  !(o.robust || o.ML) && (o.multiplier *= o._Nobs)  # will turn sum of squared errors in denom of t/z into mean
  o.sqrt && (o.multiplier = √o.multiplier)

  o.bootstrapt && (o.WREnonARubin || o.df>2 || !isnan(o.MaxMatSize)) && # unless nonWRE or df=1 or splitting weight matrix, code will create Dist element-by-element, so pre-allocate vector now
    (o.Dist = fill(T(NaN), o.B+1))
  (o.Nw>1 || o.WREnonARubin || (!o.null && o.df<=2)) && (o.numer = fill(T(NaN), o.df, o.B+1))

  if !o.WREnonARubin
    poles = anchor = zeros(T,0,0)
    o.interpolable = o.bootstrapt && o.B>0 && o.null && o.Nw==1 && (iszero(o.κ) || o.ARubin)
    if o.interpolable
      o.∂numer∂r = Vector{Matrix{T}}(undef, o.q)
      o.interpolate_u = !(o.robust || o.ML)
      o.interpolate_u && (o.∂u∂r = Vector{Matrix{T}}(undef, o.q))
      if o.robust
        o.∂denom∂r   = [Matrix{Matrix{T}}(undef, o.df, o.df) for _ in 1:o.q]
        o.∂²denom∂r² = [Matrix{Matrix{T}}(undef, o.df, o.df) for _ in 1:o.q, _ in 1:o.q]
        o.∂Jcd∂r     = [Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.df) for _ in 1:o.q]
      end
    end
  end
end

# main routine
function boottest!(o::StrBoottest{T}) where T
	if !o.initialized
    Init!(o)
  elseif !o.null
    NoNullUpdate!(o)
    return
  end

	o.Nw > 1 && MakeWildWeights!(o, last(o.WeightGrp[1])-1, first=true)

  o.WREnonARubin ? PrepWRE!(o) :
                   MakeInterpolables!(o)  # make stuff that depends linearly on r, possibly by interpolating, for first weight group

  for w ∈ 1:o.Nw  # do group 1 first because it includes col 1, which is all that might need updating in constructing CI in WCU
		w > 1 && MakeWildWeights!(o, length(o.WeightGrp[w]), first=false)

		o.WREnonARubin ? MakeWREStats!(o, w) :
			               MakeNonWREStats!(o, w)

    !o.bootstrapt && UpdateBootstrapcDenom!(o, w)
	end

  o.BFeas = isnan(o.Dist[1]) ? 0 : sum(.!(isnan.(o.Dist) .| isinf.(o.Dist))) - 1
	o.DistCDR = zeros(T,0,0)
	setdirty!(o, false, noinitialize=false)
	o.initialized = true
end

# if not imposing null and we have returned to boottest!(), then df=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
function NoNullUpdate!(o::StrBoottest{T} where T)
  if o.WREnonARubin
    o.numer[:,1] = o.R * o.DGP.Rpar * o.βs[1] - o.r
  elseif o.ARubin
    o.DGP.Estimate(o.r)
    o.numer[:,1] = o.u_sd * o.DGP.Rpar * @view o.DGP.β[o.kX₁+1:end,:] # coefficients on excluded instruments in ARubin OLS
  else
    o.numer[:,1] = o.u_sd * (o.R * (o.ML ? o.β : o.M.Rpar * o.M.β) - o.r) # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this
  end
  o.Dist[1] = isone(o.df) ? o.numer[1] / sqrt(o.statDenom[1]) : o.numer[:,1]'invsym(o.statDenom)*o.numer[:,1]
end

# compute bootstrap-c denominator from all bootstrap numerators
function UpdateBootstrapcDenom!(o::StrBoottest{T} where T, w::Integer)
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
		o.Dist = o.sqrt ? vec(o.numer ./ sqrtNaN.(o.statDenom)) : colquadform(invsym(o.statDenom), o.numer)
	end
end

# draw wild weight matrix of width _B. If first=true, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all weights
const ϕ = (1 + √5) * .5

function MakeWildWeights!(o::StrBoottest{T}, _B::Integer; first::Bool=true) where T
	if _B>0  # in scoretest or waldtest WRE, still make v a col of 1's
		if o.enumerate
			o.v = [ones(o.Nstar) count_binary(o.Nstar, -1-Int(o.WREnonARubin), 1-Int(o.WREnonARubin))]  # complete Rademacher set
		else
      if o.auxtwtype==normal
        o.v = T==Float64 ? randn(o.rng, o.Nstar, _B+first) : T.(randn(o.rng, o.Nstar, _B+first))
      elseif o.auxtwtype==gamma
        o.v = T==Float64 ? quantile.(Gamma(4,.5), rand(o.rng, o.Nstar, _B+first)) : T.(quantile.(Gamma(4,.5), rand(o.rng, o.Nstar, _B+first)))
      elseif o.auxtwtype==webb
        o.v = getindex.(Ref(T.([-√1.5, -1, -√.5, √.5, 1, √1.5])), ceil.(Int8, 6 .* rand(o.rng, o.Nstar, _B+first)))
      elseif o.auxtwtype == mammen
        o.v = getindex.(Ref(T.([1-ϕ; ϕ])), ceil.(Int8, rand(o.rng, o.Nstar, _B+first) ./ (ϕ/√5)))
      else
        o.v = 2 .* (rand(o.rng, o.Nstar, _B+first) .> .5) .- 1   # Rademacher
      end
      o.WREnonARubin && (o.v .-= 1)  # for efficiency tweak pre-subtract 1 from all weights for WRE
    end

		first && (o.v[:,1] = fill(o.WREnonARubin ? 0 : o.u_sd, o.Nstar))  # XXX unnec when enumerated and u_sd=1. keep original residuals in first entry to compute base model stat
	else
		o.v = Matrix{T}(undef,0,1)  # in places, cols(v) indicates B -- 1 for classical tests
  end
end


# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of Y[:,ind1]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,ind2] for κ constant across replications
# ind1 can be a rowvector
# (only really the Hessian when we narrow Y to Z)
# HessianFixedκ(o::StrBoottest, ind1::Vector{T} where T<:Integer, ind2::Integer, κ::Number) = vcat([_HessianFixedκ(o, ind1[i], ind2, κ) for i ∈ ind1]...)

function HessianFixedκ(o::StrBoottest{T}, ind1::Vector{S} where S<:Integer, ind2::Integer, κ::Number) where T
  retval = Matrix{T}(undef, length(ind1), cols(o.v))
  for i ∈ eachindex(ind1, axes(retval,1))
    retval[i,:] = _HessianFixedκ(o, ind1[i], ind2, κ)
  end
  retval
end

function _HessianFixedκ(o::StrBoottest, ind1::Integer, ind2::Integer, κ::Number)
  if !iszero(κ)
    T1L = iszero(ind1) ? o.Repl.Xy₁par : o.Repl.XZ[:,ind1]
    if o.Repl.Yendog[ind1+1]
		tmp = o.SstarUX[ind1+1] * o.v
		tmp .+= T1L
		T1L = tmp
	end
		T1R = iszero(ind2) ? o.Repl.invXXXy₁par : o.Repl.invXXXZ[:,ind2]
		if o.Repl.Yendog[ind2+1]
			tmp = o.SstarUXinvXX[ind2+1] * o.v
			tmp .+= T1R
			T1R = tmp  # right-side linear term
		end

    retval = coldot(T1L, T1R)  # multiply in the left-side linear term
	end

	if !isone(κ)
    if o.Repl.Yendog[ind1+1]
      T2 = o.SstarUZperpinvZperpZperp[ind1+1]'o.SstarUZperp[ind2+1]  # quadratic term
      T2[diagind(T2)] .-= ind1 <= ind2 ? o.SstarUU[ind2+1, ind1+1] : o.SstarUU[ind1+1, ind2+1]  # minus diagonal crosstab
      o.NFE>0 &&
        (T2 .+= o.CTFEU[ind1+1]'(o.invFEwt .* o.CTFEU[ind2+1]))

      retval = iszero(κ) ? o.Repl.YY[ind1+1,ind2+1] .+ colquadformminus!((                            @views o.SstaruY[ind2+1][:,ind1+1] .+ o.SstaruY[ind1+1][:,ind2+1])'o.v, T2, o.v) :
                           κ .* retval .+ (1 - κ)   .* colquadformminus!((o.Repl.YY[ind1+1,ind2+1] .+ @views o.SstaruY[ind2+1][:,ind1+1] .+ o.SstaruY[ind1+1][:,ind2+1])'o.v, T2, o.v)
    else
      retval = iszero(κ) ? hcat(o.Repl.YY[ind1+1,ind2+1]) : κ .* retval .+ (1 - κ) .* o.Repl.YY[ind1+1,ind2+1]  # hcat stores scalar as matrix for type stability
    end
	end
	o.Repl.Yendog[ind1+1] | o.Repl.Yendog[ind2+1] ? retval : fill(retval[1], 1, cols(o.v)) # if both vars exogenous, term is same for all b; this duplication is a bit inefficient, but only arises when exog vars involved in null
end


# Workhorse for WRE CRVE sandwich filling
# With reference to notional Y = [y₁ Z], given 0-based columns index within it, ind1, and a matrix βs of all the bootstrap estimates, return all bootstrap realizations of P_X * Y[:,ind1]_g ' u\hat_1g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per cap cluster, one col per bootstrap replication
function Filling(o::StrBoottest{T}, ind1::Integer, βs::AbstractMatrix) where T
	if o.granular
		if o.Nw == 1  # create or avoid NxB matrix?
			PXYstar = iszero(ind1) ? o.Repl.PXy₁ : o.Repl.PXZ[:,ind1]
			o.Repl.Yendog[ind1+1] && (PXYstar .+= o.SstarUPX[ind1+1] * o.v)

			retval = panelsum(PXYstar .* (o.Repl.y₁ .- o.SstarUMZperp[1] * o.v), o.wt, o.infoCapData)

			for ind2 ∈ 1:o.Repl.kZ
				_β = βs[ind2,:]'
				retval .-= panelsum(PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[:,ind2] * _β .- o.SstarUMZperp[ind2+1] * (o.v .* _β) :
				                                                        o.Repl.Z[:,ind2] * _β                                       ), o.wt, infoCapData)
			end
		else  # create pieces of each N x B matrix one at a time rather than whole thing at once
			retval = Matrix{T}(undef, o.clust[1].N, cols(o.v))  # XXX preallocate this & turn Filling into Filling! ?
			for ind2 ∈ 0:o.Repl.kZ
				ind2>0 && (βv = o.v .* (_β = βs[ind2,:]'))

        if o.purerobust
          for i ∈ 1:o.clust[1].N
            PXYstar = hcat(iszero(ind1) ? o.Repl.PXy₁[i] : o.Repl.PXZ[i,ind1])
            o.Repl.Yendog[ind1+1] && (PXYstar .+= o.SstarUPX[ind1+1][i,:]'o.v)

            if iszero(ind2)
              retval[i,:]   = wtsum(o.wt, PXYstar .* (o.Repl.y₁[i] .- o.SstarUMZperp[1][i,:])'o.v)
            else
              retval[i,:] .-= wtsum(o.wt, PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[i,ind2] * _β .- o.SstarUMZperp[ind2+1][i,:]'βv :
                                                      o.Repl.Z[i,ind2] * _β))
            end
          end
        else
          for i ∈ 1:o.clust[1].N
            S = o.infoCapData[i]
            PXYstar = ind1 ? o.Repl.PXZ[S,ind1] : o.Repl.PXy₁[S]
            o.Repl.Yendog[ind1+1] && (PXYstar .+= o.SstarUPX[ind1+1][S,:] * o.v)

            if iszero(ind2)
              retval[i,:]   = wtsum(o.wt, PXYstar .* (o.Repl.y₁[S] .- o.SstarUMZperp[1][S,:] * o.v))
            else
              retval[i,:] .-= wtsum(o.wt, PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[S,ind2] * _β .- o.SstarUMZperp[ind2+1][S,:] * βv :
                                                                              o.Repl.Z[S,ind2] * _β                                       ))
            end
          end
        end
			end
		end
	else  # coarse error clustering
		for ind2 ∈ 0:o.Repl.kZ
			βv = iszero(ind2) ? o.v : o.v .* (_β = -βs[ind2,:]')

			o.Repl.Yendog[ind1+1] && (T₁ = o.Repl.ScapYX[ind2+1] * o.SstarUXinvXX[ind1+1])  #  S_cap (Y_(parj).*X_par ) [S_* (U ̈_(pari).*X_par (X_par^' X_par )^(-1) )]^'

			if o.Repl.Yendog[ind2+1]  # add CT_(cap,*) (P_(X_par ) Y_(pari).*U ̈_(parj) )
				if o.NClustVar == o.NBootClustVar && iszero(o.subcluster)  # simple case of one clustering: full crosstab is diagonal
					tmp = ind1>0 ? o.Repl.XZ[:,ind1] : o.Repl.Xy₁par
          if iszero(length(T₁))
						T₁                = o.SstarUXinvXX[ind2+1]'tmp  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					else
            T₁[diagind(T₁)] .+= o.SstarUXinvXX[ind2+1]'tmp
          end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JNcapNstar)
					for i ∈ 1:o.Nstar
						T₁[o.IDCTCapstar[i], i] .+= o.SCTcapuXinvXX[ind2+1,i] * (ind1>0 ? o.Repl.XZ[:,ind1] : o.Repl.Xy₁par)
          end
				end
				length(o.Repl.Zperp) > 0 && (T₁ .-= o.Repl.ScapPXYZperp[ind1+1] * o.SstarUZperpinvZperpZperp[ind2+1])
        o.NFE                > 0 && (T₁ .-= o.Repl.CT_FEcapPY[ind1+1]'o.CTFEU[ind2+1])
			end

			if ind2>0
        if isone(cols(T₁))
          retval .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β .+ T₁ .* βv   # - x*β components
        else
          retval .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β
          matmulplus!(retval, T₁, βv)
        end
      else
        if isone(cols(T₁))
          retval = o.Repl.FillingT₀[ind1+1,1] .+ T₁ .* βv     # y component
        else
          retval = T₁ * o.v; retval .+= o.Repl.FillingT₀[ind1+1,1]
        end
      end

			if o.Repl.Yendog[ind1+1] && o.Repl.Yendog[ind2+1]
				for i ∈ 1:o.clust[1].N
					S = o.infoCapData[i]
					colquadformminus!(retval, i, cross(o.SstarUPX[ind1+1][S,:], o.haswt ? o.wt[S] : I, o.SstarUMZperp[ind2+1][S,:]), o.v, βv)
				end
      end
		end
	end
	retval
end

function PrepWRE!(o::StrBoottest)
	Estimate!(o.DGP, o.null ? [o.r₁ ; o.r] : o.r₁)
  MakeResiduals!(o.DGP)
	o.Ü₂par = o.DGP.Ü₂ * o.Repl.RparY  # XXX XB(o.DGP.Ü₂, o.Repl.RparY)

	for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		uwt = vHadw(i>0 ? o.Ü₂par[:,i] : o.DGP.u⃛₁, o.wt)

    # S_star(u .* X), S_star(u .* Zperp) for residuals u for each endog var; store transposed
    o.SstarUX[i+1]      = panelsum(o.Repl.X₁, o.Repl.X₂, uwt, o.infoBootData)'
    o.SstarUXinvXX[i+1] = o.Repl.invXX * o.SstarUX[i+1]

    if o.LIML || o.bootstrapt || !isone(o.κ)
      o.SstarUZperp[i+1]              = panelsum(o.Repl.Zperp, uwt, o.infoBootData)'
      o.SstarUZperpinvZperpZperp[i+1] = o.Repl.invZperpZperp * o.SstarUZperp[i+1]
      o.NFE>0 && (o.CTFEU[i+1] = crosstabFE(o, uwt, o.infoBootData))
    end

    if o.LIML || !o.robust || !isone(o.κ)
      o.SstaruY[i+1] = panelsum(o.Repl.y₁par, o.Repl.Z, uwt, o.infoBootData)
      for j ∈ 0:i
        o.SstarUU[i+1,j+1] = panelsum(j>0 ? o.Ü₂par[:,j] : o.DGP.u⃛₁, uwt, o.infoBootData)
      end
    end

    if o.robust && o.bootstrapt
      if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
        for g ∈ 1:o.Nstar
        	o.SCTcapuXinvXX[i+1,g] = panelsum(o.Repl.XinvXX, uwt, o.infoCTCapstar[g])
        end
      end

      o.SstarUPX[i+1]     = o.Repl.XinvXX * o.SstarUX[i+1]
      o.SstarUMZperp[i+1] = o.Repl.Zperp  * o.SstarUZperpinvZperpZperp[i+1]
      if o.Nobs == o.Nstar  # subtract "crosstab" of observation by cap-group of u
        o.SstarUMZperp[i+1][diagind(o.SstarUMZperp[i+1])] .-= i>0 ? o.Ü₂par[:,i] : o.DGP.u⃛₁  # case: bootstrapping by observation
      else
        for g ∈ 1:o.Nobs
          o.SstarUMZperp[i+1][g,o.IDBootData[g]] -= iszero(i) ? o.DGP.u⃛₁[g] : o.Ü₂par[g,i]
        end
      end
      o.NFE>0 &&
        (o.SstarUMZperp[i+1] .+= (o.invFEwt .* o.CTFEU[i+1])[o._FEID,:])  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
    end
	end
end

function MakeWREStats!(o::StrBoottest{T}, w::Integer) where T
	if isone(o.Repl.kZ)  # optimized code for 1 coefficient in bootstrap regression
		if o.LIML
      YY₁₁   = HessianFixedκ(o, [0], 0, zero(T))  # κ=0 => Y*MZperp*Y
      YY₁₂   = HessianFixedκ(o, [0], 1, zero(T))
      YY₂₂   = HessianFixedκ(o, [1], 1, zero(T))
      YPXY₁₁ = HessianFixedκ(o, [0], 0, one(T) )  # κ=1 => Y*PXpar*Y
      YPXY₁₂ = HessianFixedκ(o, [0], 1, one(T) )
      YPXY₂₂ = HessianFixedκ(o, [1], 1, one(T) )
      YY₁₂YPXY₁₂ = YY₁₂ .* YPXY₁₂
      x₁₁ = YY₂₂ .* YPXY₁₁ .- YY₁₂YPXY₁₂      # elements of YYstar^-1 * YPXYstar up to factor of det(YYstar)
      x₁₂ = YY₂₂ .* YPXY₁₂ .- YY₁₂ .* YPXY₂₂
      x₂₁ = YY₁₁ .* YPXY₁₂ .- YY₁₂ .* YPXY₁₁
      x₂₂ = YY₁₁ .* YPXY₂₂ .- YY₁₂YPXY₁₂
      κs = (x₁₁ .+ x₂₂)/2; κs = 1 ./ (1 .- (κs - sqrt.(κs.^2 .- x₁₁.*x₂₂ .+ x₁₂.*x₂₁)) ./ (YY₁₁ .* YY₂₂ .- YY₁₂ .* YY₁₂))  # solve quadratic equation for smaller eignenvalue; last term is det(YYstar)
      !iszero(o.Fuller) && (κs .-= o.Fuller / (o._Nobs - o.kX))
      o.As = κs .* (YPXY₂₂ .- YY₂₂) .+ YY₂₂
      βs = (κs .* (YPXY₁₂ .- YY₁₂) .+ YY₁₂) ./ o.As
		else
      o.As = HessianFixedκ(o, [1], 1, o.κ)
			βs = HessianFixedκ(o, [1], 0, o.κ) ./ o.As
    end

    if o.null
 			o.numerw = βs .+ (o.Repl.Rt₁ - o.r) / o.Repl.RRpar
		else
			o.numerw = βs .- o.DGP.β₀
			isone(w) && (o.numerw[1] = βs[1] + (o.Repl.Rt₁ - o.r) / o.Repl.RRpar)
		end

    @storeWtGrpResults!(o.numer, o.numerw)

		if o.bootstrapt
			if o.robust
				Jcaps = Filling(o, 1, βs) ./ o.As
				for c ∈ 1:o.NErrClustCombs  # sum sandwich over error clusterings
					length(o.clust[c].order)>0 && (Jcaps = Jcaps[o.clust[c].order,:])
					@clustAccum!(denom, c, coldot(panelsum(Jcaps, o.clust[c].info)))
        end
			else
        denom = (HessianFixedκ(o,[0],0,zero(T)) .- 2 .* βs .* HessianFixedκ(o, [0], 1, zero(T)) .+ βs.^2 .* HessianFixedκ(o, [1], 1, zero(T))) ./ o._Nobs ./ o.As  # classical error variance
      end
      @storeWtGrpResults!(o.Dist, view(o.sqrt ? o.numerw ./ sqrt.(denom) : o.numerw .^ 2 ./ denom, 1, :))
      denom *= o.Repl.RRpar[1]^2
		end

	else  # WRE bootstrap for more than 1 coefficeint in bootstrap regression

		βs = zeros(T, o.Repl.kZ, cols(o.v))
		A = Vector{Matrix{T}}(undef, cols(o.v))

		if o.LIML
      o.YYstar   = [HessianFixedκ(o, collect(0:i), i, zero(T)) for i ∈ 0:o.Repl.kZ] # κ=0 => Y*MZperp*Y
      o.YPXYstar = [HessianFixedκ(o, collect(0:i), i,  one(T)) for i ∈ 0:o.Repl.kZ] # κ=1 => Y*PXpar*Y

			for b ∈ 1:cols(o.v)
				for i ∈ 0:o.Repl.kZ
					o.YYstar_b[1:i+1,i+1]   = o.YYstar[i+1][:,b]  # fill uppper triangles, which is all that invsym() looks at
					o.YPXYstar_b[1:i+1,i+1] = o.YPXYstar[i+1][:,b]
				end
				o.κ = 1/(1 - eigvals(invsym(o.YYstar_b) * Symmetric(o.YPXYstar_b))[1])
				!iszero(o.Fuller) && (o.κ -= o.Fuller / (o._Nobs - o.kX))
				βs[:,b] = (A[b] = invsym(o.κ*o.YPXYstar_b[2:end,2:end] + (1-o.κ)*o.YYstar_b[2:end,2:end])) * (o.κ*o.YPXYstar_b[1,2:end]' + (1-o.κ)*o.YYstar_b[1,2:end]')
			end
		else
			δnumer = HessianFixedκ(o, collect(1:o.Repl.kZ), 0, o.κ)

			for i ∈ 1:o.Repl.kZ
				o.δdenom[i] = HessianFixedκ(o, collect(1:i), i, o.κ)
      end

			for b ∈ axes(o.v,2)
				for i ∈ 1:o.Repl.kZ
          o.δdenom_b[1:i,i] = o.δdenom[i][:,b] # fill uppper triangle, which is all that invsym() looks at
        end
				βs[:,b]  = (A[b] = invsym(o.δdenom_b)) * δnumer[:,b]
			end
		end

		if o.bootstrapt
      if o.robust
        Zyg = Vector{Matrix{T}}(undef, o.Repl.kZ)
        for i ∈ 1:o.Repl.kZ
          Zyg[i] = Filling(o, i, βs)  # XXX concatenate into 3-d array?
        end
      else
        o.YYstar = [HessianFixedκ(o, collect(i:o.Repl.kZ), i, zero(T)) for i ∈ 0:o.Repl.kZ] # κ=0 => Y*MZperp*Y
      end
    end

		numer_b = Vector{T}(undef,rows(o.RRpar))  # XXX move to Init()
    for b ∈ cols(o.v):-1:1
			if o.null || w==1 && b==1
        numer_b .= o.Repl.RRpar * βs[:,b] + o.Repl.Rt₁ - o.r
      else
        numer_b .= o.Repl.RRpar * (βs[:,b] - o.DGP.β₀)
      end

			if o.bootstrapt
				if o.robust  # Compute denominator for this WRE test stat
          for i ∈ 1:o.Repl.kZ  # XXX replace with 3-D array
            _Jcap[:,i] = Zyg[i][:,b]
          end
          Jcap = _Jcap * (A[b] * o.Repl.RRpar')

					for c ∈ 1:o.NErrClustCombs
						(!isone(o.NClustVar) && length(o.clust[c].order)>0) && (Jcap = Jcap[o.clust[c].order,:])
						J_b = panelsum(Jcap, o.clust[c].info)
            @clustAccum!(denom, c, J_b'J_b)
          end
				else  # non-robust
          for i ∈ 0:o.Repl.kZ
            YYstar_b[i+1,i+1:o.Repl.kZ+1] = YYstar[i+1][b,:]  # fill upper triangle
          end
          denom = (o.Repl.RRpar * A[b] * o.Repl.RRpar') * [-one(T) ; βs[:,b]]'Symmetric(YYstar_b) * [-one(T) ; βs[:,b]] / o._Nobs  # 2nd half is sig2 of errors
        end
				if o.sqrt
          o.Dist[b+first(o.WeightGrp[w])-1] .= numer_b ./ sqrt.(denom)
        else
          o.Dist[b+first(o.WeightGrp[w])-1] .= numer_b'invsym(denom)*numer_b  # hand-code for 2-dimensional?
        end
			end
			o.numer[:,b+first(o.WeightGrp[w])-1] = numer_b  # slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
		end
	end
	w==1 && o.bootstrapt && (o.statDenom = denom)  # original-sample denominator
end


# Construct stuff that depends linearly or quadratically on r, possibly by interpolation
function MakeInterpolables!(o::StrBoottest{T} where T)
  if o.interpolable
    if iszero(length(o.anchor))  # first call? save current r as permanent anchor for interpolation
      o.anchor = o.r
      _MakeInterpolables!(o, o.anchor)
      o.numer₀ = o.numer
      o.interpolate_u && (o.ü₀ = o.ü)  # XXX would need to copy() if o.ü were not allocated from scratch on each construction
      o.robust && (o.Jcd₀ = deepcopy(o.Jcd))
      return
    end

    if iszero(length(o.poles))  # second call: from anchor make set of orthogonal poles, which equal anchor except in one dimension
      o.poles = o.r - o.anchor
      o.robust && (o.denom₀ = deepcopy(o.denom))  # grab quadratic denominator from *previous* (1st) evaluation
      newPole = trues(o.q, 1)  # all poles new
    else  # been here at least twice? interpolate unless current r stretches range > 2X in some dimension(s)
      newPole = abs.(o.r - o.anchor) .> 2 * abs.(o.poles)
    end

    if any(newPole)  # prep interpolation
      for h₁ ∈ 1:o.q
        if newPole[h₁]
        	o.poles[h₁] = o.r[h₁] - o.anchor[h₁]
          thisr = copy(o.anchor); thisr[h₁] = o.r[h₁]  # if q>1 this creates anchor points that are not graphed, an inefficiency. But simpler to make the deviations from 1st point orthogonal
          _MakeInterpolables!(o, thisr)  # calculate linear stuff at new anchor

          o.∂numer∂r[h₁] = (o.numer - o.numer₀) / o.poles[h₁]
          o.interpolate_u && (o.∂u∂r[h₁] = (o.ü - o.ü₀) / o.poles[h₁])
          if o.robust  # df > 1 for an ARubin test with >1 instruments.
            for d₁ ∈ 1:o.df
              for c ∈ 1:o.NErrClustCombs
                o.∂Jcd∂r[h₁][c,d₁] = (o.Jcd[c,d₁] - o.Jcd₀[c,d₁]) / o.poles[h₁]
                for d₂ ∈ 1:d₁
                              tmp   = coldot(o.Jcd₀[c,d₁], o.∂Jcd∂r[h₁][c,d₂])
                  d₁ ≠ d₂ && (tmp .+= coldot(o.Jcd₀[c,d₂], o.∂Jcd∂r[h₁][c,d₁]))  # for diagonal items, faster to just double after the c loop
                  @clustAccum!(o.∂denom∂r[h₁][d₁,d₂], c, tmp)
                end
              end
              o.∂denom∂r[h₁][d₁,d₁] .*= 2  # double diagonal terms
            end
          end
        end
      end
      if o.robust  # quadratic interaction terms
        for h₁ ∈ 1:o.q, h₂ ∈ 1:h₁
          if newPole[h₁] || newPole[h₂]
            for d₁ ∈ 1:o.df, d₂ ∈ 1:d₁, c ∈ 1:o.NErrClustCombs
              @clustAccum!(o.∂²denom∂r²[h₁,h₂][d₁,d₂], c, coldot(o.∂Jcd∂r[h₁][c,d₁], o.∂Jcd∂r[h₂][c,d₂]))
            end
          end
        end
      end
      Δ = o.poles
      o.interpolating = true

      if o.q==2  # in this case we haven't yet actually computed interpolables at *pr, so interpolate them
        o.numerw .= o.numer₀ .+ o.∂numer∂r[1] .* Δ[1] .+ o.∂numer∂r[2] .* Δ[2]
        if o.interpolate_u
          o.ü .= o.ü₀ .+ o.∂u∂r[1] .* Δ[1] .+ o.∂u∂r[2] .* Δ[2]
        end
      end

    else  # routine linear interpolation if the anchors not moved
      Δ = o.r - o.anchor
      o.numerw = o.numer₀ + o.∂numer∂r[1] * Δ[1]
      o.q > 1 && (o.numerw .+= o.∂numer∂r[2] * Δ[2])
      if o.interpolate_u
        o.ü = o.ü₀ + o.∂u∂r * Δ[1]
        o.q > 1 && (o.ü .+= o.∂u∂r[2] * Δ[2])
      end
    end

    if o.robust  # even if an anchor was just moved, and linear components just computed from scratch, do the quadratic interpolation now, from the updated linear factors
      if isone(o.q)
        for d₁ ∈ 1:o.df, d₂ ∈ 1:d₁
          o.denom[d₁,d₂] = o.denom₀[d₁,d₂] .+ o.∂denom∂r[d₁,d₂][1,1] .* Δ .+ o.∂²denom∂r²[d₁,d₂][1,1] .* Δ.^2
        end
      else  # q==2
        for d₁ ∈ 1:o.df, d₂ ∈ 1:d₁
          o.denom[d₁,d₂] = o.denom₀[d₁,d₂] +
                            o.∂denom∂r[1][d₁,d₂] * Δ[1] +
                            o.∂denom∂r[2][d₁,d₂] * Δ[2] +
                            o.∂²denom∂r²[1,1][d₁,d₂] * (Δ[1] ^ 2) +
                            o.∂²denom∂r²[2,1][d₁,d₂] * (Δ[1] * Δ[2]) +
                            o.∂²denom∂r²[2,2][d₁,d₂] * (Δ[2] ^ 2)
        end
      end
    end
  else  # non-interpolable cases
    _MakeInterpolables!(o, o.r)
  end
end

# Construct stuff that depends linearly or quadratically on r and doesn't depend on v. No interpolation.
function _MakeInterpolables!(o::StrBoottest{T}, thisr::AbstractVector) where T
  if o.ML
		o.uXAR = o.Sc * (o.AR = o.A * o.R')
	else
    if o.ARubin
      Estimate!(o.DGP, thisr)
    elseif iszero(o.κ)  # regular OLS
    	Estimate!(o.DGP, o.null ? [o.r₁ ; thisr] : o.r₁)
    elseif o.null  # in score bootstrap for IV/GMM, if imposing null, then DGP constraints, κ, Hessian, etc. do vary with r and must be set now
      Estimate!(o.DGP, [o.r₁ ; thisr])
      InitTestDenoms!(o.DGP)
    end

    MakeResiduals!(o.DGP)
    o.ü = o.DGP.ü₁

		(o.scoreBS || (o.robust && Int(o.granular) < o.NErrClustCombs)) &&
      (o.uXAR = o.DGP.ü₁ .* o.M.XAR)
  end

  o.SuwtXA = o.scoreBS ?
	             o.B>0 ?
		             o.NClustVar ?
                   panelsum(o.uXAR, o.wt, o.infoBootData) :
					         vHadw(o.uXAR, o.wt)                    :
				         wtsum(o.wt, o.uXAR)                      :
              o.DGP.A * panelsum(o.X₁, o.X₂, vHadw(o.ü, o.wt), o.infoBootData)'  # same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X₂ empty except in Anderson-Rubin

  if o.robust && o.bootstrapt && Int(o.granular) < o.NErrClustCombs
    ustarXAR = panelsum(o.uXAR, o.wt, o.infoAllData)  # collapse data to all-boot && error-cluster-var intersections. If no collapsing needed, panelsum() will still fold in any weights
    if o.B>0
      Kd = Vector{Matrix{T}}(undef, o.df)
      if o.scoreBS
        for d ∈ 1:o.df
          Kd[d] = copy(JNcapNstar)  # inefficient, but not optimizing for the score bootstrap
        end
      else
        for d ∈ 1:o.df
          Kd[d] = (panelsum(o.X₁, o.X₂, vHadw(o.DGP.XAR[:,d], o.wt), o.infoCapData)) * o.SuwtXA  # final term in (64), for c=intersection of all error clusters
        end
      end

      o.NFE>0 && !o.FEboot && (o.CT_WE = crosstabFE(o, vHadw(o.ü, o.wt), o.infoBootData))

			for d ∈ 1:o.df  # subtract crosstab of u.*XAR wrt bootstrapping cluster combo and all-cluster-var intersections
				crosstabCapstarMinus!(o, Kd[d], ustarXAR[:,d])
        if o.NFE>0 && !o.FEboot
          Kd[d] .+= o.M.CT_XAR[d]' * (o.invFEwt .* o.CT_WE)  # middle term of (64)
        end
        o.scoreBS && (Kd[d] .-= o.ClustShare * colsum(Kd[d])) # recenter
			end

      for c ∈ 1+Int(o.granular):o.NErrClustCombs  # XXX pre-compute common iterators
        length(o.clust[c].order)>0 &&
          (Kd = [Kd[d][o.clust[c].order,:] for d ∈ 1:o.df])  # XXX replace with 3-D array?
        for d ∈ 1:o.df
          o.Kcd[c,d] = panelsum(Kd[d], o.clust[c].info)
        end
      end
    else  # B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
      o.scoreBS &&
        (ustarXAR .-= o.ClustShare * colsum(ustarXAR))  # recenter if OLS
      for c ∈ 1:o.NErrClustCombs
        length(o.clust[c].order)>0 &&
          (ustarXAR = ustarXAR[o.clust[c].order,:])
        tmp = panelsum(ustarXAR, o.clust[c].info)
        for d ∈ 1:o.df
          o.Kcd[c,d] = tmp[:,d][:,:]  # XXX a little inefficient -- fix when making Kcd's elements 3-D?
        end
      end
    end
  end
  MakeNumerAndJ!(o, 1, thisr)  # compute J = κ * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
end

# compute stuff depending linearly on v, needed to prep for interpolation
function MakeNumerAndJ!(o::StrBoottest{T}, w::Integer, r::AbstractVector=Vector{T}(undef,0)) where T  # called to *prepare* interpolation, or when w>1, in which case there is no interpolation
  o.numerw = o.scoreBS ?
               (o.B>0 ?
                 o.SuwtXA'o.v :
                 o.SuwtXA * o.u_sd    ) :
               (!o.robust || o.granular || o.purerobust ?
                  o.R * (o.βdev = o.SuwtXA * o.v) :
                 (o.R * o.SuwtXA) * o.v)

  if isone(w)
  	if o.ARubin
      o.numerw[:,1] = o.u_sd * o.DGP.Rpar * o.DGP.β[o.kX₁+1:end]  # coefficients on excluded instruments in ARubin OLS
    elseif !o.null
      o.numerw[:,1] = o.u_sd * (o.R * (o.ML ? o.β : o.M.Rpar * o.M.β) - r)  # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this.
    end
  end

  @storeWtGrpResults!(o.numer, o.numerw)

	@views if o.B>0 && o.robust && o.bootstrapt
    if o.granular || o.purerobust  # optimized treatment when bootstrapping by many/small groups
      if o.purerobust
        o.ustar = partialFE!(o, o.ü .* o.v) - X₁₂B(o.X₁, o.X₂, o.βdev)
      else  # clusters small but not all singletons
        if o.NFE>0 && !o.FEboot
          o.ustar = partialFE!(o, o.ü .* o.v[o.IDBootData,:])
          for d ∈ 1:o.df
            o.Jcd[1,d] = panelsum(o.ustar, o.M.WXAR[:,d], o.infoCapData)                                            - panelsum(o.X₁, o.X₂, o.M.WXAR[:,d], o.infoCapData) * o.βdev
          end
        else
          for d ∈ 1:o.df
            o.Jcd[1,d] = panelsum( panelsum(o.ü, o.M.WXAR[:,d], o.infoAllData) .* o.v[o.IDBootAll,:], o.infoErrAll) - panelsum(o.X₁, o.X₂, o.M.WXAR[:,d], o.infoCapData) * o.βdev
          end
        end
      end
    end
	  for c ∈ Int(o.granular)+1:o.NErrClustCombs, d ∈ eachindex(axes(o.Jcd, 2), axes(o.Kcd, 2))
      o.Jcd[c,d] = o.Kcd[c,d] * o.v
    end
  end
end

function MakeNonWREStats!(o::StrBoottest{T}, w::Integer) where T
  w > 1 && MakeNumerAndJ!(o, w)
  !o.bootstrapt && return

	if o.robust
    if !o.interpolating  # these quadratic computation needed to *prepare* for interpolation but are superseded by interpolation once it is going
      o.purerobust && (ustar2 = o.ustar .^ 2)
      for i ∈ 1:o.df, j ∈ 1:i
        o.purerobust &&
          (o.denom[i,j] = cross(o.M.WXAR[:,i], o.M.WXAR[:,j], ustar2) * (o.clust[1].even ? o.clust[1].multiplier : -o.clust[1].multiplier))
        for c ∈ Int(o.purerobust)+1:o.NErrClustCombs
          @clustAccum!(o.denom[i,j], c, j==i ? coldot(o.Jcd[c,i]) : coldot(o.Jcd[c,i],o.Jcd[c,j]))
        end
      end
    end
    if isone(o.df)
    	@storeWtGrpResults!(o.Dist, vec(o.numerw ./ sqrtNaN.(o.denom[1,1])))
      isone(w) &&
        (o.statDenom = hcat(o.denom[1,1][1]))  # original-sample denominator
		elseif o.df==2  # hand-code 2D numer'inv(denom)*numer
    	t1 = o.numerw[1,:]'; t2 = o.numerw[2,:]'; t12 = t1.*t2
			@storeWtGrpResults!(o.Dist, vec((t1.^2 .* o.denom[2,2] .- 2 .* t12 .* o.denom[2,1] .+ t2.^2 .* o.denom[1,1]) ./ (o.denom[1,1].*o.denom[2,2] .- o.denom[2,1].^2)))
      isone(w) &&
        (o.statDenom = [o.denom[1,1][1] o.denom[2,1][1] ; o.denom[2,1][1] o.denom[2,2][1]])  # original-sample denominator
    else  # build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = Matrix{T}(undef, o.df, o.df)
			for k ∈ 1:cols(o.v)  # XXX probably can simplify
				for i ∈ 1:o.df, j ∈ 1:i
					tmp[j,i] = o.denom[i,j][k]  # fill upper triangle, which is all invsym() looks at
        end
				numer_l = o.numerw[:,k]
				o.Dist[k+first(o.WeightGrp[w])-1] = numer_l'invsym(tmp)*numer_l  # in degenerate cases, cross() would turn cross(.,.) into 0
			end
			isone(w) && (o.statDenom = tmp)  # original-sample denominator
		end
	else  # non-robust
		AR = o.ML ? o.AR : o.M.AR
		if isone(o.df)  # optimize for one null constraint
			o.denom[1,1] = o.R * AR
			if !o.ML
        o.ustar = o.B>0 ? o.v .* o.ü : o.ü
				if o.scoreBS
          if o.haswt  # Center variance if interpolated
            o.ustar .-= o.ClustShare'o.ustar
          else
            o.ustar .-= colsum(o.ustar) * o.ClustShare  # Center variance if interpolated
          end
        else
          o.ustar -= X₁₂B(o.X₁, o.X₂, o.βdev)  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
        end
        if o.haswt
          o.denom[1,1] .*= o.wt'(o.ustar .^ 2)
        else
          o.denom[1,1] .*= coldot(o.ustar)
        end
      end

      @storeWtGrpResults!(o.Dist, vec(o.numerw ./ sqrtNaN.(o.denom[1,1])))
			isone(w) && (o.statDenom = o.denom[1,1])  # original-sample denominator
		else
			o.denom[1,1] = o.R * AR
			if o.ML
				for k ∈ 1:cols(o.v)
					numer_l = o.numerw[:,k]
					o.Dist[k+first(o.WeightGrp[w])-1] = numer_l'invsym(o.denom[1,1])*numer_l
				end
				isone(w) && (o.statDenom = o.denom[1,1])  # original-sample denominator
			else
        invdenom = invsym(o.denom[1,1])
				for k ∈ 1:cols(o.v)
					numer_l = o.numerw[:,k]
					o.Dist[k+first(o.WeightGrp[w])-1] = o.numer_l'invdenom*numer_l
          o.ustar = o.B>0 ? o.v[:,k] .* o.ü : o.ü
					if o.scoreBS
            if o.haswt  # Center variance if interpolated
              o.ustar .-= o.wt'o.ustar * o.ClustShare
            else
              o.ustar .-= colsum(o.ustar) * o.ClustShare  # Center variance if interpolated
            end
					else
            o.ustar .-= X₁₂B(o.X₁, o.X₂, o.βdev[:,k])  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
          end
					o.Dist[k+first(o.WeightGrp[w])-1] ./= (tmp = symcross(o.ustar, o.wt))
				end
				isone(w) && (o.statDenom = o.denom[1,1] * tmp)  # original-sample denominator
			end
		end
	end
end


# like Mata panelsetup() but can group on multiple columns, like sort(), and faster. But doesn't take minobs, maxobs arguments.
function panelsetup(X::AbstractArray{S} where S, colinds::Vector{T} where T<:Integer)
	N = rows(X)
	info = Vector{UnitRange{Int64}}(undef, N)
	lo = p = 1
	id = @view X[1, colinds]
	@inbounds for hi ∈ 2:N
		if (tmp=X[hi,colinds]) ≠ id
			info[p] = lo:hi-1
      lo = hi
			p += 1
			id = tmp
		end
	end
  info[p] = lo:N
	resize!(info, p)
  info
end
# Like above but also return standardized ID variable, starting from 1
function panelsetupID(X::AbstractArray{S} where S, colinds::Vector{T} where T<:Integer)
	N = rows(X)
	info = Vector{UnitRange{Int64}}(undef, N)
  ID = ones(Int64, N)
	lo = p = 1
	id = @view X[1, colinds]
	@inbounds for hi ∈ 2:N
		if (tmp=X[hi,colinds]) ≠ id
			info[p] = lo:hi-1
      lo = hi
			p += 1
			id = tmp
		end
		ID[hi] = p
	end
  info[p] = lo:N
	resize!(info, p)
  info, ID
end

# Do panelsum() except that a single missing value in X doesn't make all results missing and
# efficiently handles case when all groups have one row.
function panelsum(X::AbstractMatrix, info::Vector{UnitRange{T}} where T<:Integer)
  iszero(cols(X)) && return X[1:rows(info),:]
  (iszero(length(info)) || rows(info)==rows(X)) && return X

  retval = zeros(promote_type(Float32, eltype(X)), size(info,1), size(X,2))
  for g in eachindex(info)
    f, l = first(info[g]), last(info[g])
    @turbo for j ∈ axes(X,2)
      for i ∈ f:l
        retval[g,j] += X[i,j]
      end
    end
  end
  retval
end
function panelsum(X::AbstractVector, info::Vector{UnitRange{T}} where T<:Integer)
  (iszero(length(info)) || rows(info)==length(X)) && return X

  retval = zeros(promote_type(Float32, eltype(X)), size(info,1))
  for g in eachindex(info, retval)
    f, l = first(info[g]), last(info[g])
    _retval = zero(eltype(retval))
    @turbo for i ∈ f:l
      _retval += X[i]
    end
    retval[g] = _retval
  end
  retval
end

@inline panelsum(X::AbstractArray, wt::UniformScaling, info::Vector{UnitRange{T}} where T<:Integer) = panelsum(X, info)

function panelsum(X::AbstractArray, wt::AbstractVector, info::Vector{UnitRange{T}} where T<:Integer)
  iszero(cols(X)) && return X[1:rows(info),:]
  (iszero(length(info)) || rows(info)==rows(X)) && return X .* wt

  if ndims(X)==1  # X maybe a 1D view, which is not a subtype of Vector and so cannot be trapped by such a function signature
    retval = zeros(promote_type(Float32, eltype(X)), size(info,1))
    for g in eachindex(info)
      f, l = first(info[g]), last(info[g])
      @inbounds for i ∈ f:l  # @turbo doesn't seem to work with views
        retval[g] += X[i] * wt[i]
      end
    end
  else
    retval = zeros(eltype(X), size(info,1), size(X,2))
    for g in eachindex(info)
      f, l = first(info[g]), last(info[g])
      @turbo for j ∈ axes(X,2)
        for i ∈ f:l
          retval[g,j] += X[i,j] * wt[i]
        end
      end
    end
  end
  retval
end

@inline panelsum(X₁::AbstractArray, X₂::AbstractArray, wt::AbstractVector, info::Vector{UnitRange{T}} where T<:Integer) =
  iszero(length(X₁)) ? panelsum(X₂,wt,info) :
  iszero(length(X₂)) ? panelsum(X₁,wt,info) :
#                       hcat(panelsum(X₁,wt,info), panelsum(X₂,wt,info))
                       ApplyArray(hcat, panelsum(X₁,wt,info), panelsum(X₂,wt,info))

@inline panelsum(X₁::AbstractArray, X₂::AbstractArray, wt::UniformScaling, info::Vector{UnitRange{T}} where T<:Integer) =
  iszero(length(X₁)) ? panelsum(X₂,info) :
  iszero(length(X₂)) ? panelsum(X₁,info) :
#                      hcat(panelsum(X₁,info), panelsum(X₂,info))
                      ApplyArray(hcat, panelsum(X₁,info), panelsum(X₂,info))

# Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
# except use provided lo and hi values for 0 and 1
count_binary(N::Integer, lo::Number, hi::Number) = N<=1 ? [lo  hi] :
                                                          (tmp = count_binary(N-1, lo, hi);
                                                           [fill(lo, 1, cols(tmp)) fill(hi, 1, cols(tmp)) ;
                                                                             tmp                    tmp     ]
                                                          )

# cross-tab sum of a column vector w.r.t. given panel info and fixed-effect var
# one row per FE, one col per other grouping
function crosstabFE(o::StrBoottest{T}, v::AbstractVector, info::Vector{UnitRange{Int64}}) where T
	retval = zeros(T, o.NFE, rows(info))
  if length(info)>0
    for i ∈ 1:rows(info)
      FEIDi = @view o._FEID[info[i],:]
      vi    = @view       v[info[i],:]
      @inbounds for j in eachindex(vi, FEIDi)
        retval[FEIDi[j],i] += vi[j]
      end
    end
  else  # "robust" case, no clustering
    @inbounds for i in eachindex(v,o._FEID)
      retval[o._FEID[i],i] = v[i]
    end
  end
	return retval
end

# subtract crosstab of v wrt bootstrapping cluster and all-cluster-var intersections from M
# M should have one row for each all-cluster-var (including bootstrap cluster) intersection and one col for each bootstrap cluster
# *** v needs to have been panelsum'd with infoAllData
# XXX this should be reducible to a single line of the form M[inds] .-= v where inds is pre-computed
function crosstabCapstarMinus!(o::StrBoottest, M::AbstractMatrix, v::AbstractVector)
  isa(M, ApplyArray) && (M = M[:,:])  # can't edit lazy array
	if o.subcluster>0  # crosstab c,c* is wide
		for i ∈ eachindex(o.infoErrAll, axes(M,1))
			M[i, o.infoErrAll[i]] .-= v[o.infoErrAll[i]]
		end
	elseif o.NClustVar == o.NBootClustVar  # crosstab c,c* is square
		M[diagind(M)] .-= v
	else  # crosstab c,c* is tall
		for i ∈ eachindex(o.clust[o.BootClust].info, axes(M,2))
			tmp = o.clust[o.BootClust].info[i]
			M[tmp,i] .-= v[tmp]
		end
  end
end

# given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r
# used with optimize() to construct confidence intervals
# performs no error checking
function r_to_p(o::StrBoottest, r::AbstractVector)
	o.r = r
	setdirty!(o, true, noinitialize=true) # set dirty = 1, but leave initialized=0, which we want when only changing r
	return getpadj(o)
end


# Chandrupatla 1997, "A new hybrid quadratic-bisection algorithm for finding the zero of a nonlinear function without using derivatives"
# x₁, x₂ must bracket the true value, with f₁=f(x₁) and f₂=f(x₂)
function search(o::StrBoottest, α::T, f₁::T, x₁::T, f₂::T, x₂::T) where T<:Real
	t = half = T(.5)
	while true
		x = x₁ + t * (x₂ - x₁)
    fx = r_to_p(o, [x])
		((fx>f₁) == (fx>f₂)) && return x  # violation of monotonicity because of precision problems? That's as good as it gets.

    if (fx<α) == (f₁<α)
			x₃, x₁, f₃, f₁ = x₁, x, f₁, fx
		else
			x₃, x₂, x₁, f₃, f₂, f₁ = x₂, x₁, x, f₂, f₁, fx
		end

		((o.B>0 && abs(fx - α) < (1+Int(o.ptype==equaltail)) / o.BFeas * 1.000001) || ≈(x₂, x₁, rtol=o.rtol)) &&
			return abs(f₁ - α) < abs(f₂ - α) ? x₁ : x₂

		ϕ₁ = (f₁ - f₂) / (f₃ - f₂)
		ϕ₁² = ϕ₁^2
		xi1 = (x₁ - x₂) / (x₃ - x₂)
    t = ϕ₁² > xi1 || xi1 > 2 * ϕ₁ - ϕ₁² ?
          half :
          clamp(((f₃ - α) / (f₁ - f₂) + (x₃ - x₁) / ((x₂ - x₁) * (f₃ - f₁)) * (f₂ - α)) * (f₁ - α) / (f₃ - f₂), 0.000001, 0.999999)
	end
end

# derive wild bootstrap-based CI, for case of linear model with one-degree null imposed.
function plot(o::StrBoottest{T}) where T
	_quietly = o.quietly; _r = copy(o.r)
	setquietly!(o, true)
  α = one(T) - o.level
  o.gridpoints[ismissing.(o.gridpoints)] .= 25

  boottest!(o)
  if !o.ARubin
    halfwidth = T.(-1.5 * quantile(Normal(), α/2)) .* sqrtNaN.(diag(getV(o)))
    o.confpeak = getb(o) + o.r
  else
    halfwidth = abs.(o.confpeak) * T.(quantile(Normal(), getpadj(o, classical=true)/2) / quantile(Normal(), α/2))
  end

	if isone(o.q)  # 2D plot
    α<=0 && (α = T(.05))  # if level=100, no CI constructed, but we need a reasonable α to choose graphing bounds

    if α > 0 && cols(o.v)-1 <= 1/α-1e6
      setquietly!(o, _quietly)
      !o.quietly && throw(ErrorException("need at least $(ceil(1/α)) replications to resolve a $(o.level)% two-sided confidence interval."))
    end

    p_lo, p_hi = T(NaN), T(NaN)
    if ismissing(o.gridmin[1]) || ismissing(o.gridmax[1])
      if o.B>0
        lo = ismissing(o.gridmin[1]) ? o.confpeak - halfwidth : o.gridmin[1]  # initial guess based on classical distribution
        hi = ismissing(o.gridmax[1]) ? o.confpeak + halfwidth : o.gridmax[1]
      else
        tmp = vec(sqrtNaN.(o.statDenom)) * cquantile(o.small ? TDist(o.df_r) : Normal(), α/2)
        lo = ismissing(o.gridmin[1]) ? o.confpeak - tmp : o.gridmin[1]
        hi = ismissing(o.gridmax[1]) ? o.confpeak + tmp : o.gridmax[1]
        if o.scoreBS && !o.null && !o.willplot  # if doing simple Wald test with no graph, we're done
          o.CI = [lo hi]
          return
        end
      end

      if abs(lo[1] - o.r[1]) > abs(hi[1] - o.r[1])  # brute force way to ensure that first trial bound tested is the farther one from r, for better interpolation
        if ismissing(o.gridmin[1]) && o.ptype≠lower  # unless lower-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
          for i ∈ 1:10
            p_lo = r_to_p(o, lo)
            p_lo < α && break
            tmp = hi - lo
            lo .-= tmp
            ismissing(o.gridmax[1]) && o.twotailed && (hi .+= tmp)  # maintain rough symmetry unless user specified upper bound
          end
        end
        if ismissing(o.gridmax[1]) && o.ptype≠upper  # ditto for high side
          for i ∈ 1:10
            p_hi = r_to_p(o, hi)
            p_hi < α && break
            tmp = hi - lo
            ismissing(o.gridmin[1]) && o.twotailed && (lo .-= tmp)
            hi .+= tmp
          end
        end
      else
        if ismissing(o.gridmax[1]) && o.ptype≠upper  # ditto for high side
          for i ∈ 1:10
            p_hi = r_to_p(o, hi)
            p_hi < α && break
            tmp = hi - lo
            ismissing(o.gridmin[1]) && o.twotailed && (lo .-= tmp)
            hi .+= tmp
          end
        end
        if ismissing(o.gridmin[1]) && o.ptype≠lower  # unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
          for i ∈ 1:10
            p_lo = r_to_p(o, lo)
            p_lo < α && break
            tmp = hi - lo
            lo .-= tmp
            ismissing(o.gridmax[1]) && o.twotailed && (hi .+= tmp)  # maintain rough symmetry unless user specified upper bound
          end
        end
      end
    else  # both grid bounds pre-specified
      lo = [o.gridmin[1]]
      hi = [o.gridmax[1]]
    end

    o.plotX = range(lo[1], hi[1], length=o.gridpoints[1])[:,:]
    o.plotY = fill(T(NaN), length(o.plotX))
    o.plotY[1  ] = p_lo
    o.plotY[end] = p_hi
    p_confpeak = T(o.WREnonARubin ? NaN : o.twotailed ? 1 : .5)

    c = clamp((floor(Int, (o.confpeak[1] - lo[1]) / (hi[1] - lo[1]) * (o.gridpoints[1] - 1)) + 2), 1, o.gridpoints[1]+1)  # insert original point estimate into grid
    o.plotX = [o.plotX[1:c-1] ; o.confpeak ; o.plotX[c:end]][:,:]
    o.plotY = [o.plotY[1:c-1] ; p_confpeak ; o.plotY[c:end]]
  else  # 2D plot
    lo = Vector{T}(undef, 2)
    hi = Vector{T}(undef, 2)
    for d ∈ 1:o.df
      lo[d] = ismissing(o.gridmin[d]) ? o.confpeak[d] - halfwidth[d] : o.gridmin[d]
      hi[d] = ismissing(o.gridmax[d]) ? o.confpeak[d] + halfwidth[d] : o.gridmax[d]

      # stata("_natscale " + strofreal(lo[d]) + " " + strofreal(hi[d]) + " 4")  # using Stata logic for setting graph bounds ensures good-looking contour plot
      # if (o.gridmin[d]==.)
      #   stata("local min = r(min)")  # for some reason st_global("r(min)") doesn't work
      #   lo[d] = strtoreal(st_local("min"))
      # end
      # if (o.gridmax[d]==.)
      #   stata("local max = r(max)")
      #   hi[d] = strtoreal(st_local("max"))
      # end
    end
    o.plotX = [repeat(range(lo[1], hi[1], length=o.gridpoints[1]), inner=o.gridpoints[2]) repeat(range(lo[2], hi[2], length=o.gridpoints[2]), outer=o.gridpoints[1])]
    o.plotY = fill(T(NaN), rows(o.plotX))
  end

  isnan(o.plotY[1]) && (o.plotY[1] = r_to_p(o, o.plotX[1,:]))  # do in this order for widest interpolation
  @views for i ∈ length(o.plotY):-1:2
		isnan(o.plotY[i]) && (o.plotY[i] = r_to_p(o, o.plotX[i,:]))
	end

  if isone(o.q) && o.level<100  # find CI bounds
    _CI = map(x->isnan(x) ? x : T(x > α), o.plotY); _CI = _CI[2:end] - _CI[1:end-1]
    lo = findall(x->x== 1, _CI)
    hi = findall(x->x==-1, _CI)
    if length(lo)==0 && length(hi)==0
      o.CI = [NaN NaN]
    else
      if iszero(length(lo))
        lo = [NaN]
      elseif iszero(length(hi))
        hi = [NaN]
      else
         lo[1  ] >  hi[1  ] && (lo = [NaN ; lo ]) # non-rejection ranges that are not within grid range
        -lo[end] < -hi[end] && (hi = [hi  ; NaN])
      end
      o.CI = [lo hi]

      for i ∈ 1:length(lo), j ∈ 1:2
        if !isnan(o.CI[i,j])
          t = Int(o.CI[i,j])
          o.CI[i,j] = search(o, α, o.plotY[t], o.plotX[t], o.plotY[t+1], o.plotX[t+1])
        end
      end
    end
  end

  if @isdefined c  # now that it's done helping graph look good, remove peak point from returned grid for evenness, for Bayesian sampling purposes
    peak = o.plotX[c], o.plotY[c]
    o.plotX = o.plotX[[1:c-1; c+1:rows(o.plotX)],:]
    deleteat!(o.plotY, c)
  end

	setquietly!(o, _quietly); o.r = _r; o.dirty = true  # restore backups
	o.notplotted = false
end

# return matrix whose rows are all the subsets of a row of numbers. Nil is at bottom.
function combs(d::Integer)
  retval = Array{Bool}(undef, 2^d, 0)
	for i ∈ d:-1:1
		retval = [retval kron(kron(trues(2^(d-i),1), [true;false]), trues(2^(i-1),1))]
  end
	retval
end

end # module

using StatFiles, StatsModels, DataFrames, prBenchmarkTools, Plots, CategoricalArrays

df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\collapsed.dta"))
dropmissing!(df)
f = @formula(hasinsurance ~ 1 + selfemployed + post + post_self)
f = apply_schema(f, schema(f, df, Dict(:hasinsurance => ContinuousTerm)))
resp, pred = modelcols(f, df)

println("\nboottest post_self=.04, weight(webb)")
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, resp)
Boottest.setX₁!(M, pred)
Boottest.setID!(M, df.year)
Boottest.setsmall!(M, true)
Boottest.setR!(M, [0 0 0 1.], [.04])
Boottest.setB!(M, 999)
Boottest.setauxwttype!(M,Boottest.webb)
Boottest.setgrid!(M, [missing], [missing], [missing])
CI = Boottest.getCI(M)
t = Boottest.getstat(M)
p = Boottest.getp(M)
dist = Boottest.getdist(M)
histogram(dist)
println("t=$t p=$p CI=$CI")

println("\nboottest post_self=.04, weight(webb) reps(9999999) noci")
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, resp)
Boottest.setX₁!(M, pred)
Boottest.setID!(M, df[:,:year])
Boottest.setsmall!(M, true)
Boottest.setR!(M, [0 0 0 1.], [.04])
Boottest.setB!(M, 9999999)
Boottest.setauxwttype!(M,Boottest.webb)
println("t=$(Boottest.getstat(M)) p=$(Boottest.getp(M))")

println("\nregress hasinsurance selfemployed post post_self, cluster(year)")
println("boottest (post_self=.05) (post=-.02), reps(9999) weight(webb)")
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, resp)
Boottest.setX₁!(M, pred)
Boottest.setID!(M, df.year)
Boottest.setsmall!(M, true)
Boottest.setB!(M, 9999)
Boottest.setR!(M, [0 0 0 1.; 0 0 1. 0], [0.05; -0.02])
Boottest.setauxwttype!(M,Boottest.webb)
Boottest.setgrid!(M, [.02; -.03], [.08; -.01], [missing; missing])
x,y = Boottest.getplot(M)
println("F=$(Boottest.getstat(M)) p=$(Boottest.getp(M))")
gr()
plot(contour(x[25:25:625,1], x[1:25,2], reshape(y,25,25), fill=true))

df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))[:,[:wage; :tenure; :ttl_exp; :collgrad; :industry]]
dropmissing!(df)
desc = describe(df, :eltype)
for i in axes(desc, 1)  # needed only for lm()
  if desc[i,:eltype] == Float32
    sym = desc[i,:variable]
    @transform!(df, @byrow $sym=Float64($sym))
  end
end
sort!(df, :industry)
f = @formula(wage ~ 1 + tenure + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
resp, pred = modelcols(f, df)

println("\nconstraint 1 ttl_exp = .2")
println("cnsreg wage tenure ttl_exp collgrad, constr(1) cluster(industry)")
println("boottest tenure")
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, resp)
Boottest.setX₁!(M, pred)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, true)
Boottest.setR₁!(M, [0 0 1. 0], [.2])
Boottest.setR!(M, [0 1. 0 0], [.0])
Boottest.setB!(M, 999)
Boottest.setgrid!(M, [missing], [missing], [missing])
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
t = Boottest.getstat(M)
p = Boottest.getp(M)
dist = Boottest.getdist(M)
plot(x,y)
println("t=$t p=$p CI=$CI")

println("\nivregress 2sls wage ttl_exp collgrad (tenure = union), cluster(industry)")
println("boottest tenure, ptype(equaltail)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setR!(M, [0 0 0 1.], [.0])
Boottest.setB!(M, 999)
Boottest.setptype!(M, Boottest.equaltail)
Boottest.setgrid!(M, [missing], [missing], [missing])
x,y = Boottest.getplot(M)
t = Boottest.getstat(M)
p = Boottest.getp(M)
CI = Boottest.getCI(M)
println("z=$t p=$p CI=$CI")
plot(x,y)

println("\nboottest tenure, ptype(equaltail) reps(99999) weight(webb) stat(c)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setR!(M, [0 0 0 1.], [.0])
Boottest.setB!(M, 99999)
Boottest.setptype!(M, Boottest.equaltail)
Boottest.setgrid!(M, [missing], [missing], [missing])
Boottest.setauxwttype!(M,Boottest.webb)
Boottest.setstattype!(M,Boottest.c)
t = Boottest.getstat(M)
p = Boottest.getp(M)
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
dist = Boottest.getdist(M)
histogram(dist)
println("z=$t p=$p CI=$CI")
plot(x,y)

println("\nboottest, ar")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setARubin!(M, true)
Boottest.setR!(M, [0 0 0 1.], [.0])
Boottest.setB!(M, 999)
Boottest.setptype!(M, Boottest.equaltail)
Boottest.setgrid!(M, [missing], [missing], [missing])
x,y = Boottest.getplot(M)
t = Boottest.getstat(M)
p = Boottest.getp(M)
CI = Boottest.getCI(M)
dist = Boottest.getdist(M)
histogram(dist)
println("z=$t p=$p CI=$CI")
plot(x,y)

println("\nscoretest tenure")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setR!(M, [0 0 0 1.], [.0])
Boottest.setB!(M, 0)
Boottest.setscoreBS!(M,true)
Boottest.setptype!(M, Boottest.equaltail)
Boottest.setgrid!(M, Vector{Union{Missing, Float32}}([0.]), Vector{Union{Missing, Float32}}([1.5]), [missing])
t = Boottest.getstat(M)
p = Boottest.getp(M)
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
println("z=$t p=$p CI=$CI")
plot(x,y)

println("\nwaldtest tenure")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage; :tenure; :ttl_exp; :collgrad; :industry; :union]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + ttl_exp + collgrad)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setnull!(M, false)
Boottest.setR!(M, [0 0 0 1.], [.0])
Boottest.setB!(M, 0)
Boottest.setscoreBS!(M,true)
Boottest.setptype!(M, Boottest.equaltail)
Boottest.setgrid!(M, [missing], [missing], [missing])
t = Boottest.getstat(M)
p = Boottest.getp(M)
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
println("z=$t p=$p CI=$CI")
plot(x,y)

println("\nivregress liml wage (tenure = collgrad ttl_exp), cluster(industry)")
println("boottest tenure")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :tenure, :ttl_exp, :collgrad, :industry]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ collgrad + ttl_exp)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setR!(M, [0 1.], [.0])
Boottest.setB!(M, 999)
Boottest.setLIML!(M, true)
Boottest.setgrid!(M, [missing], [missing], [missing])
t = Boottest.getstat(M)
p = Boottest.getp(M)
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
println("z=$t p=$p CI=$CI")
dist = Boottest.getdist(M)
histogram(dist)
plot(x,y)

println("\nivreg2 wage collgrad smsa race age (tenure = union married), cluster(industry) fuller(1)")
println("boottest tenure, nograph weight(webb)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry]]
dropmissing!(df)
sort!(df, :industry)
f = @formula(wage ~ 1 + collgrad + smsa + race + age)
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
ivf = @formula(tenure ~ union + married)
ivf = apply_schema(ivf, schema(ivf, df))
Y2, X2 = modelcols(ivf, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setY₂!(M, Y2)
Boottest.setX₂!(M, X2)
Boottest.setID!(M, df.industry)
Boottest.setsmall!(M, false)
Boottest.setR!(M, [0 0 0 0 0 1.], [.0])
Boottest.setauxwttype!(M,Boottest.webb)
Boottest.setB!(M, 999)
Boottest.setLIML!(M, true)
Boottest.setFuller!(M, true)
Boottest.setgrid!(M, [missing], [missing], [missing])
t = Boottest.getstat(M)
p = Boottest.getp(M)
x,y = Boottest.getplot(M)
CI = Boottest.getCI(M)
println("z=$t p=$p CI=$CI")
dist = Boottest.getdist(M)
histogram(dist)
plot(x,y)

println("\nareg wage ttl_exp collgrad tenure [aw=hours] if occupation<., cluster(age) absorb(industry)")
println("boottest tenure, cluster(age occupation) bootcluster(occupation)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :ttl_exp, :collgrad, :tenure, :age, :industry, :occupation, :hours]]
dropmissing!(df)
sort!(df, [:occupation, :age])
f = @formula(wage ~ ttl_exp + collgrad + tenure)  # constant unneeded in FE model
f = apply_schema(f, schema(f, df))
Y1, X1 = modelcols(f, df)
M = Boottest.StrBoottest{Float32}()
Boottest.sety₁!(M, Y1)
Boottest.setX₁!(M, X1)
Boottest.setID!(M, Matrix(df[:, [:occupation, :age]]), NBootClustVar=1, NErrClust=2)
Boottest.setFEID!(M, df.industry)
Boottest.setsmall!(M, true)
Boottest.setobswt!(M, df.hours, false)
Boottest.setR!(M, [0 0 1.], [.0])
Boottest.setB!(M, 999)
Boottest.setgrid!(M, [missing], [missing], [missing])
x,y = Boottest.getplot(M)
t = Boottest.getstat(M)
p = Boottest.getp(M)
CI = Boottest.getCI(M)
println("t=$t p=$p CI=$CI")
dist = Boottest.getdist(M)
histogram(dist)
plot(x,y)

println("\nglobal pix lnkm pixpetro pixdia pixwaterd pixcapdist pixmal pixsead pixsuit pixelev pixbdist")
println("global geo lnwaterkm lnkm2split mean_elev mean_suit malariasuit petroleum diamondd")
println("global poly capdistance1 seadist1 borderdist1")
println("encode pixwbcode, gen(ccode)  // make numerical country identifier")
println("qui areg lnl0708s centr_tribe lnpd0 \$pix \$geo \$poly, absorb(ccode)")
println("boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(ccode)")
println("boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(pixcluster)")
println("boottest centr_tribe, nogr reps(9999) clust(ccode pixcluster) bootcluster(ccode pixcluster)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Work\Econometrics\Wild cluster\pixel-level-baseline-final.dta"))
@time begin
  pix  = [:lnkm, :pixpetro, :pixdia, :pixwaterd, :pixcapdist, :pixmal, :pixsead, :pixsuit, :pixelev, :pixbdist]
  geo  = [:lnwaterkm, :lnkm2split, :mean_elev, :mean_suit, :malariasuit, :petroleum, :diamondd]
  poly = [:capdistance1, :seadist1, :borderdist1]
  df = df[:,[pix; geo; poly; :lnl0708s; :centr_tribe; :lnpd0; :pixwbcode; :pixcluster]]
  dropmissing!(df)
  df.ccode = levelcode.(categorical(df.pixwbcode, compress=true))
  df.pixcode = levelcode.(categorical(df.pixcluster, compress=true))
  f = Term(:lnl0708s) ~ sum(term.([:centr_tribe; :lnpd0; pix; geo; poly]))
  f = apply_schema(f, schema(f, df))
  sort!(df, [:ccode, :pixcode])
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, Matrix(df[:, [:ccode, :pixcode]]), NBootClustVar=1, NErrClust=2)
  Boottest.setFEID!(M, df.ccode)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [1 zeros(1,size(pred,2)-1)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [missing])
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  M = Boottest.StrBoottest{Float32}()
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [1 zeros(1,size(pred,2)-1)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [missing])
  sort!(df, [:pixcode, :ccode])
  resp, pred = modelcols(f, df)
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, Matrix(df[:, [:pixcode, :ccode]]), NBootClustVar=1, NErrClust=2)
  Boottest.setFEID!(M, df.ccode)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  M = Boottest.StrBoottest{Float32}()
  sort!(df, [:pixcode, :ccode])
  resp, pred = modelcols(f, df)
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [1 zeros(1,size(pred,2)-1)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [missing])
  Boottest.setID!(M, Matrix(df[:, [:pixcode, :ccode]]), NBootClustVar=2, NErrClust=2)
  Boottest.setFEID!(M, df.ccode)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")
end

println("\ninfile coll merit male black asian year state chst using regm.raw, clear")
println("qui regress coll merit male black asian i.year i.state, cluster(state)	")
println("generate individual = _n  // unique ID for each observation")
println("boottest merit, nogr reps(999) gridpoints(10)  // defaults to bootcluster(state)")
println("boottest merit, nogr reps(999) gridpoints(10) nonull")
println("boottest merit, nogr reps(999) gridpoints(10) bootcluster(state year)")
println("boottest merit, nogr reps(999) gridpoints(10) nonull bootcluster(state year)")
println("boottest merit, nogr reps(999) gridpoints(10) bootcluster(individual)")
println("boottest merit, nogr reps(999) gridpoints(10) nonull bootcluster(individual)")
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Work\Econometrics\Wild cluster\regm.dta"))
@time begin
  df = DataFrame(coll=Bool.(df.coll), merit=Bool.(df.merit), male=Bool.(df.male), black=Bool.(df.black), asian=Bool.(df.asian), state=categorical(Int8.(df.state)), year=categorical(Int16.(df.year)))
  dropmissing!(df)
  df = df[df.state .∉ Ref([34,57,59,61,64,71,72,85,88]),:]  # restrict treatment group to Georgia
  f = @formula(coll ~ 1 + merit + male + black + asian + year + state)
  f = apply_schema(f, schema(f, df))
  sort!(df, [:state])
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, levelcode.(df.state))
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, levelcode.(df.state))
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  Boottest.setnull!(M, false)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  sort!(df, [:year, :state])
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, [levelcode.(df.year) levelcode.(df.state)], NBootClustVar=2, NErrClust=1)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  Boottest.setnull!(M, true)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  sort!(df, [:year, :state])
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, [levelcode.(df.year) levelcode.(df.state)], NBootClustVar=2, NErrClust=1)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  Boottest.setnull!(M, false)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  sort!(df, :state)
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, [collect(1:nrow(df)) levelcode.(df.state)], NBootClustVar=1, NErrClust=1)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")

  sort!(df, :state)
  resp, pred = modelcols(f, df)
  M = Boottest.StrBoottest{Float32}()
  Boottest.sety₁!(M, resp)
  Boottest.setX₁!(M, pred)
  Boottest.setID!(M, [collect(1:nrow(df)) levelcode.(df.state)], NBootClustVar=1, NErrClust=1)
  Boottest.setsmall!(M, true)
  Boottest.setR!(M, [0. 1 zeros(1,size(pred,2)-2)], [0.])
  Boottest.setB!(M, 9999)
  Boottest.setgrid!(M, [missing], [missing], [10])
  Boottest.setnull!(M, false)
  CI = Boottest.getCI(M)
  t = Boottest.getstat(M)
  p = Boottest.getp(M)
  println("t=$t p=$p CI=$CI")
end

# println("ivregress 2sls D.lViolentpop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year i.state, robust")
# println("boottest DL.lpris_totpop, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) gridmin(-2) gridmax(2) nogr")
# println("ivregress 2sls D.lPropertypop (DL.lpris_totpop = ibnL.stage#i(1/3)L.substage) D.(lincomepop unemp lpolicepop metrop black a*pop) i.year i.state, robust")
# println("boottest DL.lpris_totpop, cluster(state year) bootcluster(year) ptype(equaltail) reps(999) gridmin(-2) gridmax(2) nogr")

using Profile
df = DataFrame(load(raw"C:\Users\drood\OneDrive\Documents\Macros\nlsw88.dta"))
df = df[:, [:wage, :ttl_exp, :collgrad, :tenure, :age, :industry, :occupation, :hours]]
Profile.clear()
@profile begin
# @profview begin
  dropmissing!(df)
	sort!(df, [:occupation, :age])
	f = @formula(wage ~ ttl_exp + collgrad + tenure)  # constant unneeded in FE model
	f = apply_schema(f, schema(f, df))
	Y1, X1 = modelcols(f, df)
	M = Boottest.StrBoottest{Float32}()
	Boottest.sety₁!(M, Y1)
	Boottest.setX₁!(M, X1)
	Boottest.setID!(M, Matrix(df[:, [:occupation, :age]]), NBootClustVar=1, NErrClust=2)
	Boottest.setFEID!(M, df.industry)
	Boottest.setsmall!(M, true)
	Boottest.setobswt!(M, df.hours, false)
	Boottest.setR!(M, [0 0 1.], [.0])
	Boottest.setB!(M, 9999)
	Boottest.setgrid!(M, [missing], [missing], [missing])
	x,y = Boottest.getplot(M)
	t = Boottest.getstat(M)
	p = Boottest.getp(M)
	CI = Boottest.getCI(M)
end
Juno.profiler()
