# WildBootTest.jl 0.1 2 November 2021

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


module WildBootTest
export BoottestResult, wildboottest, AuxWtType, PType, MAdjType, teststat, stattype, p, padj, reps, repsfeas, NBootClust, dof, dof_r, plotpoints, peak, CI, dist, statnumer, statvar, auxweights

using LinearAlgebra, Random, Distributions, LoopVectorization, SortingAlgorithms

"Auxilliary weight types: rademacher, mammen, webb, normal, gamma"
@enum AuxWtType rademacher mammen webb normal gamma

"p value types: symmetric, equaltail, lower, upper"
@enum PType symmetric equaltail lower upper

"Multiple hypothesis adjustment types: none, bonferroni, sidak"
@enum MAdjType none bonferroni sidak

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
#   r = nrows(M)
# 	return (all(colsum(M) .== 1) && all(colsum(!M) .== r-1)?  # is M 0's but for one 1 in each col, so it's just copying/reordering ncols?
#              (all(diagonal(M)) && nrows(M)==ncols(M)?  # is M just the identity matrix?
#                 &X            :
#                 &X[:,colsum(M.*(1::r))]) :  # reorder X
#              &(X * M))
# }

# (wt of type UniformScaling <=> wt=I)
# @inline isone(X::AbstractArray) = X==ones(eltype(X),ones(Int64, ndims(X))...)
@inline sqrtNaN(x) = x<0 ? typeof(x)(NaN) : sqrt(x)
@inline invsym(X) = iszero(length(X)) ? X : inv(Symmetric(X))
@inline eigensym(X) = eigen(Symmetric(X))  # does eigen recognize symmetric matrices?
@inline symcross(X::AbstractVecOrMat, wt::Union{UniformScaling,AbstractVector}) = Symmetric(cross(X,wt,X))  # maybe bad name since it means cross product in Julia
@inline function cross(X::AbstractVecOrMat{T}, wt::AbstractVector{T}, Y::AbstractVecOrMat{T}) where T
  dest = Matrix{T}(undef, ncols(X), ncols(Y))
  mul!(dest, X', wt.*Y)
end
@inline function cross(X::AbstractVecOrMat{T}, wt::UniformScaling, Y::AbstractVecOrMat{T}) where T
  dest = Matrix{T}(undef, ncols(X), ncols(Y))
  mul!(dest, X', Y)
end
@inline function crossvec(X::AbstractMatrix{T}, wt::AbstractVector{T}, Y::AbstractVector{T}) where T
  dest = Vector{T}(undef, ncols(X))
  mul!(dest, X', wt.*Y)
end
@inline function crossvec(X::AbstractMatrix{T}, wt::UniformScaling, Y::AbstractVector{T}) where T
  dest = Vector{T}(undef, ncols(X))
  mul!(dest, X', Y)
end
@inline vHadw(v::AbstractArray, w::AbstractVector) = v .* w
@inline vHadw(v::AbstractArray, w::UniformScaling) = v
@inline nrows(X::AbstractArray) = size(X,1)
@inline ncols(X::AbstractArray) = size(X,2)
@inline colsum(X::AbstractArray) = iszero(length(X)) ? Matrix{eltype(X)}(undef,1,0) : isone(ncols(X)) ? hcat(sum(X)) : sum(X, dims=1)
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))
@inline wtsum(wt::AbstractArray, X::AbstractArray) = wt'X
@inline wtsum(wt::UniformScaling, X::AbstractArray) = sum(X,dims=1)
# checkI!(X::AbstractArray) = all(abs.(X - I) .< 10eps(eltype(X))) ? I : X
@inline X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractMatrix) = @views X₁*B[1:size(X₁,2),:] + X₂*B[size(X₁,2)+1:end,:]
@inline X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractVector) = @views X₁*B[1:size(X₁,2)  ] + X₂*B[size(X₁,2)+1:end  ]

import Base.*  # extend * to left- and right-multiply arrays by vec or mat
@inline *(A::AbstractArray, B::AbstractVecOrMat) = reshape(reshape(A,:,size(B,1)) * B, size(A)[1:end-1]..., size(B)[2:end]...)
@inline *(A::AbstractVecOrMat, B::AbstractArray) = reshape(A * reshape(B,size(A,2),:), size(A)[1:end-1]..., size(B)[2:end]...)

function coldot!(dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)  # colsum(A .* B)
  s2 = axes(A,2)
  @turbo for i ∈ s2
	  dest[i] = A[1,i] * B[1,i]
  end
  if size(A,1)>1
	  @turbo for i ∈ s2, j ∈ 2:size(A,1)
	    dest[i] += A[j,i] * B[j,i]
	  end
  end
end
function coldot!(dest::AbstractMatrix, A::AbstractMatrix)  # colsum(A .* A)
  s2 = axes(A,2)
  @turbo for i ∈ s2
	  dest[i] = A[1,i] ^ 2
  end
  if size(A,1)>1
	  @turbo for i ∈ s2, j ∈ 2:size(A,1)
	    dest[i] += A[j,i] ^ 2
	  end
  end
end
function coldot!(dest::AbstractArray, A::AbstractArray, B::AbstractArray)
  J = CartesianIndices(axes(A)[2:end])
  eachindexJ = eachindex(J)
  @turbo for j ∈ eachindexJ
	  dest[1,J[j]] = A[1,J[j]] * B[1,J[j]]
  end
  if size(A,1) > 1
	  @turbo for j ∈ eachindexJ, i ∈ 2:size(A,1)
	    dest[1,J[j]] += A[i,J[j]] * B[i,J[j]]
	  end
  end
end
function coldot!(dest::AbstractArray, A::AbstractArray)
  J = CartesianIndices(axes(A)[2:end])
  eachindexJ = eachindex(J)
  @turbo for j ∈ eachindexJ
	  dest[1,J[j]] = A[1,J[j]] ^ 2
  end
  if size(A,1) > 1
	  @turbo for j ∈ eachindexJ, i ∈ 2:size(A,1)
	    dest[1,J[j]] += A[i,J[j]] ^ 2
	  end
  end
end
function coldot(args...)
  dest = Array{promote_type(eltype.(args)...)}(undef, 1, size(args[1])[2:end]...)
  coldot!(dest, args...)
  dest
 end
# coldot(A::AbstractMatrix, B::AbstractVector) = coldot(B,A)
coldot(A::AbstractVector, B::AbstractVector) = [dot(A,B)]
function coldotplus!(dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)  # colsum(A .* B)
  @turbo for i ∈ axes(A,2), j ∈ axes(A,1)
	  dest[i] += A[j,i] * B[j,i]
  end
end

# compute the norm of each col of A using quadratic form Q
function colquadform(Q::AbstractMatrix, A::AbstractMatrix) :: AbstractVector
  dest = zeros(promote_type(eltype(Q), eltype(A)), size(A,2))
  @turbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
	  dest[i] += A[j,i] * Q[k,j] * A[k,i]
  end
  dest
end

 # From given row of given matrix, substract inner products of corresponding ncols of A & B with quadratic form Q
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

# unweighted panelsum!() along first axis of an Array
function panelsum!(dest::AbstractArray, X::AbstractArray, info::Vector{UnitRange{T}} where T<:Integer)
    iszero(length(X)) && return
    J = CartesianIndices(axes(X)[2:end])
    eachindexJ = eachindex(J)
    @inbounds for g in eachindex(info)
      f, l = first(info[g]), last(info[g])
      @turbo for j ∈ eachindexJ
        dest[g,J[j]] = X[f,J[j]]
      end
      if f<l
        @turbo for j ∈ eachindexJ, i ∈ f+1:l
          dest[g,J[j]] += X[i,J[j]]
        end
      end
    end
  end
  # single-weighted panelsum!() along first axis of an Array
  function panelsum!(dest::AbstractArray, X::AbstractArray, wt::AbstractVector, info::Vector{UnitRange{T}} where T<:Integer)
    iszero(length(X)) && return
    if iszero(length(info)) || nrows(info)==nrows(X)
      dest .= X .* wt
      return
    end
    J = CartesianIndices(axes(X)[2:end])
    eachindexJ = eachindex(J)
    @inbounds for g in eachindex(info)
      f, l = first(info[g]), last(info[g])
      _wt = wt[f]
      @turbo for j ∈ eachindexJ
        dest[g,J[j]] = X[f,J[j]] * _wt
      end
      if f<l
       	@turbo for j ∈ eachindexJ, i ∈ f+1:l
          dest[g,J[j]] += X[i,J[j]] * wt[i]
        end
      end
    end
  end
  # multiple-weighted panelsum!() along first axis of an Array
  # *2nd* dimension of resulting 3-D array corresponds to cols of v=wt;
  # this facilitates reshape() to 2-D array in which results for each col of v are stacked vertically
  function panelsum!(dest::AbstractArray, X::AbstractArray, wt::AbstractMatrix, info::Vector{UnitRange{T}} where T<:Integer)
    iszero(length(X)) && return
    if iszero(length(info)) || nrows(info)==nrows(X)
      for i ∈ axes(wt,2)
        dest[:,i,:] .= X .* view(wt,:,i)
      end
      return
    end
    J = CartesianIndices(axes(X)[2:end])
    eachindexJ = eachindex(J)
    @inbounds for g in eachindex(info)
      f, l = first(info[g]), last(info[g])
      fl = f+1:l
      for k ∈ axes(wt,2)
        _wt = wt[f,k]
        @turbo for j ∈ eachindexJ
          dest[g,k,J[j]] = X[f,J[j]] * _wt
        end
        if f<l
          @turbo for j ∈ eachindexJ, i ∈ fl
            dest[g,k,J[j]] += X[i,J[j]] * wt[i,k]
          end
        end
      end
    end
  end
  function panelsum(X::AbstractArray, wt::AbstractVecOrMat, info::Vector{UnitRange{T}} where T<:Integer)
    dest = similar(X, length(info), size(wt)[2:end]..., size(X)[2:end]...)
    panelsum!(dest, X, wt, info)
    dest
  end
  function panelsum(X::AbstractArray, info::Vector{UnitRange{T}} where T<:Integer)
    dest = similar(X, length(info), size(X)[2:end]...)
    panelsum!(dest, X, info)
    dest
  end
  function panelsum2(X₁::AbstractArray, X₂::AbstractArray, wt::AbstractVecOrMat, info::Vector{UnitRange{T}} where T<:Integer)
    if iszero(length(X₁))
      panelsum(X₂,wt,info)
    elseif iszero(length(X₂))
      panelsum(X₁,wt,info)
    else
      dest = similar(X₁, length(info), size(wt)[2:end]..., ncols(X₁)+ncols(X₂))
      panelsum!(view(dest, Vector{Colon}(undef,ndims(wt))...,           1:ncols(X₁  )), X₁, wt, info)
      panelsum!(view(dest, Vector{Colon}(undef,ndims(wt))..., ncols(X₁)+1:size(dest)[end]), X₂, wt, info)
      dest
    end
  end
  @inline panelsum(X::AbstractArray, wt::UniformScaling, info::Vector{UnitRange{T}} where T<:Integer) = panelsum(X, info)
  # macros to efficiently handle result = input
  macro panelsum(X, info)
    :( iszero(length($(esc(info)))) || length($(esc(info)))==nrows($(esc(X))) ? $(esc(X)) : panelsum($(esc(X)), $(esc(info)) ) )
  end
  macro panelsum(X, wt, info)
    :( panelsum($(esc(X)), $(esc(wt)), $(esc(info))) )
  end
  macro panelsum2(X₁, X₂, wt, info)
    :( panelsum2($(esc(X₁)), $(esc(X₂)), $(esc(wt)), $(esc(info))) )
  end

abstract type Estimator end
struct OLS<:Estimator end
struct ARubin<:Estimator end
struct IVGMM<:Estimator end

mutable struct StrEstimator{T<:AbstractFloat, E<:Estimator}
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
	WXAR::Matrix{T}; ScapPXYZperp::Vector{Matrix{T}}; ScapYX::Vector{Matrix{T}}; CT_XAR::Array{T,3}; CT_FEcapPY::Vector{Matrix{T}}

  # IV/GMM only
  ZZ::Matrix{T}; XY₂::Matrix{T}; XX::Matrix{T}; H_2SLS::Matrix{T}; V::Matrix{T}; ZY₂::Matrix{T}; X₂Y₂::Matrix{T}; X₁Y₂::Matrix{T}; ZR₁ZR₁::Matrix{T}; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZZR₁::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; ZXinvXXXZ::Matrix{T}; H_2SLSmZZ::Matrix{T}
  ZXinvXXXy₁par::Vector{T}; t₁Y::Vector{T}
  Y₂y₁::Vector{T}; twoR₁Zy₁::Vector{T}
  y₁y₁::T; y₁pary₁par::T
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}
  Y₂y₁par::Vector{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}
  kX::Integer

  StrEstimator{T,E}(parent) where T<:AbstractFloat where E<:Estimator = new(parent, true, E==IVGMM, false, T(E==IVGMM ? 1 : 0), Matrix{T}(undef,0,0), I)
end

matconvert(T::DataType, X) = !isa(X, AbstractArray) || eltype(X)==T ? X : T.(X)

mutable struct StrBootTest{T<:AbstractFloat}
  R::Matrix{T}; r::Vector{T}; R₁::Matrix{T};r₁::Vector{T}
  y₁::Vector{T}; X₁::VecOrMat{T}; Y₂::VecOrMat{T}; X₂::VecOrMat{T}
  wt::Union{Vector{T}, UniformScaling}; fweights::Bool
  LIML::Bool; Fuller::T; κ::T; ARubin::Bool
  B::Int64; auxtwtype::AuxWtType; rng::AbstractRNG; maxmatsize::Float16
  ptype::PType; null::Bool; bootstrapt::Bool
  ID::VecOrMat{Int64}; nbootclustvar::Int8; nerrclustvar::Int64; small::Bool
  FEID::Vector{Int64}; FEdfadj::Bool
  level::T; rtol::T
  madjtype::MAdjType; NH0::Int16
  ML::Bool; β::Vector{T}; A::Matrix{T}; sc::VecOrMat{T}
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
		purerobust::Bool; Nstar::Int64; Nw::Int64; enumerate::Bool; interpolable::Bool; interpolate_u::Bool; kX₂::Int64; kX::Int64
  _FEID::Vector{Int64}; AR::Matrix{T}; v::Matrix{T}; ustar::Matrix{T}; CT_WE::Matrix{T}
  infoBootData::Vector{UnitRange{Int64}}; infoBootAll::Vector{UnitRange{Int64}}; infoErrAll::Vector{UnitRange{Int64}}
  JNcapNstar::Matrix{T}; statDenom::Matrix{T}; uXAR::Matrix{T}; SuwtXA::Matrix{T}; numer₀::Matrix{T}; βdev::Matrix{T}; δdenom_b::Matrix{T}; _Jcap::Matrix{T}; YYstar_b::Matrix{T}; YPXYstar_b::Matrix{T}; numerw::Matrix{T}
	distCDR::Matrix{T}; plotX::Matrix{T}; plotY::Vector{T}; ClustShare::Vector{T}; WeightGrp::Vector{UnitRange{Int64}}
  numersum::Vector{T}; ü₀::Vector{T}; invFEwt::Vector{T}
	βs::Matrix{T}; As::Matrix{T}
	infoAllData::Vector{UnitRange{Int64}}; infoCapData::Vector{UnitRange{Int64}}; IDAll::Matrix{T}; Ü₂par::Matrix{T}
	ü::Vector{T}
	DGP::StrEstimator{T,E} where E; Repl::StrEstimator{T,E} where E; M::StrEstimator{T,E} where E
	clust::Vector{StrClust{T}}
	denom::Matrix{Matrix{T}}; Kcd::Matrix{Matrix{T}}; Jcd::Matrix{Matrix{T}}; denom₀::Matrix{Matrix{T}}; Jcd₀::Matrix{Matrix{T}}; SCTcapuXinvXX::Matrix{Matrix{T}}; SstarUU::Matrix{Vector{T}}; CTUX::Matrix{Matrix{T}}
	∂u∂r::Vector{Matrix{T}}; ∂numer∂r::Vector{Matrix{T}}; IDCTCapstar::Vector{Vector{Int64}}; infoCTCapstar::Vector{Vector{UnitRange{Int64}}}; SstarUX::Vector{Matrix{T}}; SstarUXinvXX::Vector{Matrix{T}}; SstarUZperpinvZperpZperp::Vector{Matrix{T}}; SstaruY::Vector{Matrix{T}}; SstarUMZperp::Vector{Matrix{T}}; SstarUPX::Vector{Matrix{T}}; SstarUZperp::Vector{Matrix{T}}; CTFEU::Vector{Matrix{T}}
  ∂denom∂r::Vector{Matrix{Matrix{T}}}; ∂Jcd∂r::Vector{Matrix{Matrix{T}}}
  ∂²denom∂r²::Matrix{Matrix{Matrix{T}}}
	FEs::Vector{StrFE{T}}
  T1L::Vector{Matrix{T}}; T1R::Vector{Matrix{T}}
	crosstabCapstarind::Vector{Int64}

  StrBootTest{T}(R, r, R₁, r₁, y₁, X₁, Y₂, X₂, wt, fweights, LIML, Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, scorebs, bootstrapt, ID, nbootclustvar, nerrclustvar, robust, small, FEID, FEdfadj, level, rtol, madjtype, NH0, ML, β, A, sc, willplot, gridmin, gridmax, gridpoints) where T<:Real =
	  new(matconvert(T,R), matconvert(T,r), matconvert(T,R₁), matconvert(T,r₁), matconvert(T,y₁), matconvert(T,X₁), matconvert(T,Y₂), matconvert(T,X₂), matconvert(T,wt), fweights, LIML || !iszero(Fuller), Fuller, κ, ARubin, B, auxtwtype, rng, maxmatsize, ptype, null, bootstrapt, matconvert(Int64,ID), nbootclustvar, nerrclustvar, small, FEID, FEdfadj, level, rtol, madjtype, NH0, ML, matconvert(T,β), matconvert(T,A), matconvert(T,sc), willplot, gridmin, gridmax, gridpoints,
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

function perp(A::AbstractMatrix)
  F = eigensym(A*invsym(A'A)*A')
  F.vectors[:, abs.(F.values) .< 1000eps(eltype(A))]
  # checkI!(dest)
end

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 nrows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T,E}, R₁::AbstractMatrix{T}, R::Union{UniformScaling{Bool},AbstractMatrix{T}}=Matrix{T}(undef,0,0)) where {T,E}
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
	  RR₁perp = Matrix([R ; zeros(T, o.parent.kY₂, o.parent.kX₁) I])  # nrows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

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
function InitVars!(o::StrEstimator{T,OLS}, Rperp::AbstractMatrix{T}) where T # Rperp is for replication regression--no null imposed
  o.y₁par = o.parent.y₁
  H = symcross(o.parent.X₁, o.parent.wt)
  o.invH = inv(H)

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.invH : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
  o.β₀ = Vector{T}(undef, nrows(R₁AR₁)); mul!(o.β₀, R₁AR₁, crossvec(o.parent.X₁, o.parent.wt, o.y₁par))
  o.∂β∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

  o.A = iszero(length(Rperp)) ? o.invH : Rperp * invsym(Rperp'H*Rperp) * Rperp'   # for replication regression
  o.AR = o.A * o.parent.R'
  (o.parent.scorebs || o.parent.robust) && (o.XAR = o.parent.X₁ * o.AR)
end

function InitVars!(o::StrEstimator{T,ARubin}, Rperp::AbstractMatrix{T} = Matrix{T}(undef,0,0)) where T
  X₂X₁ = cross(o.parent.X₂, o.parent.wt, o.parent.X₁)
  H = Symmetric([symcross(o.parent.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.parent.X₂, o.parent.wt)])  # XXX use LazyArrays?
  o.A = inv(H)
  o.AR = o.A * o.parent.R'
  (o.parent.scorebs || o.parent.robust) && (o.XAR = X₁₂B(o.parent.X₁, o.parent.X₂, o.AR))

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.A : o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp'
  o.β₀   = R₁AR₁ * [crossvec(o.parent.X₁, o.parent.wt, o.parent.y₁) ; crossvec(o.parent.X₂, o.parent.wt, o.parent.y₁)]  # XXX use LazyArrays?
  o.∂β∂r = R₁AR₁ * [cross(o.parent.X₁, o.parent.wt, o.parent.Y₂) ; cross(o.parent.X₂, o.parent.wt, o.parent.Y₂)]
end

function InitVars!(o::StrEstimator{T,IVGMM}, Rperp::AbstractMatrix{T}...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

  o.Zperp = o.parent.X₁ * o.RperpX  # XXX XB(o.parent.X₁, o.RperpX)
  o.invZperpZperp = iszero(length(o.Zperp)) ? Matrix{T}(undef,0,0) : inv(symcross(o.Zperp, o.parent.wt))
  o.ZperpinvZperpZperp = o.Zperp * o.invZperpZperp

  o.X₁ = o.parent.X₁ * perp(o.RperpX); o.X₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.X₁)  # FWL-process X₁ XXX XB(o.parent.X₁, perp(RperpX))
  o.X₂ = o.parent.X₂ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.X₂)                # FWL-process X₂
  X₂X₁ = cross(o.X₂, o.parent.wt, o.X₁)
  o.XX = Symmetric([symcross(o.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.X₂, o.parent.wt)])
  o.kX = ncols(o.XX)
  o.invXX = invsym(o.XX)

  o.Z   = X₁₂B(o.parent.X₁, o.parent.Y₂, o.Rpar     )  # Zpar
  o.ZR₁ = X₁₂B(o.parent.X₁, o.parent.Y₂, o.R₁invR₁R₁)

  o.Z   .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.Z  )  # partialling out
  o.ZR₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.ZR₁)
  o.Y₂ = o.parent.Y₂ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.Y₂)
  o.y₁ = o.parent.y₁ - o.ZperpinvZperpZperp * crossvec(o.Zperp, o.parent.wt, o.parent.y₁)

  o.X₁Y₂ = cross(o.X₁, o.parent.wt, o.Y₂)
  o.X₂Y₂ = cross(o.X₂, o.parent.wt, o.Y₂)
  o.XY₂ = [o.X₁Y₂ ; o.X₂Y₂]
  o.Y₂y₁ = crossvec(o.Y₂, o.parent.wt, o.y₁)
  o.X₂y₁ = crossvec(o.X₂, o.parent.wt, o.y₁)
  o.X₁y₁ = crossvec(o.X₁, o.parent.wt, o.y₁)
  o.y₁y₁ = cross(o.y₁ , o.parent.wt, o.y₁)[1]
  o.Zy₁  = crossvec(o.Z, o.parent.wt, o.y₁)
  o.XZ   = [cross(o.X₁, o.parent.wt, o.Z) ;
			      cross(o.X₂, o.parent.wt, o.Z)]
  o.ZY₂ =  cross(o.Z, o.parent.wt, o.Y₂)
  o.ZZ  =  symcross(o.Z, o.parent.wt)

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ

  if length(o.R₁invR₁R₁)>0
	  o.X₂ZR₁    = cross(o.X₂, o.parent.wt, o.ZR₁)
	  o.X₁ZR₁    = cross(o.X₁, o.parent.wt, o.ZR₁)
	  o.ZZR₁     = cross(o.Z , o.parent.wt, o.ZR₁)
	  o.twoR₁Zy₁ = 2crossvec(o.ZR₁, o.parent.wt, o.y₁)
	  o.ZR₁ZR₁   = symcross(o.ZR₁, o.parent.wt)
	  o.ZR₁Y₂    = cross(o.ZR₁, o.parent.wt, o.Y₂)
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
	  o.kZ = ncols(o.Rpar)
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
	  	    uwt = vHadw(view(o.PXZ,:,i), o.parent.wt)
	  	    o.parent.NFE>0 &&
	  	  	(o.CT_FEcapPY[i+1] = crosstabFEt(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)
	  	    tmp = @panelsum(o.Z, uwt, o.parent.infoCapData)
	  	    for j ∈ 1:o.kZ
	  	  	  o.FillingT₀[i+1,j+1] = view(tmp,:,j)
	  	    end
	  	  end

	  	  for i ∈ 1:o.kZ  # precompute various clusterwise sums
	  	    o.ScapPXYZperp[i+1] = @panelsum(o.Zperp, vHadw(view(o.PXZ,:,i), o.parent.wt), o.parent.infoCapData)  # Scap(P_(MZperpX) * Z .* Zperp)
	  	    !o.parent.granular &&
	  	  	(o.ScapYX[i+1] = @panelsum2(o.X₁, o.X₂, vHadw(view(o.Z,:,i), o.parent.wt), o.parent.infoCapData))  # Scap(M_Zperp[Z or y₁] .* P_(MZperpX)])
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
  o.y₁par = o.parent.y₁ - o.parent.Y₂ * (isone(length(r₁)) ? r₁[1] : r₁)
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
    o.y₁pary₁par = o.y₁y₁ - (o.twoR₁Zy₁'r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
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
	    tmp = @panelsum(o.Z, uwt, o.parent.infoCapData)
	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[1,i+1] = view(tmp,:,i)
	    end
	    o.ScapPXYZperp[1] = @panelsum(o.Zperp, uwt, o.parent.infoCapData)  # Scap(P_(MZperpX) * y₁ .* Zperp)

	    o.parent.NFE>0 &&
	  	(o.CT_FEcapPY[1] = crosstabFEt(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)

	    uwt = vHadw(o.y₁par, o.parent.wt)
	    tmp = @panelsum(o.PXZ, uwt, o.parent.infoCapData)
	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[i+1,1] = view(tmp,:,i)
	    end
	    o.FillingT₀[1] = @panelsum(o.PXy₁, uwt, o.parent.infoCapData)
	    !o.parent.granular &&
	  	(o.ScapYX[1] = @panelsum2(o.X₁, o.X₂  , uwt, o.parent.infoCapData))  # Scap(M_Zperp*y₁ .* P_(MZperpX)])
	  end
  end
end


@inline function MakeResiduals!(o::StrEstimator{T,<:Union{OLS,ARubin}} where T)
  o.ü₁ = o.y₁par - X₁₂B(o.parent.X₁, o.parent.X₂, o.β)
end

function MakeResiduals!(o::StrEstimator{T,IVGMM} where T)
  o.ü₁ = o.y₁par - o.Z * o.β

  if !o.parent.scorebs
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
  if o.parent.bootstrapt && (o.parent.scorebs || o.parent.robust)
	  (o.parent.granular || o.parent.purerobust) && (o.WXAR = vHadw(o.XAR, o.parent.wt))

	  if o.parent.robust && o.parent.NFE>0 && !(o.parent.FEboot || o.parent.scorebs) && o.parent.granular < o.parent.NErrClustCombs  # make first factor of second term of (64) for c=cap (c=1)
	    !isdefined(o, :WXAR) && (o.WXAR = vHadw(o.XAR, o.parent.wt))
	    o.CT_XAR = crosstabFEt(o.parent, o.WXAR, o.parent.infoCapData)
	  end
  end
end

# partial fixed effects out of a data matrix
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

function setdirty!(o::StrBootTest, _dirty::Bool; noinitialize::Bool=false)
  o.dirty = _dirty
  !noinitialize && (o.initialized = false)
end

function getdist(o::StrBootTest, diststat::String="")
  o.dirty && boottest!(o)
  if diststat == "numer"
	  _numer = isone(o.v_sd) ? o.numer : o.numer / o.v_sd
	  o.distCDR = (@view _numer[:,2:end]) .+ o.r
	  sort!(o.distCDR)
  elseif length(o.distCDR)==0
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
getNBootClust(o::StrBootTest) = o.Nstar
getreps(o::StrBootTest) = o.B  # return number of replications, possibly reduced to 2^G

function getpadj(o::StrBootTest{T}; classical::Bool=false) where T
  _p = o.dirty || classical ? getp(o, classical=classical) : o.p
  if o.madjtype==bonferroni min(one(T), o.NH0 * _p)
  elseif o.madjtype==sidak  one(T) - (one(T) - _p) ^ o.NH0
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

@inline _sortperm(X::AbstractVecOrMat) = size(X,2)==1 ? sortperm(ndims(X)==1 ? X : reshape(X,length(X)), alg=RadixSort) : sortperm(collect(eachrow(X)))  # sort a data matrix

function Init!(o::StrBootTest{T}) where T  # for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  o.Nobs = nrows(o.X₁)
  o.NClustVar = ncols(o.ID)
  o.kX = (o.kX₁ = ncols(o.X₁)) + (o.kX₂ = ncols(o.X₂))
  o.kX₂==0 && (o.X₂ = zeros(T,o.Nobs,0))
  o.kY₂ = ncols(o.Y₂)
  iszero(o.kY₂) && (o.Y₂ = zeros(T,o.Nobs,0))
  o.kZ = o.kX₁ + o.kY₂
  if o.LIML && o.kX₂==o.kY₂  # exactly identified LIML = 2SLS
  	o.κ = one(T)
  	o.LIML = false
  end
  if !(o.REst = length(o.R₁)>0)  # base model contains no restrictions?
    o.R₁ = zeros(T,0,o.kZ)
    o.r₁ = zeros(T,0)
  end
  isnan(o.κ) && (o.κ = o.kX₂>0 ? one(T) : zero(T))  # if κ in κ-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  o.WRE = !(iszero(o.κ) || o.scorebs) || o.ARubin
  o.WREnonARubin = o.WRE && !o.ARubin

  o.haswt = typeof(o.wt) <: AbstractVector
  if o.haswt
	  o.sumwt = sum(o.wt)
  else
	  o.wt = I
	  o.sumwt = 1
  end
  o._Nobs = o.haswt && o.fweights ? o.sumwt : o.Nobs

  if !iszero(o.NClustVar)
	o.subcluster = o.NClustVar - o.nerrclustvar
	p = _sortperm(view(o.ID, :, [collect(o.subcluster+1:o.NClustVar); collect(1:o.subcluster)]))  # sort by err cluster vars, then remaining boot cluster vars
	o.ID = ndims(o.ID)==1 ? o.ID[p] : o.ID[p,:]
	o.X₁ = ndims(o.X₁)==1 ? o.X₁[p] : o.X₁[p,:]
	o.X₂ = ndims(o.X₂)==1 ? o.X₂[p] : o.X₂[p,:]
	o.y₁ = o.y₁[p]
	o.Y₂ = ndims(o.Y₂)==1 ? o.Y₂[p] : o.Y₂[p,:]
	o.haswt && (o.wt = o.wt[p])
	isdefined(o, :FEID) && length(o.FEID)>0 && (o.FEID = o.FEID[p])
  end

  if o.WREnonARubin
	  if !iszero(o.NClustVar)
	    o.infoBootData, o.IDBootData = panelsetupID(o.ID, collect(1:o.nbootclustvar))
	  else
	    o.infoCapData = o.infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # no clustering, so no collapsing by cluster
	  end
  elseif !iszero(o.NClustVar)
	  o.infoBootData = panelsetup(o.ID, collect(1:min(o.NClustVar,o.nbootclustvar)))  # bootstrap cluster grouping defs rel to original data
  else
	  infoCapData = infoAllData = infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # causes no collapsing of data in panelsum() calls, only multiplying by weights if any
  end
  o.Nstar = nrows(o.infoBootData)

  if o.bootstrapt
    if !iszero(o.NClustVar)
	    minN = Inf; sumN = 0

	    combs = [x & 2^y > 0 for x in 2^o.nerrclustvar-1:-1:1, y in o.nerrclustvar-1:-1:0]  # represent all error clustering combinations. First is intersection of all error clustering vars
	    o.clust = Vector{StrClust{T}}(undef, nrows(combs))  # leave out no-cluster combination
	    o.NErrClustCombs = length(o.clust)

	    if o.NClustVar > o.nbootclustvar  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
		    if o.WREnonARubin && !o.granular
		      o.infoAllData, IDAllData = panelsetupID(o.ID, collect(1:o.NClustVar))
		    else
		      o.infoAllData            = panelsetup(o.ID, collect(1:o.NClustVar))
		    end
	    else
		    o.infoAllData = o.infoBootData  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
		    o.WREnonARubin && !o.granular && (IDAllData = o.IDBootData)
	    end

		Nall = length(o.infoAllData)

	    if o.NClustVar > o.nerrclustvar  # info for intersections of error clustering wrt data
		    if o.WREnonARubin && !o.granular
		      o.infoCapData, IDCapData = panelsetupID(o.ID, collect(o.subcluster+1:o.NClustVar))
		    else
		      o.infoCapData            = panelsetup(o.ID, collect(o.subcluster+1:o.NClustVar))
		    end
		    IDCap = length(o.infoCapData)==o.Nobs ? o.ID : @views o.ID[first.(o.infoCapData),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
		    o.IDAll = Nall==o.Nobs ? o.ID : @view o.ID[first.(o.infoAllData),:]  # version of ID matrix with one row for each all-bootstrap && error cluster-var intersection instead of 1 row for each obs
	    else
		    o.infoCapData = o.infoAllData  # info for intersections of error clustering wrt data
		    o.WREnonARubin && !o.granular && (IDCapData = IDAllData)
		    o.IDAll = IDCap = nrows(o.infoCapData)==o.Nobs ? o.ID : @views o.ID[first.(o.infoCapData),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
	    end

	    o.BootClust = 2^(o.NClustVar - o.nbootclustvar)  # location of bootstrap clustering within list of cluster combinations

	    for c ∈ 1:o.NErrClustCombs  # for each error clustering combination
		    ClustCols = o.subcluster .+ findall(@view combs[c,:])
		    even = isodd(length(ClustCols))  # not a typo

		    if isone(c)
		      if iszero(o.subcluster)
		    	  order = Vector{Int64}(undef,0)
		    	  info  = Vector{UnitRange{Int64}}(undef, Nall)  # causes no collapsing of data in panelsum() calls
		      else
		    	  order = _sortperm(@view IDCap[:,ClustCols])
		    	  IDCap = @view IDCap[order, :]
		    	  info  = panelsetup(IDCap, ClustCols)
		      end
		    else
		      if any(combs[c, min(findall(view(combs,c,:) .≠ view(combs,c-1,:))...):end])  # if this sort ordering same as last to some point and missing thereafter, no need to re-sort
		    	  order = _sortperm(@view IDCap[:,ClustCols])
		    	  IDCap = @view IDCap[order,:]
		      else
		    	  order = Vector{Int64}(undef,0)
		      end
		      info = panelsetup(IDCap, ClustCols)
		    end

		    N = nrows(info)
		    sumN += N

			if o.small
			  	multiplier = T(N / (N-1))
			  	N < minN && (minN = N)
			else
			  	multiplier = one(T)
			end
			o.clust[c] = StrClust{T}(N, multiplier, even, order, info)
	    end

	    (o.scorebs || !o.WREnonARubin) &&
		  	(o.ClustShare = o.haswt ? @panelsum(o.wt, o.infoCapData)/o.sumwt : length.(o.infoCapData)./o.Nobs) # share of observations by group

    else  # if no clustering, cast "robust" as clustering by observation
      o.clust = StrClust{T}(Nobs, small ? _Nobs / (_Nobs - 1) : 1, true, Vector{Int64}(undef,0), Vector{UnitRange{Int64}}(undef,0))
      sumN = o.Nobs
      o.NErrClustCombs = 1
      (o.scorebs || !o.WREnonARubin) &&
    		(o.ClustShare = o.haswt ? o.wt/o.sumwt : 1/o._Nobs)
    end

		o.purerobust = o.robust && !o.scorebs && iszero(o.subcluster) && o.Nstar==o.Nobs  # do we ever error-cluster *and* bootstrap-cluster by individual?
		o.granular   = o.WREnonARubin ? 2*o.Nobs*o.B*(2*o.Nstar+1) < o.Nstar*(o.Nstar*o.Nobs+o.clust[1].N*o.B*(o.Nstar+1)) :
											o.NClustVar>0 && !o.scorebs && (o.purerobust || (o.clust[1].N+o.Nstar)*o.kZ*o.B + (o.clust[1].N-o.Nstar)*o.B + o.kZ*o.B < o.clust[1].N*o.kZ^2 + o.Nobs*o.kZ + o.clust[1].N * o.Nstar * o.kZ + o.clust[1].N * o.Nstar)

		if o.robust && !o.purerobust
			(o.subcluster>0 || o.granular) &&
				(o.infoErrAll = panelsetup(o.IDAll, collect(o.subcluster+1:o.NClustVar)))  # info for error clusters wrt data collapsed to intersections of all bootstrapping && error clusters; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusterings
			((o.scorebs && o.B>0) || (o.WREnonARubin && !o.granular && o.bootstrapt)) &&
				(o.JNcapNstar = zeros(T, o.clust[1].N, o.Nstar))
		end

		if o.WREnonARubin && o.robust && o.bootstrapt && !o.granular
			if iszero(length(IDAllData))
				_, IDAllData = panelsetupID(o.ID, collect(             1:o.NClustVar))
				_, IDCapData = panelsetupID(o.ID, collect(o.subcluster+1:o.NClustVar))
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
	  minN = nrows(o.infoBootData)
		end

		if isdefined(o, :FEID) && length(o.FEID)>0
			p = _sortperm(o.FEID)
			sortID = o.FEID[p]
			i_FE = 1; o.FEboot = o.B>0 && !o.WREnonARubin && o.NClustVar>0; j = o.Nobs; o._FEID = ones(Int64, o.Nobs)
			o.invFEwt = zeros(T, o.NFE>0 ? o.NFE : o.Nobs)
			o.FEs = Vector{StrFE{T}}(undef, o.NFE>0 ? o.NFE : o.Nobs)
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
				o.FEs[i_FE] = StrFE{T}(is, wt)
				if (o.B>0 && o.robust && o.granular < o.nerrclustvar) || (o.WREnonARubin && o.robust && o.granular && o.bootstrapt)
					o.invFEwt[i_FE] = 1 / sumFEwt
				end

			j = i

			if o.FEboot  # are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
				tmp = o.ID[is, 1:o.nbootclustvar]
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
		o.FEs[i_FE] = StrFE{T}(is, wt)
		o.robust && ((o.B>0 && o.granular < o.nerrclustvar) || (o.WREnonARubin && o.granular && o.bootstrapt)) &&
			(o.invFEwt[i_FE] = 1 / sumFEwt)
		o.NFE = i_FE
		resize!(o.invFEwt, o.NFE)
		resize!(o.FEs, o.NFE)
		if o.FEboot  # are all of this FE's obs in same bootstrapping cluster?
			tmp = o.ID[is, 1:o.nbootclustvar]
			o.FEboot = all(tmp .== @view tmp[1,:])
		end

		if o.robust && o.B>0 && o.bootstrapt && !o.FEboot && o.granular < o.nerrclustvar
			o.infoBootAll = panelsetup(o.IDAll, collect(1:o.nbootclustvar))  # info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping && error clusters
		end

		o.X₁ = partialFE(o, o.X₁)  # don't overwrite caller's data
		o.X₂ = partialFE(o, o.X₂)
		o.y₁ = partialFE(o, o.y₁)
		o.Y₂ = partialFE(o, o.Y₂)
	end

	if o.B>0 && o.robust && o.granular && !o.purerobust && o.bootstrapt && !o.WREnonARubin
		if o.NFE>0 && !o.FEboot
			_, o.IDBootData = panelsetupID(o.ID   , collect(1:o.nbootclustvar))
		else
			_, o.IDBootAll  = panelsetupID(o.IDAll, collect(1:o.nbootclustvar))
		end
	end

	o.enumerate = o.B>0 && o.auxtwtype==rademacher && o.Nstar*log(2) < log(o.B)+1e-6  # generate full Rademacher set?
	o.enumerate && (o.maxmatsize = 0)

	o.Nw = iszero(o.maxmatsize) ? 1 : ceil((o.B+1) * Float64(max(nrows(o.IDBootData), length(o.IDBootAll), o.Nstar) * sizeof(T)) / o.maxmatsize / 1073741824) # 1073741824 = giga(byte)
	if isone(o.Nw)
		MakeWildWeights!(o, o.B, first=true)  # make all wild weights, once
		o.enumerate && (o.B = ncols(o.v) - 1)  # replications reduced to 2^G
		o.WeightGrp = [1:ncols(o.v)]
	else
		_B = ceil(Int64, (o.B+1) / o.Nw)
		o.Nw = ceil(Int64, (o.B+1) / _B)
		o.WeightGrp = [(i-1)*_B+1:i*_B for i ∈ 1:o.Nw]
		o.WeightGrp[end] = first(o.WeightGrp[end]):o.B+1
	end

	if o.ML
		o.dof = nrows(o.R)
	else
		if o.ARubin
			o.R = hcat(zeros(o.kX₂,o.kX₁), Matrix(I(o.kX₂)))  # attack surface is all endog vars
			o.R₁ = o.kX₁>0 && nrows(o.R₁)>0 ? hcat(o.R₁[:,1:kX₁], zeros(nrows(o.R₁),o.kX₂)) : zeros(0, o.kX)  # and convert model constraints from referring to X₁, Y₂ to X₁, X₂
		end
		o.dof = nrows(o.R)

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
			o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kX, ncols(o.v))] :
									[Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[end]))]
			o.T1R = deepcopy(o.T1L)

			if o.bootstrapt
				o.δdenom_b = zeros(o.Repl.kZ, o.Repl.kZ)
				o.SstarUMZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
				o.SstarUPX     = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
				o._Jcap = zeros(o.clust[1].N, o.Repl.kZ)
				!o.granular && (o.SCTcapuXinvXX = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Nstar))
				if o.LIML || !o.robust
					o.YYstar_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
					o.YPXYstar_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
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
			o.M = o.Repl  # StrEstimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scorebs for IV/GMM mixes the two
			if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
				Estimate!(o.DGP, o.r₁)
				MakeResiduals!(o.DGP)
			end
		end
  end

  if !o.WREnonARubin && o.bootstrapt
		o.denom = [Matrix{T}(undef,0,0) for _ in 1:o.dof, _ in 1:o.dof]
		if o.robust
			o.Kcd =                       Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof)
			o.Jcd = iszero(o.B) ? o.Kcd : Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof)  # if B = 0, Kcd will be multiplied by v, which is all 1's, and will constitute Jcd
		end

		if o.robust && o.granular < o.NErrClustCombs && o.B>0
			if o.subcluster>0  # crosstab c,c* is wide
				o.crosstabCapstarind = Vector{Int64}(undef, Nall)
				for i ∈ axes(o.infoErrAll, 1)
					o.crosstabCapstarind[o.infoErrAll[i]]              = collect(i                            .+ o.clust[1].N * o.dof * (o.infoErrAll[i] .- 1))
				end
			elseif o.NClustVar == o.nbootclustvar  # crosstab c,c* is square
				o.crosstabCapstarind = collect(1 : 1+Nall*o.dof : 1+(Nall-1)*(1+Nall*o.dof))  # ~diagind()
			else  # crosstab c,c* is tall
				o.crosstabCapstarind = Vector{Int64}(undef, Nall)
				for i ∈ axes(o.clust[o.BootClust].info, 1)
					o.crosstabCapstarind[o.clust[o.BootClust].info[i]] = collect(o.clust[o.BootClust].info[i] .+ o.clust[1].N * o.dof * (              i  - 1))
				end
			end
			o.crosstabCapstarind = vcat([o.crosstabCapstarind .+ i*o.clust[1].N for i in 0:o.dof-1]...)
		end
  end

  o.small && (o.dof_r = o.NClustVar>0 ? minN - 1 : o._Nobs - o.kZ - o.NFE)

  o.sqrt = isone(o.dof)  # work with t/z stats instead of F/chi2

  if o.small
		o.multiplier = (o.smallsample = (o._Nobs - o.kZ - o.FEdfadj * o.NFE) / (o._Nobs - o.robust)) / o.dof  # divide by # of constraints because F stat is so defined
  else
		o.multiplier = o.smallsample = 1
  end

  !(o.robust || o.ML) && (o.multiplier *= o._Nobs)  # will turn sum of squared errors in denom of t/z into mean
  o.sqrt && (o.multiplier = √o.multiplier)

  ((!o.bootstrapt && o.dof==1) || o.bootstrapt && (o.WREnonARubin || o.dof>1+o.robust || !isnan(o.maxmatsize))) &&  # unless nonWRE or dof=1 or splitting weight matrix, code will create dist element-by-element, so pre-allocate vector now
	(o.dist = fill(T(NaN), o.B+1))
  (o.Nw>1 || o.WREnonARubin || (!o.null && o.dof<=2)) && (o.numer = fill(T(NaN), o.dof, o.B+1))

  if !o.WREnonARubin
		o.poles = o.anchor = zeros(T,0)
		o.interpolable = o.bootstrapt && o.B>0 && o.null && o.Nw==1 && (iszero(o.κ) || o.ARubin)
		if o.interpolable
			o.∂numer∂r = Vector{Matrix{T}}(undef, o.q)
			o.interpolate_u = !(o.robust || o.ML)
			o.interpolate_u && (o.∂u∂r = Vector{Matrix{T}}(undef, o.q))
			if o.robust
				o.∂denom∂r   = [Matrix{Matrix{T}}(undef, o.dof, o.dof) for _ in 1:o.q]
				o.∂²denom∂r² = [Matrix{Matrix{T}}(undef, o.dof, o.dof) for _ in 1:o.q, _ in 1:o.q]
				o.∂Jcd∂r     = [Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof) for _ in 1:o.q]
			end
		end
  end
end

# main routine
function boottest!(o::StrBootTest{T}) where T
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

  o.BFeas = isnan(o.dist[1]) ? 0 : sum(.!(isnan.(o.dist) .| isinf.(o.dist))) - 1
  o.distCDR = zeros(T,0,0)
  o.dirty = false
  o.initialized = true
end

# if not imposing null and we have returned to boottest!(), then dof=1 or 2; we're plotting or finding CI, and only test stat, not distribution, changes with r
function NoNullUpdate!(o::StrBootTest{T} where T)
  if o.WREnonARubin
		o.numer[:,1] = o.R * o.DGP.Rpar * o.βs[1] - o.r
  elseif o.ARubin
		o.DGP.Estimate(o.r)
		o.numer[:,1] = o.v_sd * o.DGP.Rpar * @view o.DGP.β[o.kX₁+1:end,:] # coefficients on excluded instruments in ARubin OLS
  else
		o.numer[:,1] = o.v_sd * (o.R * (o.ML ? o.β : o.M.Rpar * o.M.β) - o.r) # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this
  end
  o.dist[1] = isone(o.dof) ? o.numer[1] / sqrt(o.statDenom[1]) : o.numer[:,1]'invsym(o.statDenom)*o.numer[:,1]
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
			o.dist[:,:] .= o.numer' ./ sqrtNaN.(o.statDenom)
		else
			o.dist = colquadform(invsym(o.statDenom), o.numer)
		end
  end
end

# draw wild weight matrix of width _B. If first=true, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all
const ϕ = (1 + √5)/2

function MakeWildWeights!(o::StrBootTest{T}, _B::Integer; first::Bool=true) where T
  if _B>0  # in scoretest or waldtest WRE, still make v a col of 1's
    if o.enumerate
			o.v = o.WREnonARubin ? [zeros(o.Nstar) count_binary(o.Nstar, -2, 0)] :  # complete Rademacher set
								[ones( o.Nstar) count_binary(o.Nstar, -1, 1)]
		elseif o.auxtwtype==normal
			o.v = randn(o.rng, T, o.Nstar, _B+first)
			o.WREnonARubin && (o.v .-= one(T))
		elseif o.auxtwtype==gamma
			tmp = quantile.(Gamma(4,.5), rand(o.rng, o.Nstar, _B+first))
			o.v = T==Float64 ? tmp : T.(tmp)
			o.WREnonARubin && (o.v .-= one(T))
		elseif o.auxtwtype==webb
			o.v = rand(o.rng, T.([-√1.5, -1, -√.5, √.5, 1, √1.5] .- o.WREnonARubin), o.Nstar, _B+first)
		elseif o.auxtwtype == mammen
			o.v = getindex.(Ref(T.([1-ϕ; ϕ] .- o.WREnonARubin)), ceil.(Int16, rand(o.rng, o.Nstar, _B+first) ./ (ϕ/√5)))
		elseif o.WREnonARubin  # Rademacher
			o.v = -2rand(o.rng, Bool, o.Nstar, _B+first)
		else
			o.v = rand(o.rng, Bool, o.Nstar, _B+first) .- T(.5)  # rand(o.rng, Bool, o.Nstar, _B+first) .- T(.5)
			o.v_sd = .5
		end

		first && !(o.enumerate && isone(o.v_sd)) && (o.v[:,1] .= o.WREnonARubin ? zero(T) : o.v_sd)  # keep original residuals in first entry to compute base model stat
  else
		o.v = Matrix{T}(undef,0,1)  # in places, ncols(v) indicates B -- 1 for classical tests
  end
end


# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of Y[:,ind1]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,ind2] for κ constant across replications
# ind1 can be a rowvector
# (only really the Hessian when we narrow Y to Z)
function HessianFixedκ(o::StrBootTest{T}, ind1::Vector{S} where S<:Integer, ind2::Integer, κ::Number, w::Integer) where T
  dest = Matrix{T}(undef, length(ind1), ncols(o.v))
  for i ∈ eachindex(ind1, axes(dest,1))
		_HessianFixedκ!(o, view(dest,i:i,:), ind1[i], ind2, κ, w)
  end
  dest
end

function _HessianFixedκ!(o::StrBootTest, dest::AbstractMatrix, ind1::Integer, ind2::Integer, κ::Number, w::Integer)
  if !(o.Repl.Yendog[ind1+1] || o.Repl.Yendog[ind2+1])  # if both vars exog, result = order-0 term only, same for all draws
		!iszero(κ) && (@views coldot!(dest, o.Repl.XZ[:,ind1], o.Repl.invXXXZ[:,ind2]))
		if !isone(κ)
			if iszero(κ)
				dest = o.Repl.YY[ind1+1,ind2+1]
			else
				dest .= κ .* dest .+ (1 - κ) .* o.Repl.YY[ind1+1,ind2+1]
			end
		end
		dest .= dest
	else
		if !iszero(κ)
			_T1L = iszero(ind1) ? o.Repl.Xy₁par : @view o.Repl.XZ[:,ind1]
			if o.Repl.Yendog[ind1+1]
				T1L = o.T1L[isone(o.Nw) || w<o.Nw ? 1 : 2]
				T1L .= o.SstarUX[ind1+1] * o.v
				T1L .+= _T1L
			else
				T1L = _T1L
			end
			_T1R = iszero(ind2) ? o.Repl.invXXXy₁par : @view o.Repl.invXXXZ[:,ind2]
			if o.Repl.Yendog[ind2+1]
				T1R = o.T1R[isone(o.Nw) || w<o.Nw ? 1 : 2]
				T1R .= o.SstarUXinvXX[ind2+1] * o.v
				T1R .+= _T1R
			else
				T1R = _T1R  # right-side linear term
			end
			coldot!(dest, T1L, T1R)
		end
		if !isone(κ)
			if o.Repl.Yendog[ind1+1]
				T2 = o.SstarUZperpinvZperpZperp[ind1+1]'o.SstarUZperp[ind2+1]  # quadratic term
				T2[diagind(T2)] .-= ind1 <= ind2 ? o.SstarUU[ind2+1, ind1+1] : o.SstarUU[ind1+1, ind2+1]  # minus diagonal crosstab
				o.NFE>0 &&
					(T2 .+= o.CTFEU[ind1+1]'(o.invFEwt .* o.CTFEU[ind2+1]))

				if iszero(κ)
					dest .= o.Repl.YY[ind1+1,ind2+1] .+ colquadformminus!((                            @views o.SstaruY[ind2+1][:,ind1+1] .+ o.SstaruY[ind1+1][:,ind2+1])'o.v, T2, o.v)
				else
					dest .=   κ .* dest .+ (1 - κ)   .* colquadformminus!((o.Repl.YY[ind1+1,ind2+1] .+ @views o.SstaruY[ind2+1][:,ind1+1] .+ o.SstaruY[ind1+1][:,ind2+1])'o.v, T2, o.v)
				end
			elseif iszero(κ)
				dest .= o.Repl.YY[ind1+1,ind2+1]
			else
				dest .= κ .* dest .+ (1 - κ) .* o.Repl.YY[ind1+1,ind2+1]
			end
		end
  end
end


# Workhorse for WRE CRVE sandwich filling
# With reference to notional Y = [y₁ Z], given 0-based columns index within it, ind1, and a matrix βs of all the bootstrap estimates, return all bootstrap realizations of P_X * Y[:,ind1]_g ' u\hat_1g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per cap cluster, one col per bootstrap replication
function Filling(o::StrBootTest{T}, ind1::Integer, βs::AbstractMatrix) where T
  	if o.granular
    	if o.Nw == 1  # create or avoid NxB matrix?
			PXYstar = iszero(ind1) ? o.Repl.PXy₁ : o.Repl.PXZ[:,ind1]
			o.Repl.Yendog[ind1+1] && (PXYstar = PXYstar .+ o.SstarUPX[ind1+1] * o.v)

			dest = @panelsum(PXYstar .* (o.Repl.y₁ .- o.SstarUMZperp[1] * o.v), o.wt, o.infoCapData)

			for ind2 ∈ 1:o.Repl.kZ
				_β = view(βs,ind2,:)'
				dest .-= @panelsum(PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[:,ind2] * _β .- o.SstarUMZperp[ind2+1] * (o.v .* _β) :
															           o.Repl.Z[:,ind2] * _β                                            ), o.wt, o.infoCapData)
			end
		else  # create pieces of each N x B matrix one at a time rather than whole thing at once
			dest = Matrix{T}(undef, o.clust[1].N, ncols(o.v))  # XXX preallocate this & turn Filling into Filling! ?
			for ind2 ∈ 0:o.Repl.kZ
				ind2>0 && (βv = o.v .* (_β = view(βs,ind2,:)'))

				if o.purerobust
					for i ∈ 1:o.clust[1].N
						PXYstar = hcat(iszero(ind1) ? o.Repl.PXy₁[i] : o.Repl.PXZ[i,ind1])
						o.Repl.Yendog[ind1+1] && (PXYstar = PXYstar .+ view(o.SstarUPX[ind1+1],i,:)'o.v)

						if iszero(ind2)
							dest[i,:]   = wtsum(o.wt, PXYstar .* (o.Repl.y₁[i] .- view(o.SstarUMZperp[1],i,:))'o.v)
						else
							dest[i,:] .-= wtsum(o.wt, PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[i,ind2] * _β .- view(o.SstarUMZperp[ind2+1],i,:)'βv :
																						  o.Repl.Z[i,ind2] * _β))
						end
					end
				else
					for i ∈ 1:o.clust[1].N
						S = o.infoCapData[i]
						PXYstar = ind1 ? o.Repl.PXZ[S,ind1] : o.Repl.PXy₁[S]
						o.Repl.Yendog[ind1+1] && (PXYstar = PXYstar .+ view(o.SstarUPX[ind1+1],S,:) * o.v)

						if iszero(ind2)
							dest[i,:]   = wtsum(o.wt, PXYstar .* (o.Repl.y₁[S] .- view(o.SstarUMZperp[1],S,:) * o.v))
						else
							dest[i,:] .-= wtsum(o.wt, PXYstar .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[S,ind2] * _β .- view(o.SstarUMZperp[ind2+1],S,:) * βv :
																						  o.Repl.Z[S,ind2] * _β                                       ))
						end
					end
				end
			end
    	end
  	else  # coarse error clustering
		for ind2 ∈ 0:o.Repl.kZ
			βv = iszero(ind2) ? o.v : o.v .* (_β = -view(βs,ind2,:)')

			# T1 * o.v will be 1st-order terms
			T₁ = o.Repl.Yendog[ind1+1] ? o.Repl.ScapYX[ind2+1] * o.SstarUXinvXX[ind1+1] : Matrix{T}(undef,0,0)  #  S_∩ (Y_(∥j):*X_∥ ) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

			if o.Repl.Yendog[ind2+1]  # add CT_(cap,*) (P_(X_par ) Y_(pari).*U ̈_(parj) )
				if o.NClustVar == o.nbootclustvar && iszero(o.subcluster)  # simple case of one clustering: full crosstab is diagonal
					tmp = ind1>0 ? o.Repl.XZ[:,ind1] : o.Repl.Xy₁par
					if length(T₁)>0
						T₁[diagind(T₁)] .+= o.SstarUXinvXX[ind2+1]'tmp
					else
						T₁                = o.SstarUXinvXX[ind2+1]'tmp  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JNcapNstar)
					for i ∈ 1:o.Nstar
						T₁[o.IDCTCapstar[i], i] .+= o.SCTcapuXinvXX[ind2+1,i] * (ind1>0 ? view(o.Repl.XZ,:,ind1) : o.Repl.Xy₁par)
					end
				end
				length(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.Repl.ScapPXYZperp[ind1+1] * o.SstarUZperpinvZperpZperp[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
					o.NFE              > 0 && (T₁ = T₁ .- o.Repl.CT_FEcapPY[ind1+1] * o.CTFEU[ind2+1])
			end

			if ind2>0  # order-0 and -1 terms
				if iszero(ncols(T₁))  # - x*β components
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β
				elseif isone(ncols(T₁))
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β .+ T₁ .* βv
				else
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β
					matmulplus!(dest, T₁, βv)
				end
			else  # y component
				if iszero(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1]
				elseif isone(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1] .+ T₁ .* βv
				else
					dest = T₁ * o.v; dest .+= o.Repl.FillingT₀[ind1+1,1]
				end
			end

			if o.Repl.Yendog[ind1+1] && o.Repl.Yendog[ind2+1]
				for i ∈ 1:o.clust[1].N
					S = o.infoCapData[i]
					colquadformminus!(dest, i, cross(view(o.SstarUPX[ind1+1],S,:), o.haswt ? o.wt[S] : I, view(o.SstarUMZperp[ind2+1],S,:)), o.v, βv)
				end
			end
		end
  end
  dest
end

function PrepWRE!(o::StrBootTest)
  Estimate!(o.DGP, o.null ? [o.r₁ ; o.r] : o.r₁)
  MakeResiduals!(o.DGP)
  o.Ü₂par = o.DGP.Ü₂ * o.Repl.RparY  # XXX XB(o.DGP.Ü₂, o.Repl.RparY)

  for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		uwt = vHadw(i>0 ? o.Ü₂par[:,i] : o.DGP.u⃛₁, o.wt)

		# S_star(u .* X), S_star(u .* Zperp) for residuals u for each endog var; store transposed
		o.SstarUX[i+1]      = @panelsum2(o.Repl.X₁, o.Repl.X₂, uwt, o.infoBootData)'
		o.SstarUXinvXX[i+1] = o.Repl.invXX * o.SstarUX[i+1]

		if o.LIML || o.bootstrapt || !isone(o.κ)
			o.SstarUZperp[i+1]              = @panelsum(o.Repl.Zperp, uwt, o.infoBootData)'
			o.SstarUZperpinvZperpZperp[i+1] = o.Repl.invZperpZperp * o.SstarUZperp[i+1]
			o.NFE>0 && (o.CTFEU[i+1] = crosstabFE(o, uwt, o.infoBootData))
		end

		if o.LIML || !o.robust || !isone(o.κ)
			o.SstaruY[i+1] = @panelsum2(o.Repl.y₁par, o.Repl.Z, uwt, o.infoBootData)
			for j ∈ 0:i
				o.SstarUU[i+1,j+1] = @panelsum(j>0 ? o.Ü₂par[:,j] : o.DGP.u⃛₁, uwt, o.infoBootData)
			end
		end

		if o.robust && o.bootstrapt
			if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				for g ∈ 1:o.Nstar
					o.SCTcapuXinvXX[i+1,g] = @panelsum(o.Repl.XinvXX, uwt, o.infoCTCapstar[g])
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
				(o.SstarUMZperp[i+1] .+= (o.invFEwt .* view(o.CTFEU[i+1]),o._FEID,:))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
		end
  end
end

function MakeWREStats!(o::StrBootTest{T}, w::Integer) where T
    if isone(o.Repl.kZ)  # optimized code for 1 coefficient in bootstrap regression
		if o.LIML
			YY₁₁   = HessianFixedκ(o, [0], 0, zero(T), w)  # κ=0 => Y*MZperp*Y
			YY₁₂   = HessianFixedκ(o, [0], 1, zero(T), w)
			YY₂₂   = HessianFixedκ(o, [1], 1, zero(T), w)
			YPXY₁₁ = HessianFixedκ(o, [0], 0, one(T) , w)  # κ=1 => Y*PXpar*Y
			YPXY₁₂ = HessianFixedκ(o, [0], 1, one(T) , w)
			YPXY₂₂ = HessianFixedκ(o, [1], 1, one(T) , w)
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
			o.As = HessianFixedκ(o, [1], 1, o.κ, w)
			βs = HessianFixedκ(o, [1], 0, o.κ, w) ./ o.As
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
					length(o.clust[c].order)>0 && 
						(Jcaps = view(Jcaps,o.clust[c].order,:))
					@clustAccum!(denom, c, coldot(@panelsum(Jcaps, o.clust[c].info)))
				end
			else
				denom = (HessianFixedκ(o, [0], 0, zero(T), w) .- 2 .* βs .* HessianFixedκ(o, [0], 1, zero(T), w) .+ βs.^2 .* HessianFixedκ(o, [1], 1, zero(T), w)) ./ o._Nobs ./ o.As  # classical error variance
			end
			@storeWtGrpResults!(o.dist, view(o.sqrt ? o.numerw ./ sqrt.(denom) : o.numerw .^ 2 ./ denom, 1, :))
			denom *= o.Repl.RRpar[1]^2
		end

  else  # WRE bootstrap for more than 1 coefficeint in bootstrap regression

		βs = zeros(T, o.Repl.kZ, ncols(o.v))
		A = Vector{Matrix{T}}(undef, ncols(o.v))

		if o.LIML
			o.YYstar   = [HessianFixedκ(o, collect(0:i), i, zero(T), w) for i ∈ 0:o.Repl.kZ] # κ=0 => Y*MZperp*Y
			o.YPXYstar = [HessianFixedκ(o, collect(0:i), i,  one(T), w) for i ∈ 0:o.Repl.kZ] # κ=1 => Y*PXpar*Y

			for b ∈ axes(o.v,2)
				for i ∈ 0:o.Repl.kZ
					o.YYstar_b[1:i+1,i+1]   = o.YYstar[i+1][:,b]  # fill uppper triangles, which is all that invsym() looks at
					o.YPXYstar_b[1:i+1,i+1] = o.YPXYstar[i+1][:,b]
				end
				o.κ = 1/(1 - eigvals(invsym(o.YYstar_b) * Symmetric(o.YPXYstar_b))[1])
				!iszero(o.Fuller) && (o.κ -= o.Fuller / (o._Nobs - o.kX))
				βs[:,b] = (A[b] = invsym(o.κ*o.YPXYstar_b[2:end,2:end] + (1-o.κ)*o.YYstar_b[2:end,2:end])) * (o.κ*o.YPXYstar_b[1,2:end]' + (1-o.κ)*o.YYstar_b[1,2:end]')
			end
		else
			δnumer = HessianFixedκ(o, collect(1:o.Repl.kZ), 0, o.κ, w)
			δdenom = [HessianFixedκ(o, collect(1:i), i, o.κ, w) for i ∈ 1:o.Repl.kZ]
			
			for b ∈ axes(o.v,2)
				for i ∈ 1:o.Repl.kZ
					o.δdenom_b[1:i,i] = view(δdenom[i],:,b)  # fill uppper triangle, which is all that invsym() looks at
				end
				βs[:,b] = (A[b] = invsym(o.δdenom_b)) * view(δnumer,:,b)
			end
		end

		if o.bootstrapt
			if o.robust
				Zyg = [Filling(o, i, βs) for i ∈ 1:o.Repl.kZ]  # XXX concatenate into 3-d array?
			else
				o.YYstar = [HessianFixedκ(o, collect(i:o.Repl.kZ), i, zero(T), w) for i ∈ 0:o.Repl.kZ]  # κ=0 => Y*MZperp*Y
			end
		end

		numer_b = Vector{T}(undef,nrows(o.Repl.RRpar))  # XXX move to Init()
		for b ∈ reverse(axes(o.v,2))
			if o.null || w==1 && b==1
				numer_b .= o.Repl.RRpar * βs[:,b] + o.Repl.Rt₁ - o.r
			else
				numer_b .= o.Repl.RRpar * (βs[:,b] - o.DGP.β₀)
			end

			if o.bootstrapt
				if o.robust  # Compute denominator for this WRE test stat
					for i ∈ 1:o.Repl.kZ  # XXX replace with 3-D array
						o._Jcap[:,i] = view(Zyg[i],:,b)
					end
					Jcap = o._Jcap * (A[b] * o.Repl.RRpar')

					for c ∈ 1:o.NErrClustCombs
						(!isone(o.NClustVar) && length(o.clust[c].order)>0) && (Jcap = view(Jcap,o.clust[c].order,:))
						J_b = @panelsum(Jcap, o.clust[c].info)
						@clustAccum!(denom, c, J_b'J_b)
					end
				else  # non-robust
					for i ∈ 0:o.Repl.kZ
						YYstar_b[i+1,i+1:o.Repl.kZ+1] = view(YYstar[i+1],b,:)  # fill upper triangle
					end
					denom = (o.Repl.RRpar * A[b] * o.Repl.RRpar') * [-one(T) ; βs[:,b]]'Symmetric(YYstar_b) * [-one(T) ; βs[:,b]] / o._Nobs  # 2nd half is sig2 of errors
				end
				if o.sqrt
					o.dist[b+first(o.WeightGrp[w])-1] = numer_b[1] / sqrt(denom[1])
				else
					o.dist[b+first(o.WeightGrp[w])-1] = numer_b'invsym(denom)*numer_b  # hand-code for 2-dimensional?
				end
			end
			o.numer[:,b+first(o.WeightGrp[w])-1] = numer_b  # slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
		end
	end
	w==1 && o.bootstrapt && (o.statDenom = denom)  # original-sample denominator
end


# Construct stuff that depends linearly or quadratically on r, possibly by interpolation
function MakeInterpolables!(o::StrBootTest{T} where T)
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
					if o.robust  # dof > 1 for an ARubin test with >1 instruments.
						for d₁ ∈ 1:o.dof
							for c ∈ 1:o.NErrClustCombs
								o.∂Jcd∂r[h₁][c,d₁] = (o.Jcd[c,d₁] - o.Jcd₀[c,d₁]) / o.poles[h₁]
								for d₂ ∈ 1:d₁
									tmp = coldot(o.Jcd₀[c,d₁], o.∂Jcd∂r[h₁][c,d₂])
									d₁ ≠ d₂ && (coldotplus!(tmp, o.Jcd₀[c,d₂], o.∂Jcd∂r[h₁][c,d₁]))  # for diagonal items, faster to just double after the c loop
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
						for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁, c ∈ 1:o.NErrClustCombs
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
				for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁
					o.denom[d₁,d₂] = o.denom₀[d₁,d₂] .+ o.∂denom∂r[d₁,d₂][1,1] .* Δ .+ o.∂²denom∂r²[d₁,d₂][1,1] .* Δ.^2
				end
			else  # q==2
				for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁
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
function _MakeInterpolables!(o::StrBootTest{T}, thisr::AbstractVector) where T
	if o.ML
		o.uXAR = o.sc * (o.AR = o.A * o.R')
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

		(o.scorebs || (o.robust && o.granular < o.NErrClustCombs)) &&
			(o.uXAR = o.DGP.ü₁ .* o.M.XAR)
	end

	o.SuwtXA = o.scorebs ?
				 o.B>0 ?
					 o.NClustVar ?
				          @panelsum(o.uXAR, o.wt, o.infoBootData) :
					      vHadw(o.uXAR, o.wt)                    :
				        wtsum(o.wt, o.uXAR)                      :
			  o.DGP.A * @panelsum2(o.X₁, o.X₂, vHadw(o.ü, o.wt), o.infoBootData)'  # same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X₂ empty except in Anderson-Rubin

	if o.robust && o.bootstrapt && o.granular < o.NErrClustCombs
		ustarXAR = @panelsum(o.uXAR, o.wt, o.infoAllData)  # collapse data to all-boot && error-cluster-var intersections. If no collapsing needed, panelsum() will still fold in any weights
		if o.B>0
			if o.scorebs
				Kd = zeros(T, o.clust[1].N, o.dof, o.Nstar)  # inefficient, but not optimizing for the score bootstrap
			else
				Kd = @panelsum2(o.X₁, o.X₂, vHadw(o.DGP.XAR, o.wt), o.infoCapData) * o.SuwtXA  # overloaded def of * for >2D arrays
			end

			o.NFE>0 && !o.FEboot && (o.CT_WE = crosstabFE(o, vHadw(o.ü, o.wt), o.infoBootData))

			o.NFE>0 && !o.FEboot &&
				(Kd .+= o.M.CT_XAR * (o.invFEwt .* o.CT_WE))  # overloaded def of * for >2D arrays
			Kd[o.crosstabCapstarind] .-= reshape(ustarXAR,:)  # subtract crosstab of ustarXAR wrt bootstrapping cluster and all-cluster-var intersections from M
			o.scorebs && (Kd .-= o.ClustShare * colsum(Kd))  # recenter

			for c ∈ 1+o.granular:o.NErrClustCombs  # XXX pre-compute common iterators
				length(o.clust[c].order)>0 &&
					(Kd = view(Kd, o.clust[c].order,:,:))
				for d ∈ 1:o.dof
					o.Kcd[c,d] = @panelsum(view(Kd,:,d,:), o.clust[c].info)
				end
			end
		else  # B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
			o.scorebs &&
				(ustarXAR .-= o.ClustShare * colsum(ustarXAR))  # recenter if OLS
			for c ∈ 1:o.NErrClustCombs
				length(o.clust[c].order)>0 &&
					(ustarXAR = view(ustarXAR, o.clust[c].order,:))
				tmp = @panelsum(ustarXAR, o.clust[c].info)
				for d ∈ 1:o.dof
					o.Kcd[c,d] = reshape(view(tmp,:,d),:,1)
				end
			end
		end
	end
	MakeNumerAndJ!(o, 1, thisr)  # compute J = κ * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
end

# compute stuff depending linearly on v, needed to prep for interpolation
function MakeNumerAndJ!(o::StrBootTest{T}, w::Integer, r::AbstractVector=Vector{T}(undef,0)) where T  # called to *prepare* interpolation, or when w>1, in which case there is no interpolation
	o.numerw = o.scorebs ?
			   (o.B>0 ?
				 o.SuwtXA'o.v :
				 o.SuwtXA * o.v_sd    ) :
			   (!o.robust || o.granular || o.purerobust ?
				  o.R * (o.βdev = o.SuwtXA * o.v) :
				 (o.R * o.SuwtXA) * o.v)

	if isone(w)
		if o.ARubin
			o.numerw[:,1] = o.v_sd * o.DGP.Rpar * o.DGP.β[o.kX₁+1:end]  # coefficients on excluded instruments in ARubin OLS
		elseif !o.null
			o.numerw[:,1] = o.v_sd * (o.R * (o.ML ? o.β : o.M.Rpar * o.M.β) - r)  # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this.
		end
	end

	@storeWtGrpResults!(o.numer, o.numerw)

	@views if o.B>0 && o.robust && o.bootstrapt
		if o.granular || o.purerobust  # optimized treatment when bootstrapping by many/small groups
			if o.purerobust
				o.ustar = o.ü .* o.v
				partialFE!(o, o.ustar)
				o.ustar -= X₁₂B(o.X₁, o.X₂, o.βdev)  # XXX make X₁₂Bminus
			else  # clusters small but not all singletons
				if o.NFE>0 && !o.FEboot
					o.ustar = o.ü .* view(o.v, o.IDBootData, :)
					partialFE!(o, o.ustar)
					for d ∈ 1:o.dof
						o.Jcd[1,d] = @panelsum(o.ustar, o.M.WXAR[:,d], o.infoCapData)                           - @panelsum2(o.X₁, o.X₂, o.M.WXAR[:,d], o.infoCapData) * o.βdev
					end
				else
					_v = view(o.v,o.IDBootAll,:)
					for d ∈ 1:o.dof
						o.Jcd[1,d] = panelsum( panelsum(o.ü, o.M.WXAR[:,d], o.infoAllData) .* _v, o.infoErrAll) - @panelsum2(o.X₁, o.X₂, o.M.WXAR[:,d], o.infoCapData) * o.βdev
					end
				end
			end
		end
		for c ∈ o.granular+1:o.NErrClustCombs, d ∈ eachindex(axes(o.Jcd, 2), axes(o.Kcd, 2))
			o.Jcd[c,d] = o.Kcd[c,d] * o.v
		end
	end
end

function MakeNonWREStats!(o::StrBootTest{T}, w::Integer) where T
	w > 1 && MakeNumerAndJ!(o, w)
	!o.bootstrapt && return

	if o.robust
    	if !o.interpolating  # these quadratic computation needed to *prepare* for interpolation but are superseded by interpolation once it is going
      		o.purerobust && (ustar2 = o.ustar .^ 2)
      		for i ∈ 1:o.dof, j ∈ 1:i
    			o.purerobust &&
  	      			(o.denom[i,j] = cross(view(o.M.WXAR,:,i), view(o.M.WXAR,:,j), ustar2) * (o.clust[1].even ? o.clust[1].multiplier : -o.clust[1].multiplier))
				for c ∈ o.purerobust+1:o.NErrClustCombs
					@clustAccum!(o.denom[i,j], c, j==i ? coldot(o.Jcd[c,i]) : coldot(o.Jcd[c,i],o.Jcd[c,j]))
				end
  	 		 end
   		end

		if isone(o.dof)
			@storeWtGrpResults!(o.dist, vec(o.numerw ./ sqrtNaN.(o.denom[1,1])))
			isone(w) &&
				(o.statDenom = hcat(o.denom[1,1][1]))  # original-sample denominator
		elseif o.dof==2  # hand-code 2D numer'inv(denom)*numer
			t1 = view(o.numerw,1,:)'; t2 = view(o.numerw,2,:)'; t12 = t1.*t2
			@storeWtGrpResults!(o.dist, vec((t1.^2 .* o.denom[2,2] .- 2 .* t12 .* o.denom[2,1] .+ t2.^2 .* o.denom[1,1]) ./ (o.denom[1,1].*o.denom[2,2] .- o.denom[2,1].^2)))
			isone(w) &&
				(o.statDenom = [o.denom[1,1][1] o.denom[2,1][1] ; o.denom[2,1][1] o.denom[2,2][1]])  # original-sample denominator
		else  # build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = Matrix{T}(undef, o.dof, o.dof)
			for k ∈ 1:ncols(o.v)  # XXX probably can simplify
				for i ∈ 1:o.dof, j ∈ 1:i
					tmp[j,i] = o.denom[i,j][k]  # fill upper triangle, which is all invsym() looks at
				end
				numer_l = view(o.numerw,:,k)
				o.dist[k+first(o.WeightGrp[w])-1] = numer_l'invsym(tmp)*numer_l  # in degenerate cases, cross() would turn cross(.,.) into 0
			end
			isone(w) && (o.statDenom = tmp)  # original-sample denominator
		end
	else  # non-robust
		AR = o.ML ? o.AR : o.M.AR
		if isone(o.dof)  # optimize for one null constraint
			o.denom[1,1] = o.R * AR
			if !o.ML
				o.ustar = o.B>0 ? o.v .* o.ü : o.ü
				if o.scorebs
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

			@storeWtGrpResults!(o.dist, vec(o.numerw ./ sqrtNaN.(o.denom[1,1])))
			isone(w) && (o.statDenom = o.denom[1,1])  # original-sample denominator
		else
			o.denom[1,1] = o.R * AR
			if o.ML
				for k ∈ 1:ncols(o.v)
					numer_l = view(o.numerw,:,k)
					o.dist[k+first(o.WeightGrp[w])-1] = numer_l'invsym(o.denom[1,1])*numer_l
				end
				isone(w) && (o.statDenom = o.denom[1,1])  # original-sample denominator
			else
				invdenom = invsym(o.denom[1,1])
				for k ∈ 1:ncols(o.v)
					numer_l =view(o.numerw,:,k)
					o.dist[k+first(o.WeightGrp[w])-1] = o.numer_l'invdenom*numer_l
					o.ustar = o.B>0 ? view(o.v,:,k) .* o.ü : o.ü
					if o.scorebs
						if o.haswt  # Center variance if interpolated
							o.ustar .-= o.wt'o.ustar * o.ClustShare
						else
							o.ustar .-= colsum(o.ustar) * o.ClustShare  # Center variance if interpolated
						end
					else
						o.ustar .-= X₁₂B(o.X₁, o.X₂, view(o.βdev,:,k))  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
					end
					o.dist[k+first(o.WeightGrp[w])-1] ./= (tmp = symcross(o.ustar, o.wt))
				end
				isone(w) && (o.statDenom = o.denom[1,1] * tmp)  # original-sample denominator
			end
		end
	end
end


# like Mata panelsetup() but can group on multiple columns, like sort(), and faster. But doesn't take minobs, maxobs arguments.
function panelsetup(X::AbstractArray{S} where S, colinds::Vector{T} where T<:Integer)
  N = nrows(X)
  info = Vector{UnitRange{Int64}}(undef, N)
  lo = p = 1
  id = view(X,1,colinds)
  @inbounds for hi ∈ 2:N
  	if (tmp=view(X,hi,colinds)) ≠ id
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
  N = nrows(X)
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


# Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
# except use provided lo and hi values for 0 and 1
count_binary(N::Integer, lo::Number, hi::Number) = N<=1 ? [lo  hi] :
														  (tmp = count_binary(N-1, lo, hi);
														   [fill(lo, 1, ncols(tmp)) fill(hi, 1, ncols(tmp)) ;
																			  tmp                     tmp     ])

# cross-tab sum of a column vector w.r.t. given panel info and fixed-effect var
# one row per FE, one col per other grouping
function crosstabFE(o::StrBootTest{T}, v::AbstractVector, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, o.NFE, nrows(info))
  if length(info)>0
	  for i ∈ 1:nrows(info)
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
	  for i ∈ 1:nrows(info)
	    FEIDi = @view o._FEID[info[i]]
	    vi    = @view       v[info[i]]
	    @inbounds @simd for j ∈ eachindex(vi, FEIDi)
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
	  for i ∈ 1:nrows(info)
	    FEIDi = @view o._FEID[info[i]]
	    vi    = @view       v[info[i],:]
	    @inbounds @simd for j ∈ 1:length(FEIDi)
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


# given a pre-configured boottest linear model with one-degree null imposed, compute distance from target p value of boostrapped one associated with given value of r
# used with optimize() to construct confidence intervals
# performs no error checking
function r_to_p(o::StrBootTest{T}, r::AbstractVector{T}) where T
  o.r = r
  o.dirty = true
  getpadj(o)
end


# Chandrupatla 1997, "A new hybrid quadratic-bisection algorithm for finding the zero of a nonlinear function without using derivatives"
# x₁, x₂ must bracket the true value, with f₁=f(x₁) and f₂=f(x₂)
function search(o::StrBootTest{T}, α::T, f₁::T, x₁::T, f₂::T, x₂::T) where T<:Real
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

		((o.B>0 && abs(fx - α) < (1+(o.ptype==equaltail)) / o.BFeas * 1.000001) || ≈(x₂, x₁, rtol=o.rtol)) &&
			return abs(f₁ - α) < abs(f₂ - α) ? x₁ : x₂

		ϕ₁ = (f₁ - f₂) / (f₃ - f₂)
		ϕ₁² = ϕ₁^2
		xi1 = (x₁ - x₂) / (x₃ - x₂)
		t = ϕ₁² > xi1 || xi1 > 2 * ϕ₁ - ϕ₁² ?
					half :
					clamp(((f₃ - α) / (f₁ - f₂) + (x₃ - x₁) / ((x₂ - x₁) * (f₃ - f₁)) * (f₂ - α)) * (f₁ - α) / (f₃ - f₂), T(0.000001), T(0.999999))
  end
end

# derive wild bootstrap-based CI, for case of linear model with one-degree null imposed
# and generate plot data
function plot(o::StrBootTest{T}) where T
  _r = copy(o.r)
  α = one(T) - o.level
  o.gridpoints[ismissing.(o.gridpoints)] .= 25

  boottest!(o)
  if !o.ARubin
		halfwidth = T.(-1.5 * quantile(Normal(), α/2)) .* sqrtNaN.(diag(getV(o)))
		o.confpeak = getb(o) + o.r
  else
		halfwidth = abs.(o.confpeak) * T.(quantile(Normal(), getpadj(o, classical=true)/2) / quantile(Normal(), α/2))
  end

	if isone(o.q)  # 1D plot
		α<=0 && (α = T(.05))  # if level=100, no CI constructed, but we need a reasonable α to choose graphing bounds

		if α > 0 && ncols(o.v)-1 <= 1/α-1e6
			throw(ErrorException("need at least $(ceil(1/α)) replications to resolve a $(o.level)% two-sided confidence interval."))
		end

		p_lo, p_hi = T(NaN), T(NaN)
		if ismissing(o.gridmin[1]) || ismissing(o.gridmax[1])
			if o.B>0  # initial guess based on classical distribution
				lo = Vector{T}(ismissing(o.gridmin[1]) ? o.confpeak - halfwidth : o.gridmin)  # signal compiler that lo and hi cannot be missing now
				hi = Vector{T}(ismissing(o.gridmax[1]) ? o.confpeak + halfwidth : o.gridmax)
			else
				tmp = vec(sqrtNaN.(o.statDenom)) * cquantile(o.small ? TDist(o.dof_r) : Normal(), α/2)
				lo = Vector{T}(ismissing(o.gridmin[1]) ? o.confpeak - tmp : o.gridmin)
				hi = Vector{T}(ismissing(o.gridmax[1]) ? o.confpeak + tmp : o.gridmax)
				if o.scorebs && !o.null && !o.willplot  # if doing simple Wald test with no graph, we're done
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

		o.plotX = reshape(range(lo[1], hi[1], length=o.gridpoints[1]),:,1)
		o.plotY = fill(T(NaN), length(o.plotX))
		o.plotY[1  ] = p_lo
		o.plotY[end] = p_hi
		p_confpeak = o.WREnonARubin ? T(NaN) : o.twotailed ? one(T) : T(.5)

		c = clamp((floor(Int, (o.confpeak[1] - lo[1]) / (hi[1] - lo[1]) * (o.gridpoints[1] - 1)) + 2), 1, o.gridpoints[1]+1)  # insert original point estimate into grid
		o.plotX = [@view o.plotX[1:c-1,:] ; o.confpeak ; @view o.plotX[c:end,:]]
		insert!(o.plotY, c, p_confpeak)
	else  # 2D plot
		lo = Vector{T}(undef, 2)
		hi = Vector{T}(undef, 2)
		for d ∈ 1:o.dof
			lo[d] = ismissing(o.gridmin[d]) ? o.confpeak[d] - halfwidth[d] : o.gridmin[d]
			hi[d] = ismissing(o.gridmax[d]) ? o.confpeak[d] + halfwidth[d] : o.gridmax[d]
		end
		o.plotX = [repeat(range(lo[1], hi[1], length=o.gridpoints[1]), inner=o.gridpoints[2]) repeat(range(lo[2], hi[2], length=o.gridpoints[2]), outer=o.gridpoints[1])]
		o.plotY = fill(T(NaN), nrows(o.plotX))
	end

	isnan(o.plotY[1]) && (o.plotY[1] = r_to_p(o, view(o.plotX,1,:)))  # do in this order for widest interpolation
	@views for i ∈ length(o.plotY):-1:2
		isnan(o.plotY[i]) && (o.plotY[i] = r_to_p(o, view(o.plotX,i,:)))
	end

	if any(isnan.(o.plotY))
		o.CI = [T(-Inf) T(Inf)]
	elseif isone(o.q) && o.level<100 # find CI bounds
		_CI = Vector{T}(undef, nrows(o.plotY))
		for i in eachindex(_CI)  # map() version hampers type inference in Julia 1.6.2
			_CI[i] = isnan(o.plotY[i]) ? o.plotY[i] : T(o.plotY[i] > α)
		end
		_CI = _CI[2:end] - _CI[1:end-1]
		lo = T.(findall(x->x== 1, _CI))
		hi = T.(findall(x->x==-1, _CI))
		if iszero(length(lo)) && iszero(length(hi))
			o.CI = [T(-Inf) T(Inf)]
		else
			if iszero(length(lo))
				lo = [T(-Inf)]
			elseif iszero(length(hi))
				hi = [T(Inf)]
			else
				lo[1  ] > hi[1  ] && (lo = [T(-Inf) ; lo    ]) # non-rejection ranges that are not within grid range
				lo[end] > hi[end] && (hi = [hi      ; T(Inf)])
			end
			o.CI = [lo hi]

			for i ∈ 1:length(lo), j ∈ 1:2
				if !isinf(o.CI[i,j])
					t = Int(o.CI[i,j])
					o.CI[i,j] = search(o, α, o.plotY[t], o.plotX[t], o.plotY[t+1], o.plotX[t+1])
				end
			end
		end
	end

  if @isdefined c  # now that it's done helping graph look good, remove peak point from returned grid for evenness, for Bayesian sampling purposes
		o.peak = (X = view(o.plotX,c,:), p = o.plotY[c])
		o.plotX = view(o.plotX,[1:c-1; c+1:nrows(o.plotX)],:)
		deleteat!(o.plotY, c)
  end

	o.r = _r; o.dirty = true  # restore backups
	o.notplotted = false
end

"Structure to store wild bootstrap test results"
struct BoottestResult{T}
  stat::T; stattype::String
  p::T; padj::T
  reps::Int64; repsfeas::Int64
  NBootClust::Int64
  dof::Int64; dof_r::Int64
  plot::Union{Nothing, NamedTuple{(:X, :p), Tuple{Matrix{T},Vector{T}}}}
  peak::Union{Nothing, NamedTuple{(:X, :p), Tuple{Vector{T}, T}}}
  CI::Union{Nothing, Matrix{T}}
  dist::Matrix{T}
  b::Vector{T}
  V::Matrix{T}
  auxweights::Union{Nothing,Matrix{T}}
  M::StrBootTest
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

"Return actual  number of replications in wild bootstrap test, subject to enumeration of Rademacher draws"
repsfeas(o::BoottestResult) = o.repsfeas

"Return number of bootstrapping clusters in wild bootstrap test"
NBootClust(o::BoottestResult) = o.NBootClust

"Return degrees of freedom wild bootstrap test"
dof(o::BoottestResult) = o.dof

"Return residual degrees of freedom wild bootstrap test"
dof_r(o::BoottestResult) = o.dof_r

"Return data for confidence plot of wild bootstrap test"
plotpoints(o::BoottestResult) = o.plot

"Return parameter value with peak p value in wild bootstrap test"
peak(o::BoottestResult) = o.peak

"Return confidence interval matrix from wild bootstrap test, one row per disjoint piece"
CI(o::BoottestResult) = o.CI

"Return bootstrap distribution of statistic or statistic numerator in wild bootstrap test"
dist(o::BoottestResult) = o.dist

"Return auxilliary weight matrix for wild bootstrap"
auxweights(o::BoottestResult) = o.auxweights


"""
wildboottest([T=Float32], H₀::Tuple{AbstractMatrix, AbstractVector}; resp, <optional keyword arguments>) <: WildBootTest.BoottestResult

Function to perform wild-bootstrap-based hypothesis test

# Positional arguments
* `T::DataType`: data type for inputs, results, and computations: Float32 (default) or Float64
* `H₀::Tuple{AbstractMatrix, AbstractVector}`: required matrix-vector tuple (R,r) expressesing the null Rβ=r

# Required keyword argument
* `resp::AbstractVector`: response/dependent variable (y/y₁)

# Optional keyword arguments
* `predexog::AbstractVecOrMat`: exogenous predictors (X/X₁)
* `predendog::AbstractVecOrMat`: endogenous predictors (Y₂)
* `inst::AbstractVecOrMat`: instruments (X₂)
* `H₁::Tuple{AbstractMatrix, AbstractVector`: model constraints; same format as for H₀
* `clustid::AbstractVecOrMat{<:Integer}`: data vector/matrix of error and bootstrapping cluster identifiers; see Notes 
* `nbootclustvar::Integer=1`: number of bootstrap-clustering variables
* `nerrclustvar::Integer=nbootclustvar`: number of error-clustering variables
* `hetrobust::Bool=false`: true unless errors are treated as iid
* `feid::AbstractVector{<:Integer}`: data vector for fixed effect group identifier
* `fedfadj::Bool=true`: true if small-sample adjustment should reflect number of fixed effects
* `obswt`: observation weight vector
* `fweights::Bool=false`: true for frequency weights
* `maxmatsize::Number`: maximum size of auxilliary weight matrix (v), in gigabytes
* `ptype::PType=symmetric`: p value type (symmetric, equaltail, lower, upper)
* `bootstrapc::Bool=false`: true for bootstrap-c
* `LIML::Bool=false`: true for LIML or Fuller LIML
* `Fuller::Number`: Fuller factor
* `κ::Number`: fixed κ for _k_-class estimation
* `ARubin::Bool=false`: true for Anderson-Rubin test
* `small::Bool=true`: true for small-sample corrections
* `scorebs::Bool=false`: true for score bootstrap instead of wild bootstrap
* `reps::Integer=999`: number of bootstrap replications
* `imposenull::Bool=true`: true to impose null
* `auxwttype::AuxWtType=rademacher`: auxilliary weight type (rademacher mammen webb normal gamma)
* `rng::AbstractRNG=MersenneTwister()`: randon number generator
* `level::Number=.95`: significance level (0-1)
* `rtol::Number=1e-6`: tolerance for CI bound determination
* `madjtype::MAdjType=none`: multiple hypothesis adjustment (none, bonferroni, sidak)
* `NH0::Integer=1`: number of hypotheses tested, including one being tested now
* `ML::Bool=false`: true for (nonlinear) ML estimation
* `scores::AbstractVecOrMat`: for ML, pre-computed scores
* `β::AbstractVector`: for ML, parameter estimates
* `A::AbstractMatrix`: for ML, covariance estimates
* `gridmin`: vector of graph lower bounds, max length 2, entries may be `missing`
* `gridmax`: vector of graph upper bounds
* `gridpoints`: vector of number of sampling points
* `diststat::String=""`: "numer" or "t" to save associated bootstrap distribution
* `getCI::Bool=true`: whether to return CI
* `getplot::Bool=getCI`: whether to generate plot data
* `getauxweights::Bool=false`: whether to save auxilliary weight matrix (v)

# Notes
Order the columns of `clustid` this way:
1. Variables only used to define bootstrapping clusters, as in the subcluster bootstrap.
2. Variables used to define both bootstrapping and error clusters.
3. Variables only used to define error clusters.
In the most common case, `clustid` is a single column of type 2.

The code does not handle missing data values: all data and identifier matrices must 
match the estimation sample.

"""
function wildboottest(T::DataType,
					  H₀::Tuple{AbstractMatrix, AbstractVector};
					  resp::AbstractVector,
					  predexog::AbstractVecOrMat=zeros(T,0,0),
					  predendog::AbstractVecOrMat=zeros(T,0,0),
					  inst::AbstractVecOrMat=zeros(T,0,0),
					  H₁::Tuple{AbstractMatrix, AbstractVector}=(zeros(T,0,0), zeros(T,0)),
					  clustid::AbstractVecOrMat=zeros(Int,0,0),  # bootstrap-only clust vars, then boot&err clust vars, then err-only clust vars
					  nbootclustvar::Integer=1,
					  nerrclustvar::Integer=nbootclustvar,
					  hetrobust::Bool=false,
					  feid::AbstractVector=Vector(undef,0),
					  fedfadj::Bool=true,
					  obswt::Union{AbstractVector,UniformScaling{Bool}}=I,
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
					  madjtype::MAdjType=none,
					  NH0::Integer=1,
					  ML::Bool=false,
					  scores::AbstractVecOrMat=Matrix{Float32}(undef,0,0),
					  β::AbstractVector=zeros(T,0),
					  A::AbstractMatrix=zeros(T,0,0),
					  gridmin::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = [missing],
					  gridmax::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Number = [missing],
					  gridpoints::Union{Vector{S},Vector{Union{S,Missing}}} where S<:Integer = [missing],
					  diststat::String="",
					  getCI::Bool=true,
					  getplot::Bool=getCI,
					  getauxweights::Bool=false)

  if nrows(H₀[1])==2
		_gridmin    = ismissing.(gridmin)==[true] ? [missing; missing] : gridmin
		_gridmax    = ismissing.(gridmax)==[true] ? [missing; missing] : gridmax
		_gridpoints = ismissing.(gridpoints)==[true] ? [missing; missing] : gridpoints
  else
		_gridmin, _gridmax, _gridpoints = gridmin, gridmax, gridpoints
  end
  M = StrBootTest{T}(H₀[1], H₀[2], H₁[1], H₁[2], resp, predexog, predendog, inst, obswt, fweights, LIML, Fuller, κ, ARubin,
	reps, auxwttype, rng, maxmatsize, ptype, imposenull, scorebs, !bootstrapc, clustid, nbootclustvar, nerrclustvar, hetrobust, small,
	feid, fedfadj, level, rtol, madjtype, NH0, ML, β, A, scores, getplot,
	map(x->ismissing(x) ? missing : T(x), _gridmin),
	map(x->ismissing(x) ? missing : T(x), _gridmax),
	map(x->ismissing(x) ? missing : Int32(x), _gridpoints))

  if getplot || (level<1 && getCI)
		plot = getplot ? _getplot(M) : nothing
		peak = getpeak(M)
		CI = level<1 & getCI ? _getCI(M) : nothing
  else
		CI = plot = peak = nothing
  end

  BoottestResult{T}(getstat(M),
					isone(nrows(H₀[1])) ? (small ? "t" : "z") : (small ? "F" : "χ²"),
					getp(M), getpadj(M), getreps(M), getrepsfeas(M), getNBootClust(M), getdf(M), getdf_r(M), plot, peak, CI,
					getdist(M,diststat),
					getb(M), getV(M),
					getauxweights && reps>0 ? getauxweights(M) : nothing, M)
end

wildboottest(H₀::Tuple{AbstractMatrix, AbstractVector}; args...) = wildboottest(Float32, H₀; args...)

end # module