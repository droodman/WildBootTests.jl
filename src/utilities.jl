# core functions not referencing StrBootest or StrEstimator "classes"

@inline sqrtNaN(x) = x<0 ? typeof(x)(NaN) : sqrt(x)
@inline invsym(X) = iszero(nrows(X)) ? Symmetric(X) : inv(Symmetric(X))

@inline symcross(X::AbstractVecOrMat, wt::AbstractVector) = Symmetric(cross(X,wt,X))  # maybe bad name since it means cross product in Julia
function cross(X::AbstractVecOrMat{T}, wt::AbstractVector{T}, Y::AbstractVecOrMat{T}) where T
	dest = Matrix{T}(undef, ncols(X), ncols(Y))
	if iszero(nrows(wt))
		mul!(dest, X', Y)
	elseif ncols(X)>ncols(Y)
		mul!(dest, X', wt.*Y)
	else
		mul!(dest, (X.*wt)', Y)
	end
end
function crossvec(X::AbstractMatrix{T}, wt::AbstractVector{T}, Y::AbstractVector{T}) where T
  dest = Vector{T}(undef, ncols(X))
  if iszero(nrows(wt))
		mul!(dest, X', Y)
	elseif ncols(X)>ncols(Y)
		mul!(dest, X', wt.*Y)
	else
		mul!(dest, (X.*wt)', Y)
	end
end

@inline vHadw(v::AbstractArray, w::AbstractVector) = iszero(nrows(w)) ? v : v .* w  # not type-stable

@inline nrows(X::AbstractArray) = size(X,1)
@inline ncols(X::AbstractArray) = size(X,2)

@inline colsum(X::AbstractArray) = iszero(length(X)) ? similar(X, 1, size(X)[2:end]...) : sum(X, dims=1)
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))

function X₁₂B(X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
	dest = X₁ * view(B,1:size(X₁,2),:)
	length(dest)>0 && length(X₂)>0 && matmulplus!(dest, X₂, B[size(X₁,2)+1:end,:])
	dest
end
function X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractVector)
	dest = X₁ * view(B,1:size(X₁,2))
	length(dest)>0 && length(X₂)>0 && matmulplus!(dest, X₂, B[size(X₁,2)+1:end])
	dest
end

function coldot!(dest::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T  # colsum(A .* B)
  fill!(dest, zero(T))
	@tturbo for i ∈ eachindex(axes(A,2),axes(B,2)), j ∈ eachindex(axes(A,1),axes(B,1))
		dest[i] += A[j,i] * B[j,i]
  end
end
function coldot(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
  dest = Matrix{T}(undef, 1, size(A,2))
  coldot!(dest, A, B)
  dest
end
coldot(A::AbstractMatrix) = coldot(A, A)
coldot(A::AbstractVector, B::AbstractVector) = [dot(A,B)]

# colsum(A .* B); dest should be a one-row matrix
function coldotplus!(dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1)
	  dest[i] += A[j,i] * B[j,i]
  end
end

# compute the norm of each col of A using quadratic form Q; dest should be a one-row matrix
function colquadform!(dest::AbstractMatrix{T}, Q::AbstractMatrix{T}, A::AbstractMatrix{T}) where T
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
	  dest[i] += A[j,i] * Q[k,j] * A[k,i]
  end
end


 # From given row of given matrix, substract inner products of corresponding cols of A & B with quadratic form Q
function colquadformminus!(X::AbstractMatrix, row::Integer, Q::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
    X[row,i] -= A[j,i] * Q[k,j] * B[k,i]
  end
  X
end
colquadformminus!(X::AbstractMatrix, Q::AbstractMatrix, A::AbstractMatrix) = colquadformminus!(X, 1, Q, A, A)

function matmulplus!(A::Matrix, B::Matrix, C::Matrix)  # add B*C to A in place
	@tturbo for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2)), j ∈ eachindex(axes(B,2),axes(C,1))
		A[i,k] += B[i,j] * C[j,k]
	end
end
function matmulplus!(A::Vector, B::Matrix, C::Vector)  # add B*C to A in place
	@tturbo for j ∈ eachindex(axes(B,2),C), i ∈ eachindex(axes(A,1),axes(B,1))
		A[i] += B[i,j] * C[j]
	end
end
# like Mata panelsetup() but can group on multiple columns, like sort(). But doesn't take minobs, maxobs arguments.
function panelsetup(X::AbstractArray{S} where S, colinds::AbstractVector{T} where T<:Integer)
  N = nrows(X)
  info = Vector{UnitRange{Int64}}(undef, N)
  lo = p = 1
  @inbounds for hi ∈ 2:N
    for j ∈ colinds
      if X[hi,j] ≠ X[lo,j]
        info[p] = lo:hi-1
        lo = hi
        p += 1
        break
      end
  	end
  end
  info[p] = lo:N
  resize!(info, p)
  info
end
# Like above but also return standardized ID variable, starting from 1
function panelsetupID(X::AbstractArray{S} where S, colinds::UnitRange{T} where T<:Integer)
  N = nrows(X)
  info = Vector{UnitRange{Int64}}(undef, N)
  ID = ones(Int64, N)
  lo = p = 1
  @inbounds for hi ∈ 2:N
    for j ∈ colinds
      if X[hi,j] ≠ X[lo,j]
        info[p] = lo:hi-1
        lo = hi
        p += 1
        break
      end
  	end
	  ID[hi] = p
  end
  info[p] = lo:N
  resize!(info, p)
  info, ID
end


# Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
# except use provided lo and hi values for 0 and 1
count_binary(N::Integer, lo::Number, hi::Number) = N≤1 ? [lo  hi] :
														  (tmp = count_binary(N-1, lo, hi);
														   [fill(lo, 1, ncols(tmp)) fill(hi, 1, ncols(tmp)) ;
																			            tmp                     tmp     ])

# unweighted panelsum!() along first axis of a VecOrMat
function panelsum!(dest::AbstractVecOrMat, X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	iszero(length(X)) && return
	J = CartesianIndices(axes(X)[2:end])
	eachindexJ = eachindex(J)
	@inbounds for g in eachindex(info)
		f, l = first(info[g]), last(info[g])
		fl = f+1:l
		if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj]
				for i ∈ fl
					tmp += X[i,Jj]
				end
				dest[g,Jj] = tmp
			end
		else
			@simd for j ∈ eachindexJ
				dest[g,J[j]] = X[f,J[j]]
			end
		end
	end
end
# single-weighted panelsum!() along first axis of a VecOrMat
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
    fl = f+1:l
		_wt = wt[f]
    if f<l
      for j ∈ eachindexJ
        Jj = J[j]
        tmp = X[f,Jj] * _wt
        for i ∈ fl
          tmp += X[i,Jj] * wt[i]
        end
        dest[g,Jj] = tmp
      end
    else
      for j ∈ eachindexJ
        dest[g,J[j]] = X[f,J[j]] * _wt
      end
    end
  end
end
function panelsum(X::AbstractVector{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(nrows(wt))
		panelsum(X, info)
	else
		dest = Vector{T}(undef, length(info))
		panelsum!(dest, X, wt, info)
		dest
	end
end
function panelsum(X::AbstractMatrix{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(nrows(wt))
		panelsum(X, info)
	else
		dest = Matrix{T}(undef, length(info), size(X,2))
		panelsum!(dest, X, wt, info)
		dest
	end
end
function panelsum(X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	dest = similar(X, length(info), size(X)[2:end]...)
	panelsum!(dest, X, info)
	dest
end
function panelsum2(X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(length(X₁))
		panelsum(X₂,wt,info)
	elseif iszero(length(X₂))
		panelsum(X₁,wt,info)
	else
		dest = Matrix{T}(undef, length(info), ncols(X₁)+ncols(X₂))
		panelsum!(view(dest, :,           1:ncols(X₁   )), X₁, wt, info)
		panelsum!(view(dest, :, ncols(X₁)+1:size(dest,2)), X₂, wt, info)
		dest
	end
end

# macros to efficiently handle result = input
macro panelsum(X, info)
	:(local _X = $(esc(X)); iszero(length($(esc(info)))) || length($(esc(info)))==nrows(_X) ? _X : panelsum(_X, $(esc(info)) ) )
end
macro panelsum(X, wt, info)
	:( panelsum($(esc(X)), $(esc(wt)), $(esc(info))) )
end
macro panelsum2(X₁, X₂, wt, info)
	:( panelsum2($(esc(X₁)), $(esc(X₂)), $(esc(wt)), $(esc(info))) )
end

import Base.size
struct FakeArray{N} <: AbstractArray{Bool,N} size::Tuple{Vararg{Int64,N}} end # AbstractArray with almost no storage, just for LinearIndices() conversion         
FakeArray(size...) = FakeArray{length(size)}(size)
size(X::FakeArray) = X.size