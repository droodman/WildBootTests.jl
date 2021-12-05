# core functions not referencing StrBootest or StrEstimator "classes"

@inline sqrtNaN(x) = x<0 ? typeof(x)(NaN) : sqrt(x)
@inline invsym(X) = iszero(length(X)) ? X : inv(Symmetric(X))
@inline eigensym(X) = eigen(Symmetric(X))  # does eigen recognize symmetric matrices?

@inline symcross(X::AbstractVecOrMat, wt::AbstractVector) = Symmetric(cross(X,wt,X))  # maybe bad name since it means cross product in Julia
@inline function cross(X::AbstractVecOrMat{T}, wt::AbstractVector{T}, Y::AbstractVecOrMat{T}) where T
  dest = Matrix{T}(undef, ncols(X), ncols(Y))
  if iszero(nrows(wt))
		mul!(dest, X', Y)
	elseif ncols(X)>ncols(Y)
		mul!(dest, X', wt.*Y)
	else
		mul!(dest, (X.*wt)', Y)
	end
end
@inline function crossvec(X::AbstractMatrix{T}, wt::AbstractVector{T}, Y::AbstractVector{T}) where T
  dest = Vector{T}(undef, ncols(X))
  if iszero(nrows(wt))
		mul!(dest, X', Y)
	elseif ncols(X)>ncols(Y)
		mul!(dest, X', wt.*Y)
	else
		mul!(dest, (X.*wt)', Y)
	end
end

@inline vHadw(v::AbstractArray, w::AbstractVector) = iszero(nrows(w)) ? v : v .* w

@inline nrows(X::AbstractArray) = size(X,1)
@inline ncols(X::AbstractArray) = size(X,2)

@inline colsum(X::AbstractArray) = iszero(length(X)) ? similar(X, 1, size(X)[2:end]...) : sum(X, dims=1)
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))

@inline wtsum(wt::AbstractVector, X::AbstractMatrix) = iszero(nrows(wt)) ? sum(X,dims=1) : reshape(X'wt,1,:)  # return 1xN Matrix

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

function minusX₁₂B!(dest::AbstractVecOrMat, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
  matmulminus!(dest, X₁, view(B,           1:size(X₁,2),:))
  matmulminus!(dest, X₂, view(B,size(X₁,2)+1:size(B ,1),:))
end

import Base.*  # extend * to left- and right-multiply arrays by vec or mat
@inline *(A::AbstractArray, B::AbstractVecOrMat) = reshape(reshape(A,:,size(B,1)) * B, size(A)[1:end-1]..., size(B)[2:end]...)
@inline *(A::AbstractVecOrMat, B::AbstractArray) = reshape(A * reshape(B,size(A,2),:), size(A)[1:end-1]..., size(B)[2:end]...)

function coldot!(dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)  # colsum(A .* B)
  s2 = axes(A,2)
  @tturbo for i ∈ s2
	  dest[i] = A[1,i] * B[1,i]
  end
  if size(A,1)>1
	  @tturbo for i ∈ s2, j ∈ 2:size(A,1)
	    dest[i] += A[j,i] * B[j,i]
	  end
  end
end
function coldot!(dest::AbstractMatrix, A::AbstractMatrix)  # colsum(A .* A)
  s2 = axes(A,2)
  @tturbo for i ∈ s2
	  dest[i] = A[1,i] ^ 2
  end
  if size(A,1)>1
	  @tturbo for i ∈ s2, j ∈ 2:size(A,1)
	    dest[i] += A[j,i] ^ 2
	  end
  end
end
function coldot!(dest::AbstractArray, A::AbstractArray, B::AbstractArray)
  J = CartesianIndices(axes(A)[2:end])
  eachindexJ = eachindex(J)
  @tturbo for j ∈ eachindexJ
	  dest[1,J[j]] = A[1,J[j]] * B[1,J[j]]
  end
  if size(A,1) > 1
	  @tturbo for j ∈ eachindexJ, i ∈ 2:size(A,1)
	    dest[1,J[j]] += A[i,J[j]] * B[i,J[j]]
	  end
  end
end
function coldot!(dest::AbstractArray, A::AbstractArray)
  J = CartesianIndices(axes(A)[2:end])
  eachindexJ = eachindex(J)
  @tturbo for j ∈ eachindexJ
	  dest[1,J[j]] = A[1,J[j]] ^ 2
  end
  if size(A,1) > 1
	  @tturbo for j ∈ eachindexJ, i ∈ 2:size(A,1)
	    dest[1,J[j]] += A[i,J[j]] ^ 2
	  end
  end
end
function coldot(args...)
  dest = Array{promote_type(eltype.(args)...)}(undef, 1, size(args[1])[2:end]...)
  coldot!(dest, args...)
  dest
 end
coldot(A::AbstractVector, B::AbstractVector) = [dot(A,B)]

function coldotplus!(dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)  # colsum(A .* B)
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1)
	  dest[i] += A[j,i] * B[j,i]
  end
end

# compute the norm of each col of A using quadratic form Q
function colquadform(Q::AbstractMatrix, A::AbstractMatrix) :: AbstractVector
  dest = zeros(promote_type(eltype(Q), eltype(A)), size(A,2))
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
	  dest[i] += A[j,i] * Q[k,j] * A[k,i]
  end
  dest
end

 # From given row of given matrix, substract inner products of corresponding cols of A & B with quadratic form Q
function colquadformminus!(X::AbstractMatrix, row::Integer, Q::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
    X[row,i] -= A[j,i] * Q[k,j] * B[k,i]
  end
  X
end
colquadformminus!(X::AbstractMatrix, Q::AbstractMatrix, A::AbstractMatrix) = colquadformminus!(X, 1, Q, A, A)

function matmulminus!(A::Matrix, B::Matrix, C::AbstractMatrix)  # add B*C to A in place
	@tturbo for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2)), j ∈ eachindex(axes(B,2),axes(C,1))
		A[i,k] -= B[i,j] * C[j,k]
	end
end

function matmulplus!(A::Matrix, B::Vector, C::Matrix)  # add B*C to A in place
	@tturbo for k ∈ eachindex(axes(A,2), axes(C,2)), i ∈ eachindex(axes(A,1),axes(B,1))
		A[i,k] += B[i] * C[k]
	end
end
function matmulplus!(A::Matrix, B::Matrix, C::Matrix)  # add B*C to A in place
	@tturbo for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2)), j ∈ eachindex(axes(B,2),axes(C,1))
		A[i,k] += B[i,j] * C[j,k]
	end
end
function matmulplus!(A::Vector, B::Matrix, C::Vector)  # add B*C to A in place
	@tturbo for j ∈ eachindex(axes(B,2),axes(C,1)), i ∈ eachindex(axes(A,1),axes(B,1))
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

# unweighted panelsum!() along first axis of an Array
function panelsum!(dest::AbstractArray, X::AbstractArray, info::Vector{UnitRange{T}} where T<:Integer)
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
				@tturbo for i ∈ fl
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
    fl = f+1:l
		_wt = wt[f]
		if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj] * _wt
				@tturbo for i ∈ fl
					tmp += X[i,Jj] * wt[i]
				end
				dest[g,Jj] = tmp
			end
		else
			@simd for j ∈ eachindexJ
				dest[g,J[j]] = X[f,J[j]] * _wt
			end
		end
	end
end
# multiple-weighted panelsum!() along first axis of an Array
# *2nd* dimension of resulting 3-D array corresponds to cols of wt;
# this facilitates reshape() to 2-D array in which results for each col of wt are stacked vertically
function panelsum!(dest::AbstractArray, X::AbstractArray, wt::AbstractMatrix, info::Vector{UnitRange{T}} where T<:Integer)
	iszero(length(X)) && return
	if iszero(length(info)) || nrows(info)==nrows(X)
		@inbounds @simd for i ∈ axes(wt,2)
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
			if f<l
				for j ∈ eachindexJ
					Jj = J[j]
					tmp = X[f,Jj] * _wt
					@tturbo for i ∈ fl
						tmp += X[i,Jj] * wt[i,k]
					end
					dest[g,k,Jj] = tmp
				end
			else
				@simd for j ∈ eachindexJ
					dest[g,k,J[j]] = X[f,J[j]] * _wt
				end
			end
		end
	end
end
function panelsum(X::AbstractArray, wt::AbstractVecOrMat, info::Vector{UnitRange{T}} where T<:Integer)
	if iszero(nrows(wt))
		panelsum(X, info)
	else
		dest = similar(X, length(info), size(wt)[2:end]..., size(X)[2:end]...)
		panelsum!(dest, X, wt, info)
		dest
	end
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
		panelsum!(view(dest, Vector{Colon}(undef,ndims(wt))...,           1:ncols(X₁      )), X₁, wt, info)
		panelsum!(view(dest, Vector{Colon}(undef,ndims(wt))..., ncols(X₁)+1:size(dest)[end]), X₂, wt, info)
		dest
	end
end

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

# SelectionMatrix type to efficiently represent selection matrices
import Base.size, Base.getindex, Base.*
struct SelectionMatrix{T} <: AbstractMatrix{T}
	p::Vector{Int64}
	size::Tuple{Int64,Int64}
end
SelectionMatrix(X::AbstractMatrix{T}) where T = SelectionMatrix{T}(X'collect(1:size(X,1)), size(X))
selectify(X) = X==I ? I :
                      isa(X,AbstractMatrix) &&
											  length(X) > 0 &&
											  all(sum(isone.(X), dims=1) .== 1) &&
												all(sum(iszero.(X), dims=1) .== size(X,1)-1) ? SelectionMatrix(X) :
											                                                 X
size(X::SelectionMatrix) = X.size
getindex(X::SelectionMatrix, i, j) = i==X.p[j]
*(X::AbstractVector{T} where T, Y::SelectionMatrix)::SubArray{T, 1, Vector{T}} = view(X,Y.p)
*(X::AbstractMatrix{T} where T, Y::SelectionMatrix)::SubArray{T, 2, Matrix{T}} = view(X,:,Y.p)
*(X::SelectionMatrix, Y::AbstractMatrix{T} where T)::SubArray{T, 2, Matrix{T}} = view(Y, X.p, :)
*(X::SelectionMatrix, Y::AbstractVector{T} where T)::SubArray{T, 1, Vector{T}} = view(Y, X.p)

struct FakeArray{N} <: AbstractArray{Bool,N}  # AbstractArray with almost no storage, just for LinearIndices() conversion         
	size::Tuple{Vararg{Int64,N}}
end
FakeArray(size...) = FakeArray{length(size)}(size)
size(X::FakeArray) = X.size