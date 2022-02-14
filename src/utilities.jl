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

@inline nrows(X::AbstractArray) = size(X,1)
@inline ncols(X::AbstractArray) = size(X,2)

@inline colsum(X::AbstractArray) = iszero(length(X)) ? similar(X, 1, size(X)[2:end]...) : sum(X, dims=1)
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))

function X₁₂B(o::StrBootTest, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
	dest = X₁ * view(B,1:size(X₁,2),:)
	length(dest)>0 && length(X₂)>0 && matmulplus!(Val(o.turbo), dest, X₂, B[size(X₁,2)+1:end,:])
	dest
end
function X₁₂B(o::StrBootTest, X₁::AbstractArray, X₂::AbstractArray, B::AbstractVector)
	dest = X₁ * view(B,1:size(X₁,2))
	length(dest)>0 && length(X₂)>0 && matmulplus!(Val(o.turbo), dest, X₂, B[size(X₁,2)+1:end])
	dest
end

function coldot!(o::StrBootTest, dest::AbstractMatrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T  # colsum(A .* B)
  fill!(dest, zero(T))
	coldotplus!(Val(o.turbo), dest, A, B)
	nothing
end
function coldot(o::StrBootTest, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
  dest = Matrix{T}(undef, 1, size(A,2))
  coldot!(o, dest, A, B)
  dest
end
coldot(o::StrBootTest, A::AbstractMatrix) = coldot(o, A, A)
coldot(o::StrBootTest, A::AbstractVector, B::AbstractVector) = [dot(A,B)]

# colsum(A .* B); dest should be a one-row matrix
function coldotplus!(turbo::Val{true}, dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
  @tturbo for i ∈ eachindex(axes(A,2),axes(B,2)), j ∈ eachindex(axes(A,1),axes(B,1))
	  dest[i] += A[j,i] * B[j,i]
  end
	nothing
end
function coldotplus!(turbo::Val{false}, dest::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
	@inbounds Threads.@threads for i ∈ eachindex(axes(A,2),axes(B,2))
		@inbounds for j ∈ eachindex(axes(A,1),axes(B,1))
			dest[i] += A[j,i] * B[j,i]
		end
	end
	nothing
end

 # From given row of given matrix, substract inner products of corresponding cols of A & B with quadratic form Q; despite "!", puts result in return value too
function colquadformminus!(turbo::Val{true}, dest::AbstractMatrix, row::Integer, A::AbstractMatrix, Q::AbstractMatrix, B::AbstractMatrix)
  @tturbo for i ∈ axes(A,2), j ∈ axes(A,1), k ∈ axes(A,1)
    dest[row,i] -= A[j,i] * Q[k,j] * B[k,i]
  end
  dest
end
function colquadformminus!(turbo::Val{false}, dest::AbstractMatrix{T}, row::Integer, A::AbstractMatrix, Q::AbstractMatrix, B::AbstractMatrix) where T
  nt = Threads.nthreads()
  cs = [round(Int, size(A,2)/nt*i) for i ∈ 0:nt]
  if A===B 
		@inbounds Threads.@threads for t ∈ 1:nt
			tmp = Vector{T}(undef, size(Q,1))  # thread-safe scratchpad to minimize allocations
			@inbounds for i ∈ cs[t]+1:cs[t+1]
				v = view(A,:,i)
				mul!(tmp, Q, v)
				dest[row,i] -= v'tmp
			end
		end
	else
		@inbounds Threads.@threads for t ∈ 1:nt
			tmp = Vector{T}(undef, size(Q,1))  # thread-safe scratchpad to minimize allocations
			@inbounds for i ∈ cs[t]+1:cs[t+1]
				mul!(tmp, Q, view(B,:,i))
				dest[row,i] -= view(A,:,i)'tmp
			end
		end
	end
	dest
end
# function colquadformminus!(turbo::Val{false}, dest::AbstractMatrix{T}, A::AbstractMatrix, Q::AbstractArray{T,3}, B::AbstractMatrix) where T
#   nt = Threads.nthreads()
#   cs = [round(Int, size(A,2)/nt*i) for i ∈ 0:nt]
#   if A===B 
# 		@inbounds Threads.@threads for t ∈ 1:nt
# 			tmp = Matrix{T}(undef, size(Q,1), size(Q,2))  # thread-safe scratchpad to minimize allocations
# 			@inbounds for i ∈ cs[t]+1:cs[t+1]
# 				v = view(A,:,i)
# 				mul!(tmp, Q, v)
# 				mul!(dest[:,i], tmp', v)
# 			end
# 		end
# 	else
# 		@inbounds Threads.@threads for t ∈ 1:nt
# 			tmp = Matrix{T}(undef, size(Q,1), size(Q,2))  # thread-safe scratchpad to minimize allocations
# 			@inbounds for i ∈ cs[t]+1:cs[t+1]
# 				mul!(tmp, Q, view(B,:,i))
# 				mul!(dest[:,i], tmp', view(A,:,i))
# 			end
# 		end
# 	end
# 	dest
# end
colquadformminus!(o::StrBootTest, dest::AbstractMatrix, Q::AbstractMatrix, A::AbstractMatrix) = colquadformminus!(Val(o.turbo),dest,1,A,Q,A)
# compute negative of the norm of each col of A using quadratic form Q; dest should be a one-row matrix
function negcolquadform!(o::StrBootTest, dest::AbstractMatrix{T}, Q::AbstractMatrix{T}, A::AbstractMatrix{T}) where T
  fill!(dest, zero(T))
	colquadformminus!(Val(o.turbo),dest,1,A,Q,A)
	nothing
end

function matmulplus!(turbo::Val{true}, A::Matrix, B::Matrix, C::Matrix)  # add B*C to A in place
	@tturbo for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2)), j ∈ eachindex(axes(B,2),axes(C,1))
		A[i,k] += B[i,j] * C[j,k]
	end
	nothing
end
function matmulplus!(turbo::Val{true}, A::Vector, B::Matrix, C::Vector)  # add B*C to A in place
	@tturbo for j ∈ eachindex(axes(B,2),C), i ∈ eachindex(axes(A,1),axes(B,1))
		A[i] += B[i,j] * C[j]
	end
	nothing
end
function matmulplus!(turbo::Val{false}, A::VecOrMat{T}, B::Matrix{T}, C::VecOrMat{T}) where T  # add B*C to A in place
	BLAS.gemm!('N','N',one(T),B,C,one(T),A)
	nothing
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
function panelsum!(turbo::Val{true}, dest::AbstractVecOrMat, X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
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
function panelsum!(turbo::Val{false}, dest::AbstractVecOrMat, X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	iszero(length(X)) && return
	J = CartesianIndices(axes(X)[2:end])
	eachindexJ = eachindex(J)
	@inbounds Threads.@threads for g in eachindex(info)
		f, l = first(info[g]), last(info[g])
		fl = f+1:l
		@inbounds if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj]
				for i ∈ fl
					tmp += X[i,Jj]
				end
				dest[g,Jj] = tmp
			end
		else
			for j ∈ eachindexJ
				dest[g,J[j]] = X[f,J[j]]
			end
		end
	end
end
# single-weighted panelsum!() along first axis of a VecOrMat
function panelsum!(turbo::Val{true}, dest::AbstractVecOrMat{T}, X::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
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
			for j ∈ eachindexJ
				dest[g,J[j]] = X[f,J[j]] * _wt
			end
		end
	end
end
function panelsum!(turbo::Val{true}, dest::AbstractArray, X::AbstractArray{T,3} where T, info::Vector{UnitRange{S}} where S<:Integer)
  iszero(length(X)) && return
  @inbounds for g in eachindex(info)
    f, l = first(info[g]), last(info[g])
    fl = f+1:l
    if f<l
      for i ∈ axes(X,1), k ∈ axes(X,3)
        tmp = X[i,f,k]
        @tturbo for j ∈ fl
          tmp += X[i,j,k]
        end
        dest[i,g,k] = tmp
      end
    else
      for i ∈ axes(X,1), k ∈ axes(X,3)
        dest[i,g,k] = X[i,f,k]
      end
    end
  end
end
function panelsum!(turbo::Val{false}, dest::AbstractArray, X::AbstractArray, wt::AbstractVector, info::Vector{UnitRange{T}} where T<:Integer)
  iszero(length(X)) && return
  if iszero(length(info)) || nrows(info)==nrows(X)
    dest .= X .* wt
    return
  end
  J = CartesianIndices(axes(X)[2:end])
  eachindexJ = eachindex(J)
  @inbounds Threads.@threads for g in eachindex(info)
    f, l = first(info[g]), last(info[g])
    fl = f+1:l
		_wt = wt[f]
    @inbounds if f<l
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
function panelsum!(turbo::Val{false}, dest::AbstractArray, X::AbstractArray{T,3} where T, info::Vector{UnitRange{S}} where S<:Integer)
  iszero(length(X)) && return
  @inbounds Threads.@threads for g in eachindex(info)
    f, l = first(info[g]), last(info[g])
    fl = f+1:l
    @inbounds if f<l
      for i ∈ axes(X,1), k ∈ axes(X,3)
        tmp = X[i,f,k]
        for j ∈ fl
          tmp += X[i,j,k]
        end
        dest[i,g,k] = tmp
      end
    else
      for i ∈ axes(X,1), k ∈ axes(X,3)
        dest[i,g,k] = X[i,f,k]
      end
    end
  end
end
# groupwise inner product of two two data matrices
# 1st dimension of result corresponds to columns of X, second to rows of both, third to columns of Y
function panelcross!(turbo::Val{false}, dest::AbstractArray{T,3}, X::AbstractMatrix{T}, Y::AbstractMatrix{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	iszero(length(X)) && return
	if iszero(length(info)) || nrows(info)==nrows(X)
		@inbounds for i ∈ axes(Y,2)
			dest[:,:,i] .= X .* view(Y,:,i)
		end
		return
	elseif X===Y
    @inbounds Threads.@threads for g in eachindex(info)
      v = view(X,info[g],:)
      dest[:,g,:] = v'v
    end
  else
    @inbounds Threads.@threads for g in eachindex(info)
      infog = info[g]
      dest[:,g,:] = view(X,infog,:)'view(Y,infog,:)
    end
  end
end
function panelcross!(turbo::Val{false}, dest::AbstractMatrix{T}, X::AbstractVecOrMat{T}, Y::AbstractVector{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	iszero(length(X)) && return
	if iszero(length(info)) || nrows(info)==nrows(X)
		dest .= X .* Y
		return
	end
	if X===Y
    @inbounds Threads.@threads for g in eachindex(info)
      v = view(X,info[g])
      dest[1,g] = dot(v,v)
    end
  else
    @inbounds Threads.@threads for g in eachindex(info)
      infog = info[g]
			dest[:,g] .= view(X,infog,:)'view(Y,infog)
    end
  end
end

function panelsum(o::StrBootTest, X::AbstractVector{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(nrows(wt))
		panelsum(o, X, info)
	else
		dest = Vector{T}(undef, length(info))
		panelsum!(Val(o.turbo), dest, X, wt, info)
		dest
	end
end
function panelsum(o::StrBootTest, X::AbstractMatrix{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(nrows(wt))
		panelsum(o, X, info)
	else
		dest = Matrix{T}(undef, length(info), size(X,2))
		panelsum!(Val(o.turbo), dest, X, wt, info)
		dest
	end
end
function panelcross11(o::StrBootTest, X::AbstractVecOrMat{T}, Y::AbstractMatrix{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X,2), length(info), size(Y,2))
	panelcross!(Val(false), dest, X, Y, info)
	dest
end
function panelcross11(o::StrBootTest, X::AbstractVecOrMat{T}, Y::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = Matrix{T}(undef, size(X,2), length(info))
	panelcross!(Val(false), dest, X, Y, info)
	dest
end
function panelsum(o::StrBootTest, X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	dest = similar(X, length(info), size(X)[2:end]...)
	panelsum!(Val(o.turbo), dest, X, info)
	dest
end
function panelsum(o::StrBootTest{T}, X::AbstractArray{T,3}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X,1), size(info,1), size(X,3))
	panelsum!(Val(o.turbo), dest, X, info)
	dest
end
function panelsum2(o::StrBootTest, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(length(X₁))
		panelsum(o, X₂,wt,info)
	elseif iszero(length(X₂))
		panelsum(o, X₁,wt,info)
	else
		dest = Matrix{T}(undef, length(info), ncols(X₁)+ncols(X₂))
		panelsum!(Val(o.turbo), view(dest, :,           1:ncols(X₁   )), X₁, wt, info)
		panelsum!(Val(o.turbo), view(dest, :, ncols(X₁)+1:size(dest,2)), X₂, wt, info)
		dest
	end
end
function panelcross21(o::StrBootTest, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y::AbstractMatrix{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(ncols(X₁))
		panelcross11(o,X₂,Y,info)
	elseif iszero(ncols(X₂))
		panelcross11(o,X₁,Y,info)
	else
		dest = Array{T,3}(undef, ncols(X₁)+ncols(X₂), length(info), ncols(Y))
		panelcross!(Val(false), view(dest,           1:ncols(X₁   ), :, :), X₁, Y, info)
		panelcross!(Val(false), view(dest, ncols(X₁)+1:size(dest,1), :, :), X₂, Y, info)
		dest
	end
end
function panelcross21(o::StrBootTest, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(ncols(X₁))
		panelcross11(o,X₂,Y,info)
	elseif iszero(ncols(X₂))
		panelcross11(o,X₁,Y,info)
	else
		dest = Matrix{T}(undef, ncols(X₁)+ncols(X₂), length(info))
		panelcross!(Val(false), view(dest,           1:ncols(X₁   ), :), X₁, Y, info)
		panelcross!(Val(false), view(dest, ncols(X₁)+1:size(dest,1), :), X₂, Y, info)
		dest
	end
end
function panelcross12(o::StrBootTest, X::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(ncols(Y₁))
		panelcross11(o,X,Y₂,info)
	elseif iszero(ncols(Y₂))
		panelcross11(o,X,Y₁,info)
	else
		dest = Array{T,3}(undef, ncols(X), length(info), ncols(Y₁)+ncols(Y₂))
		panelcross!(Val(false), view(dest, :, :,           1:ncols(Y₁   )), X, Y₁, info)
		panelcross!(Val(false), view(dest, :, :, ncols(Y₁)+1:size(dest,3)), X, Y₂, info)
		dest
	end
end
function panelcross22(o::StrBootTest, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(ncols(Y₁))
		panelcross21(o,X₁,X₂,Y₂,info)
	elseif iszero(ncols(Y₂))
		panelcross21(o,X₁,X₂,Y₁,info)
	else
		dest = Array{T,3}(undef, ncols(X₁)+ncols(X₂), length(info), ncols(Y₁)+ncols(Y₂))
		panelcross!(Val(false), view(dest, 1:ncols(X₁)             , :, 1:ncols(Y₁)             ), X₁, Y₁, info)
		panelcross!(Val(false), view(dest, ncols(X₁)+1:size(dest,1), :, 1:ncols(Y₁)             ), X₂, Y₁, info)
		panelcross!(Val(false), view(dest, 1:ncols(X₁)             , :, ncols(Y₁)+1:size(dest,3)), X₁, Y₂, info)
		panelcross!(Val(false), view(dest, ncols(X₁)+1:size(dest,1), :, ncols(Y₁)+1:size(dest,3)), X₂, Y₂, info)
		dest
	end
end

# macros to efficiently handle result = input
macro panelsum(o, X, info)
	:(local _X = $(esc(X)); iszero(length($(esc(info)))) || length($(esc(info)))==nrows(_X) ? _X : panelsum($(esc(o)), _X, $(esc(info)) ) )
end
macro panelsum(o, X, wt, info)
	:( panelsum($(esc(o)), $(esc(X)), $(esc(wt)), $(esc(info))) )
end

@inline sumpanelcross(X::Array{T} where T) = dropdims(sum(X, dims=2), dims=2)

# cross-tab sum of a column vector w.r.t. given panel info and fixed-effect var
# one row per FE, one col per other grouping
function crosstabFE(o::StrBootTest{T}, v::AbstractVector{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, o.NFE, nrows(info))
	if o.haswt
		vw = v .* o.sqrtwt
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = @view o._FEID[info[i]]
				vi    = @view      vw[info[i]]
				@inbounds for j in eachindex(vi, FEIDi)
					dest[FEIDi[j],i] += vi[j]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i in eachindex(v,o._FEID)
				dest[o._FEID[i],i] = vw[i]
			end
		end
	else
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = @view o._FEID[info[i]]
				vi    = @view       v[info[i]]
				@inbounds for j in eachindex(vi, FEIDi)
					dest[FEIDi[j],i] += vi[j]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i in eachindex(v,o._FEID)
				dest[o._FEID[i],i] = v[i]
			end
		end
	end		
  dest
end
# same, transposed
function crosstabFEt(o::StrBootTest{T}, v::AbstractVector{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, nrows(info), o.NFE)
	if o.haswt
		vw = v .* o.sqrtwt
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = @view o._FEID[info[i]]
				vi    = @view      vw[info[i]]
				@inbounds for j ∈ eachindex(vi, FEIDi)
					dest[i,FEIDi[j]] += vi[j]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i ∈ eachindex(v,o._FEID)
				dest[i,o._FEID[i]] = vw[i]
			end
		end
	else
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = @view o._FEID[info[i]]
				vi    = @view       v[info[i]]
				@inbounds for j ∈ eachindex(vi, FEIDi)
					dest[i,FEIDi[j]] += vi[j]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i ∈ eachindex(v,o._FEID)
				dest[i,o._FEID[i]] = v[i]
			end
		end
	end	
  dest
end
# crossab handling multiple columns in v
# dimensions: (FEs,entries of info, cols of v)
# this facilitates reshape() to 2-D array in which results for each col of v are stacked vertically
function crosstabFE(o::StrBootTest{T}, v::AbstractMatrix{T}, info::Vector{UnitRange{Int64}}) where T
  dest = zeros(T, o.NFE, nrows(info), ncols(v))
	if o.haswt
		vw = v .* o.sqrtwt
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = view(o._FEID, info[i])
				vi    = @view      vw[info[i],:]
				@inbounds for j ∈ axes(FEIDi,1)
					dest[FEIDi[j],i,:] += @view vi[j,:]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i ∈ axes(o._FEID,1)
				dest[o._FEID[i],i,:] .= @view vw[i,:]
			end
		end
	else
		if nrows(info)>0
			@inbounds Threads.@threads for i ∈ axes(info,1)
				FEIDi = view(o._FEID, info[i])
				vi    = @view       v[info[i],:]
				@inbounds for j ∈ axes(FEIDi,1)
					dest[FEIDi[j],i,:] += @view vi[j,:]
				end
			end
		else  # "robust" case, no clustering
			@inbounds Threads.@threads for i ∈ axes(o._FEID,1)
				dest[o._FEID[i],i,:] .= @view v[i,:]
			end
		end
	end		
  dest
end


# partial any fixed effects out of a data matrix
function partialFE!(o::StrBootTest, In::AbstractArray)
  if length(In)>0
		if o.haswt
			Threads.@threads for f ∈ o.FEs
				tmp = @view In[f.is,:]
				tmp .-= f.sqrtwt * (f.wt'tmp)
			end
		else
			Threads.@threads for f ∈ o.FEs
				tmp = @view In[f.is,:]
				tmp .-= f.wt[1] * sum(tmp; dims=1)
			end
		end
  end
	nothing
end
function partialFE(o::StrBootTest, In::AbstractArray)
  Out = similar(In)
  if length(In)>0
		if o.haswt
			for f ∈ o.FEs
				tmp = @view In[f.is,:]
				Out[f.is,:] .= tmp .- f.sqrtwt *  (f.wt'tmp)
			end
		else
			for f ∈ o.FEs
				tmp = @view In[f.is,:]
				Out[f.is,:] .= tmp .- f.wt[1] * sum(tmp; dims=1)
			end
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

import Base.size
struct FakeArray{N} <: AbstractArray{Bool,N} size::Tuple{Vararg{Int64,N}} end # AbstractArray with almost no storage, just for LinearIndices() conversion         
FakeArray(size...) = FakeArray{length(size)}(size)
size(X::FakeArray) = X.size

import Base.*, Base.adjoint, Base.hcat  # extend * to left- and right-multiply 3-arrays by vec or mat, 2nd index of 3-array corresponds to left and 3rd index to right
@inline *(A::AbstractArray{T,3}, B::AbstractMatrix{T}) where T = reshape(reshape(A, size(A,1) * size(A,2), size(A,3)) * B, size(A,1), size(A,2), size(B,2))
@inline *(A::AbstractArray{T,3}, B::AbstractVector{T}) where T = reshape(reshape(A, size(A,1) * size(A,2), size(A,3)) * B, size(A,1), size(A,2))
@inline *(A::AbstractVecOrMat, B::AbstractArray) = reshape(A * reshape(B, size(B,1), size(B,2) * size(B,3)), size(A,1), size(B,2), size(B,3))
@inline adjoint(A::AbstractArray{T,3} where T) = permutedims(A,(3,2,1))
@inline hcat(A::Array{T,3}, B::Array{T,3}) where T = cat(A,B; dims=3)::Array{T,3}

# 5-argument version of mul!() for 3-arrays but forces α=1, β=0
# import LinearAlgebra.mul!
# @inline mul!(dest::AbstractMatrix{T}, A::AbstractArray{T,3}, B::AbstractVector{T}; α=one(T), β=zero(T)) where T = (mul!(reshape(dest, size(A,1) * size(A,2)), reshape(A, size(A,1) * size(A,2), size(A,3)), B))
