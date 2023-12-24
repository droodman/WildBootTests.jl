@inline nrows(X::AbstractArray) = size(X,1)
@inline ncols(X::AbstractArray) = size(X,2)
@inline sqrtNaN(x) = x<0 ? typeof(x)(NaN) : sqrt(x)

# given an a x b scratchpad matrix X and r≤a, create an r x b view onto an _adjacent_ subsection of the data
# avoids allocations by reusing X and speeds memory access by compatifying
@inline rowsubview(X,r) = reshape(view(vec(X),1:r*ncols(X)),r,:)

# similar but for square sub-blocks
@inline squaresubview(X,r) = reshape(view(vec(X),1:r*r),r,:)

# Apply .* to a matrix and row of another matrix
function matbyrow!(dest::Matrix{T}, A::Matrix{T}, B::AbstractMatrix{T}, row::Int) where T
	@tturbo warn_check_args=false for j ∈ eachindex(axes(dest,2), axes(A,2), axes(B,2)), i ∈ eachindex(axes(dest,1), axes(A,1))
		dest[i,j] = A[i,j] * B[row,j]
	end
end

# do X .= Y[:,k], hopefully faster
@inline function fillcols!(X::Matrix{T}, Y::Matrix{T}, k::Int) where T
	@tturbo warn_check_args=false for i ∈ indices((X,Y),1)
		Yᵢₖ = Y[i,k]
		for j ∈ indices(X,2)
			X[i,j] = Yᵢₖ
		end
	end
	nothing
end
@inline function fillcols!(X::Matrix{T}, Y::Array{T,3}, k::Int, l::Int) where T
	@tturbo warn_check_args=false for j ∈ indices(X,2), i ∈ indices((X,Y),(1,3))
		X[i,j] = Y[k,l,i]
	end
	nothing
end


# iszero(nrows(X)) && (return Symmetric(X))
# X, ipiv, info = LinearAlgebra.LAPACK.sytrf!('U', Matrix(X))
# iszero(info) && LinearAlgebra.LAPACK.sytri!('U', X, ipiv)
@inline invsym(X) = inv(_cholesky(X))
@inline invsym!(X) = inv(_cholesky!(X))

eigvalsNaN(X) =	try eigvals(X) catch _ fill(eltype(X)(NaN), size(X)) end

@inline colsum(X::AbstractArray) = iszero(length(X)) ? similar(X, 1, size(X)[2:end]...) : sum(X, dims=1)
@inline colsum(X::AbstractArray{Bool}) = iszero(length(X)) ? Array{Int}(undef, 1, size(X)[2:end]...) : sum(X, dims=1)  # type-stable
@inline rowsum(X::AbstractArray) = vec(sum(X, dims=2))

function X₁₂Bminus!(dest::AbstractVector, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractVector)
	t✻minus!(dest, X₁, view(B,1:size(X₁,2)))
	t✻minus!(dest, X₂, @view B[size(X₁,2)+1:end])
	nothing
end
function X₁₂Bminus!(dest::AbstractVecOrMat, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
	t✻minus!(dest, X₁, view(B,1:size(X₁,2),:))
	t✻minus!(dest, X₂, @view B[size(X₁,2)+1:end,:])
	nothing
end
function X₁₂Bplus!(dest::AbstractVecOrMat, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractVector)
	t✻plus!(dest, X₁, view(B,1:size(X₁,2)))
	t✻plus!(dest, X₂, @view B[size(X₁,2)+1:end])
	nothing
end
function X₁₂Bplus!(dest::AbstractVecOrMat, X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
	t✻plus!(dest, X₁, view(B,1:size(X₁,2),:))
	t✻plus!(dest, X₂, @view B[size(X₁,2)+1:end,:])
	nothing
end
function X₁₂B!(dest::AbstractVecOrMat{T}, X₁::AbstractVecOrMat{T}, X₂::AbstractArray{T}, B::AbstractVecOrMat{T}) where T
	fill!(dest, zero(T))
	X₁₂Bplus!(dest, X₁, X₂, B)
	nothing
end
function X₁₂B(X₁::AbstractVecOrMat, X₂::AbstractArray, B::AbstractMatrix)
	dest = X₁ * view(B,1:size(X₁,2),:)
	length(dest)>0 && length(X₂)>0 && t✻plus!(dest, X₂, B[size(X₁,2)+1:end,:])
	dest
end
function X₁₂B(X₁::AbstractArray, X₂::AbstractArray, B::AbstractVector)
	dest = X₁ * view(B,1:size(X₁,2))
	length(dest)>0 && length(X₂)>0 && t✻plus!(dest, X₂, B[size(X₁,2)+1:end])
	dest
end

function coldot!(dest::AbstractMatrix{T}, row::Integer, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T  # colsum(A .* B)
  @tturbo warn_check_args=false for i ∈ indices((A,B),2)
		destᵢ = zero(T)
		for j ∈ indices((A,B),1)
	  	destᵢ += A[j,i] * B[j,i]
		end
		dest[row,i] = destᵢ
  end
	nothing
end
function coldot(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
  dest = Matrix{T}(undef, 1, size(A,2))
  coldot!(dest, 1, A, B)
  dest
end
coldot(A::AbstractMatrix) = coldot(A, A)
coldot(A::AbstractVector, B::AbstractVector) = [dot(A,B)]

function coldotplus!(dest::AbstractMatrix{T}, row::Integer, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	if !iszero(nrows(A))
	  @tturbo warn_check_args=false for i ∈ indices((A,B),2)
			destᵢ = zero(T)
			for j ∈ indices((A,B),1)
		  	destᵢ += A[j,i] * B[j,i]
			end
			dest[row,i] += destᵢ
	  end
	end
	nothing
end
function coldotplus!(dest::AbstractMatrix{T}, c::T, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	if !iszero(nrows(A))
		@tturbo warn_check_args=false for i ∈ indices((A,B),2)
			destᵢ = zero(T)
			for j ∈ indices((A,B),1)
		  	destᵢ += A[j,i] * B[j,i]
			end
			dest[1,i] += c * destᵢ
	  end
	end
	nothing
end
function coldotplus!(dest::AbstractMatrix{T}, c::T, A::AbstractMatrix{T}) where T
	if !iszero(nrows(A))
		@tturbo warn_check_args=false for i ∈ indices((dest,A),2)
			destᵢ = zero(T)
			for j ∈ indices(A,1)
		  	destᵢ += A[j,i]^2
			end
			dest[1,i] += c * destᵢ
	  end
	end
	nothing
end
function coldot!(dest::AbstractMatrix{T}, c::T, A::AbstractMatrix{T}) where T
	fill!(dest, zero(T))
	coldotplus!(dest, c, A)
	nothing
end
function coldotplus!(dest::AbstractMatrix{T}, row::Integer, A::AbstractMatrix{T}, v::AbstractVector{T}, B::AbstractMatrix{T}) where T
  if !iszero(nrows(A))
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,2),axes(B,2))
			destᵢ = zero(T)
			for j ∈ eachindex(axes(A,1),axes(B,1))
		  	destᵢ += A[j,i] * v[j] * B[j,i]
			end
			dest[row,i] += destᵢ
	  end
	end
	nothing
end

function coldotminus!(dest::AbstractVecOrMat{T}, row::Integer, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
  if !iszero(nrows(A))
		@tturbo warn_check_args=false for i ∈ indices((dest,A,B),2)
			destᵢ = zero(T)
			for j ∈ indices((A,B),1)
		  	destᵢ += A[j,i] * B[j,i]
			end
			dest[row,i] -= destᵢ
	  end
	end
	nothing
end

# [A[i,:]'Q*B[i,:] for i] -> dest. Q doesn't have to be square or symmetric
function rowquadformplus!(dest::AbstractVector{T}, A::AbstractMatrix{T}, Q::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	if !iszero(nrows(A))
	  @tturbo warn_check_args=false for i ∈ indices((A,B),1)
			destᵢ = zero(T)
			for j ∈ indices((A,Q),(2,1)), k ∈ indices((Q,B),2)
	    	destᵢ += A[i,j] * Q[j,k] * B[i,k]
			end
			dest[i] += destᵢ
	  end
	end
  dest
end
function rowquadform(A::AbstractMatrix{T}, Q::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
	dest = zeros(T, size(A,1))
	rowquadformplus!(dest, A, Q, B)
	dest
end

# compute norm of each col of A using quadratic form Q; returns one-row matrix
function colquadform(Q::AbstractMatrix{T}, A::AbstractMatrix{T}) where T
  dest = Matrix{T}(undef, 1, size(A,2))
	if !iszero(nrows(A))
	  @tturbo warn_check_args=false for i ∈ eachindex(axes(A,2),axes(dest,2))
			destᵢ = zero(T)
			for j ∈ eachindex(axes(A,1),axes(Q,2)), k ∈ eachindex(axes(A,1),axes(Q,1))
	    	destᵢ += A[j,i] * Q[k,j] * A[k,i]
			end
			dest[1,i] = destᵢ
	  end
	end
	dest
end

# @tturbo-based matrix multiplication
# accepts views as destination
# no error checking
@inline function t✻!(A::AbstractVecOrMat{T}, B::AbstractVecOrMat{T}, C::AbstractVecOrMat{T}) where T
	mul!(A, B, C)
end
@inline function t✻!(A::AbstractVecOrMat{T}, b::T, C::AbstractVecOrMat{T}) where T
	if !iszero(length(A))
		@tturbo for i ∈ indices((A,C),1), j ∈ indices((A,C),2)
			A[i,j] = b * C[i,j]
		end
	end
	nothing
end
function t✻(A::AbstractVecOrMat{T}, B::AbstractVector{T}) where T
	dest = Vector{T}(undef, size(A,1))
	mul!(dest, A, B)
	dest
end
function t✻(A::AbstractVecOrMat{T}, B::AbstractMatrix{T}) where T
	dest = Matrix{T}(undef, size(A,1), size(B,2))
	mul!(dest, A, B)
	dest
end
function t✻!(A::AbstractMatrix{T}, c::T, B::AbstractVecOrMat{T}, C::AbstractMatrix{T}) where T  # c*B*C -> A
	fill!(A, zero(T))
	t✻plus!(A, c, B, C)
end

function t✻plus!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T}, C::AbstractVector{T}) where T  # add B*C to A in place
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1))
			Aᵢₖ = zero(T)
			for j ∈ eachindex(axes(B,2),axes(C,1))
				Aᵢₖ += B[i,j] * C[j]
			end
			A[i] += Aᵢₖ
		end
	end
	nothing
end
function t✻plus!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T}, C::AbstractMatrix{T}) where T  # add B*C to A in place
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2))
			Aᵢₖ = zero(T)
			for j ∈ eachindex(axes(B,2),axes(C,1))
				Aᵢₖ += B[i,j] * C[j,k]
			end
			A[i,k] += Aᵢₖ
		end
	end
	nothing
end
function t✻plus!(A::AbstractMatrix{T}, c::T, B::AbstractVecOrMat{T}, C::AbstractMatrix{T}) where T  # A += c*B*C
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1)), k ∈ eachindex(axes(A,2), axes(C,2))
			Aᵢₖ = zero(T)
			for j ∈ eachindex(axes(B,2),axes(C,1))
				Aᵢₖ += B[i,j] * C[j,k]
			end
			A[i,k] += c * Aᵢₖ
		end
	end
	nothing
end
function t✻plus!(A::AbstractVector{T}, B::AbstractMatrix{T}, C::AbstractVector{T}) where T  # A += B*C
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1))
			Aᵢ = zero(T)
			for j ∈ eachindex(axes(B,2),C)
				Aᵢ += B[i,j] * C[j]
			end
			A[i] += Aᵢ
		end
	end
	nothing
end
function t✻plus!(A::AbstractVector{T}, c::T, B::AbstractMatrix{T}, C::AbstractVector{T}) where T  # A += B*C
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1))
			Aᵢ = zero(T)
			for j ∈ eachindex(axes(B,2),C)
				Aᵢ += B[i,j] * C[j]
			end
			A[i] += c * Aᵢ
		end
	end
	nothing
end
function t✻minus!(A::AbstractMatrix{T}, B::AbstractVecOrMat{T}, C::AbstractMatrix{T}) where T  # A -= B*C
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ indices((A,B),1), k ∈ indices((A,C),2)
			Aᵢₖ = zero(T)
			for j ∈ indices((B,C),(2,1))
				Aᵢₖ += B[i,j] * C[j,k]
			end
			A[i,k] -= Aᵢₖ	
		end
	end
	nothing
end
function t✻minus!(A::AbstractVecOrMat{T}, c::T, B::AbstractVecOrMat{T}, C::AbstractVecOrMat{T}) where T  # A -= c*B*C
	if length(B)>0 && length(C)>0
		@tturbo warn_check_args=false for i ∈ indices((A,B),1), k ∈ indices((A,C),2)
			Aᵢₖ = zero(T)
			for j ∈ indices((B,C),(2,1))
				Aᵢₖ += B[i,j] * C[j,k]
			end
			A[i,k] -= c * Aᵢₖ	
		end
	end
	nothing
end
function t✻minus!(A::AbstractVector{T}, B::AbstractVecOrMat{T}, C::AbstractVector{T}) where T  # add B*C to A in place
	if length(B)>0 && length(C)>0
			@tturbo warn_check_args=false for i ∈ eachindex(axes(A,1),axes(B,1))
				Aᵢ = zero(T)
				for j ∈ eachindex(axes(B,2),axes(C,1))
					Aᵢ += B[i,j] * C[j]
				end
				A[i] -= Aᵢ
		end
	end
	nothing
end
function t✻minus!(A::AbstractVecOrMat{T}, B::T, C::AbstractVecOrMat{T}) where T  # add B*C to A in place
	if length(C)>0
		@tturbo warn_check_args=false for j ∈ indices((A,C),2), i ∈ indices((A,C),1)
			A[i,j] -= B * C[i,j]
		end
	end
	nothing
end

_cholesky!(X) = cholesky!(Symmetric(X), RowMaximum(), check=false)
_cholesky(X)  = cholesky( Symmetric(X), RowMaximum(), check=false)

_cholesky(X::Array{T,3}) where T = [_cholesky(view(X,:,g,:)) for g ∈ eachindex(axes(X,2))]
_cholesky!(X::Array{T,3}) where T = [_cholesky!(view(X,:,g,:)) for g ∈ eachindex(axes(X,2))]


# ch \ Y => Y
function cholldiv!(ch::CholeskyPivoted{T}, Y::AbstractVecOrMat{T}) where T  # adapted from GLM, https://github.com/JuliaStats/GLM.jl/blob/afbb5130ab2773c4b72a3efb4737cf6c6f0c1b09/src/linpred.jl#L134C1-L134C30
  if !iszero(length(Y))
		rnk = rank(ch)
	  len = size(Y,1)
	  if rnk == len
	    ldiv!(ch, Y)
	  else
			invpiv = invperm(ch.piv)
	    for v ∈ eachcol(Y)
	      permute!(v, ch.piv)
	      v[rnk+1:len] .= zero(T)
	      LAPACK.potrs!(ch.uplo, view(ch.factors, 1:rnk, 1:rnk), view(v, 1:rnk, :))
	      permute!(v, invpiv)
	    end
	  end
	end
  Y
end
cholldiv!(ch::CholeskyPivoted{T}, Y::Array{T,3}) where T = cholldiv!(ch, reshape(Y,size(Y,1),:))

# ch \ Y => A
cholldiv!(A, ch::CholeskyPivoted{T}, Y) where T = cholldiv!(ch, (A .= Y))
function cholldiv!(A::AbstractArray{T,3}, vch::Vector{S} where S<:CholeskyPivoted{T}, Y::AbstractMatrix{T}) where T
	@inbounds for g ∈ eachindex(axes(A,2))
		A[:,g,:] .= Y
		cholldiv!(vch[g], view(A,:,g,:))
	end
	A
end
# ch \ Y => Y
function cholldiv!(vch::Vector{S} where S<:CholeskyPivoted{T}, Y::AbstractArray{T,3}) where T
	@inbounds for g ∈ eachindex(axes(Y,2))
		cholldiv!(vch[g], view(Y,:,g,:))
	end
	Y
end


# ch \ Y => copy(Y)
cholldiv(ch::CholeskyPivoted{T}, Y::AbstractVecOrMat{T}) where T = cholldiv!(ch, copy(Y))
cholldiv(ch::CholeskyPivoted{T}, Y::AbstractArray{T,3} ) where T = reshape(cholldiv(ch, reshape(Y,size(Y,1),:)), size(Y))
function cholldiv(vch::Vector{S} where S<:CholeskyPivoted{T}, Y::AbstractArray{T,3}) where T
	dest = copy(Y)
	cholldiv!(vch, dest)
end
function cholldiv(vch::Vector{S} where S<:CholeskyPivoted{T}, Y::AbstractVecOrMat{T}) where T
	dest = Array{T,3}(undef, size(Y,1), length(vch), size(Y,2))
	cholldiv!(dest, vch, Y)
end

cholldiv!!(X::AbstractMatrix{T}, Y::AbstractVecOrMat{T}) where T = cholldiv!(_cholesky!(X),Y)  # overwrites both args


# like Mata panelsetup() but can group on multiple columns, like sort(). But doesn't take minobs, maxobs arguments.
function panelsetup(X::AbstractArray{S} where S, colinds::AbstractVector{T} where T<:Integer)
  N = nrows(X)
	iszero(N) && return(Vector{UnitRange{Int64}}(undef,0))
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
  ID   = Vector{Int64}(undef,N)
	ID[1] = lo = p = 1
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

function panelsum!(dest::AbstractVecOrMat, X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	iszero(length(dest)) && return
	J = CartesianIndices(axes(X)[2:end])
	eachindexJ = eachindex(J)
	@inbounds for g in eachindex(info)
		f, l = first(info[g]), last(info[g])
		fl = f+1:l
		if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj]
				@tturbo warn_check_args=false for i ∈ fl
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
function panelsum!(dest::AbstractVecOrMat{T}, X::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	iszero(length(dest)) && return
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
				@tturbo warn_check_args=false for i ∈ fl
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
# like above but subtract from destination
function panelsumminus!(dest::AbstractVecOrMat{T}, X::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	iszero(length(dest)) && return
	if iszero(length(info)) || nrows(info)==nrows(X)
		dest .-= X .* wt
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
				@tturbo warn_check_args=false for i ∈ fl
					tmp += X[i,Jj] * wt[i]
				end
				dest[g,Jj] -= tmp
			end
		else
			for j ∈ eachindexJ
				dest[g,J[j]] -= X[f,J[j]] * _wt
			end
		end
	end
end
function panelsum!(dest::AbstractArray, X::AbstractArray{T,3} where T, info::Vector{UnitRange{S}} where S<:Integer)
  iszero(length(dest)) && return
  @inbounds for g in eachindex(info)
    f, l = first(info[g]), last(info[g])
    fl = f+1:l
    if f<l
      for i ∈ axes(X,1), k ∈ axes(X,3)
        tmp = X[i,f,k]
        @tturbo warn_check_args=false for j ∈ fl
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
# 1st dimension of result corresponds to columns of X, second to rows of info, third to columns of Y
function panelcross!(dest::AbstractArray{T,3}, X::AbstractVecOrMat{T}, Y::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	iszero(length(dest)) && return
	if iszero(nrows(info)) || nrows(info)==nrows(X)
		if isone(size(X,2))
			if isone(size(Y,2))
				dest[1,:,1] .= X .* Y
			else
				@tturbo warn_check_args=false for k ∈ eachindex(axes(dest,3),axes(Y,2)), j ∈ eachindex(axes(dest,2),axes(X,1),axes(Y,1))
					dest[1,j,k] = X[j] * Y[j,k]
				end
			end
		else
			indicesᵢ = eachindex(axes(dest,1),axes(X,2))
			for k ∈ eachindex(axes(dest,3),axes(Y,2)), j ∈ eachindex(axes(dest,2),axes(X,1),axes(Y,1))
				Yⱼₖ = Y[j,k]
				@tturbo warn_check_args=false for i ∈ indicesᵢ
					dest[i,j,k] = X[j,i] * Yⱼₖ
				end
			end
		end
	else
		if isone(size(Y,2))  # handle special cases to prevent @tturbo crash
			if isone(size(X,2))
				@inbounds for (g, infog) ∈ enumerate(info)
						dest_jgk = zero(T)
					@tturbo warn_check_args=false for i ∈ infog
						dest_jgk += X[i] * Y[i]
					end
					dest[1,g,1] = dest_jgk
				end
			else
				indicesⱼ = eachindex(axes(dest,1),axes(X,2))
				@inbounds for (g, infog) ∈ enumerate(info)
					@tturbo warn_check_args=false for j ∈ indicesⱼ
						dest_jgk = zero(T)
						for i ∈ infog
							dest_jgk += X[i,j] * Y[i]
						end
						dest[j,g,1] = dest_jgk
					end
				end
			end
		elseif isone(size(X,2))
			indicesₖ = eachindex(axes(dest,3),axes(Y,2))
			@inbounds for (g, infog) ∈ enumerate(info)
				@tturbo warn_check_args=false for k ∈ indicesₖ
					dest_jgk = zero(T)
					for i ∈ infog
						dest_jgk += X[i] * Y[i,k]
					end
					dest[1,g,k] = dest_jgk
				end
			end
		else
			indicesⱼ = eachindex(axes(dest,1),axes(X,2)); indiciesₖ = eachindex(axes(dest,3),axes(Y,2))
			@inbounds for (g, infog) ∈ enumerate(info)
				@tturbo warn_check_args=false for j ∈ indicesⱼ, k ∈ indiciesₖ
					dest_jgk = zero(T)
					for i ∈ infog
						dest_jgk += X[i,j] * Y[i,k]
					end
					dest[j,g,k] = dest_jgk
				end
			end
		end
  end
	nothing
end
# version for two matrices on left
function panelcross21!(dest::AbstractArray{T,3}, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	panelcross!(view(dest,            1:size(X₁  ,2),:,:), X₁, Y, info)
	panelcross!(view(dest, size(X₁,2)+1:size(dest,1),:,:), X₂, Y, info)
end
# version for two matrices on right
function panelcross12!(dest::AbstractArray{T,3}, X::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	panelcross!((@view dest[:,:,1:size(Y₁,2)    ]), X, Y₁, info)
	panelcross!((@view dest[:,:,1+size(Y₁,2):end]), X, Y₂, info)
end
function panelcross12(X::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X,2), iszero(length(info)) ? nrows(X) : length(info), size(Y₁,2)+size(Y₂,2))
	panelcross12!(dest, X, Y₁, Y₂, info)
	dest
end

# version for two on left and two on right
function panelcross22!(dest::AbstractArray{T,3}, X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	panelcross!((@view dest[           1:size(X₁,2),:,1:size(Y₁,2)    ]), X₁, Y₁, info)
	panelcross!((@view dest[size(X₁,2)+1:end       ,:,1:size(Y₁,2)    ]), X₂, Y₁, info)
	panelcross!((@view dest[           1:size(X₁,2),:,1+size(Y₁,2):end]), X₁, Y₂, info)
	panelcross!((@view dest[size(X₁,2)+1:end       ,:,1+size(Y₁,2):end]), X₂, Y₂, info)
end
function panelcross(X::AbstractVecOrMat{T}, Y::AbstractVecOrMat{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X,2), iszero(length(info)) ? nrows(X) : length(info), size(Y,2))
	panelcross!(dest, X, Y, info)
	dest
end
function panelcross22(X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, Y₁::AbstractVecOrMat{T}, Y₂::AbstractVecOrMat{T}, info::Vector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X₁,2)+size(X₂,2), iszero(length(info)) ? nrows(X) : length(info), size(Y₁,2)+size(Y₂,2))
	panelcross22!(dest, X₁, X₂, Y₁, Y₂, info)
	dest
end

function panelsum(X::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = isa(X, AbstractVector{T}) ? Vector{T}(undef, iszero(length(info)) ? nrows(X) : length(info)           ) :
		                                 Matrix{T}(undef, iszero(length(info)) ? nrows(X) : length(info), size(X,2))
	if iszero(length(info)) || length(info)==length(X)
		dest .= X .* wt
	else
		panelsum!(dest, X, wt, info)
	end
	dest
end
function panelsum(X::AbstractVecOrMat, info::AbstractVector{UnitRange{T}} where T<:Integer)
	dest = similar(X, iszero(length(info)) ? nrows(X) : length(info), size(X)[2:end]...)
	panelsum!(dest, X, info)
	dest
end
function panelsum(X::AbstractArray{T,3}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	dest = Array{T,3}(undef, size(X,1), size(info,1), size(X,3))
	panelsum!(dest, X, info)
	dest
end
function panelsum2(X₁::AbstractVecOrMat{T}, X₂::AbstractVecOrMat{T}, wt::AbstractVector{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	if iszero(length(X₂))
		panelsum(X₁,wt,info)
	else
		dest = Matrix{T}(undef, iszero(length(info)) ? nrows(X₁) : length(info), ncols(X₁)+ncols(X₂))
		panelsum!(view(dest, :,           1:ncols(X₁   )), X₁, wt, info)
		panelsum!(view(dest, :, ncols(X₁)+1:size(dest,2)), X₂, wt, info)
		dest
	end
end


# macros to efficiently handle result = input
macro panelsum(X, info)
	:(local _X = $(esc(X)); iszero(length($(esc(info)))) || length($(esc(info)))==size(_X,ndims(_X)==3 ? 2 : 1) ? _X : panelsum(_X, $(esc(info)) ) )
end
macro panelsum!(dest, X, info)
	:(local _X = $(esc(X)); iszero(length($(esc(info)))) || length($(esc(info)))==size(_X,ndims(_X)==3 ? 2 : 1) ? $(esc(dest)) .= _X : panelsum!($(esc(dest)), _X, $(esc(info)) ) )
end

@inline sumpanelcross(X::Array{T} where T) = dropdims(sum(X, dims=2); dims=2)

# given two similar matrices, compute their column- and panel-wise dot products
function panelcoldot!(dest::AbstractMatrix{T}, X::AbstractMatrix{T}, Y::AbstractMatrix{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	iszero(length(X)) && return
	if length(info)==nrows(X)
		dest .= X .* Y
		return
	end

	J = CartesianIndices(axes(X)[2:end])
	eachindexJ = eachindex(J)
	@inbounds for g ∈ eachindex(info)
		f, l = first(info[g]), last(info[g])
		fl = f+1:l
		if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj] * Y[f,Jj]
				@tturbo warn_check_args=false for i ∈ fl
					tmp += X[i,Jj] * Y[i,Jj]
				end
				dest[g,Jj] = tmp
			end
		else
			@simd for j ∈ eachindexJ
				dest[g,J[j]] = X[f,J[j]] * Y[f,J[j]]
			end
		end
	end
end
# like panelcoldot! but subtract from instead of overwriting destination
function panelcoldotminus!(dest::AbstractMatrix{T}, X::AbstractMatrix{T}, Y::AbstractMatrix{T}, info::AbstractVector{UnitRange{S}} where S<:Integer) where T
	iszero(length(X)) && return
	if length(info)==nrows(X)
		dest .-= X .* Y
		return
	end

	J = CartesianIndices(axes(X)[2:end])
	eachindexJ = eachindex(J)
	@inbounds for g in eachindex(info)
		f, l = first(info[g]), last(info[g])
		fl = f+1:l
		if f<l
			for j ∈ eachindexJ
				Jj = J[j]
				tmp = X[f,Jj] * Y[f,Jj]
				@tturbo warn_check_args=false for i ∈ fl
					tmp += X[i,Jj] * Y[i,Jj]
				end
				dest[g,Jj] -= tmp
			end
		else
			@simd for j ∈ eachindexJ
				dest[g,J[j]] -= X[f,J[j]] * Y[f,J[j]]
			end
		end
	end
end


# crosstab of  the columns of a vector of matrices given panel info and fixed-effect var
# returns a vector of sparse matrices, one for each col of each matrix
# each returned sparse metric has one row per FE, one col per other grouping
function crosstabFE!(o::StrBootTest{T}, dest::AbstractVector{SparseMatrixCSC{T}}, v::Vector{S} where S<:AbstractVecOrMat{T}, ID::Vector{Int64}, NID::Int64) where T
	i = 1
	if o.haswt
		for M ∈ v
			for Mj ∈ eachcol(M)
				o.Mjw .= Mj .* o.sqrtwt
				dest[i] = sparse(o._FEID, ID, o.Mjw, o.NFE, NID, +)
				i += 1
			end
		end	
	else
		for M ∈ v
			for Mj ∈ eachcol(M)
				dest[i] = sparse(o._FEID, ID, Mj, o.NFE, NID, +)
				i += 1
			end
		end	
	end
	# if iszero(nrows(info)) || nrows(info)==o.Nobs  # "robust" case, no clustering
	# 	@inbounds Threads.@threads for i ∈ eachindex(axes(o._FEID,1),axes(dest,2),axes(vw,1))
	# 		FEIDi = o._FEID[i]
	# 		@inbounds for k ∈ eachindex(axes(dest,3),axes(vw,2))
	# 			dest[FEIDi,i,k] = vw[i,k]
	# 		end
	# 	end
	# else
	# 	fill!(dest, zero(T))
	# 	@inbounds #=Threads.@threads=# for i ∈ eachindex(axes(info,1))
	# 		infoi = info[i]
	# 		@inbounds @fastmath for infoij ∈ infoi
	# 			FEIDij = o._FEID[infoij]
	# 			for k ∈ eachindex(axes(vw,2))
	# 				dest[FEIDij,i,k] += vw[infoij,k]
	# 			end
	# 		end
	# 	end
	# end
  dest
end
function crosstabFE(o::StrBootTest{T}, v::AbstractVecOrMat{T}, ID::Vector{Int64}, NID::Int64) where T
  dest = Vector{SparseMatrixCSC{T}}(undef, ncols(v))
	crosstabFE!(o, dest, [v], ID, NID)
end
function crosstabFE(o::StrBootTest{T}, v₁::AbstractVecOrMat{T}, v₂::AbstractVecOrMat{T}, ID::Vector{Int64}, NID::Int64) where T
  dest = Vector{SparseMatrixCSC{T}}(undef, ncols(v₁)+ncols(v₂))
	crosstabFE!(o, dest, [v₁, v₂], ID, NID)
end


# partial any fixed effects out of a data matrix
function partialFE!(o::StrBootTest{T}, In::AbstractVecOrMat{T}) where T
	if length(In)>0
		is = eachindex(axes(In,1))
		if o.haswt
			@inbounds Threads.@threads for j ∈ eachindex(axes(In,2))
				S = zeros(T,o.NFE)
				for i ∈ is
					S[o._FEID[i]] += In[i,j] * o.FEwt[i]
				end
				for i ∈ is
					In[i,j] -= S[o._FEID[i]] * o.sqrtwt[i]
				end
			end
		else
			@inbounds Threads.@threads for j ∈ eachindex(axes(In,2))
				S = zeros(T,o.NFE)
				for i ∈ is
					S[o._FEID[i]] += In[i,j]
				end
				S .*= o.invsumFEwt
				for i ∈ is
					In[i,j] -= S[o._FEID[i]]
				end
			end
		end
  end
	nothing
end
function partialFE(o::StrBootTest{T}, In::AbstractVecOrMat{T}) where T
  Out = copy(In)
	partialFE!(o, Out)
  Out
end

macro storeWtGrpResults!(dest, content)  # poor hygiene in referencing caller's o and w
  if dest == :(o.dist)
		return quote
			if isone($(esc(:o)).Nw)
				$(esc(dest)) = $(esc(content))
			else
				$(esc(dest))[$(esc(:o)).WeightGrp[$(esc(:w))]] = reshape($(esc(content)),:)
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
	  	  $(esc(X)) .= $(esc(:o)).clust[1].even ? _Y : -_Y
	    else
	  	  $(esc(X)) .= _Y * $(esc(:o)).clust[1].multiplier
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
	    $(esc(X)) .+= _Y .* $(esc(:o)).clust[$(esc(c))].multiplier
	  end
  end
end

import Base.*, Base.adjoint, Base.hcat, Base.vcat, Base.-, Base.size

# use 3-arrays to hold single-indexed sets of matrices. Index in _middle_ dimension.
@inline each(A::Array{T,3}) where T = [view(A,:,i,:) for i ∈ 1:size(A,2)]  #	eachslice(A; dims=2) more elegant but type-unstable
@inline adjoint(A::AbstractArray{T,3} where T) = permutedims(A,(3,2,1))
@inline hcat(A::Array{T,3}, B::Array{T,3}) where T = cat(A,B; dims=3)::Array{T,3}
@inline vcat(A::Array{T,3}, B::Array{T,3}) where T = cat(A,B; dims=1)::Array{T,3}

# copy all the upper triangles to lower
function symmetrize!(A::AbstractArray{T,3}) where T
	d = size(A,1)
	@inbounds for i ∈ 1:d-1, k ∈ axes(A,2), j ∈ i+1:d
		A[j,k,i] = A[i,k,j]
	end
	nothing
end

function *(A::AbstractArray{T,3}, B::AbstractVecOrMat{T}) where T
	dest = zeros(T, size(A,1), size(A,2), size(B,2))
	t✻plus!(dest, A, B)
	dest
end
function *(A::AbstractArray{T,3}, B::AbstractArray{T,3}) where T
	dest = zeros(T, size(A,1), size(A,2), size(B,3))
	t✻plus!(dest, A, B)
	dest
end
function *(A::AbstractVecOrMat{T}, B::AbstractArray{T,3}) where T
	dest = zeros(T, size(A,1), size(B,2), size(B,3))
	t✻plus!(dest, A, B)
	dest
end
function t✻!(dest::AbstractArray{T}, A::AbstractArray{T}, B::AbstractArray{T}) where T
	fill!(dest, zero(T))
	t✻plus!(dest, A, B)
	nothing
end
function t✻!(dest::AbstractArray{T}, c::T, A::AbstractArray{T}, B::AbstractArray{T}) where T
	fill!(dest, zero(T))
	t✻plus!(dest, c, A, B)
	nothing
end
function t✻plus!(dest::AbstractArray{T,3}, A::AbstractArray{T,3}, B::AbstractVecOrMat{T}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,2),axes(dest,3)), j ∈ eachindex(axes(A,1),axes(dest,1)), g ∈ eachindex(axes(A,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,3),axes(B,1))
				dest_jgi += A[j,g,k] * B[k,i]
			end
			dest[j,g,i] += dest_jgi
		end
	end
	nothing
end
function t✻plus!(dest::AbstractArray{T,3}, A::AbstractVecOrMat{T}, B::AbstractArray{T,3}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,3),axes(dest,3)), j ∈ eachindex(axes(A,1), axes(dest,1)), g ∈ eachindex(axes(B,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,2),axes(B,1))
				dest_jgi += A[j,k] * B[k,g,i]
			end
			dest[j,g,i] += dest_jgi
		end
	end
	nothing
end
function t✻plus!(dest::AbstractArray{T,3}, A::AbstractArray{T,3}, B::AbstractArray{T,3}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,3),axes(dest,3)), j ∈ eachindex(axes(A,1),axes(dest,1)), g ∈ eachindex(axes(A,2),axes(B,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,3),axes(B,1))
				dest_jgi += A[j,g,k] * B[k,g,i]
			end
			dest[j,g,i] += dest_jgi
		end
	end
	nothing
end
function t✻plus!(dest::AbstractArray{T,3}, c::T, A::AbstractArray{T,3}, B::AbstractArray{T,3}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,3),axes(dest,3)), j ∈ eachindex(axes(A,1),axes(dest,1)), g ∈ eachindex(axes(A,2),axes(B,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,3),axes(B,1))
				dest_jgi += A[j,g,k] * B[k,g,i]
			end
			dest[j,g,i] += c * dest_jgi
		end
	end
	nothing
end
function t✻minus!(dest::AbstractArray{T,3}, A::AbstractArray{T,3}, B::AbstractVecOrMat{T}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ indices((B,dest), (2,3)), j ∈ indices((A,dest),1), g ∈ indices((A,dest),2)
			dest_jgi = zero(T)
			for k ∈ indices((A,B),(3,1))
				dest_jgi += A[j,g,k] * B[k,i]
			end
			dest[j,g,i] -= dest_jgi
		end
	end
	nothing
end
function t✻minus!(dest::AbstractArray{T,3}, A::AbstractVecOrMat{T}, B::AbstractArray{T,3}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,3),axes(dest,3)), j ∈ eachindex(axes(A,1), axes(dest,1)), g ∈ eachindex(axes(B,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,2),axes(B,1))
				dest_jgi += A[j,k] * B[k,g,i]
			end
			dest[j,g,i] -= dest_jgi
		end
	end
	nothing
end
function t✻minus!(dest::AbstractArray{T,3}, A::AbstractArray{T,3}, B::AbstractArray{T,3}) where T
	if length(dest)>0 && !iszero(length(A)) && length(B)>0
		@tturbo warn_check_args=false for i ∈ eachindex(axes(B,3),axes(dest,3)), j ∈ eachindex(axes(A,1), axes(dest,1)), g ∈ eachindex(axes(A,2),axes(B,2),axes(dest,2))
			dest_jgi = zero(T)
			for k ∈ eachindex(axes(A,3),axes(B,1))
				dest_jgi += A[j,g,k] * B[k,g,i]
			end
			dest[j,g,i] -= dest_jgi
		end
	end
	nothing
end

# in-place inverse of a set of symmetric matrices
function invsym(A::Array{T,3}) where T
	dest = similar(A)
	@inbounds for g ∈ eachindex(axes(A,2))
		dest[:,g,:] = invsym(@view A[:,g,:])
	end
	dest
end

@inline (-)(A::AbstractMatrix{T}, B::Array{T,3}) where T = reshape(A, (size(A,1),1,size(A,2))) .- B  # would be better to overload .-, but more complicated
@inline (-)(B::Array{T,3}, A::AbstractMatrix{T}) where T = B .- reshape(A, (size(A,1),1,size(A,2)))

# delete-g inner products of two vector/matrices; returns full inner product too 
function crossjk(A::VecOrMat{T}, B::AbstractMatrix{T}, info::Vector{UnitRange{Int64}}) where T
	t = panelcross(A,B,info)
	sumt = sum(t; dims=2)  # A'B
	t .= sumt .- t
	(dropdims(sumt; dims=2), t)
end
# version with two data matrices on right
function crossjk(A::VecOrMat{T}, B₁::AbstractMatrix{T}, B₂::AbstractMatrix{T}, info::Vector{UnitRange{Int64}}) where T
	t = panelcross12(A,B₁,B₂,info)
	sumt = sum(t; dims=2)  # A'[B₁ B₂]
	t .= sumt .- t
	(dropdims(sumt; dims=2), t)
end
function crossjk(A::VecOrMat{T}, B::Vector{T}, info::Vector{UnitRange{Int64}}) where T
	(sumt, t) = crossjk(A, view(B,:,:), info)
	(vec(sumt), t)
end

# given data matrix X and cluster-defining info vector, compute X'X, cholesky(X'X), and delete-g invsym(X'X)'s efficiently
function invsymcrossjk(X::Matrix{T}, info::Vector{UnitRange{Int64}}) where T
	SXX =  panelcross(X,X,info)
	XX = sumpanelcross(SXX)
	cholXX = _cholesky(XX)
	invXX = inv(cholXX)
	for (g,S) ∈ enumerate(info)
		Xg = view(X,S,:)
		if size(Xg,1) > size(Xg,2)
			SXX[:,g,:] = invsym(XX - Xg'Xg)
		else
			invXXXg = cholldiv(cholXX, Xg')
			tmp = Xg * invXXXg; tmp -= I
			SXX[:,g,:] = invXX; t✻minus!(view(SXX,:,g,:), invXXXg, cholldiv(_cholesky!(tmp), invXXXg'))
		end
	end
	(XX, cholXX, SXX)
end

# Partial Zperp from A, jackknifed. A and Z are data matrices/vectors. ZZZA is a 3-array
# Returns {A_g - Z_g * ZZZA_g} stacked in A
function partialjk!(A::AbstractVecOrMat{T}, Z::AbstractMatrix{T}, ZZZA::AbstractArray{T}, info::Vector{UnitRange{Int64}}) where T
	if !iszero(length(A))
		indicesᵢ = indices((A,ZZZA),(2,3))
		for (g,G) ∈ enumerate(info)
	    @tturbo warn_check_args=false for i ∈ indicesᵢ, j ∈ G
				destⱼᵢ = zero(T)
				for k ∈ indices((Z,ZZZA), (2,1))
			  	destⱼᵢ += Z[j,k] * ZZZA[k,g,i]
				end
				A[j,i] -= destⱼᵢ
	    end
		end
	end
  A
end
function partialjk(A::AbstractVecOrMat{T}, Z::AbstractMatrix{T}, ZZZA::AbstractArray{T}, info::Vector{UnitRange{Int64}}) where T
	dest = copy(A)
	partialjk!(dest,Z,ZZZA,info)
end
function partialjk(A₁::AbstractVecOrMat{T}, A₂::AbstractVecOrMat{T}, Z::AbstractMatrix{T}, ZZZA::AbstractArray{T}, info::Vector{UnitRange{Int64}}) where T
	dest = [A₁ A₂]
	partialjk!(dest,Z,ZZZA,info)
end

# multiply an a-vector of b x c sparse matrices, with shared sparsity pattern, by an a x d matrix, producing a d-vector of b x c matrices, still same sparsity pattern
function t✻plus!(dest::AbstractVector{S}, A::AbstractVector{S}, B::AbstractVecOrMat{T}) where S<:SparseMatrixCSC{T} where T
	@inbounds @fastmath for i ∈ eachindex(axes(A,1),axes(B,1)), j ∈ eachindex(axes(dest,1), axes(B,2))
		dest[j].nzval .+= A[i].nzval .* B[i,j]
	end
	dest
end
function t✻minus!(dest::AbstractVector{S}, A::AbstractVector{S}, B::AbstractVecOrMat{T}) where S<:SparseMatrixCSC{T} where T
	@inbounds @fastmath for i ∈ eachindex(axes(A,1),axes(B,1)), j ∈ eachindex(axes(dest,1), axes(B,2))
		dest[j].nzval .-= A[i].nzval .* B[i,j]
	end
	dest
end
function t✻!(dest::AbstractVector{S}, A::AbstractVector{S}, B::AbstractVecOrMat{T}) where S<:SparseMatrixCSC{T} where T
	for desti ∈ dest
		fill!(desti.nzval, zero(T))
	end
	t✻plus!(dest, A, B)
end
# import Base.*
# function (*)(A::Vector{SparseMatrixCSC{Tv,Ti}}, B::Matrix{Tv}) where Tv where Ti
# 	dest = [SparseMatrixCSC{Tv,Ti}(A[1].m, A[1].n, A[1].colptr, A[1].rowval, zeros(Tv,size(A[1].nzval))) for _ ∈ 1:size(B,2)]
# 	t✻plus!(dest, A, B)
# end

# do B - A -> A efficiently for sparse matrices with same sparsity pattern
function lsub!(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}) where Tv where Ti
	A.nzval .= B.nzval .- A.nzval
	A
end
function lsub!(A::AbstractVector{S}, B::Vector{S}) where S<:SparseMatrixCSC
	for i ∈ eachindex(axes(A,1), axes(B,1))
		A[i].nzval .= B[i].nzval .- A[i].nzval
	end
	A
end

# multiply a sparse matrix by a scalar in place without changing sparsity pattern even if some results are 0
function rmul!(A::SparseMatrixCSC, b::Real)
	A.nzval .*= b
	A
end
