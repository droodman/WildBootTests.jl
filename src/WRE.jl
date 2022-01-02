
# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of 
# Y[:,ind1]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,ind2] for κ constant across replications
# ind1 can be a rowvector
# (only really the Hessian when we narrow Y to Z)
function HessianFixedkappa(o::StrBootTest{T}, ind1::Vector{S} where S<:Integer, ind2::Integer, κ::Number, w::Integer) where T
  dest = Matrix{T}(undef, length(ind1), ncols(o.v))
  @inbounds for i ∈ eachindex(ind1, axes(dest,1))
		_HessianFixedkappa!(o, view(dest,i:i,:), ind1[i], ind2, κ, w)
  end
  dest
end

function _HessianFixedkappa!(o::StrBootTest, dest::AbstractMatrix, ind1::Integer, ind2::Integer, κ::Number, w::Integer)
  if !(o.Repl.Yendog[ind1+1] || o.Repl.Yendog[ind2+1])  # if both vars exog, result = order-0 term only, same for all draws
		!iszero(κ) && 
			coldot!(o, dest, view(o.Repl.XZ,:,ind1), view(o.Repl.invXXXZ,:,ind2))
		if !isone(κ)
			if iszero(κ)
				fill!(dest, o.Repl.YY[ind1+1,ind2+1])
			else
				dest .= κ .* dest .+ (1 - κ) .* o.Repl.YY[ind1+1,ind2+1]
			end
		end
	else
		if !iszero(κ)  # repititiveness in this section to maintain type stability
			if o.Repl.Yendog[ind1+1]
				T1L = o.T1L[isone(o.Nw) || w<o.Nw ? 1 : 2]  # use preallocated destinations
				mul!(T1L, o.S✻UX[ind1+1], o.v)
				if iszero(ind1)
					T1L .+= o.Repl.Xy₁par
				else
					T1L .+= view(o.Repl.XZ,:,ind1)
				end
				if o.Repl.Yendog[ind2+1]
					T1R = o.T1R[isone(o.Nw) || w<o.Nw ? 1 : 2]  # use preallocated destinations
					mul!(T1R, o.S✻UXinvXX[ind2+1], o.v)
					if iszero(ind2)
						T1R .+=  o.Repl.invXXXy₁par
					else
						T1R .+= view(o.Repl.invXXXZ,:,ind2)
					end
					coldot!(o, dest, T1L, T1R)
				else
					coldot!(o, dest, T1L, view(o.Repl.invXXXZ,:,ind2))
				end
			else
				if o.Repl.Yendog[ind2+1]
					T1R = o.T1R[isone(o.Nw) || w<o.Nw ? 1 : 2]  # use preallocated destinations
					mul!(T1R, o.S✻UXinvXX[ind2+1], o.v)
					if iszero(ind2)
						T1R .+=  o.Repl.invXXXy₁par
					else
						T1R .+= view(o.Repl.invXXXZ,:,ind2)
					end
					coldot!(o, dest, view(o.Repl.XZ,:,ind1), T1R)
				else
					dest .= coldot(o, view(o.Repl.XZ,:,ind1), view(o.Repl.invXXXZ,:,ind2))
				end
			end
		end
		if !isone(κ)
			if o.Repl.Yendog[ind1+1]
				T2 = o.S✻UZperpinvZperpZperp[ind1+1]'o.S✻UZperp[ind2+1]  # quadratic term
				T2[diagind(T2)] .-= ind1 ≤ ind2 ? o.S✻UU[ind2+1, ind1+1] : o.S✻UU[ind1+1, ind2+1]  # minus diagonal crosstab
				o.NFE>0 &&
					(T2 .+= o.CTFEU[ind1+1]'(o.invFEwt .* o.CTFEU[ind2+1]))
				if iszero(κ)
					dest .= o.Repl.YY[ind1+1,ind2+1] .+ colquadformminus!(o, (                            @views o.S✻uY[ind2+1][:,ind1+1] .+ o.S✻uY[ind1+1][:,ind2+1])'o.v, T2, o.v)
				else
					dest .=   κ .* dest .+ (1 - κ)   .* colquadformminus!(o, (o.Repl.YY[ind1+1,ind2+1] .+ @views o.S✻uY[ind2+1][:,ind1+1] .+ o.S✻uY[ind1+1][:,ind2+1])'o.v, T2, o.v)
				end
			elseif iszero(κ)
				dest .= o.Repl.YY[ind1+1,ind2+1]
			else
				dest .= κ .* dest .+ (1 - κ) .* o.Repl.YY[ind1+1,ind2+1]
			end
		end
  end
	nothing
end

# put threaded loops in functions to prevent type instability https://discourse.julialang.org/t/type-inference-with-threads/2004/3
function FillingLoop1!(o::StrBootTest{T}, dest::Matrix{T}, ind1::Integer, ind2::Integer, _β̂::AbstractMatrix{T}) where T
	Threads.@threads for i ∈ 1:o.clust[1].N
		PXY✻ = reshape(o.Repl.PXZ[i,ind1], :, 1)
		o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX[ind1+1],i,:)'o.v)

		if iszero(ind2)
			dest[i,:]   = colsum(PXY✻ .* (o.Repl.y₁[i] .- view(o.S✻UMZperp[1],i,:))'o.v)
		elseif o.Repl.Yendog[ind2+1]
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Z[i,ind2] * _β̂ .- view(o.S✻UMZperp[ind2+1],i,:)'β̂v))
		else
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Z[i,ind2] * _β̂))
		end
	end
	nothing
end
function FillingLoop2!(o::StrBootTest{T}, dest::Matrix{T}, ind1::Integer, ind2::Integer, _β̂::AbstractMatrix{T}) where T
	Threads.@threads for i ∈ 1:o.clust[1].N
		S = o.info⋂Data[i]
		PXY✻ = o.Repl.Yendog[ind1+1] ? view(o.Repl.PXZ,S,ind1) .+ view(o.S✻UPX[ind1+1],S,:) * o.v :
																	reshape(o.Repl.PXZ[S,ind1], :, 1)

		if iszero(ind2)
			dest[i,:]   = colsum(PXY✻ .* (o.Repl.y₁[S] .- view(o.S✻UMZperp[1],S,:) * o.v))
		else
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[S,ind2] * _β̂ .- view(o.S✻UMZperp[ind2+1],S,:) * β̂v :
																																o.Repl.Z[S,ind2] * _β̂                                       ))
		end
	end
	nothing
end
function FillingLoop3!(o::StrBootTest{T}, T₁::Matrix{T}, ind1::Integer, ind2::Integer) where T
	Threads.@threads for i ∈ 1:o.N✻
		T₁[o.IDCT⋂✻[i], i] .+= o.SCT⋂uXinvXX[ind2+1,i] * view(o.Repl.XZ,:,ind1)
	end
	nothing
end

# Workhorse for WRE CRVE sandwich filling
# Given a column index within it, ind1, and a matrix β̂s of all the bootstrap estimates, 
# return all bootstrap realizations of P_X * Z[:,ind1]_g ' û₁g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per ⋂ cluster, one col per bootstrap replication
function Filling(o::StrBootTest{T}, ind1::Integer, β̂s::AbstractMatrix) where T
	if o.granular
   	if o.Nw == 1  # create or avoid NxB matrix?
			PXY✻ = reshape(o.Repl.PXZ[:,ind1], :, 1)  # store as matrix to reduce compiler confusion
			o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ o.S✻UPX[ind1+1] * o.v)

			dest = @panelsum(o, PXY✻ .* (o.Repl.y₁ .- o.S✻UMZperp[1] * o.v), o.info⋂Data)

			@inbounds for ind2 ∈ 1:o.Repl.kZ
				_β̂ = view(β̂s,ind2,:)'
				dest .-= @panelsum(o, PXY✻ .* (o.Repl.Yendog[ind2+1] ?  view(o.Repl.Z,:,ind2) * _β̂ .- o.S✻UMZperp[ind2+1] * (o.v .* _β̂) :
															                              (view(o.Repl.Z,:,ind2) * _β̂)                                      ), o.info⋂Data)
			end
		else  # create pieces of each N x B matrix one at a time rather than whole thing at once
			dest = Matrix{T}(undef, o.clust[1].N, ncols(o.v))  # XXX preallocate this & turn Filling into Filling! ?
			@inbounds for ind2 ∈ 0:o.Repl.kZ
				ind2>0 && (β̂v = o.v .* (_β̂ = view(β̂s,ind2,:)'))

				if o.purerobust
					FillingLoop1!(o, dest, ind1, ind2, _β̂)
				else
					FillingLoop2!(o, dest, ind1, ind2, _β̂)
				end
			end
    end
  else  # coarse error clustering
		@inbounds for ind2 ∈ 0:o.Repl.kZ
			β̂v = iszero(ind2) ? o.v : o.v .* (_β̂ = -view(β̂s,ind2,:)')

			# T1 * o.v will be 1st-order terms
			T₁ = o.Repl.Yendog[ind1+1] ? o.Repl.S⋂YX[ind2+1]'o.S✻UXinvXX[ind1+1] : Matrix{T}(undef,0,0)  #  S_∩ (Y_(∥j):*X_∥) (X_∥^' X_∥)^(-1) [S_* (U ̈_(∥i):*X_∥)]^' 

			if o.Repl.Yendog[ind2+1]  # add CT_(⋂,*) (P_(X_par ) Y_(pari).*U ̈_(parj) )
				if o.NClustVar == o.nbootclustvar && iszero(o.subcluster)  # simple case of one clustering: full crosstab is diagonal
					tmp = view(o.Repl.XZ,:,ind1)
					if length(T₁)>0
						T₁[diagind(T₁)] .+=         o.S✻UXinvXX[ind2+1]'tmp
					else
						T₁                = reshape(o.S✻UXinvXX[ind2+1]'tmp, :, 1)  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JN⋂N✻)
					FillingLoop3(o, T₁, ind1, ind2)
				end
				ncols(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.Repl.S⋂PXYZperp[ind1+1] * o.S✻UZperpinvZperpZperp[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
				o.NFE               > 0 && (T₁ = T₁ .- o.Repl.CT_FE⋂PY[ind1+1] * o.CTFEU[ind2+1])
			end

			if iszero(ind2)  # order-0 and -1 terms
				if iszero(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1]
				elseif isone(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1] .+ T₁ .* β̂v
				else
					dest = T₁ * o.v; dest .+= o.Repl.FillingT₀[ind1+1,1]
				end
			else  # y component
				if iszero(ncols(T₁))  # - x*β̂ components
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̂
				elseif isone(ncols(T₁))
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̂ .+ T₁ .* β̂v
				else
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̂
					o.matmulplus!(dest, T₁, β̂v)
				end
			end

			if o.Repl.Yendog[ind1+1] && o.Repl.Yendog[ind2+1]
				for i ∈ 1:o.clust[1].N
					S = o.info⋂Data[i]
					o.colquadformminus!(dest, i, o.v, view(o.S✻UPX[ind1+1],S,:)'view(o.S✻UMZperp[ind2+1],S,:), β̂v)
				end
			end
		end
  end
  dest
end


function PrepWRE!(o::StrBootTest{T}) where T
  EstimateIV!(o.DGP, o, o.null ? [o.r₁ ; o.r] : o.r₁)
  MakeResidualsIV!(o.DGP, o)
  Ü₂par = view(o.DGP.Ü₂ * o.Repl.RparY,:,:)

  @inbounds for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		u = i>0 ? view(Ü₂par,:,i) : view(o.DGP.u⃛₁,:)
		uwt = vHadw(u, o.wt)::Union{Vector{T}, SubArray{T, 1}}

		# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
		o.S✻UX[i+1]      = @panelsum2(o, o.Repl.X₁, o.Repl.X₂, uwt, o.infoBootData)'
		o.S✻UXinvXX[i+1] = o.Repl.invXX * o.S✻UX[i+1]

		if o.LIML || o.bootstrapt || !isone(o.κ)
			o.S✻UZperp[i+1]              = @panelsum(o, o.Repl.Zperp, uwt, o.infoBootData)'
			o.S✻UZperpinvZperpZperp[i+1] = o.Repl.invZperpZperp * o.S✻UZperp[i+1]
			o.NFE>0 && (o.CTFEU[i+1] = crosstabFE(o, uwt, o.infoBootData))
		end

		if o.LIML || !o.robust || !isone(o.κ)
			o.S✻uY[i+1] = @panelsum2(o, o.Repl.y₁par, o.Repl.Z, uwt, o.infoBootData)
			for j ∈ 0:i
				o.S✻UU[i+1,j+1] = vec(@panelsum(o, j>0 ? view(Ü₂par,:,j) : view(o.DGP.u⃛₁,:), uwt, o.infoBootData))
			end
		end

		if o.robust && o.bootstrapt
			if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				for g ∈ 1:o.N✻
					o.SCT⋂uXinvXX[i+1,g] = @panelsum(o, o.Repl.XinvXX, u, o.infoCT⋂✻[g])
				end
			end

			i>0 && (o.S✻UPX[i+1] = o.Repl.XinvXX * o.S✻UX[i+1])
			o.S✻UMZperp[i+1] = o.Repl.Zperp * o.S✻UZperpinvZperpZperp[i+1]

			if iszero(i)  # subtract crosstab of observation by ∩-group of u
				o.S✻UMZperp[  1][o.crosstabBootind] .-= o.DGP.u⃛₁
			else
				o.S✻UMZperp[i+1][o.crosstabBootind] .-= view(Ü₂par,:,i)
			end

			o.NFE>0 &&
				(o.S✻UMZperp[i+1] .+= view(o.invFEwt .* o.CTFEU[i+1], o._FEID, :))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
		end
  end
	nothing
end

function MakeWREStats!(o::StrBootTest{T}, w::Integer) where T
  if isone(o.Repl.kZ)  # optimized code for 1 coefficient in bootstrap regression
		if o.LIML
			YY₁₁   = HessianFixedkappa(o, [0], 0, zero(T), w)  # κ=0 => Y*MZperp*Y
			YY₁₂   = HessianFixedkappa(o, [0], 1, zero(T), w)
			YY₂₂   = HessianFixedkappa(o, [1], 1, zero(T), w)
			YPXY₁₁ = HessianFixedkappa(o, [0], 0, one(T) , w)  # κ=1 => Y*PXpar*Y
			YPXY₁₂ = HessianFixedkappa(o, [0], 1, one(T) , w)
			YPXY₂₂ = HessianFixedkappa(o, [1], 1, one(T) , w)
			YY₁₂YPXY₁₂ = YY₁₂ .* YPXY₁₂
			x₁₁ = YY₂₂ .* YPXY₁₁ .- YY₁₂YPXY₁₂      # elements of YY✻^-1 * YPXY✻ up to factor of det(YY✻)
			x₁₂ = YY₂₂ .* YPXY₁₂ .- YY₁₂ .* YPXY₂₂
			x₂₁ = YY₁₁ .* YPXY₁₂ .- YY₁₂ .* YPXY₁₁
			x₂₂ = YY₁₁ .* YPXY₂₂ .- YY₁₂YPXY₁₂
			κs = (x₁₁ .+ x₂₂)./2; κs = 1 ./ (1 .- (κs .- sqrt.(κs.^2 .- x₁₁ .* x₂₂ .+ x₁₂ .* x₂₁)) ./ (YY₁₁ .* YY₂₂ .- YY₁₂ .* YY₁₂))  # solve quadratic equation for smaller eignenvalue; last term is det(YY✻)
			!iszero(o.Fuller) && (κs .-= o.Fuller / (o._Nobs - o.kX))
			o.As = κs .* (YPXY₂₂ .- YY₂₂) .+ YY₂₂
			_β̂s = (κs .* (YPXY₁₂ .- YY₁₂) .+ YY₁₂) ./ o.As
		else
			o.As = HessianFixedkappa(o, [1], 1, o.κ, w)
			_β̂s  = HessianFixedkappa(o, [1], 0, o.κ, w) ./ o.As
		end

		if o.null
			o.numerw = _β̂s .+ (o.Repl.Rt₁ - o.r) / o.Repl.RRpar
		else
			o.numerw = _β̂s .- o.DGP.β̂₀
			isone(w) && (o.numerw[1] = _β̂s[1] + (o.Repl.Rt₁ - o.r) / o.Repl.RRpar)
		end

		@storeWtGrpResults!(o.numer, o.numerw)

		if o.bootstrapt
			if o.robust
				J⋂s = Filling(o, 1, _β̂s) ./ o.As
				@inbounds for c ∈ 1:o.NErrClustCombs  # sum sandwich over error clusterings
					nrows(o.clust[c].order)>0 && 
						(J⋂s = J⋂s[o.clust[c].order,:])
					@clustAccum!(denom, c, coldot(o, @panelsum(o, J⋂s, o.clust[c].info)))
				end
			else
				denom = (HessianFixedkappa(o, [0], 0, zero(T), w) .- 2 .* _β̂s .* HessianFixedkappa(o, [0], 1, zero(T), w) .+ _β̂s.^2 .* HessianFixedkappa(o, [1], 1, zero(T), w)) ./ o._Nobs ./ o.As  # classical error variance
			end
			@storeWtGrpResults!(o.dist, o.sqrt ? o.numerw ./ sqrt.(denom) : o.numerw .^ 2 ./ denom)
			denom *= o.Repl.RRpar[1]^2
		end

  else  # WRE bootstrap for more than 1 retained coefficient in bootstrap regression

		β̂s = zeros(T, o.Repl.kZ, ncols(o.v))
		A = Vector{Matrix{T}}(undef, ncols(o.v))

		if o.LIML
			YY✻   = [HessianFixedkappa(o, collect(0:i), i, zero(T), w) for i ∈ 0:o.Repl.kZ] # κ=0 => Y*MZperp*Y
			o.YPXY✻ = [HessianFixedkappa(o, collect(0:i), i,  one(T), w) for i ∈ 0:o.Repl.kZ] # κ=1 => Y*PXpar*Y

			@inbounds for b ∈ axes(o.v,2)
				for i ∈ 0:o.Repl.kZ
					o.YY✻_b[1:i+1,i+1]   = YY✻[i+1][:,b]  # fill uppper triangles, which is all that invsym() looks at
					o.YPXY✻_b[1:i+1,i+1] = o.YPXY✻[i+1][:,b]
				end
				o.κ = 1/(1 - eigvals(invsym(o.YY✻_b) * Symmetric(o.YPXY✻_b))[1])
				!iszero(o.Fuller) && (o.κ -= o.Fuller / (o._Nobs - o.kX))
				β̂s[:,b] = (A[b] = invsym(o.κ*o.YPXY✻_b[2:end,2:end] + (1-o.κ)*o.YY✻_b[2:end,2:end])) * (o.κ*o.YPXY✻_b[1,2:end]' + (1-o.κ)*YY✻_b[1,2:end]')
			end
		else
			δnumer =  HessianFixedkappa(o, collect(1:o.Repl.kZ), 0, o.κ, w)
			δdenom = [HessianFixedkappa(o, collect(1:i), i, o.κ, w) for i ∈ 1:o.Repl.kZ]
			
			@inbounds Threads.@threads for b ∈ axes(o.v,2)
				for i ∈ 1:o.Repl.kZ
					o.δdenom_b[1:i,i] = view(δdenom[i],:,b)  # fill uppper triangle
				end
				β̂s[:,b] = (A[b] = invsym(o.δdenom_b)) * view(δnumer,:,b)
			end
		end

		if o.bootstrapt
			if o.robust
				@inbounds for i ∈ 1:o.Repl.kZ  # avoid list comprehension construction because of compiler type detection issue
					o.Zyg[i] = Filling(o, i, β̂s)
				end
			else
				YY✻ = [HessianFixedkappa(o, collect(i:o.Repl.kZ), i, zero(T), w) for i ∈ 0:o.Repl.kZ]  # κ=0 => Y*MZperp*Y
			end
		end

		@inbounds for b ∈ reverse(axes(o.v,2))
			if o.null || w==1 && b==1
				o.numer_b .= o.Repl.RRpar * view(β̂s,:,b) + o.Repl.Rt₁ - o.r
			else
				o.numer_b .= o.Repl.RRpar * (view(β̂s,:,b) - o.DGP.β̂₀)
			end

			if o.bootstrapt
				if o.robust  # Compute denominator for this WRE test stat
					for i ∈ 1:o.Repl.kZ  # XXX replace with 3-D array?
						o._J⋂[:,i] = view(o.Zyg[i],:,b)
					end
					J⋂ = o._J⋂ * (A[b] * o.Repl.RRpar')

					for c ∈ 1:o.NErrClustCombs
						(!isone(o.NClustVar) && nrows(o.clust[c].order)>0) &&
							(J⋂ = J⋂[o.clust[c].order,:])
						J_b = @panelsum(o, J⋂, o.clust[c].info)
						@clustAccum!(denom, c, J_b'J_b)
					end
				else  # non-robust
					for i ∈ 0:o.Repl.kZ
						o.YY✻_b[i+1,i+1:o.Repl.kZ+1] = view(YY✻[i+1],:,b)  # fill upper triangle
					end
					denom = (o.Repl.RRpar * A[b] * o.Repl.RRpar') * [-one(T) ; β̂s[:,b]]'Symmetric(o.YY✻_b) * [-one(T) ; β̂s[:,b]] / o._Nobs  # 2nd half is sig2 of errors
				end
				if o.sqrt
					o.dist[b+first(o.WeightGrp[w])-1] = o.numer_b[1] / sqrt(denom[1])
				else
					o.dist[b+first(o.WeightGrp[w])-1] = o.numer_b'invsym(denom)*o.numer_b  # hand-code for 2-dimensional?
				end
			end
			o.numer[:,b+first(o.WeightGrp[w])-1] = o.numer_b  # slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
		end
	end
	w==1 && o.bootstrapt && (o.statDenom = denom)  # original-sample denominator
	nothing
end
