
# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, ind1, ind2, return all bootstrap realizations of 
# Y[:,ind1]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,ind2] for κ constant across replications
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
				T1L .= o.S✻UX[ind1+1] * o.v
				T1L .+= _T1L
			else
				T1L = _T1L
			end
			_T1R = iszero(ind2) ? o.Repl.invXXXy₁par : @view o.Repl.invXXXZ[:,ind2]
			if o.Repl.Yendog[ind2+1]
				T1R = o.T1R[isone(o.Nw) || w<o.Nw ? 1 : 2]
				T1R .= o.S✻UXinvXX[ind2+1] * o.v
				T1R .+= _T1R
			else
				T1R = _T1R  # right-side linear term
			end
			coldot!(dest, T1L, T1R)
		end
		if !isone(κ)
			if o.Repl.Yendog[ind1+1]
				T2 = o.S✻UZperpinvZperpZperp[ind1+1]'o.S✻UZperp[ind2+1]  # quadratic term
				T2[diagind(T2)] .-= ind1 ≤ ind2 ? o.S✻UU[ind2+1, ind1+1] : o.S✻UU[ind1+1, ind2+1]  # minus diagonal crosstab
				o.NFE>0 &&
					(T2 .+= o.CTFEU[ind1+1]'(o.invFEwt .* o.CTFEU[ind2+1]))

				if iszero(κ)
					dest .= o.Repl.YY[ind1+1,ind2+1] .+ colquadformminus!((                            @views o.S✻uY[ind2+1][:,ind1+1] .+ o.S✻uY[ind1+1][:,ind2+1])'o.v, T2, o.v)
				else
					dest .=   κ .* dest .+ (1 - κ)   .* colquadformminus!((o.Repl.YY[ind1+1,ind2+1] .+ @views o.S✻uY[ind2+1][:,ind1+1] .+ o.S✻uY[ind1+1][:,ind2+1])'o.v, T2, o.v)
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
			PXY✻ = iszero(ind1) ? o.Repl.PXy₁ : view(o.Repl.PXZ,:,ind1)
			o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ o.S✻UPX[ind1+1] * o.v)

			dest = @panelsum(PXY✻ .* (o.Repl.y₁ .- o.S✻UMZperp[1] * o.v), o.wt, o.infoCapData)

			for ind2 ∈ 1:o.Repl.kZ
				_β = view(βs,ind2,:)'
				dest .-= @panelsum(PXY✻ .* (o.Repl.Yendog[ind2+1] ? view(o.Repl.Z,:,ind2) * _β .- o.S✻UMZperp[ind2+1] * (o.v .* _β) :
															                              view(o.Repl.Z,:,ind2) * _β                                        ), o.wt, o.infoCapData)
			end
		else  # create pieces of each N x B matrix one at a time rather than whole thing at once
			dest = Matrix{T}(undef, o.clust[1].N, ncols(o.v))  # XXX preallocate this & turn Filling into Filling! ?
			for ind2 ∈ 0:o.Repl.kZ
				ind2>0 && (βv = o.v .* (_β = view(βs,ind2,:)'))

				if o.purerobust
					for i ∈ 1:o.clust[1].N
						PXY✻ = hcat(iszero(ind1) ? o.Repl.PXy₁[i] : o.Repl.PXZ[i,ind1])
						o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX[ind1+1],i,:)'o.v)

						if iszero(ind2)
							dest[i,:]   = wtsum(o.wt, PXY✻ .* (o.Repl.y₁[i] .- view(o.S✻UMZperp[1],i,:))'o.v)
						else
							dest[i,:] .-= wtsum(o.wt, PXY✻ .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[i,ind2] * _β .- view(o.S✻UMZperp[ind2+1],i,:)'βv :
																						                             o.Repl.Z[i,ind2] * _β))
						end
					end
				else
					for i ∈ 1:o.clust[1].N
						S = o.infoCapData[i]
						PXY✻ = ind1 ? o.Repl.PXZ[S,ind1] : o.Repl.PXy₁[S]
						o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX[ind1+1],S,:) * o.v)

						if iszero(ind2)
							dest[i,:]   = wtsum(o.wt, PXY✻ .* (o.Repl.y₁[S] .- view(o.S✻UMZperp[1],S,:) * o.v))
						else
							dest[i,:] .-= wtsum(o.wt, PXY✻ .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[S,ind2] * _β .- view(o.S✻UMZperp[ind2+1],S,:) * βv :
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
			T₁ = o.Repl.Yendog[ind1+1] ? o.Repl.ScapYX[ind2+1] * o.S✻UXinvXX[ind1+1] : Matrix{T}(undef,0,0)  #  S_∩ (Y_(∥j):*X_∥ ) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

			if o.Repl.Yendog[ind2+1]  # add CT_(cap,*) (P_(X_par ) Y_(pari).*U ̈_(parj) )
				if o.NClustVar == o.nbootclustvar && iszero(o.subcluster)  # simple case of one clustering: full crosstab is diagonal
					tmp = ind1>0 ? view(o.Repl.XZ,:,ind1) : o.Repl.Xy₁par
					if length(T₁)>0
						T₁[diagind(T₁)] .+= o.S✻UXinvXX[ind2+1]'tmp
					else
						T₁                = o.S✻UXinvXX[ind2+1]'tmp  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JNcapN✻)
					for i ∈ 1:o.N✻
						T₁[o.IDCTCap✻[i], i] .+= o.SCTcapuXinvXX[ind2+1,i] * (ind1>0 ? view(o.Repl.XZ,:,ind1) : o.Repl.Xy₁par)
					end
				end
				length(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.Repl.ScapPXYZperp[ind1+1] * o.S✻UZperpinvZperpZperp[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
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
					colquadformminus!(dest, i, cross(view(o.S✻UPX[ind1+1],S,:), o.haswt ? o.wt[S] : I, view(o.S✻UMZperp[ind2+1],S,:)), o.v, βv)
				end
			end
		end
  end
  dest
end

function PrepWRE!(o::StrBootTest)
  Estimate!(o.DGP, o.null ? [o.r₁ ; o.r] : o.r₁)
  MakeResiduals!(o.DGP)
  o.Ü₂par = o.DGP.Ü₂ * o.Repl.RparY

  for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		uwt = vHadw(i>0 ? view(o.Ü₂par,:,i) : o.DGP.u⃛₁, o.wt)

		# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
		o.S✻UX[i+1]      = @panelsum2(o.Repl.X₁, o.Repl.X₂, uwt, o.infoBootData)'
		o.S✻UXinvXX[i+1] = o.Repl.invXX * o.S✻UX[i+1]

		if o.LIML || o.bootstrapt || !isone(o.κ)
			o.S✻UZperp[i+1]              = @panelsum(o.Repl.Zperp, uwt, o.infoBootData)'
			o.S✻UZperpinvZperpZperp[i+1] = o.Repl.invZperpZperp * o.S✻UZperp[i+1]
			o.NFE>0 && (o.CTFEU[i+1] = crosstabFE(o, uwt, o.infoBootData))
		end

		if o.LIML || !o.robust || !isone(o.κ)
			o.S✻uY[i+1] = @panelsum2(o.Repl.y₁par, o.Repl.Z, uwt, o.infoBootData)
			for j ∈ 0:i
				o.S✻UU[i+1,j+1] = @panelsum(j>0 ? view(o.Ü₂par,:,j) : o.DGP.u⃛₁, uwt, o.infoBootData)
			end
		end

		if o.robust && o.bootstrapt
			if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				for g ∈ 1:o.N✻
					o.SCTcapuXinvXX[i+1,g] = @panelsum(o.Repl.XinvXX, uwt, o.infoCTCap✻[g])
				end
			end

			o.S✻UPX[i+1]     = o.Repl.XinvXX * o.S✻UX[i+1]
			o.S✻UMZperp[i+1] = o.Repl.Zperp * o.S✻UZperpinvZperpZperp[i+1]
			if o.Nobs == o.N✻  # subtract "crosstab" of observation by cap-group of u
				o.S✻UMZperp[i+1][diagind(o.S✻UMZperp[i+1])] .-= i>0 ? view(o.Ü₂par,:,i) : o.DGP.u⃛₁  # case: bootstrapping by observation
			else
				if iszero(i)
					for g ∈ 1:o.Nobs
						o.S✻UMZperp[i+1][g,o.IDBootData[g]] -= o.DGP.u⃛₁[g]
					end
				else
					for g ∈ 1:o.Nobs
						o.S✻UMZperp[i+1][g,o.IDBootData[g]] -= o.Ü₂par[g,i]
					end
				end
			end
			o.NFE>0 &&
				(o.S✻UMZperp[i+1] .+= view(o.invFEwt .* o.CTFEU[i+1], o._FEID,:))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
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
			x₁₁ = YY₂₂ .* YPXY₁₁ .- YY₁₂YPXY₁₂      # elements of YY✻^-1 * YPXY✻ up to factor of det(YY✻)
			x₁₂ = YY₂₂ .* YPXY₁₂ .- YY₁₂ .* YPXY₂₂
			x₂₁ = YY₁₁ .* YPXY₁₂ .- YY₁₂ .* YPXY₁₁
			x₂₂ = YY₁₁ .* YPXY₂₂ .- YY₁₂YPXY₁₂
			κs = (x₁₁ .+ x₂₂)/2; κs = 1 ./ (1 .- (κs - sqrt.(κs.^2 .- x₁₁.*x₂₂ .+ x₁₂.*x₂₁)) ./ (YY₁₁ .* YY₂₂ .- YY₁₂ .* YY₁₂))  # solve quadratic equation for smaller eignenvalue; last term is det(YY✻)
			!iszero(o.Fuller) && (κs .-= o.Fuller / (o._Nobs - o.kX))
			o.As = κs .* (YPXY₂₂ .- YY₂₂) .+ YY₂₂
			βs = (κs .* (YPXY₁₂ .- YY₁₂) .+ YY₁₂) ./ o.As
		else
			o.As = HessianFixedκ(o, [1], 1, o.κ, w)
			βs   = HessianFixedκ(o, [1], 0, o.κ, w) ./ o.As
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
						(Jcaps = Jcaps[o.clust[c].order,:])  # physically reorder rows to keep @turbo happy
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
			YY✻   = [HessianFixedκ(o, collect(0:i), i, zero(T), w) for i ∈ 0:o.Repl.kZ] # κ=0 => Y*MZperp*Y
			o.YPXY✻ = [HessianFixedκ(o, collect(0:i), i,  one(T), w) for i ∈ 0:o.Repl.kZ] # κ=1 => Y*PXpar*Y

			for b ∈ axes(o.v,2)
				for i ∈ 0:o.Repl.kZ
					o.YY✻_b[1:i+1,i+1]   = YY✻[i+1][:,b]  # fill uppper triangles, which is all that invsym() looks at
					o.YPXY✻_b[1:i+1,i+1] = o.YPXY✻[i+1][:,b]
				end
				o.κ = 1/(1 - eigvals(invsym(o.YY✻_b) * Symmetric(o.YPXY✻_b))[1])
				!iszero(o.Fuller) && (o.κ -= o.Fuller / (o._Nobs - o.kX))
				βs[:,b] = (A[b] = invsym(o.κ*o.YPXY✻_b[2:end,2:end] + (1-o.κ)*o.YY✻_b[2:end,2:end])) * (o.κ*o.YPXY✻_b[1,2:end]' + (1-o.κ)*YY✻_b[1,2:end]')
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
				Zyg = Vector{Matrix{T}}(undef, o.Repl.kZ)  # XXX move to Init()
				for i ∈ 1:o.Repl.kZ  # avoid list comprehension construction because of compiler type detection issue
					Zyg[i] = Filling(o, i, βs)
				end
			else
				o.YY✻ = [HessianFixedκ(o, collect(i:o.Repl.kZ), i, zero(T), w) for i ∈ 0:o.Repl.kZ]  # κ=0 => Y*MZperp*Y
			end
		end

		numer_b = Vector{T}(undef,nrows(o.Repl.RRpar))  # XXX move to Init()
		for b ∈ reverse(axes(o.v,2))
			if o.null || w==1 && b==1
				numer_b .= o.Repl.RRpar * view(βs,:,b) + o.Repl.Rt₁ - o.r
			else
				numer_b .= o.Repl.RRpar * (view(βs,:,b) - o.DGP.β₀)
			end

			if o.bootstrapt
				if o.robust  # Compute denominator for this WRE test stat
					for i ∈ 1:o.Repl.kZ  # XXX replace with 3-D array
						o._Jcap[:,i] = view(Zyg[i],:,b)
					end
					Jcap = o._Jcap * (A[b] * o.Repl.RRpar')

					for c ∈ 1:o.NErrClustCombs
						(!isone(o.NClustVar) && length(o.clust[c].order)>0) &&
							(Jcap = Jcap[o.clust[c].order,:])  # physically reorder rows to keep @turbo happy
						J_b = @panelsum(Jcap, o.clust[c].info)
						@clustAccum!(denom, c, J_b'J_b)
					end
				else  # non-robust
					for i ∈ 0:o.Repl.kZ
						o.YY✻_b[i+1,i+1:o.Repl.kZ+1] = view(YY✻[i+1],b,:)  # fill upper triangle
					end
					denom = (o.Repl.RRpar * A[b] * o.Repl.RRpar') * [-one(T) ; βs[:,b]]'Symmetric(o.YY✻_b) * [-one(T) ; βs[:,b]] / o._Nobs  # 2nd half is sig2 of errors
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
