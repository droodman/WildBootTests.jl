# Workhorse for WRE CRVE sandwich filling
# Return all bootstrap realizations of P_X * Z_g ' û₁g^*b
# for all groups in the intersection of all error clusterings
# return value is Ncap x kZ x ncols(o.v)
function Filling2(o::StrBootTest{T}, βs::AbstractMatrix) where T
	uwt = o.wt==I ? [o.DGP.u⃛₁ o.Ü₂par] : [o.DGP.u⃛₁ o.Ü₂par] .* o.wt

	# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
	S✻UX      = @panelsum2(o.Repl.X₁, o.Repl.X₂, uwt, o.infoBootData)  # *** un-transposed
	S✻UXinvXX = S✻UX * o.Repl.invXX    # *** un-transposed

	if o.LIML || o.bootstrapt || !isone(o.κ)
		S✻UZperp              = @panelsum(o.Repl.Zperp, uwt, o.infoBootData)  # *** un-transposed
		S✻UZperpinvZperpZperp[i+1] = S✻UZperp * o.Repl.invZperpZperp   # *** un-transposed
		o.NFE>0 && (CTFEU = crosstabFEt(o, uwt, o.infoBootData))  # *** newly transposed for lack of 3-D crosstabFE
	end

	if o.LIML || !o.robust || !isone(o.κ)
		S✻uY = @panelsum2(o.Repl.y₁par, o.Repl.Z, uwt, o.infoBootData)
		for j ∈ 0:i
			S✻UU[i+1,j+1] = @panelsum(j>0 ? view(o.Ü₂par,:,j) : o.DGP.u⃛₁, uwt, o.infoBootData)
		end
	end

	if o.robust && o.bootstrapt
		if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
			SCTcapuXinvXX = [@panelsum(o.Repl.XinvXX, uwt, o.infoCTCap✻[g]) for g ∈ 1:o.N✻]
		end

		i>0 && (o.S✻UPX[i+1] = o.Repl.XinvXX * o.S✻UX[i+1])
		o.S✻UMZperp[i+1] = o.Repl.Zperp * o.S✻UZperpinvZperpZperp[i+1]

		if iszero(i)  # subtract crosstab of observation by ∩-group of u
			o.S✻UMZperp[  1][o.crosstabBootind] .-= o.DGP.u⃛₁
		else
			o.S✻UMZperp[i+1][o.crosstabBootind] .-= view(o.Ü₂par,:,i)
		end

		o.NFE>0 &&
			(o.S✻UMZperp[i+1] .+= view(o.invFEwt .* o.CTFEU[i+1], o._FEID, :))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
	end

  if o.granular
   	if o.Nw == 1  # create or avoid NxB matrix?
			PXY✻ = reshape(o.Repl.PXZ, size(infoCapData)..., 1)
			PXY✻ = PXY✻ .+ o.S✻UPX[ind1+1] * o.v

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
						PXY✻ = hcat(o.Repl.PXZ[i,ind1])
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
						PXY✻ = view(o.Repl.PXZ,S,ind1)
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
					tmp = view(o.Repl.XZ,:,ind1)
					if length(T₁)>0
						T₁[diagind(T₁)] .+= o.S✻UXinvXX[ind2+1]'tmp
					else
						T₁                = o.S✻UXinvXX[ind2+1]'tmp  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JNcapN✻)
					for i ∈ 1:o.N✻
						T₁[o.IDCTCap✻[i], i] .+= o.SCTcapuXinvXX[ind2+1,i] * view(o.Repl.XZ,:,ind1)
					end
				end
				ncols(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.Repl.ScapPXYZperp[ind1+1] * o.S✻UZperpinvZperpZperp[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
				o.NFE               > 0 && (T₁ = T₁ .- o.Repl.CT_FEcapPY[ind1+1] * o.CTFEU[ind2+1])
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
