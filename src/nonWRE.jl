# OLS, ML (score bootstrap), and Anderson-Rubin tests

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
		u✻XAR = @panelsum(o.uXAR, o.wt, o.infoAllData)  # collapse data to all-boot && error-cluster-var intersections. If no collapsing needed, panelsum() will still fold in any weights
		if o.B>0
			if o.scorebs
				Kd = zeros(T, o.clust[1].N, o.dof, o.N✻)  # inefficient, but not optimizing for the score bootstrap
			else
				Kd = @panelsum2(o.X₁, o.X₂, vHadw(o.DGP.XAR, o.wt), o.infoCapData) * o.SuwtXA  # overloaded def of * for >2D arrays
			end

			o.NFE>0 && !o.FEboot && (o.CT_WE = crosstabFE(o, vHadw(o.ü, o.wt), o.infoBootData))

			o.NFE>0 && !o.FEboot &&
				(Kd .+= o.M.CT_XAR * (o.invFEwt .* o.CT_WE))  # overloaded def of * for >2D arrays
			Kd[o.crosstabCap✻ind] .-= reshape(u✻XAR,:)  # subtract crosstab of u✻XAR wrt bootstrapping cluster and all-cluster-var intersections from M
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
				(u✻XAR .-= o.ClustShare * colsum(u✻XAR))  # recenter if OLS
			for c ∈ 1:o.NErrClustCombs
				length(o.clust[c].order)>0 &&
					(u✻XAR = view(u✻XAR, o.clust[c].order,:))
				tmp = @panelsum(u✻XAR, o.clust[c].info)
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
				o.u✻ = o.ü .* o.v
				partialFE!(o, o.u✻)
				o.u✻ -= X₁₂B(o.X₁, o.X₂, o.βdev)  # XXX make X₁₂Bminus
			else  # clusters small but not all singletons
				if o.NFE>0 && !o.FEboot
					o.u✻ = o.ü .* view(o.v, o.IDBootData, :)
					partialFE!(o, o.u✻)
					for d ∈ 1:o.dof
						o.Jcd[1,d] = @panelsum(o.u✻, o.M.WXAR[:,d], o.infoCapData)                           - @panelsum2(o.X₁, o.X₂, o.M.WXAR[:,d], o.infoCapData) * o.βdev
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
      		o.purerobust && (u✻2 = o.u✻ .^ 2)
      		for i ∈ 1:o.dof, j ∈ 1:i
    			o.purerobust &&
  	      			(o.denom[i,j] = cross(view(o.M.WXAR,:,i), view(o.M.WXAR,:,j), u✻2) * (o.clust[1].even ? o.clust[1].multiplier : -o.clust[1].multiplier))
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
				o.u✻ = o.B>0 ? o.v .* o.ü : o.ü
				if o.scorebs
					if o.haswt  # Center variance if interpolated
						o.u✻ .-= o.ClustShare'o.u✻
					else
						o.u✻ .-= colsum(o.u✻) * o.ClustShare  # Center variance if interpolated
					end
				else
					o.u✻ -= X₁₂B(o.X₁, o.X₂, o.βdev)  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
				end
				if o.haswt
					o.denom[1,1] .*= o.wt'(o.u✻ .^ 2)
				else
					o.denom[1,1] .*= coldot(o.u✻)
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
					o.u✻ = o.B>0 ? view(o.v,:,k) .* o.ü : o.ü
					if o.scorebs
						if o.haswt  # Center variance if interpolated
							o.u✻ .-= o.wt'o.u✻ * o.ClustShare
						else
							o.u✻ .-= colsum(o.u✻) * o.ClustShare  # Center variance if interpolated
						end
					else
						o.u✻ .-= X₁₂B(o.X₁, o.X₂, view(o.βdev,:,k))  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
					end
					o.dist[k+first(o.WeightGrp[w])-1] ./= (tmp = symcross(o.u✻, o.wt))
				end
				isone(w) && (o.statDenom = o.denom[1,1] * tmp)  # original-sample denominator
			end
		end
	end
end