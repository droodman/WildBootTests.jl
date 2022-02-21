# OLS, ML (score bootstrap), and Anderson-Rubin tests

# Construct stuff that depends linearly or quadratically on r, possibly by interpolation
function MakeInterpolables!(o::StrBootTest{T}) where T
	if o.interpolable
		if iszero(nrows(o.anchor))  # first call? save current r as anchor for interpolation
			o.anchor = o.r
			_MakeInterpolables!(o, o.anchor)
			o.numer₀ = o.numer
			o.interpolate_u && (o.ü₀ = copy(o.ü))
			o.robust && (o.Jcd₀ = deepcopy(o.Jcd))
			return
		end

		if iszero(nrows(o.poles))  # second call: from anchor make set of orthogonal poles, which equal anchor except in one dimension
			o.poles = o.r - o.anchor
			o.robust && (o.denom₀ = deepcopy(o.denom))  # grab quadratic denominator from *previous* (1st) evaluation
			newPole = trues(o.q)  # all poles new
		else  # been here at least twice? interpolate unless current r stretches range > 2X in some dimension(s)
			newPole = abs.(o.r .- o.anchor) .> T(2) .* abs.(o.poles)
		end

		if any(newPole)  # prep interpolation
			@inbounds for h₁ ∈ 1:o.q
				if newPole[h₁]
					o.poles[h₁] = o.r[h₁] - o.anchor[h₁]
					thisr = copy(o.anchor); thisr[h₁] = o.r[h₁]  # if q>1 this creates anchor points that are not graphed, an inefficiency. But simpler to make the deviations from 1st point orthogonal
					_MakeInterpolables!(o, thisr)  # calculate linear stuff at new anchor

					o.∂numer∂r[h₁] = (o.numer .- o.numer₀) ./ o.poles[h₁]
					o.interpolate_u && (o.∂u∂r[h₁] = (o.ü .- o.ü₀) ./ o.poles[h₁])
					if o.robust  # dof > 1 for an ARubin test with >1 instruments.
						for d₁ ∈ 1:o.dof
							for c ∈ 1:o.NErrClustCombs
								o.∂Jcd∂r[h₁][c,d₁] = (o.Jcd[c,d₁] .- o.Jcd₀[c,d₁]) ./ o.poles[h₁]
								for d₂ ∈ 1:d₁
									tmp = coldot(o, o.Jcd₀[c,d₁], o.∂Jcd∂r[h₁][c,d₂])
									d₁ ≠ d₂ && (coldotplus!(Val(o.turbo), tmp, o.Jcd₀[c,d₂], o.∂Jcd∂r[h₁][c,d₁]))  # for diagonal items, faster to just double after the c loop
									@clustAccum!(o.∂denom∂r[h₁][d₁,d₂], c, tmp)
								end
							end
							o.∂denom∂r[h₁][d₁,d₁] .*= T(2)  # double diagonal terms
						end
					end
				end
			end
			if o.robust  # quadratic interaction terms
				@inbounds for h₁ ∈ 1:o.q, h₂ ∈ 1:h₁
					if newPole[h₁] || newPole[h₂]
						for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁, c ∈ 1:o.NErrClustCombs
							@clustAccum!(o.∂²denom∂r²[h₁,h₂][d₁,d₂], c, coldot(o, o.∂Jcd∂r[h₁][c,d₁], o.∂Jcd∂r[h₂][c,d₂]))
						end
					end
				end
			end
			Δ = o.poles
			o.interpolating = true

			if o.q==2  # in this case we haven't yet actually computed interpolables at r, so interpolate them
				o.numerw .= o.numer₀ .+ o.∂numer∂r[1] .* Δ[1] .+ o.∂numer∂r[2] .* Δ[2]
				if o.interpolate_u
					o.ü = o.ü₀ .+ o.∂u∂r[1] .* Δ[1] .+ o.∂u∂r[2] .* Δ[2]
				end
			end

		else  # routine linear interpolation if the anchors not moved
			Δ = o.r - o.anchor
			o.numerw .= o.numer₀ .+ o.∂numer∂r[1] .* Δ[1]
			o.q > 1 && (o.numerw .+= o.∂numer∂r[2] .* Δ[2])
			if o.interpolate_u
				o.ü .= o.ü₀ .+ o.∂u∂r .* Δ[1]
				o.q > 1 && (o.ü .+= o.∂u∂r[2] .* Δ[2])
			end
		end

		if o.robust  # even if an anchor was just moved, and linear components just computed from scratch, do the quadratic interpolation now, from the updated linear factors
			if isone(o.q)
				@inbounds for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁
					o.denom[d₁,d₂] .= o.denom₀[d₁,d₂] .+ o.∂denom∂r[d₁,d₂][1,1] .* Δ .+ o.∂²denom∂r²[d₁,d₂][1,1] .* Δ.^2
				end
			else  # q==2
				@inbounds for d₁ ∈ 1:o.dof, d₂ ∈ 1:d₁
					o.denom[d₁,d₂] .= o.denom₀[d₁,d₂] .+
														o.∂denom∂r[1][d₁,d₂] .* Δ[1] .+
														o.∂denom∂r[2][d₁,d₂] .* Δ[2] .+
														o.∂²denom∂r²[1,1][d₁,d₂] .* (Δ[1] .^ 2) .+
														o.∂²denom∂r²[2,1][d₁,d₂] .* (Δ[1] .* Δ[2]) .+
														o.∂²denom∂r²[2,2][d₁,d₂] .* (Δ[2] .^ 2)
				end
			end
		end
	else  # non-interpolable cases
		_MakeInterpolables!(o, o.r)
	end
	nothing
end

# Construct stuff that depends linearly or quadratically on r and doesn't depend on v. No interpolation.
function _MakeInterpolables!(o::StrBootTest{T}, thisr::AbstractVector) where T
	if o.ML
		o.uXAR = o.sc * (o.AR = o.A * o.R')
	else
		if o.ARubin
			EstimateARubin!(o.DGP, o, thisr)
			MakeResidualsOLSARubin!(o.DGP, o)
		elseif iszero(o.κ)  # regular OLS
			EstimateOLS!(o.DGP, o.null ? [o.r₁ ; thisr] : o.r₁)
			MakeResidualsOLSARubin!(o.DGP, o)
		elseif o.null  # in score bootstrap for IV/GMM, if imposing null, then DGP constraints, κ, Hessian, etc. do vary with r and must be set now
			EstimateIV!(o.DGP, o, [o.r₁ ; thisr])
			InitTestDenoms!(o.DGP, o)
			MakeResidualsIV!(o.DGP, o)
		end

		o.ü = o.DGP.ü₁

		(o.scorebs || (o.robust && o.granular < o.NErrClustCombs)) &&
			(o.uXAR = o.DGP.ü₁ .* o.M.XAR)
	end

	o.SuwtXA = o.scorebs ?
				o.B>0 ?
					 o.NClustVar ?
				      	@panelsum(o, o.uXAR, o.info✻) :
					      o.uXAR                    :
				        sum(o.uXAR,dims=1)  :
			  o.DGP.A * panelsum2(o, o.X₁, o.X₂, o.ü, o.info✻)'  # same calc as in score BS but broken apart to grab intermediate stuff, and assuming residuals defined; X₂ empty except in Anderson-Rubin

	if o.robust && o.bootstrapt && o.granular < o.NErrClustCombs
		u✻XAR = @panelsum(o, o.uXAR, o.info✻⋂)  # collapse data to all-boot && error-cluster-var intersections. If no collapsing needed, panelsum() will still fold in any weights
		if o.B>0
			if o.scorebs
				K = [zeros(T, o.N⋂, o.N✻) for _ in o.dof]::Vector{Matrix{T}}  # inefficient, but not optimizing for the score bootstrap
			else
				K = [panelsum2(o, o.X₁, o.X₂, view(o.DGP.XAR,:,d), o.info⋂) * o.SuwtXA for d ∈ 1:o.dof]::Vector{Matrix{T}}
			end

			o.NFE>0 && !o.FEboot && (o.CT_WE = crosstabFE(o, o.ü, o.info✻))

			if o.NFE>0 && !o.FEboot
				tmp = o.invFEwt .* o.CT_WE
				@inbounds for d ∈ 1:o.dof
					K[d] .+= o.M.CT_XAR[d] * tmp
				end
			end
			@inbounds for d ∈ 1:o.dof
				K[d][o.crosstab⋂✻ind] .-= view(u✻XAR,:,d)  # subtract crosstab of u✻XAR wrt bootstrapping cluster and all-cluster-var intersections from M
				o.scorebs && (K[d] .-= o.ClustShare * colsum(K[d]))  # recenter
			end

			@inbounds for c ∈ 1+o.granular:o.NErrClustCombs
				for d ∈ 1:o.dof
					nrows(o.clust[c].order)>0 &&
						(K[d] = K[d][o.clust[c].order,:])
					o.Kcd[c,d] = @panelsum(o, K[d], o.clust[c].info)
				end
			end
		else  # B = 0. In this case, only 1st term of (64) is non-zero after multiplying by v* (= all 1's), and it is then a one-way sum by c
			o.scorebs &&
				(u✻XAR .-= o.ClustShare * colsum(u✻XAR))  # recenter if OLS
			@inbounds for c ∈ 1:o.NErrClustCombs
				nrows(o.clust[c].order)>0 &&
					(u✻XAR = u✻XAR[o.clust[c].order,:])
				tmp = @panelsum(o, u✻XAR, o.clust[c].info)
				for d ∈ 1:o.dof
					o.Kcd[c,d] = reshape(view(tmp,:,d),:,1)
				end
			end
		end
	end
	MakeNumerAndJ!(o, 1, thisr)  # compute J = κ * v; if Nw > 1, then this is for 1st group; if interpolating, it is only group, and may be needed now to prep interpolation
	nothing
end

# compute stuff depending linearly on v, needed to prep for interpolation
function MakeNumerAndJ!(o::StrBootTest{T}, w::Integer, r::AbstractVector=Vector{T}(undef,0)) where T  # called to *prepare* interpolation, or when w>1, in which case there is no interpolation
	o.numerw = o.scorebs ?
			   (o.B>0 ?
				 	o.SuwtXA'o.v :
				 	o.SuwtXA * o.v_sd    ) :
			   (!o.robust || o.granular || o.purerobust ?
				  	 o.R * (o.β̈dev = o.SuwtXA * o.v) :
				 		(o.R * o.SuwtXA) * o.v)

	if isone(w)
		if o.ARubin
			o.numerw[:,1] = o.v_sd * o.DGP.β̈[o.kX₁+1:end]  # coefficients on excluded instruments in ARubin OLS
		elseif !o.null
			o.numerw[:,1] = o.v_sd * (o.R * (o.ML ? o.β̈ : o.M.β̈) - r)  # Analytical Wald numerator; if imposing null then numer[:,1] already equals this. If not, then it's 0 before this.
		end
	end

	@storeWtGrpResults!(o.numer, o.numerw)

	if o.B>0 && o.robust && o.bootstrapt
		if o.granular || o.purerobust  # optimized treatment when bootstrapping by many/small groups
			if o.purerobust
				o.u✻ = o.ü .* o.v
				o.NFE>0 && partialFE!(o, o.u✻)
				minusX₁₂B(o.u✻, o.X₁, o.X₂, o.β̈dev)
			else  # clusters small but not all singletons
				if o.NFE>0 && !o.FEboot
					o.u✻ = o.ü .* view(o.v, o.ID✻, :)
					partialFE!(o, o.u✻)
					@inbounds for d ∈ 1:o.dof
						o.Jcd[1,d] = @panelsum(o, o.u✻, view(o.M.WXAR,:,d), o.info⋂)                                - panelsum2(o, o.X₁, o.X₂, view(o.M.WXAR,:,d), o.info⋂) * o.β̈dev
					end
				else
					_v = view(o.v,o.ID✻_✻⋂,:)
					@inbounds for d ∈ 1:o.dof
						o.Jcd[1,d] = panelsum(o, panelsum(o, o.ü, view(o.M.WXAR,:,d), o.info✻⋂) .* _v, o.info⋂_✻⋂) - panelsum2(o, o.X₁, o.X₂, view(o.M.WXAR,:,d), o.info⋂) * o.β̈dev
					end
				end
			end
		end
		@inbounds	for c ∈ o.granular+1:o.NErrClustCombs, d ∈ eachindex(axes(o.Jcd, 2), axes(o.Kcd, 2))
			o.Jcd[c,d] = o.Kcd[c,d] * o.v
		end
	end
	nothing
end

function MakeNonWRELoop1!(o::StrBootTest, tmp::Matrix, w::Integer)
	@inbounds Threads.@threads for k ∈ 1:ncols(o.v)
		@inbounds for i ∈ 1:o.dof
			for j ∈ 1:i
				tmp[j,i] = o.denom[i,j][k]  # fill upper triangle, which is all invsym() looks at
			end
		end
		numer_l = view(o.numerw,:,k)
		o.dist[k+first(o.WeightGrp[w])-1] = numer_l'invsym(tmp)*numer_l  # in degenerate cases, cross() would turn cross(.,.) into 0
	end
	nothing
end

function MakeNonWREStats!(o::StrBootTest{T}, w::Integer) where T
	w > 1 && MakeNumerAndJ!(o, w)
	!o.bootstrapt && return

	if o.robust
    if !o.interpolating  # these quadratic computations needed to *prepare* for interpolation but are superseded by interpolation once it is going
    	o.purerobust && (u✻2 = o.u✻ .^ 2)
    	@inbounds for i ∈ 1:o.dof, j ∈ 1:i
    		o.purerobust &&
  	   		(o.denom[i,j] = cross(view(o.M.WXAR,:,i), view(o.M.WXAR,:,j), u✻2) * (o.clust[1].even ? o.clust[1].multiplier : -o.clust[1].multiplier))
				for c ∈ o.purerobust+1:o.NErrClustCombs
					@clustAccum!(o.denom[i,j], c, j==i ? coldot(o, o.Jcd[c,i]) : coldot(o, o.Jcd[c,i],o.Jcd[c,j]))
				end
  		end
   	end

		if isone(o.dof)
			@storeWtGrpResults!(o.dist, o.numerw ./ sqrtNaN.(o.denom[1,1]))
			isone(w) &&
				(o.statDenom = hcat(o.denom[1,1][1]))  # original-sample denominator
		elseif o.dof==2  # hand-code 2D numer'inv(denom)*numer
			t1 = view(o.numerw,1,:)'; t2 = view(o.numerw,2,:)'; t12 = t1.*t2
			@storeWtGrpResults!(o.dist, (t1.^2 .* o.denom[2,2] .- 2 .* t12 .* o.denom[2,1] .+ t2.^2 .* o.denom[1,1]) ./ (o.denom[1,1].*o.denom[2,2] .- o.denom[2,1].^2))
			isone(w) &&
				(o.statDenom = [o.denom[1,1][1] o.denom[2,1][1] ; o.denom[2,1][1] o.denom[2,2][1]])  # original-sample denominator
		else  # build each replication's denominator from vectors that hold values for each position in denominator, all replications
			tmp = Matrix{T}(undef, o.dof, o.dof)
			MakeNonWRELoop1!(o,tmp, w)
			isone(w) && (o.statDenom = tmp)  # original-sample denominator
		end
	else  # non-robust
		AR = o.ML ? o.AR : o.M.AR
		if isone(o.dof)  # optimize for one null constraint
			o.denom[1,1] = o.R * AR
			if !o.ML
				o.u✻ = o.B>0 ? o.v .* o.ü : reshape(o.ü,:,1)  # reshape for type stability
				if o.scorebs
					if o.haswt  # Center variance if interpolated
						o.u✻ .-= o.ClustShare'o.u✻
					else
						o.u✻ .-= colsum(o.u✻) * o.ClustShare  # Center variance if interpolated
					end
				else
					minusX₁₂B(o.u✻, o.X₁, o.X₂, o.β̈dev)  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
				end
				o.denom[1,1] .*= coldot(o,o.u✻)
			end

			@storeWtGrpResults!(o.dist, o.numerw ./ sqrtNaN.(o.denom[1,1]))
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
				if o.B>0
					for k ∈ 1:ncols(o.v)
						numer_l = view(o.numerw,:,k)
						o.dist[k+first(o.WeightGrp[w])-1] = o.numer_l'invdenom*numer_l
						o.u✻ = view(o.v,:,k) .* o.ü
						if o.scorebs
							o.u✻ .-= colsum(o.u✻) * o.ClustShare
						else
							minusX₁₂B(o.u✻, o.X₁, o.X₂, view(o.β̈dev,:,k))  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
						end
						o.dist[k+first(o.WeightGrp[w])-1] ./= (tmp = Symetric(o.u✻'o.u✻))
					end
				else
					o.dist[1] = o.numerw'invdenom*o.numerw
					o.u✻ = o.ü
					if o.scorebs
						o.u✻ .-= colsum(o.u✻) * o.ClustShare  # Center variance if interpolated
					else
						minusX₁₂B(o.u✻, o.X₁, o.X₂, o.β̈dev)  # residuals of wild bootstrap regression are the wildized residuals after partialling out X (or XS) (Kline && Santos eq (11))
					end
					o.dist[1] /= (tmp = Symmetric(o.u✻'o.u✻))
				end
				isone(w) && (o.statDenom = o.denom[1,1] * tmp)  # original-sample denominator
			end
		end
	end
	nothing
end
