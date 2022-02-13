# stuff done once per exucution--not depending on r
function InitWRE!(o::StrBootTest{T}) where T
	o.LIML && o.Repl.kZ==1 && o.Nw==1 && (o.As = o.β̈s = zeros(1, o.B+1))

	o.S✻ZperpU              = [Matrix{T}(undef, ncols(o.Repl.Zperp), o.N✻) for _ ∈ 0:o.Repl.kZ]
	o.invZperpZperpS✻ZperpU = [Matrix{T}(undef, ncols(o.Repl.Zperp), o.N✻) for _ ∈ 0:o.Repl.kZ]
	o.S✻YU                  = [Vector{T}(undef, o.N✻) for _ ∈ 0:o.Repl.kZ, _ ∈ 0:o.Repl.kZ]
	o.S✻XU                  = [Matrix{T}(undef, o.DGP.kX, o.N✻) for _ ∈ 0:o.Repl.kZ]
	o.invXXS✻XU             = [Matrix{T}(undef, o.DGP.kX, o.N✻) for _ ∈ 0:o.Repl.kZ]
	o.S✻UU                  = [Vector{T}(undef, o.N✻) for _ ∈ 0:o.Repl.kZ, _ ∈ 0:o.Repl.kZ]
	o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kX, ncols(o.v))] :
												[Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[end]))]
	o.T1R = deepcopy(o.T1L)
	o.robust && o.bootstrapt && iszero(o.granular) &&
		(o.S✻UPX_S✻UMZperp = [Array{T,3}(undef, o.N✻, nrows(o.info⋂), o.N✻) for _ ∈ 0:o.Repl.kZ, _ ∈ 0:o.Repl.kZ])

	if o.bootstrapt
		o.δdenom_b = zeros(o.Repl.kZ, o.Repl.kZ)
		if o.robust && o.granular
			o.S✻UMZperp = [Matrix{T}(undef, o.Nobs, o.N✻) for _ ∈ 0:o.Repl.kZ]
			o.S✻UPX     = [Matrix{T}(undef, o.Nobs, o.N✻) for _ ∈ 0:o.Repl.kZ]
		end
		o._J⋂ = zeros(o.clust[1].N, o.Repl.kZ)
		!o.granular && !(o.NClustVar == o.nbootclustvar && iszero(o.subcluster)) &&
			(o.S✻⋂u₁XinvXX = [Matrix{T}(undef, length(o.info✻[g]), o.Repl.kX) for _ ∈ 0:o.Repl.kZ, g ∈ 1:o.N✻])
		if o.LIML || !o.robust
			o.YY✻_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
			o.YPXY✻_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
		end
		o.NFE>0 && (o.bootstrapt || !isone(o.κ) || o.LIML) && (o.CTFEU = Vector{Matrix{T}}(undef, o.Repl.kZ+1))
	end

	o.S✻⋂XY₂      = o.Repl.S✻⋂XY₂     - o.Repl.S✻⋂XZperp     * o.Repl.invZperpZperpZperpY₂  - o.Repl.invZperpZperpZperpX' * (o.Repl.S✻⋂ZperpY₂  - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpY₂ )
	o.S✻⋂XX       = o.Repl.S✻⋂XX      - o.Repl.S✻⋂XZperp     * o.Repl.invZperpZperpZperpX   - o.Repl.invZperpZperpZperpX' * (o.Repl.S✻⋂XZperp'  - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpX  )
	o.S✻⋂XDGPZ    = o.DGP.S✻⋂XZpar    - o.Repl.S✻⋂XZperp     * o.DGP.invZperpZperpZperpZpar - o.Repl.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZpar - o.Repl.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZpar)
	o.S✻⋂Xy₁      = o.Repl.S✻⋂Xy₁     - o.Repl.S✻⋂XZperp     * o.Repl.invZperpZperpZperpy₁  - o.Repl.invZperpZperpZperpX' * (o.Repl.S✻⋂Zperpy₁  - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpy₁ )
	o.S✻⋂X_DGPZR₁ = o.DGP.S✻⋂XZR₁     - o.Repl.S✻⋂XZperp     * o.DGP.invZperpZperpZperpZR₁  - o.Repl.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZR₁  - o.Repl.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZR₁ )
	  S✻⋂ZperpX   = o.Repl.S✻⋂XZperp' - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpX

	o.invXXS✻⋂XY₂     = o.Repl.invXX * o.S✻⋂XY₂ 
	o.invXXS✻⋂XX      = o.Repl.invXX * o.S✻⋂XX  
	o.invXXS✻⋂XDGPZ   = o.Repl.invXX * o.S✻⋂XDGPZ 
	o.invXXS✻⋂Xy₁     = o.Repl.invXX * o.S✻⋂Xy₁ 
	o.invXXS✻⋂XDGPZR₁ = o.Repl.invXX * o.S✻⋂X_DGPZR₁

	_S✻ZperpY₂      = @panelsum(o, o.Repl.S✻⋂ZperpY₂ , o.info✻_✻⋂)  # moments of variables _before_ FWL processing
	_S✻Zperpy₁      = dropdims(@panelsum(o, reshape(o.Repl.S✻⋂Zperpy₁,Val(3)), o.info✻_✻⋂); dims=3)
	_S✻ZperpDGPZpar = @panelsum(o, o.DGP.S✻⋂ZperpZpar, o.info✻_✻⋂)
	_S✻ZperpDGPZR₁  = @panelsum(o, o.DGP.S✻⋂ZperpZR₁, o.info✻_✻⋂)

	S✻ZperpZperp    = @panelsum(o, o.Repl.S✻⋂ZperpZperp, o.info✻_✻⋂)
	o.S✻XY₂         = @panelsum(o, o.S✻⋂XY₂   , o.info✻_✻⋂)
	o.S✻XX          = @panelsum(o, o.S✻⋂XX    , o.info✻_✻⋂)
	o.S✻XDGPZ       = @panelsum(o, o.S✻⋂XDGPZ, o.info✻_✻⋂)
	o.S✻Xy₁         = dropdims(@panelsum(o, reshape(o.S✻⋂Xy₁,Val(3)), o.info✻_✻⋂); dims=3)
	o.S✻XZR₁        = @panelsum(o, o.S✻⋂X_DGPZR₁, o.info✻_✻⋂)
	o.S✻ZperpX      = @panelsum(o, S✻⋂ZperpX, o.info✻_✻⋂)
	o.S✻ZperpY₂     = _S✻ZperpY₂ - S✻ZperpZperp * o.Repl.invZperpZperpZperpY₂
	o.S✻ZperpDGPZ   = _S✻ZperpDGPZpar - S✻ZperpZperp * o.DGP.invZperpZperpZperpZpar
	o.S✻ZperpDGPZR₁ = @panelsum(o, o.DGP.S✻⋂ZperpZR₁ , o.info✻_✻⋂) - S✻ZperpZperp * o.DGP.invZperpZperpZperpZR₁
	o.S✻Zperpy₁     = _S✻Zperpy₁ - S✻ZperpZperp * o.Repl.invZperpZperpZperpy₁

	if o.NFE>0 && (o.LIML || !isone(o.κ) || o.bootstrapt)
		  CT✻⋂FEX  = [crosstabFE(o, o.Repl.X₁, o.info✻⋂) crosstabFE(o, o.Repl.X₂, o.info✻⋂)]
		o.CT✻FEX   = @panelsum(o, CT✻⋂FEX, o.info✻_✻⋂)
		o.CT✻FEY₂  = crosstabFE(o, o.DGP.Y₂, o.info✻)
		o.CT✻FEZ   = crosstabFE(o, o.DGP.Z, o.info✻)
		o.CT✻FEy₁  = crosstabFE(o, o.DGP.y₁, o.info✻)
		o.CT✻FEZR₁ = crosstabFE(o, o.DGP.ZR₁, o.info✻)
	end

	if ((o.robust && o.bootstrapt) || o.LIML || !o.robust || !isone(o.κ))
		S✻⋂ReplZX = (o.Repl.S✻⋂XZpar - o.Repl.S✻⋂XZperp * o.Repl.invZperpZperpZperpZpar - o.Repl.invZperpZperpZperpX' * (o.Repl.S✻⋂ZperpZpar - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpZpar))'
	end

	if o.bootstrapt & o.robust
		o.info⋂_✻⋂ = panelsetup(o.ID✻⋂, o.subcluster+1:o.NClustVar)
		o.S⋂YX       = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
		o.S⋂PXYZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
		o.FillingT₀ = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)  # fixed component of groupwise term in sandwich filling

		if o.NFE>0
			o.CT_FE⋂PY = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			_CT_FE⋂PY  = o.invFEwt .* @panelsum(o, CT✻⋂FEX, o.info⋂_✻⋂) * o.Repl.invXXXZ
		end

		S⋂ReplZ_X = @panelsum(o, S✻⋂ReplZX, o.info⋂_✻⋂)
		S⋂ZperpX  = @panelsum(o, S✻⋂ZperpX, o.info⋂_✻⋂)
		S⋂Xy₁  = dropdims(@panelsum(o, reshape(o.S✻⋂Xy₁,Val(3)), o.info⋂_✻⋂); dims=3)
		@inbounds for i ∈ 1:o.Repl.kZ
			o.NFE>0 &&
				(o.CT_FE⋂PY[i+1] = view(_CT_FE⋂PY,:,:,i))
			tmp = (S⋂ReplZ_X * view(o.Repl.invXXXZ,:,i))'
			for j ∈ 1:o.Repl.kZ
				o.FillingT₀[i+1,j+1] = reshape(view(tmp,:,j),:,1)
			end

			o.S⋂PXYZperp[i+1] = (S⋂ZperpX * view(o.Repl.invXXXZ,:,i))'  # S⋂(P_(MZperpX) * Z .* Zperp)
			iszero(o.granular) &&
				(o.S⋂YX[i+1] = view(S⋂ReplZ_X,i,:,:))  # S⋂(M_Zperp[Z or y₁] .* P_(MZperpX)])
		end

		tmp = (o.Repl.invXXXZ' * S⋂Xy₁)'
		if o.Repl.restricted
			S⋂ReplXZR₁r₁  = dropdims(@panelsum(o, reshape((o.Repl.S✻⋂XZR₁ - 
			                              o.Repl.S✻⋂XZperp * o.Repl.invZperpZperpZperpZR₁ - 
																		o.Repl.invZperpZperpZperpX' * (o.Repl.S✻⋂ZperpZR₁ - o.Repl.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpZR₁)) * o.r₁, Val(3)), o.info⋂_✻⋂); dims=3)
			tmp .-= S⋂ReplXZR₁r₁'o.Repl.invXXXZ
		end

		for i ∈ 1:o.Repl.kZ
			o.FillingT₀[i+1,1] = reshape(view(tmp,:,i),:,1)
		end

		if o.granular
			o.XinvXX = X₁₂B(o, o.Repl.X₁, o.Repl.X₂, o.invXX)
			o.PXZ    = X₁₂B(o, o.Repl.X₁, o.Repl.X₂, o.Repl.invXXXZ)
		else
			o.S⋂YX[1] = dropdims(@panelsum(o, reshape(o.S✻⋂Xy₁,Val(3)), o.info⋂_✻⋂); dims=3)'  # S⋂(M_Zperp*y₁ .* P_(MZperpX)])
			o.Repl.restricted &&
				(o.S⋂YX[1] .-= S⋂ReplXZR₁r₁')

			inds = o.subcluster>0 ?
							[CartesianIndex(i,j) for (j,v) ∈ enumerate(o.info⋂_✻⋂) for i ∈ v] :  # crosstab *,∩ is tall
							o.NClustVar == o.nbootclustvar ?
									[CartesianIndex(i,i) for i ∈ 1:o.N✻⋂] :  # crosstab *,∩ is square
									[CartesianIndex(j,i) for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab *,∩ is wide
			inds = [CartesianIndex(k,I) for I ∈ inds for k ∈ 1:o.Repl.kX]
			o.crosstab⋂✻ind = LinearIndices(FakeArray(Tuple(max(inds...))...))[inds]

			o.S⋂XZperpinvZperZperp = S⋂ZperpX' * o.Repl.invZperpZperp

			o.NFE>0 && (o.CT⋂FEX = o.invFEwt .* @panelsum(o, CT✻⋂FEX, o.info⋂_✻⋂))
		end
	end

	if o.LIML || !o.robust || !isone(o.κ)  # cluster-wise moments after FWL
		o.S✻Y₂Y₂     = o.Repl.S✻Y₂Y₂    - _S✻ZperpY₂'   * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpY₂'   * o.S✻ZperpY₂
		o.S✻DGPZDGPZ = o.DGP.S✻ZparZpar - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpZpar' * o.S✻ZperpDGPZ
		o.S✻DGPZY₂   = o.DGP.S✻ZparY₂   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpZpar' * o.S✻ZperpY₂
		o.S✻DGPZy₁   = o.DGP.S✻Zpary₁   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpZpar' * o.S✻Zperpy₁   
		o.S✻Y₂y₁     = o.Repl.S✻Y₂y₁    - _S✻ZperpY₂'   * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpY₂'   * o.S✻Zperpy₁
		o.S✻y₁y₁     = o.Repl.S✻y₁y₁    - _S✻Zperpy₁'   * o.DGP.invZperpZperpZperpy₁   - o.S✻Zperpy₁' * o.DGP.invZperpZperpZperpy₁
		o.DGP.restricted && (o.S✻DGPZR₁y₁ = o.DGP.S✻ZR₁y₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpy₁ - o.DGP.invZperpZperpZperpZR₁' * o.S✻Zperpy₁)

		if o.Repl.restricted
			_S✻ZperpReplZR₁ = @panelsum(o, o.Repl.S✻⋂ZperpZR₁, o.info✻_✻⋂)
			_S✻⋂XReplZR₁    = @panelsum(o, o.Repl.S✻⋂XZR₁    , o.info✻_✻⋂)
			
			o.r₁S✻ReplZR₁Y₂     = o.r₁' * (o.Repl.S✻ZR₁Y₂ - _S✻ZperpReplZR₁' * o.Repl.invZperpZperpZperpY₂ - o.Repl.invZperpZperpZperpZR₁' * o.S✻ZperpY₂)  
			o.r₁S✻ReplZR₁X      = o.r₁' * (_S✻⋂XReplZR₁'  - _S✻ZperpReplZR₁' * o.Repl.invZperpZperpZperpX  - o.Repl.invZperpZperpZperpZR₁' * o.S✻ZperpX )
			o.r₁S✻ReplZR₁DGPZ   = o.r₁' * panelcross11(o, o.Repl.ZR₁, o.DGP.Z, o.info✻)   
			o.r₁S✻ReplZR₁y₁     = o.r₁' * (o.Repl.S✻ZR₁y₁ - _S✻ZperpReplZR₁' * o.Repl.invZperpZperpZperpy₁ - o.Repl.invZperpZperpZperpZR₁' * o.S✻Zperpy₁)
			o.DGP.restricted &&
				(o.r₁S✻ReplZR₁DGPZR₁ = o.r₁' * panelcross11(o, o.Repl.ZR₁, o.DGP.ZR₁, o.info✻))
		end

		_S✻ZperpReplZpar = @panelsum(o, o.Repl.S✻⋂ZperpZpar, o.info✻_✻⋂)
		_S✻ReplXZ        = @panelsum(o, o.Repl.S✻⋂XZpar    , o.info✻_✻⋂)

		o.S✻ReplZY₂      = o.Repl.S✻ZparY₂ - _S✻ZperpReplZpar' * o.Repl.invZperpZperpZperpY₂ - o.Repl.invZperpZperpZperpZpar' * o.S✻ZperpY₂
		o.S✻ReplZX       = _S✻ReplXZ'      - _S✻ZperpReplZpar' * o.Repl.invZperpZperpZperpX  - o.Repl.invZperpZperpZperpZpar' * o.S✻ZperpX
		o.S✻ReplZDGPZ    = panelcross11(o, o.Repl.Z, o.DGP.Z, o.info✻)   
		o.S✻ReplZy₁      = o.Repl.S✻Zpary₁ - _S✻ZperpReplZpar' * o.Repl.invZperpZperpZperpy₁ - o.Repl.invZperpZperpZperpZpar' * o.S✻Zperpy₁
		
		if o.DGP.restricted
			_S✻⋂XDGPZR₁ = @panelsum(o, o.DGP.S✻⋂XZR₁, o.info✻_✻⋂)

			o.S✻ReplZDGPZR₁  = panelcross11(o, o.Repl.Z, o.DGP.ZR₁, o.info✻)   
			o.S✻DGPZR₁Y₂     = o.DGP.S✻ZR₁Y₂  - _S✻ZperpDGPZR₁' * o.Repl.invZperpZperpZperpY₂  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpY₂
			o.S✻DGPZR₁DGPZR₁ = o.DGP.S✻ZR₁ZR₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZR₁  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZR₁
			o.S✻DGPZR₁DGPZ   = o.DGP.S✻ZR₁Z   - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZ
			o.S✻DGPZR₁X      = _S✻⋂XDGPZR₁'   - _S✻ZperpDGPZR₁' * o.Repl.invZperpZperpZperpX   - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpX
		end
	end

	o.invXXS✻XY₂   = @panelsum(o, o.invXXS✻⋂XY₂  , o.info✻_✻⋂)
	o.invXXS✻XX    = @panelsum(o, o.invXXS✻⋂XX   , o.info✻_✻⋂)
	o.invXXS✻XDGPZ = @panelsum(o, o.invXXS✻⋂XDGPZ, o.info✻_✻⋂)
	o.invXXS✻Xy₁   = dropdims(@panelsum(o, reshape(o.invXXS✻⋂Xy₁,Val(3)), o.info✻_✻⋂); dims=3)
	o.invZperpZperpS✻ZperpY₂   = o.Repl.invZperpZperp * o.S✻ZperpY₂ 
	o.invZperpZperpS✻ZperpX    = o.Repl.invZperpZperp * o.S✻ZperpX  
	o.invZperpZperpS✻Zperpy₁   = o.Repl.invZperpZperp * o.S✻Zperpy₁ 
	o.invZperpZperpS✻ZperpDGPZ = o.Repl.invZperpZperp * o.S✻ZperpDGPZ

	if o.DGP.restricted
		o.invXXS✻XDGPZR₁ = @panelsum(o, o.invXXS✻⋂XDGPZR₁, o.info✻_✻⋂)
		o.invZperpZperpS✻ZperpDGPZR₁ = o.Repl.invZperpZperp * o.S✻ZperpDGPZR₁
	end
end

function PrepWRE!(o::StrBootTest{T}) where T
	r₁ = o.null ? [o.r₁ ; o.r] : o.r₁
  EstimateIV!(o.DGP, o, r₁)
  MakeResidualsIV!(o.DGP, o)
  o.robust && o.bootstrapt && o.granular && (Ü₂par = view(o.DGP.Ü₂ * o.Repl.RparY,:,:))

	S✻⋂XU₂ = o.S✻⋂XY₂ - o.S✻⋂XX * o.DGP.Π̂  # XXX preallocate these
	S✻⋂XU₂RparY = S✻⋂XU₂ * o.Repl.RparY
	S✻XU₂ = o.S✻XY₂ - o.S✻XX * o.DGP.Π̂
	S✻XU₂RparY = S✻XU₂ * o.Repl.RparY
	
	S✻ZperpU₂ = o.S✻ZperpY₂ - o.S✻ZperpX * o.DGP.Π̂
	S✻ZperpU₂RparY = S✻ZperpU₂ * o.Repl.RparY
	invZperpZperpS✻ZperpU₂ = o.invZperpZperpS✻ZperpY₂ - o.invZperpZperpS✻ZperpX * o.DGP.Π̂
	invZperpZperpS✻ZperpU₂RparY = invZperpZperpS✻ZperpU₂ * o.Repl.RparY

	invXXS✻XU₂ = o.Repl.invXX * S✻XU₂
	invXXS✻XU₂RparY = invXXS✻XU₂ * o.Repl.RparY
	if iszero(o.granular)
		invXXS✻⋂XU₂ = o.Repl.invXX * S✻⋂XU₂
		invXXS✻⋂XU₂RparY = invXXS✻⋂XU₂ * o.Repl.RparY
	end

	if o.LIML || !o.robust || !isone(o.κ)
		S✻U₂y₁ = o.S✻Y₂y₁ - o.DGP.Π̂' * o.S✻Xy₁
		S✻ZU₂ = o.S✻ReplZY₂ - o.S✻ReplZX * o.DGP.Π̂
		S✻ZU₂RparY = S✻ZU₂ * o.Repl.RparY
		Π̂S✻XÜ₂γ̈ = o.DGP.Π̂' * S✻XU₂ * o.DGP.γ̈
		S✻Ü₂Y₂ = o.S✻Y₂Y₂ - o.DGP.Π̂' * o.S✻XY₂
		S✻Y₂Ü₂γ̈ = S✻Ü₂Y₂' * o.DGP.γ̈
		S✻Ü₂parÜ₂par = o.Repl.RparY' * (S✻Ü₂Y₂ - S✻XU₂' * o.DGP.Π̂) * o.Repl.RparY

		if o.Repl.restricted
			S✻ReplZR₁r₁U₂ = o.r₁S✻ReplZR₁Y₂ - o.r₁S✻ReplZR₁X * o.DGP.Π̂
			S✻ReplZR₁r₁U₂RparY = S✻ReplZR₁r₁U₂ * o.Repl.RparY
		end	
		if o.DGP.restricted
			r₁S✻DGPZR₁y₁ = r₁' * o.S✻DGPZR₁y₁
			S✻ReplZDGPZR₁r₁ = o.S✻ReplZDGPZR₁ * r₁
			o.Repl.restricted && (S✻ReplZR₁r₁DGPZR₁r₁ = o.r₁S✻ReplZR₁DGPZR₁ * r₁)
		end

		if o.NFE>0
			CT✻FEU = o.CT✻FEY₂ - o.CT✻FEX * o.DGP.Π̂
			CT✻FEURparY = CT✻FEU * o.Repl.RparY
		end
	end

  @inbounds for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
		if iszero(i)
			o.S✻XU[1]      .= o.S✻Xy₁      - o.S✻XDGPZ      * o.DGP.β̈ + S✻XU₂      * o.DGP.γ̈
			o.invXXS✻XU[1] .= o.invXXS✻Xy₁ - o.invXXS✻XDGPZ * o.DGP.β̈ + invXXS✻XU₂ * o.DGP.γ̈
			if o.DGP.restricted
				o.S✻XU[1]      .-= o.S✻XZR₁      * r₁
				o.invXXS✻XU[1] .-= o.invXXS✻XDGPZR₁ * r₁
			end
		else
			o.S✻XU[i+1]      .= view(S✻XU₂RparY,:,:,i)
			o.invXXS✻XU[i+1] .= view(invXXS✻XU₂RparY,:,:,i)
		end

		if o.LIML || !isone(o.κ) || o.bootstrapt
			if iszero(i)
				o.S✻ZperpU[1]              .= o.S✻Zperpy₁              - o.S✻ZperpDGPZ              * o.DGP.β̈ + S✻ZperpU₂              * o.DGP.γ̈
				o.invZperpZperpS✻ZperpU[1] .= o.invZperpZperpS✻Zperpy₁ - o.invZperpZperpS✻ZperpDGPZ * o.DGP.β̈ + invZperpZperpS✻ZperpU₂ * o.DGP.γ̈
				if o.DGP.restricted
					o.S✻ZperpU[1]              .-= o.S✻ZperpDGPZR₁              * r₁
					o.invZperpZperpS✻ZperpU[1] .-= o.invZperpZperpS✻ZperpDGPZR₁ * r₁
				end
			else
				o.S✻ZperpU[i+1]              .= view(S✻ZperpU₂RparY,:,:,i)
				o.invZperpZperpS✻ZperpU[i+1] .= view(invZperpZperpS✻ZperpU₂RparY,:,:,i)
			end
	
			if o.NFE>0
				if iszero(i)
					o.CTFEU[1] = o.CT✻FEy₁ - o.CT✻FEZ * o.DGP.β̈ + CT✻FEU * o.DGP.γ̈
					o.DGP.restricted &&
						(o.CTFEU[1] .-= o.CT✻FEZR₁ * r₁)
				else
					o.CTFEU[i+1] = @view CT✻FEURparY[:,:,i]
				end
			end
		end

		if o.LIML || !isone(o.κ) || !o.robust
			if iszero(i)  # panelsum2(o, o.Repl.y₁par, o.Repl.Z, uwt, o.info✻)
				o.S✻YU[1,1] .= o.S✻y₁y₁ - o.S✻DGPZy₁'o.DGP.β̈ + S✻U₂y₁'o.DGP.γ̈
				o.DGP.restricted &&
					(o.S✻YU[1,1] .-= view(r₁S✻DGPZR₁y₁, 1,:))

				if o.Repl.restricted
					o.S✻YU[1,1] .-= view(o.r₁S✻ReplZR₁y₁ - o.r₁S✻ReplZR₁DGPZ * o.DGP.β̈ + S✻ReplZR₁r₁U₂ * o.DGP.γ̈, 1,:,1)
					o.DGP.restricted &&
						(o.S✻YU[1,1] .+= view(S✻ReplZR₁r₁DGPZR₁r₁, 1,:,1))
				end

				for j ∈ 1:o.Repl.kZ
					o.S✻YU[j+1,1] .= view(o.S✻ReplZy₁ - o.S✻ReplZDGPZ * o.DGP.β̈ + S✻ZU₂ * o.DGP.γ̈, j,:,1)
					o.DGP.restricted &&
						(o.S✻YU[j+1,1] .-= view(S✻ReplZDGPZR₁r₁, j,:,1))
				end
      else
				o.S✻YU[1,i+1] .= @view S✻U₂y₁[i,:,1]
				o.Repl.restricted &&
					(o.S✻YU[1,i+1] .-= @view S✻ReplZR₁r₁U₂RparY[1,:,i])
				for j ∈ 1:o.Repl.kZ
					o.S✻YU[j+1,i+1] .= @view S✻ZU₂RparY[j,:,i]
				end
			end
			XXX = o.S✻Y₂y₁ - o.S✻DGPZY₂' * o.DGP.β̈ - o.DGP.Π̂' * (o.S✻Xy₁  - o.S✻XDGPZ * o.DGP.β̈)
			if iszero(i)
		 		o.S✻UU[1,i+1] .= o.S✻y₁y₁ - (2 * o.S✻DGPZy₁ + o.S✻DGPZDGPZ * o.DGP.β̈)'o.DGP.β̈ + (2 * XXX - Π̂S✻XÜ₂γ̈ + S✻Y₂Ü₂γ̈)'o.DGP.γ̈
				o.DGP.restricted
					(o.S✻UU[1,i+1] .+= dropdims(-2 * r₁S✻DGPZR₁y₁ + r₁' * (o.S✻DGPZR₁DGPZR₁ * r₁) + 2 * r₁' * (o.S✻DGPZR₁DGPZ * o.DGP.β̈ + (o.S✻DGPZR₁X * o.DGP.Π̂ - o.S✻DGPZR₁Y₂) * o.DGP.γ̈); dims=1))  # XXX dropdims vs view...
			else
				o.S✻UU[1,i+1] .= view(o.Repl.RparY' * (XXX + S✻Y₂Ü₂γ̈ - Π̂S✻XÜ₂γ̈), i,:)

				if o.DGP.restricted
					o.S✻UU[1,i+1] .-= view(r₁' * (o.S✻DGPZR₁Y₂ - o.S✻DGPZR₁X * o.DGP.Π̂) * o.Repl.RparY, 1,:,i)
				end
			end
			for j ∈ 1:i
				o.S✻UU[j+1,i+1] .= view(S✻Ü₂parÜ₂par, j,:,i)
			end
		end

		if o.robust && o.bootstrapt
			if iszero(o.granular)  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				if iszero(i)
					_S✻⋂UXinvXX = o.invXXS✻⋂Xy₁ - o.invXXS✻⋂XDGPZ * o.DGP.β̈ + invXXS✻⋂XU₂ * o.DGP.γ̈
					o.DGP.restricted &&
						(_S✻⋂UXinvXX .-= o.invXXS✻⋂XDGPZR₁ * r₁)
					if !(o.NClustVar == o.nbootclustvar && iszero(o.subcluster))
						for g ∈ 1:o.N✻
							o.S✻⋂u₁XinvXX[i+1,g] .= view(_S✻⋂UXinvXX,:,o._ID✻⋂[o.info✻[g]],1)'  # panelsum(o, o.XinvXX, u, o.infoCT⋂✻[g])
						end
					end
				else
					if !(o.NClustVar == o.nbootclustvar && iszero(o.subcluster))
						for g ∈ 1:o.N✻
							o.S✻⋂u₁XinvXX[i+1,g] .= view(invXXS✻⋂XU₂RparY,:,o._ID✻⋂[o.info✻[g]],1,i)'  # panelsum(o, o.XinvXX, u, o.infoCT⋂✻[g])
						end
					end
				end
			end

			if o.granular
				i>0 && (o.S✻UPX[i+1] .= o.XinvXX * o.S✻XU[i+1])

				o.S✻UMZperp[i+1] .= o.Repl.Zperp * o.invZperpZperpS✻ZperpU[i+1] 
				if iszero(i)  # subtract crosstab of observation by ∩-group of u
					o.S✻UMZperp[  1][o.crosstab✻ind] .-= o.DGP.u⃛₁
				else
					o.S✻UMZperp[i+1][o.crosstab✻ind] .-= view(Ü₂par,:,i)
				end

				o.NFE>0 &&
					(o.S✻UMZperp[i+1] .-= view(o.invFEwt .* o.CTFEU[i+1], o._FEID, :))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
			end
		end
  end

	if o.robust && o.bootstrapt && iszero(o.granular)
		@inbounds for ind2 ∈ 0:o.Repl.kZ
			if o.Repl.Yendog[ind2+1]
				S✻UMZperpX = o.S⋂XZperpinvZperZperp * o.S✻ZperpU[ind2+1]  # S_*  diag⁡(U ̈_(∥j) ) Z_⊥ (Z_⊥^' Z_⊥ )^(-1) Z_(⊥g)^' X_(∥g)
				if iszero(ind2)  # - S_*  diag⁡(U ̈_(∥j) ) I_g^' X_(∥g)
					S✻UMZperpX[o.crosstab⋂✻ind] .-= vec(o.S✻⋂Xy₁ - o.S✻⋂XDGPZ * o.DGP.β̈ + S✻⋂XU₂ * o.DGP.γ̈)
					o.DGP.restricted &&
						(S✻UMZperpX[o.crosstab⋂✻ind] .+= vec(o.S✻⋂X_DGPZR₁ * r₁))
				else
					S✻UMZperpX[o.crosstab⋂✻ind] .-= vec(view(S✻⋂XU₂RparY,:,:,ind2))
				end
				o.NFE>0 &&
					(S✻UMZperpX .-= o.CT⋂FEX'o.CTFEU[ind2+1])  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
 				S✻UMZperpX = o.Repl.invXX * S✻UMZperpX
				for ind1 ∈ 1:o.Repl.kZ
					o.Repl.Yendog[ind1+1] &&
						(o.S✻UPX_S✻UMZperp[ind1+1,ind2+1] .= o.S✻XU[ind1+1]'S✻UMZperpX)
				end
			end
		end
	end
	nothing
end

# stuff done once per exucution--not depending on r
function InitWRE!(o::StrBootTest{T}) where T
	o.LIML && o.Repl.kZ==1 && o.Nw==1 && (o.As = o.β̂s = zeros(1, o.B+1))
	o.S✻UZperpinvZperpZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻UZperp              = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻UY                  = Matrix{Vector{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
	o.S✻UXinvXX             = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻UX                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻UU                  = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
	o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kX, ncols(o.v))] :
												[Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[end]))]
	o.T1R = deepcopy(o.T1L)

	S✻⋂XZperp = panelcross21(o, o.Repl._X₁, o.X₂, o.DGP.Zperp, o.info✻⋂)
	S✻⋂ZperpZperp = panelcross11(o, o.Repl.Zperp, o.DGP.Zperp, o.info✻⋂)
	S✻⋂XY₂ = panelcross21(o, o.Repl._X₁, o.X₂, o.Y₂, o.info✻⋂)
	S✻⋂ZperpY₂ = panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻⋂)
	S✻⋂XX = panelcross22(o, o.Repl._X₁, o.X₂, o.DGP._X₁, o.X₂, o.info✻⋂)  # XXX X₂-X₂ part same for DGP & Repl
	S✻⋂ZperpX = panelcross12(o, o.Repl.Zperp, o.DGP._X₁, o.X₂, o.info✻⋂)
	S✻⋂Zperpy₁ = panelcross11(o, o.Repl.Zperp, o.y₁, o.info✻⋂)
	S✻⋂ZperpZR₁ = panelcross11(o, o.Repl.Zperp, o.DGP._ZR₁, o.info✻⋂)
	S✻⋂ZperpZ = panelcross11(o, o.Repl.Zperp, o.DGP._Z, o.info✻⋂)
	S✻⋂Xy₁ = panelcross21(o, o.Repl._X₁, o.X₂, o.y₁, o.info✻⋂)
	S✻⋂XZ = panelcross21(o, o.Repl._X₁, o.X₂, o.DGP._Z, o.info✻⋂)
	S✻⋂XZR₁ = panelcross21(o, o.Repl._X₁, o.X₂, o.DGP._ZR₁, o.info✻⋂)

	if o.robust && o.bootstrapt && !o.granular
		o.S✻⋂u₁XinvXX0    = o.Repl.invXX * (S✻⋂Xy₁  - S✻⋂XZperp * o.DGP.invZperpZperpZperpy₁  - o.Repl.invZperpZperpZperpX' * (S✻⋂Zperpy₁  - S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpy₁))
		o.∂S✻⋂u₁XinvXX∂β̂  = o.Repl.invXX * (S✻⋂XZ   - S✻⋂XZperp * o.DGP.invZperpZperpZperpZ   - o.Repl.invZperpZperpZperpX' * (S✻⋂ZperpZ   - S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZ))
		o.∂S✻⋂u₁XinvXX∂γ̈  = o.Repl.invXX * (S✻⋂XY₂  - S✻⋂XZperp * o.DGP.invZperpZperpZperpY₂  - o.Repl.invZperpZperpZperpX' * (S✻⋂ZperpY₂  - S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpY₂))
		o.∂S✻⋂u₁XinvXX∂γ̈Π̂ = o.Repl.invXX * (S✻⋂XX   - S✻⋂XZperp * o.DGP.invZperpZperpZperpX   - o.Repl.invZperpZperpZperpX' * (S✻⋂ZperpX   - S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpX))
		o.∂S✻⋂u₁XinvXX∂r  = o.Repl.invXX * (S✻⋂XZR₁ - S✻⋂XZperp * o.DGP.invZperpZperpZperpZR₁ - o.Repl.invZperpZperpZperpX' * (S✻⋂ZperpZR₁ - S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZR₁))
		o.S✻⋂U₂XinvXX0 = o.∂S✻⋂u₁XinvXX∂γ̈
		o.∂S✻⋂U₂XinvXX0∂Π̂ = o.∂S✻⋂u₁XinvXX∂γ̈Π̂
	end

	info✻_✻⋂, _ = panelsetupID(o.ID✻⋂, 1:o.nbootclustvar)

	S✻XZperp = @panelsum(o, S✻⋂XZperp, info✻_✻⋂)
	S✻ZperpZperp = @panelsum(o, S✻⋂ZperpZperp, info✻_✻⋂)
	S✻XY₂ = @panelsum(o, S✻⋂XY₂, info✻_✻⋂)
	S✻ZperpY₂ = @panelsum(o, S✻⋂ZperpY₂, info✻_✻⋂)
	S✻XX = @panelsum(o, S✻⋂XX, info✻_✻⋂)
	S✻ZperpX = @panelsum(o, S✻⋂ZperpX, info✻_✻⋂)
	S✻Zperpy₁ = @panelsum(o, S✻⋂Zperpy₁, info✻_✻⋂)
	S✻ZperpZR₁ = @panelsum(o, S✻⋂ZperpZR₁, info✻_✻⋂)
	S✻ZperpZ = @panelsum(o, S✻⋂ZperpZ, info✻_✻⋂)
	S✻Xy₁ = @panelsum(o, S✻⋂Xy₁, info✻_✻⋂)
	S✻XZ = @panelsum(o, S✻⋂XZ, info✻_✻⋂)
	S✻XZR₁ = @panelsum(o, S✻⋂XZR₁, info✻_✻⋂)

	o.S✻Zperp_DGPy₁ = S✻Zperpy₁  - S✻ZperpZperp * o.DGP.invZperpZperpZperpy₁
	o.S✻Zperp_DGPZ  = S✻ZperpZ   - S✻ZperpZperp * o.DGP.invZperpZperpZperpZ
	o.S✻Zperp_DGPY₂ = S✻ZperpY₂  - S✻ZperpZperp * o.DGP.invZperpZperpZperpY₂
	o.S✻Zperp_DGPX = S✻ZperpX   - S✻ZperpZperp * o.DGP.invZperpZperpZperpX
	o.S✻Zperp_DGPZR₁ = S✻ZperpZR₁ - S✻ZperpZperp * o.DGP.invZperpZperpZperpZR₁

	o.S✻u₁X0    = S✻Xy₁  - S✻XZperp * o.DGP.invZperpZperpZperpy₁ - o.Repl.invZperpZperpZperpX' * o.S✻Zperp_DGPy₁
	o.∂S✻u₁X∂β̂  = S✻XZ   - S✻XZperp * o.DGP.invZperpZperpZperpZ  - o.Repl.invZperpZperpZperpX' * o.S✻Zperp_DGPZ
	o.∂S✻u₁X∂γ̈  = S✻XY₂  - S✻XZperp * o.DGP.invZperpZperpZperpY₂ - o.Repl.invZperpZperpZperpX' * o.S✻Zperp_DGPY₂
	o.∂S✻u₁X∂γ̈Π̂ = S✻XX   - S✻XZperp * o.DGP.invZperpZperpZperpX  - o.Repl.invZperpZperpZperpX' * o.S✻Zperp_DGPX
	o.∂S✻u₁X∂r  = S✻XZR₁ - S✻XZperp * o.DGP.invZperpZperp * o.DGP.Zperp_ZR₁ - o.Repl.invZperpZperpZperpX' * o.S✻Zperp_DGPZR₁
	o.S✻U₂X0 = o.∂S✻u₁X∂γ̈ * o.Repl.RparY
	o.∂S✻U₂X∂Π̂ = o.∂S✻u₁X∂γ̈Π̂

	if o.LIML || !o.robust || !isone(o.κ)
		# S✻y₁y₁ = panelcross11(o, o.y₁, o.y₁, o.info✻)
		# S✻y₁Y₂ = panelcross11(o, o.y₁, o.Y₂, o.info✻)
		# S✻y₁Zperp = panelcross11(o, o.y₁, o.DGP.Zperp, o.info✻)
		# S✻ZZ = panelcross11(o, o.Repl._Z, o.DGP._Z, o.info✻)
		# S✻ZZperp = panelcross11(o, o.Repl._Z, o.DGP.Zperp, o.info✻)
		# S✻y₁Z = panelcross11(o, o.y₁, o.DGP._Z, o.info✻)
		# S✻ZY₂ = panelcross11(o, o.Repl._Z, o.Y₂, o.info✻)
		# S✻ZX = panelcross12(o, o.Repl._Z, o.DGP._X₁, o.X₂, o.info✻)
		# S✻ZZR₁ = panelcross11(o, o.Repl._Z, o.DGP._ZR₁, o.info✻)
		# S✻y₁X = panelcross12(o, o.y₁, o.DGP._X₁, o.X₂, o.info✻)
		# S✻y₁ZR₁ = panelcross11(o, o.y₁, o.DGP._ZR₁, o.info✻)

		o.Reply₁_DGPy₁  = panelcross11(o, o.Repl.y₁, o.DGP.y₁, o.info✻)
		o.Reply₁_DGPZR₁ = panelcross11(o, o.Repl.y₁, o.DGP.ZR₁, o.info✻)
		o.Reply₁_DGPZ   = panelcross11(o, o.Repl.y₁, o.DGP.Z, o.info✻)
		o.Reply₁_DGPY₂  = panelcross11(o, o.Repl.y₁, o.DGP.Y₂, o.info✻)
		o.Reply₁_DGPX   = panelcross11(o, o.Repl.y₁, [o.DGP.X₁ o.DGP.X₂], o.info✻)

		o.ReplZ_DGPy₁  = panelcross11(o, o.Repl.Z, o.DGP.y₁, o.info✻)
		o.ReplZ_DGPZR₁ = panelcross11(o, o.Repl.Z, o.DGP.ZR₁, o.info✻)
		o.ReplZ_DGPZ   = panelcross11(o, o.Repl.Z, o.DGP.Z, o.info✻)
		o.ReplZ_DGPY₂  = panelcross11(o, o.Repl.Z, o.DGP.Y₂, o.info✻)
		o.ReplZ_DGPX   = panelcross11(o, o.Repl.Z, [o.DGP.X₁ o.DGP.X₂], o.info✻)

		if o.REst
			ZR₁r₁ = o.Repl.ZR₁ * o.r₁
			o.Reply₁_DGPy₁ .-= panelcross11(o, ZR₁r₁, o.DGP.y₁, o.info✻)
			o.Reply₁_DGPZR₁ .-= panelcross11(o, ZR₁r₁, o.DGP.ZR₁, o.info✻)
			o.Reply₁_DGPZ .-= panelcross11(o, ZR₁r₁, o.DGP.Z, o.info✻)
			o.Reply₁_DGPY₂ .-= panelcross11(o, ZR₁r₁, o.DGP.Y₂, o.info✻)
			o.Reply₁_DGPX .-= panelcross11(o, ZR₁r₁, [o.DGP.X₁ o.DGP.X₂], o.info✻)
		end

		# o.∂S✻y₁U₂∂γ̈ = S✻y₁Y₂ - o.Repl.invZperpZperpZperpy₁' * S✻ZperpY₂ - (S✻y₁Zperp - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂
		# o.S✻y₁u₁0 = -o.∂S✻y₁U₂∂γ̈ * o.Repl.RparY
		# o.∂S✻y₁u₁∂β̂ = S✻ZZ - S✻ZZperp * o.DGP.invZperpZperpZperpZ - (S✻y₁Z - S✻y₁Zperp * o.DGP.invZperpZperpZperpZ)
		# o.∂S✻y₁u₁∂γ̈ = S✻ZY₂ - o.Repl.invZperpZperpZperpy₁' * S✻ZperpY₂ - (S✻ZZperp - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂
		# o.∂S✻y₁u₁∂γ̈Π̂ = S✻ZX   - S✻ZZperp * o.DGP.invZperpZperpZperpX - (S✻y₁X - S✻y₁Zperp * o.DGP.invZperpZperpZperpX)
		# o.∂S✻y₁u₁∂r = S✻ZZR₁  - S✻ZZperp * o.DGP.invZperpZperpZperpZR₁ - (S✻y₁ZR₁ - S✻y₁Zperp * o.DGP.invZperpZperpZperpZR₁)
		# o.S✻y₁U₂0 = S✻y₁y₁ - o.Repl.invZperpZperpZperpy₁' * S✻Zperpy₁ - o.∂S✻y₁U₂∂γ̈
		# o.∂S✻y₁U₂∂β̂ = S✻y₁Z - S✻y₁Zperp * o.DGP.invZperpZperpZperpZ - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZ
		# o.∂S✻y₁U₂∂γ̈Π̂ = S✻y₁X - o.Repl.invZperpZperpZperpy₁' * S✻ZperpX - (S✻y₁Zperp - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpX
		# o.∂S✻y₁U₂∂Π̂ =  o.∂S✻y₁U₂∂γ̈Π̂
		# o.∂S✻y₁U₂∂r = S✻y₁ZR₁ - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZR₁ - (S✻y₁Zperp - o.Repl.invZperpZperpZperpy₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpZR₁

		# if o.REst
		# 	o.∂S✻y₁u₁∂β̂ .-= o.r₁' * (panelcross11(o, o.Repl._ZR₁, o.DGP._Z, o.info✻) - panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) * o.DGP.invZperpZperpZperpZ -
		# 											o.Repl.invZperpZperpZperpZR₁' * (S✻ZperpZ - 
		# 																												S✻ZperpZperp * o.DGP.invZperpZperpZperpZ) -
		# 																												(panelcross11(o, o.Repl._ZR₁, o.DGP._Z, o.info✻) -
		# 																												panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) * o.DGP.invZperpZperpZperpZ -
		# 																												o.Repl.invZperpZperpZperpZR₁' * (S✻ZperpZ -
		# 																																												 S✻ZperpZperp * o.DGP.invZperpZperpZperpZ)))
		# 	o.∂S✻y₁u₁∂γ̈ .-= o.r₁' * ((panelcross11(o, o.Repl._ZR₁, o.Y₂, o.info✻) -
		# 										o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻)) -
		# 									(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 										o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂ -
		# 										((panelcross11(o, o.Repl._ZR₁, o.Y₂, o.info✻) -
		# 										 o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻)) - 
		# 											(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 										 o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂))
		# 	o.∂S✻y₁u₁∂γ̈Π̂ .-= o.r₁' * ((panelcross12(o, o.Repl._ZR₁, o.DGP._X₁, o.X₂, o.info✻) -
		# 										panelcross12(o, o.Repl.Zperp * o.Repl.invZperpZperpZperpZR₁, o.DGP._X₁, o.X₂, o.info✻)) -
		# 										(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 										o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpX - 
		# 										((panelcross12(o, o.Repl._ZR₁, o.DGP._X₁, o.X₂, o.info✻) -
		# 									 panelcross12(o, o.Repl.Zperp * o.Repl.invZperpZperpZperpZR₁, o.DGP._X₁, o.X₂, o.info✻)) -
		# 										 (panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 									 o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpX))
		# 	o.∂S✻y₁u₁∂r .-= o.r₁' * ((panelcross11(o, o.Repl._ZR₁, o.DGP._ZR₁, o.info✻) -
		# 									o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZR₁) -
		# 									(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 									o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpZR₁ -
		# 									((panelcross11(o, o.Repl._ZR₁, o.DGP._ZR₁, o.info✻) -
		# 									o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZR₁) -
		# 									 (panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 									o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpZR₁))
		# 	o.S✻y₁U₂0 .-= o.r₁' * (panelcross11(o, o.Repl._ZR₁, o.y₁, o.info✻) -
		# 	                        o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.y₁, o.info✻) +
		# 													(panelcross11(o, o.Repl._ZR₁, o.Y₂, o.info✻) -
		# 													o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻)) -
		# 													(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 													o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂) * o.Repl.RparY
		# 	o.∂S✻y₁U₂∂β̂ .-= o.r₁' * (panelcross11(o, o.Repl._ZR₁, o.DGP._Z, o.info✻) - panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) * o.DGP.invZperpZperpZperpZ -
		# 														o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZ)
		# 	o.∂S✻y₁U₂∂γ̈ .-= o.r₁' * ((panelcross11(o, o.Repl._ZR₁, o.Y₂, o.info✻) -
    #                       o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻)) -
    #                      (panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
    #                       o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂)
		# 	o.∂S✻y₁U₂∂γ̈Π̂ .-= o.r₁' * ((panelcross12(o, o.Repl._ZR₁, o.DGP._X₁, o.X₂, o.info✻) -
		# 										 panelcross12(o, o.Repl.Zperp * o.Repl.invZperpZperpZperpZR₁, o.DGP._X₁, o.X₂, o.info✻)) -
    #                       (panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
    #                       o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpX)
		# 	o.∂S✻y₁U₂∂Π̂ .+= o.r₁' * ((panelcross11(o, o.Repl._ZR₁, o.Y₂, o.info✻) -
		# 														o.Repl.invZperpZperpZperpZR₁' * panelcross11(o, o.Repl.Zperp, o.Y₂, o.info✻)) -
		# 														(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 														o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpY₂)
		# 	o.∂S✻y₁U₂∂r .-= o.r₁' * ((panelcross11(o, o.Repl._ZR₁, o.DGP._ZR₁, o.info✻) -
		# 														o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZR₁) -
		# 													(panelcross11(o, o.Repl._ZR₁, o.DGP.Zperp, o.info✻) -
		# 														o.Repl.invZperpZperpZperpZR₁' * S✻ZperpZperp) * o.DGP.invZperpZperpZperpZR₁)
		# end
	end
end

function PrepWRE!(o::StrBootTest{T}) where T
  r₁ = o.null ? [o.r₁ ; o.r] : o.r₁
	EstimateIV!(o.DGP, o, r₁)
  MakeResidualsIV!(o.DGP, o)
  Ü₂par = view(o.DGP.Ü₂ * o.Repl.RparY,:,:)

	_S✻U₂X = o.S✻U₂X0 - o.∂S✻U₂X∂Π̂ * o.DGP.Π̂ * o.Repl.RparY  # panelsum2(o, o.Repl.X₁, o.Repl.X₂, uwt, o.info✻)'
	(o.LIML || o.bootstrapt || !isone(o.κ)) && 
		(_S✻U₂Zperp = (o.S✻Zperp_DGPY₂ - o.S✻Zperp_DGPX * o.DGP.Π̂) * o.Repl.RparY)
	o.robust && o.bootstrapt && !o.granular &&
		(_S✻⋂U₂XinvXX = o.S✻⋂U₂XinvXX0 - o.∂S✻⋂U₂XinvXX0∂Π̂ * o.DGP.Π̂ * o.Repl.RparY)
	if o.LIML || !o.robust || !isone(o.κ)
		Reply₁_DGPU₂ = o.Reply₁_DGPY₂ - o.Reply₁_DGPX * o.DGP.Π̂
		ReplZ_DGPU₂  = o.ReplZ_DGPY₂  - o.ReplZ_DGPX  * o.DGP.Π̂
		Reply₁_DGPU₂par = Reply₁_DGPU₂ * o.Repl.RparY
		ReplZ_DGPU₂par = ReplZ_DGPU₂ * o.Repl.RparY
	end

  @inbounds for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		uwt = i>0 ? view(Ü₂par,:,i) : view(o.DGP.u⃛₁,:)::Union{Vector{T}, SubArray{T, 1}}

		# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
		if iszero(i)  # panelsum2(o, o.Repl.X₁, o.Repl.X₂, uwt, o.info✻)'
			o.S✻UX[1] = panelsum2(o, o.Repl.X₁, o.Repl.X₂, uwt, o.info✻)'  #dropdims(o.S✻u₁X0 - o.∂S✻u₁X∂β̂ * o.DGP.β̈ + (o.∂S✻u₁X∂γ̈ - o.∂S✻u₁X∂γ̈Π̂ * o.DGP.Π̂) * o.DGP.γ̈; dims=3)
			# o.DGP.restricted &&
			# 	(o.S✻UX[1] .-= dropdims(o.∂S✻u₁X∂r * r₁; dims=3))
		else
			o.S✻UX[i+1] = panelsum2(o, o.Repl.X₁, o.Repl.X₂, uwt, o.info✻)'  # view(_S✻U₂X,:,:,i)
		end

		o.S✻UXinvXX[i+1] = o.Repl.invXX * o.S✻UX[i+1]

		if o.LIML || !isone(o.κ) || o.bootstrapt
			if iszero(i)  # panelsum(o, o.Repl.Zperp, uwt, o.infoBootData)'
				o.S✻UZperp[i+1] = view(o.S✻Zperp_DGPy₁ - o.S✻Zperp_DGPZ * o.DGP.β̈ + (o.S✻Zperp_DGPY₂ - o.S✻Zperp_DGPX * o.DGP.Π̂) * o.DGP.γ̈,:,:,1)
				o.DGP.restricted &&
					(o.S✻UZperp[i+1] .-= view(o.S✻Zperp_DGPZR₁ * r₁,:,:,1))
			else
				o.S✻UZperp[i+1] = view(_S✻U₂Zperp,:,:,i)
			end
			o.S✻UZperpinvZperpZperp[i+1] = o.Repl.invZperpZperp * o.S✻UZperp[i+1]
			o.NFE>0 && (o.CTFEU[i+1] = crosstabFE(o, uwt, o.info✻))
		end

		if o.LIML || !isone(o.κ) || !o.robust
			if iszero(i)  # panelsum2(o, o.Repl.y₁par, o.Repl.Z, uwt, o.info✻)
				o.S✻UY[1,1] = dropdims(o.Reply₁_DGPy₁ - o.Reply₁_DGPZR₁ * r₁ - o.Reply₁_DGPZ * o.DGP.β̈ - Reply₁_DGPU₂ * o.DGP.γ̈; dims=(1,3))
				for j ∈ 1:o.Repl.kZ
					o.S✻UY[1,j+1] = @view (o.ReplZ_DGPy₁ - o.ReplZ_DGPZR₁ * r₁ - o.ReplZ_DGPZ * o.DGP.β̈ - ReplZ_DGPU₂ * o.DGP.γ̈)[j,:,i]
				end
      else
				o.S✻UY[i+1,1] = @view Reply₁_DGPU₂par[1,:,i]
				for j ∈ 1:o.Repl.kZ
					o.S✻UY[i+1,j+1] = @view ReplZ_DGPU₂par[j,:,i]
				end
			end
			for j ∈ 0:i
				o.S✻UU[i+1,j+1] = vec(@panelsum(o, j>0 ? view(Ü₂par,:,j) : view(o.DGP.u⃛₁,:), uwt, o.info✻))
			end
		end

		if o.robust && o.bootstrapt
			if !o.granular  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				_S✻⋂U₂XinvXX = view(o.S✻⋂U₂XinvXX0 - o.∂S✻⋂U₂XinvXX0∂Π̂ * o.DGP.Π̂ * o.Repl.RparY,:,:,i)
				if iszero(i)
					_S✻⋂UXinvXX = dropdims(o.S✻⋂u₁XinvXX0 - o.∂S✻⋂u₁XinvXX∂β̂ * o.DGP.β̈ + (o.∂S✻⋂u₁XinvXX∂γ̈  - o.∂S✻⋂u₁XinvXX∂γ̈Π̂ * o.DGP.Π̂) * o.DGP.γ̈; dims=3)
					o.DGP.restricted
						(_S✻⋂UXinvXX .-= dropdims(o.∂S✻⋂u₁XinvXX∂r * r₁; dims=3))
				else
					_S✻⋂UXinvXX = view(_S✻⋂U₂XinvXX,:,:,i)
				end
				for g ∈ 1:o.N✻
					o.S✻⋂u₁XinvXX[i+1,g] = view(_S✻⋂UXinvXX,:,o._ID✻⋂[o.info✻[g]],1)'  # panelsum(o, o.Repl.XinvXX, u, o.infoCT⋂✻[g])
				end
			end

			i>0 && (o.S✻UPX[i+1] = o.Repl.XinvXX * o.S✻UX[i+1])
			o.S✻UMZperp[i+1] = o.Repl.Zperp * o.S✻UZperpinvZperpZperp[i+1]  # S_* diag⁡(U ̈_(∥j) ) Z_⊥ (Z_⊥^' Z_⊥ )^(-1) Z_(⊥g)^'  over all g

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

function _HessianFixedkappa!(o::StrBootTest{T}, dest::AbstractMatrix, ind1::Integer, ind2::Integer, κ::Number, w::Integer) where T
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
				mul!(T1L, o.S✻XU[ind1+1], o.v)
				if iszero(ind1)
					T1L .+= o.Repl.Xy₁par
				else
					T1L .+= view(o.Repl.XZ,:,ind1)
				end
				if o.Repl.Yendog[ind2+1]
					T1R = o.T1R[isone(o.Nw) || w<o.Nw ? 1 : 2]  # use preallocated destinations
					mul!(T1R, o.invXXS✻XU[ind2+1], o.v)
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
					mul!(T1R, o.invXXS✻XU[ind2+1], o.v)
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
				T2 = o.invZperpZperpS✻ZperpU[ind1+1]'o.S✻ZperpU[ind2+1]  # quadratic term
				T2[diagind(T2)] .-= ind1 ≤ ind2 ? o.S✻UU[ind1+1, ind2+1] : o.S✻UU[ind2+1, ind1+1]  # minus diagonal crosstab
				o.NFE>0 &&
					(T2 .+= o.CTFEU[ind1+1]' * (o.invFEwt .* o.CTFEU[ind2+1]))
				if iszero(κ)
					dest .= o.Repl.YY[ind1+1,ind2+1] .+ colquadformminus!(o, (                            o.S✻YU[ind1+1,ind2+1] .+ o.S✻YU[ind2+1,ind1+1])'o.v, T2, o.v)
				else
					dest .=   κ .* dest .+ (1 - κ)   .* colquadformminus!(o, (o.Repl.YY[ind1+1,ind2+1] .+ o.S✻YU[ind1+1,ind2+1] .+ o.S✻YU[ind2+1,ind1+1])'o.v, T2, o.v)
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
function FillingLoop1!(o::StrBootTest{T}, dest::Matrix{T}, ind1::Integer, ind2::Integer, _β̈::AbstractMatrix{T}) where T
	#=Threads.@threads=# for i ∈ 1:o.clust[1].N
		PXY✻ = hcat(o.PXZ[i,ind1])
		o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX[ind1+1],i,:)'o.v)

		if iszero(ind2)
			dest[i,:]   = colsum(PXY✻ .* (o.Repl.y₁[i] .- view(o.S✻UMZperp[1],i,:))'o.v)
		elseif o.Repl.Yendog[ind2+1]
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Z[i,ind2] * _β̈ .- view(o.S✻UMZperp[ind2+1],i,:)'β̈v))
		else
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Z[i,ind2] * _β̈))
		end
	end
	nothing
end

function FillingLoop2!(o::StrBootTest{T}, dest::Matrix{T}, ind1::Integer, ind2::Integer, _β̈::AbstractMatrix{T}) where T
	#=Threads.@threads=# for i ∈ 1:o.clust[1].N
		S = o.info⋂[i]
		PXY✻ = o.Repl.Yendog[ind1+1] ? view(o.PXZ,S,ind1) .+ view(o.S✻UPX[ind1+1],S,:) * o.v :
																	reshape(o.PXZ[S,ind1], :, 1)
		if iszero(ind2)
			dest[i,:]   = colsum(PXY✻ .* (o.Repl.y₁[S] .- view(o.S✻UMZperp[1],S,:) * o.v))
		else
			dest[i,:] .-= colsum(PXY✻ .* (o.Repl.Yendog[ind2+1] ? o.Repl.Z[S,ind2] * _β̈ .- view(o.S✻UMZperp[ind2+1],S,:) * β̈v :
																																o.Repl.Z[S,ind2] * _β̈                                       ))
		end
	end
	nothing
end

function FillingLoop3!(o::StrBootTest{T}, T₁::Matrix{T}, ind1::Integer, ind2::Integer) where T
	Threads.@threads for i ∈ 1:o.N✻
		T₁[o.IDCT⋂✻[i], i] .+= o.S✻⋂u₁XinvXX[ind2+1,i] * view(o.Repl.XZ,:,ind1)
	end
	nothing
end

# Workhorse for WRE CRVE sandwich filling
# Given a column index within it, ind1, and a matrix β̈s of all the bootstrap estimates, 
# return all bootstrap realizations of P_X * Z[:,ind1]_g ' û₁g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per ⋂ cluster, one col per bootstrap replication
function Filling(o::StrBootTest{T}, ind1::Integer, β̈s::AbstractMatrix) where T
	if o.granular
   	if o.Nw == 1  # create or avoid NxB matrix?
			PXY✻ = reshape(o.PXZ[:,ind1], :, 1)  # store as matrix to reduce compiler confusion
			o.Repl.Yendog[ind1+1] && (PXY✻ = PXY✻ .+ o.S✻UPX[ind1+1] * o.v)

			dest = @panelsum(o, PXY✻ .* (o.Repl.y₁ .- o.S✻UMZperp[1] * o.v), o.info⋂)

			@inbounds for ind2 ∈ 1:o.Repl.kZ
				_β̈ = view(β̈s,ind2,:)'
				dest .-= @panelsum(o, PXY✻ .* (o.Repl.Yendog[ind2+1] ? view(o.Repl.Z,:,ind2) * _β̈ .- o.S✻UMZperp[ind2+1] * (o.v .* _β̈) :
															                                (view(o.Repl.Z,:,ind2) * _β̈)                                      ), o.info⋂)
			end
		else  # create pieces of each N x B matrix one at a time rather than whole thing at once
			dest = Matrix{T}(undef, o.clust[1].N, ncols(o.v))  # XXX preallocate this & turn Filling into Filling! ?
			@inbounds for ind2 ∈ 0:o.Repl.kZ
				ind2>0 && (β̈v = o.v .* (_β̈ = view(β̈s,ind2,:)'))

				if o.purerobust
					FillingLoop1!(o, dest, ind1, ind2, _β̈)
				else
					FillingLoop2!(o, dest, ind1, ind2, _β̈)
				end
			end
    end
  else  # coarse error clustering
		@inbounds for ind2 ∈ 0:o.Repl.kZ
			β̈v = iszero(ind2) ? o.v : o.v .* (_β̈ = -view(β̈s,ind2,:)')

			# T1 * o.v will be 1st-order terms
			T₁ = o.Repl.Yendog[ind1+1] ? o.S⋂YX[ind2+1] * o.invXXS✻XU[ind1+1] : Matrix{T}(undef,0,0)  #  S_∩ (Y_(∥j):*X_∥ ) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

			if o.Repl.Yendog[ind2+1]  # add CT_(⋂,*) (P_(X_par ) Y_(pari).*U ̈_(parj) )
				if o.NClustVar == o.nbootclustvar && iszero(o.subcluster)  # simple case of one clustering: full crosstab is diagonal
					tmp = view(o.Repl.XZ,:,ind1)
					if length(T₁)>0
						T₁[diagind(T₁)] .+= o.invXXS✻XU[ind2+1]'tmp
					else
						T₁                = reshape(o.invXXS✻XU[ind2+1]'tmp, :, 1)  # keep T₁ as vector rather than Diagonal matrix; probably better for fusion loop
					end
				else
					!o.Repl.Yendog[ind1+1] && (T₁ = o.JN⋂N✻)
					FillingLoop3!(o, T₁, ind1, ind2)
				end
				ncols(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.S⋂PXYZperp[ind1+1] * o.invZperpZperpS✻ZperpU[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
				o.NFE               > 0 && (T₁ = T₁ .- o.CT_FE⋂PY[ind1+1]'o.CTFEU[ind2+1])
			end

			if iszero(ind2)  # order-0 and -1 terms
				if iszero(ncols(T₁))
					dest = o.FillingT₀[ind1+1,1]
				elseif isone(ncols(T₁))
					dest = o.FillingT₀[ind1+1,1] .+ T₁ .* β̈v
				else
					dest = T₁ * o.v; dest .+= o.FillingT₀[ind1+1,1]
				end
			else  # y component
				if iszero(ncols(T₁))  # - x*β̈ components
					dest .+= o.FillingT₀[ind1+1,ind2+1] * _β̈
				elseif isone(ncols(T₁))
					dest .+= o.FillingT₀[ind1+1,ind2+1] * _β̈ .+ T₁ .* β̈v
				else
					dest .+= o.FillingT₀[ind1+1,ind2+1] * _β̈
					matmulplus!(Val(o.turbo), dest, T₁, β̈v)
				end
			end

			if o.Repl.Yendog[ind1+1] && o.Repl.Yendog[ind2+1]  # (S_*:*U ̈_(∥j)^'-S_* (U ̈_(∥j):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) Z_⊥^' )^' v^*
				for i ∈ 1:o.clust[1].N  # equivalent to this loop but slower: colquadformminus!(Val(o.turbo), dest, o.v, o.S✻UPX_S✻UMZperp[ind1+1,ind2+1], β̈v)
					colquadformminus!(Val(o.turbo), dest, i, o.v, o.S✻UPX_S✻UMZperp[ind1+1,ind2+1][:,i,:], β̈v)  # colquadformminus runs faster on array than view; convert to Array{3} of vectors
				end
			end
		end
  end
  dest
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
			_β̈s = (κs .* (YPXY₁₂ .- YY₁₂) .+ YY₁₂) ./ o.As
		else
			o.As = HessianFixedkappa(o, [1], 1, o.κ, w)
			_β̈s  = HessianFixedkappa(o, [1], 0, o.κ, w) ./ o.As
		end

		if o.null
			o.numerw = _β̈s .+ (o.Repl.Rt₁ - o.r) / o.Repl.RRpar
		else
			o.numerw = _β̈s .- o.DGP.β̈
			isone(w) && (o.numerw[1] = _β̈s[1] + (o.Repl.Rt₁[1] - o.r[1]) / o.Repl.RRpar[1])
		end

		@storeWtGrpResults!(o.numer, o.numerw)

		if o.bootstrapt
			if o.robust
				J⋂s = Filling(o, 1, _β̈s); J⋂s ./= o.As  # preallocate J⋂s
				@inbounds for c ∈ 1:o.NErrClustCombs  # sum sandwich over error clusterings
					nrows(o.clust[c].order)>0 && 
						(J⋂s = J⋂s[o.clust[c].order,:])
					@clustAccum!(denom, c, coldot(o, @panelsum(o, J⋂s, o.clust[c].info)))
				end
			else
				denom = (HessianFixedkappa(o, [0], 0, zero(T), w) .- 2 .* _β̈s .* HessianFixedkappa(o, [0], 1, zero(T), w) .+ _β̈s.^2 .* HessianFixedkappa(o, [1], 1, zero(T), w)) ./ o._Nobs ./ o.As  # classical error variance
			end
			@storeWtGrpResults!(o.dist, o.sqrt ? o.numerw ./ sqrt.(denom) : o.numerw .^ 2 ./ denom)
			denom *= o.Repl.RRpar[1]^2
		end

  else  # WRE bootstrap for more than 1 retained coefficient in bootstrap regression

		β̈s = zeros(T, o.Repl.kZ, ncols(o.v))
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
				β̈s[:,b] = (A[b] = invsym(o.κ*o.YPXY✻_b[2:end,2:end] + (1-o.κ)*o.YY✻_b[2:end,2:end])) * (o.κ*o.YPXY✻_b[1,2:end]' + (1-o.κ)*YY✻_b[1,2:end]')
			end
		else
			δnumer =  HessianFixedkappa(o, collect(1:o.Repl.kZ), 0, o.κ, w)
			δdenom = [HessianFixedkappa(o, collect(1:i), i, o.κ, w) for i ∈ 1:o.Repl.kZ]
			
			@inbounds #=Threads.@threads=# for b ∈ axes(o.v,2)
				for i ∈ 1:o.Repl.kZ
					o.δdenom_b[1:i,i] = view(δdenom[i],:,b)  # fill uppper triangle
				end
				β̈s[:,b] = (A[b] = invsym(o.δdenom_b)) * view(δnumer,:,b)
			end
		end

		if o.bootstrapt
			if o.robust
				@inbounds for i ∈ 1:o.Repl.kZ  # avoid list comprehension construction because of compiler type detection issue
					o.Zyg[i] = Filling(o, i, β̈s)
				end
			else
				YY✻ = [HessianFixedkappa(o, collect(i:o.Repl.kZ), i, zero(T), w) for i ∈ 0:o.Repl.kZ]  # κ=0 => Y*MZperp*Y
			end
		end

		@inbounds for b ∈ reverse(axes(o.v,2))
			if o.null || w==1 && b==1
				o.numer_b .= o.Repl.RRpar * view(β̈s,:,b) + o.Repl.Rt₁ - o.r
			else
				o.numer_b .= o.Repl.RRpar * (view(β̈s,:,b) - o.DGP.β̈₀)
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
					denom = (o.Repl.RRpar * A[b] * o.Repl.RRpar') * [-one(T) ; β̈s[:,b]]'Symmetric(o.YY✻_b) * [-one(T) ; β̈s[:,b]] / o._Nobs  # 2nd half is sig2 of errors
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
