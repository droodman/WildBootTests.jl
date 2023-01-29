# stuff done once per exucution--not depending on r
function InitWRE!(o::StrBootTest{T}) where T
	o.Repl.kZ>1 && (o.numer_b = Vector{T}(undef,nrows(o.Repl.RRpar)))

	o.J⋂s = isone(o.Nw) ? [Array{T,3}(undef, o.N⋂, ncols(o.v), o.Repl.kZ)] :
	                       [Array{T,3}(undef, o.N⋂, length(o.WeightGrp[1]), o.Repl.kZ), Array{T,3}(undef, o.N⋂, length(o.WeightGrp[end]), o.Repl.kZ)]


	o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.DGP.kX, ncols(o.v))] :
	                      [Matrix{T}(undef, o.DGP.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.DGP.kX, length(o.WeightGrp[end]))]
  o.T1R = deepcopy(o.T1L)
	o.β̈s = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kZ, ncols(o.v))] :
	                     [Matrix{T}(undef, o.Repl.kZ, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kZ, length(o.WeightGrp[end]))]
	o.As = isone(o.Nw) ? [Array{T,3}(undef, o.Repl.kZ, o.Repl.kZ, ncols(o.v))] :
	                     [Array{T,3}(undef, o.Repl.kZ, o.Repl.kZ, length(o.WeightGrp[1])), Array{T,3}(undef, o.Repl.kZ, o.Repl.kZ, length(o.WeightGrp[end]))]
	o.numerWRE = isone(o.Nw) ? [Matrix{T}(undef, o.dof, ncols(o.v))] :
			                        [Matrix{T}(undef, o.dof, length(o.WeightGrp[1])), Matrix{T}(undef, o.dof, length(o.WeightGrp[end]))]
	if isone(o.Repl.kZ)
		if o.liml
			o.YY₁₁ = deepcopy(o.β̈s)
			o.YY₁₂ = deepcopy(o.β̈s)
			o.YY₂₂   = deepcopy(o.β̈s)
			o.YPXY₁₁ = deepcopy(o.β̈s)
			o.YPXY₁₂ = deepcopy(o.β̈s)
			o.YPXY₂₂ = deepcopy(o.β̈s)
			o.YY₁₂YPXY₁₂ = deepcopy(o.β̈s)
			o.x₁₁ = deepcopy(o.β̈s)
			o.x₁₂ = deepcopy(o.β̈s)
			o.x₂₁ = deepcopy(o.β̈s)
			o.x₂₂ = deepcopy(o.β̈s)
			o.κs = deepcopy(o.β̈s)
		end
	else
		if o.liml
			o.YY✻ = isone(o.Nw) ? [[Matrix{T}(undef, i+1, ncols(o.v)) for i ∈ 0:o.Repl.kZ]] :
		                         [[Matrix{T}(undef, i+1, length(o.WeightGrp[1])) for i ∈ 0:o.Repl.kZ], 
														  [Matrix{T}(undef, i+1, length(o.WeightGrp[end])) for i ∈ 0:o.Repl.kZ]]
			o.YPXY✻ = deepcopy(o.YY✻)
		else
			o.δnumer = deepcopy(o.β̈s)
			o.δdenom = isone(o.Nw) ? [[Matrix{T}(undef, i, ncols(o.v))               for i ∈ 1:o.Repl.kZ]] :
		                           [[Matrix{T}(undef, i, length(o.WeightGrp[1]))   for i ∈ 1:o.Repl.kZ], 
															  [Matrix{T}(undef, i, length(o.WeightGrp[end])) for i ∈ 1:o.Repl.kZ]]
		end
	end

	if o.bootstrapt
		if o.liml || !o.robust
			o.YY✻_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
			o.YPXY✻_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
		end
		if o.NFE>0 && !o.FEboot && (o.bootstrapt || !isone(o.κ) || o.liml)
			o.CT✻FEu₁           = Array{T,3}(undef, o.NFE, o.N✻, 1)
			o.CT✻FEU₂par        = Array{T,3}(undef, o.NFE, o.N✻, o.Repl.kZ)
			o.invFEwtCT✻FEu₁    = Array{T,3}(undef, o.NFE, o.N✻, 1)
			o.invFEwtCT✻FEU₂par = Array{T,3}(undef, o.NFE, o.N✻, o.Repl.kZ)
			o.CT✻FEU        = [i>0 ? view(o.CT✻FEU₂par       ,:,:,i) : view(o.CT✻FEu₁       ,:,:,1) for i ∈ 0:o.Repl.kZ]
			o.invFEwtCT✻FEU = [i>0 ? view(o.invFEwtCT✻FEU₂par,:,:,i) : view(o.invFEwtCT✻FEu₁,:,:,1) for i ∈ 0:o.Repl.kZ]
		end
	end

	o.S✻Xu₁         = Array{T,3}(undef, o.DGP.kX, o.N✻, 1)
	o.S✻XU₂         = Array{T,3}(undef, o.DGP.kX, o.N✻, o.kY₂)
	o.S✻XU₂par      = Array{T,3}(undef, o.DGP.kX, o.N✻, o.Repl.kZ) 
	o.invXXS✻Xu₁    = Array{T,3}(undef, o.DGP.kX, o.N✻, 1)
	o.invXXS✻XU₂par = Array{T,3}(undef, o.DGP.kX, o.N✻, o.Repl.kZ)
	o.S✻XU      = [i>0 ? view(o.S✻XU₂par     ,:,:,i) : view(o.S✻Xu₁     ,:,:,1) for i ∈ 0:o.Repl.kZ]
	o.invXXS✻XU = [i>0 ? view(o.invXXS✻XU₂par,:,:,i) : view(o.invXXS✻Xu₁,:,:,1) for i ∈ 0:o.Repl.kZ]

	if o.bootstrapt || o.liml || !isone(o.κ)
		o.S✻Zperpu₁                 = Array{T,3}(undef, o.DGP.kZperp, o.N✻, 1)
		o.S✻ZperpU₂par              = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.Repl.kZ)
		o.invZperpZperpS✻Zperpu₁    = Array{T,3}(undef, o.DGP.kZperp, o.N✻, 1)
		o.invZperpZperpS✻ZperpU₂par = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.Repl.kZ)
		o.S✻ZperpU              = [i>0 ? view(o.S✻ZperpU₂par             ,:,:,i) : view(o.S✻Zperpu₁             ,:,:,1) for i ∈ 0:o.Repl.kZ]
		o.invZperpZperpS✻ZperpU = [i>0 ? view(o.invZperpZperpS✻ZperpU₂par,:,:,i) : view(o.invZperpZperpS✻Zperpu₁,:,:,1) for i ∈ 0:o.Repl.kZ]

		o.S✻u₁u₁       = Array{T,3}(undef, 1        , o.N✻, 1        )
		o.S✻U₂paru₁    = Array{T,3}(undef, o.Repl.kZ, o.N✻, 1        )
		o.S✻U₂parU₂par = Array{T,3}(undef, o.Repl.kZ, o.N✻, o.Repl.kZ)
		o.S✻UU = [i>0 ? j>0 ? view(o.S✻U₂parU₂par,i,:,j) : view(o.S✻U₂paru₁,i,:,1) : j>0 ? view(o.S✻U₂paru₁,j,:,1) : view(o.S✻u₁u₁,1,:,1) for i ∈ 0:o.Repl.kZ, j ∈ 0:o.Repl.kZ]

		o.S✻y₁paru₁     = Array{T,3}(undef, 1        , o.N✻, 1        )
		o.S✻Zu₁         = Array{T,3}(undef, o.Repl.kZ, o.N✻, 1        )
		o.S✻y₁parU₂par  = Array{T,3}(undef, 1        , o.N✻, o.Repl.kZ)
		o.S✻ZU₂par      = Array{T,3}(undef, o.Repl.kZ, o.N✻, o.Repl.kZ)
		o.S✻YU = [i>0 ? j>0 ? view(o.S✻ZU₂par,i,:,j) : view(o.S✻Zu₁,i,:,1) : j>0 ? view(o.S✻y₁parU₂par,1,:,j) : view(o.S✻y₁paru₁,1,:,1) for i ∈ 0:o.Repl.kZ, j ∈ 0:o.Repl.kZ]
		o.S✻YUfold = Array{T,3}(undef, o.Repl.kZ+1, o.N✻, o.Repl.kZ+1)
	end

	if o.granular
		if o.bootstrapt && o.robust
			o.S✻UMZperp = Array{T,3}(undef, o.Nobs, o.N✻, o.Repl.kZ+1)
			o.S✻UPX     = Array{T,3}(undef, o.Nobs, o.N✻, o.Repl.kZ  )
			o.crosstab✻ind = o.Nobs==o.N✻ ? Vector(diagind(FakeArray(o.N✻,o.N✻))) : LinearIndices(FakeArray(o.Nobs,o.N✻))[CartesianIndex.(1:o.Nobs, o.ID✻)]
			o.XinvXX = X₁₂B(o.Repl.X₁, o.Repl.X₂, o.Repl.invXX)
			o.PXZ    = X₁₂B(o.Repl.X₁, o.Repl.X₂, o.Repl.V)
		end
	else
		o.Repl.Zperp = o.DGP.Zperp = Matrix{T}(undef,0,0)  # drop this potentially large array

		o.S✻⋂XU₂     = Array{T,3}(undef, o.DGP.kX, o.N✻⋂, o.kY₂)
		o.S✻⋂XU₂par  = Array{T,3}(undef, o.DGP.kX, o.N✻⋂, o.Repl.kZ)
		o.invXXS✻XU₂ = Array{T,3}(undef, o.DGP.kX, o.N✻, o.kY₂)

		if o.bootstrapt || o.liml || !isone(o.κ)  
			o.S✻ZperpU₂ = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.kY₂)
			o.invZperpZperpS✻ZperpU₂ = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.kY₂)
		end

		o.bootstrapt && o.robust &&
			(o.negS✻UMZperpX = [Array{T,3}(undef, o.DGP.kX, o.N⋂, o.N✻) for _ in 0:o.Repl.kZ])

		o.S✻⋂XY₂      = o.DGP.S✻⋂XY₂     - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpY₂   - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpY₂ )
		o.S✻⋂XX       = o.DGP.S✻⋂XX      - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpX    - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂XZperp'   - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpX  )
		o.S✻⋂XDGPZ    = o.DGP.S✻⋂XZpar   - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZpar - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZpar)
		o.S✻⋂Xy₁      = o.DGP.S✻⋂Xy₁     - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂Zperpy₁   - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpy₁ )
		  S✻⋂ZperpX   = o.DGP.S✻⋂XZperp' - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpX
		o.DGP.restricted &&
			(o.S✻⋂X_DGPZR₁ = o.DGP.S✻⋂XZR₁ - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpZR₁  - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZR₁  - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZR₁ ))

		o.invXXS✻⋂XY₂    = o.DGP.invXX * o.S✻⋂XY₂ 
		o.invXXS✻⋂XX     = o.DGP.invXX * o.S✻⋂XX  
		o.invXXS✻⋂XDGPZ  = o.DGP.invXX * o.S✻⋂XDGPZ 
		o.invXXS✻⋂Xy₁    = o.DGP.invXX * o.S✻⋂Xy₁ 
		o.DGP.restricted &&
			(o.invXXS✻⋂XDGPZR₁ = o.DGP.invXX * o.S✻⋂X_DGPZR₁)

		_S✻ZperpY₂      = @panelsum(o.DGP.S✻⋂ZperpY₂, o.info✻_✻⋂)  # moments of variables _before_ FWL processing
		_S✻Zperpy₁      = @panelsum(o.DGP.S✻⋂Zperpy₁, o.info✻_✻⋂)
		_S✻ZperpDGPZpar = @panelsum(o.DGP.S✻⋂ZperpZpar, o.info✻_✻⋂)
		o.DGP.restricted &&
			(_S✻ZperpDGPZR₁  = @panelsum(o.DGP.S✻⋂ZperpZR₁, o.info✻_✻⋂))

		S✻ZperpZperp    = @panelsum(o.DGP.S✻⋂ZperpZperp, o.info✻_✻⋂)
		o.S✻XY₂         = @panelsum(o.S✻⋂XY₂   , o.info✻_✻⋂)
		o.S✻XX          = @panelsum(o.S✻⋂XX    , o.info✻_✻⋂)
		o.S✻XDGPZ       = @panelsum(o.S✻⋂XDGPZ, o.info✻_✻⋂)
		o.S✻Xy₁         = @panelsum(o.S✻⋂Xy₁, o.info✻_✻⋂)
		o.S✻ZperpX      = @panelsum(S✻⋂ZperpX, o.info✻_✻⋂)
		o.S✻ZperpY₂     = _S✻ZperpY₂ - S✻ZperpZperp * o.DGP.invZperpZperpZperpY₂
		o.S✻ZperpDGPZ   = _S✻ZperpDGPZpar - S✻ZperpZperp * o.DGP.invZperpZperpZperpZpar
		o.S✻Zperpy₁     = _S✻Zperpy₁ - S✻ZperpZperp * o.DGP.invZperpZperpZperpy₁
		if o.DGP.restricted
			o.S✻XZR₁        = @panelsum(o.S✻⋂X_DGPZR₁, o.info✻_✻⋂)
			o.S✻ZperpDGPZR₁ = @panelsum(o.DGP.S✻⋂ZperpZR₁ , o.info✻_✻⋂) - S✻ZperpZperp * o.DGP.invZperpZperpZperpZR₁
		end

		if o.NFE>0 && !o.FEboot && (o.liml || !isone(o.κ) || o.bootstrapt)
			  CT✻⋂FEX  = [crosstabFE(o, o.DGP.X₁, o.info✻⋂) crosstabFE(o, o.DGP.X₂, o.info✻⋂)]
			o.CT✻FEX   = @panelsum(CT✻⋂FEX, o.info✻_✻⋂)
			o.CT✻FEY₂  = crosstabFE(o, o.DGP.Y₂, o.info✻)
			o.CT✻FEZ   = crosstabFE(o, o.DGP.Z, o.info✻)
			o.CT✻FEy₁  = crosstabFE(o, o.DGP.y₁, o.info✻)
			o.DGP.restricted &&
				(o.CT✻FEZR₁ = crosstabFE(o, o.DGP.ZR₁, o.info✻))
		end

		((o.robust && o.bootstrapt) || o.liml || !o.robust || !isone(o.κ)) &&
			(S✻⋂ReplZX = (o.Repl.S✻⋂XZpar - o.DGP.S✻⋂XZperp * o.Repl.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpX' * (o.Repl.S✻⋂ZperpZpar - o.DGP.S✻⋂ZperpZperp * o.Repl.invZperpZperpZperpZpar))')

		if o.bootstrapt && o.robust
			o.info⋂_✻⋂ = panelsetup(o.ID✻⋂, o.subcluster+1:o.NClustVar)

			o.S⋂ReplZX = @panelsum(S✻⋂ReplZX, o.info⋂_✻⋂)
			S⋂ZperpX   = @panelsum(S✻⋂ZperpX, o.info⋂_✻⋂)
			o.S⋂Xy₁    = @panelsum(o.S✻⋂Xy₁, o.info⋂_✻⋂)

			_inds = o.subcluster>0 ?
							[CartesianIndex(j,i) for (j,v) ∈ enumerate(o.info⋂_✻⋂) for i ∈ v] :  # crosstab ∩,* is wide
							o.NClustVar == o.NBootClustVar ?
									[CartesianIndex(i,i) for i ∈ 1:o.N✻⋂] :  # crosstab *,∩ is square
									[CartesianIndex(i,j) for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab ∩,* is tall
			inds = [CartesianIndex(k,I) for I ∈ _inds for k ∈ 1:o.DGP.kX]
			o.crosstab⋂✻ind = LinearIndices(FakeArray(Tuple(max(inds...))...))[inds]

			o.S⋂XZperpinvZperpZperp = S⋂ZperpX' * o.DGP.invZperpZperp

			o.Q    = Array{T,3}(undef, o.N✻, o.N⋂, o.N✻)
			o.NFE>0 && !o.FEboot && (o.CT⋂FEX = o.invFEwt .* @panelsum(CT✻⋂FEX, o.info⋂_✻⋂))

			o.β̈v = isone(o.Nw) ? [Matrix{T}(undef, o.N✻, ncols(o.v))] :
			                     [Matrix{T}(undef, o.N✻, length(o.WeightGrp[1])), Matrix{T}(undef, o.N✻, length(o.WeightGrp[end]))]
		end

		if o.liml || !o.robust || !isone(o.κ)  # cluster-wise moments after FWL
			o.S✻Y₂Y₂     = o.DGP.S✻Y₂Y₂     - _S✻ZperpY₂'      * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpY₂'   * o.S✻ZperpY₂
			o.S✻DGPZDGPZ = o.DGP.S✻ZparZpar - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpZpar' * o.S✻ZperpDGPZ
			o.S✻DGPZY₂   = o.DGP.S✻ZparY₂   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpZpar' * o.S✻ZperpY₂
			o.S✻DGPZy₁   = o.DGP.S✻Zpary₁   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpZpar' * o.S✻Zperpy₁   
			o.S✻Y₂y₁     = o.DGP.S✻Y₂y₁     - _S✻ZperpY₂'      * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpY₂'   * o.S✻Zperpy₁
			o.S✻y₁y₁     = o.DGP.S✻y₁y₁     - _S✻Zperpy₁'      * o.DGP.invZperpZperpZperpy₁   - o.S✻Zperpy₁'                  * o.DGP.invZperpZperpZperpy₁
			o.DGP.restricted && 
				(o.S✻DGPZR₁y₁ = o.DGP.S✻ZR₁y₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpy₁ - o.DGP.invZperpZperpZperpZR₁' * o.S✻Zperpy₁)

			if o.Repl.restricted
				_S✻ZperpReplZR₁ = @panelsum(o.Repl.S✻⋂ZperpZR₁, o.info✻_✻⋂)
				_S✻⋂XReplZR₁    = @panelsum(o.Repl.S✻⋂XZR₁    , o.info✻_✻⋂)
				
				o.r₁S✻ReplZR₁Y₂     = o.r₁' * (o.Repl.S✻ZR₁Y₂ - _S✻ZperpReplZR₁' * o.DGP.invZperpZperpZperpY₂ - o.Repl.invZperpZperpZperpZR₁' * o.S✻ZperpY₂)  
				o.r₁S✻ReplZR₁X      = o.r₁' * (_S✻⋂XReplZR₁'  - _S✻ZperpReplZR₁' * o.DGP.invZperpZperpZperpX  - o.Repl.invZperpZperpZperpZR₁' * o.S✻ZperpX )
				o.r₁S✻ReplZR₁DGPZ   = o.r₁' * panelcross(o.Repl.ZR₁, o.DGP.Z, o.info✻)   
				o.r₁S✻ReplZR₁y₁     = o.r₁' * (o.Repl.S✻ZR₁y₁ - _S✻ZperpReplZR₁' * o.DGP.invZperpZperpZperpy₁ - o.Repl.invZperpZperpZperpZR₁' * o.S✻Zperpy₁)
				o.DGP.restricted &&
					(o.r₁S✻ReplZR₁DGPZR₁ = o.r₁' * panelcross(o.Repl.ZR₁, o.DGP.ZR₁, o.info✻))
			end

			_S✻ZperpReplZpar = @panelsum(o.Repl.S✻⋂ZperpZpar, o.info✻_✻⋂)
			_S✻ReplXZ        = @panelsum(o.Repl.S✻⋂XZpar    , o.info✻_✻⋂)

			o.S✻ReplZY₂      = o.Repl.S✻ZparY₂ - _S✻ZperpReplZpar' * o.DGP.invZperpZperpZperpY₂ - o.Repl.invZperpZperpZperpZpar' * o.S✻ZperpY₂
			o.S✻ReplZX       = _S✻ReplXZ'      - _S✻ZperpReplZpar' * o.DGP.invZperpZperpZperpX  - o.Repl.invZperpZperpZperpZpar' * o.S✻ZperpX
			o.S✻ReplZDGPZ    = panelcross(o.Repl.Z, o.DGP.Z, o.info✻)   
			o.S✻ReplZy₁      = o.Repl.S✻Zpary₁ - _S✻ZperpReplZpar' * o.DGP.invZperpZperpZperpy₁ - o.Repl.invZperpZperpZperpZpar' * o.S✻Zperpy₁
			
			if o.DGP.restricted
				_S✻⋂XDGPZR₁ = @panelsum(o.DGP.S✻⋂XZR₁, o.info✻_✻⋂)

				o.S✻ReplZDGPZR₁  = panelcross(o.Repl.Z, o.DGP.ZR₁, o.info✻)   
				o.S✻DGPZR₁Y₂     = o.DGP.S✻ZR₁Y₂  - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpY₂  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpY₂
				o.S✻DGPZR₁DGPZR₁ = o.DGP.S✻ZR₁ZR₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZR₁  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZR₁
				o.S✻DGPZR₁DGPZ   = o.DGP.S✻ZR₁Z   - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZpar - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZ
				o.S✻DGPZR₁X      = _S✻⋂XDGPZR₁'   - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpX   - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpX
			end
		end

		o.invXXS✻XY₂   = @panelsum(o.invXXS✻⋂XY₂  , o.info✻_✻⋂)
		o.invXXS✻XX    = @panelsum(o.invXXS✻⋂XX   , o.info✻_✻⋂)
		o.invXXS✻XDGPZ = @panelsum(o.invXXS✻⋂XDGPZ, o.info✻_✻⋂)
		o.invXXS✻Xy₁   = @panelsum(o.invXXS✻⋂Xy₁, o.info✻_✻⋂)
		o.invZperpZperpS✻ZperpY₂   = o.DGP.invZperpZperp * o.S✻ZperpY₂ 
		o.invZperpZperpS✻ZperpX    = o.DGP.invZperpZperp * o.S✻ZperpX  
		o.invZperpZperpS✻Zperpy₁   = o.DGP.invZperpZperp * o.S✻Zperpy₁ 
		o.invZperpZperpS✻ZperpDGPZ = o.DGP.invZperpZperp * o.S✻ZperpDGPZ

		if o.DGP.restricted
			o.invXXS✻XDGPZR₁ = @panelsum(o.invXXS✻⋂XDGPZR₁, o.info✻_✻⋂)
			o.invZperpZperpS✻ZperpDGPZR₁ = o.DGP.invZperpZperp * o.S✻ZperpDGPZR₁
		end
	end
	nothing
end

function PrepWRE!(o::StrBootTest{T}) where T
	if o.null
		r₁ = [o.r₁ ; o.r]
	  EstimateIV!(o.DGP, o, o.jk, r₁)
	  MakeResidualsIV!(o.DGP, o)
	else
		r₁ = o.r₁
	end

	_β̈  = view(o.DGP.β̈  ,:,1)

	tmul!(o.S✻⋂XU₂, o.S✻⋂XX, o.DGP.Π̂ ) ; o.S✻⋂XU₂ .= o.S✻⋂XY₂ .- o.S✻⋂XU₂
	tmul!(o.S✻⋂XU₂par, o.S✻⋂XU₂, o.Repl.RparY)
	tmul!(o.S✻XU₂, o.S✻XX, o.DGP.Π̂ ); o.S✻XU₂ .= o.S✻XY₂ .- o.S✻XU₂
	tmul!(o.S✻XU₂par, o.S✻XU₂, o.Repl.RparY)
	tmul!(o.invXXS✻XU₂, o.DGP.invXX, o.S✻XU₂)
	tmul!(o.invXXS✻XU₂par, o.invXXS✻XU₂, o.Repl.RparY)
	if o.bootstrapt || o.liml || !isone(o.κ)
		tmul!(o.S✻ZperpU₂, o.S✻ZperpX, o.DGP.Π̂ ); o.S✻ZperpU₂ .= o.S✻ZperpY₂ .- o.S✻ZperpU₂
		o.invZperpZperpS✻ZperpU₂ .= o.invZperpZperpS✻ZperpY₂ - o.invZperpZperpS✻ZperpX * o.DGP.Π̂
		tmul!(o.S✻ZperpU₂par, o.S✻ZperpU₂, o.Repl.RparY)
		tmul!(o.invZperpZperpS✻ZperpU₂par, o.invZperpZperpS✻ZperpU₂, o.Repl.RparY)
	end

	o.S✻Xu₁ .= o.S✻Xy₁ .- o.S✻XDGPZ * _β̈  .+ o.S✻XU₂ * o.DGP.γ̈ 
	o.DGP.restricted &&
		(o.S✻Xu₁ .-= o.S✻XZR₁ * r₁)
	@panelsum!(o.S✻XU₂par, o.S✻⋂XU₂par, o.info✻_✻⋂)

	o.invXXS✻Xu₁ .= o.invXXS✻Xy₁ .- o.invXXS✻XDGPZ * _β̈  .+ o.invXXS✻XU₂ * o.DGP.γ̈ 
	o.DGP.restricted &&
		(o.invXXS✻Xu₁ .-= o.invXXS✻XDGPZR₁ * r₁)

	if o.liml || !isone(o.κ) || !o.robust
		S✻U₂y₁ = o.S✻Y₂y₁ - o.DGP.Π̂' * o.S✻Xy₁
		S✻ZU₂ = o.S✻ReplZY₂ - o.S✻ReplZX * o.DGP.Π̂
		o.S✻ZU₂par .= S✻ZU₂ * o.Repl.RparY
		Π̂S✻XÜ₂γ̈ = o.DGP.Π̂' * o.S✻XU₂ * o.DGP.γ̈
		S✻Ü₂Y₂ = o.S✻Y₂Y₂ - o.DGP.Π̂' * o.S✻XY₂
		S✻Y₂Ü₂γ̈ = S✻Ü₂Y₂' * o.DGP.γ̈

		S✻UUterm = o.S✻Y₂y₁ - o.S✻DGPZY₂' * view(o.DGP.β̈ ,:,1) - o.DGP.Π̂' * (o.S✻Xy₁  - o.S✻XDGPZ * o.DGP.β̈)
		o.S✻u₁u₁ .= o.S✻y₁y₁ .- (2 * o.S✻DGPZy₁ - o.S✻DGPZDGPZ * _β̈ )'_β̈  .+ (2 * S✻UUterm - Π̂S✻XÜ₂γ̈ + S✻Y₂Ü₂γ̈ )'o.DGP.γ̈ 
		o.S✻U₂paru₁ .= o.Repl.RparY' * (S✻UUterm + S✻Y₂Ü₂γ̈ - Π̂S✻XÜ₂γ̈ )
		o.S✻U₂parU₂par .= o.Repl.RparY' * (S✻Ü₂Y₂ - o.S✻XU₂' * o.DGP.Π̂) * o.Repl.RparY

		o.S✻y₁paru₁ .= o.S✻y₁y₁ - o.S✻DGPZy₁'_β̈ + S✻U₂y₁'o.DGP.γ̈

		o.S✻Zu₁ .= o.S✻ReplZy₁ - o.S✻ReplZDGPZ * _β̈ + S✻ZU₂ * o.DGP.γ̈

		o.S✻y₁parU₂par .= S✻U₂y₁' * o.Repl.RparY

		if o.DGP.restricted
			r₁S✻DGPZR₁y₁ = r₁' * o.S✻DGPZR₁y₁
			o.S✻u₁u₁ .+= -2 .* r₁S✻DGPZR₁y₁ .+ r₁' * (o.S✻DGPZR₁DGPZR₁ * r₁) .+ 2 .* r₁' * (o.S✻DGPZR₁DGPZ * _β̈ + (o.S✻DGPZR₁X * o.DGP.Π̂ - o.S✻DGPZR₁Y₂) * o.DGP.γ̈ )
			o.S✻U₂paru₁ .-= o.Repl.RparY' * (o.S✻DGPZR₁Y₂ - o.S✻DGPZR₁X * o.DGP.Π̂ )' * r₁
			o.S✻y₁paru₁ .-= r₁S✻DGPZR₁y₁
			o.S✻Zu₁ .-= o.S✻ReplZDGPZR₁ * r₁
		end

		if o.Repl.restricted
			S✻ReplZR₁r₁U₂ = o.r₁S✻ReplZR₁Y₂ - o.r₁S✻ReplZR₁X * o.DGP.Π̂
			o.S✻y₁paru₁ .-= o.r₁S✻ReplZR₁y₁ - o.r₁S✻ReplZR₁DGPZ * _β̈ + S✻ReplZR₁r₁U₂ * o.DGP.γ̈ 
			o.S✻y₁parU₂par .-= S✻ReplZR₁r₁U₂ * o.Repl.RparY
			o.DGP.restricted &&
				(o.S✻y₁paru₁ .+= o.r₁S✻ReplZR₁DGPZR₁ * r₁)
		end

		@inbounds for i ∈ 0:o.Repl.kZ, j ∈ 0:i
			o.S✻YUfold[i+1,:,j+1] .= o.S✻YU[i+1,j+1] .+ o.S✻YU[j+1,i+1]
			o.S✻YUfold[j+1,:,i+1] .= o.S✻YUfold[i+1,:,j+1]
		end
	end

	if o.liml || !isone(o.κ) || o.bootstrapt
		o.S✻Zperpu₁              .= o.S✻Zperpy₁              - o.S✻ZperpDGPZ              * _β̈ + o.S✻ZperpU₂              * o.DGP.γ̈ 
		o.invZperpZperpS✻Zperpu₁ .= o.invZperpZperpS✻Zperpy₁ - o.invZperpZperpS✻ZperpDGPZ * _β̈ + o.invZperpZperpS✻ZperpU₂ * o.DGP.γ̈ 
		if o.DGP.restricted
			o.S✻Zperpu₁              .-= o.S✻ZperpDGPZR₁              * r₁
			o.invZperpZperpS✻Zperpu₁ .-= o.invZperpZperpS✻ZperpDGPZR₁ * r₁
		end

		if o.NFE>0 && !o.FEboot
			CT✻FEU₂par = o.CT✻FEY₂ - o.CT✻FEX * o.DGP.Π̂
			o.CT✻FEu₁ .= o.CT✻FEy₁ - o.CT✻FEZ * _β̈ + CT✻FEU₂par * o.DGP.γ̈ 
			o.CT✻FEU₂par .= CT✻FEU₂par * o.Repl.RparY
			o.DGP.restricted &&
				(o.CT✻FEu₁ .-= o.CT✻FEZR₁ * r₁)
			o.invFEwtCT✻FEu₁    .= o.invFEwt .* o.CT✻FEu₁
			o.invFEwtCT✻FEU₂par .= o.invFEwt .* o.CT✻FEU₂par
		end
	end

	if o.robust && o.bootstrapt
		@inbounds for j ∈ 0:o.Repl.kZ
			if o.Repl.Yendog[j+1]
				o.negS✻UMZperpX[j+1] .= o.S⋂XZperpinvZperpZperp * o.S✻ZperpU[j+1]  # S_* diag⁡(U ̈_(∥j) ) Z_⊥ (Z_⊥^' Z_⊥ )^(-1) Z_(⊥g)^' X_(∥g)
				if iszero(j)  # - S_*  diag⁡(U ̈_(∥j) ) I_g^' X_(∥g)
					o.negS✻UMZperpX[j+1][o.crosstab⋂✻ind] .-= vec(o.S✻⋂Xy₁ - o.S✻⋂XDGPZ * _β̈ + o.S✻⋂XU₂ * o.DGP.γ̈)
					o.DGP.restricted &&
						(o.negS✻UMZperpX[j+1][o.crosstab⋂✻ind] .+= vec(o.S✻⋂X_DGPZR₁ * r₁))
				else
					o.negS✻UMZperpX[j+1][o.crosstab⋂✻ind] .-= vec(view(o.S✻⋂XU₂par,:,:,j))
				end
				o.NFE>0 && !o.FEboot &&
					(o.negS✻UMZperpX[j+1] .+= o.CT⋂FEX' * o.CT✻FEU[j+1])  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
			end
		end
	end
	nothing
end

function PrepWREGranular!(o::StrBootTest{T}) where T
	if o.null
		r₁ = [o.r₁ ; o.r]
	  EstimateIV!(o.DGP, o, o.jk, r₁)
	  MakeResidualsIV!(o.DGP, o)
  	o.Ü₂par = o.DGP.Ü₂[1] * o.Repl.RparY
	else
		r₁ = o.r₁
	end

	panelcross!(o.S✻Xu₁, o.DGP.X₁, o.DGP.X₂, o.DGP.u⃛₁[1], o.info✻)
	panelcross!(o.S✻XU₂par, o.DGP.X₁, o.DGP.X₂, o.Ü₂par, o.info✻)
	tmul!(o.invXXS✻Xu₁   , o.DGP.invXX, o.S✻Xu₁   )
	tmul!(o.invXXS✻XU₂par, o.DGP.invXX, o.S✻XU₂par)

	if o.bootstrapt || o.liml || !isone(o.κ)
		panelcross!(o.S✻Zperpu₁, o.DGP.Zperp, o.DGP.u⃛₁[1], o.info✻)
		panelcross!(o.S✻ZperpU₂par, o.DGP.Zperp, o.Ü₂par, o.info✻)
		tmul!(o.invZperpZperpS✻Zperpu₁, o.DGP.invZperpZperp, o.S✻Zperpu₁)
		tmul!(o.invZperpZperpS✻ZperpU₂par, o.DGP.invZperpZperp, o.S✻ZperpU₂par)
		if o.NFE>0 && !o.FEboot
			crosstabFE!(o, o.CT✻FEu₁    , o.DGP.u⃛₁[1], o.info✻)
			crosstabFE!(o, o.CT✻FEU₂par, o.Ü₂par     , o.info✻)
			o.invFEwtCT✻FEu₁ .= o.invFEwt .* o.CT✻FEu₁
			o.invFEwtCT✻FEU₂par .= o.invFEwt .* o.CT✻FEU₂par
		end
	end

	if o.liml || !o.robust || !isone(o.κ)
		panelcross!(o.S✻u₁u₁, o.DGP.u⃛₁[1], o.DGP.u⃛₁[1], o.info✻)
		panelcross!(o.S✻U₂paru₁, o.Ü₂par, o.DGP.u⃛₁[1], o.info✻)
		panelcross!(o.S✻U₂parU₂par, o.Ü₂par, o.Ü₂par, o.info✻)

		panelcross!(o.S✻y₁paru₁, o.Repl.y₁par, o.DGP.u⃛₁[1], o.info✻)
		panelcross!(o.S✻Zu₁, o.Repl.Z, o.DGP.u⃛₁[1], o.info✻)
		panelcross!(o.S✻y₁parU₂par, o.Repl.y₁par, o.Ü₂par, o.info✻)
		panelcross!(o.S✻ZU₂par, o.Repl.Z, o.Ü₂par, o.info✻)

		@inbounds for i ∈ 0:o.Repl.kZ, j ∈ 0:i
			o.S✻YUfold[i+1,:,j+1] .= o.S✻YU[i+1,j+1] .+ o.S✻YU[j+1,i+1]
			o.S✻YUfold[j+1,:,i+1] .= o.S✻YUfold[i+1,:,j+1]
		end
	end

	if o.robust && o.bootstrapt
		tmul!(o.S✻UMZperp[:,:,1    ][:,:,:], o.DGP.Zperp, o.invZperpZperpS✻Zperpu₁   )
		tmul!(o.S✻UMZperp[:,:,2:end], o.DGP.Zperp, o.invZperpZperpS✻ZperpU₂par)
		tmul!(o.S✻UPX, o.XinvXX, o.S✻XU₂par)
		@inbounds for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
			o.S✻UMZperp[:,:,i+1][o.crosstab✻ind] .-= i>0 ? view(o.Ü₂par,:,i) : view(o.DGP.u⃛₁[1], :)  # subtract crosstab of observation by ∩-group of u
			o.NFE>0 && !o.FEboot &&
				(o.S✻UMZperp[:,:,i+1] .-= view(o.invFEwtCT✻FEU[i+1], o._FEID, :))  # CT_(*,FE) (U ̈_(parj) ) (S_FE S_FE^' )^(-1) S_FE
		end
  end
	nothing
end


# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, i, j, return all bootstrap realizations of 
# Y[:,i]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,j] for κ constant across replications
# i can be a rowvector
# (only really the Hessian when we narrow Y to Z)
function HessianFixedkappa(o::StrBootTest{T}, is::Vector{S} where S<:Integer, j::Integer, κ::Number, _w::Integer) where T
  dest = Matrix{T}(undef, length(is), ncols(o.v))
  @inbounds for i ∈ eachindex(is, axes(dest,1))
		_HessianFixedkappa!(o, dest, i, is[i], j, κ, _w)
  end
  dest
end
function HessianFixedkappa!(o::StrBootTest{T}, dest::AbstractMatrix{T}, is::Vector{S} where S<:Integer, j::Integer, κ::Number, _w::Integer) where T
  @inbounds for i ∈ eachindex(is, axes(dest,1))
		_HessianFixedkappa!(o, dest, i, is[i], j, κ, _w)
  end
  dest
end

function _HessianFixedkappa!(o::StrBootTest, dest::AbstractMatrix, row::Integer, i::Integer, j::Integer, κ::Number, _w::Integer)
  if !(o.Repl.Yendog[i+1] || o.Repl.Yendog[j+1])  # if both vars exog, result = order-0 term only, same for all draws
		!iszero(κ) && 
			(dest[row,:] .= dot(view(o.Repl.XZ,:,i), view(o.Repl.V,:,j)))
		if !isone(κ)
			if iszero(κ)
				dest[row,:] .= o.Repl.YY[i+1,j+1]
			else
				dest[row,:] .= κ .* dest[row,:] .+ (1 - κ) .* o.Repl.YY[i+1,j+1]
			end
		end
	else
		if !iszero(κ)  # repetitiveness in this section to maintain type stability
			if o.Repl.Yendog[i+1]
				T1L = o.T1L[_w]  # preallocated destinations
				tmul!(T1L, o.S✻XU[i+1], o.v)
				if iszero(i)
					T1L .+= o.Repl.Xy₁par
				else
					T1L .+= view(o.Repl.XZ,:,i)
				end
				if o.Repl.Yendog[j+1]
					T1R = o.T1R[_w]  # preallocated destinations
					tmul!(T1R, o.invXXS✻XU[j+1], o.v)
					if iszero(j)
						T1R .+=  o.Repl.invXXXy₁par
					else
						T1R .+= view(o.Repl.V,:,j)
					end
					coldot!(dest, row, T1L, T1R)
				else
					dest[row,:] .= T1L'view(o.Repl.V,:,j)  # coldot!(dest, row, T1L, view(o.Repl.V,:,j))
				end
			else
				if o.Repl.Yendog[j+1]
					T1R = o.T1R[_w]  # use preallocated destinations
					tmul!(T1R, o.invXXS✻XU[j+1], o.v)
					if iszero(j)
						T1R .+=  o.Repl.invXXXy₁par
					else
						T1R .+= view(o.Repl.V,:,j)
					end
					dest[row,:] .= T1R'view(o.Repl.XZ,:,i)
				else
					dest[row,:] .= dot(view(o.Repl.V,:,j), view(o.Repl.XZ,:,i))
				end
			end
		end
		if !isone(κ)
			if o.Repl.Yendog[i+1]
				if iszero(κ)
					dest[row,:] .= o.Repl.YY[i+1,j+1] .+ o.v'view(o.S✻YUfold,i+1,:,j+1)
					coldotminus!(dest, row, o.invZperpZperpS✻ZperpU[i+1] * o.v, o.S✻ZperpU[j+1] * o.v)  # when is this term 0??
					coldotplus!(dest, row, o.v, o.S✻UU[i+1,j+1], o.v)
					o.NFE>0 && !o.FEboot &&
						coldotminus!(dest, row, o.CT✻FEU[i+1] * o.v, o.invFEwtCT✻FEU[j+1] * o.v)
				else
					_dest = o.Repl.YY[i+1,j+1] .+ view(o.S✻YUfold,i+1,:,j+1)'o.v
					coldotminus!(_dest, 1, o.invZperpZperpS✻ZperpU[i+1] * o.v, o.S✻ZperpU[j+1] * o.v)
					coldotplus!(_dest, 1, o.v, o.S✻UU[i+1, j+1], o.v)
					o.NFE>0 && !o.FEboot &&
						coldotminus!(_dest, 1, o.CT✻FEU[i+1] * o.v, o.invFEwtCT✻FEU[j+1] * o.v)
					dest[row,:] .= κ .* dest[row,:] .+ (1 - κ) .* _dest
				end
			elseif iszero(κ)
				dest[row,:] .= o.Repl.YY[i+1,j+1]
			else
				dest[row,:] .= κ .* dest[row,:] .+ (1 - κ) .* o.Repl.YY[i+1,j+1]
			end
		end
  end
	nothing
end

# put threaded loops in functions to prevent compiler-perceived type instability https://discourse.julialang.org/t/type-inference-with-threads/2004/3
function FillingLoop1!(o::StrBootTest{T}, dest::Matrix{T}, i::Integer, j::Integer, _β̈::AbstractMatrix{T}, β̈v::AbstractMatrix{T}) where T
	#=Threads.@threads=# for g ∈ 1:o.N⋂
		PXY✻ = [o.PXZ[g,i]]
		o.Repl.Yendog[i+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX,g,:,i)'o.v)

		if iszero(j)
			dest[g,:]   = dropdims(colsum(PXY✻ .* (o.DGP.y₁[g] .- view(o.S✻UMZperp,g,:,1)'o.v)); dims=1)
		elseif o.Repl.Yendog[j+1]
			dest[g,:] .-= dropdims(colsum(PXY✻ .* (o.Repl.Z[g,j] * _β̈ .- view(o.S✻UMZperp,g,:,j+1)'β̈v)); dims=1)
		else
			dest[g,:] .-= dropdims(colsum(PXY✻ .* (o.Repl.Z[g,j] * _β̈)); dims=1)
		end
	end
	nothing
end
function FillingLoop2!(o::StrBootTest{T}, dest::Matrix{T}, i::Integer, j::Integer, _β̈::AbstractMatrix{T}, β̈v::AbstractMatrix{T}) where T
	#=Threads.@threads=# for g ∈ 1:o.N⋂
		S = o.info⋂[g]
		PXY✻ = o.Repl.Yendog[i+1] ? view(o.PXZ,S,i) .+ view(o.S✻UPX[:,:,i],S,:) * o.v :
												         reshape(view(o.PXZ,S,i), :, 1)

		if iszero(j)
			dest[g,:]   = dropdims(colsum(PXY✻ .* (o.DGP.y₁[S] .- view(o.S✻UMZperp,S,:,1) * o.v)); dims=1)
		else
			dest[g,:] .-= dropdims(colsum(PXY✻ .* (o.Repl.Yendog[j+1] ? o.Repl.Z[S,j] * _β̈ .- view(o.S✻UMZperp,S,:,j+1) * β̈v :
																														       o.Repl.Z[S,j] * _β̈                                       )); dims=1)
		end
	end
	nothing
end

# Workhorse for WRE CRVE sandwich filling
# Given a zero-based column index, i>0, and a matrix β̈s of all the bootstrap estimates, 
# return all bootstrap realizations of P_X * Z[:,i]_g ' û₁g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per ⋂ cluster, one col per bootstrap replication
# that is, given i, β̈s = δ ̂_CRκ^(*), return, over all g, b (P_(X_∥ g) Z_(∥i)^(*b) )^' (M_(Z_⊥ ) y_(1∥)^(*b) )_g-(P_(X_∥ g) Z_(∥i)^(*b) )^' (M_(Z_⊥ ) Z_∥^(*b) )_g δ ̂_CRκ^(*b)
function Filling!(o::StrBootTest{T}, dest::AbstractMatrix{T}, i::Int64, β̈s::AbstractMatrix, _w::Integer) where T
	if o.granular
		β̈v = _β̈ = Matrix{T}(undef,0,0)
   	if o.Nw == 1  # create or avoid NxB matrix?
			PXY✻ = reshape(view(o.PXZ,:,i), :, 1)
			o.Repl.Yendog[i+1] && (PXY✻ = PXY✻ .+ view(o.S✻UPX,:,:,i) * o.v)
			dest .= @panelsum(PXY✻ .* (o.DGP.y₁ .- view(o.S✻UMZperp,:,:,1) * o.v), o.info⋂)
			@inbounds for j ∈ 1:o.Repl.kZ
				_β̈ = view(β̈s,j,:)'
				dest .-= @panelsum(PXY✻ .* (o.Repl.Yendog[j+1] ? view(o.Repl.Z,:,j) * _β̈ .- view(o.S✻UMZperp,:,:,j+1) * (o.v .* _β̈) :
															                           view(o.Repl.Z,:,j) * _β̈                                    ), o.info⋂)
			end
		else  # create pieces of each N x B matrix one at a time
			@inbounds for j ∈ 0:o.Repl.kZ
				j>0 && (β̈v = o.v .* (_β̈ = view(β̈s,j,:)'))
				(o.purerobust ? FillingLoop1! : FillingLoop2!)(o, dest, i, j, _β̈, β̈v)
			end
    end
  else  # coarse error clustering
		# (P_(X_∥ g) Z_∥^* )^' (M_(Z_⊥ ) y_(1∥)^* )_g
		F1₀ = view(o.Repl.V,:,i)
		F1₁ = view(o.invXXS✻XU₂par,:,:,i)  # zero when FEboot??
		F2₀ = o.S⋂Xy₁
		F2₁ = o.negS✻UMZperpX[1]
		if o.Repl.Yendog[i+1]  # add terms that are zero only if Zpar[i] is exogenous, i.e. if a null refers only to exogenous variables
			dest .= vec(F1₀'F2₀) .- dropdims(F1₀'F2₁ - F2₀'F1₁; dims=1) * o.v  # 0th- & 1st-order terms
			o.Q .= F1₁'F2₁
			@inbounds for g ∈ 1:o.N⋂
				colquadformminus!(dest, g, o.v, o.Q[:,g,:], o.v)
			end
		else
			dest .= vec(F2₀'F1₀) .- dropdims(F1₀'F2₁; dims=1) * o.v  # 0th- & 1st-order terms
		end

		# -(P_(X_∥ g) Z_∥^* )^' (M_(Z_⊥ ) Z_∥^(*b) )_g δ ̂_CRκ^*
		@inbounds for j ∈ 1:o.Repl.kZ
			F2₀ = view(o.S⋂ReplZX,j,:,:)'
			β̈v = o.β̈v[_w]; β̈v .= o.v .* (_β̈ = -view(β̈s,j,:)')
			if o.Repl.Yendog[j+1]
				F2₁ = o.negS✻UMZperpX[j+1]
				dest .+= F2₀'F1₀ .* _β̈ .- (dropdims(F1₀'F2₁; dims=1) - F2₀'F1₁) * β̈v  # "-" because S✻UMZperpX is stored negated as F2₁=negS✻UMZperpX[j+1]
				o.Q .= F1₁'F2₁
				for g ∈ 1:o.N⋂
					colquadformminus!(dest, g, o.v, o.Q[:,g,:], β̈v)
				end
			elseif o.Repl.Yendog[i+1]
				dest .+= F2₀'F1₀ .* _β̈  .+ F2₀'F1₁ * β̈v
			else
				dest .+= F2₀'F1₀ .* _β̈
			end
		end
  end
  nothing
end

function MakeWREStats!(o::StrBootTest{T}, w::Integer) where T
	_w = isone(o.Nw) || w<o.Nw ? 1 : 2

  if isone(o.Repl.kZ)  # optimized code for 1 retained coefficient in bootstrap regression
		_As = view(o.As[_w],1,:,:)
		if o.liml
			HessianFixedkappa!(o, o.YY₁₁[_w]  , [0], 0, zero(T), _w)  # κ=0 => Y*MZperp*Y
			HessianFixedkappa!(o, o.YY₁₂[_w]  , [0], 1, zero(T), _w)
			HessianFixedkappa!(o, o.YY₂₂[_w]  , [1], 1, zero(T), _w)
			HessianFixedkappa!(o, o.YPXY₁₁[_w], [0], 0, one(T) , _w)  # κ=1 => Y*PXpar*Y
			HessianFixedkappa!(o, o.YPXY₁₂[_w], [0], 1, one(T) , _w)
			HessianFixedkappa!(o, o.YPXY₂₂[_w], [1], 1, one(T) , _w)
			o.YY₁₂YPXY₁₂[_w] .= o.YY₁₂[_w] .* o.YPXY₁₂[_w]
			o.x₁₁[_w] .= o.YY₂₂[_w] .* o.YPXY₁₁[_w] .- o.YY₁₂YPXY₁₂[_w]      # elements of YY✻^-1 * YPXY✻ up to factor of det(YY✻)
			o.x₁₂[_w] .= o.YY₂₂[_w] .* o.YPXY₁₂[_w] .- o.YY₁₂[_w] .* o.YPXY₂₂[_w]
			o.x₂₁[_w] .= o.YY₁₁[_w] .* o.YPXY₁₂[_w] .- o.YY₁₂[_w] .* o.YPXY₁₁[_w]
			o.x₂₂[_w] .= o.YY₁₁[_w] .* o.YPXY₂₂[_w] .- o.YY₁₂YPXY₁₂[_w]
			o.κs[_w] .= (o.x₁₁[_w] .+ o.x₂₂[_w])./2
			o.κs[_w] .= 1 ./ (1 .- (o.κs[_w] .- sqrtNaN.(o.κs[_w].^2 .- o.x₁₁[_w] .* o.x₂₂[_w] .+ o.x₁₂[_w] .* o.x₂₁[_w])) ./ 
			                 (o.YY₁₁[_w] .* o.YY₂₂[_w] .- o.YY₁₂[_w] .* o.YY₁₂[_w]))  # solve quadratic equation for smaller eignenvalue; last term is det(YY✻)
			!iszero(o.fuller) && (o.κs[_w] .-= o.fuller / (o._Nobs - o.kX))
			_As .= o.κs[_w] .* (o.YPXY₂₂[_w] .- o.YY₂₂[_w]) .+ o.YY₂₂[_w]
			o.β̈s[_w] .= (o.κs[_w] .* (o.YPXY₁₂[_w] .- o.YY₁₂[_w]) .+ o.YY₁₂[_w]) ./ _As
		else
			HessianFixedkappa!(o, _As, [1], 1, o.κ, _w)
			HessianFixedkappa!(o, o.β̈s[_w], [1], 0, o.κ, _w); o.β̈s[_w] ./= _As
		end

		if o.null
			o.numerWRE[_w] .= o.β̈s[_w] .+ (o.Repl.Rt₁ - o.r) / o.Repl.RRpar
		else
			o.numerWRE[_w] .= o.β̈s[_w] .- view(o.DGP.β̈ ,:,1)
			isone(w) && (o.numerWRE[_w][1] = o.β̈s[_w][1] + (o.Repl.Rt₁[1] - o.r[1]) / o.Repl.RRpar[1])
		end

		@storeWtGrpResults!(o.numer, o.numerWRE[_w])
		if o.bootstrapt
			if o.robust
				J⋂s1 = dropdims(o.J⋂s[_w]; dims=3)
				Filling!(o, J⋂s1, 1, o.β̈s[_w], _w); 
				J⋂s1 ./= _As
				@inbounds for c ∈ 1:o.NErrClustCombs  # sum sandwich over error clusterings
					nrows(o.clust[c].order)>0 && 
						(J⋂s1 .= J⋂s1[o.clust[c].order,:])
					@clustAccum!(denom, c, coldot(@panelsum(J⋂s1, o.clust[c].info)))  # XXX allocation in here where c=1
				end
			else
				denom = (HessianFixedkappa(o, [0], 0, zero(T), _w) .-   # XXX rewrite to avoid allocations
				         2 .* o.β̈s[_w] .* HessianFixedkappa(o, [0], 1, zero(T), _w) .+ 
								 o.β̈s[_w].^2 .* HessianFixedkappa(o, [1], 1, zero(T), _w)) ./ o._Nobs ./ _As  # classical error variance
			end
			@storeWtGrpResults!(o.dist, o.sqrt ? o.numerWRE[_w] ./ sqrtNaN.(denom) : o.numerWRE[_w] .^ 2 ./ denom)
			denom *= o.Repl.RRpar[1]^2
		end
		w==1 && o.bootstrapt && (o.statDenom = fill(denom[1],1,1))  # original-sample denominator

  else  # WRE bootstrap for more than 1 retained coefficient in bootstrap regression

		if o.liml
			for i ∈ 0:o.Repl.kZ
				HessianFixedkappa!(o, o.YY✻[_w][i+1]  , collect(0:i), i, zero(T), _w)  # κ=0 => Y*MZperp*Y
				HessianFixedkappa!(o, o.YPXY✻[_w][i+1], collect(0:i), i,  one(T), _w)  # κ=1 => Y*PXpar*Y
			end

			@inbounds for b ∈ axes(o.v,2)
				for i ∈ 0:o.Repl.kZ
					o.YY✻_b[1:i+1,i+1]   = o.YY✻[_w][i+1][:,b]  # fill uppper triangles, which is all that invsym() looks at
					o.YPXY✻_b[1:i+1,i+1] = o.YPXY✻[_w][i+1][:,b]
				end
				o.κ = 1/(1 - real(eigvalsNaN(invsym(o.YY✻_b) * Symmetric(o.YPXY✻_b))[1]))
				!iszero(o.fuller) && (o.κ -= o.fuller / (o._Nobs - o.kX))
				o.β̈s[_w][:,b] = (o.As[_w][:,:,b] = invsym(o.κ*o.YPXY✻_b[2:end,2:end] + (1-o.κ)*o.YY✻_b[2:end,2:end])) * (o.κ*o.YPXY✻_b[1,2:end] + (1-o.κ)*o.YY✻_b[1,2:end])
			end
		else
			HessianFixedkappa!(o, o.δnumer[_w], collect(1:o.Repl.kZ), 0, o.κ, _w)
			for i ∈ 1:o.Repl.kZ
				HessianFixedkappa!(o, o.δdenom[_w][i], collect(1:i), i, o.κ, _w)
			end

			nt = Threads.nthreads()
  		cs = [round(Int, size(o.v,2)/nt*i) for i ∈ 0:nt]
			@inbounds #=Threads.@threads=# for t ∈ 1:nt
				δdenom_b = zeros(T, o.Repl.kZ, o.Repl.kZ)  # thread-safe scratch pad
				@inbounds for b ∈ cs[t]+1:cs[t+1]
					for i ∈ 1:o.Repl.kZ
						δdenom_b[1:i,i] = view(o.δdenom[_w][i],:,b)  # fill uppper triangle
					end
					o.β̈s[_w][:,b] = (o.As[_w][:,:,b] = invsym(δdenom_b)) * view(o.δnumer[_w],:,b)  # XXX move b to middle index and use 3-array operators?
				end
			end
		end

		if o.bootstrapt
			if o.robust
				J⋂s = o.J⋂s[_w]
				@inbounds for i ∈ 1:o.Repl.kZ  # avoid list comprehension construction for compiler-perceived type stability
					Filling!(o, view(J⋂s,:,:,i), i, o.β̈s[_w], _w)
				end
			else
				YY✻ = [HessianFixedkappa(o, collect(i:o.Repl.kZ), i, zero(T), _w) for i ∈ 0:o.Repl.kZ]  # κ=0 => Y*MZperp*Y
			end
		end

		let denom_b
			@inbounds for b ∈ reverse(axes(o.v,2))
				o.numer_b .= o.null || w==1 && b==1 ? o.Repl.RRpar * view(o.β̈s[_w],:,b) + o.Repl.Rt₁ - o.r : o.Repl.RRpar * (view(o.β̈s[_w],:,b) - view(o.DGP.β̈ ,:,1))
				if o.bootstrapt
					if o.robust  # Compute denominator for this WRE test stat
						J⋂ = view(J⋂s,:,b,:) * (view(o.As[_w],:,:,b) * o.Repl.RRpar')
						for c ∈ 1:o.NErrClustCombs
							(!isone(o.NClustVar) && nrows(o.clust[c].order)>0) &&
								(J⋂ = J⋂[o.clust[c].order,:])
							J_b = @panelsum(J⋂, o.clust[c].info)
							@clustAccum!(denom_b, c, J_b'J_b)
						end
					else  # non-robust
						for i ∈ 0:o.Repl.kZ
							o.YY✻_b[i+1,i+1:o.Repl.kZ+1] = view(YY✻[i+1],:,b)  # fill upper triangle
						end
						tmp = [-one(T) ; o.β̈s[_w][:,b]]
						denom_b = (o.Repl.RRpar * view(o.As[_w],:,:,b) * o.Repl.RRpar') * tmp'Symmetric(o.YY✻_b) * tmp / o._Nobs  # 2nd half is sig2 of errors
					end
					o.dist[b+first(o.WeightGrp[w])-1] = o.sqrt ? o.numer_b[1] / sqrtNaN(denom_b[1]) : o.numer_b'invsym(denom_b)*o.numer_b  # hand-code for 2-dimensional?
				end
				o.numer[:,b+first(o.WeightGrp[w])-1] = o.numer_b  # slight inefficiency: in usual bootstrap-t case, only need to save numerators in numer if getdist("numer") is coming because of svmat(numer)
			end
			w==1 && o.bootstrapt && (o.statDenom = denom_b)  # original-sample denominator
		end
	end
	nothing
end
