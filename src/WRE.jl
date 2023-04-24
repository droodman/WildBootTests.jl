# stuff done once per exucution--not depending on r
function InitWRE!(o::StrBootTest{T}) where T
	if o.null
		if o.granular || o.jk
			o.Ü₂par = Matrix{T}(undef, o.Nobs, o.Repl.kZ)
			o.Z̄ = similar(o.Ü₂par)
		end
	else  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
		EstimateIV!(o.DGP, o, o.jk, o.r₁)
		MakeResidualsIV!(o.DGP, o)
		if o.granular || o.jk
			o.Ü₂par = view(o.DGP.Ü₂ * o.Repl.RparY,:,:)
			o.Z̄ .= t✻(o.DGP.Ȳ₂, o.Repl.RparY); o.Z̄ .+= o.Repl.ZparX
		end
	end

	o.Repl.kZ>1 && (o.numer_b = Vector{T}(undef,nrows(o.Repl.RRpar)))

	if o.willfill
		if o.Repl.kZ==1 
			o.J⋂s1 = Vector{Matrix{T}}(undef, o.NErrClustCombs)
			for c ∈ 1:o.NErrClustCombs
				o.J⋂s1[c] = c==1 || nrows(o.clust[c].order)>0 ? Matrix{T}(undef, o.N⋂, o.ncolsv) : o.J⋂s1[c-1]
			end
		else
			o.J⋂s = Array{T,3}(undef, o.N⋂, o.ncolsv, o.Repl.kZ)
			o.J⋂ARpars = Vector{Matrix{T}}(undef, o.NErrClustCombs)
			for c ∈ 1:o.NErrClustCombs
				o.J⋂ARpars[c] = c==1 || nrows(o.clust[c].order)>0 ? Array{T,3}(undef, o.N⋂, o.ncolsv, o.q) : o.J⋂ARpars[c-1]
			end
			o.ARpars = Array{T,3}(undef, o.Repl.kZ, o.ncolsv, o.q)
		end
		o.Jc = [c==1 ?  Array{T,3}(undef,0,0,0) : Array{T,3}(undef, o.clust[c].N, o.ncolsv, o.q) for c ∈ 1:o.NErrClustCombs]

		o.jk && (o.uⱼₖ = Vector{T}(undef,o.Nobs))
	end

	o.T1L = Matrix{T}(undef, o.DGP.kX, o.ncolsv)
	o.T1R = similar(o.T1L)
	if o.Repl.kZ==1
		o.β̈sAs = Matrix{T}(undef, 2, o.ncolsv)
	else
		o.β̈s = Matrix{T}(undef, o.Repl.kZ, o.ncolsv)
		o.As = Array{T,3}(undef, o.Repl.kZ, o.ncolsv, o.Repl.kZ)
	end
	o.numerWRE = Matrix{T}(undef, o.dof, o.ncolsv)
	o.invXXXZ̄ = Matrix{T}(undef, o.DGP.kX, o.Repl.kZ)
	o.XȲ = Matrix{T}(undef, o.DGP.kX, o.Repl.kZ+1) 
	o.ZÜ₂par = Matrix{T}(undef, o.Repl.kZ, o.Repl.kZ)
	o.ȲȲ = Matrix{T}(undef, o.Repl.kZ+1, o.Repl.kZ+1)

	if isone(o.Repl.kZ)
		if o.liml
			o.YY₁₁ = Matrix{T}(undef, o.Repl.kZ, o.ncolsv)
			o.YY₁₂ = similar(o.YY₁₁)
			o.YY₂₂   = similar(o.YY₁₁)
			o.YPXY₁₁ = similar(o.YY₁₁)
			o.YPXY₁₂ = similar(o.YY₁₁)
			o.YPXY₂₂ = similar(o.YY₁₁)
			o.YY₁₂YPXY₁₂ = similar(o.YY₁₁)
			o.x₁₁ = similar(o.YY₁₁)
			o.x₁₂ = similar(o.YY₁₁)
			o.x₂₁ = similar(o.YY₁₁)
			o.x₂₂ = similar(o.YY₁₁)
			o.κs = similar(o.YY₁₁)
		end
	else
		o.denomWRE = Array{T,3}(undef, o.q, o.ncolsv, o.q)
		if o.liml
			o.YY✻ = Array{T,3}(undef, o.Repl.kZ+1, o.ncolsv, o.Repl.kZ+1)
			o.YPXY✻ = similar(o.YY✻)
			o.κWRE = Array{T,3}(undef,1,o.ncolsv,1)
		else
			o.δnumer = Matrix{T}(undef, o.Repl.kZ, o.ncolsv)
		end
		if o.bootstrapt && !o.robust
			o.YY✻ = Array{T,3}(undef, o.Repl.kZ+1, o.ncolsv, o.Repl.kZ+1)
		end
	end

	if o.bootstrapt
		if o.NFE>0 && !o.FEboot && (o.not2SLS || o.willfill)
			t = crosstabFE(o, o.y₁, o.ID✻, o.N✻)[]
			o.CT✻FEU₂       = fill(t, o.kY₂      )  # placeholder vectors of sparse matrices; sparsity pattern will never change
			o.CT✻FEU        = fill(t, o.Repl.kZ+1)
			o.invFEwtCT✻FEU = fill(t, o.Repl.kZ+1)
		end

		o.willfill &&
			(o.β̈v = Matrix{T}(undef, o.N✻, o.ncolsv))
	end

	o.S✻Xu₁         = Array{T,3}(undef, o.DGP.kX, o.N✻, 1)
	o.S✻XU₂         = Array{T,3}(undef, o.DGP.kX, o.N✻, o.kY₂)
	o.S✻XU₂par      = Array{T,3}(undef, o.DGP.kX, o.N✻, o.Repl.kZ) 
	o.invXXS✻Xu₁    = Array{T,3}(undef, o.DGP.kX, o.N✻, 1)
	o.invXXS✻XU₂par = Array{T,3}(undef, o.DGP.kX, o.N✻, o.Repl.kZ)
	o.S✻XU      = [i>0 ? view(o.S✻XU₂par     ,:,:,i) : view(o.S✻Xu₁     ,:,:,1) for i ∈ 0:o.Repl.kZ]
	o.invXXS✻XU = [i>0 ? view(o.invXXS✻XU₂par,:,:,i) : view(o.invXXS✻Xu₁,:,:,1) for i ∈ 0:o.Repl.kZ]

	if o.willfill || o.not2SLS
		o.S✻Zperpu₁                 = Array{T,3}(undef, o.DGP.kZperp, o.N✻, 1)
		o.S✻ZperpU₂par              = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.Repl.kZ)
		o.invZperpZperpS✻Zperpu₁    = Array{T,3}(undef, o.DGP.kZperp, o.N✻, 1)
		o.invZperpZperpS✻ZperpU₂par = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.Repl.kZ)
		o.S✻ZperpU              = [i>0 ? view(o.S✻ZperpU₂par             ,:,:,i) : view(o.S✻Zperpu₁             ,:,:,1) for i ∈ 0:o.Repl.kZ]
		o.invZperpZperpS✻ZperpU = [i>0 ? view(o.invZperpZperpS✻ZperpU₂par,:,:,i) : view(o.invZperpZperpS✻Zperpu₁,:,:,1) for i ∈ 0:o.Repl.kZ]

		o.S✻ZperpUv              = Matrix{T}(undef, o.DGP.kZperp, o.ncolsv)
		o.invZperpZperpS✻ZperpUv = Matrix{T}(undef, o.DGP.kZperp, o.ncolsv)

		if o.NFE>0 & !o.FEboot
			o.CT✻FEUv = Matrix{T}(undef, o.NFE, o.ncolsv)
			o.invFEwtCT✻FEUv = Matrix{T}(undef, o.NFE, o.ncolsv)
		end

		o.S✻u₁u₁       = Array{T,3}(undef, 1        , o.N✻, 1        )
		o.S✻U₂paru₁    = Array{T,3}(undef, o.Repl.kZ, o.N✻, 1        )
		o.S✻U₂parU₂par = Array{T,3}(undef, o.Repl.kZ, o.N✻, o.Repl.kZ)
		o.S✻UU = [i>0 ? j>0 ? view(o.S✻U₂parU₂par,i,:,j) : view(o.S✻U₂paru₁,i,:,1) : j>0 ? view(o.S✻U₂paru₁,j,:,1) : view(o.S✻u₁u₁,1,:,1) for i ∈ 0:o.Repl.kZ, j ∈ 0:o.Repl.kZ]

		o.S✻ȲUfold = Array{T,3}(undef, o.Repl.kZ+1, o.N✻, o.Repl.kZ+1)
	end

	if o.not2SLS
		o.S✻ȳ₁u₁   = Array{T,3}(undef, 1, o.N✻, 1)
		o.S✻Z̄u₁    = Array{T,3}(undef, o.Repl.kZ, o.N✻, 1)
		o.S✻ȳ₁Ü₂par = Array{T,3}(undef, 1, o.N✻, o.Repl.kZ)
		o.S✻Z̄Ü₂par = Array{T,3}(undef, o.Repl.kZ, o.N✻, o.Repl.kZ)
		o.S✻ȲU     = [i>0 ? j>0 ? view(o.S✻Z̄Ü₂par,i,:,j) : view(o.S✻Z̄u₁,i,:,1) : j>0 ? view(o.S✻ȳ₁Ü₂par,1,:,j) : view(o.S✻ȳ₁u₁,1,:,1) for i ∈ 0:o.Repl.kZ, j ∈ 0:o.Repl.kZ]
	end

	if o.granular
		if o.willfill
			o.XinvXX = X₁₂B(o.Repl.X₁, o.Repl.X₂, o.Repl.invXX)
			o.PXZ̄ = Matrix{T}(undef, o.Nobs, o.Repl.kZ)
			o.S✻XUv = Matrix{T}(undef, o.DGP.kX, o.ncolsv)
			o.PXY✻ = Matrix{T}(undef, o.purerobust ? 1 : mapreduce(length, max, o.info⋂), o.ncolsv)
			o.S✻UMZperp = similar(o.PXY✻)
			o.S✻UZperpinvZperpZperpv = Matrix{T}(undef, o.DGP.kZperp, o.ncolsv)
		end
	else
		o.Π̈Rpar = Matrix{T}(undef, o.DGP.kX, o.Repl.kZ)

		if o.willfill || !o.jk
			S✻⋂ZperpX = o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpX; S✻⋂ZperpX .= o.DGP.S✻⋂XZperp' .- S✻⋂ZperpX
			o.S✻⋂XX   = o.DGP.S✻⋂XZperp * o.DGP.invZperpZperpZperpX; o.S✻⋂XX .= o.DGP.S✻⋂XX .- o.S✻⋂XX; t✻minus!(o.S✻⋂XX, o.DGP.invZperpZperpZperpX', S✻⋂ZperpX)
		end

		if o.willfill
			o.S⋂XX = @panelsum(o.S✻⋂XX, o.info⋂_✻⋂)
			S⋂ZperpX   = @panelsum(S✻⋂ZperpX, o.info⋂_✻⋂)
			o.S⋂ȳ₁X 	  = Array{T,3}(undef, 1, o.N⋂, o.DGP.kX)
			o.S⋂ReplZ̄X = Array{T,3}(undef, o.Repl.kZ, o.N⋂, o.DGP.kX)

			o.S⋂XZperpinvZperpZperp = S⋂ZperpX' * o.DGP.invZperpZperp
			o.negS✻UMZperpX = [Array{T,3}(undef, o.DGP.kX, o.N⋂, o.N✻) for _ in 0:o.Repl.kZ]

			_inds = o.subcluster>0 ?
							[CartesianIndex(j,i) for (j,v) ∈ enumerate(o.info⋂_✻⋂) for i ∈ v] :  # crosstab ∩,* is wide
							o.NClustVar == o.NBootClustVar ?
									[CartesianIndex(i,i) for i ∈ 1:o.N✻⋂] :  # crosstab *,∩ is square
									[CartesianIndex(i,j) for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab ∩,* is tall
			inds = [CartesianIndex(k,I) for I ∈ _inds for k ∈ 1:o.DGP.kX]
			o.crosstab⋂✻ind = LinearIndices(o.negS✻UMZperpX[1])[inds]

			o.F₁ = Matrix{T}(undef, o.DGP.kX, o.ncolsv)
			o.F₁β = similar(o.F₁)
			o.F₂ = similar(o.F₁)
	
			if o.NFE>0 && !o.FEboot
				o.CT⋂FEX = crosstabFE(o, o.X₁, o.X₂, o.ID⋂, o.N⋂)
				broadcast!(.*, o.CT⋂FEX, Ref(o.invsumFEwt), o.CT⋂FEX)
			end
			o.S✻⋂Xu₁ = Array{T,3}(undef, o.DGP.kX, o.N✻⋂, 1)
			o.S✻⋂XÜ₂par = Array{T,3}(undef, o.DGP.kX, o.N✻⋂, o.Repl.kZ)
		end

		if !o.jk
			o.invXXS✻XU₂ = Array{T,3}(undef, o.DGP.kX, o.N✻, o.kY₂)
			if o.willfill || o.not2SLS  
				o.S✻ZperpU₂ = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.kY₂)
				o.invZperpZperpS✻ZperpU₂ = Array{T,3}(undef, o.DGP.kZperp, o.N✻, o.kY₂)
			end
	
			o.S✻⋂XY₂        = o.DGP.S✻⋂XY₂     - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpY₂   - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpY₂   - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpY₂ )
			o.S✻⋂XDGPZ      = o.DGP.S✻⋂XZpar   - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpZ - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZpar - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZ)
			o.S✻⋂Xy₁        = o.DGP.S✻⋂Xy₁     - o.DGP.S✻⋂XZperp     * o.DGP.invZperpZperpZperpy₁   - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂Zperpy₁   - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpy₁ )
			_S✻ZperpY₂      = @panelsum(o.DGP.S✻⋂ZperpY₂, o.info✻_✻⋂)  # moments of variables _before_ FWL processing
			_S✻Zperpy₁      = @panelsum(o.DGP.S✻⋂Zperpy₁, o.info✻_✻⋂)
			_S✻ZperpDGPZpar = @panelsum(o.DGP.S✻⋂ZperpZpar, o.info✻_✻⋂)

			S✻ZperpZperp  = @panelsum(o.DGP.S✻⋂ZperpZperp, o.info✻_✻⋂)
			o.S✻XY₂       = @panelsum(o.S✻⋂XY₂  , o.info✻_✻⋂)
			o.S✻XX        = @panelsum(o.S✻⋂XX   , o.info✻_✻⋂)
			o.S✻XDGPZ     = @panelsum(o.S✻⋂XDGPZ, o.info✻_✻⋂)
			o.S✻Xy₁       = @panelsum(o.S✻⋂Xy₁  , o.info✻_✻⋂)
			o.S✻ZperpX    = @panelsum(S✻⋂ZperpX , o.info✻_✻⋂)
			o.S✻ZperpY₂   = S✻ZperpZperp * o.DGP.invZperpZperpZperpY₂; o.S✻ZperpY₂   .= _S✻ZperpY₂      .- o.S✻ZperpY₂
			o.S✻ZperpDGPZ = S✻ZperpZperp * o.DGP.invZperpZperpZperpZ ; o.S✻ZperpDGPZ .= _S✻ZperpDGPZpar .- o.S✻ZperpDGPZ
			o.S✻Zperpy₁   = S✻ZperpZperp * o.DGP.invZperpZperpZperpy₁; o.S✻Zperpy₁   .= _S✻Zperpy₁      .- o.S✻Zperpy₁

			if o.NFE>0 && !o.FEboot && (o.willfill || o.not2SLS)
				o.CT✻FEX   = crosstabFE(o, o.X₁, o.X₂, o.ID✻, o.N✻)
				o.CT✻FEY₂  = crosstabFE(o, o.DGP.Y₂, o.ID✻, o.N✻)
				o.CT✻FEZ   = crosstabFE(o, o.DGP.Zpar, o.ID✻, o.N✻)
				o.CT✻FEy₁  = crosstabFE(o, o.DGP.y₁, o.ID✻, o.N✻)
				o.DGP.restricted &&
					(o.CT✻FEZR₁ = crosstabFE(o, o.DGP.ZR₁, o.ID✻, o.N✻))  #  XXX just do o.CT✻FEZ * R₁ ?
			end

			o.willfill &&
				(o.S✻⋂XU₂ = Array{T,3}(undef, o.DGP.kX, o.N✻⋂, o.kY₂))

			if o.DGP.restricted
				o.S✻⋂X_DGPZR₁ = o.DGP.S✻⋂XZR₁ - o.DGP.S✻⋂XZperp * o.DGP.invZperpZperpZperpZR₁ - o.DGP.invZperpZperpZperpX' * (o.DGP.S✻⋂ZperpZR₁  - o.DGP.S✻⋂ZperpZperp * o.DGP.invZperpZperpZperpZR₁ )
				_S✻ZperpDGPZR₁ = @panelsum(o.DGP.S✻⋂ZperpZR₁, o.info✻_✻⋂)
				o.S✻XZR₁ = @panelsum(o.S✻⋂X_DGPZR₁, o.info✻_✻⋂)
				o.S✻ZperpDGPZR₁ = @panelsum(o.DGP.S✻⋂ZperpZR₁, o.info✻_✻⋂) - S✻ZperpZperp * o.DGP.invZperpZperpZperpZR₁
				# o.S✻ZperpDGPZR₁ = @panelsum(o.DGP.S✻⋂ZperpZR₁, o.info✻_✻⋂); t✻minus!(o.S✻ZperpDGPZR₁, S✻ZperpZperp, o.DGP.invZperpZperpZperpZR₁)
			end

			if o.not2SLS  # cluster-wise moments after FWL
				o.S✻Y₂Y₂     = o.DGP.S✻Y₂Y₂     - _S✻ZperpY₂'      * o.DGP.invZperpZperpZperpY₂ - o.DGP.invZperpZperpZperpY₂' * o.S✻ZperpY₂
				o.S✻DGPZDGPZ = o.DGP.S✻ZparZpar - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpZ  - o.DGP.invZperpZperpZperpZ'  * o.S✻ZperpDGPZ
				o.S✻DGPZY₂   = o.DGP.S✻ZparY₂   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpY₂ - o.DGP.invZperpZperpZperpZ'  * o.S✻ZperpY₂
				o.S✻DGPZy₁   = o.DGP.S✻Zpary₁   - _S✻ZperpDGPZpar' * o.DGP.invZperpZperpZperpy₁ - o.DGP.invZperpZperpZperpZ'  * o.S✻Zperpy₁   
				o.S✻Y₂y₁     = o.DGP.S✻Y₂y₁     - _S✻ZperpY₂'      * o.DGP.invZperpZperpZperpy₁ - o.DGP.invZperpZperpZperpY₂' * o.S✻Zperpy₁
				o.S✻y₁y₁     = o.DGP.S✻y₁y₁     - _S✻Zperpy₁'      * o.DGP.invZperpZperpZperpy₁ - o.S✻Zperpy₁'               * o.DGP.invZperpZperpZperpy₁
				o.DGP.restricted && 
					(o.S✻DGPZR₁y₁ = o.DGP.S✻ZR₁y₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpy₁ - o.DGP.invZperpZperpZperpZR₁' * o.S✻Zperpy₁)

				if o.DGP.restricted
					o.S✻DGPZR₁Y₂     = o.DGP.S✻ZR₁Y₂  - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpY₂  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpY₂
					o.S✻DGPZR₁DGPZR₁ = o.DGP.S✻ZR₁ZR₁ - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZR₁  - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZR₁
					o.S✻DGPZR₁DGPZ   = o.DGP.S✻ZR₁Z   - _S✻ZperpDGPZR₁' * o.DGP.invZperpZperpZperpZ - o.DGP.invZperpZperpZperpZR₁' * o.S✻ZperpDGPZ
					o.S✻DGPZR₁X      = (@panelsum(o.DGP.S✻⋂XZR₁, o.info✻_✻⋂))'; t✻minus!(o.S✻DGPZR₁X, _S✻ZperpDGPZR₁', o.DGP.invZperpZperpZperpX); t✻minus!(o.S✻DGPZR₁X, o.DGP.invZperpZperpZperpZR₁', o.S✻ZperpX)
				end
			end

			o.invXXS✻XDGPZ = @panelsum(o.DGP.invXX * o.S✻⋂XDGPZ, o.info✻_✻⋂)
			o.invXXS✻Xy₁   = o.DGP.invXX * o.S✻Xy₁
			o.invZperpZperpS✻ZperpY₂   = o.DGP.invZperpZperp * o.S✻ZperpY₂ 
			o.invZperpZperpS✻ZperpX    = o.DGP.invZperpZperp * o.S✻ZperpX  
			o.invZperpZperpS✻Zperpy₁   = o.DGP.invZperpZperp * o.S✻Zperpy₁ 
			o.invZperpZperpS✻ZperpDGPZ = o.DGP.invZperpZperp * o.S✻ZperpDGPZ

			if o.DGP.restricted
				o.invXXS✻XDGPZR₁ = o.DGP.invXX * o.S✻XZR₁ 
				o.invZperpZperpS✻ZperpDGPZR₁ = o.DGP.invZperpZperp * o.S✻ZperpDGPZR₁
			end
		end
	end
	nothing
end

# stuff done when r changes, but outside looping in HessianFixedkappa() and Filling!()
function PrepWRE!(o::StrBootTest{T}) where T
	if o.null
		r₁ = [o.r₁ ; o.r]
	  EstimateIV!(o.DGP, o, o.jk, r₁)
	  MakeResidualsIV!(o.DGP, o)
	else
		r₁ = o.r₁
	end

	o.invXXXZ̄ .= o.Repl.XZ - o.DGP.XÜ₂ * o.Repl.RparY
	o.XȲ .= [o.DGP.Xȳ₁ o.invXXXZ̄]
  o.invXXXZ̄ .= o.Repl.invXX * o.invXXXZ̄
  o.ZÜ₂par .= (o.Repl.ZY₂	 - o.Repl.XZ'o.DGP.Π̈ ) * o.Repl.RparY
  _ȲȲ = o.DGP.γ⃛'o.Repl.XZ - o.DGP.ȳ₁Ü₂ * o.Repl.RparY 
  o.ȲȲ .= [o.DGP.ȳ₁ȳ₁ _ȲȲ
           _ȲȲ'       o.Repl.ZZ - o.ZÜ₂par' - o.ZÜ₂par + o.Repl.RparY'o.DGP.Ü₂Ü₂*o.Repl.RparY]

	if o.granular || o.jk
		if o.null
	  	mul!(o.Ü₂par, o.DGP.Ü₂, o.Repl.RparY)
			mul!(o.Z̄, o.DGP.Ȳ₂, o.Repl.RparY); o.Z̄ .+= o.Repl.ZparX
		end

		panelcross21!(o.S✻Xu₁, o.DGP.X₁, o.DGP.X₂, o.DGP.u⃛₁, o.info✻)
		panelcross21!(o.S✻XU₂par, o.DGP.X₁, o.DGP.X₂, o.Ü₂par, o.info✻)
		t✻!(o.invXXS✻Xu₁   , o.DGP.invXX, o.S✻Xu₁   )
		t✻!(o.invXXS✻XU₂par, o.DGP.invXX, o.S✻XU₂par)

		if o.willfill || o.not2SLS
			panelcross!(o.S✻Zperpu₁, o.DGP.Zperp, o.DGP.u⃛₁, o.info✻)
			panelcross!(o.S✻ZperpU₂par, o.DGP.Zperp, o.Ü₂par, o.info✻)
			t✻!(o.invZperpZperpS✻Zperpu₁, o.DGP.invZperpZperp, o.S✻Zperpu₁)
			t✻!(o.invZperpZperpS✻ZperpU₂par, o.DGP.invZperpZperp, o.S✻ZperpU₂par)
			if o.NFE>0 && !o.FEboot
				crosstabFE!(o, (@view o.CT✻FEU[1:1   ]), [o.DGP.u⃛₁], o.ID✻, o.N✻)
				crosstabFE!(o, (@view o.CT✻FEU[2:end]), [o.Ü₂par  ], o.ID✻, o.N✻)
				broadcast!(.*, o.invFEwtCT✻FEU, Ref(o.invsumFEwt), o.CT✻FEU)
			end
		end

		if o.not2SLS
			panelcross!(o.S✻u₁u₁, o.DGP.u⃛₁, o.DGP.u⃛₁, o.info✻)
			panelcross!(o.S✻U₂paru₁, o.Ü₂par, o.DGP.u⃛₁, o.info✻)
			panelcross!(o.S✻U₂parU₂par, o.Ü₂par, o.Ü₂par, o.info✻)

			panelcross!(o.S✻ȳ₁u₁, o.DGP.ȳ₁, o.DGP.u⃛₁, o.info✻)
			panelcross!(o.S✻Z̄u₁, o.Z̄, o.DGP.u⃛₁, o.info✻)
			panelcross!(o.S✻ȳ₁Ü₂par, o.DGP.ȳ₁, o.Ü₂par, o.info✻)
			panelcross!(o.S✻Z̄Ü₂par, o.Z̄, o.Ü₂par, o.info✻)
		end
	end

	o.granular && o.willfill &&
		X₁₂B!(o.PXZ̄, o.Repl.X₁, o.Repl.X₂, o.invXXXZ̄)

	if !o.granular
		if o.willfill || o.not2SLS
			t✻!(o.Π̈Rpar, o.DGP.Π̈ , o.Repl.RparY)
			!iszero(o.DGP.kX₁) && (o.Π̈Rpar[1:o.DGP.kX₁,:] += o.Repl.Xpar₁toZparX)
		end

		if !o.jk  # in coarse case, if not jackknifing, construct things while avoiding O(N) operations
			Π⃛y = [o.DGP.RperpXperp'o.DGP.γ̈X ; zeros(T, o.kX₂, 1)] + o.DGP.Π̈ * o.DGP.γ̈Y

			t✻!(o.S✻XU₂, o.S✻XX, o.DGP.Π̈); o.S✻XU₂ .= o.S✻XY₂ .- o.S✻XU₂
			o.S✻XU₂par .= o.S✻XU₂ * o.Repl.RparY  # use this syntax for 3-array x DesignerMatrix
			t✻!(o.invXXS✻XU₂, o.DGP.invXX, o.S✻XU₂)
			o.invXXS✻XU₂par .= o.invXXS✻XU₂ * o.Repl.RparY  # use this syntax for 3-array x DesignerMatrix
			if o.willfill || o.not2SLS
				t✻!(o.S✻ZperpU₂, o.S✻ZperpX, o.DGP.Π̈); o.S✻ZperpU₂ .= o.S✻ZperpY₂ .- o.S✻ZperpU₂
				o.invZperpZperpS✻ZperpU₂ .= o.invZperpZperpS✻ZperpY₂; t✻minus!(o.invZperpZperpS✻ZperpU₂, o.invZperpZperpS✻ZperpX, o.DGP.Π̈)
				o.S✻ZperpU₂par .= o.S✻ZperpU₂ * o.Repl.RparY  # use this syntax for 3-array x DesignerMatrix
				o.invZperpZperpS✻ZperpU₂par .= o.invZperpZperpS✻ZperpU₂ * o.Repl.RparY
			end

			o.S✻Xu₁ .= o.S✻Xy₁; t✻minus!(o.S✻Xu₁, o.S✻XDGPZ, o.DGP.β̈ ); t✻plus!(o.S✻Xu₁, o.S✻XU₂, o.DGP.γ̈Y )
			o.DGP.restricted &&
				t✻minus!(o.S✻Xu₁, o.S✻XZR₁, r₁)

			o.invXXS✻Xu₁ .= o.invXXS✻Xy₁; t✻minus!(o.invXXS✻Xu₁, o.invXXS✻XDGPZ, o.DGP.β̈ ); t✻plus!(o.invXXS✻Xu₁, o.invXXS✻XU₂, o.DGP.γ̈Y )
			o.DGP.restricted &&
				t✻minus!(o.invXXS✻Xu₁, o.invXXS✻XDGPZR₁, r₁)

			if o.not2SLS
				Π̂S✻XÜ₂γ̈Y = (o.DGP.Π̈ )' * o.S✻XU₂ * o.DGP.γ̈Y
				S✻Ü₂Y₂ = o.S✻Y₂Y₂ - (o.DGP.Π̈ )' * o.S✻XY₂
				S✻Y₂Ü₂γ̈Y = S✻Ü₂Y₂' * o.DGP.γ̈Y

				S✻UUterm = o.S✻Y₂y₁ - o.S✻DGPZY₂' * view(o.DGP.β̈ ,:,1) - (o.DGP.Π̈ )'* (o.S✻Xy₁ - o.S✻XDGPZ * o.DGP.β̈)
				o.S✻u₁u₁ .= o.S✻y₁y₁ .- (2 * o.S✻DGPZy₁ - o.S✻DGPZDGPZ * o.DGP.β̈ )'o.DGP.β̈  .+ (2 * S✻UUterm - Π̂S✻XÜ₂γ̈Y + S✻Y₂Ü₂γ̈Y )'o.DGP.γ̈Y 
				o.S✻U₂paru₁ .= o.Repl.RparY' * (S✻UUterm + S✻Y₂Ü₂γ̈Y - Π̂S✻XÜ₂γ̈Y )
				o.S✻U₂parU₂par .= o.Repl.RparY' * (S✻Ü₂Y₂ - o.S✻XU₂' * o.DGP.Π̈ ) * o.Repl.RparY

				if o.DGP.restricted
					r₁S✻DGPZR₁y₁ = r₁' * o.S✻DGPZR₁y₁
					o.S✻u₁u₁ .+= -2 .* r₁S✻DGPZR₁y₁ .+ r₁' * (o.S✻DGPZR₁DGPZR₁ * r₁) .+ 2 .* r₁' * (o.S✻DGPZR₁DGPZ * o.DGP.β̈ + (o.S✻DGPZR₁X * o.DGP.Π̈- o.S✻DGPZR₁Y₂) * o.DGP.γ̈Y )
					o.S✻U₂paru₁ .-= o.Repl.RparY' * (o.S✻DGPZR₁Y₂ - o.S✻DGPZR₁X * o.DGP.Π̈)' * r₁
				end

				o.S✻ȳ₁u₁ .= Π⃛y'o.S✻Xu₁
				o.S✻ȳ₁Ü₂par .= Π⃛y'o.S✻XU₂par
				o.S✻Z̄u₁ .= o.Π̈Rpar'o.S✻Xu₁
				o.S✻Z̄Ü₂par .= o.Π̈Rpar'o.S✻XU₂par
			end

			if o.willfill || o.not2SLS  # make Z⟂U
				o.S✻Zperpu₁              .= o.S✻Zperpy₁; t✻minus!(o.S✻Zperpu₁, o.S✻ZperpDGPZ, o.DGP.β̈ ); t✻plus!(o.S✻Zperpu₁, o.S✻ZperpU₂, o.DGP.γ̈Y )
				o.invZperpZperpS✻Zperpu₁ .= o.invZperpZperpS✻Zperpy₁; t✻minus!(o.invZperpZperpS✻Zperpu₁, o.invZperpZperpS✻ZperpDGPZ, o.DGP.β̈ ); t✻plus!(o.invZperpZperpS✻Zperpu₁, o.invZperpZperpS✻ZperpU₂, o.DGP.γ̈Y )
				if o.DGP.restricted
					t✻minus!(o.S✻Zperpu₁, o.S✻ZperpDGPZR₁, r₁)
					t✻minus!(o.invZperpZperpS✻Zperpu₁, o.invZperpZperpS✻ZperpDGPZR₁, r₁)
				end

				if o.NFE>0 && !o.FEboot
					t✻!(o.CT✻FEU₂, o.CT✻FEX, o.DGP.Π̈ ); lsub!(o.CT✻FEU₂, o.CT✻FEY₂)
					t✻!((@view o.CT✻FEU[1:1]), o.CT✻FEZ, o.DGP.β̈ ); lsub!((@view o.CT✻FEU[1:1]), o.CT✻FEy₁); t✻plus!((@view o.CT✻FEU[1:1]), o.CT✻FEU₂, o.DGP.γ̈Y)
					t✻!((@view o.CT✻FEU[2:end]), o.CT✻FEU₂, o.Repl.RparY)
					o.DGP.restricted &&
						t✻minus!((@view o.CT✻FEU[1:1]), o.CT✻FEZR₁, r₁)
					broadcast!(.*, o.invFEwtCT✻FEU, Ref(o.invsumFEwt), o.CT✻FEU)  # optimize by summing invsumFEwt to conform to zvars
				end
			end

			if o.willfill
				t✻!(o.S✻⋂XU₂, o.S✻⋂XX, o.DGP.Π̈) ; o.S✻⋂XU₂ .= o.S✻⋂XY₂ .- o.S✻⋂XU₂
				o.S✻⋂Xu₁ .= o.S✻⋂Xy₁; t✻minus!(o.S✻⋂Xu₁, o.S✻⋂XDGPZ, o.DGP.β̈ ); t✻plus!(o.S✻⋂Xu₁, o.S✻⋂XU₂, o.DGP.γ̈Y )
				o.DGP.restricted &&
					t✻minus!(o.S✻⋂Xu₁, o.S✻⋂X_DGPZR₁, r₁)
				o.S✻⋂XÜ₂par .= o.S✻⋂XU₂ * o.Repl.RparY
			end
		elseif o.willfill  # for coarse, jk, construct this input in granular fashion rather than in for coarse, non-jk above
			panelcross21!(o.S✻⋂Xu₁, o.DGP.X₁, o.DGP.X₂, o.DGP.u⃛₁, o.info✻⋂)
			panelcross21!(o.S✻⋂XÜ₂par, o.DGP.X₁, o.DGP.X₂, o.Ü₂par, o.info✻⋂)
		end

		if o.willfill
			t✻!(o.S⋂ȳ₁X   , o.DGP.γ⃛', o.S⋂XX)
			t✻!(o.S⋂ReplZ̄X, o.Π̈Rpar' , o.S⋂XX)

			@inbounds for j ∈ 0:o.Repl.kZ
				if o.Repl.Yendog[j+1]
					t✻!(o.negS✻UMZperpX[j+1], o.S⋂XZperpinvZperpZperp, o.S✻ZperpU[j+1])  # S_* diag⁡(U ̈_(∥j) ) Z_⊥ (Z_⊥^' Z_⊥ )^(-1) Z_(⊥g)^' X_(∥g)
					o.negS✻UMZperpX[j+1][o.crosstab⋂✻ind] .-= vec(j>0 ? view(o.S✻⋂XÜ₂par,:,:,j) : view(o.S✻⋂Xu₁,:,:,1))
					if o.NFE>0 && !o.FEboot
						for i ∈ 1:o.DGP.kX
							t✻!(view(o.negS✻UMZperpX[j+1],i,:,:), o.CT⋂FEX[i]', o.CT✻FEU[j+1])  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
						end
					end
				end
			end
		end
	end

	if o.not2SLS
		@inbounds for i ∈ 0:o.Repl.kZ, j ∈ 0:i
			o.S✻ȲUfold[i+1,:,j+1] .= o.S✻ȲU[i+1,j+1] .+ o.S✻ȲU[j+1,i+1]
			o.S✻ȲUfold[j+1,:,i+1] .= o.S✻ȲUfold[i+1,:,j+1]
		end
	end
	nothing
end

# For WRE, and with reference to Y = [y₁ Z], given 0-based columns indexes within it, i, j, return all bootstrap realizations of 
# Y[:,i]'((1-κ)*M_Zperp-κ*M_Xpar)*Y[:,j] for κ constant across replications
# i can be a rowvector
# (only really the Hessian when we narrow Y to Z)
function HessianFixedkappa(o::StrBootTest{T}, is::Vector{S} where S<:Integer, j::Integer, κ::Number, _jk::Bool) where T
  dest = Matrix{T}(undef, length(is), o.ncolsv)
  HessianFixedkappa!(o, dest, is, j, κ, _jk)
  dest
end
function HessianFixedkappa!(o::StrBootTest{T}, dest::AbstractMatrix{T}, is::Vector{S} where S<:Integer, j::Integer, κ::Number, _jk::Bool) where T
	if !iszero(κ)
		if iszero(j)
			o.T1R .= o.DGP.γ⃛
		else
			fillcols!(o.T1R, o.invXXXZ̄, j)
		end
	end
	if o.Repl.Yendog[j+1] && any(o.Repl.Yendog[is.+1])
		if !iszero(κ)
			t✻plus!(o.T1R, o.invXXS✻XU[j+1], o.v)
		end
		if !isone(κ)
			mul!(o.S✻ZperpUv, o.S✻ZperpU[j+1], o.v)
			o.NFE>0 && !o.FEboot &&
				mul!(o.invFEwtCT✻FEUv, o.invFEwtCT✻FEU[j+1], o.v)
		end
	end

	@inbounds for (row,i) ∈ enumerate(is)
		if iszero(κ)
			dest[row,:] .= o.ȲȲ[i+1,j+1]; t✻plus!(view(dest,row:row,:), view(o.S✻ȲUfold,i+1,:,j+1)', o.v); 
			coldotplus!(dest, row, o.v, o.S✻UU[i+1, j+1], o.v)
			if o.Repl.Yendog[i+1] && o.Repl.Yendog[j+1]
				mul!(o.invZperpZperpS✻ZperpUv, o.invZperpZperpS✻ZperpU[i+1], o.v)
				coldotminus!(dest, row, o.invZperpZperpS✻ZperpUv, o.S✻ZperpUv)
			end

			if o.NFE>0 && !o.FEboot
				mul!(o.CT✻FEUv, o.CT✻FEU[i+1], o.v)
				coldotminus!(dest, row, o.CT✻FEUv, o.invFEwtCT✻FEUv)
			end
		else
	    fillcols!(o.T1L, o.XȲ,i+1)  # X_∥^' Y_(∥i)
	    o.Repl.Yendog[i+1] &&
				t✻plus!(o.T1L, o.S✻XU[i+1], o.v)
			coldot!(dest, row, o.T1L, o.T1R)  # multiply in the left-side linear term

		  if !isone(κ)
		    _dest = t✻(view(o.S✻ȲUfold,i+1,:,j+1)', o.v); _dest .+= o.ȲȲ[i+1,j+1]
				coldotplus!(_dest, 1, o.v, o.S✻UU[i+1, j+1], o.v)
		    if o.Repl.Yendog[i+1] && o.Repl.Yendog[j+1]
					mul!(o.invZperpZperpS✻ZperpUv, o.invZperpZperpS✻ZperpU[i+1], o.v)
		      coldotminus!(_dest, 1, o.invZperpZperpS✻ZperpUv, o.S✻ZperpUv)
				end

		    if o.NFE>0 && !o.FEboot
					mul!(o.CT✻FEUv, o.CT✻FEU[i+1], o.v)
					coldotminus!(_dest, 1, o.CT✻FEUv, o.invFEwtCT✻FEUv)
				end

		    view(dest,row,:) .*= κ;  view(dest,row:row,:) .+= (1-κ) .* _dest  # "view(dest,row:row,:)" avoids allocations where dest[row:row,:] doesn't
			end
		end

		if _jk
			!iszero(κ) &&
				(dest[row,1] = dot(i>0 ? o.Repl.XZ[:,i] : o.Repl.Xy₁par, 
				                   j>0 ? o.Repl.V[ :,j] : o.Repl.invXXXy₁par))
			!isone(κ) &&
				(dest[row,1] = iszero(κ) ? o.Repl.YY[i+1,j+1] : κ * dest[row,1] + (1 - κ) * o.Repl.YY[i+1,j+1])
		end
	end
  nothing
end

# Workhorse for WRE CRVE sandwich filling
# Given a zero-based column index, i>0, and a matrix β̈s of all the bootstrap estimates, 
# return all bootstrap realizations of P_X * Z̄[:,i]_g ' û₁g^*b
# for all groups in the intersection of all error clusterings
# return value has one row per ⋂ cluster, one col per bootstrap replication
# that is, given i, β̈s = δ ̂_CRκ^(*), return, over all g, b (P_(X_∥ g) Z_(∥i)^(*b) )^' (M_(Z_⊥ ) y_(1∥)^(*b) )_g-(P_(X_∥ g) Z_(∥i)^(*b) )^' (M_(Z_⊥ ) Z_∥^(*b) )_g δ ̂_CRκ^(*b)
function Filling!(o::StrBootTest{T}, dest::AbstractMatrix{T}, i::Int64, β̈s::AbstractMatrix{T}, _jk::Bool) where T
	if o.granular
    t✻!(o.S✻XUv, o.S✻XU[i+1], o.v)
		t✻!(o.S✻UZperpinvZperpZperpv, o.invZperpZperpS✻ZperpU[1], o.v)
		o.NFE>0 && !o.FEboot &&
			t✻!(o.invFEwtCT✻FEUv, o.invFEwtCT✻FEU[1], o.v)
		if o.purerobust
			@inbounds for g ∈ 1:o.N⋂
				o.PXY✻ .= o.PXZ̄[g,i]
				o.Repl.Yendog[i+1] &&
					t✻plus!(o.PXY✻, view(o.XinvXX,g:g,:), o.S✻XUv)

				t✻!(o.S✻UMZperp, view(o.Repl.Zperp,g:g,:), o.S✻UZperpinvZperpZperpv); t✻minus!(o.S✻UMZperp, o.DGP.u⃛₁[g], view(o.v,g:g,:))
				o.NFE>0 && !o.FEboot &&
					(o.S✻UMZperp .+= view(o.invFEwtCT✻FEUv, o._FEID[g]:o._FEID[g], :))  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE

				t✻!(view(dest,g:g,:), o.DGP.ȳ₁[g], o.PXY✻)
				coldotminus!(dest, g, o.PXY✻, o.S✻UMZperp)
			end
		else
			@inbounds for (g,S) ∈ enumerate(o.info⋂)
				S✻UMZperpg = rowsubview(o.S✻UMZperp, length(S))
				PXY✻g  = rowsubview(o.PXY✻, length(S))
				PXY✻g .= view(o.PXZ̄,S,i)
				o.Repl.Yendog[i+1] &&
					t✻plus!(PXY✻g, view(o.XinvXX,S,:), o.S✻XUv)

				t✻!(S✻UMZperpg, view(o.Repl.Zperp,S,:), o.S✻UZperpinvZperpZperpv); S✻UMZperpg .-= view(o.DGP.u⃛₁, S) .*  view(o.v, view(o.ID✻, S),:)
				o.NFE>0 && !o.FEboot &&
					(S✻UMZperpg .+= view(o.invFEwtCT✻FEUv, view(o._FEID,S), :))  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE

				t✻!(view(dest,g,:), PXY✻g', view(o.DGP.ȳ₁, S))
				coldotminus!(dest, g, PXY✻g, S✻UMZperpg)
			end
		end
		@inbounds for j ∈ 1:o.Repl.kZ
			_β̈  = view(β̈s,j:j,:)
			if o.Repl.Yendog[j+1]
				matbyrow!(o.β̈v, o.v, β̈s, j)
				t✻!(o.S✻UZperpinvZperpZperpv, o.invZperpZperpS✻ZperpU[j+1], o.β̈v)
				o.NFE>0 && !o.FEboot &&
					t✻!(o.invFEwtCT✻FEUv, o.invFEwtCT✻FEU[j+1], o.β̈v)
			end
			if o.purerobust
				@inbounds for g ∈ 1:o.N⋂
					o.PXY✻ .= o.PXZ̄[g,i]
					o.Repl.Yendog[i+1] &&
						t✻plus!(o.PXY✻, view(o.XinvXX,g:g,:), o.S✻XUv)
	
					if o.Repl.Yendog[j+1]
						t✻!(o.S✻UMZperp, view(o.Repl.Zperp,g:g,:), o.S✻UZperpinvZperpZperpv); t✻minus!(o.S✻UMZperp, o.Ü₂par[g,j], view(o.β̈v,g:g,:))
						o.NFE>0 && !o.FEboot &&
							(o.S✻UMZperp .+= view(o.invFEwtCT✻FEUv, o._FEID[g]:o._FEID[g], :))  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
	
						# dest[g:g,:] .+= o.PXY✻ .* (o.S✻UMZperp .- o.Z̄[g,j] .* _β̈ )
						@tturbo for c ∈ indices((dest, o.PXY✻, _β̈ , o.S✻UMZperp), 2)
							dest[g,c] += o.PXY✻[c] * (o.S✻UMZperp[c] - o.Z̄[g,j] * _β̈[c])
						end
					else
						# dest[g:g,:] .-= o.Z̄[g,j] .* o.PXY✻ .* _β̈ 
						@tturbo for c ∈ indices((dest, o.PXY✻, _β̈ ), 2)
							dest[g,c] -= o.Z̄[g,j] * o.PXY✻[c] * _β̈[c]
						end
					end
				end
			else
				for (g,S) ∈ enumerate(o.info⋂)
					PXY✻g  = rowsubview(o.PXY✻, length(S))
					PXY✻g .= view(o.PXZ̄,S,i)
					o.Repl.Yendog[i+1] &&
						t✻plus!(PXY✻g, view(o.XinvXX,S,:), o.S✻XUv)
	
					if o.Repl.Yendog[j+1]
						S✻UMZperpg = rowsubview(o.S✻UMZperp, length(S))
						t✻!(S✻UMZperpg, view(o.Repl.Zperp,  S,:), o.S✻UZperpinvZperpZperpv)
						S✻UMZperpg .-= view(o.Ü₂par, S, j) .* view(o.β̈v,view(o.ID✻, S),:)
						# IDS = view(o.ID✻,S)
						# @tturbo for c ∈ indices((S✻UMZperpg,o.β̈v),2), r ∈ indices((S✻UMZperpg,S),1)
						# 	S✻UMZperpg[r,c] -= o.Ü₂par[S[r], j] * o.β̈v[IDS[r],c]
						# end

						t✻minus!(S✻UMZperpg, view(o.Z̄,S,j), _β̈ )
	
						o.NFE>0 && !o.FEboot &&
							(S✻UMZperpg .+= view(o.invFEwtCT✻FEUv, view(o._FEID,S), :))  # CT_(*,FE) (U ̈_(∥j) ) (S_FE S_FE^' )^(-1) S_FE
	
						coldotplus!(dest, g, PXY✻g, S✻UMZperpg)
					else
						coldotminus!(dest, g, PXY✻g, view(o.Z̄,S,j) * _β̈)
					end
				end
			end
		end
	else  # coarse error clustering
		fillcols!(o.F₁, o.invXXXZ̄, i)
		o.Repl.Yendog[i+1] && t✻plus!(o.F₁, o.invXXS✻XU[i+1], o.v)
    @inbounds for g ∈ 1:o.N⋂
			fillcols!(o.F₂, o.S⋂ȳ₁X, 1, g)
			t✻minus!(o.F₂, view(o.negS✻UMZperpX[1],:,g,:), o.v)
      coldot!(dest, g, o.F₁, o.F₂)
		end
		@inbounds for j ∈ 1:o.Repl.kZ
	    matbyrow!(o.F₁β, o.F₁, β̈s, j)
	    for g ∈ 1:o.N⋂
	      fillcols!(o.F₂, o.S⋂ReplZ̄X, j, g)
				o.Repl.Yendog[j+1] && t✻minus!(o.F₂, view(o.negS✻UMZperpX[j+1],:,g,:), o.v)
	      coldotminus!(dest, g, o.F₁β, o.F₂)
			end
		end
	end
	if _jk
		o.uⱼₖ .= o.Repl.y₁par; t✻minus!(o.uⱼₖ, o.Repl.Zpar, β̈s[:,1])
		panelsum!(view(dest,:,1), view(o.Repl.PXZ,:,i), o.uⱼₖ, o.info⋂)
	end
  nothing
end

function MakeWREStats!(o::StrBootTest{T}, w::Integer) where T
	_jk = o.jk & w==1

  if isone(o.Repl.kZ)  # optimized code for 1 retained coefficient in bootstrap regression
		if o.liml
			HessianFixedkappa!(o, o.YY₁₁  , [0], 0, zero(T), _jk)  # κ=0 => Y*MZperp*Y
			HessianFixedkappa!(o, o.YY₁₂  , [0], 1, zero(T), _jk)
			HessianFixedkappa!(o, o.YY₂₂  , [1], 1, zero(T), _jk)
			HessianFixedkappa!(o, o.YPXY₁₁, [0], 0, one(T) , _jk)  # κ=1 => Y*PXpar*Y
			HessianFixedkappa!(o, o.YPXY₁₂, [0], 1, one(T) , _jk)
			HessianFixedkappa!(o, o.YPXY₂₂, [1], 1, one(T) , _jk)
			o.YY₁₂YPXY₁₂ .= o.YY₁₂ .* o.YPXY₁₂
			o.x₁₁ .= o.YY₂₂ .* o.YPXY₁₁ .- o.YY₁₂YPXY₁₂      # elements of o.YY✻^-1 * o.YPXY✻ up to factor of det(o.YY✻)
			o.x₁₂ .= o.YY₂₂ .* o.YPXY₁₂ .- o.YY₁₂ .* o.YPXY₂₂
			o.x₂₁ .= o.YY₁₁ .* o.YPXY₁₂ .- o.YY₁₂ .* o.YPXY₁₁
			o.x₂₂ .= o.YY₁₁ .* o.YPXY₂₂ .- o.YY₁₂YPXY₁₂
			o.κs .= (o.x₁₁ .+ o.x₂₂)./2
			o.κs .= 1 ./ (1 .- (o.κs .- sqrtNaN.(o.κs.^2 .- o.x₁₁ .* o.x₂₂ .+ o.x₁₂ .* o.x₂₁)) ./ 
			                   (o.YY₁₁ .* o.YY₂₂ .- o.YY₁₂ .* o.YY₁₂))  # solve quadratic equation for smaller eignenvalue; last term is det(o.YY✻)
			!iszero(o.fuller) && (o.κs .-= o.fuller / (o._Nobs - o.kX))
			o.β̈sAs[2:2,:] .= o.κs .* (o.YPXY₂₂ .- o.YY₂₂) .+ o.YY₂₂
			o.β̈sAs[1:1,:] .= (o.κs .* (o.YPXY₁₂ .- o.YY₁₂) .+ o.YY₁₂) ./ view(o.β̈sAs,2:2,:)
		else
			HessianFixedkappa!(o, o.β̈sAs, [0; 1], 1, o.κ, _jk)
			o.β̈sAs[1,:] ./= view(o.β̈sAs,2,:)  # o.β̈s ./= _As
		end

		if o.null
			o.numerWRE .= view(o.β̈sAs,1:1,:) .+ (o.Repl.Rt₁ - o.r) / o.Repl.RRpar
		else
			o.numerWRE .= view(o.β̈sAs,1:1,:) .- view(o.DGP.β̈ ,:,1)
			isone(w) && (o.numerWRE[1] = o.β̈sAs[1] + (o.Repl.Rt₁[] - o.r[]) / o.Repl.RRpar[])
		end

		@storeWtGrpResults!(o.numer, o.numerWRE)
		if o.bootstrapt
			if o.robust
				Filling!(o, o.J⋂s1[1], 1, view(o.β̈sAs,1:1,:), _jk)

				_J = o.J⋂s1[1]
				@tturbo for j ∈ indices((o.J⋂s1[1],o.β̈sAs),2)  # ./= As
					t = o.β̈sAs[2,j]
					for i ∈ indices(o.J⋂s1[1], 1)
						_J[i,j] /= t
					end
				end

				coldot!(o.denom[1,1], o.clust[1].multiplier, o.J⋂s1[1])
				@inbounds for c ∈ 2:o.NErrClustCombs  # sum sandwich over error clusteringssrc/WRE.jl
					nrows(o.clust[c].order)>0 && 
						(o.J⋂s1[c] .= o.J⋂s1[c-1][o.clust[c].order,:])
					panelsum!(dropdims(o.Jc[c]; dims=3), o.J⋂s1[c], o.clust[c].info)
		    	coldotplus!(o.denom[1,1], o.clust[c].multiplier, dropdims(o.Jc[c]; dims=3))
				end
			else
				o.denom[1,1] .= (HessianFixedkappa(o, [0], 0, zero(T), _jk) .-   # XXX rewrite to avoid allocations
				                     2 .* view(o.β̈sAs,1:1,:) .* HessianFixedkappa(o, [0], 1, zero(T), _jk) .+ 
								             view(o.β̈sAs,1:1,:).^2 .* HessianFixedkappa(o, [1], 1, zero(T), _jk)) ./ view(o.β̈sAs,2:2,:)  # classical error variance
			end
			@storeWtGrpResults!(o.dist, o.sqrt ? o.numerWRE ./ sqrtNaN.(o.denom[1,1]) : o.numerWRE .^ 2 ./ o.denom[1,1])
			o.denom[1,1] .*= o.Repl.RRpar[]^2
		end
		w==1 && o.bootstrapt && (o.statDenom = [o.denom[1,1][1];;])  # original-sample denominator

  else  # WRE bootstrap for more than 1 retained coefficient in bootstrap regression

		if o.liml
			@inbounds for i ∈ 0:o.Repl.kZ
				HessianFixedkappa!(o, view(o.YY✻  , 1:i+1, :, i+1), collect(0:i), i, zero(T), _jk)  # κ=0 => Y*MZperp*Y
				HessianFixedkappa!(o, view(o.YPXY✻, 1:i+1, :, i+1), collect(0:i), i,  one(T), _jk)  # κ=1 => Y*PXpar*Y
			end
			symmetrize!(o.YY✻)
			symmetrize!(o.YPXY✻)

			M = Matrix{T}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
			@inbounds for b ∈ eachindex(axes(o.κWRE,2))
				ldiv!(M, bunchkaufman(view(o.YY✻,:,b,:)), view(o.YPXY✻,:,b,:))
				o.κWRE[b] = one(T)/(one(T) - real(eigvalsNaN(M)[1]))
			end
			# view(o.κWRE,1,:,1) .= one(T) ./ (one(T) .- getindex.(real.(eigvalsNaN.(each(invsym(o.YY✻) * o.YPXY✻))), 1))
			!iszero(o.fuller) && (o.κWRE .-= o.fuller / (o._Nobs - o.kX))

			o.As .= o.κWRE .* view(o.YPXY✻, 2:o.Repl.kZ+1, :, 2:o.Repl.kZ+1) .+ (1 .- o.κWRE) .* view(o.YY✻, 2:o.Repl.kZ+1, :, 2:o.Repl.kZ+1)
			invsym!(o.As)
			t✻!(view(o.β̈s,:,:,1:1), o.As, o.κWRE .* view(o.YPXY✻, 2:o.Repl.kZ+1, :, 1) .+ (1 .- o.κWRE) .* view(o.YY✻, 2:o.Repl.kZ+1, :,  1))
		else
			HessianFixedkappa!(o, o.δnumer, collect(1:o.Repl.kZ), 0, o.κ, _jk)
			@inbounds for i ∈ 1:o.Repl.kZ
				HessianFixedkappa!(o, view(o.As, 1:i, :, i), collect(1:i), i, o.κ, _jk)
			end
			symmetrize!(o.As)
			invsym!(o.As)
			t✻!(view(o.β̈s,:,:,1:1), o.As, view(o.δnumer,:,:,1:1))
		end

		if o.bootstrapt
			if o.robust
				@inbounds for i ∈ 1:o.Repl.kZ
					Filling!(o, view(o.J⋂s,:,:,i), i, o.β̈s, _jk)
				end
			else
				@inbounds for i ∈ 0:o.Repl.kZ
					HessianFixedkappa!(o, view(o.YY✻, 1:i+1, :, i+1), collect(0:i), i, zero(T), _jk)
				end
				symmetrize!(o.YY✻)
			end
		end

		if o.null
			o.numerWRE .= o.Repl.Rt₁ - o.r                ; t✻plus!(o.numerWRE, o.Repl.RRpar, o.β̈s)
		else
			o.numerWRE .= o.Repl.RRpar * (-o.DGP.β̈[:,1:1]); t✻plus!(o.numerWRE, o.Repl.RRpar, o.β̈s); 
			w==1 && (o.numerWRE[:,1:1] .= o.Repl.RRpar * o.β̈s[:,1:1] .+ o.Repl.Rt₁ .- o.r)
		end

		if o.bootstrapt
			if o.robust  # Compute denominator for this WRE test stat
				t✻!(o.ARpars, o.As, o.Repl.RRpar')
				t✻!(o.J⋂ARpars[1], o.J⋂s, o.ARpars)
				t✻!(o.denomWRE, o.clust[1].multiplier, o.J⋂ARpars[1]', o.J⋂ARpars[1])
				for c ∈ 2:o.NErrClustCombs
					(!isone(o.NClustVar) && nrows(o.clust[c].order)>0) &&
						(o.J⋂ARpars[c] .= o.J⋂ARpars[c-1][o.clust[c].order,:,:])
					panelsum!(reshape(o.Jc[c], o.clust[c].N, :), reshape(o.J⋂ARpars[c], o.N⋂, :), o.clust[c].info)
					t✻plus!(o.denomWRE, o.clust[c].multiplier, o.Jc[c]', o.Jc[c])
				end
			else  # non-robust
				tmp = view([fill(T(-1), 1, o.ncolsv) ; o.β̈s], :, :, 1:1)
				o.denomWRE .= (o.Repl.RRpar * o.As * o.Repl.RRpar') .* (tmp'o.YY✻ * tmp)  # 2nd half is sig2 of errors
			end
			if w==1
				o.statDenom = o.denomWRE[:,1,:]
				o.numer = o.numerWRE[:,1,1:1]  # just save full-sample numerator
			end
			if o.sqrt
				@storeWtGrpResults!(o.dist, o.numerWRE ./ sqrtNaN.(dropdims(o.denomWRE; dims=3)))
			else
				invsym!(o.denomWRE)
				_numer = view(o.numerWRE,:,:,1:1)
				@storeWtGrpResults!(o.dist, dropdims(_numer'o.denomWRE*_numer; dims=3))  # hand-code for 2-dimensional?  XXX allocations
			end
		else
			@storeWtGrpResults!(o.numer, o.numerWRE)
		end
	end
	nothing
end
