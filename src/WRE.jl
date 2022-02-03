# stuff done once per exucution--not depending on r
function InitWRE!(o::StrBootTest{T}) where T
	o.LIML && o.Repl.kZ==1 && o.Nw==1 && (o.As = o.β̈s = zeros(1, o.B+1))

	o.invZperpZperpS✻ZperpU = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻ZperpU              = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻YU                  = Matrix{Vector{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
	o.invXXS✻XU             = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻XU                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
	o.S✻UU                  = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
	o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kX, ncols(o.v))] :
												[Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[end]))]
	o.T1R = deepcopy(o.T1L)
	o.S✻UPX_S✻UMZperp = Matrix{Array{T,3}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)

	if o.bootstrapt
		o.δdenom_b = zeros(o.Repl.kZ, o.Repl.kZ)
		o.S✻UMZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
		o.granular && (o.S✻UPX     = Vector{Matrix{T}}(undef, o.Repl.kZ+1))
		o._J⋂ = zeros(o.clust[1].N, o.Repl.kZ)
		!o.granular && (o.S✻⋂u₁XinvXX = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.N✻))
		if o.LIML || !o.robust
			o.YY✻_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
			o.YPXY✻_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
		end
		o.NFE>0 && (o.bootstrapt || !isone(o.κ) || o.LIML) && (o.CTFEU = Vector{Matrix{T}}(undef, o.Repl.kZ+1))
	end

	info✻_✻⋂ = panelsetup(o.ID✻⋂, 1:o.nbootclustvar)

	o.S✻⋂XY₂  = panelcross21(o, o.Repl.X₁, o.Repl.X₂, o.DGP.Y₂, o.info✻⋂)
	o.S✻⋂XX   = panelcross22(o, o.Repl.X₁, o.Repl.X₂, o.Repl.X₁, o.Repl.X₂, o.info✻⋂)
	o.S✻⋂X_DGPZ   = panelcross21(o, o.Repl.X₁, o.Repl.X₂, o.DGP.Z, o.info✻⋂)
	o.S✻⋂Xy₁  = panelcross21(o, o.Repl.X₁, o.Repl.X₂, o.DGP.y₁, o.info✻⋂)
	o.S✻⋂X_DGPZR₁ = panelcross21(o, o.Repl.X₁, o.Repl.X₂, o.DGP.ZR₁, o.info✻⋂)
	  S✻⋂ZperpX   = panelcross12(o, o.Repl.Zperp, o.Repl.X₁, o.Repl.X₂, o.info✻⋂)

	o.invXXS✻⋂XY₂  = o.Repl.invXX * o.S✻⋂XY₂ 
	o.invXXS✻⋂XX   = o.Repl.invXX * o.S✻⋂XX  
	o.invXXS✻⋂X_DGPZ   = o.Repl.invXX * o.S✻⋂X_DGPZ 
	o.invXXS✻⋂Xy₁  = o.Repl.invXX * o.S✻⋂Xy₁ 
	o.invXXS✻⋂X_DGPZR₁ = o.Repl.invXX * o.S✻⋂X_DGPZR₁

	o.S✻XY₂  = @panelsum(o, o.S✻⋂XY₂ , info✻_✻⋂)
	o.S✻XX   = @panelsum(o, o.S✻⋂XX  , info✻_✻⋂)
	o.S✻X_DGPZ   = @panelsum(o, o.S✻⋂X_DGPZ  , info✻_✻⋂)
	o.S✻Xy₁  = dropdims(@panelsum(o, reshape(o.S✻⋂Xy₁,Val(3)), info✻_✻⋂); dims=3)
	o.S✻XZR₁ = @panelsum(o, o.S✻⋂X_DGPZR₁, info✻_✻⋂)
	o.S✻ZperpY₂  = panelcross11(o, o.Repl.Zperp, o.DGP.Y₂, o.info✻)
	o.S✻ZperpX   = @panelsum(o, S✻⋂ZperpX, info✻_✻⋂)
	o.S✻Zperp_DGPZ   = panelcross11(o, o.Repl.Zperp, o.DGP.Z, o.info✻)
	o.S✻Zperpy₁  = panelcross11(o, o.Repl.Zperp, o.DGP.y₁, o.info✻)
	o.S✻Zperp_DGPZR₁ = panelcross11(o, o.Repl.Zperp, o.DGP.ZR₁, o.info✻)

	if o.NFE>0 && (o.LIML || !isone(o.κ) || o.bootstrapt)
		  CT✻⋂FEX = cat(crosstabFE(o, o.Repl.X₁, o.info✻⋂), crosstabFE(o, o.Repl.X₂, o.info✻⋂); dims=3)
		o.CT✻FEX  = @panelsum(o, CT✻⋂FEX, info✻_✻⋂)
		o.CT✻FEY₂ = crosstabFE(o, o.DGP.Y₂, o.info✻)
		o.CT✻FEZ  = crosstabFE(o, o.DGP.Z, o.info✻)
		o.CT✻FEy₁ = crosstabFE(o, o.DGP.y₁, o.info✻)
		o.CT✻FEZR₁ = crosstabFE(o, o.DGP.ZR₁, o.info✻)
	end

	if o.bootstrapt & o.robust & iszero(o.granular)
		o.info⋂_✻⋂ = panelsetup(o.ID✻⋂, o.subcluster+1:o.NClustVar)

		inds = o.subcluster>0 ?
						[CartesianIndex(i,j) for (j,v) ∈ enumerate(o.info⋂_✻⋂) for i ∈ v] :  # crosstab *,∩ is tall
						o.NClustVar == o.nbootclustvar ?
								[CartesianIndex(i,i) for i ∈ 1:o.N✻⋂] :  # crosstab *,∩ is square
								[CartesianIndex(j,i) for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab *,∩ is wide
		inds = [CartesianIndex(k,I) for I ∈ inds for k ∈ 1:o.Repl.kX]
		o.crosstab⋂✻ind = LinearIndices(FakeArray(Tuple(max(inds...))...))[inds]

		o.S⋂XZperpinvZperZperp = (@panelsum(o, S✻⋂ZperpX, o.info⋂_✻⋂))' * o.Repl.invZperpZperp

		o.NFE>0 && (o.CT⋂FEX = o.invFEwt .* @panelsum(o, CT✻⋂FEX, o.info⋂_✻⋂))
	end

	if o.LIML || !o.robust || !isone(o.κ)
		o.S✻Y₂Y₂ = panelcross11(o, o.DGP.Y₂, o.DGP.Y₂, o.info✻)
		o.S✻DGPZ_DGPZ = panelcross11(o, o.DGP.Z, o.DGP.Z, o.info✻)
		o.S✻DGPZ_Y₂ = panelcross11(o, o.DGP.Z, o.DGP.Y₂, o.info✻)
		o.S✻y₁Y₂  = panelcross11(o, o.Repl.y₁, o.DGP.Y₂, o.info✻)   
		o.S✻y₁X   = panelcross12(o, o.Repl.y₁, o.Repl.X₁, o.Repl.X₂, o.info✻)   
		o.S✻y₁_DGPZ   = panelcross11(o, o.Repl.y₁, o.DGP.Z, o.info✻)   
		o.S✻y₁y₁  = panelcross11(o, o.Repl.y₁, o.DGP.y₁, o.info✻)   
		o.S✻y₁_DGPZR₁ = panelcross11(o, o.Repl.y₁, o.DGP.ZR₁, o.info✻)
		ReplZR₁r₁ = o.Repl.ZR₁ * o.r₁
		o.S✻ReaplZR₁r₁_Y₂  = panelcross11(o, ReplZR₁r₁, o.DGP.Y₂, o.info✻)   
		o.S✻ReplZR₁r₁_X   = panelcross12(o, ReplZR₁r₁, o.Repl.X₁, o.Repl.X₂, o.info✻)   
		o.S✻ReplZR₁r₁_DGPZ   = panelcross11(o, ReplZR₁r₁, o.DGP.Z, o.info✻)   
		o.S✻ReplZR₁r₁_y₁  = panelcross11(o, ReplZR₁r₁, o.DGP.y₁, o.info✻)   
		o.S✻ReplZR₁r₁_DGPZR₁ = panelcross11(o, ReplZR₁r₁, o.DGP.ZR₁, o.info✻)   
		o.S✻ReplZ_Y₂  = panelcross11(o, o.Repl.Z, o.DGP.Y₂, o.info✻)   
		o.S✻ReplZ_X   = panelcross12(o, o.Repl.Z, o.Repl.X₁, o.Repl.X₂, o.info✻)   
		o.S✻ReplZ_DGPZ   = panelcross11(o, o.Repl.Z, o.DGP.Z, o.info✻)   
		o.S✻ReplZ_y₁  = panelcross11(o, o.Repl.Z, o.DGP.y₁, o.info✻)   
		o.S✻ReplZ_DGPZR₁ = panelcross11(o, o.Repl.Z, o.DGP.ZR₁, o.info✻)   
		o.S✻Y₂_DGPZR₁ = panelcross11(o, o.DGP.Y₂, o.DGP.ZR₁, o.info✻)
		o.S✻DGPZR₁_DGPZR₁ = panelcross11(o, o.DGP.ZR₁, o.DGP.ZR₁, o.info✻)
		o.S✻DGPZR₁_DGPZ = panelcross11(o, o.DGP.ZR₁, o.DGP.Z, o.info✻)
		o.S✻X_DGPZR₁ = @panelsum(o, o.S✻⋂X_DGPZR₁, info✻_✻⋂)
	end

	o.invXXS✻XY₂  = @panelsum(o, o.invXXS✻⋂XY₂ , info✻_✻⋂)
	o.invXXS✻XX   = @panelsum(o, o.invXXS✻⋂XX  , info✻_✻⋂)
	o.invXXS✻X_DGPZ   = @panelsum(o, o.invXXS✻⋂X_DGPZ  , info✻_✻⋂)
	o.invXXS✻Xy₁  = dropdims(@panelsum(o, reshape(o.invXXS✻⋂Xy₁,Val(3)), info✻_✻⋂); dims=3)
	o.invXXS✻X_ReplZR₁ = @panelsum(o, o.invXXS✻⋂X_DGPZR₁, info✻_✻⋂)
	o.invZperpZperpS✻ZperpY₂  = o.Repl.invZperpZperp * o.S✻ZperpY₂ 
	o.invZperpZperpS✻ZperpX   = o.Repl.invZperpZperp * o.S✻ZperpX  
	o.invZperpZperpS✻Zperp_DGPZ   = o.Repl.invZperpZperp * o.S✻Zperp_DGPZ  
	o.invZperpZperpS✻Zperpy₁  = o.Repl.invZperpZperp * o.S✻Zperpy₁ 
	o.invZperpZperpS✻Zperp_DGPZR₁ = o.Repl.invZperpZperp * o.S✻Zperp_DGPZR₁
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
		S✻y₁U₂ = o.S✻y₁Y₂ - o.S✻y₁X * o.DGP.Π̂
		S✻ZR₁r₁U₂ = o.S✻ReaplZR₁r₁_Y₂ - o.S✻ReplZR₁r₁_X * o.DGP.Π̂
		S✻ZR₁r₁U₂RparY = S✻ZR₁r₁U₂ * o.Repl.RparY
		S✻ZU₂ = o.S✻ReplZ_Y₂ - o.S✻ReplZ_X * o.DGP.Π̂
		S✻ZU₂RparY = S✻ZU₂ * o.Repl.RparY
		S✻y₁ZR₁r₁ = o.S✻y₁_DGPZR₁ * r₁
		S✻ZR₁r₁ZR₁r₁ = o.S✻ReplZR₁r₁_DGPZR₁ * r₁
		S✻ZZR₁r₁ = o.S✻ReplZ_DGPZR₁ * r₁
		γ̈S✻Ü₂XΠ̂ = o.DGP.γ̈' * S✻XU₂' * o.DGP.Π̂
		S✻Ü₂Y₂ = o.S✻Y₂Y₂ - o.DGP.Π̂' * o.S✻XY₂
		γ̈S✻Ü₂Y₂ = o.DGP.γ̈' * S✻Ü₂Y₂
	S✻Ü₂parÜ₂par = o.Repl.RparY' * (S✻Ü₂Y₂ - S✻XU₂' * o.DGP.Π̂) * o.Repl.RparY

		if o.NFE>0
			CT✻FEU = o.CT✻FEY₂ - o.CT✻FEX * o.DGP.Π̂
			CT✻FEURparY = CT✻FEU * o.Repl.RparY
		end
	end

  @inbounds for i ∈ 0:o.Repl.kZ  # precompute various clusterwise sums
		# S_✻(u .* X), S_✻(u .* Zperp) for residuals u for each endog var; store transposed
		if iszero(i)
			o.S✻XU[1]      = o.S✻Xy₁      - o.S✻X_DGPZ      * o.DGP.β̈ + S✻XU₂      * o.DGP.γ̈
			o.invXXS✻XU[1] = o.invXXS✻Xy₁ - o.invXXS✻X_DGPZ * o.DGP.β̈ + invXXS✻XU₂ * o.DGP.γ̈
			if o.DGP.restricted
				o.S✻XU[1]      .-= o.S✻XZR₁      * r₁
				o.invXXS✻XU[1] .-= o.invXXS✻X_ReplZR₁ * r₁
			end
		else
			o.S✻XU[i+1]      = view(S✻XU₂RparY,:,:,i)
			o.invXXS✻XU[i+1] = view(invXXS✻XU₂RparY,:,:,i)
		end

		if o.LIML || !isone(o.κ) || o.bootstrapt
			if iszero(i)
				o.S✻ZperpU[1]              = o.S✻Zperpy₁              - o.S✻Zperp_DGPZ              * o.DGP.β̈ + S✻ZperpU₂              * o.DGP.γ̈
				o.invZperpZperpS✻ZperpU[1] = o.invZperpZperpS✻Zperpy₁ - o.invZperpZperpS✻Zperp_DGPZ * o.DGP.β̈ + invZperpZperpS✻ZperpU₂ * o.DGP.γ̈
				if o.DGP.restricted
					o.S✻ZperpU[1]              .-= o.S✻Zperp_DGPZR₁              * r₁
					o.invZperpZperpS✻ZperpU[1] .-= o.invZperpZperpS✻Zperp_DGPZR₁ * r₁
				end
			else
				o.S✻ZperpU[i+1]              = view(S✻ZperpU₂RparY,:,:,i)
				o.invZperpZperpS✻ZperpU[i+1] = view(invZperpZperpS✻ZperpU₂RparY,:,:,i)
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
				o.S✻YU[1,1] = view(o.S✻y₁y₁ - o.S✻y₁_DGPZ * o.DGP.β̈ + S✻y₁U₂ * o.DGP.γ̈, 1,:,1)
				o.DGP.restricted &&
					(o.S✻YU[1,1] .-= view(S✻y₁ZR₁r₁, 1,:,1))

				if o.Repl.restricted
					o.S✻YU[1,1] .-= view(o.S✻ReplZR₁r₁_y₁ - o.S✻ReplZR₁r₁_DGPZ * o.DGP.β̈ + S✻ZR₁r₁U₂ * o.DGP.γ̈, 1,:,1)
					o.DGP.restricted &&
						(o.S✻YU[1,1] .+= view(S✻ZR₁r₁ZR₁r₁, 1,:,1))
				end

				for j ∈ 1:o.Repl.kZ
					o.S✻YU[j+1,1] = view(o.S✻ReplZ_y₁ - o.S✻ReplZ_DGPZ * o.DGP.β̈ + S✻ZU₂ * o.DGP.γ̈, j,:,1)
					o.DGP.restricted &&
						(o.S✻YU[j+1,1] .-= view(S✻ZZR₁r₁, j,:,1))
				end
      else
				o.S✻YU[1,i+1] = @view S✻y₁U₂[1,:,i]
				o.Repl.restricted &&
					(o.S✻YU[1,i+1] .-= @view S✻ZR₁r₁U₂RparY[1,:,i])
				for j ∈ 1:o.Repl.kZ
					o.S✻YU[j+1,i+1] = @view S✻ZU₂RparY[j,:,i]
				end
			end
			XXX = o.S✻y₁Y₂ - o.DGP.β̈' * o.S✻DGPZ_Y₂ - (o.S✻y₁X  - o.DGP.β̈' * o.S✻X_DGPZ') * o.DGP.Π̂
			# XXX = o.S✻y₁U₂ - o.DGP.β̈' * (o.S✻DGPZ_Y₂ + o.S✻X_DGPZ' * o.DGP.Π̂)
			if iszero(i)
		 		o.S✻UU[1,i+1] = dropdims(o.S✻y₁y₁ - (2 * o.S✻y₁_DGPZ + o.DGP.β̈' * o.S✻DGPZ_DGPZ) * o.DGP.β̈ +
																(2 * XXX - γ̈S✻Ü₂XΠ̂ + permutedims(γ̈S✻Ü₂Y₂,(3,2,1))) * o.DGP.γ̈; dims=1)
				o.DGP.restricted
					(o.S✻UU[1,i+1] .+= dropdims(-2 * S✻y₁ZR₁r₁ + r₁' * o.S✻DGPZR₁_DGPZR₁ * r₁ + 2 * r₁' * (o.S✻DGPZR₁_DGPZ * o.DGP.β̈ + (o.S✻X_DGPZR₁' * o.DGP.Π̂ - o.S✻Y₂_DGPZR₁') * o.DGP.γ̈); dims=1))
			else
				o.S✻UU[1,i+1] = view((XXX + γ̈S✻Ü₂Y₂ - γ̈S✻Ü₂XΠ̂) * o.Repl.RparY, 1,:,i)

				if o.DGP.restricted
					o.S✻UU[1,i+1] .-= view(r₁' * (o.S✻Y₂_DGPZR₁' - o.S✻X_DGPZR₁' * o.DGP.Π̂) * o.Repl.RparY, 1,:,i)
				end
			end
			for j ∈ 1:i
				o.S✻UU[j+1,i+1] = view(S✻Ü₂parÜ₂par, j,:,i)
			end
		end

		if o.robust && o.bootstrapt
			if iszero(o.granular)  # Within each bootstrap cluster, groupwise sum by all-error-cluster-intersections of u.*X and u.Zperp (and times invXX or invZperpZperp)
				if iszero(i)
					_S✻⋂UXinvXX = o.invXXS✻⋂Xy₁ - o.invXXS✻⋂X_DGPZ * o.DGP.β̈ + invXXS✻⋂XU₂ * o.DGP.γ̈
					o.DGP.restricted &&
						(_S✻⋂UXinvXX .-= o.invXXS✻⋂X_DGPZR₁ * r₁)
				else
					_S✻⋂UXinvXX = view(invXXS✻⋂XU₂RparY,:,:,i)
				end
				for g ∈ 1:o.N✻
					o.S✻⋂u₁XinvXX[i+1,g] = view(_S✻⋂UXinvXX,:,o._ID✻⋂[o.info✻[g]],1)'  # panelsum(o, o.Repl.XinvXX, u, o.infoCT⋂✻[g])
				end
			end

			if o.granular
				i>0 && (o.S✻UPX[i+1] = o.Repl.XinvXX * o.S✻XU[i+1])

				o.S✻UMZperp[i+1] = o.Repl.Zperp * o.invZperpZperpS✻ZperpU[i+1] 
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
					S✻UMZperpX[o.crosstab⋂✻ind] .-= vec(o.S✻⋂Xy₁ - o.S✻⋂X_DGPZ * o.DGP.β̈ + S✻⋂XU₂ * o.DGP.γ̈)
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
						(o.S✻UPX_S✻UMZperp[ind1+1,ind2+1] = o.S✻XU[ind1+1]'S✻UMZperpX) 
				end
			end
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
	Threads.@threads for i ∈ 1:o.clust[1].N
		PXY✻ = hcat(o.Repl.PXZ[i,ind1])
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
	Threads.@threads for i ∈ 1:o.clust[1].N
		S = o.info⋂[i]
		PXY✻ = o.Repl.Yendog[ind1+1] ? view(o.Repl.PXZ,S,ind1) .+ view(o.S✻UPX[ind1+1],S,:) * o.v :
																	reshape(o.Repl.PXZ[S,ind1], :, 1)

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
			PXY✻ = reshape(o.Repl.PXZ[:,ind1], :, 1)  # store as matrix to reduce compiler confusion
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
			T₁ = o.Repl.Yendog[ind1+1] ? o.Repl.S⋂YX[ind2+1] * o.invXXS✻XU[ind1+1] : Matrix{T}(undef,0,0)  #  S_∩ (Y_(∥j):*X_∥ ) (X_∥^' X_∥ )^(-1) [S_* (U ̈_(∥i):*X_∥ )]^'

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
				ncols(o.Repl.Zperp) > 0 && (T₁ = T₁ .- o.Repl.S⋂PXYZperp[ind1+1] * o.invZperpZperpS✻ZperpU[ind2+1])  # subtract S_∩ (P_(X_∥ ) Y_(∥i):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) [S_* (U ̈_(∥j):*Z_⊥ )]^'
				o.NFE               > 0 && (T₁ = T₁ .- o.Repl.CT_FE⋂PY[ind1+1] * o.CTFEU[ind2+1])
			end

			if iszero(ind2)  # order-0 and -1 terms
				if iszero(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1]
				elseif isone(ncols(T₁))
					dest = o.Repl.FillingT₀[ind1+1,1] .+ T₁ .* β̈v
				else
					dest = T₁ * o.v; dest .+= o.Repl.FillingT₀[ind1+1,1]
				end
			else  # y component
				if iszero(ncols(T₁))  # - x*β̈ components
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̈
				elseif isone(ncols(T₁))
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̈ .+ T₁ .* β̈v
				else
					dest .+= o.Repl.FillingT₀[ind1+1,ind2+1] * _β̈
					o.matmulplus!(dest, T₁, β̈v)
				end
			end

			if o.Repl.Yendog[ind1+1] && o.Repl.Yendog[ind2+1]  # (S_*:*U ̈_(∥j)^'-S_* (U ̈_(∥j):*Z_⊥ ) (Z_⊥^' Z_⊥ )^(-1) Z_⊥^' )^' v^*
				for i ∈ 1:o.clust[1].N  # put this logic inside colquadform(...Array{T,3}...)
					o.colquadformminus!(dest, i, o.v, view(o.S✻UPX_S✻UMZperp[ind1+1,ind2+1] ,:,i,:), β̈v)
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
				J⋂s = Filling(o, 1, _β̈s) ./ o.As
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
			
			@inbounds Threads.@threads for b ∈ axes(o.v,2)
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
