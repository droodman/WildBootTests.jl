# Logically, wild bootstrap tests perform estimation at two stages, once as part of the bootstrap DGP, once in each bootstrap replication
# The StrEstimator "class" and its three "children" hold the estimation logic for the OLS, Anderson-Rubin, and IV/GMM cases

function perp(A::AbstractMatrix)
  F = eigen(Symmetric(A*invsym(A'A)*A'))
  F.vectors[:, abs.(F.values) .< 1000eps(eltype(A))]
end

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 nrows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T}, parent::StrBootTest{T}, R₁::AbstractMatrix{T}, R::Union{UniformScaling{Bool},AbstractMatrix{T}}=Matrix{T}(undef,0,0)) where T
  o.restricted = nrows(R₁) > 0
	if o.restricted
		singular, invR₁R₁ = invsymsingcheck(R₁ * R₁')
		singular && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
	  o.R₁invR₁R₁ = R₁'invR₁R₁
	  F = eigen(Symmetric(o.R₁invR₁R₁ * R₁))
	  o.R₁perp = F.vectors[:, abs.(F.values) .< 1000*eps(T)]  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
  else
	  o.R₁invR₁R₁ = Matrix{T}(undef, parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
	  RR₁perp = Matrix{T}([R ; zeros(T, parent.kY₂, parent.kX₁) I])  # rows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

	  o.restricted && (RR₁perp *= o.R₁perp)

	  F = eigen(Symmetric(RR₁perp'pinv(RR₁perp * RR₁perp')*RR₁perp)); val = abs.(F.values) .> 1000*eps(T)
	  o.Rpar   = F.vectors[:, val]  # defines attack surface for null
	  o.restricted && (o.Rpar = o.R₁perp * o.Rpar)  # fold model constraint factors into Rpar, RperpX
	  o.RRpar = R * o.Rpar
		o.RparX = o.Rpar[1:parent.kX₁,:]  # part of Rpar that refers to Y₂
		o.RparY = o.Rpar[parent.kX₁+1:end,:]  # part of Rpar that refers to Y₂
	  o.R₁invR₁R₁Y = o.R₁invR₁R₁[parent.kX₁+1:end,:]
	  o.RR₁invR₁R₁ = R * o.R₁invR₁R₁

		if o.isDGP || !parent.WREnonARubin  # DGP regressions will just copy the replication-regression Zperp, X₂par for speed
			o.RperpX = F.vectors[:, abs.(F.values) .< 1000*eps(T)]
			o.restricted && (o.RperpX = o.R₁perp * o.RperpX)  # fold model constraint factors into Rpar, RperpX
			o.RperpX = o.RperpX[1:parent.kX₁,:]  # Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X₁, since all endogenous ones will be retained
			o.RperpXperp = perp(o.RperpX)
		end
  end
	nothing
end

# stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
function InitVarsOLS!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}) where T # Rperp is for replication regression--no null imposed
  o.y₁par = parent.y₁
	o.ü₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk]

	if parent.granular
		X₁y₁ = parent.X₁'parent.y₁
		H    = parent.X₁'parent.X₁
	else	
		S✻X₁y₁ = panelcross(parent.X₁, parent.y₁, parent.info✻)
		X₁y₁ = vec(sumpanelcross(S✻X₁y₁))
		o.S✻XX = panelcross(parent.X₁, parent.X₁, parent.info✻)
	  H = sumpanelcross(o.S✻XX)
	end

	o.invH = #=Symmetric=#(pinv(H))
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.invH : #=Symmetric=#(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
	o.β̈₀ = R₁AR₁ * X₁y₁
	o.∂β̈∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

	if parent.jk
		if parent.purerobust
			o.invMjkv = rowquadform(o.invH, parent.X₁)
			o.invMjkv .= 1 ./ (1 .- o.invMjkv)  # standard hc3 multipliers
		elseif parent.granularjk
			o.invMjk = Vector{Matrix{T}}(undef, parent.N✻)
			for g ∈ 1:parent.N✻
				S = parent.info✻[g]
				v = view(parent.X₁, S,:)
				o.invMjk[g] = - v * R₁AR₁ * v'
				o.invMjk[g][1:length(S)+1:length(S)^2] .+= one(T)  # add I
				o.invMjk[g] .= invsym(o.invMjk[g])
			end
		else
			!isdefined(o, :S✻XX) && (o.S✻XX = panelcross(parent.X₁, parent.X₁, parent.info✻))
			_H = reshape(H, (parent.kX, 1, parent.kX)) .- o.S✻XX
			_invH = iszero(nrows(o.R₁perp)) ? invsym(_H) : o.R₁perp * invsym(o.R₁perp' * _H * o.R₁perp) * o.R₁perp'
			o.XinvHjk = [view(parent.X₁, parent.info✻[g],:) * view(_invH,:,g,:) for g ∈ 1:parent.N✻]
		end
	end

  o.A = iszero(nrows(Rperp)) ? o.invH : #=Symmetric=#(Rperp * invsym(Rperp'H*Rperp) * Rperp')  # for replication regression
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = parent.X₁ * o.AR)
	nothing
end

function InitVarsARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
	o.y₁par = Vector{T}(undef, parent.Nobs)
	o.ü₁    = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk]

	X₁y₁ = parent.X₁'parent.y₁
	X₂y₁ = parent.X₂'parent.y₁
	X₁Y₂ = parent.X₁'parent.Y₂
	X₂Y₂ = parent.X₂'parent.Y₂

	if !parent.jk || !(parent.granularjk || parent.purerobust)
		X₂X₁ = parent.X₂'parent.X₁
		X₁X₁ = parent.X₁'parent.X₁
		X₂X₂ = parent.X₂'parent.X₂
		if !(parent.granularjk || parent.purerobust)
			S✻X₁X₁ = panelcross(parent.X₁, parent.X₁, parent.info✻)
			S✻X₂X₁ = panelcross(parent.X₂, parent.X₁, parent.info✻)
			S✻X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻)
		end
	else
		S✻X₂X₁ = panelcross(parent.X₂, parent.X₁, parent.info✻)
		X₂X₁ = sumpanelcross(S✻X₂X₁)
		S✻X₁X₁ = panelcross(parent.X₁, parent.X₁, parent.info✻)
		X₁X₁ = sumpanelcross(S✻X₁X₁)
		S✻X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻)
		X₂X₂ = sumpanelcross(S✻X₂X₂)
	end

  H = #=Symmetric=#([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
  o.A = invsym(H)
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.A : #=Symmetric=#(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')
	o.β̈₀   = R₁AR₁ * [X₁y₁ ; X₂y₁]
	o.∂β̈∂r = R₁AR₁ * [X₁Y₂ ; X₂Y₂]

	if parent.jk
		if parent.purerobust
			o.invMjkv =     rowquadform(parent.X₁, (@view o.invH[1:parent.kX₁,     1:parent.kX₁    ]), parent.X₁) +
			            2 * rowquadform(parent.X₁, (@view o.invH[1:parent.kX₁,     parent.kX₁+1:end]), parent.X₂) +
									    rowquadform(parent.X₂, (@view o.invH[parent.kX₁+1:end, parent.kX₁+1:end]), parent.X₂)
			o.invMjkv .= 1 ./ (1 .- o.invMjkv)  # standard hc3 multipliers
		elseif parent.granularjk
			o.invMjk = Vector{Matrix{T}}(undef, parent.N✻)
			negR₁AR₁ = -R₁AR₁
			for g ∈ 1:parent.N✻
				S = parent.info✻[g] 
				v₁ = view(parent.X₁, S, :); v₂ = view(parent.X₂, S, :)
				o.invMjk[g] = v₁ * (@view negR₁AR₁[1:parent.kX₁    , 1:parent.kX₁    ]) * v₁' +
											v₁ * (@view negR₁AR₁[1:parent.kX₁    , parent.kX₁+1:end]) * v₂' +
											v₂ * (@view negR₁AR₁[parent.kX₁+1:end, 1:parent.kX₁    ]) * v₁' +
											v₂ * (@view negR₁AR₁[parent.kX₁+1:end, parent.kX₁+1:end]) * v₂'
				o.invMjk[g][1:length(S)+1:length(S)^2] .+= one(T)  # add I
				o.invMjk[g] .= invsym(o.invMjk[g])
			end
		else
			!isdefined(o, :S✻XX) && (o.S✻XX = [[S✻X₁X₁ S✻X₂X₁'] ; [S✻X₂X₁ S✻X₂X₂]])
			_H = reshape(H, (parent.kX, 1, parent.kX)) .- o.S✻XX
			_invH = iszero(nrows(o.R₁perp)) ? invsym(_H) : o.R₁perp * invsym(o.R₁perp' * _H * o.R₁perp) * o.R₁perp'
			o.XinvHjk = [(S = parent.info✻[g]; X₁₂B(view(parent.X₁, S,:), view(parent.X₂, S,:), view(_invH,:,g,:))) for g ∈ 1:parent.N✻]
		end
	end

  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = X₁₂B(parent.X₁, parent.X₂, o.AR))
	nothing
end

function InitVarsIV!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}...) where T
	!isempty(Rperp) && (o.Rperp = Rperp[1])
	o.kZ = ncols(o.Rpar)

	if parent.granular
		if !o.isDGP && parent.WREnonARubin
			o.X₁ = parent.DGP.X₁
			o.Zperp = parent.DGP.Zperp
			o.invZperpZperp = parent.DGP.invZperpZperp
		else
			o.kZperp = ncols(o.RperpX)
			o.Zperp = parent.X₁ * o.RperpX
			o.invZperpZperp = iszero(length(o.Zperp)) ? Matrix{T}(undef,0,0) : invsym(o.Zperp'o.Zperp)
			o.X₁ = parent.X₁ * o.RperpXperp
			o.kX = ncols(o.X₁) + parent.kX₂
		end
	
		o.Z   = X₁₂B(parent.X₁, parent.Y₂, o.Rpar     )  # Z∥
		o.ZR₁ = X₁₂B(parent.X₁, parent.Y₂, o.R₁invR₁R₁)
		ZperpZ   = o.Zperp'o.Z
		ZperpZR₁ = o.Zperp'o.ZR₁

		if parent.jk && o.isDGP && parent.WREnonARubin
			o.β̈ⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻,1)
			o.YYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
			o.invXXXy₁parⱼₖ = Array{T,3}(undef, o.kX, parent.N✻, 1)
			o.ZXinvXXXy₁parⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻, 1)
			if o.liml
				o.κⱼₖ = Array{T,3}(undef, 1, parent.N✻, 1)
				o.YPXYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
			end

			ZperpX₁   , _ZperpX₁    = crossjk(o.Zperp, o.X₁     , parent.info✻)  # full-sample and delete-g cross-moments
			ZperpX₂   , _ZperpX₂    = crossjk(o.Zperp, parent.X₂, parent.info✻)
			Zperpy₁   , _Zperpy₁    = crossjk(o.Zperp, parent.y₁, parent.info✻)
			ZperpY₂   , _ZperpY₂    = crossjk(o.Zperp, parent.Y₂, parent.info✻)
			ZperpZ    , _ZperpZ     = crossjk(o.Zperp, o.Z      , parent.info✻)
			ZperpZR₁  , _ZperpZR₁   = crossjk(o.Zperp, o.ZR₁    , parent.info✻)
			  ZperpZperp, _ZperpZperp = crossjk(o.Zperp, o.Zperp, parent.info✻)
			_invZperpZperp = invsym(_ZperpZperp)  # inverses of all delete-g Zperp'Zperp

      _invZperpZperpZperpX₁  = _invZperpZperp * _ZperpX₁
      _invZperpZperpZperpX₂  = _invZperpZperp * _ZperpX₂
      _invZperpZperpZperpZ   = _invZperpZperp * _ZperpZ
      _invZperpZperpZperpZR₁ = _invZperpZperp * _ZperpZR₁
      _invZperpZperpZperpY₂  = _invZperpZperp * _ZperpY₂
      _invZperpZperpZperpy₁  = _invZperpZperp * _Zperpy₁

			o.X₁ⱼₖ  = partialjk(     o.X₁ , o.Zperp, _invZperpZperpZperpX₁ , parent.info✻)    # FWL-process
			o.X₂ⱼₖ  = partialjk(parent.X₂ , o.Zperp, _invZperpZperpZperpX₂ , parent.info✻)
			o.y₁ⱼₖ  = partialjk(parent.y₁ , o.Zperp, _invZperpZperpZperpy₁ , parent.info✻)
			o.Y₂ⱼₖ  = partialjk(parent.Y₂ , o.Zperp, _invZperpZperpZperpY₂ , parent.info✻)
			o.Zⱼₖ   = partialjk(     o.Z  , o.Zperp, _invZperpZperpZperpZ  , parent.info✻)
			o.ZR₁ⱼₖ = partialjk(     o.ZR₁, o.Zperp, _invZperpZperpZperpZR₁, parent.info✻)

			tX₁ = ZperpX₁ - ZperpZperp * _invZperpZperpZperpX₁
      tX₂ = ZperpX₂ - ZperpZperp * _invZperpZperpZperpX₂
      tZ  = ZperpZ  - ZperpZperp * _invZperpZperpZperpZ 
      tY₂ = ZperpY₂ - ZperpZperp * _invZperpZperpZperpY₂
      ty₁ = Zperpy₁ .- ZperpZperp * _invZperpZperpZperpy₁

			_y₁ = view(parent.y₁,:,:)
			_X₂X₁   = parent.X₂'o.X₁      - ZperpX₂'_invZperpZperpZperpX₁ - _invZperpZperpZperpX₂'tX₁ - panelcross(o.X₂ⱼₖ, o.X₁ⱼₖ, parent.info✻)
      _X₁X₁   = o.X₁'o.X₁           - ZperpX₁'_invZperpZperpZperpX₁ - _invZperpZperpZperpX₁'tX₁ - panelcross(o.X₁ⱼₖ, o.X₁ⱼₖ, parent.info✻)
      _X₂X₂   = parent.X₂'parent.X₂ - ZperpX₂'_invZperpZperpZperpX₂ - _invZperpZperpZperpX₂'tX₂ - panelcross(o.X₂ⱼₖ, o.X₂ⱼₖ, parent.info✻)
      _X₁Y₂   = o.X₁'parent.Y₂      - ZperpX₁'_invZperpZperpZperpY₂ - _invZperpZperpZperpX₁'tY₂ - panelcross(o.X₁ⱼₖ, o.Y₂ⱼₖ, parent.info✻)
      _X₂Y₂   = parent.X₂'parent.Y₂ - ZperpX₂'_invZperpZperpZperpY₂ - _invZperpZperpZperpX₂'tY₂ - panelcross(o.X₂ⱼₖ, o.Y₂ⱼₖ, parent.info✻)
      o.Y₂y₁ⱼₖ = parent.Y₂'_y₁ - ZperpY₂'_invZperpZperpZperpy₁ - _invZperpZperpZperpY₂'ty₁ - panelcross(o.Y₂ⱼₖ, o.y₁ⱼₖ, parent.info✻)
      o.X₂y₁ⱼₖ = parent.X₂'_y₁ - ZperpX₂'_invZperpZperpZperpy₁ - _invZperpZperpZperpX₂'ty₁ - panelcross(o.X₂ⱼₖ, o.y₁ⱼₖ, parent.info✻)
      o.X₁y₁ⱼₖ = o.X₁'_y₁      - ZperpX₁'_invZperpZperpZperpy₁ - _invZperpZperpZperpX₁'ty₁ - panelcross(o.X₁ⱼₖ, o.y₁ⱼₖ, parent.info✻)
      o.Zy₁ⱼₖ  = o.Z'_y₁       - ZperpZ'_invZperpZperpZperpy₁  - _invZperpZperpZperpZ'ty₁  - panelcross(o.Zⱼₖ,  o.y₁ⱼₖ, parent.info✻)
      o.XZⱼₖ   = [o.X₁'o.Z           - ZperpX₁'_invZperpZperpZperpZ  - _invZperpZperpZperpX₁'tZ  - panelcross(o.X₁ⱼₖ, o.Zⱼₖ,  parent.info✻)
                  parent.X₂'o.Z      - ZperpX₂'_invZperpZperpZperpZ  - _invZperpZperpZperpX₂'tZ  - panelcross(o.X₂ⱼₖ, o.Zⱼₖ, parent.info✻)]
      o.ZZⱼₖ   = o.Z'o.Z             - ZperpZ'_invZperpZperpZperpZ   - _invZperpZperpZperpZ'tZ   - panelcross(o.Zⱼₖ,  o.Zⱼₖ,  parent.info✻)
      o.ZY₂ⱼₖ  = o.Z'parent.Y₂       - ZperpZ'_invZperpZperpZperpY₂  - _invZperpZperpZperpZ'tY₂  - panelcross(o.Zⱼₖ,  o.Y₂ⱼₖ, parent.info✻)
      o.y₁y₁ⱼₖ = _y₁'_y₁ -2 * (Zperpy₁'_invZperpZperpZperpy₁) + _invZperpZperpZperpy₁'ZperpZperp*_invZperpZperpZperpy₁ - panelcross(o.y₁ⱼₖ, o.y₁ⱼₖ, parent.info✻)

      o.XY₂ⱼₖ = [_X₁Y₂ ; _X₂Y₂]
      o.XXⱼₖ = [_X₁X₁ ; _X₂X₁ ;;; _X₂X₁' ;  _X₂X₂]
			o.invXXⱼₖ = invsym(o.XXⱼₖ)
      o.H_2SLSⱼₖ = o.XZⱼₖ'o.invXXⱼₖ * o.XZⱼₖ
      (!isone(o.κ) || o.liml) && (o.H_2SLSmZZⱼₖ = o.H_2SLSⱼₖ - o.ZZⱼₖ)

			if o.restricted
        tZR₁ = ZperpZR₁ - ZperpZperp * _invZperpZperpZperpZR₁
        o.X₁ZR₁ⱼₖ    =  o.X₁'o.ZR₁      - ZperpX₁'_invZperpZperpZperpZR₁  - _invZperpZperpZperpX₁'tZR₁  - panelcross(o.X₁ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻)
        o.X₂ZR₁ⱼₖ    =  parent.X₂'o.ZR₁ - ZperpX₂'_invZperpZperpZperpZR₁  - _invZperpZperpZperpX₂'tZR₁  - panelcross(o.X₂ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻)
        o.ZZR₁ⱼₖ     =  o.Z'o.ZR₁       - ZperpZ'_invZperpZperpZperpZR₁   - _invZperpZperpZperpZ'tZR₁   - panelcross(o.Zⱼₖ,   o.ZR₁ⱼₖ, parent.info✻)
        o.twoZR₁y₁ⱼₖ =2(o.ZR₁'_y₁       - ZperpZR₁'_invZperpZperpZperpy₁  - _invZperpZperpZperpZR₁'ty₁  - panelcross(o.ZR₁ⱼₖ,  o.y₁ⱼₖ, parent.info✻))
        o.ZR₁ZR₁ⱼₖ   =  o.ZR₁'o.ZR₁     - ZperpZR₁'_invZperpZperpZperpZR₁ - _invZperpZperpZperpZR₁'tZR₁ - panelcross(o.ZR₁ⱼₖ, o.ZR₁ⱼₖ, parent.info✻)
        o.ZR₁Y₂ⱼₖ    =  o.ZR₁'parent.Y₂ - ZperpZR₁'_invZperpZperpZperpY₂  - _invZperpZperpZperpZR₁'tY₂  - panelcross(o.ZR₁ⱼₖ, o.Y₂ⱼₖ,  parent.info✻)
			else
				o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ 
				o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ
				o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ 
				o.Xy₁parⱼₖ     = [o.X₁y₁ⱼₖ ; o.X₂y₁ⱼₖ]
				o.y₁parⱼₖ      = o.y₁ⱼₖ
			end

			!iszero(o.fuller) &&
				(o.Nobsⱼₖ = parent._Nobs .- (parent.fweights ? @panelsum(parent.wt, parent.info✻) : length.(parent.info✻)))
		else
			ZperpX₁    = o.Zperp'o.X₁
			ZperpX₂    = o.Zperp'parent.X₂
			Zperpy₁    = o.Zperp'parent.y₁
			ZperpY₂    = o.Zperp'parent.Y₂
			ZperpZ     = o.Zperp'o.Z
			ZperpZR₁   = o.Zperp'o.ZR₁
			ZperpZperp = o.Zperp'o.Zperp
		end

		if !o.isDGP && parent.WREnonARubin
			o.X₂ = parent.DGP.X₂
			o.y₁ = parent.DGP.y₁
			o.Y₂ = parent.DGP.Y₂
			o.XX = parent.DGP.XX
			o.invXX = parent.DGP.invXX
			o.XY₂ = parent.DGP.XY₂
			o.Y₂y₁ = parent.DGP.Y₂y₁
			o.X₂y₁ = parent.DGP.X₂y₁
			o.X₁y₁ = parent.DGP.X₁y₁
			o.y₁y₁ = parent.DGP.y₁y₁
		else
			o.X₁ .-= o.Zperp * (o.invZperpZperp * ZperpX₁)  # FWL-process X₁
			o.X₂ = o.Zperp * (o.invZperpZperp * ZperpX₂); o.X₂ .= parent.X₂ .- o.X₂   # FWL-process X₂
			o.y₁ = parent.y₁ - o.Zperp * (o.invZperpZperp * Zperpy₁)
			o.Y₂ = parent.Y₂ - o.Zperp * (o.invZperpZperp * ZperpY₂)
		
			X₂X₁ = o.X₂'o.X₁
			o.XX = Symmetric([o.X₁'o.X₁ X₂X₁' ; X₂X₁ o.X₂'o.X₂])
			o.invXX = invsym(o.XX)
			X₁Y₂ = o.X₁'o.Y₂
			X₂Y₂ = o.X₂'o.Y₂
			o.XY₂ = [X₁Y₂ ; X₂Y₂]
			o.Y₂y₁ = o.Y₂'o.y₁
			o.X₂y₁ = o.X₂'o.y₁
			o.X₁y₁ = o.X₁'o.y₁
			o.y₁y₁ = dot(o.y₁, o.y₁)
		end

		o.Z   .-= o.Zperp * (o.invZperpZperp * ZperpZ  )  # partialling out
		o.ZR₁ .-= o.Zperp * (o.invZperpZperp * ZperpZR₁)

		o.Zy₁ = o.Z'o.y₁
		o.ZY₂ = o.Z'o.Y₂
		o.ZZ  = Symmetric(o.Z'o.Z)
		o.XZ  = [o.X₁'o.Z ; o.X₂'o.Z]
	
		if o.restricted
			o.X₂ZR₁    = o.X₂'o.ZR₁
			o.X₁ZR₁    = o.X₁'o.ZR₁
			o.ZR₁Z     = o.ZR₁'o.Z
			o.twoZR₁y₁ = 2*(o.ZR₁'o.y₁)
			o.ZR₁ZR₁   = Symmetric(o.ZR₁'o.ZR₁)
			o.ZR₁Y₂    = o.ZR₁'o.Y₂
		else
			o.Y₂y₁par    = o.Y₂y₁
			o.X₂y₁par    = o.X₂y₁
			o.X₁y₁par    = o.X₁y₁
			o.Zy₁par     = o.Zy₁
			o.y₁pary₁par = o.y₁y₁
			o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
		end

		o.y₁par = copy(o.y₁)

		if o.isDGP && !parent.scorebs
			o.Ü₂ = [Matrix{T}(undef, parent.Nobs, parent.kY₂) for _ in 0:parent.jk]
			o.u⃛₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk]
		end

	else  # coarse

		if !o.isDGP && parent.WREnonARubin
			o.Zperp = parent.DGP.Zperp
			o.Xpar₁ = parent.DGP.Xpar₁
			isdefined(parent.DGP, :X₁) && (o.X₁ = parent.DGP.X₁)
			isdefined(parent.DGP, :X₂) && (o.X₂ = parent.DGP.X₂)
			isdefined(parent.DGP, :Y₂) && (o.Y₂ = parent.DGP.Y₂)
			isdefined(parent.DGP, :y₁) && (o.y₁ = parent.DGP.y₁)
			o.invZperpZperp = parent.DGP.invZperpZperp
			o.XX = parent.DGP.XX
			o.invXX = parent.DGP.invXX
			o.XY₂ = parent.DGP.XY₂
			o.Y₂Y₂ = parent.DGP.Y₂Y₂
			o.S✻Y₂Y₂ = parent.DGP.S✻Y₂Y₂
			o.Y₂y₁ = parent.DGP.Y₂y₁
			o.X₂y₁ = parent.DGP.X₂y₁
			o.X₁y₁ = parent.DGP.X₁y₁
			o.y₁y₁ = parent.DGP.y₁y₁
			o.invZperpZperpZperpX₁ = parent.DGP.invZperpZperpZperpX₁
			o.invZperpZperpZperpX₂ = parent.DGP.invZperpZperpZperpX₂
			o.invZperpZperpZperpy₁ = parent.DGP.invZperpZperpZperpy₁
			o.invZperpZperpZperpY₂ = parent.DGP.invZperpZperpZperpY₂
			o.S✻⋂X₁Y₂ = parent.DGP.S✻⋂X₁Y₂
			o.S✻⋂X₂Y₂ = parent.DGP.S✻⋂X₂Y₂
			o.S✻⋂ZperpY₂ = parent.DGP.S✻⋂ZperpY₂
			o.S✻Y₂y₁ = parent.DGP.S✻Y₂y₁
		else
			o.kZperp = ncols(o.RperpX)
			o.Zperp = t✻(parent.X₁, o.RperpX)
			o.S✻⋂ZperpZperp = panelcross(o.Zperp, o.Zperp, parent.info✻⋂)
			o.ZperpZperp = iszero(ncols(o.RperpX)) ? #=Symmetric=#(Matrix{T}(undef,0,0)) : sumpanelcross(o.S✻⋂ZperpZperp)
			o.invZperpZperp = invsym(o.ZperpZperp)

			o.Xpar₁ = parent.X₁ * o.RperpXperp  # X∥ := [Xpar₁ X₂]
			S✻⋂X₁Zperp = panelcross(o.Xpar₁, o.Zperp, parent.info✻⋂)
			ZperpX₁ = sumpanelcross(S✻⋂X₁Zperp)'
			o.invZperpZperpZperpX₁ = o.invZperpZperp * ZperpX₁

			S✻⋂X₂Zperp = panelcross(parent.X₂, o.Zperp, parent.info✻⋂)
			ZperpX₂ = sumpanelcross(S✻⋂X₂Zperp)'
			o.ZperpX = [ZperpX₁ ZperpX₂]
			o.S✻⋂XZperp = [S✻⋂X₁Zperp; S✻⋂X₂Zperp]

			o.invZperpZperpZperpX₂ = o.invZperpZperp * ZperpX₂
			o.invZperpZperpZperpX = [o.invZperpZperpZperpX₁ o.invZperpZperpZperpX₂]

			if parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt) || !o.liml && !isempty(Rperp)
				o.X₁ = o.Xpar₁ - o.Zperp * o.invZperpZperpZperpX₁  # shrink and FWL-process X₁; do it as an O(N) operation because it can be so size-reducing
				o.X₂ = o.Zperp * o.invZperpZperpZperpX₂; o.X₂ .= parent.X₂ .- o.X₂  # FWL-process X₂
			end
			
			S✻⋂X₁X₁ = panelcross(o.Xpar₁, o.Xpar₁, parent.info✻⋂)
			S✻⋂X₂X₁ = panelcross(parent.X₂, o.Xpar₁, parent.info✻⋂)
			S✻⋂X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻⋂)
			X₂X₁ = sumpanelcross(S✻⋂X₂X₁)
			X₁X₁ = #=Symmetric=#(sumpanelcross(S✻⋂X₁X₁))
			X₂X₂ = #=Symmetric=#(sumpanelcross(S✻⋂X₂X₂))
			o.S✻⋂XX = [[S✻⋂X₁X₁ S✻⋂X₂X₁'] ; [S✻⋂X₂X₁ S✻⋂X₂X₂]]  # [a b; c d] syntax would call hvcat() to concatenate horizontally along dim 2 rather than 3

			X₂X₁ -= ZperpX₂'o.invZperpZperp * ZperpX₁
			X₁X₁ -= #=Symmetric=#(ZperpX₁'o.invZperpZperp * ZperpX₁)
			X₂X₂ -= #=Symmetric=#(ZperpX₂'o.invZperpZperp * ZperpX₂)
			o.XX = #=Symmetric=#([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
			o.kX = ncols(o.XX)
			o.invXX = invsym(o.XX)

			o.S✻⋂ZperpY₂ = panelcross(o.Zperp, parent.Y₂, parent.info✻⋂)
			ZperpY₂ = sumpanelcross(o.S✻⋂ZperpY₂)
			o.invZperpZperpZperpY₂ = o.invZperpZperp * ZperpY₂
			(parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt)) &&
				(o.Y₂ = parent.Y₂ - o.Zperp * o.invZperpZperpZperpY₂)
			o.S✻⋂Zperpy₁ = panelcross(o.Zperp, parent.y₁, parent.info✻⋂)
			Zperpy₁ = vec(sumpanelcross(o.S✻⋂Zperpy₁))
			o.invZperpZperpZperpy₁ = o.invZperpZperp * Zperpy₁
			((parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt)) || parent.scorebs) &&
				(o.y₁ = parent.y₁ - o.Zperp * o.invZperpZperpZperpy₁)
		  o.S✻⋂X₁Y₂ = panelcross(o.Xpar₁, parent.Y₂, parent.info✻⋂)
		  o.S✻⋂X₂Y₂ = panelcross(parent.X₂, parent.Y₂, parent.info✻⋂)
			o.S✻⋂XY₂ = [o.S✻⋂X₁Y₂; o.S✻⋂X₂Y₂]

			o.XY₂ = sumpanelcross(o.S✻⋂XY₂)
			o.S✻Y₂y₁ = panelcross(parent.Y₂, parent.y₁, parent.info✻)
		  o.Y₂y₁ = vec(sumpanelcross(o.S✻Y₂y₁))
			o.S✻Y₂Y₂ = panelcross(parent.Y₂, parent.Y₂, parent.info✻)
			o.Y₂Y₂ = #=Symmetric=#(sumpanelcross(o.S✻Y₂Y₂))
			S✻⋂X₂y₁ = panelcross(parent.X₂, parent.y₁, parent.info✻⋂)
			o.X₂y₁ = vec(sumpanelcross(S✻⋂X₂y₁))
			S✻⋂X₁y₁ = panelcross(o.Xpar₁, parent.y₁, parent.info✻⋂)
			o.S✻⋂Xy₁ = [S✻⋂X₁y₁; S✻⋂X₂y₁]
			o.X₁y₁ = vec(sumpanelcross(S✻⋂X₁y₁))
			o.S✻y₁y₁ = panelcross(parent.y₁, parent.y₁, parent.info✻)
			o.y₁y₁ = sum(o.S✻y₁y₁)
			o.XY₂ -= o.invZperpZperpZperpX'ZperpY₂
		  o.Y₂y₁ -= ZperpY₂'o.invZperpZperpZperpy₁
			o.Y₂Y₂ -= #=Symmetric=#(ZperpY₂'o.invZperpZperpZperpY₂)
			o.X₂y₁ -= ZperpX₂'o.invZperpZperpZperpy₁
			o.X₁y₁ -= ZperpX₁'o.invZperpZperpZperpy₁
			o.y₁y₁ -= Zperpy₁'o.invZperpZperpZperpy₁
		end

	  o.Z = X₁₂B(parent.X₁, parent.Y₂, o.Rpar)  # Z∥

		X₁par = parent.X₁ * o.RparX  # XXX expressible as a linear combination of Xpar₁??
		S✻⋂X₁Zpar = panelcross(o.Xpar₁  , X₁par, parent.info✻⋂) + o.S✻⋂X₁Y₂ * o.RparY
		S✻⋂X₂Zpar = panelcross(parent.X₂, X₁par, parent.info✻⋂) + o.S✻⋂X₂Y₂ * o.RparY
		o.S✻⋂XZpar = [S✻⋂X₁Zpar; S✻⋂X₂Zpar]
		X₁Zpar = sumpanelcross(S✻⋂X₁Zpar)
		X₂Zpar = sumpanelcross(S✻⋂X₂Zpar)
		S✻⋂ZperpX₁par = panelcross(o.Zperp, X₁par, parent.info✻⋂)
		ZperpX₁par = sumpanelcross(S✻⋂ZperpX₁par)
		o.S✻⋂ZperpZpar = S✻⋂ZperpX₁par + o.S✻⋂ZperpY₂ * o.RparY
	  ZperpZpar = ZperpX₁par + sumpanelcross(o.S✻⋂ZperpY₂) * o.RparY
		o.invZperpZperpZperpZpar = o.invZperpZperp * ZperpZpar

		S✻X₁pary₁ = panelcross(X₁par, parent.y₁, parent.info✻)
		o.S✻Zpary₁ = S✻X₁pary₁ + o.RparY' * o.S✻Y₂y₁
		o.Zy₁ = vec(sumpanelcross(S✻X₁pary₁)) + o.RparY' * o.Y₂y₁

		S✻X₁parY₂ = panelcross(X₁par, parent.Y₂, parent.info✻)
		o.XZ = [X₁Zpar - o.invZperpZperpZperpX₁'ZperpZpar ; X₂Zpar - o.invZperpZperpZperpX₂'ZperpZpar]
	  o.S✻ZparY₂ = S✻X₁parY₂ + o.RparY' * o.S✻Y₂Y₂
	  tmp = S✻X₁parY₂ * o.RparY; o.S✻ZparZpar = panelcross(X₁par, X₁par, parent.info✻) + tmp + tmp' + o.RparY' * o.S✻Y₂Y₂ * o.RparY
	  o.ZZ = #=Symmetric=#(sumpanelcross(o.S✻ZparZpar)); o.ZZ -= #=Symmetric=#(ZperpZpar'o.invZperpZperpZperpZpar)
		o.Zy₁ -= ZperpX₁par'o.invZperpZperpZperpy₁

	  if o.restricted
			o.ZR₁ = X₁₂B(parent.X₁, parent.Y₂, o.R₁invR₁R₁)
			S✻⋂X₁ZR₁ = panelcross(o.Xpar₁, o.ZR₁, parent.info✻⋂)
			S✻⋂X₂ZR₁ = panelcross(parent.X₂, o.ZR₁, parent.info✻⋂)
			o.S✻⋂XZR₁ = [S✻⋂X₁ZR₁ ; S✻⋂X₂ZR₁]
			o.S✻⋂ZperpZR₁ = panelcross(o.Zperp, o.ZR₁, parent.info✻⋂)
			o.ZperpZR₁ = sumpanelcross(o.S✻⋂ZperpZR₁)
			o.invZperpZperpZperpZR₁ = o.invZperpZperp * o.ZperpZR₁
		  o.X₁ZR₁    = sumpanelcross(S✻⋂X₁ZR₁)
		  o.X₂ZR₁    = sumpanelcross(S✻⋂X₂ZR₁)
			o.S✻ZR₁Z   = panelcross(o.ZR₁, o.Z, parent.info✻)
		  o.ZR₁Z     = sumpanelcross(o.S✻ZR₁Z)
			o.S✻ZR₁Y₂  = panelcross(o.ZR₁, parent.Y₂, parent.info✻)
		  o.ZR₁Y₂    = sumpanelcross(o.S✻ZR₁Y₂) 
			o.S✻ZR₁y₁  = panelcross(o.ZR₁, parent.y₁, parent.info✻)
		  o.twoZR₁y₁ = 2 * vec(sumpanelcross(o.S✻ZR₁y₁))
			o.S✻ZR₁ZR₁ = panelcross(o.ZR₁, o.ZR₁, parent.info✻)
		  o.ZR₁ZR₁   = #=Symmetric=#(sumpanelcross(o.S✻ZR₁ZR₁))
			o.ZR₁ -= o.Zperp * o.invZperpZperpZperpZR₁
		  t✻minus!(o.X₁ZR₁, o.invZperpZperpZperpX₁', o.ZperpZR₁)
		  t✻minus!(o.X₂ZR₁, o.invZperpZperpZperpX₂', o.ZperpZR₁)
		  o.ZR₁Z    -= o.ZperpZR₁'o.invZperpZperp * ZperpZpar
		  t✻minus!(o.ZR₁Y₂, o.ZperpZR₁', o.invZperpZperpZperpY₂)
		  o.twoZR₁y₁-= 2 * o.ZperpZR₁'o.invZperpZperpZperpy₁
		  o.ZR₁ZR₁  -= #=Symmetric=#(o.ZperpZR₁'o.invZperpZperp * o.ZperpZR₁)
		else
		  o.Y₂y₁par    = o.Y₂y₁
		  o.X₂y₁par    = o.X₂y₁
		  o.X₁y₁par    = o.X₁y₁
		  o.Zy₁par     = o.Zy₁
		  o.y₁pary₁par = o.y₁y₁
		  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
		end

		parent.scorebs && (o.y₁par = copy(o.y₁))

		o.isDGP && parent.WREnonARubin && 
			(o.ZY₂ = sumpanelcross(o.S✻ZparY₂) - ZperpZpar'o.invZperpZperpZperpY₂)
		
		t✻minus!(o.Z, o.Zperp, o.invZperpZperpZperpZpar)

	end  # end coarse

	o.V =  o.invXX * o.XZ  # in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and liml robust VCV by Stata convention
	o.H_2SLS = o.XZ'o.V  # Hessian

	if o.isDGP
		if parent.jk && parent.WREnonARubin
			if o.liml
				o.invHⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻, o.kZ)
			else
				o.invHⱼₖ = invsym(isone(o.κ) ? o.H_2SLSⱼₖ : o.ZZⱼₖ + o.κ * o.H_2SLSmZZⱼₖ)
			end
		end

		parent.scorebs &&
			(o.ü₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk])
		if o.liml
			o.H_2SLSmZZ = o.H_2SLS - o.ZZ
		else
			MakeH!(o, parent, !isempty(Rperp)) # DGP is liml except possibly when getting confidence peak for A-R plot; but liml=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
		end
	else
		o.Yendog = [true colsum(o.RparY .!= zero(T)).!=0]
	end

	!o.restricted && (o.t₁Y = zeros(T, parent.kY₂))

	nothing
end


# do most of estimation; for liml r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β̈₀ (β̈ when r=0) and ∂β̈∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function EstimateOLS!(o::StrEstimator, _jk::Bool, r₁::AbstractVector)
	o.β̈ = o.β̈₀ - o.∂β̈∂r * r₁
	_jk && !iszero(nrows(o.R₁perp)) &&
		(o.t₁ = o.R₁invR₁R₁ * r₁)
	nothing
end

function EstimateARubin!(o::StrEstimator{T}, parent::StrBootTest{T}, _jk::Bool, r₁::AbstractVector) where T
	EstimateOLS!(o, _jk, r₁)
  o.y₁par .= parent.y₁ .- parent.Y₂ * r₁
	nothing
end

function MakeH!(o::StrEstimator{T}, parent::StrBootTest{T}, makeXAR::Bool=false) where T
  H = isone(o.κ) ? o.H_2SLS : o.ZZ + o.κ * o.H_2SLSmZZ
  o.invH = invsym(H)

	# parent.jk && (_invH = inv(isone(o.κ) ? _H_2SLS : _ZZ + o.κ * _H_2SLSmZZ))

  if makeXAR  # for replication regression in score bootstrap of IV/GMM
	  o.A = ncols(o.Rperp)>0 ? #=Symmetric=#(o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp') : o.invH
	  o.AR = o.A * (o.Rpar'parent.R')
	  o.XAR = X₁₂B(o.X₁, o.X₂, o.V * o.AR)
  end
	nothing
end


function EstimateIV!(o::StrEstimator{T}, parent::StrBootTest{T}, _jk::Bool, r₁::AbstractVector) where T
  if o.restricted
	  o.t₁Y = o.R₁invR₁R₁Y * r₁

    o.y₁pary₁par = o.y₁y₁ - (o.twoZR₁y₁'r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
	  o.Y₂y₁par = o.Y₂y₁  - o.ZR₁Y₂'r₁
	  o.X₂y₁par = o.X₂y₁ - o.X₂ZR₁ * r₁
	  o.X₁y₁par = o.X₁y₁ - o.X₁ZR₁ * r₁
	  o.Zy₁par  = o.Zy₁ -  o.ZR₁Z'r₁
	  o.Xy₁par  = [o.X₁y₁par ; o.X₂y₁par]
	  (parent.scorebs || parent.granular && o.isDGP) && 
			(o.y₁par .= o.y₁ .- o.ZR₁ * r₁)

		if _jk
			o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ - o.twoZR₁y₁ⱼₖ'r₁ + r₁'o.ZR₁ZR₁ⱼₖ * r₁
			o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ - o.ZR₁Y₂ⱼₖ'r₁
			o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ  - o.ZZR₁ⱼₖ * r₁
			  X₂y₁parⱼₖ    = o.X₂y₁ⱼₖ - o.X₂ZR₁ⱼₖ * r₁
			  X₁y₁parⱼₖ    = o.X₁y₁ⱼₖ - o.X₁ZR₁ⱼₖ * r₁
		  o.Xy₁parⱼₖ     = [X₁y₁parⱼₖ ; X₂y₁parⱼₖ]
		end
  end

  o.invXXXy₁par = o.invXX * o.Xy₁par
  o.YY = #=Symmetric=#([[o.y₁pary₁par] o.Zy₁par'; o.Zy₁par o.ZZ])

  if o.isDGP
		o.ZXinvXXXy₁par = o.XZ'o.invXXXy₁par

		if o.liml
			o.YPXY = #=Symmetric=#([[o.invXXXy₁par'o.Xy₁par] o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par  o.H_2SLS])
	    o.κ = 1/(1 - real(eigvalsNaN(invsym(o.YY) * o.YPXY)[1]))  # like Fast & Wild (81), but more stable, at least in Mata
	    !iszero(o.fuller) && (o.κ -= o.fuller / (parent._Nobs - parent.kX))
	    MakeH!(o, parent)
		end
		o.β̈ = o.invH * (isone(o.κ) ? o.ZXinvXXXy₁par : o.κ * (o.ZXinvXXXy₁par - o.Zy₁par) + o.Zy₁par)

		if _jk
			o.invXXXy₁parⱼₖ .= o.invXXⱼₖ * o.Xy₁parⱼₖ
			o.ZXinvXXXy₁parⱼₖ .= o.XZⱼₖ'o.invXXXy₁parⱼₖ
			o.YYⱼₖ .= [o.y₁pary₁parⱼₖ; o.Zy₁parⱼₖ ;;; o.Zy₁parⱼₖ'; o.ZZⱼₖ]

			if o.liml
				o.YPXYⱼₖ .= [o.invXXXy₁parⱼₖ'o.Xy₁parⱼₖ ; o.ZXinvXXXy₁parⱼₖ ;;; o.ZXinvXXXy₁parⱼₖ' ; o.H_2SLSⱼₖ]
				o.κⱼₖ .= reshape(one(T) ./ (one(T) .- real(getindex.(eigvalsNaN.(each(invsym(o.YYⱼₖ) * o.YPXYⱼₖ)), 1))), (1,:,1))
				!iszero(o.fuller) && (o.κⱼₖ .-= reshape(o.fuller ./ (o.Nobsⱼₖ .- parent.kX)), (1,:,1))
				o.invHⱼₖ .= o.ZZⱼₖ .+ o.κⱼₖ .* o.H_2SLSmZZⱼₖ; invsym!(o.invHⱼₖ)
				o.β̈ⱼₖ .= o.κⱼₖ .* (o.ZXinvXXXy₁parⱼₖ .- o.Zy₁parⱼₖ) .+ o.Zy₁parⱼₖ
			else
				if isone(o.κ)
					o.β̈ⱼₖ.= o.ZXinvXXXy₁parⱼₖ
				else
					o.β̈ⱼₖ .= o.κ .* (o.ZXinvXXXy₁parⱼₖ .- o.Zy₁parⱼₖ) .+ o.Zy₁parⱼₖ
				end
			end
			o.β̈ⱼₖ .= o.invHⱼₖ * o.β̈ⱼₖ
		end
  elseif parent.WREnonARubin  # if not score bootstrap of IV/GMM...
		o.Rt₁ = o.RR₁invR₁R₁ * r₁
	end
	nothing
end

function MakeResidualsOLS!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
	o.ü₁[1] .= o.y₁par .- X₁₂B(parent.X₁, parent.X₂, view(o.β̈ , :,1))   # ordinary non-jk residuals

	if parent.jk
		m = parent.small ? sqrt((parent.N✻ - 1) / T(parent.N✻)) : one(T)
		if parent.purerobust
			if iszero(nrows(o.R₁perp))
				o.ü₁[2] .= m .* o.invMjkv .* o.ü₁[1]
			else
				Xt₁	= X₁₂B(parent.X₁, parent.X₂, o.t₁)
				o.ü₁[2] .= m .* (o.invMjkv .* (o.ü₁[1] .+ Xt₁) .- Xt₁)
			end
		elseif parent.granularjk
	    for g ∈ 1:parent.N✻
				S = parent.info✻[g]
				if nrows(o.R₁perp)>0
					Xt₁	= X₁₂B(view(parent.X₁,S,:), view(parent.X₂,S,:), o.t₁)
					o.ü₁[2][S] .= m .* (o.invMjk[g] * (view(o.ü₁[1], S) + Xt₁) .- Xt₁)
				else
					o.ü₁[2][S] .= m .*  o.invMjk[g] *  view(o.ü₁[1], S)
				end
			end
		else
	    for g ∈ 1:parent.N✻
				S = parent.info✻[g]
				ü₁g = view(o.ü₁[1],S)
				tmp = [view(parent.X₁,S,:)'ü₁g ; view(parent.X₂,S,:)'ü₁g]
				o.ü₁[2][S] .= m .* (ü₁g .+ o.XinvHjk[g] * (iszero(nrows(o.R₁perp)) ? tmp : tmp + view(o.S✻XX,:,g,:) * o.t₁))
			end
		end
	end
	nothing
end

function MakeResidualsIV!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  if parent.scorebs
		o.ü₁[1] .= o.y₁par .- o.Z * o.β̈
	else
		_β = [1 ; -o.β̈ ]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β̈  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  o.Π̂ = invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β̈)' + o.XY₂)
		o.γ̈ = o.RparY * view(o.β̈ ,:,1) + o.t₁Y - parent.Repl.t₁Y
	  if parent.granular
			o.Ü₂[1] .= o.Y₂ .- X₁₂B(o.X₁, o.X₂, o.Π̂ )
			o.u⃛₁[1] .= o.y₁par .- o.Z * view(o.β̈ ,:,1) .+ o.Ü₂[1] * o.γ̈
		end
  end
	nothing
end

# non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
# since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
# depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then κ from DGP liml affects Hessian, H.
function InitTestDenoms!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  if parent.bootstrapt && (parent.scorebs || parent.robust)
	  (parent.granular || parent.purerobust) && (o.WXAR = o.XAR)  # XXX simplify

	  if parent.robust && parent.NFE>0 && !(parent.FEboot || parent.scorebs) && parent.granular < parent.NErrClustCombs  # make first factor of second term of (64) for c=⋂ (c=1)
	    !isdefined(o, :WXAR) && (o.WXAR = o.XAR)  # XXX simplify
	    o.CT_XAR = [crosstabFEt(parent, view(o.WXAR,:,d), parent.info⋂) for d ∈ 1:parent.dof]
	  end
  end
	nothing
end