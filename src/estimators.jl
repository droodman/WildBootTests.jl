# Logically, wild bootstrap tests perform estimation at two stages, once as part of the bootstrap DGP, once in each bootstrap replication
# The StrEstimator "class" and its three "children" hold the estimation logic for the OLS, Anderson-Rubin, and IV/GMM cases

@inline denegate(X) = X .* [any(map(<(0), x)) ? -1 : 1 for x ∈ eachcol(X)]'  # try to turn -1's into 1's for proper selection matrix
@inline identify(X::AbstractMatrix{T}) where T = size(X,1)==size(X,2) && DesignerMatrix(X).type==selection ? Matrix{T}(I(size(X,1))) : X  # try to turn square selection matrix into idnentity: same projection space

function par(X::AbstractMatrix{T}) where T
	F = eigen(Symmetric(X'pinv(X * X')*X))
	return identify(denegate(F.vectors[:, abs.(F.values) .> 1000*eps(T)]))
end
@inline perp(X) = identify(denegate(nullspace(X')))

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 nrows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T}, parent::StrBootTest{T}, R₁::AbstractMatrix{T}, R::Union{UniformScaling{Bool},AbstractMatrix{T}}=Matrix{T}(undef,0,0)) where T
  o.restricted = nrows(R₁) > 0
	if o.restricted
		singular, invR₁R₁ = invsymsingcheck(R₁ * R₁')
		singular && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
	  o.R₁invR₁R₁ = R₁'invR₁R₁
	  o.R₁perp = perp(o.R₁invR₁R₁ * R₁)  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
  else
	  o.R₁invR₁R₁ = Matrix{T}(undef, parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
		o.R₁invR₁R₁X = DesignerMatrix(o.R₁invR₁R₁[1:parent.kX₁,:])
		o.R₁invR₁R₁Y = DesignerMatrix(o.R₁invR₁R₁[parent.kX₁+1:end,:])
		
	  RR₁perp = Matrix{T}([R ; zeros(T, parent.kY₂, parent.kX₁) I])  # rows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

	  o.restricted && (RR₁perp *= o.R₁perp)

	  o.Rpar   = par(RR₁perp)  # defines attack surface for null
	  o.restricted && (o.Rpar = o.R₁perp * o.Rpar)  # fold model constraint factors into Rpar, RperpX
	  o.RRpar = R * o.Rpar
		o.RparX = o.Rpar[1:parent.kX₁,:]  # part of Rpar that refers to Y₂
		o.RparY = DesignerMatrix(o.Rpar[parent.kX₁+1:end,:])  # part of Rpar that refers to Y₂
	  o.RR₁invR₁R₁ = R * o.R₁invR₁R₁

		if o.isDGP || !parent.WREnonARubin  # DGP regressions will just copy the replication-regression Zperp, X₂par for speed
			_RperpX = perp(par(RR₁perp))
			o.restricted && (_RperpX = o.R₁perp * _RperpX)  # fold model constraint factors into Rpar, RperpX
			_RperpX = _RperpX[1:parent.kX₁,:]
			o.RperpX = DesignerMatrix(identify(_RperpX))  # Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X₁, since all endogenous ones will be retained
			o.RperpXperp = DesignerMatrix(perp(_RperpX))
		end
  end
	nothing
end

# stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
function InitVarsOLS!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}) where T # Rperp is for replication regression--no null imposed
  o.y₁par = parent.y₁
	o.ü₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk]

	if parent.jk
		o.S✻XX = panelcross(parent.X₁, parent.X₁, parent.info✻)
	  H = sumpanelcross(o.S✻XX)
	else	
		H    = parent.X₁'parent.X₁
	end

	o.invH = (pinv(H))
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.invH : (o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
	o.β̈₀ = R₁AR₁ * (parent.X₁'parent.y₁)
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

  o.A = iszero(nrows(Rperp)) ? o.invH : (Rperp * invsym(Rperp'H*Rperp) * Rperp')  # for replication regression
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = parent.X₁ * o.AR)
	nothing
end

function InitVarsARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
	o.y₁par = Vector{T}(undef, parent.Nobs)
	o.ü₁    = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk]

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

  H = ([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
  o.A = invsym(H)
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.A : (o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')
	o.β̈₀   = R₁AR₁ * [parent.X₁'parent.y₁ ; parent.X₂'parent.y₁]
	o.∂β̈∂r = R₁AR₁ * [parent.X₁'parent.Y₂ ; parent.X₂'parent.Y₂]

	if parent.jk
		if parent.purerobust
			o.invMjkv =     rowquadform(parent.X₁, (@view o.A[1:parent.kX₁,     1:parent.kX₁    ]), parent.X₁) +
			            2 * rowquadform(parent.X₁, (@view o.A[1:parent.kX₁,     parent.kX₁+1:end]), parent.X₂) +
									    rowquadform(parent.X₂, (@view o.A[parent.kX₁+1:end, parent.kX₁+1:end]), parent.X₂)
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
	prepjk = parent.jk && o.isDGP && parent.WREnonARubin
	
	!isempty(Rperp) && (o.Rperp = Rperp[1])

	o.kZ = ncols(o.Rpar)
	o.kZperp = ncols(parent.DGP.RperpX)

	if !o.isDGP && parent.WREnonARubin
		o.X₁noFWL = parent.DGP.X₁noFWL
		o.Zperp = parent.DGP.Zperp
		o.invZperpZperp = parent.DGP.invZperpZperp
	else
		o.X₁noFWL = parent.X₁ * o.RperpXperp
		o.kX = (o.kX₁ = ncols(o.X₁noFWL)) + parent.kX₂
		o.Zperp = parent.X₁ * o.RperpX
	end
	if parent.WREnonARubin  # faster, doesn't work for score test because only do FWL in DGP, not Repl
		o.Xpar₁toZparX = DesignerMatrix(parent.DGP.RperpXperp'parent.DGP.RperpXperp \ parent.DGP.RperpXperp'o.RparX)
		o.Zpar = o.X₁noFWL * o.Xpar₁toZparX
	else
		o.Zpar = parent.X₁ * o.RparX
	end
	!o.isDGP && !parent.scorebs && (parent.jk || parent.granular) &&
		(o.ZparX = copy(o.Zpar))

	o.Zpar .+= parent.Y₂ * o.RparY
	o.restricted && (o.ZR₁ = parent.X₁ * o.R₁invR₁R₁X; t✻plus!(o.ZR₁, parent.Y₂, o.R₁invR₁R₁Y))

	if prepjk
		o.β̈ⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻,1)
		o.YYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
		o.invXXXy₁parⱼₖ = Array{T,3}(undef, o.kX, parent.N✻, 1)
		o.ZXinvXXXy₁parⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻, 1)
		if o.liml
			o.κⱼₖ = Array{T,3}(undef, 1, parent.N✻, 1)
			o.YPXYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
		end

		ZperpX₁   , _ZperpX₁    = crossjk(o.Zperp, o.X₁noFWL, parent.info✻)  # full-sample and delete-g cross-moments
		ZperpX₂   , _ZperpX₂    = crossjk(o.Zperp, parent.X₂ , parent.info✻)
		Zperpy₁   , _Zperpy₁    = crossjk(o.Zperp, parent.y₁ , parent.info✻)
		ZperpY₂   , _ZperpY₂    = crossjk(o.Zperp, parent.Y₂ , parent.info✻)
		ZperpZ    , _ZperpZ     = crossjk(o.Zperp, o.Zpar    , parent.info✻)
		ZperpZR₁  , _ZperpZR₁   = crossjk(o.Zperp, o.ZR₁     , parent.info✻)

		ZperpZperp, o.invZperpZperp, _invZperpZperp = invsymcrossjk(o.Zperp, parent.info✻)

		_invZperpZperpZperpX₁  = _invZperpZperp * _ZperpX₁
    _invZperpZperpZperpX₂  = _invZperpZperp * _ZperpX₂
    _invZperpZperpZperpZ   = _invZperpZperp * _ZperpZ
    _invZperpZperpZperpZR₁ = _invZperpZperp * _ZperpZR₁
    _invZperpZperpZperpY₂  = _invZperpZperp * _ZperpY₂
    _invZperpZperpZperpy₁  = _invZperpZperp * _Zperpy₁

		o.X₁ⱼₖ  = partialjk(o.X₁noFWL , o.Zperp, _invZperpZperpZperpX₁, parent.info✻)    # FWL-process
		o.X₂ⱼₖ  = partialjk(parent.X₂, o.Zperp, _invZperpZperpZperpX₂ , parent.info✻)
		o.y₁ⱼₖ  = partialjk(parent.y₁, o.Zperp, _invZperpZperpZperpy₁ , parent.info✻)
		o.Y₂ⱼₖ  = partialjk(parent.Y₂, o.Zperp, _invZperpZperpZperpY₂ , parent.info✻)
		o.Zⱼₖ   = partialjk(  o.Zpar , o.Zperp, _invZperpZperpZperpZ  , parent.info✻)
		o.ZR₁ⱼₖ = partialjk(    o.ZR₁, o.Zperp, _invZperpZperpZperpZR₁, parent.info✻)

		tX₁ = ZperpX₁ - ZperpZperp * _invZperpZperpZperpX₁
    tX₂ = ZperpX₂ - ZperpZperp * _invZperpZperpZperpX₂
    tZ  = ZperpZ  - ZperpZperp * _invZperpZperpZperpZ 
    tY₂ = ZperpY₂ - ZperpZperp * _invZperpZperpZperpY₂
    ty₁ = Zperpy₁ .- ZperpZperp * _invZperpZperpZperpy₁

		_y₁ = view(parent.y₁,:,:)
		_X₂X₁   = parent.X₂'o.X₁noFWL      - ZperpX₂'_invZperpZperpZperpX₁ - _invZperpZperpZperpX₂'tX₁ - panelcross(o.X₂ⱼₖ, o.X₁ⱼₖ, parent.info✻)
    _X₁X₁   = o.X₁noFWL'o.X₁noFWL           - ZperpX₁'_invZperpZperpZperpX₁ - _invZperpZperpZperpX₁'tX₁ - panelcross(o.X₁ⱼₖ, o.X₁ⱼₖ, parent.info✻)
    _X₂X₂   = parent.X₂'parent.X₂ - ZperpX₂'_invZperpZperpZperpX₂ - _invZperpZperpZperpX₂'tX₂ - panelcross(o.X₂ⱼₖ, o.X₂ⱼₖ, parent.info✻)
    _X₁Y₂   = o.X₁noFWL'parent.Y₂      - ZperpX₁'_invZperpZperpZperpY₂ - _invZperpZperpZperpX₁'tY₂ - panelcross(o.X₁ⱼₖ, o.Y₂ⱼₖ, parent.info✻)
    _X₂Y₂   = parent.X₂'parent.Y₂ - ZperpX₂'_invZperpZperpZperpY₂ - _invZperpZperpZperpX₂'tY₂ - panelcross(o.X₂ⱼₖ, o.Y₂ⱼₖ, parent.info✻)
    o.Y₂y₁ⱼₖ = parent.Y₂'_y₁       - ZperpY₂'_invZperpZperpZperpy₁ - _invZperpZperpZperpY₂'ty₁ - panelcross(o.Y₂ⱼₖ, o.y₁ⱼₖ, parent.info✻)
    o.X₂y₁ⱼₖ = parent.X₂'_y₁       - ZperpX₂'_invZperpZperpZperpy₁ - _invZperpZperpZperpX₂'ty₁ - panelcross(o.X₂ⱼₖ, o.y₁ⱼₖ, parent.info✻)
    o.X₁y₁ⱼₖ = o.X₁noFWL'_y₁            - ZperpX₁'_invZperpZperpZperpy₁ - _invZperpZperpZperpX₁'ty₁ - panelcross(o.X₁ⱼₖ, o.y₁ⱼₖ, parent.info✻)
    o.Zy₁ⱼₖ  = o.Zpar'_y₁             - ZperpZ'_invZperpZperpZperpy₁  - _invZperpZperpZperpZ'ty₁  - panelcross(o.Zⱼₖ,  o.y₁ⱼₖ, parent.info✻)
    o.XZⱼₖ   = [o.X₁noFWL'o.Zpar           - ZperpX₁'_invZperpZperpZperpZ  - _invZperpZperpZperpX₁'tZ  - panelcross(o.X₁ⱼₖ, o.Zⱼₖ,  parent.info✻)
                parent.X₂'o.Zpar      - ZperpX₂'_invZperpZperpZperpZ  - _invZperpZperpZperpX₂'tZ  - panelcross(o.X₂ⱼₖ, o.Zⱼₖ, parent.info✻)]
    o.ZZⱼₖ   = o.Zpar'o.Zpar             - ZperpZ'_invZperpZperpZperpZ   - _invZperpZperpZperpZ'tZ   - panelcross(o.Zⱼₖ,  o.Zⱼₖ,  parent.info✻)
    o.ZY₂ⱼₖ  = o.Zpar'parent.Y₂       - ZperpZ'_invZperpZperpZperpY₂  - _invZperpZperpZperpZ'tY₂  - panelcross(o.Zⱼₖ,  o.Y₂ⱼₖ, parent.info✻)
    o.y₁y₁ⱼₖ = _y₁'_y₁ -2 * (Zperpy₁'_invZperpZperpZperpy₁) + _invZperpZperpZperpy₁'ZperpZperp*_invZperpZperpZperpy₁ - panelcross(o.y₁ⱼₖ, o.y₁ⱼₖ, parent.info✻)

    o.XY₂ⱼₖ = [_X₁Y₂ ; _X₂Y₂]
    o.XXⱼₖ  = [_X₁X₁ ; _X₂X₁ ;;; _X₂X₁' ;  _X₂X₂]
		o.invXXⱼₖ = invsym(o.XXⱼₖ)
    o.H_2SLSⱼₖ = o.XZⱼₖ'o.invXXⱼₖ * o.XZⱼₖ
    (!isone(o.κ) || o.liml) && (o.H_2SLSmZZⱼₖ = o.H_2SLSⱼₖ - o.ZZⱼₖ)

		if o.restricted
      tZR₁ = ZperpZR₁ - ZperpZperp * _invZperpZperpZperpZR₁
      o.X₁ZR₁ⱼₖ    =  o.X₁noFWL'o.ZR₁      - ZperpX₁'_invZperpZperpZperpZR₁  - _invZperpZperpZperpX₁'tZR₁  - panelcross(o.X₁ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻)
      o.X₂ZR₁ⱼₖ    =  parent.X₂'o.ZR₁ - ZperpX₂'_invZperpZperpZperpZR₁  - _invZperpZperpZperpX₂'tZR₁  - panelcross(o.X₂ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻)
      o.ZZR₁ⱼₖ     =  o.Zpar'o.ZR₁       - ZperpZ'_invZperpZperpZperpZR₁   - _invZperpZperpZperpZ'tZR₁   - panelcross(o.Zⱼₖ,   o.ZR₁ⱼₖ, parent.info✻)
      o.twoZR₁y₁ⱼₖ =2(o.ZR₁'_y₁       - ZperpZR₁'_invZperpZperpZperpy₁  - _invZperpZperpZperpZR₁'ty₁  - panelcross(o.ZR₁ⱼₖ,  o.y₁ⱼₖ, parent.info✻))
      o.ZR₁ZR₁ⱼₖ   =  o.ZR₁'o.ZR₁     - ZperpZR₁'_invZperpZperpZperpZR₁ - _invZperpZperpZperpZR₁'tZR₁ - panelcross(o.ZR₁ⱼₖ, o.ZR₁ⱼₖ, parent.info✻)
      o.ZR₁Y₂ⱼₖ    =  o.ZR₁'parent.Y₂ - ZperpZR₁'_invZperpZperpZperpY₂  - _invZperpZperpZperpZR₁'tY₂  - panelcross(o.ZR₁ⱼₖ, o.Y₂ⱼₖ,  parent.info✻)
			o.y₁parⱼₖ    = Vector{T}(undef, parent.Nobs)
		else
			o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ 
			o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ
			o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ 
			o.Xy₁parⱼₖ     = [o.X₁y₁ⱼₖ ; o.X₂y₁ⱼₖ]
			o.y₁parⱼₖ      = o.y₁ⱼₖ
		end

		!iszero(o.fuller) &&
			(o.Nobsⱼₖ = parent._Nobs .- (parent.fweights ? @panelsum(parent.wt, parent.info✻) : T.(length.(parent.info✻))))
	elseif parent.granular  # if no jk prep, leave behind same directly computed cross-products; except coarse treatment will sum them from panelwise cross-products for speed
		if o.isDGP
			ZperpX₁ = o.Zperp'o.X₁noFWL
			ZperpX₂ = o.Zperp'parent.X₂
			Zperpy₁ = o.Zperp'parent.y₁
			ZperpY₂ = o.Zperp'parent.Y₂
		end
		ZperpZ   = o.Zperp'o.Zpar
		o.restricted && (ZperpZR₁ = o.Zperp'o.ZR₁)
		o.invZperpZperp = invsym(o.Zperp'o.Zperp)
	end

	if !o.isDGP && parent.WREnonARubin
		o.Zperp = parent.DGP.Zperp
		o.X₁noFWL = parent.DGP.X₁noFWL
		isdefined(parent.DGP, :X₁) && (o.X₁ = parent.DGP.X₁)
		isdefined(parent.DGP, :X₂) && (o.X₂ = parent.DGP.X₂)
		isdefined(parent.DGP, :Y₂) && (o.Y₂ = parent.DGP.Y₂)
		isdefined(parent.DGP, :y₁) && (o.y₁ = parent.DGP.y₁)
		o.invZperpZperp = parent.DGP.invZperpZperp
		o.XX = parent.DGP.XX
		o.invXX = parent.DGP.invXX
		o.XY₂ = parent.DGP.XY₂
		o.Y₂Y₂ = parent.DGP.Y₂Y₂
		o.Y₂y₁ = parent.DGP.Y₂y₁
		o.X₂y₁ = parent.DGP.X₂y₁
		o.X₁y₁ = parent.DGP.X₁y₁
		o.y₁y₁ = parent.DGP.y₁y₁
	end

	# copy or construct objects that are same in DGP and Repl, including O(N) ones
	if parent.granular || prepjk
		if o.isDGP || !parent.WREnonARubin
			o.X₁ = o.Zperp * (o.invZperpZperp * ZperpX₁); o.X₁ .= o.X₁noFWL .- o.X₁  # FWL-processing
			o.X₂ = o.Zperp * (o.invZperpZperp * ZperpX₂); o.X₂ .= parent.X₂  .- o.X₂
			o.y₁ = o.Zperp * (o.invZperpZperp * Zperpy₁); o.y₁ .= parent.y₁  .- o.y₁
			o.Y₂ = o.Zperp * (o.invZperpZperp * ZperpY₂); o.Y₂ .= parent.Y₂  .- o.Y₂

			X₂X₁ = o.X₂'o.X₁
			o.XX = [o.X₁'o.X₁ X₂X₁' ; X₂X₁ o.X₂'o.X₂]
			o.invXX = invsym(o.XX)
			X₁Y₂ = o.X₁'o.Y₂
			X₂Y₂ = o.X₂'o.Y₂
			o.XY₂ = [X₁Y₂ ; X₂Y₂]
			o.Y₂y₁ = o.Y₂'o.y₁
			o.Y₂Y₂ = o.Y₂'o.Y₂
			o.X₂y₁ = o.X₂'o.y₁
			o.X₁y₁ = o.X₁'o.y₁
			o.y₁y₁ = dot(o.y₁, o.y₁)
		end

		t✻minus!(o.Zpar, o.Zperp, o.invZperpZperp * ZperpZ  )
		o.restricted && t✻minus!(o.ZR₁ , o.Zperp, o.invZperpZperp * ZperpZR₁)

		o.Zy₁ = o.Zpar'o.y₁
		o.ZY₂ = o.Zpar'o.Y₂
		o.ZZ  = o.Zpar'o.Zpar
		o.XZ  = [o.X₁'o.Zpar ; o.X₂'o.Zpar]
	
		if o.restricted
			o.X₂ZR₁    = o.X₂'o.ZR₁
			o.X₁ZR₁    = o.X₁'o.ZR₁
			o.ZR₁Z     = o.ZR₁'o.Zpar
			o.twoZR₁y₁ = o.ZR₁'o.y₁; o.twoZR₁y₁ .*= 2
			o.ZR₁ZR₁   = o.ZR₁'o.ZR₁
			o.ZR₁Y₂    = o.ZR₁'o.Y₂
		end

		if o.isDGP && !parent.scorebs
			o.u⃛₁ = Vector{T}(undef, parent.Nobs)
			o.Ü₂ = Matrix{T}(undef, parent.Nobs, parent.kY₂)
		end
		o.ȳ₁ = Vector{T}(undef, parent.Nobs)
		o.Ȳ₂ = Matrix{T}(undef, parent.Nobs, parent.kY₂)
	end

	if !parent.granular
		if !o.isDGP && parent.WREnonARubin
			o.S✻Y₂Y₂ = parent.DGP.S✻Y₂Y₂
			o.invZperpZperpZperpX₁ = parent.DGP.invZperpZperpZperpX₁
			o.invZperpZperpZperpX₂ = parent.DGP.invZperpZperpZperpX₂
			o.invZperpZperpZperpy₁ = parent.DGP.invZperpZperpZperpy₁
			o.invZperpZperpZperpY₂ = parent.DGP.invZperpZperpZperpY₂
			o.S✻⋂X₁Y₂ = parent.DGP.S✻⋂X₁Y₂
			o.S✻⋂X₂Y₂ = parent.DGP.S✻⋂X₂Y₂
			o.S✻⋂ZperpY₂ = parent.DGP.S✻⋂ZperpY₂
			o.S✻Y₂y₁ = parent.DGP.S✻Y₂y₁
		else  # XXX how much of this stuff is already done if jk, especially if no subcluster?
			o.X₁noFWL = parent.X₁ * o.RperpXperp  # X∥ := [X₁noFWL X₂]

			o.S✻⋂ZperpZperp = panelcross(o.Zperp, o.Zperp, parent.info✻⋂)
			S✻⋂X₁Zperp = panelcross(o.X₁noFWL, o.Zperp, parent.info✻⋂)
			S✻⋂X₂Zperp = panelcross(parent.X₂, o.Zperp, parent.info✻⋂)

			if !(parent.jk && parent.WREnonARubin)
				o.invZperpZperp = iszero(o.kZperp) ? Matrix{T}(undef,0,0) : invsym(sumpanelcross(o.S✻⋂ZperpZperp))
				ZperpX₁ = sumpanelcross(S✻⋂X₁Zperp)'
				ZperpX₂ = sumpanelcross(S✻⋂X₂Zperp)'
			end
			
			o.invZperpZperpZperpX₁ = o.invZperpZperp * ZperpX₁
			o.S✻⋂XZperp = [S✻⋂X₁Zperp; S✻⋂X₂Zperp]

			o.invZperpZperpZperpX₂ = o.invZperpZperp * ZperpX₂
			o.invZperpZperpZperpX = [o.invZperpZperpZperpX₁ o.invZperpZperpZperpX₂]

			S✻⋂X₁X₁ = panelcross(o.X₁noFWL, o.X₁noFWL, parent.info✻⋂)
			S✻⋂X₂X₁ = panelcross(parent.X₂, o.X₁noFWL, parent.info✻⋂)
			S✻⋂X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻⋂)
			X₂X₁ = sumpanelcross(S✻⋂X₂X₁)
			X₁X₁ = sumpanelcross(S✻⋂X₁X₁)
			X₂X₂ = sumpanelcross(S✻⋂X₂X₂)
			o.S✻⋂XX = [[S✻⋂X₁X₁ S✻⋂X₂X₁'] ; [S✻⋂X₂X₁ S✻⋂X₂X₂]]  # [a b; c d] syntax would call hvcat() to concatenate horizontally along dim 2 rather than 3

			X₂X₁ .-= ZperpX₂'o.invZperpZperp * ZperpX₁
			X₁X₁ .-= ZperpX₁'o.invZperpZperp * ZperpX₁
			X₂X₂ .-= ZperpX₂'o.invZperpZperp * ZperpX₂
			o.XX = ([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
			o.invXX = invsym(o.XX)
			o.kX = ncols(o.XX)
			o.kX₁ = ncols(X₂X₁)

			o.S✻⋂ZperpY₂ = panelcross(o.Zperp, parent.Y₂, parent.info✻⋂)
			!(parent.jk && parent.WREnonARubin) && (ZperpY₂ = sumpanelcross(o.S✻⋂ZperpY₂))
			o.invZperpZperpZperpY₂ = o.invZperpZperp * ZperpY₂
			(parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt)) &&
				(o.Y₂ = parent.Y₂ - o.Zperp * o.invZperpZperpZperpY₂)
			o.S✻⋂Zperpy₁ = panelcross(o.Zperp, parent.y₁, parent.info✻⋂)
			!(parent.jk && parent.WREnonARubin) && (Zperpy₁ = vec(sumpanelcross(o.S✻⋂Zperpy₁)))
			o.invZperpZperpZperpy₁ = o.invZperpZperp * Zperpy₁
			((parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt)) || parent.scorebs) &&
				(o.y₁ = parent.y₁ - o.Zperp * o.invZperpZperpZperpy₁)
		  o.S✻⋂X₁Y₂ = panelcross(o.X₁noFWL, parent.Y₂, parent.info✻⋂)
		  o.S✻⋂X₂Y₂ = panelcross(parent.X₂, parent.Y₂, parent.info✻⋂)
			o.S✻⋂XY₂ = [o.S✻⋂X₁Y₂; o.S✻⋂X₂Y₂]

			o.XY₂ = sumpanelcross(o.S✻⋂XY₂)
			o.S✻Y₂y₁ = panelcross(parent.Y₂, parent.y₁, parent.info✻)
		  o.Y₂y₁ = vec(sumpanelcross(o.S✻Y₂y₁))
			o.S✻Y₂Y₂ = panelcross(parent.Y₂, parent.Y₂, parent.info✻)
			o.Y₂Y₂ = (sumpanelcross(o.S✻Y₂Y₂))
			S✻⋂X₂y₁ = panelcross(parent.X₂, parent.y₁, parent.info✻⋂)
			o.X₂y₁ = vec(sumpanelcross(S✻⋂X₂y₁))
			S✻⋂X₁y₁ = panelcross(o.X₁noFWL, parent.y₁, parent.info✻⋂)
			o.S✻⋂Xy₁ = [S✻⋂X₁y₁; S✻⋂X₂y₁]
			o.X₁y₁ = vec(sumpanelcross(S✻⋂X₁y₁))
			o.S✻y₁y₁ = panelcross(parent.y₁, parent.y₁, parent.info✻)
			o.y₁y₁ = sum(o.S✻y₁y₁)
			o.XY₂  .-= o.invZperpZperpZperpX'ZperpY₂
		  o.Y₂y₁ .-= ZperpY₂'o.invZperpZperpZperpy₁
			o.Y₂Y₂ .-= ZperpY₂'o.invZperpZperpZperpY₂
			o.X₂y₁ .-= ZperpX₂'o.invZperpZperpZperpy₁

			o.X₁y₁ .-= ZperpX₁'o.invZperpZperpZperpy₁
			o.y₁y₁  -= Zperpy₁'o.invZperpZperpZperpy₁

			if !isdefined(o, :X₁) && parent.NFE>0 && (parent.liml || !isone(parent.κ) || parent.bootstrapt) || !o.liml && !isempty(Rperp)
				o.X₁ = o.X₁noFWL - o.Zperp * o.invZperpZperpZperpX₁  # shrink and FWL-process X₁; do it as an O(N) operation because it can be so size-reducing
				o.X₂ = o.Zperp * o.invZperpZperpZperpX₂; o.X₂ .= parent.X₂ .- o.X₂  # FWL-process X₂
			end
		end

		S✻⋂X₁Zpar = panelcross(o.X₁noFWL, o.Zpar, parent.info✻⋂)
		S✻⋂X₂Zpar = panelcross(parent.X₂, o.Zpar, parent.info✻⋂)
		o.S✻⋂XZpar = [S✻⋂X₁Zpar; S✻⋂X₂Zpar]
		X₁Zpar = sumpanelcross(S✻⋂X₁Zpar)
		X₂Zpar = sumpanelcross(S✻⋂X₂Zpar)
		o.S✻⋂ZperpZpar = panelcross(o.Zperp, o.Zpar, parent.info✻⋂)
		ZperpZpar = sumpanelcross(o.S✻⋂ZperpZpar)
		o.invZperpZperpZperpZpar = o.invZperpZperp * ZperpZpar
		o.S✻Zpary₁ = panelcross(o.Zpar, parent.y₁, parent.info✻)
		o.Zy₁ = vec(sumpanelcross(o.S✻Zpary₁)); o.Zy₁ .-= ZperpZpar'o.invZperpZperpZperpy₁
		o.S✻ZparY₂ = panelcross(o.Zpar, parent.Y₂, parent.info✻)
		o.XZ = [X₁Zpar - o.invZperpZperpZperpX₁'ZperpZpar ; X₂Zpar - o.invZperpZperpZperpX₂'ZperpZpar]
	 	o.S✻ZparZpar = panelcross(o.Zpar, o.Zpar, parent.info✻)
		o.ZZ = sumpanelcross(o.S✻ZparZpar); o.ZZ .-= ZperpZpar'o.invZperpZperpZperpZpar

	  if o.restricted
			S✻⋂X₁ZR₁ = panelcross(o.X₁noFWL, o.ZR₁, parent.info✻⋂)
			S✻⋂X₂ZR₁ = panelcross(parent.X₂, o.ZR₁, parent.info✻⋂)
			o.S✻⋂XZR₁ = [S✻⋂X₁ZR₁ ; S✻⋂X₂ZR₁]
			o.S✻⋂ZperpZR₁ = panelcross(o.Zperp, o.ZR₁, parent.info✻⋂)
			o.ZperpZR₁ = sumpanelcross(o.S✻⋂ZperpZR₁)
			o.invZperpZperpZperpZR₁ = o.invZperpZperp * o.ZperpZR₁
		  o.X₁ZR₁    = sumpanelcross(S✻⋂X₁ZR₁)
		  o.X₂ZR₁    = sumpanelcross(S✻⋂X₂ZR₁)
			o.S✻ZR₁Z   = panelcross(o.ZR₁, o.Zpar, parent.info✻)
		  o.ZR₁Z     = sumpanelcross(o.S✻ZR₁Z)
			o.S✻ZR₁Y₂  = panelcross(o.ZR₁, parent.Y₂, parent.info✻)
		  o.ZR₁Y₂    = sumpanelcross(o.S✻ZR₁Y₂) 
			o.S✻ZR₁y₁  = panelcross(o.ZR₁, parent.y₁, parent.info✻)
		  o.twoZR₁y₁ = 2 * vec(sumpanelcross(o.S✻ZR₁y₁))
			o.S✻ZR₁ZR₁ = panelcross(o.ZR₁, o.ZR₁, parent.info✻)
		  o.ZR₁ZR₁   = sumpanelcross(o.S✻ZR₁ZR₁)
			(parent.scorebs || (parent.NFE>0 && !parent.FEboot && o.isDGP && o.restricted && (parent.willfill || parent.not2SLS))) &&
				(o.ZR₁ .-= o.Zperp * o.invZperpZperpZperpZR₁)
		  t✻minus!(o.X₁ZR₁, o.invZperpZperpZperpX₁', o.ZperpZR₁)
		  t✻minus!(o.X₂ZR₁, o.invZperpZperpZperpX₂', o.ZperpZR₁)
		  o.ZR₁Z    -= o.ZperpZR₁'o.invZperpZperp * ZperpZpar
		  t✻minus!(o.ZR₁Y₂, o.ZperpZR₁', o.invZperpZperpZperpY₂)
		  o.twoZR₁y₁ .-= 2 * o.ZperpZR₁'o.invZperpZperpZperpy₁
		  o.ZR₁ZR₁   .-= o.ZperpZR₁'o.invZperpZperp * o.ZperpZR₁
		else
		  o.Y₂y₁par    = o.Y₂y₁
		  o.X₂y₁par    = o.X₂y₁
		  o.X₁y₁par    = o.X₁y₁
		  o.Zy₁par     = o.Zy₁
		  o.y₁pary₁par = o.y₁y₁
		  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
		end

		parent.scorebs && (o.y₁par = copy(o.y₁))

		parent.WREnonARubin && !(o.isDGP && parent.scorebs) && 
			(o.ZY₂ = sumpanelcross(o.S✻ZparY₂) - ZperpZpar'o.invZperpZperpZperpY₂)
		
		t✻minus!(o.Zpar, o.Zperp, o.invZperpZperpZperpZpar)
	end  # end coarse

	o.Y₂y₁par    = copy(o.Y₂y₁)  # RHS objects are === across DGP and Repl, but LHS objects can be modified later; break link
	o.X₂y₁par    = copy(o.X₂y₁)
	o.X₁y₁par    = copy(o.X₁y₁)
	o.Zy₁par     = copy(o.Zy₁)
	o.y₁pary₁par = copy(o.y₁y₁)
	o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
	((parent.jk && !o.isDGP) || (parent.granular || parent.scorebs)) &&
		(o.y₁par   = copy(o.y₁))

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

		(parent.scorebs || parent.granular || parent.jk) &&
			(o.ü₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk])
		if o.liml
			o.H_2SLSmZZ = o.H_2SLS - o.ZZ
		else
			MakeH!(o, parent, !isempty(Rperp)) # DGP is liml except possibly when getting confidence peak for A-R plot; but liml=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
		end
	else
		o.Yendog = [true colsum(o.RparY .!= zero(T)).!=0]
		parent.jk && parent.willfill && (o.PXZ = X₁₂B(o.X₁, o.X₂, o.V))
		!parent.scorebs && (parent.jk || parent.granular) &&
			t✻minus!(o.ZparX, o.Zperp, o.invZperpZperp * (o.Zperp'o.ZparX))
	end

	if !o.restricted
		o.t₁ = zeros(T, parent.kZ)
		parent.jk && (o.t₁Y = o.t₁[parent.kX₁+1:end])
	end

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

  if makeXAR  # for replication regression in score bootstrap of IV/GMM
	  o.A = ncols(o.Rperp)>0 ? (o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp') : o.invH
	  o.AR = o.A * (o.Rpar'parent.R')
	  o.XAR = X₁₂B(o.X₁, o.X₂, o.V * o.AR)
  end
	nothing
end


function EstimateIV!(o::StrEstimator{T}, parent::StrBootTest{T}, _jk::Bool, r₁::AbstractVector) where T
  if o.restricted
	  o.t₁ = o.R₁invR₁R₁ * r₁

    o.y₁pary₁par = o.y₁y₁ - (o.twoZR₁y₁'r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
	  o.Y₂y₁par .= o.Y₂y₁  - o.ZR₁Y₂'r₁
	  o.X₂y₁par .= o.X₂y₁ - o.X₂ZR₁ * r₁
	  o.X₁y₁par .= o.X₁y₁ - o.X₁ZR₁ * r₁
	  o.Zy₁par  .= o.Zy₁ -  o.ZR₁Z'r₁
	  o.Xy₁par  .= [o.X₁y₁par ; o.X₂y₁par]
	  (parent.scorebs || parent.granular && o.isDGP) && 
			(o.y₁par .= o.y₁ .- o.ZR₁ * r₁)

		if _jk
			o.t₁Y = o.t₁[parent.kX₁+1:end]
			o.y₁parⱼₖ .= o.y₁ⱼₖ; t✻minus!(o.y₁parⱼₖ, o.ZR₁ⱼₖ, r₁)
			o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ - o.twoZR₁y₁ⱼₖ'r₁ + r₁'o.ZR₁ZR₁ⱼₖ * r₁
			o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ - o.ZR₁Y₂ⱼₖ'r₁
			o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ  - o.ZZR₁ⱼₖ * r₁
			  X₂y₁parⱼₖ    = o.X₂y₁ⱼₖ - o.X₂ZR₁ⱼₖ * r₁
			  X₁y₁parⱼₖ    = o.X₁y₁ⱼₖ - o.X₁ZR₁ⱼₖ * r₁
		  o.Xy₁parⱼₖ     = [X₁y₁parⱼₖ ; X₂y₁parⱼₖ]
		end
  end

  o.invXXXy₁par = o.invXX * o.Xy₁par
  o.YY = ([[o.y₁pary₁par] o.Zy₁par'; o.Zy₁par o.ZZ])

  if o.isDGP
		o.ZXinvXXXy₁par = o.XZ'o.invXXXy₁par

		if o.liml
			o.YPXY = ([[o.invXXXy₁par'o.Xy₁par] o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par  o.H_2SLS])
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
				o.κⱼₖ .= reshape(one(T) ./ (one(T) .- getindex.(real.(eigvalsNaN.(each(invsym(o.YYⱼₖ) * o.YPXYⱼₖ))), 1)), (1,:,1))
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
	if parent.jk
		o.ü₁[1] .= o.y₁parⱼₖ
		for (g,S) ∈ enumerate(parent.info✻)
			t✻minus!(view(o.ü₁[1],S), view(o.Zⱼₖ,S,:), view(o.β̈ⱼₖ,:,g,1))
		end
	elseif parent.granular || parent.scorebs
		o.ü₁[1] .= o.y₁par .- o.Zpar * o.β̈  
	end

  if !parent.scorebs
		_β = [1 ; -o.β̈ ]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β̈  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  o.Π̈ = invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β̈ )' + o.XY₂)
		o.γ̈ = o.Rpar * o.β̈  + o.t₁ - parent.Repl.t₁
    o.Ü₂Ü₂ = o.Y₂Y₂ - (o.Π̈)'o.XY₂ - o.XY₂'o.Π̈ + (o.Π̈)'o.XX*o.Π̈

    o.γ̈X = o.γ̈[1:parent.kX₁]
    o.γ̈Y = o.γ̈[parent.kX₁+1:end]

    tmp = o.RperpXperp'o.γ̈X
    o.γ⃛ = o.Π̈ * o.γ̈Y; o.γ⃛[1:o.kX₁] += tmp  # (X_∥'X_∥)^-1 * X_∥'y1bar
    o.Xȳ₁ = o.XX * o.γ⃛
    o.ȳ₁ȳ₁ = o.γ⃛'o.Xȳ₁
    o.XÜ₂ = o.XY₂ - o.XX * o.Π̈	
    o.ȳ₁Ü₂ = o.γ⃛'o.XÜ₂

    if parent.granular || parent.jk
      X₁₂B!(o.Ȳ₂, o.X₁, o.X₂, o.Π̈ )  
      mul!(o.ȳ₁, o.X₁, tmp); t✻plus!(o.ȳ₁, o.Ȳ₂, o.γ̈Y)
		end

		if parent.jk
			o.Ü₂ .= o.Y₂ⱼₖ
			o.u⃛₁ .= o.ü₁[1]
			_t₁Y = o.t₁Y - parent.Repl.t₁Y
			Π̈ⱼₖ = Matrix{T}(undef, o.kX, parent.kY₂)
			for (g,S) ∈ enumerate(parent.info✻)
				_β .= [1 ; -view(o.β̈ⱼₖ,:,g,1)]
				uu = _β'view(o.YYⱼₖ,:,g,:) * _β
		
				Xu .= view(o.Xy₁parⱼₖ,:,g,:); t✻minus!(Xu, view(o.XZⱼₖ,:,g,:), view(o.β̈ⱼₖ,:,g,1))  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
				negXuinvuu .= Xu ./ -uu
				Π̈ⱼₖ .= invsym(view(o.XXⱼₖ,:,g,:) + negXuinvuu * Xu') * (negXuinvuu * (view(o.Y₂y₁parⱼₖ,:,g,:) - view(o.ZY₂ⱼₖ,:,g,:)'view(o.β̈ⱼₖ,:,g,1))' + view(o.XY₂ⱼₖ,:,g,:))
				X₁₂Bminus!(view(o.Ü₂,S,:), view(o.X₁ⱼₖ,S,:), view(o.X₂ⱼₖ,S,:), Π̈ⱼₖ)
				t✻plus!(view(o.u⃛₁,S), view(o.Ü₂,S,:), o.RparY * view(o.β̈ⱼₖ,:,g,1) + _t₁Y)
			end
		elseif parent.granular
			o.Ü₂ .= o.Y₂; X₁₂Bminus!(o.Ü₂, o.X₁, o.X₂, o.Π̈ )
			o.u⃛₁ .= o.ü₁[1]; t✻plus!(o.u⃛₁, o.Ü₂, o.γ̈Y)
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