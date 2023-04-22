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
	o.copyfromDGP = !o.isDGP && parent.WREnonARubin
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

		if !o.copyfromDGP  # DGP regressions will just copy the replication-regression Zperp, X₂par for speed
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

	if parent.jk && !(parent.purerobust || parent.granularjk)
		o.S✻XX = panelcross(parent.X₁, parent.X₁, parent.info✻)
	  H = sumpanelcross(o.S✻XX)
	else	
		H = parent.X₁'parent.X₁
	end

	o.invH = invsym(H)
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.invH : (o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
	o.β̈₀ = R₁AR₁ * (parent.X₁'parent.y₁)
	o.∂β̈∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

	if parent.jk
		if parent.purerobust
			o.invMjkv = rowquadform(o.invH, parent.X₁)
			o.invMjkv .= 1 ./ (1 .- o.invMjkv)  # standard hc3 multipliers
		elseif parent.granularjk
			negR₁AR₁ = -R₁AR₁  # N.B.: likely that R₁AR₁===A, so don't overwrite it
			o.invMjk = Vector{Matrix{T}}(undef, parent.N✻)
			M       = Matrix{T}(undef, parent.maxNg, parent.maxNg)
			X₁R₁AR₁ = Matrix{T}(undef, parent.maxNg, parent.kX₁  )
			for (g,S) ∈ enumerate(parent.info✻)
				# compute X₁_g R₁AR₁ X₁_g' while minimizing allocations
				Ng = length(S)
				Mⱼₖ = squaresubview(M, Ng)
				_XR₁AR₁ = rowsubview(X₁R₁AR₁, Ng)
				_X₁ = view(parent.X₁, S, :)
				t✻!(_XR₁AR₁, _X₁, negR₁AR₁)
				t✻!(Mⱼₖ, _XR₁AR₁, _X₁')
				view(Mⱼₖ, diagind(Mⱼₖ)) .+= one(T)  # add I
				o.invMjk[g] = invsym(Mⱼₖ)
			end
			if o.restricted  # scratch matrices for MakeResidualsOLS!()
				o.Xt₁    = Matrix{T}(undef, parent.maxNg, parent.dof)
				o.Xt₁pu  = Matrix{T}(undef, parent.maxNg, parent.dof)
				o.MXt₁pu = Matrix{T}(undef, parent.maxNg, parent.dof)
			end
		else
			o.Xu = Vector{T}(undef, parent.kX)
			_H = reshape(H, (parent.kX, 1, parent.kX)) .- o.S✻XX
			_invH = iszero(nrows(o.R₁perp)) ? invsym(_H) : o.R₁perp * invsym(o.R₁perp' * _H * o.R₁perp) * o.R₁perp'
			o.XinvHⱼₖ = [view(parent.X₁, parent.info✻[g],:) * view(_invH,:,g,:) for g ∈ 1:parent.N✻]
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

	if parent.jk && !(parent.granularjk || parent.purerobust)
		S✻X₁X₁ = panelcross(parent.X₁, parent.X₁, parent.info✻)
		S✻X₂X₁ = panelcross(parent.X₂, parent.X₁, parent.info✻)
		S✻X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻)
		X₁X₁ = sumpanelcross(S✻X₁X₁)
		X₂X₁ = sumpanelcross(S✻X₂X₁)
		X₂X₂ = sumpanelcross(S✻X₂X₂)
	else
		X₁X₁ = parent.X₁'parent.X₁
		X₂X₁ = parent.X₂'parent.X₁
		X₂X₂ = parent.X₂'parent.X₂	
	end

  H = [X₁X₁ X₂X₁' ; X₂X₁ X₂X₂]
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
			negR₁AR₁ = -R₁AR₁  # N.B.: likely that R₁AR₁===A, so don't overwrite it
			negR₁AR₁_₁₁ = @view negR₁AR₁[1:parent.kX₁    , 1:parent.kX₁    ]
			negR₁AR₁_₁₂ = @view negR₁AR₁[1:parent.kX₁    , parent.kX₁+1:end]
			negR₁AR₁_₂₁ = @view negR₁AR₁[parent.kX₁+1:end, 1:parent.kX₁    ]
			negR₁AR₁_₂₂ = @view negR₁AR₁[parent.kX₁+1:end, parent.kX₁+1:end]
			M = Matrix{T}(undef, parent.maxNg, parent.maxNg)
			X₁R₁AR₁ = Matrix{T}(undef, parent.maxNg, parent.kX₁)
			X₂R₁AR₁ = Matrix{T}(undef, parent.maxNg, parent.kX₂)
			for (g,S) ∈ enumerate(parent.info✻)
				# compute -[X₁ X₂]_g R₁AR₁ [X₁ X₂]_g' while minimizing allocations
				Ng = length(S)
				_X₁ = view(parent.X₁, S, :); _XR₁AR₁_₁ = rowsubview(X₁R₁AR₁, Ng)
				_X₂ = view(parent.X₂, S, :); _XR₁AR₁_₂ = rowsubview(X₂R₁AR₁, Ng)
				t✻!(_XR₁AR₁_₁, _X₁, negR₁AR₁_₁₁); t✻plus!(_XR₁AR₁_₁, _X₂, negR₁AR₁_₂₁)
				t✻!(_XR₁AR₁_₂, _X₁, negR₁AR₁_₁₂); t✻plus!(_XR₁AR₁_₂, _X₂, negR₁AR₁_₂₂)
				Mⱼₖ = squaresubview(M, Ng)
				t✻!(Mⱼₖ, _XR₁AR₁_₁, _X₁'); t✻plus!(Mⱼₖ, _XR₁AR₁_₂, _X₂')

				view(Mⱼₖ, diagind(Mⱼₖ)) .+= one(T)  # add I
				o.invMjk[g] = invsym(Mⱼₖ)
			end
		else
			o.Xu = Vector{T}(undef, parent.kX)
			o.S✻XX = [[S✻X₁X₁ S✻X₂X₁'] ; [S✻X₂X₁ S✻X₂X₂]]
			_H = reshape(H, (parent.kX, 1, parent.kX)) .- o.S✻XX
			_invH = iszero(nrows(o.R₁perp)) ? invsym(_H) : o.R₁perp * invsym(o.R₁perp' * _H * o.R₁perp) * o.R₁perp'
			o.XinvHⱼₖ = [(S = parent.info✻[g]; X₁₂B(view(parent.X₁, S,:), view(parent.X₂, S,:), view(_invH,:,g,:))) for g ∈ 1:parent.N✻]
		end
	end

  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = X₁₂B(parent.X₁, parent.X₂, o.AR))
	nothing
end

function PrepJKIV!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
	o.β̈ⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻,1)
	o.YYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
	o.invXXXy₁parⱼₖ = Array{T,3}(undef, o.kX, parent.N✻, 1)
	o.ZXinvXXXy₁parⱼₖ = Array{T,3}(undef, o.kZ, parent.N✻, 1)
	if o.liml
		o.κⱼₖ = Array{T,3}(undef, 1, parent.N✻, 1)
		o.YPXYⱼₖ = Array{T,3}(undef, o.kZ+1, parent.N✻, o.kZ+1)
	end

	ZperpX₁  , _ZperpX₁   = crossjk(o.Zperp, o.X₁, parent.info✻)  # full-sample and delete-g cross-moments
	ZperpX₂  , _ZperpX₂   = crossjk(o.Zperp, parent.X₂, parent.info✻)
	Zperpy₁  , _Zperpy₁   = crossjk(o.Zperp, parent.y₁, parent.info✻)
	ZperpY₂  , _ZperpY₂   = crossjk(o.Zperp, parent.Y₂, parent.info✻)
	ZperpZpar, _ZperpZpar = crossjk(o.Zperp, o.Zpar   , parent.info✻)

	ZperpZperp, o.invZperpZperp, _invZperpZperp = invsymcrossjk(o.Zperp, parent.info✻)

	_invZperpZperpZperpX₁  = _invZperpZperp * _ZperpX₁
	_invZperpZperpZperpX₂  = _invZperpZperp * _ZperpX₂
	_invZperpZperpZperpZ   = _invZperpZperp * _ZperpZpar
	_invZperpZperpZperpY₂  = _invZperpZperp * _ZperpY₂
	_invZperpZperpZperpy₁  = _invZperpZperp * _Zperpy₁

	o.X₁ⱼₖ  = partialjk(o.X₁, o.Zperp, _invZperpZperpZperpX₁ , parent.info✻)  # FWL-process
	o.X₂ⱼₖ  = partialjk(parent.X₂, o.Zperp, _invZperpZperpZperpX₂ , parent.info✻)
	o.y₁ⱼₖ  = partialjk(parent.y₁, o.Zperp, _invZperpZperpZperpy₁ , parent.info✻)
	o.Y₂ⱼₖ  = partialjk(parent.Y₂, o.Zperp, _invZperpZperpZperpY₂ , parent.info✻)
	o.Zⱼₖ   = partialjk(  o.Zpar , o.Zperp, _invZperpZperpZperpZ  , parent.info✻)

	tX₁ = ZperpZperp * _invZperpZperpZperpX₁; tX₁ .= ZperpX₁   - tX₁
	tX₂ = ZperpZperp * _invZperpZperpZperpX₂; tX₂ .= ZperpX₂   - tX₂
	tZ  = ZperpZperp * _invZperpZperpZperpZ ; tZ  .= ZperpZpar - tZ
	tY₂ = ZperpZperp * _invZperpZperpZperpY₂; tY₂ .= ZperpY₂   - tY₂
	ty₁ = ZperpZperp * _invZperpZperpZperpy₁; ty₁ .= Zperpy₁  .- ty₁

	_X₁X₁   = panelcross(o.X₁ⱼₖ, o.X₁ⱼₖ, parent.info✻); _X₁X₁    .= o.X₁X₁           - _X₁X₁   ; t✻minus!(_X₁X₁   , ZperpX₁'   , _invZperpZperpZperpX₁); t✻minus!(_X₁X₁    , _invZperpZperpZperpX₁', tX₁)
	_X₂X₁   = panelcross(o.X₂ⱼₖ, o.X₁ⱼₖ, parent.info✻); _X₂X₁    .= o.X₂X₁           - _X₂X₁   ; t✻minus!(_X₂X₁   , ZperpX₂'   , _invZperpZperpZperpX₁); t✻minus!(_X₂X₁    , _invZperpZperpZperpX₂', tX₁)
	_X₂X₂   = panelcross(o.X₂ⱼₖ, o.X₂ⱼₖ, parent.info✻); _X₂X₂    .= o.X₂X₂           - _X₂X₂   ; t✻minus!(_X₂X₂   , ZperpX₂'   , _invZperpZperpZperpX₂); t✻minus!(_X₂X₂    , _invZperpZperpZperpX₂', tX₂)
	_X₁Y₂   = panelcross(o.X₁ⱼₖ, o.Y₂ⱼₖ, parent.info✻); _X₁Y₂    .= o.X₁Y₂           - _X₁Y₂   ; t✻minus!(_X₁Y₂   , ZperpX₁'   , _invZperpZperpZperpY₂); t✻minus!(_X₁Y₂    , _invZperpZperpZperpX₁', tY₂)
	_X₂Y₂   = panelcross(o.X₂ⱼₖ, o.Y₂ⱼₖ, parent.info✻); _X₂Y₂    .= o.X₂Y₂           - _X₂Y₂   ; t✻minus!(_X₂Y₂   , ZperpX₂'   , _invZperpZperpZperpY₂); t✻minus!(_X₂Y₂    , _invZperpZperpZperpX₂', tY₂)
	o.Y₂y₁ⱼₖ = panelcross(o.Y₂ⱼₖ, o.y₁ⱼₖ, parent.info✻); o.Y₂y₁ⱼₖ .= view(o.Y₂y₁,:,:) - o.Y₂y₁ⱼₖ; t✻minus!( o.Y₂y₁ⱼₖ, ZperpY₂'   , _invZperpZperpZperpy₁); t✻minus!( o.Y₂y₁ⱼₖ, _invZperpZperpZperpY₂', ty₁)
	o.X₂y₁ⱼₖ = panelcross(o.X₂ⱼₖ, o.y₁ⱼₖ, parent.info✻); o.X₂y₁ⱼₖ .= view(o.X₂y₁,:,:) - o.X₂y₁ⱼₖ; t✻minus!(o.X₂y₁ⱼₖ , ZperpX₂'   , _invZperpZperpZperpy₁); t✻minus!(o.X₂y₁ⱼₖ , _invZperpZperpZperpX₂', ty₁)
	o.X₁y₁ⱼₖ = panelcross(o.X₁ⱼₖ, o.y₁ⱼₖ, parent.info✻); o.X₁y₁ⱼₖ .= view(o.X₁y₁,:,:) - o.X₁y₁ⱼₖ; t✻minus!(o.X₁y₁ⱼₖ , ZperpX₁'   , _invZperpZperpZperpy₁); t✻minus!(o.X₁y₁ⱼₖ , _invZperpZperpZperpX₁', ty₁)
	o.Zy₁ⱼₖ  = panelcross(o.Zⱼₖ,  o.y₁ⱼₖ, parent.info✻); o.Zy₁ⱼₖ  .= view( o.Zy₁,:,:) - o.Zy₁ⱼₖ ; t✻minus!(o.Zy₁ⱼₖ  , ZperpZpar' , _invZperpZperpZperpy₁); t✻minus!(o.Zy₁ⱼₖ  , _invZperpZperpZperpZ' , ty₁)
	X₁Zⱼₖ    = panelcross(o.X₁ⱼₖ, o.Zⱼₖ , parent.info✻); X₁Zⱼₖ    .= o.X₁Z            - X₁Zⱼₖ   ; t✻minus!(X₁Zⱼₖ    , ZperpX₁'   , _invZperpZperpZperpZ ); t✻minus!(X₁Zⱼₖ    , _invZperpZperpZperpX₁', tZ )
	X₂Zⱼₖ    = panelcross(o.X₂ⱼₖ, o.Zⱼₖ , parent.info✻); X₂Zⱼₖ    .= o.X₂Z            - X₂Zⱼₖ   ; t✻minus!(X₂Zⱼₖ    , ZperpX₂'   , _invZperpZperpZperpZ ); t✻minus!(X₂Zⱼₖ    , _invZperpZperpZperpX₂', tZ )
	o.ZZⱼₖ   = panelcross(o.Zⱼₖ,  o.Zⱼₖ , parent.info✻); o.ZZⱼₖ   .= o.ZZ             - o.ZZⱼₖ  ; t✻minus!(o.ZZⱼₖ   , ZperpZpar' , _invZperpZperpZperpZ ); t✻minus!(o.ZZⱼₖ   , _invZperpZperpZperpZ' , tZ )
	o.ZY₂ⱼₖ  = panelcross(o.Zⱼₖ,  o.Y₂ⱼₖ, parent.info✻); o.ZY₂ⱼₖ  .= o.ZY₂            - o.ZY₂ⱼₖ ; t✻minus!(o.ZY₂ⱼₖ  , ZperpZpar' , _invZperpZperpZperpY₂); t✻minus!(o.ZY₂ⱼₖ  , _invZperpZperpZperpZ' , tY₂)
	o.y₁y₁ⱼₖ = panelcross(o.y₁ⱼₖ, o.y₁ⱼₖ, parent.info✻); o.y₁y₁ⱼₖ .= [o.y₁y₁;;;]     .- o.y₁y₁ⱼₖ; t✻minus!(o.y₁y₁ⱼₖ, 2 * Zperpy₁', _invZperpZperpZperpy₁); o.y₁y₁ⱼₖ .+= _invZperpZperpZperpy₁'ZperpZperp*_invZperpZperpZperpy₁
	o.XZⱼₖ   = [X₁Zⱼₖ; X₂Zⱼₖ]

	o.XY₂ⱼₖ = [_X₁Y₂ ; _X₂Y₂]
	o.XXⱼₖ  = [_X₁X₁ ; _X₂X₁ ;;; _X₂X₁' ;  _X₂X₂]
	o.invXXⱼₖ = invsym(o.XXⱼₖ)
	o.H_2SLSⱼₖ = o.XZⱼₖ'o.invXXⱼₖ * o.XZⱼₖ
	(!isone(o.κ) || o.liml) && (o.H_2SLSmZZⱼₖ = o.H_2SLSⱼₖ - o.ZZⱼₖ)

	if o.restricted
		ZperpZR₁, _ZperpZR₁ = crossjk(o.Zperp, o.ZR₁, parent.info✻)
		_invZperpZperpZperpZR₁ = _invZperpZperp * _ZperpZR₁
		o.ZR₁ⱼₖ = partialjk(o.ZR₁, o.Zperp, _invZperpZperpZperpZR₁, parent.info✻)
		tZR₁ = ZperpZperp * _invZperpZperpZperpZR₁; tZR₁ .= ZperpZR₁ - tZR₁
		o.X₁ZR₁ⱼₖ    = panelcross(o.X₁ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻); o.X₁ZR₁ⱼₖ    .= o.X₁ZR₁ - o.X₁ZR₁ⱼₖ   ; t✻minus!(o.X₁ZR₁ⱼₖ   , ZperpX₁'  , _invZperpZperpZperpZR₁); t✻minus!(o.X₁ZR₁ⱼₖ   , _invZperpZperpZperpX₁' , tZR₁)
		o.X₂ZR₁ⱼₖ    = panelcross(o.X₂ⱼₖ,  o.ZR₁ⱼₖ, parent.info✻); o.X₂ZR₁ⱼₖ    .= o.X₂ZR₁ - o.X₂ZR₁ⱼₖ   ; t✻minus!(o.X₂ZR₁ⱼₖ   , ZperpX₂'  , _invZperpZperpZperpZR₁); t✻minus!(o.X₂ZR₁ⱼₖ   , _invZperpZperpZperpX₂' , tZR₁)
		o.ZR₁Zⱼₖ     = panelcross(o.ZR₁ⱼₖ,   o.Zⱼₖ, parent.info✻); o.ZR₁Zⱼₖ     .= o.ZR₁Z  - o.ZR₁Zⱼₖ    ; t✻minus!( o.ZR₁Zⱼₖ   , ZperpZR₁' , _invZperpZperpZperpZ  ); t✻minus!( o.ZR₁Zⱼₖ   , _invZperpZperpZperpZR₁', tZ) 
		o.twoZR₁y₁ⱼₖ = panelcross(o.ZR₁ⱼₖ,  o.y₁ⱼₖ, parent.info✻); o.twoZR₁y₁ⱼₖ .= reshape(o.ZR₁'parent.y₁,:,1) - o.twoZR₁y₁ⱼₖ; t✻minus!(o.twoZR₁y₁ⱼₖ, ZperpZR₁' , _invZperpZperpZperpy₁ ); t✻minus!(o.twoZR₁y₁ⱼₖ, _invZperpZperpZperpZR₁', ty₁ );  o.twoZR₁y₁ⱼₖ .*= 2
		o.ZR₁ZR₁ⱼₖ   = panelcross(o.ZR₁ⱼₖ, o.ZR₁ⱼₖ, parent.info✻); o.ZR₁ZR₁ⱼₖ   .= o.ZR₁'o.ZR₁     - o.ZR₁ZR₁ⱼₖ  ; t✻minus!(o.ZR₁ZR₁ⱼₖ  , ZperpZR₁' , _invZperpZperpZperpZR₁); t✻minus!(o.ZR₁ZR₁ⱼₖ  , _invZperpZperpZperpZR₁', tZR₁)
		o.ZR₁Y₂ⱼₖ    = panelcross(o.ZR₁ⱼₖ, o.Y₂ⱼₖ,  parent.info✻); o.ZR₁Y₂ⱼₖ    .= o.ZR₁'parent.Y₂ - o.ZR₁Y₂ⱼₖ   ; t✻minus!(o.ZR₁Y₂ⱼₖ   , ZperpZR₁' , _invZperpZperpZperpY₂ ); t✻minus!(o.ZR₁Y₂ⱼₖ   , _invZperpZperpZperpZR₁', tY₂ )
		o.y₁parⱼₖ = Vector{T}(undef, parent.Nobs)
	else
		o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ 
		o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ
		o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ 
		o.Xy₁parⱼₖ     = [o.X₁y₁ⱼₖ ; o.X₂y₁ⱼₖ]
		o.y₁parⱼₖ      = o.y₁ⱼₖ
		ZperpZR₁ = Matrix{T}(undef,0,0)
	end

	!iszero(o.fuller) &&
		(o.Nobsⱼₖ = parent._Nobs .- (parent.fweights ? @panelsum(parent.wt, parent.info✻) : T.(length.(parent.info✻))))
	
	return ZperpX₁, ZperpX₂, Zperpy₁, ZperpY₂, ZperpZpar, ZperpZR₁
end

function InitVarsIV!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}...) where T
	dojkprep   = parent.jk        && o.isDGP && parent.WREnonARubin
	minimizeON = !parent.granular && o.isDGP && parent.WREnonARubin  # form panel-level cross-products to minimize O(N) operations (MacKinnon 2021)?
	needZY₂ = parent.WREnonARubin && !(parent.scorebs && o.isDGP)

	!isempty(Rperp) && (o.Rperp = Rperp[1])

	o.kZ = ncols(o.Rpar)
	o.kZperp = ncols(parent.DGP.RperpX)

	# subtlety for WREnonARubin: in DGP these matrices have not yet been FWL'd; in Repl, most except X₁ have, either because copied from DGP, or because DGP prep overwrote parent data matrices
	if o.copyfromDGP
		o.X₁ = parent.DGP.X₁
		o.X₂ = parent.DGP.X₂
		o.Y₂ = parent.DGP.Y₂
		o.y₁ = parent.DGP.y₁
		o.Zperp = parent.DGP.Zperp
		o.invZperpZperp = parent.DGP.invZperpZperp
	else
		o.X₁ = parent.X₁ * o.RperpXperp
		o.Zperp = parent.X₁ * o.RperpX
	end
	o.kX = (o.kX₁ = ncols(o.X₁)) + parent.kX₂

	if parent.WREnonARubin  # faster, doesn't work for score test because only do FWL in DGP, not Repl
		o.Xpar₁toZparX = DesignerMatrix(parent.DGP.RperpXperp'parent.DGP.RperpXperp \ parent.DGP.RperpXperp'o.RparX)
		o.Zpar = o.X₁ * o.Xpar₁toZparX
		!o.isDGP && t✻minus!(o.Zpar, o.Zperp, parent.DGP.invZperpZperp * (parent.DGP.Zperp'o.Zpar))   # parent.X₁ has not been FWL'd
	else
		o.Zpar = parent.X₁ * o.RparX
	end
	!o.isDGP && !parent.scorebs && (parent.granular || parent.jk) &&
		(o.ZparX = copy(o.Zpar))
	o.Zpar .+= parent.Y₂ * o.RparY
	if o.restricted
		o.ZR₁ = parent.X₁ * o.R₁invR₁R₁X
		!o.isDGP && t✻minus!(o.ZR₁, o.Zperp, parent.DGP.invZperpZperp * (parent.DGP.Zperp'o.ZR₁))   # if Repl restricted (rare) parent.X₁ has not been FWL'd but parent.Y₁ has
		t✻plus!(o.ZR₁, parent.Y₂, o.R₁invR₁R₁Y)
	end

	if minimizeON
		S✻⋂X₁Zpar = panelcross(o.X₁, o.Zpar, parent.info✻⋂)
		S✻⋂X₂Zpar = panelcross(parent.X₂, o.Zpar, parent.info✻⋂)
		o.S✻⋂XZpar = [S✻⋂X₁Zpar; S✻⋂X₂Zpar]
		o.X₁Z = sumpanelcross(S✻⋂X₁Zpar)
		o.X₂Z = sumpanelcross(S✻⋂X₂Zpar)
		o.S✻Zpary₁ = panelcross(o.Zpar, parent.y₁, parent.info✻)
		o.Zy₁ = vec(sumpanelcross(o.S✻Zpary₁))
		o.S✻ZparY₂ = panelcross(o.Zpar, parent.Y₂, parent.info✻)
		o.S✻ZparZpar = panelcross(o.Zpar, o.Zpar, parent.info✻)
		o.ZZ = sumpanelcross(o.S✻ZparZpar)
		needZY₂ && 
			(o.ZY₂ = sumpanelcross(o.S✻ZparY₂))
	
		if o.restricted
			S✻⋂X₁ZR₁ = panelcross(o.X₁, o.ZR₁, parent.info✻⋂)
			S✻⋂X₂ZR₁ = panelcross(parent.X₂, o.ZR₁, parent.info✻⋂)
			o.S✻⋂XZR₁ = [S✻⋂X₁ZR₁ ; S✻⋂X₂ZR₁]
			o.X₁ZR₁    = sumpanelcross(S✻⋂X₁ZR₁)
			o.X₂ZR₁    = sumpanelcross(S✻⋂X₂ZR₁)
			o.S✻ZR₁Z  = panelcross(o.ZR₁, o.Zpar, parent.info✻)
			o.ZR₁Z     = sumpanelcross(o.S✻ZR₁Z)
			o.S✻ZR₁Y₂ = panelcross(o.ZR₁, parent.Y₂, parent.info✻)
			o.ZR₁Y₂    = sumpanelcross(o.S✻ZR₁Y₂) 
			o.S✻ZR₁y₁ = panelcross(o.ZR₁, parent.y₁, parent.info✻)
			o.twoZR₁y₁ = 2 * vec(sumpanelcross(o.S✻ZR₁y₁))
			o.S✻ZR₁ZR₁ = panelcross(o.ZR₁, o.ZR₁, parent.info✻)
			o.ZR₁ZR₁   = sumpanelcross(o.S✻ZR₁ZR₁)
		end
	else
		o.Zy₁ = o.Zpar'parent.y₁
		o.X₁Z = o.X₁'o.Zpar
		o.X₂Z = parent.X₂'o.Zpar
		o.ZZ  = o.Zpar'o.Zpar
		needZY₂ && 
			(o.ZY₂ = o.Zpar'parent.Y₂)
		if o.restricted
			o.X₁ZR₁    = o.X₁'o.ZR₁
			o.X₂ZR₁    = parent.X₂'o.ZR₁
			o.ZR₁Z     = o.ZR₁'o.Zpar
			o.twoZR₁y₁ = o.ZR₁'parent.y₁; lmul!(2, o.twoZR₁y₁)
			o.ZR₁ZR₁   = o.ZR₁'o.ZR₁    
			o.ZR₁Y₂    = o.ZR₁'parent.Y₂
		end
	end

	if o.copyfromDGP
		o.XX = parent.DGP.XX
		o.invXX = parent.DGP.invXX
		o.XY₂  = parent.DGP.XY₂
		o.Y₂Y₂ = parent.DGP.Y₂Y₂
		o.Y₂y₁ = parent.DGP.Y₂y₁
		o.X₂y₁ = parent.DGP.X₂y₁
		o.X₁y₁ = parent.DGP.X₁y₁
		o.y₁y₁ = parent.DGP.y₁y₁
		if !parent.granular
			o.S✻Y₂Y₂     = parent.DGP.S✻Y₂Y₂
			o.S✻⋂X₁Y₂    = parent.DGP.S✻⋂X₁Y₂
			o.S✻⋂X₂Y₂    = parent.DGP.S✻⋂X₂Y₂
			o.S✻⋂ZperpY₂ = parent.DGP.S✻⋂ZperpY₂
			o.S✻Y₂y₁     = parent.DGP.S✻Y₂y₁
		end	
	else
		if minimizeON
			o.S✻⋂X₁Y₂ = panelcross(o.X₁, parent.Y₂, parent.info✻⋂)
			o.S✻⋂X₂Y₂ = panelcross(parent.X₂, parent.Y₂, parent.info✻⋂)
			o.S✻⋂XY₂ = [o.S✻⋂X₁Y₂; o.S✻⋂X₂Y₂]
			S✻⋂X₁X₁ = panelcross(o.X₁, o.X₁, parent.info✻⋂)
			S✻⋂X₂X₁ = panelcross(parent.X₂, o.X₁, parent.info✻⋂)
			S✻⋂X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻⋂)
			o.X₁X₁ = sumpanelcross(S✻⋂X₁X₁)
			o.X₂X₁ = sumpanelcross(S✻⋂X₂X₁)
			o.X₂X₂ = sumpanelcross(S✻⋂X₂X₂)
			o.S✻⋂XX = [[S✻⋂X₁X₁ S✻⋂X₂X₁'] ; [S✻⋂X₂X₁ S✻⋂X₂X₂]]
			o.X₁Y₂ = sumpanelcross(o.S✻⋂X₁Y₂)
			o.X₂Y₂ = sumpanelcross(o.S✻⋂X₂Y₂)
			o.S✻Y₂y₁ = panelcross(parent.Y₂, parent.y₁, parent.info✻)
			o.Y₂y₁ = vec(sumpanelcross(o.S✻Y₂y₁))
			o.S✻Y₂Y₂ = panelcross(parent.Y₂, parent.Y₂, parent.info✻)
			o.Y₂Y₂ = sumpanelcross(o.S✻Y₂Y₂)
			S✻⋂X₂y₁ = panelcross(parent.X₂, parent.y₁, parent.info✻⋂)
			o.X₂y₁ = vec(sumpanelcross(S✻⋂X₂y₁))
			S✻⋂X₁y₁ = panelcross(o.X₁, parent.y₁, parent.info✻⋂)
			o.S✻⋂Xy₁ = [S✻⋂X₁y₁; S✻⋂X₂y₁]
			o.X₁y₁ = vec(sumpanelcross(S✻⋂X₁y₁))
			o.S✻y₁y₁ = panelcross(parent.y₁, parent.y₁, parent.info✻)
			o.y₁y₁ = sum(o.S✻y₁y₁)
			o.S✻⋂ZperpZperp = panelcross(o.Zperp, o.Zperp, parent.info✻⋂)
			S✻⋂X₁Zperp = panelcross(o.X₁, o.Zperp, parent.info✻⋂)
			S✻⋂X₂Zperp = panelcross(parent.X₂, o.Zperp, parent.info✻⋂)
			o.S✻⋂XZperp = [S✻⋂X₁Zperp; S✻⋂X₂Zperp]
			o.S✻⋂ZperpY₂ = panelcross(o.Zperp, parent.Y₂, parent.info✻⋂)
			o.S✻⋂Zperpy₁ = panelcross(o.Zperp, parent.y₁, parent.info✻⋂)
		else
			o.X₁X₁ = o.X₁'o.X₁
			o.X₂X₁ = parent.X₂'o.X₁
			o.X₂X₂ = parent.X₂'parent.X₂
			o.X₁Y₂ = o.X₁'parent.Y₂
			o.X₂Y₂ = parent.X₂'parent.Y₂
			o.Y₂y₁ = parent.Y₂'parent.y₁
			o.Y₂Y₂ = parent.Y₂'parent.Y₂
			o.X₂y₁ = parent.X₂'parent.y₁       
			o.X₁y₁ = o.X₁'parent.y₁
			o.y₁y₁ = dot(parent.y₁, parent.y₁)
		end
	end

	if dojkprep
		o.ZperpX₁, o.ZperpX₂, o.Zperpy₁, o.ZperpY₂, o.ZperpZpar, o.ZperpZR₁ = PrepJKIV!(o, parent)
	else
		if minimizeON
			o.S✻⋂ZperpZpar = panelcross(o.Zperp, o.Zpar, parent.info✻⋂)
			!(parent.jk && parent.WREnonARubin) &&
				(o.ZperpZpar = sumpanelcross(o.S✻⋂ZperpZpar))
			if o.restricted
				o.S✻⋂ZperpZR₁ = panelcross(o.Zperp, o.ZR₁, parent.info✻⋂)
				!(parent.jk && parent.WREnonARubin) &&
					(o.ZperpZR₁ = sumpanelcross(o.S✻⋂ZperpZR₁))
			end
		else
			o.ZperpZpar = o.Zperp'o.Zpar
			o.restricted && (o.ZperpZR₁ = o.Zperp'o.ZR₁)
		end

		if o.copyfromDGP
			o.ZperpX₁ = parent.DGP.ZperpX₁
			o.ZperpX₂ = parent.DGP.ZperpX₂
			o.Zperpy₁ = parent.DGP.Zperpy₁
			o.ZperpY₂ = parent.DGP.ZperpY₂
		elseif minimizeON
			o.ZperpY₂ = sumpanelcross(o.S✻⋂ZperpY₂)
			o.Zperpy₁ = vec(sumpanelcross(o.S✻⋂Zperpy₁))
			o.ZperpX₁ = sumpanelcross(S✻⋂X₁Zperp)'
			o.ZperpX₂ = sumpanelcross(S✻⋂X₂Zperp)'
			o.invZperpZperp = iszero(o.kZperp) ? Matrix{T}(undef,0,0) : invsym(sumpanelcross(o.S✻⋂ZperpZperp))
		else
			o.ZperpX₁ = o.Zperp'o.X₁
			o.ZperpX₂ = o.Zperp'parent.X₂
			o.Zperpy₁ = o.Zperp'parent.y₁
			o.ZperpY₂ = o.Zperp'parent.Y₂
			o.invZperpZperp = invsym(o.Zperp'o.Zperp)
		end
	end

  # now that coarse treatment, if any, is done with data matrices, potentially overwrite them with FWL
	if o.copyfromDGP
		o.invZperpZperpZperpX₁ = parent.DGP.invZperpZperpZperpX₁
		o.invZperpZperpZperpX₂ = parent.DGP.invZperpZperpZperpX₂
		o.invZperpZperpZperpX  = parent.DGP.invZperpZperpZperpX 
		o.invZperpZperpZperpY₂ = parent.DGP.invZperpZperpZperpY₂
		o.invZperpZperpZperpy₁ = parent.DGP.invZperpZperpZperpy₁
	else
		o.invZperpZperpZperpX₁ = o.invZperpZperp * o.ZperpX₁
		o.invZperpZperpZperpX₂ = o.invZperpZperp * o.ZperpX₂
		o.invZperpZperpZperpX  = [o.invZperpZperpZperpX₁ o.invZperpZperpZperpX₂]
		o.invZperpZperpZperpY₂ = o.invZperpZperp * o.ZperpY₂
		o.invZperpZperpZperpy₁ = o.invZperpZperp * o.Zperpy₁

		if !parent.overwrite && parent.granular
			parent.X₂ = copy(parent.X₂)
			parent.y₁ = copy(parent.y₁)
			parent.Y₂ = copy(parent.Y₂)
			parent.overwrite = true
		end
		o.X₂ = parent.X₂
		o.y₁ = parent.y₁
		o.Y₂ = parent.Y₂

		if parent.granular || dojkprep || (parent.NFE>0 && !parent.FEboot && (parent.willfill || parent.not2SLS))
			t✻minus!(o.X₁, o.Zperp, o.invZperpZperpZperpX₁)
			t✻minus!(o.X₂, o.Zperp, o.invZperpZperpZperpX₂)
			t✻minus!(o.y₁, o.Zperp, o.invZperpZperpZperpy₁)
			t✻minus!(o.Y₂, o.Zperp, o.invZperpZperpZperpY₂)
		elseif parent.scorebs
			t✻minus!(o.y₁, o.Zperp, o.invZperpZperpZperpy₁)
		end

		o.X₂X₁ .-= o.ZperpX₂'o.invZperpZperpZperpX₁
		o.X₁X₁ .-= o.ZperpX₁'o.invZperpZperpZperpX₁
		o.X₂X₂ .-= o.ZperpX₂'o.invZperpZperpZperpX₂
		o.X₁Y₂ .-= o.invZperpZperpZperpX₁'o.ZperpY₂
		o.X₂Y₂ .-= o.invZperpZperpZperpX₂'o.ZperpY₂
	  o.Y₂y₁ .-= o.ZperpY₂'o.invZperpZperpZperpy₁
		o.Y₂Y₂ .-= o.ZperpY₂'o.invZperpZperpZperpY₂
		o.X₂y₁ .-= o.ZperpX₂'o.invZperpZperpZperpy₁
		o.X₁y₁ .-= o.ZperpX₁'o.invZperpZperpZperpy₁
		o.y₁y₁  -= o.Zperpy₁'o.invZperpZperpZperpy₁

		o.XY₂ = [o.X₁Y₂; o.X₂Y₂]
		o.XX = [o.X₁X₁ o.X₂X₁' ; o.X₂X₁ o.X₂X₂]
		o.invXX = invsym(o.XX)
	end

	o.invZperpZperpZperpZ = o.invZperpZperp * o.ZperpZpar
	t✻minus!(o.Zpar, o.Zperp, o.invZperpZperpZperpZ )
	o.restricted && t✻minus!(o.ZR₁, o.Zperp, o.invZperpZperp * o.ZperpZR₁)

	t✻minus!(o.Zy₁, o.ZperpZpar', o.invZperpZperpZperpy₁)
	t✻minus!(o.ZZ , o.ZperpZpar', o.invZperpZperpZperpZ)
	t✻minus!(o.X₁Z, o.ZperpX₁'  , o.invZperpZperpZperpZ)
	t✻minus!(o.X₂Z, o.ZperpX₂'  ,o.invZperpZperpZperpZ)
	o.XZ = [o.X₁Z; o.X₂Z]
	needZY₂ && 
		(o.ZY₂ .-= o.ZperpZpar'o.invZperpZperpZperpY₂)

	if o.restricted
		o.invZperpZperpZperpZR₁ = o.invZperpZperp * o.ZperpZR₁
		t✻minus!(o.X₂ZR₁   ,       o.invZperpZperpZperpX₂' , o.ZperpZR₁ )
		t✻minus!(o.X₁ZR₁   ,       o.invZperpZperpZperpX₁' , o.ZperpZR₁ )
		t✻minus!(o.ZR₁Z    ,       o.invZperpZperpZperpZR₁', o.ZperpZpar)
		t✻minus!(o.twoZR₁y₁, T(2), o.invZperpZperpZperpZR₁', o.Zperpy₁  )
		t✻minus!(o.ZR₁ZR₁  ,       o.invZperpZperpZperpZR₁', o.ZperpZR₁ )
		t✻minus!(o.ZR₁Y₂   ,       o.invZperpZperpZperpZR₁', o.ZperpY₂  )
	else
		o.t₁ = zeros(T, parent.kZ)
		parent.jk && (o.t₁Y = o.t₁[parent.kX₁+1:end])
	end

	o.Y₂y₁par    = copy(o.Y₂y₁)  # RHS objects are === across DGP and Repl, but LHS objects can be modified later; break link
	o.X₂y₁par    = copy(o.X₂y₁)
	o.X₁y₁par    = copy(o.X₁y₁)
	o.Zy₁par     = copy(o.Zy₁)
	o.y₁pary₁par = copy(o.y₁y₁)
	o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]

	o.y₁par = copy(o.y₁)

	o.V = o.invXX * o.XZ  # in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and liml robust VCV by Stata convention
	o.H_2SLS = o.XZ'o.V  # Hessian

	if parent.granular || dojkprep
		o.ȳ₁ = Vector{T}(undef, parent.Nobs)
		o.Ȳ₂ = Matrix{T}(undef, parent.Nobs, parent.kY₂)
	end

	if o.isDGP
		(parent.scorebs || parent.jk || parent.granular) &&
			(o.ü₁ = [Vector{T}(undef, parent.Nobs) for _ in 0:parent.jk])
		if !parent.scorebs && (dojkprep || parent.granular)
			o.u⃛₁ = Vector{T}(undef, parent.Nobs)
			o.Ü₂ = Matrix{T}(undef, parent.Nobs, parent.kY₂)
		end

		parent.jk && parent.WREnonARubin &&
			(o.invHⱼₖ = o.liml ? Array{T,3}(undef, o.kZ, parent.N✻, o.kZ) : invsym(isone(o.κ) ? o.H_2SLSⱼₖ : o.ZZⱼₖ + o.κ * o.H_2SLSmZZⱼₖ))

		if o.liml
			o.H_2SLSmZZ = o.H_2SLS - o.ZZ
		else
			MakeH!(o, parent, !isempty(Rperp)) # DGP is liml except possibly when getting confidence peak for A-R plot; but liml=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
		end
	else
		o.Yendog = [true colsum(o.RparY .!= zero(T)).!=0]
		parent.jk && parent.willfill && (o.PXZ = X₁₂B(o.X₁, o.X₂, o.V))
		!parent.scorebs && (parent.granular || parent.jk) &&
			t✻minus!(o.ZparX, o.Zperp, o.invZperpZperp * (o.Zperp'o.ZparX))
	end
	nothing
end


# do most of estimation; for liml r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β̈₀ (β̈ when r=0) and ∂β̈∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function EstimateOLS!(o::StrEstimator, _jk::Bool, r₁::AbstractVector)
	o.β̈ = o.β̈₀ - o.∂β̈∂r * r₁
	_jk && o.restricted &&
		(o.t₁ = o.R₁invR₁R₁ * r₁)
	nothing
end

function EstimateARubin!(o::StrEstimator{T}, parent::StrBootTest{T}, _jk::Bool, r₁::AbstractVector) where T
	EstimateOLS!(o, _jk, r₁)
  o.y₁par .= parent.y₁; t✻minus!(o.y₁par, parent.Y₂, r₁)
	nothing
end

function MakeH!(o::StrEstimator{T}, parent::StrBootTest{T}, makeXAR::Bool=false) where T
  H = isone(o.κ) ? o.H_2SLS : o.ZZ + o.κ * o.H_2SLSmZZ
  o.invH = invsym(H)

  if makeXAR  # for replication regression in score bootstrap of IV/GMM
	  o.A = ncols(o.Rperp)>0 ? (o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp') : o.invH
	  o.AR = o.A * o.Rpar'parent.R'
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
	  if parent.scorebs || parent.granular && o.isDGP || parent.jk && !o.isDGP
			o.y₁par .= o.y₁;
			t✻minus!(o.y₁par, o.ZR₁, r₁)
		end

		if _jk
			o.t₁Y = o.t₁[parent.kX₁+1:end]
			o.y₁parⱼₖ .= o.y₁ⱼₖ; t✻minus!(o.y₁parⱼₖ, o.ZR₁ⱼₖ, r₁)
			o.y₁pary₁parⱼₖ = o.y₁y₁ⱼₖ - o.twoZR₁y₁ⱼₖ'r₁ + r₁'o.ZR₁ZR₁ⱼₖ * r₁
			o.Y₂y₁parⱼₖ    = o.Y₂y₁ⱼₖ - o.ZR₁Y₂ⱼₖ'r₁
			o.Zy₁parⱼₖ     = o.Zy₁ⱼₖ  - o.ZR₁Zⱼₖ'r₁
			  X₁y₁parⱼₖ    = o.X₁y₁ⱼₖ - o.X₁ZR₁ⱼₖ * r₁
			  X₂y₁parⱼₖ    = o.X₂y₁ⱼₖ - o.X₂ZR₁ⱼₖ * r₁
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
	o.ü₁[1] .= o.y₁par; X₁₂Bminus!(o.ü₁[1], parent.X₁, parent.X₂, view(o.β̈ ,:,1))   # ordinary non-jk residuals

	if parent.jk
		m = parent.small ? sqrt((parent.N✻ - 1) / T(parent.N✻)) : one(T)
		if parent.purerobust
			if o.restricted
				X₁₂B!(o.ü₁[2], parent.X₁, parent.X₂, o.t₁)
				o.ü₁[2] .= m .* (o.invMjkv .* (o.ü₁[1] .+ o.ü₁[2]) .- o.ü₁[2])
			else
				o.ü₁[2] .= m .* o.invMjkv .* o.ü₁[1]
			end
		elseif parent.granularjk
	    for (g,S) ∈ enumerate(parent.info✻)
				if o.restricted
					_Xt₁    = rowsubview(o.Xt₁   , length(S))
					_Xt₁pu  = rowsubview(o.Xt₁pu , length(S))
					_MXt₁pu = rowsubview(o.MXt₁pu, length(S))
					X₁₂B!(_Xt₁, view(parent.X₁,S,:), view(parent.X₂,S,:), o.t₁)
					_Xt₁pu .= view(o.ü₁[1], S) .+ _Xt₁
					mul!(_MXt₁pu, o.invMjk[g], _Xt₁pu)
					o.ü₁[2][S] .= m .* (_MXt₁pu .- _Xt₁)
				else
					t✻!(view(o.ü₁[2],S), m, o.invMjk[g], view(o.ü₁[1], S))
				end
			end
		else
	    for (g,S) ∈ enumerate(parent.info✻)
				ü₁g = view(o.ü₁[1],S)
				t✻!(view(o.Xu, 1:parent.kX₁     ), view(parent.X₁,S,:)', ü₁g)
				t✻!(view(o.Xu, parent.kX₁+1:parent.kX), view(parent.X₂,S,:)', ü₁g)
				o.restricted && t✻plus!(o.Xu, view(o.S✻XX,:,g,:), o.t₁)
				ü₂g = view(o.ü₁[2],S)
				mul!(ü₂g, o.XinvHⱼₖ[g], o.Xu)
				ü₂g .= m .* (ü₁g .+ ü₂g)
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
  parent.bootstrapt && parent.robust && parent.NFE>0 && !parent.FEboot && parent.granular < parent.NErrClustCombs &&  # make first factor of second term of (64) for c=⋂ (c=1)
	  (o.CT_XAR = crosstabFE(parent, o.XAR, parent.ID⋂, parent.N⋂))
	nothing
end