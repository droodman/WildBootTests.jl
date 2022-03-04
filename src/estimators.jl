# Logically, wild bootstrap tests perform estimation at two stages, once as part of the bootstrap DGP, once in each bootstrap replication
# The StrEstimator "class" and its three "children" hold the estimation logic for the OLS, Anderson-Rubin, and IV/GMM cases

function perp(A::AbstractMatrix)
  F = eigen(Symmetric(A*invsym(A'A)*A'))
  F.vectors[:, abs.(F.values) .< 1000eps(eltype(A))]
end

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 nrows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T}, parent::StrBootTest{T}, R₁::AbstractMatrix{T}, R::Union{UniformScaling{Bool},AbstractMatrix{T}}=Matrix{T}(undef,0,0)) where {T,E}
  o.restricted = nrows(R₁) > 0
	if o.restricted > 0
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

		if !(o.isDGP && parent.WREnonARubin)  # DGP regressions will just copy the replication-regression Zperp, X2par for speed
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
	o.ü₁ = Vector{T}(undef, parent.Nobs)
  H = Symmetric(parent.X₁'parent.X₁)
  o.invH = Symmetric(inv(H))
  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.invH : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
  o.β̈₀ = R₁AR₁ * (parent.X₁'o.y₁par)
  o.∂β̈∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

  o.A = iszero(nrows(Rperp)) ? o.invH : Symmetric(Rperp * invsym(Rperp'H*Rperp) * Rperp')  # for replication regression
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = parent.X₁ * o.AR)
	nothing
end

function InitVarsARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
	o.y₁par = Vector{T}(undef, parent.Nobs)
	o.ü₁    = Vector{T}(undef, parent.Nobs)

	X₂X₁ = parent.X₂'parent.X₁
  H = Symmetric([parent.X₁'parent.X₁ X₂X₁' ; X₂X₁ parent.X₂'parent.X₂])
  o.A = inv(H)
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = X₁₂B(parent, parent.X₁, parent.X₂, o.AR))

  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.A : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')
  o.β̈₀   = R₁AR₁ * [parent.X₁'parent.y₁ ; parent.X₂'parent.y₁]
  o.∂β̈∂r = R₁AR₁ * [parent.X₁'parent.Y₂ ; parent.X₂'parent.Y₂]
	nothing
end

function InitVarsIV!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

	if o.isDGP && parent.WREnonARubin
		o.Zperp = parent.Repl.Zperp
		o.Xpar₁ = parent.Repl.Xpar₁
		isdefined(parent.Repl, :X₁) && (o.X₁ = parent.Repl.X₁)
		isdefined(parent.Repl, :X₂) && (o.X₂ = parent.Repl.X₂)
		isdefined(parent.Repl, :Y₂) && (o.Y₂ = parent.Repl.Y₂)
		isdefined(parent.Repl, :y₁) && (o.y₁ = parent.Repl.y₁)
		o.invZperpZperp = parent.Repl.invZperpZperp
		o.XX = parent.Repl.XX
		o.kX = parent.Repl.kX
		o.invXX = parent.Repl.invXX
		o.XY₂ = parent.Repl.XY₂
		o.Y₂Y₂ = parent.Repl.Y₂Y₂
		o.S✻Y₂Y₂ = parent.Repl.S✻Y₂Y₂
		o.Y₂y₁ = parent.Repl.Y₂y₁
		o.X₂y₁ = parent.Repl.X₂y₁
		o.X₁y₁ = parent.Repl.X₁y₁
		o.y₁y₁ = parent.Repl.y₁y₁
		o.invZperpZperpZperpX₁ = parent.Repl.invZperpZperpZperpX₁
		o.invZperpZperpZperpX₂ = parent.Repl.invZperpZperpZperpX₂
		o.invZperpZperpZperpy₁ = parent.Repl.invZperpZperpZperpy₁
		o.invZperpZperpZperpY₂ = parent.Repl.invZperpZperpZperpY₂
		o.S✻⋂X₁Y₂ = parent.Repl.S✻⋂X₁Y₂
		o.S✻⋂X₂Y₂ = parent.Repl.S✻⋂X₂Y₂
		o.S✻⋂ZperpY₂ = parent.Repl.S✻⋂ZperpY₂
		o.S✻Y₂y₁ = parent.Repl.S✻Y₂y₁
	else
		o.kZperp = ncols(o.RperpX)
		o.Zperp = parent.X₁ * o.RperpX
		o.S✻⋂ZperpZperp = panelcross(o.Zperp, o.Zperp, parent.info✻⋂)
		o.invZperpZperp = iszero(ncols(o.RperpX)) ? Symmetric(Matrix{T}(undef,0,0)) : invsym(sumpanelcross(o.S✻⋂ZperpZperp))

		o.Xpar₁ = parent.X₁ * o.RperpXperp
		S✻⋂X₁Zperp = panelcross(o.Xpar₁, o.Zperp, parent.info✻⋂)
		ZperpX₁ = sumpanelcross(S✻⋂X₁Zperp)'
		o.invZperpZperpZperpX₁ = o.invZperpZperp * ZperpX₁

		S✻⋂X₂Zperp = panelcross(parent.X₂, o.Zperp, parent.info✻⋂)
		ZperpX₂ = sumpanelcross(S✻⋂X₂Zperp)'
		o.ZperpX = [ZperpX₁ ZperpX₂]
		o.S✻⋂XZperp = [S✻⋂X₁Zperp; S✻⋂X₂Zperp]

		o.invZperpZperpZperpX₂ = o.invZperpZperp * ZperpX₂
		o.invZperpZperpZperpX = [o.invZperpZperpZperpX₁ o.invZperpZperpZperpX₂]

		if parent.NFE>0 && (parent.LIML || !isone(parent.κ) || parent.bootstrapt)  ||  parent.granular && parent.robust && parent.bootstrapt  ||  !o.LIML && !isempty(Rperp)
			o.X₁ = o.Xpar₁ - o.Zperp * o.invZperpZperpZperpX₁  # shrink and FWL-process X₁; do it as an O(N) operation because it can be so size-reducing
			o.X₂ = o.Zperp * o.invZperpZperpZperpX₂; o.X₂ .= parent.X₂ .- o.X₂  # FWL-process X₂
		end
		
		S✻⋂X₁X₁ = panelcross(o.Xpar₁, o.Xpar₁, parent.info✻⋂)  # S⋂(X₂:*M_(Z⟂)X₁∥)
		S✻⋂X₂X₁ = panelcross(parent.X₂, o.Xpar₁, parent.info✻⋂)
		S✻⋂X₂X₂ = panelcross(parent.X₂, parent.X₂, parent.info✻⋂)
		X₂X₁ = sumpanelcross(S✻⋂X₂X₁) - ZperpX₂'o.invZperpZperp * ZperpX₁
		X₁X₁ = Symmetric(sumpanelcross(S✻⋂X₁X₁)) - Symmetric(ZperpX₁'o.invZperpZperp * ZperpX₁)
		X₂X₂ = Symmetric(sumpanelcross(S✻⋂X₂X₂)) - Symmetric(ZperpX₂'o.invZperpZperp * ZperpX₂)
		o.S✻⋂XX = [[S✻⋂X₁X₁ S✻⋂X₂X₁'] ; [S✻⋂X₂X₁ S✻⋂X₂X₂]]  # [a b; c d] syntax would call hvcat() to concatenate horizontally along dim 2 rather than 3
		o.XX = Symmetric([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
		o.kX = ncols(o.XX)
		o.invXX = invsym(o.XX)
	
		o.S✻⋂ZperpY₂ = panelcross(o.Zperp, parent.Y₂, parent.info✻⋂)
		ZperpY₂ = sumpanelcross(o.S✻⋂ZperpY₂)
		o.invZperpZperpZperpY₂ = o.invZperpZperp * ZperpY₂
		((parent.NFE>0 && (parent.LIML || !isone(parent.κ) || parent.bootstrapt)) || (parent.robust && parent.bootstrapt && parent.granular)) &&
			(o.Y₂ = parent.Y₂ - o.Zperp * o.invZperpZperpZperpY₂)
		o.S✻⋂Zperpy₁ = panelcross(o.Zperp, parent.y₁, parent.info✻⋂)
		Zperpy₁ = sumpanelcross(o.S✻⋂Zperpy₁)
		o.invZperpZperpZperpy₁ = o.invZperpZperp * Zperpy₁
		((parent.NFE>0 && (parent.LIML || !isone(parent.κ) || parent.bootstrapt)) || (parent.scorebs || parent.robust && parent.bootstrapt && parent.granular)) &&
			(o.y₁ = parent.y₁ - o.Zperp * o.invZperpZperpZperpy₁)

	  o.S✻⋂X₁Y₂ = panelcross(o.Xpar₁, parent.Y₂, parent.info✻⋂)
	  o.S✻⋂X₂Y₂ = panelcross(parent.X₂, parent.Y₂, parent.info✻⋂)
		o.S✻⋂XY₂ = [o.S✻⋂X₁Y₂; o.S✻⋂X₂Y₂]
		o.XY₂ = sumpanelcross(o.S✻⋂XY₂) - o.invZperpZperpZperpX'ZperpY₂
		o.S✻Y₂y₁ = panelcross(parent.Y₂, parent.y₁, parent.info✻)
	  o.Y₂y₁ = sumpanelcross(o.S✻Y₂y₁) - ZperpY₂'o.invZperpZperpZperpy₁
		o.S✻Y₂Y₂ = panelcross(parent.Y₂, parent.Y₂, parent.info✻)
		o.Y₂Y₂ = Symmetric(sumpanelcross(o.S✻Y₂Y₂)) - Symmetric(ZperpY₂'o.invZperpZperpZperpY₂)
		S✻⋂X₂y₁ = panelcross(parent.X₂, parent.y₁, parent.info✻⋂)
		o.X₂y₁ = reshape(sumpanelcross(S✻⋂X₂y₁), :) - ZperpX₂'o.invZperpZperpZperpy₁
		S✻⋂X₁y₁ = panelcross(o.Xpar₁, parent.y₁, parent.info✻⋂)
		o.S✻⋂Xy₁ = [S✻⋂X₁y₁; S✻⋂X₂y₁]
		o.X₁y₁ = reshape(sumpanelcross(S✻⋂X₁y₁), :) - ZperpX₁'o.invZperpZperpZperpy₁
		o.S✻y₁y₁ = reshape(panelcross(parent.y₁, parent.y₁, parent.info✻), :)
		o.y₁y₁ = sum(o.S✻y₁y₁) - Zperpy₁'o.invZperpZperpZperpy₁
	end

  o.Z = X₁₂B(parent, parent.X₁, parent.Y₂, o.Rpar)  # Z∥

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
	o.Zy₁ = sumpanelcross(S✻X₁pary₁) + o.RparY' * o.Y₂y₁ - ZperpX₁par'o.invZperpZperpZperpy₁

	S✻X₁parY₂ = panelcross(X₁par, parent.Y₂, parent.info✻)
	o.XZ = [X₁Zpar - o.invZperpZperpZperpX₁'ZperpZpar ; X₂Zpar - o.invZperpZperpZperpX₂'ZperpZpar]
  o.S✻ZparY₂ = S✻X₁parY₂ + o.RparY' * o.S✻Y₂Y₂
  o.ZY₂ = sumpanelcross(S✻X₁parY₂) - ZperpX₁par'o.invZperpZperpZperpY₂
  tmp = S✻X₁parY₂ * o.RparY; o.S✻ZparZpar = panelcross(X₁par, X₁par, parent.info✻) + tmp + tmp' + o.RparY' * o.S✻Y₂Y₂ * o.RparY
  o.ZZ = Symmetric(sumpanelcross(o.S✻ZparZpar)) - Symmetric(ZperpZpar'o.invZperpZperpZperpZpar)

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ  # symmetric but converting to Symmetric() hampers type inference in the one place it's used

  if o.restricted
		_ZR₁ = X₁₂B(parent, parent.X₁, parent.Y₂, o.R₁invR₁R₁)
		S✻⋂X₁ZR₁ = panelcross(o.Xpar₁, _ZR₁, parent.info✻⋂)
		S✻⋂X₂ZR₁ = panelcross(parent.X₂, _ZR₁, parent.info⋂)
		o.S✻⋂XZR₁ = [S✻⋂X₁ZR₁ ; S✻⋂X₂ZR₁]
		o.S✻⋂ZperpZR₁ = panelcross(o.Zperp, _ZR₁, parent.info✻⋂)
		o.ZperpZR₁ = sumpanelcross(o.S✻⋂ZperpZR₁)
		o.invZperpZperpZperpZR₁ = o.invZperpZperp * o.ZperpZR₁
		o.ZR₁ = _ZR₁ - o.Zperp * o.invZperpZperpZperpZR₁
	  o.X₁ZR₁    = sumpanelcross(S✻⋂X₁ZR₁) - o.invZperpZperpZperpX₁'o.ZperpZR₁
	  o.X₂ZR₁    = sumpanelcross(S✻⋂X₂ZR₁) - o.invZperpZperpZperpX₂'o.ZperpZR₁
		o.S✻ZR₁Z   = panelcross(_ZR₁, o.Z, parent.info✻)
	  o.ZR₁Z     = sumpanelcross(o.S✻ZR₁Z) - o.ZperpZR₁'o.invZperpZperp * ZperpZpar
		o.S✻ZR₁Y₂  = panelcross(_ZR₁, parent.Y₂, parent.info✻)
	  o.ZR₁Y₂    = sumpanelcross(o.S✻ZR₁Y₂) - o.ZperpZR₁'o.invZperpZperpZperpY₂
		o.S✻ZR₁y₁  = panelcross(_ZR₁, parent.y₁, parent.info✻)
	  o.twoZR₁y₁ = 2 * (sumpanelcross(o.S✻ZR₁y₁) - o.ZperpZR₁'o.invZperpZperpZperpy₁)
		o.S✻ZR₁ZR₁ = panelcross(_ZR₁, _ZR₁, parent.info✻)
	  o.ZR₁ZR₁   = Symmetric(sumpanelcross(o.S✻ZR₁ZR₁)) - Symmetric(o.ZperpZR₁'o.invZperpZperp * o.ZperpZR₁)
  else
	  o.Y₂y₁par    = o.Y₂y₁
	  o.X₂y₁par    = o.X₂y₁
	  o.X₁y₁par    = o.X₁y₁
	  o.Zy₁par     = o.Zy₁
	  o.y₁pary₁par = o.y₁y₁
	  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
  end
	(parent.scorebs || parent.robust && parent.bootstrapt && parent.granular) && 
		(o.y₁par = copy(o.y₁))

	if o.isDGP
		if parent.scorebs
			o.ü₁ = Vector{T}(undef, parent.Nobs)
		elseif parent.robust && parent.bootstrapt && parent.granular
			o.Ü₂ = Matrix{T}(undef, parent.Nobs, parent.kY₂)
			o.u⃛₁ = Vector{T}(undef, parent.Nobs)
		end
	end
	
	o.Z .-= o.Zperp * o.invZperpZperpZperpZpar

  o.V =  o.invXX * o.XZ  # in 2SLS case, estimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and LIML robust VCV by Stata convention
  o.H_2SLS = Symmetric(o.V'o.XZ)  # Hessian
  (o.LIML || !isone(o.κ)) && (o.H_2SLSmZZ = o.H_2SLS - o.ZZ)

  if o.isDGP
		!o.LIML && MakeH!(o, parent, !isempty(Rperp)) # DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
	else
	  o.kZ = ncols(o.Rpar)
	  o.Yendog = [true colsum(o.RparY .!= zero(T)).!=0]  # columns of Y = [y₁par Zpar] that are endogenous (normally all)
	end
	nothing
end


# do most of estimation; for LIML r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β̈₀ (β̈ when r=0) and ∂β̈∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function EstimateOLS!(o::StrEstimator{T} where T, r₁::AbstractVector)
  o.β̈ = o.β̈₀ - o.∂β̈∂r * r₁
	nothing
end

function EstimateARubin!(o::StrEstimator{T}, parent::StrBootTest{T}, r₁::AbstractVector) where T
  o.β̈ = o.β̈₀ - o.∂β̈∂r * r₁
  o.y₁par .= parent.y₁ .- parent.Y₂ * r₁
	nothing
end

function MakeH!(o::StrEstimator{T}, parent::StrBootTest{T}, makeXAR::Bool=false) where T
  H = isone(o.κ) ? o.H_2SLS : o.ZZ + o.κ * o.H_2SLSmZZ
  o.invH = invsym(H)
  if makeXAR  # for replication regression in score bootstrap of IV/GMM
	  o.A = ncols(o.Rperp)>0 ? Symmetric(o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp') : o.invH
	  o.AR = o.A * (o.Rpar'parent.R')
	  o.XAR = X₁₂B(parent, o.X₁, o.X₂, o.V * o.AR)
  end
	nothing
end

function EstimateIV!(o::StrEstimator{T}, parent::StrBootTest{T}, r₁::AbstractVector) where T
  if o.restricted
    o.y₁pary₁par = o.y₁y₁ - (o.twoZR₁y₁'r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
	  o.Y₂y₁par = o.Y₂y₁  - o.ZR₁Y₂'r₁
	  o.X₂y₁par = o.X₂y₁ - o.X₂ZR₁ * r₁
	  o.X₁y₁par = o.X₁y₁ - o.X₁ZR₁ * r₁
	  o.Zy₁par  = o.Zy₁ -  o.ZR₁Z'r₁
	  o.Xy₁par  = [o.X₁y₁par ; o.X₂y₁par]
	  (parent.scorebs || parent.robust && parent.bootstrapt && parent.granular) && 
			(o.y₁par .= o.y₁ .- o.ZR₁ * r₁)
  end

  o.invXXXy₁par = o.invXX * o.Xy₁par
  o.ZXinvXXXy₁par = o.XZ'o.invXXXy₁par
  o.YY   = Symmetric([[o.y₁pary₁par          ] o.Zy₁par'        ; o.Zy₁par        Matrix(o.ZZ)])
  o.YPXY = Symmetric([[o.invXXXy₁par'o.Xy₁par] o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par  o.ZXinvXXXZ])

  if o.isDGP
	  if o.LIML
	    o.κ = 1/(1 - real(eigvals(invsym(o.YY) * o.YPXY)[1]))  # like Fast & Wild (81), but more stable, at least in Mata
	    !iszero(o.Fuller) && (o.κ -= o.Fuller / (parent._Nobs - parent.kX))
	    MakeH!(o, parent)
	  end

	  o.β̈ = o.invH * (isone(o.κ) ? o.ZXinvXXXy₁par : o.κ * (o.ZXinvXXXy₁par - o.Zy₁par) + o.Zy₁par)
	  o.t₁Y = o.R₁invR₁R₁Y * r₁
  elseif parent.WREnonARubin  # if not score bootstrap of IV/GMM
	  o.Rt₁ = o.RR₁invR₁R₁ * r₁
  end
	nothing
end


@inline function MakeResidualsOLSARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  o.ü₁ .= o.y₁par .- X₁₂B(parent, parent.X₁, parent.X₂, o.β̈)
	nothing
end

function MakeResidualsIV!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  if parent.scorebs
		o.ü₁ .= o.y₁par .- o.Z * o.β̈
	else
		_β = [1 ; -o.β̈]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β̈  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  o.Π̂ = invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β̈)' + o.XY₂)
		o.γ̈ = o.RparY * o.β̈ + o.t₁Y
	  if parent.robust && parent.bootstrapt && parent.granular
			o.Ü₂ .= o.Y₂ .- X₁₂B(parent, o.X₁, o.X₂, o.Π̂)
			o.u⃛₁ .= o.y₁par .- o.Z * o.β̈ .+ o.Ü₂ * o.γ̈
		end
  end
	nothing
end


# non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
# since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
# depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then κ from DGP LIML affects Hessian, H.
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
