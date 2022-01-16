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
	if o.restricted
	  invR₁R₁ = invsym(R₁ * R₁')
	  all(iszero.(diag(invR₁R₁))) && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
	  o.R₁invR₁R₁ = R₁'invR₁R₁
	  F = eigen(Symmetric(o.R₁invR₁R₁ * R₁))
	  o.R₁perp = F.vectors[:, abs.(F.values) .< 1000*eps(T)]  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
  else
	  o.R₁invR₁R₁ = Matrix{T}(undef, parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
	  RR₁perp = Matrix([R
		                  zeros(T, parent.kY₂, parent.kX₁) I])  # rows to prevent partialling out of endogenous regressors
	  nrows(R₁)>0 && (RR₁perp *= o.R₁perp)
	  F = eigen(Symmetric(RR₁perp'pinv(RR₁perp * RR₁perp')*RR₁perp)); val = abs.(F.values) .> 1000*eps(T)
	  o.Rpar   = F.vectors[:, val]  # perp and par of RR₁perp
    o.RperpX = F.vectors[:, abs.(F.values) .< 1000*eps(T)]

	  if nrows(R₁) > 0  # fold model constraint factors into Rpar, RperpX
	    o.Rpar   = o.R₁perp * o.Rpar
	    o.RperpX = o.R₁perp * o.RperpX
	  end

	  o.RRpar = R * o.Rpar
	  o.RperpX = o.RperpX[1:parent.kX₁,:]  # Zperp=X₁*RperpX
		o.RperpXperp = perp(o.RperpX)
		o.RparY = o.Rpar[parent.kX₁+1:end,:]  # part of Rpar that refers to Y₂
	  o.R₁invR₁R₁Y = o.R₁invR₁R₁[parent.kX₁+1:end,:]
	  o.RR₁invR₁R₁ = R * o.R₁invR₁R₁
  end
	nothing
end

# stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
function InitVarsOLS!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}) where T # Rperp is for replication regression--no null imposed
  o.y₁par = parent.y₁
  H = symcross(parent.X₁, parent.wt)
  o.invH = Symmetric(inv(H))

  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.invH : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
  o.β̂₀ = R₁AR₁ * crossvec(parent.X₁, parent.wt, o.y₁par)
  o.∂β̂∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

  o.A = iszero(nrows(Rperp)) ? o.invH : Symmetric(Rperp * invsym(Rperp'H*Rperp) * Rperp')  # for replication regression
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = parent.X₁ * o.AR)
	nothing
end

function InitVarsARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  X₂X₁ = cross(parent.X₂, parent.wt, parent.X₁)
  H = Symmetric([symcross(parent.X₁, parent.wt) X₂X₁' ; X₂X₁ symcross(parent.X₂, parent.wt)])
  o.A = inv(H)
  o.AR = o.A * parent.R'
  (parent.scorebs || parent.robust) && (o.XAR = X₁₂B(parent, parent.X₁, parent.X₂, o.AR))

  R₁AR₁ = iszero(nrows(o.R₁perp)) ? o.A : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')
  o.β̂₀   = R₁AR₁ * [crossvec(parent.X₁, parent.wt, parent.y₁) ; crossvec(parent.X₂, parent.wt, parent.y₁)]
  o.∂β̂∂r = R₁AR₁ * [cross(parent.X₁, parent.wt, parent.Y₂) ; cross(parent.X₂, parent.wt, parent.Y₂)]
	nothing
end

@inline sumpanelsum(X::Array{T,3} where T) = dropdims(sum(X, dims=2), dims=2)

function InitVarsIV!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

  o.Zperp = parent.X₁ * o.RperpX
	o.S⋂ZperpZperp = panelcross11(parent, o.Zperp, o.Zperp, parent.info⋂)  #  XXX exploit symmetry?  XXX missing weights
  o.invZperpZperp = iszero(length(o.Zperp)) ? Symmetric(Matrix{T}(undef,0,0)) : inv(Symmetric(sumpanelsum(o.S⋂ZperpZperp)))
	#=o.X₁=# _X₁ = parent.X₁ * o.RperpXperp
	S⋂X₁Zperp = panelcross11(parent, _X₁, o.Zperp, parent.info⋂)
	ZperpX₁ = sumpanelsum(S⋂X₁Zperp)'
	
	invZperpZperpZperpX₁ = o.invZperpZperp * ZperpX₁
	o.X₁ = _X₁ - o.Zperp * invZperpZperpZperpX₁  # shrink and FWL-process X₁; do it as an O(N) operation because it can be so size-reducing

	S⋂X₂X₁ = panelcross11(parent, parent.X₂, _X₁, parent.info⋂)  # S⋂(X₂:*M_(Z⟂)X₁∥)
	S⋂X₂Zperp = panelcross11(parent, parent.X₂, o.Zperp, parent.info⋂)
	ZperpX₂ = sumpanelsum(S⋂X₂Zperp)'
	o.ZperpX = [ZperpX₁ ZperpX₂]

	invZperpZperpZperpX₂ = o.invZperpZperp * ZperpX₂
	o.X₂ = o.Zperp * invZperpZperpZperpX₂; o.X₂ .= parent.X₂ .- o.X₂  # FWL-process X₂
	o.invZperpZperpZperpX = [invZperpZperpZperpX₁ invZperpZperpZperpX₂]

	X₂X₁ = sumpanelsum(S⋂X₂X₁) - ZperpX₂'o.invZperpZperp * ZperpX₁  # cross(o.X₂, parent.wt, o.X₁)
  X₁X₁ = Symmetric(cross(_X₁, parent.wt, _X₁) - ZperpX₁'o.invZperpZperp * ZperpX₁)
	X₂X₂ = symcross(parent.X₂, parent.wt) - Symmetric(ZperpX₂'o.invZperpZperp * ZperpX₂)
	o.XX = Symmetric([X₁X₁ X₂X₁' ; X₂X₁ X₂X₂])
  o.kX = ncols(o.XX)
  o.invXX = invsym(o.XX)

  #=o.Z=# _Z   = X₁₂B(parent, parent.X₁, parent.Y₂, o.Rpar     )  # Z∥

	S⋂X₁Z = panelcross11(parent,       _X₁, _Z, parent.info⋂)
	S⋂X₂Z = panelcross11(parent, parent.X₂, _Z, parent.info⋂)
	S⋂XZ  = [S⋂X₁Z; S⋂X₂Z]
	X₁Z = sumpanelsum(S⋂X₁Z)
	X₂Z = sumpanelsum(S⋂X₂Z)
	S⋂X₂Zperp = panelcross11(parent, parent.X₂, o.Zperp, parent.info⋂)
	ZperpX₂ = sumpanelsum(S⋂X₂Zperp)'
	o.S⋂XZperp = [S⋂X₁Zperp; S⋂X₂Zperp]
	S⋂ZperpZ = panelcross11(parent, o.Zperp, _Z, parent.info⋂)
  o.ZperpZ = sumpanelsum(S⋂ZperpZ)
	o.invZperpZperpZperpZ = o.invZperpZperp * o.ZperpZ

#	!o.isDGP && ((o.LIML || !parent.robust || !isone(o.κ)) || parent.granular) &&
		(o.Z = _Z - o.Zperp * o.invZperpZperpZperpZ)
	o.ZperpY₂ = cross(o.Zperp, parent.wt, parent.Y₂)
	o.invZperpZperpZperpY₂ = o.invZperpZperp * o.ZperpY₂
  o.Y₂ = parent.Y₂ - o.Zperp * o.invZperpZperpZperpY₂
	o.S⋂Zperpy₁ = panelcross11(parent, o.Zperp, parent.y₁, parent.info⋂)
  Zperpy₁ = reshape(sumpanelsum(o.S⋂Zperpy₁), :)
	o.invZperpZperpZperpy₁ = o.invZperpZperp * Zperpy₁
	o.y₁ = parent.y₁ - o.Zperp * o.invZperpZperpZperpy₁

  o.X₁Y₂ = cross(_X₁, parent.wt, parent.Y₂) - ZperpX₁'o.invZperpZperp * o.ZperpY₂  # cross(o.X₁, parent.wt, o.Y₂)
  o.X₂Y₂ = cross(parent.X₂, parent.wt, parent.Y₂) .- ZperpX₂'o.invZperpZperp * o.ZperpY₂  # cross(o.X₂, parent.wt, o.Y₂)
  o.XY₂ = [o.X₁Y₂ ; o.X₂Y₂]
  o.Y₂y₁ = crossvec(parent.Y₂, parent.wt, parent.y₁) - o.ZperpY₂'o.invZperpZperpZperpy₁  # crossvec(o.Y₂, parent.wt, o.y₁)
	S⋂X₂y₁ = panelcross11(parent, parent.X₂, parent.y₁, parent.info⋂)
  o.X₂y₁ = reshape(sumpanelsum(S⋂X₂y₁), :) - ZperpX₂'o.invZperpZperpZperpy₁  # crossvec(o.X₂, parent.wt, o.y₁)
	S⋂X₁y₁ = panelcross11(parent, _X₁, parent.y₁, parent.info⋂)
  o.X₁y₁ = reshape(sumpanelsum(S⋂X₁y₁), :) - ZperpX₁'o.invZperpZperpZperpy₁  # crossvec(o.X₁, parent.wt, o.y₁)
  o.y₁y₁ = cross(parent.y₁, parent.wt, parent.y₁)[1] -  Zperpy₁'o.invZperpZperpZperpy₁  # cross(o.y₁, parent.wt, o.y₁)[1]
  o.Zy₁  = crossvec(_Z, parent.wt, parent.y₁) - o.ZperpZ'o.invZperpZperpZperpy₁  # crossvec(o.Z, parent.wt, o.y₁)
	o.XZ   = [X₁Z - ZperpX₁'o.invZperpZperpZperpZ ; X₂Z - ZperpX₂'o.invZperpZperpZperpZ]
  o.ZY₂ =  cross(_Z, parent.wt, parent.Y₂) - o.ZperpZ'o.invZperpZperp * o.ZperpY₂  # cross(o.Z, parent.wt, o.Y₂)
  o.ZZ  =  symcross(_Z, parent.wt) - Symmetric(o.ZperpZ'o.invZperpZperpZperpZ)

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ  # this is symmetric but converting to Symmetric() only hampers type inference in the one place it's used

  if o.restricted
		#=o.ZR₁=# _ZR₁ = X₁₂B(parent, parent.X₁, parent.Y₂, o.R₁invR₁R₁)
		S⋂X₁_ZR₁ = panelcross11(parent, _X₁, _ZR₁, parent.info⋂)
		S⋂X₂_ZR₁ = panelcross11(parent, parent.X₂, _ZR₁, parent.info⋂)
		o.S⋂X_ZR₁ = [S⋂X₁_ZR₁; S⋂X₂_ZR₁]
		o.S⋂Zperp_ZR₁ = panelcross11(parent, o.Zperp, _ZR₁, parent.info⋂)
		o.Zperp_ZR₁ = sumpanelsum(o.S⋂Zperp_ZR₁)
		o.invZperpZperpZperpZR₁ = o.invZperpZperp * o.Zperp_ZR₁
		#=o.ZR₁ .-= =# o.ZR₁ = _ZR₁ - o.Zperp * o.invZperpZperpZperpZR₁
	  o.X₁ZR₁    = sumpanelsum(S⋂X₁_ZR₁) - ZperpX₁'o.invZperpZperp * o.Zperp_ZR₁  # cross(o.X₁, parent.wt, o.ZR₁)
	  o.X₂ZR₁    = sumpanelsum(S⋂X₂_ZR₁) - ZperpX₂'o.invZperpZperp * o.Zperp_ZR₁  # cross(o.X₂, parent.wt, o.ZR₁)
	  o.ZZR₁     = cross(       _Z, parent.wt, _ZR₁) - o.ZperpZ'o.invZperpZperp * o.Zperp_ZR₁  # cross(o.Z , parent.wt, o.ZR₁)
	  o.ZR₁Y₂    = cross(_ZR₁, parent.wt, parent.y₁) - o.Zperp_ZR₁'o.invZperpZperpZperpy₁  # cross(o.ZR₁, parent.wt, o.Y₂)
	  o.twoR₁Zy₁ = 2 * (crossvec(_ZR₁, parent.wt, parent.y₁) - o.Zperp_ZR₁'o.invZperpZperpZperpy₁)  # 2crossvec(o.ZR₁, parent.wt, o.y₁)
	  o.ZR₁ZR₁   = symcross(_ZR₁, parent.wt) - Symmetric(o.Zperp_ZR₁'o.invZperpZperp * o.Zperp_ZR₁)  # symcross(o.ZR₁, parent.wt)
o._ZR₁=_ZR₁
  else
	  o.Y₂y₁par    = o.Y₂y₁
	  o.X₂y₁par    = o.X₂y₁
	  o.X₁y₁par    = o.X₁y₁
	  o.Zy₁par     = o.Zy₁
	  o.y₁pary₁par = o.y₁y₁
	  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
	  o.y₁par      = o.y₁
  end
o._X₁=_X₁
o._Z=_Z
  o.V =  o.invXX * o.XZ  # in 2SLS case, StrEstimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and LIML robust VCV by Stata convention
  o.H_2SLS = Symmetric(o.V'o.XZ)  # Hessian
  (o.LIML || o.κ ≠ 1) && (o.H_2SLSmZZ = o.H_2SLS - o.ZZ)

  if o.isDGP
	  !o.LIML && MakeH!(o, parent, !isempty(Rperp)) # DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
  else
	  o.kZ = ncols(o.Rpar)
	  o.Yendog = [true; #=o.RparY.type==identity ? fill(true, o.kZ) :=# mapreduce(v->v.≠0, .|, eachcol(o.RparY))]  # columns of Y = [y₁par Zpar] that are endogenous (normally all)

	  if parent.robust && parent.bootstrapt  # for WRE replication regression, prepare for CRVE
			S⋂X₁y₁ = panelcross11(parent, _X₁, parent.y₁, parent.info⋂)  # S⋂(M_(Z⟂)X₁∥ :* y₁)  XXX missing weights
			S⋂X₂y₁ = panelcross11(parent, parent.X₂, parent.y₁, parent.info⋂)  # S⋂(X₂∥:*y₁)
			o.S⋂Xy₁ = [S⋂X₁y₁; S⋂X₂y₁]

			o.S⋂YX       = Vector{Matrix{T}}(undef, o.kZ+1)
			o.S⋂PXYZperp = Vector{Matrix{T}}(undef, o.kZ+1)

			o.XinvXX = X₁₂B(parent, o.X₁, o.X₂, o.invXX); o.PXZ = X₁₂B(parent, o.X₁, o.X₂, o.invXXXZ)
			if parent.haswt
				o.PXZ .*= parent.wt
				o.XinvXX .*= parent.wt
			end

			if parent.NFE>0
				o.CT_FE⋂PY = Vector{Matrix{T}}(undef, o.kZ+1)  # XXX would array comprehension break type stability?
				@inbounds for i ∈ 1:o.kZ
					o.CT_FE⋂PY[i+1] = crosstabFEt(parent, view(o.PXZ,:,i), parent.info⋂) .* parent.invFEwt
				end
			end

			invZperpZperpZperpX = o.invZperpZperp * [ZperpX₁ ZperpX₂]
			S⋂YX = S⋂XZ - o.S⋂XZperp * o.invZperpZperpZperpZ - invZperpZperpZperpX' * (S⋂ZperpZ - o.S⋂ZperpZperp * o.invZperpZperpZperpZ)  # S⋂(M_Zperp[Z or y₁] .* P_(MZperpX)])

			FillingT₀ = o.invXXXZ' * S⋂YX
			o.FillingT₀ = Matrix{Matrix{T}}(undef, o.kZ+1, o.kZ+1)  # fixed component of groupwise term in sandwich filling  XXX would array comprehension break type stability?
			@inbounds for i ∈ 1:o.kZ, j ∈ 1:o.kZ
				o.FillingT₀[i+1,j+1] = reshape(view(FillingT₀,i,:,j),:,1)
			end

			S⋂PXYZperp = o.invXXXZ' * (o.S⋂XZperp - invZperpZperpZperpX' * o.S⋂ZperpZperp)
			@inbounds for i ∈ 1:o.kZ  # precompute various clusterwise sums
				o.S⋂PXYZperp[i+1] = S⋂PXYZperp[i,:,:]  # S⋂(P_(MZperpX) * Z .* Zperp)
				if !parent.granular
					o.S⋂YX[i+1] = S⋂YX[:,:,i]
				end
	    end

			o.FillingT₀₀ = o.invXXXZ' * (o.S⋂Xy₁ - o.ZperpX' * o.invZperpZperp * o.S⋂Zperpy₁ - (o.S⋂XZperp - o.ZperpX' * o.invZperpZperp * o.S⋂ZperpZperp) * o.invZperpZperpZperpy₁)
			o.restricted &&
				(o.∂FillingT₀∂r = o.invXXXZ' * (o.S⋂X_ZR₁ - o.S⋂XZperp * o.invZperpZperp * o.Zperp_ZR₁ - o.ZperpX' * o.invZperpZperp * o.S⋂Zperp_ZR₁))
			if !parent.granular
				S⋂Xy₁ = [S⋂X₁y₁; S⋂X₂y₁]
				o.S⋂y₁X₀ = dropdims(S⋂Xy₁ - o.S⋂XZperp * o.invZperpZperpZperpy₁ - o.ZperpX' * o.invZperpZperp * (o.S⋂Zperpy₁ - o.S⋂ZperpZperp * o.invZperpZperpZperpy₁); dims=3)
				o.restricted &&
					(o.∂S⋂y₁X∂r = o.S⋂X_ZR₁ - o.S⋂XZperp * o.invZperpZperp * o.Zperp_ZR₁ - o.ZperpX' * o.invZperpZperp * (o.S⋂Zperp_ZR₁ - o.S⋂ZperpZperp * o.invZperpZperp * o.Zperp_ZR₁))
			end
	  end
  end
	nothing
end


# do most of estimation; for LIML r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β̂₀ (β̈ when r=0) and ∂β̂∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function EstimateOLS!(o::StrEstimator{T} where T, r₁::AbstractVector)
  o.β̈ = o.β̂₀ - o.∂β̂∂r * r₁
	nothing
end

function EstimateARubin!(o::StrEstimator{T}, parent::StrBootTest{T}, r₁::AbstractVector) where T
  o.β̈ = o.β̂₀ - o.∂β̂∂r * r₁
  o.y₁par = parent.y₁ - parent.Y₂ * r₁
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
    o.y₁pary₁par = o.y₁y₁ - (o.twoR₁Zy₁'r₁)[1] + r₁'o.ZR₁ZR₁ * r₁
	  o.y₁par   = o.y₁ - o.ZR₁ * r₁
	  o.Y₂y₁par = o.Y₂y₁  - o.ZR₁Y₂'r₁
	  o.X₂y₁par = o.X₂y₁ - o.X₂ZR₁ * r₁
	  o.X₁y₁par = o.X₁y₁ - o.X₁ZR₁ * r₁
	  o.Zy₁par  = o.Zy₁ -  o.ZZR₁ * r₁
	  o.Xy₁par  = [o.X₁y₁par ; o.X₂y₁par]
  end

  o.invXXXy₁par = o.invXX * o.Xy₁par
  o.ZXinvXXXy₁par = o.XZ'o.invXXXy₁par
  o.YY   = Symmetric([o.y₁pary₁par           o.Zy₁par'        ; o.Zy₁par        Matrix(o.ZZ)])
  o.YPXY = Symmetric([o.invXXXy₁par'o.Xy₁par o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par  o.ZXinvXXXZ])

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

	  if parent.robust && parent.bootstrapt  # prepare WRE replication regressions
			tmp = o.restricted ? o.FillingT₀₀ - o.∂FillingT₀∂r * r₁ : o.FillingT₀₀ # @panelsum(parent, o.PXZ, o.y₁par, parent.info⋂)

	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[i+1,1] = reshape(view(tmp,i,:,1),:,1)
	    end
	    if !parent.granular
	  		o.S⋂YX[1] = o.S⋂y₁X₀
				o.restricted && (o.S⋂YX[1] .-= dropdims(o.∂S⋂y₁X∂r * r₁, dims=3))   # S⋂(M_Zperp*y₁ .* P_(MZperpX)])
			end
	  end
  end
	nothing
end

@inline function MakeResidualsOLSARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  o.ü₁ = o.y₁par - X₁₂B(parent, parent.X₁, parent.X₂, o.β̈)
	nothing
end

function MakeResidualsIV!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  o.ü₁ = o.y₁par - o.Z * o.β̈

  if !parent.scorebs
	  _β = [1 ; -o.β̈]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β̈  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  o.Π̂ = invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β̈)' + o.XY₂)
		o.Ü₂ = X₁₂B(parent, o.X₁, o.X₂, o.Π̂); o.Ü₂ .= o.Y₂ .- o.Ü₂

	  o.γ̈ = o.RparY * o.β̈ + o.t₁Y
		o.u⃛₁ = o.Ü₂ * o.γ̈; o.u⃛₁ .+= o.ü₁
  end
	nothing
end


# non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
# since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
# depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then κ from DGP LIML affects Hessian, H.
function InitTestDenoms!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  if parent.bootstrapt && (parent.scorebs || parent.robust)
	  (parent.granular || parent.purerobust) && (o.WXAR = vHadw(o.XAR, parent.wt))

	  if parent.robust && parent.NFE>0 && !(parent.FEboot || parent.scorebs) && parent.granular < parent.NErrClustCombs  # make first factor of second term of (64) for c=⋂ (c=1)
	    !isdefined(o, :WXAR) && (o.WXAR = vHadw(o.XAR, parent.wt))
	    o.CT_XAR = [crosstabFEt(parent, view(o.WXAR,:,d), parent.info⋂) for d ∈ 1:parent.dof]
	  end
  end
	nothing
end
