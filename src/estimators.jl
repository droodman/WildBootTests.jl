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
  if nrows(R₁) > 0
	  invR₁R₁ = invsym(R₁ * R₁')
	  all(iszero.(diag(invR₁R₁))) && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
	  o.R₁invR₁R₁ = R₁'invR₁R₁
	  F = eigen(Symmetric(o.R₁invR₁R₁ * R₁))
	  o.R₁perp = F.vectors[:, abs.(F.values) .< 1000*eps(T)]  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
  else
	  o.R₁invR₁R₁ = Matrix{T}(undef, parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
	  RR₁perp = Matrix([R ; zeros(T, parent.kY₂, parent.kX₁) I])  # nrows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

	  length(R₁)>0 && (RR₁perp *= o.R₁perp)

	  F = eigen(Symmetric(RR₁perp'pinv(RR₁perp * RR₁perp')*RR₁perp)); val = abs.(F.values) .> 1000*eps(T)
	  o.Rpar   = F.vectors[:,   val]  # perp and par of RR₁perp
    o.RperpX = F.vectors[:, .!val]

	  if nrows(R₁) > 0  # fold model constraint factors into Rpar, RperpX
	    o.Rpar   = o.R₁perp * o.Rpar
	    o.RperpX = o.R₁perp * o.RperpX
	  end

	  o.RRpar = R * o.Rpar
	  o.RperpX = o.RperpX[1:parent.kX₁,:]  # Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X₁, since all endogenous ones will be retained
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

function InitVarsIV!(o::StrEstimator{T}, parent::StrBootTest{T}, Rperp::AbstractMatrix{T}...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

  o.Zperp = parent.X₁ * o.RperpX
  o.invZperpZperp = iszero(length(o.Zperp)) ? Symmetric(Matrix{T}(undef,0,0)) : inv(symcross(o.Zperp, parent.wt))

  o.X₁ = parent.X₁ * o.RperpXperp; o.X₁ .-= o.Zperp * (o.invZperpZperp * cross(o.Zperp, parent.wt, o.X₁))  # FWL-process X₁
	o.X₂ = o.Zperp * (o.invZperpZperp * cross(o.Zperp, parent.wt, parent.X₂)); o.X₂ .= parent.X₂ .- o.X₂   # FWL-process X₂

	X₂X₁ = cross(o.X₂, parent.wt, o.X₁)
  o.XX = Symmetric([symcross(o.X₁, parent.wt) X₂X₁' ; X₂X₁ symcross(o.X₂, parent.wt)])
  o.kX = ncols(o.XX)
  o.invXX = invsym(o.XX)
  o.Z   = X₁₂B(parent, parent.X₁, parent.Y₂, o.Rpar     )  # Z∥
  o.ZR₁ = X₁₂B(parent, parent.X₁, parent.Y₂, o.R₁invR₁R₁)

  o.Z   .-= o.Zperp * (o.invZperpZperp * cross(o.Zperp, parent.wt, o.Z))  # partialling out
  o.ZR₁ .-= o.Zperp * (o.invZperpZperp * cross(o.Zperp, parent.wt, o.ZR₁))
  o.Y₂ = parent.Y₂ - o.Zperp * (o.invZperpZperp * cross(o.Zperp, parent.wt, parent.Y₂))
  o.y₁ = parent.y₁ - o.Zperp * o.invZperpZperp * (crossvec(o.Zperp, parent.wt, parent.y₁))

  o.X₁Y₂ = cross(o.X₁, parent.wt, o.Y₂)
  o.X₂Y₂ = cross(o.X₂, parent.wt, o.Y₂)
  o.XY₂ = [o.X₁Y₂ ; o.X₂Y₂]
  o.Y₂y₁ = crossvec(o.Y₂, parent.wt, o.y₁)
  o.X₂y₁ = crossvec(o.X₂, parent.wt, o.y₁)
  o.X₁y₁ = crossvec(o.X₁, parent.wt, o.y₁)
  o.y₁y₁ = cross(o.y₁ , parent.wt, o.y₁)[1]
  o.Zy₁  = crossvec(o.Z, parent.wt, o.y₁)
  o.XZ   = [cross(o.X₁, parent.wt, o.Z) ;
			      cross(o.X₂, parent.wt, o.Z)]
  o.ZY₂ =  cross(o.Z, parent.wt, o.Y₂)
  o.ZZ  =  symcross(o.Z, parent.wt)

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ  # this is symmetric but converting to Symmetric() only hampers type inference in the one place it's used

  if ncols(o.R₁invR₁R₁)>0
	  o.X₂ZR₁    = cross(o.X₂, parent.wt, o.ZR₁)
	  o.X₁ZR₁    = cross(o.X₁, parent.wt, o.ZR₁)
	  o.ZZR₁     = cross(o.Z , parent.wt, o.ZR₁)
	  o.twoR₁Zy₁ = 2crossvec(o.ZR₁, parent.wt, o.y₁)
	  o.ZR₁ZR₁   = symcross(o.ZR₁, parent.wt)
	  o.ZR₁Y₂    = cross(o.ZR₁, parent.wt, o.Y₂)
  else
	  o.Y₂y₁par    = o.Y₂y₁
	  o.X₂y₁par    = o.X₂y₁
	  o.X₁y₁par    = o.X₁y₁
	  o.Zy₁par     = o.Zy₁
	  o.y₁pary₁par = o.y₁y₁
	  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
	  o.y₁par      = o.y₁
  end

  o.V =  o.invXX * o.XZ  # in 2SLS case, StrEstimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and LIML robust VCV by Stata convention
  o.H_2SLS = Symmetric(o.V'o.XZ)  # Hessian
  (o.LIML || o.κ ≠ 1) && (o.H_2SLSmZZ = o.H_2SLS - o.ZZ)

  if o.isDGP
	  !o.LIML && MakeH!(o, parent, !isempty(Rperp)) # DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
  else
	  o.kZ = ncols(o.Rpar)
	  o.Yendog = [true; #=o.RparY.type==identity ? fill(true, o.kZ) :=# mapreduce(v->v.≠0, .|, eachcol(o.RparY))]  # columns of Y = [y₁par Zpar] that are endogenous (normally all)

	  if parent.robust && parent.bootstrapt  # for WRE replication regression, prepare for CRVE
			o.S⋂YX       = Vector{Matrix{T}}(undef, o.kZ+1)
			o.S⋂PXYZperp = Vector{Matrix{T}}(undef, o.kZ+1)

			o.XinvXX = X₁₂B(parent, o.X₁, o.X₂, o.invXX); o.PXZ = X₁₂B(parent, o.X₁, o.X₂, o.invXXXZ)
			if parent.haswt
				o.PXZ .*= parent.wt
				o.XinvXX .*= parent.wt
			end

			o.FillingT₀ = Matrix{Matrix{T}}(undef, o.kZ+1, o.kZ+1)  # fixed component of groupwise term in sandwich filling
			parent.NFE>0 &&
				(o.CT_FE⋂PY = Vector{Matrix{T}}(undef, o.kZ+1))
			@inbounds for i ∈ 1:o.kZ
				u = view(o.PXZ,:,i)
				parent.NFE>0 &&
					(o.CT_FE⋂PY[i+1] = crosstabFEt(parent, u, parent.info⋂Data) .* parent.invFEwt)
				tmp = @panelsum(parent, o.Z, u, parent.info⋂Data)
				for j ∈ 1:o.kZ
					o.FillingT₀[i+1,j+1] = reshape(view(tmp,:,j),:,1)
				end
			end
			@inbounds for i ∈ 1:o.kZ  # precompute various clusterwise sums
				o.S⋂PXYZperp[i+1] = @panelsum(parent, o.Zperp, view(o.PXZ,:,i), parent.info⋂Data)  # S⋂(P_(MZperpX) * Z .* Zperp)
				if !parent.granular
					if parent.haswt
						(o.S⋂YX[i+1] = @panelsum2(parent, o.X₁, o.X₂, view(o.Z,:,i) .* parent.wt, parent.info⋂Data))  # S⋂(M_Zperp[Z or y₁] .* P_(MZperpX)])
					else
						(o.S⋂YX[i+1] = @panelsum2(parent, o.X₁, o.X₂, view(o.Z,:,i), parent.info⋂Data))  # S⋂(M_Zperp[Z or y₁] .* P_(MZperpX)])
					end
				end
	    end
	  end
  end
	nothing
end


# do most of estimation; for LIML r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β̂₀ (β̂ when r=0) and ∂β̂∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function EstimateOLS!(o::StrEstimator{T} where T, r₁::AbstractVector)
  o.β̂ = o.β̂₀ - o.∂β̂∂r * r₁
	nothing
end

function EstimateARubin!(o::StrEstimator{T}, parent::StrBootTest{T}, r₁::AbstractVector) where T
  o.β̂ = o.β̂₀ - o.∂β̂∂r * r₁
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
  if ncols(o.R₁invR₁R₁)>0
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
  o.YY   = Symmetric([[o.y₁pary₁par          ] o.Zy₁par'        ; o.Zy₁par        Matrix(o.ZZ)])
  o.YPXY = Symmetric([[o.invXXXy₁par'o.Xy₁par] o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par  o.ZXinvXXXZ])

  if o.isDGP
	  if o.LIML
	    o.κ = 1/(1 - real(eigvals(invsym(o.YY) * o.YPXY)[1]))  # like Fast & Wild (81), but more stable, at least in Mata
	    !iszero(o.Fuller) && (o.κ -= o.Fuller / (parent._Nobs - parent.kX))
	    MakeH!(o, parent)
	  end

	  o.β̂ = o.invH * (isone(o.κ) ? o.ZXinvXXXy₁par : o.κ * (o.ZXinvXXXy₁par - o.Zy₁par) + o.Zy₁par)
	  o.t₁Y = o.R₁invR₁R₁Y * r₁
  elseif parent.WREnonARubin  # if not score bootstrap of IV/GMM
	  o.Rt₁ = o.RR₁invR₁R₁ * r₁

	  if parent.robust && parent.bootstrapt  # prepare WRE replication regressions
	    tmp = @panelsum(parent, o.PXZ, o.y₁par, parent.info⋂Data)
	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[i+1,1] = reshape(view(tmp,:,i),:,1)
	    end
	    !parent.granular &&
	  		(o.S⋂YX[1] = @panelsum2(parent, o.X₁, o.X₂, vHadw(o.y₁par, parent.wt), parent.info⋂Data))  # S⋂(M_Zperp*y₁ .* P_(MZperpX)])
	  end
  end
	nothing
end


@inline function MakeResidualsOLSARubin!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  o.ü₁ = o.y₁par - X₁₂B(parent, parent.X₁, parent.X₂, o.β̂)
	nothing
end

function MakeResidualsIV!(o::StrEstimator{T}, parent::StrBootTest{T}) where T
  o.ü₁ = o.y₁par - o.Z * o.β̂

  if !parent.scorebs
	  _β = [1 ; -o.β̂]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β̂  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  Π̂ = invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β̂)' + o.XY₂)
		o.Ü₂ = X₁₂B(parent, o.X₁, o.X₂, Π̂); o.Ü₂ .= o.Y₂ .- o.Ü₂

	  o.u⃛₁ = o.Ü₂ * (o.t₁Y + o.RparY * o.β̂); o.u⃛₁ .+= o.ü₁
# Ü₂par = [o.y₁par - o.Z * o.β̂; o.Y₂ - X₁₂B(parent, o.X₁, o.X₂, Π̂)] * [1                        o.DGP.Ü₂ * o.Repl.RparY
# 			                                                       (o.t₁Y + o.RparY * o.β̂)  0                      ]
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
	    o.CT_XAR = [crosstabFEt(parent, view(o.WXAR,:,d), parent.info⋂Data) for d ∈ 1:parent.dof]
	  end
  end
	nothing
end
