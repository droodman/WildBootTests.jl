# Logically, wild bootstrap tests perform estimation at two stages, once as part of the bootstrap DGP, once in each bootstrap replication
# The StrEstimator "class" and its three "children" hold the estimation logic for the OLS, Anderson-Rubin, and IV/GMM cases

abstract type Estimator end
struct OLS<:Estimator end
struct ARubin<:Estimator end
struct IVGMM<:Estimator end

mutable struct StrEstimator{T<:AbstractFloat, E<:Estimator}
  parent
  isDGP::Bool; LIML::Bool; Fuller::T; κ::T
  R₁perp::Union{Matrix{T},UniformScaling{Bool}}; Rpar::Union{Matrix{T},UniformScaling}

  kZ::Int64
  y₁::Vector{T}; ü₁::Vector{T}; u⃛₁::Vector{T}; β::Vector{T}; β₀::Vector{T}; PXy₁::Vector{T}; invXXXy₁par::Vector{T}
  Yendog::Vector{Bool}
  invZperpZperp::Matrix{T}; ZperpinvZperpZperp::Matrix{T}; XZ::Matrix{T}; PXZ::Matrix{T}; YPXY::Matrix{T}; R₁invR₁R₁::Union{Matrix{T},UniformScaling}
	RperpX::Union{Matrix{T},UniformScaling,SelectionMatrix}; RperpXperp::Union{Matrix{T},UniformScaling,SelectionMatrix}; RRpar::Matrix{T}; RparY::Union{Matrix{T},UniformScaling,SelectionMatrix}; RR₁invR₁R₁::Matrix{T}
	∂β∂r::Matrix{T}; YY::Matrix{T}; AR::Matrix{T}; XAR::Matrix{T}; R₁invR₁R₁Y::Matrix{T}; invXXXZ::Matrix{T}; Ü₂::Matrix{T}; XinvXX::Matrix{T}; Rt₁::Vector{T}
	invXX::Matrix{T}; Y₂::Matrix{T}; X₂::Matrix{T}; invH
	y₁par::Vector{T}; Xy₁par::Vector{T}
	A::Matrix{T}; Z::Matrix{T}; Zperp::Matrix{T}; X₁::Matrix{T}
	FillingT₀::Matrix{Vector{T}}
	WXAR::Matrix{T}; ScapPXYZperp::Vector{Matrix{T}}; ScapYX::Vector{Matrix{T}}; CT_XAR::Array{T,3}; CT_FEcapPY::Vector{Matrix{T}}

  # IV/GMM only
  ZZ::Matrix{T}; XY₂::Matrix{T}; XX::Matrix{T}; H_2SLS::Matrix{T}; V::Matrix{T}; ZY₂::Matrix{T}; X₂Y₂::Matrix{T}; X₁Y₂::Matrix{T}; ZR₁ZR₁::Matrix{T}; X₂ZR₁::Matrix{T}; ZR₁Y₂::Matrix{T}; X₁ZR₁::Matrix{T}
  ZZR₁::Matrix{T}; X₂y₁::Vector{T}; X₁y₁::Vector{T}; Zy₁::Vector{T}; ZXinvXXXZ::Matrix{T}; H_2SLSmZZ::Matrix{T}
  ZXinvXXXy₁par::Vector{T}; t₁Y::Vector{T}
  Y₂y₁::Vector{T}; twoR₁Zy₁::Vector{T}
  y₁y₁::T; y₁pary₁par::T
  X₂y₁par::Vector{T}; X₁y₁par::Vector{T}; Zy₁par::Vector{T}
  Y₂y₁par::Vector{T}
  Rperp::Matrix{T}; ZR₁::Matrix{T}
  kX::Integer

  StrEstimator{T,E}(parent) where T<:AbstractFloat where E<:Estimator = new(parent, true, E==IVGMM, false, T(E==IVGMM ? 1 : 0), Matrix{T}(undef,0,0), I)
end

function perp(A::AbstractMatrix)
  F = eigensym(A*invsym(A'A)*A')
  F.vectors[:, abs.(F.values) .< 1000eps(eltype(A))]
  # checkI!(dest)
end

# R₁ is constraints. R is attack surface for null; only needed when using FWL for WRE
# for DGP regression, R₁ is maintained constraints + null if imposed while R should have 0 nrows
# for replication regressions R₁ is maintained constraints, R is null
function setR!(o::StrEstimator{T,E}, R₁::AbstractMatrix{T}, R::Union{UniformScaling{Bool},AbstractMatrix{T}}=Matrix{T}(undef,0,0)) where {T,E}
  if length(R₁) > 0
	  invR₁R₁ = invsym(R₁ * R₁')
	  all(iszero.(diag(invR₁R₁))) && throw(ErrorException("Null hypothesis or model constraints are inconsistent or redundant."))
	  o.R₁invR₁R₁ = R₁'invR₁R₁
	  F = eigensym(o.R₁invR₁R₁ * R₁)
	  o.R₁perp = F.vectors[:, abs.(F.values) .< 1000*eps(T)]  # eigenvectors orthogonal to span of R₁; foundation for parameterizing subspace compatible with constraints
  else
	  o.R₁invR₁R₁ = Matrix{T}(undef, o.parent.kZ, 0)  # and R₁perp = I
  end

  if !iszero(o.κ)
	  RR₁perp = Matrix([R ; zeros(T, o.parent.kY₂, o.parent.kX₁) I])  # nrows to prevent partialling out of endogenous regressors; convert sparse matrix produced by constructor to dense

	  length(R₁)>0 && (RR₁perp *= o.R₁perp)

	  F = eigensym(RR₁perp'pinv(RR₁perp * RR₁perp')*RR₁perp); val = abs.(F.values) .> 1000*eps(T)
	  o.Rpar   = F.vectors[:,   val]  # perp and par of RR₁perp
    o.RperpX = F.vectors[:, .!val]

	  if length(R₁) > 0  # fold model constraint factors into Rpar, RperpX
	    o.Rpar   = o.R₁perp * o.Rpar
	    o.RperpX = o.R₁perp * o.RperpX
	  end

	  o.RRpar = R * o.Rpar
	  o.RperpX = o.RperpX[1:o.parent.kX₁,:]  # Zperp=Z*RperpX; though formally a multiplier on Z, it will only extract exogenous components, in X₁, since all endogenous ones will be retained
		o.RperpXperp = selectify(perp(o.RperpX))
		o.RperpX = selectify(o.RperpX)
	  o.RparY = selectify(o.Rpar[o.parent.kX₁+1:end,:])  # part of Rpar that refers to Y₂
	  o.R₁invR₁R₁Y = o.R₁invR₁R₁[o.parent.kX₁+1:end,:]
	  o.RR₁invR₁R₁ = R * o.R₁invR₁R₁
  end
end

# stuff that can be done before r set, and depends only on exogenous variables, which are fixed throughout all bootstrap methods
function InitVars!(o::StrEstimator{T,OLS}, Rperp::AbstractMatrix{T}) where T # Rperp is for replication regression--no null imposed
  o.y₁par = o.parent.y₁
  H = symcross(o.parent.X₁, o.parent.wt)
  o.invH = inv(H)

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.invH : Symmetric(o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp')  # for DGP regression
  o.β₀ = R₁AR₁ * crossvec(o.parent.X₁, o.parent.wt, o.y₁par)
  o.∂β∂r = R₁AR₁ * H * o.R₁invR₁R₁ - o.R₁invR₁R₁

  o.A = iszero(length(Rperp)) ? o.invH : Rperp * invsym(Rperp'H*Rperp) * Rperp'  # for replication regression
  o.AR = o.A * o.parent.R'
  (o.parent.scorebs || o.parent.robust) && (o.XAR = o.parent.X₁ * o.AR)
end

function InitVars!(o::StrEstimator{T,ARubin}, Rperp::AbstractMatrix{T} = Matrix{T}(undef,0,0)) where T
  X₂X₁ = cross(o.parent.X₂, o.parent.wt, o.parent.X₁)
  H = Symmetric([symcross(o.parent.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.parent.X₂, o.parent.wt)])  # XXX use LazyArrays?
  o.A = inv(H)
  o.AR = o.A * o.parent.R'
  (o.parent.scorebs || o.parent.robust) && (o.XAR = X₁₂B(o.parent.X₁, o.parent.X₂, o.AR))

  R₁AR₁ = iszero(length(o.R₁perp)) ? o.A : o.R₁perp * invsym(o.R₁perp'H*o.R₁perp) * o.R₁perp'
  o.β₀   = R₁AR₁ * [crossvec(o.parent.X₁, o.parent.wt, o.parent.y₁) ; crossvec(o.parent.X₂, o.parent.wt, o.parent.y₁)]
  o.∂β∂r = R₁AR₁ * [cross(o.parent.X₁, o.parent.wt, o.parent.Y₂) ; cross(o.parent.X₂, o.parent.wt, o.parent.Y₂)]
end

function InitVars!(o::StrEstimator{T,IVGMM}, Rperp::AbstractMatrix{T}...) where T
  !isempty(Rperp) && (o.Rperp = Rperp[1])

  o.Zperp = o.parent.X₁ * o.RperpX
  o.invZperpZperp = iszero(length(o.Zperp)) ? Matrix{T}(undef,0,0) : inv(symcross(o.Zperp, o.parent.wt))
  o.ZperpinvZperpZperp = o.Zperp * o.invZperpZperp

  o.X₁ = o.parent.X₁ * o.RperpXperp; o.X₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.X₁)  # FWL-process X₁
  o.X₂ = o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.X₂)
		o.X₂ .= o.parent.X₂ .- o.X₂                 # FWL-process X₂
  X₂X₁ = cross(o.X₂, o.parent.wt, o.X₁)
  o.XX = Symmetric([symcross(o.X₁, o.parent.wt) X₂X₁' ; X₂X₁ symcross(o.X₂, o.parent.wt)])
  o.kX = ncols(o.XX)
  o.invXX = invsym(o.XX)
  o.Z   = X₁₂B(o.parent.X₁, o.parent.Y₂, o.Rpar     )  # Zpar
  o.ZR₁ = X₁₂B(o.parent.X₁, o.parent.Y₂, o.R₁invR₁R₁)

  o.Z   .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.Z  )  # partialling out
  o.ZR₁ .-= o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.ZR₁)
  o.Y₂ = o.parent.Y₂ - o.ZperpinvZperpZperp * cross(o.Zperp, o.parent.wt, o.parent.Y₂)
  o.y₁ = o.parent.y₁ - o.ZperpinvZperpZperp * crossvec(o.Zperp, o.parent.wt, o.parent.y₁)

  o.X₁Y₂ = cross(o.X₁, o.parent.wt, o.Y₂)
  o.X₂Y₂ = cross(o.X₂, o.parent.wt, o.Y₂)
  o.XY₂ = [o.X₁Y₂ ; o.X₂Y₂]
  o.Y₂y₁ = crossvec(o.Y₂, o.parent.wt, o.y₁)
  o.X₂y₁ = crossvec(o.X₂, o.parent.wt, o.y₁)
  o.X₁y₁ = crossvec(o.X₁, o.parent.wt, o.y₁)
  o.y₁y₁ = cross(o.y₁ , o.parent.wt, o.y₁)[1]
  o.Zy₁  = crossvec(o.Z, o.parent.wt, o.y₁)
  o.XZ   = [cross(o.X₁, o.parent.wt, o.Z) ;
			      cross(o.X₂, o.parent.wt, o.Z)]
  o.ZY₂ =  cross(o.Z, o.parent.wt, o.Y₂)
  o.ZZ  =  symcross(o.Z, o.parent.wt)

  o.invXXXZ = o.invXX * o.XZ
  o.ZXinvXXXZ = o.XZ'o.invXXXZ

  if length(o.R₁invR₁R₁)>0
	  o.X₂ZR₁    = cross(o.X₂, o.parent.wt, o.ZR₁)
	  o.X₁ZR₁    = cross(o.X₁, o.parent.wt, o.ZR₁)
	  o.ZZR₁     = cross(o.Z , o.parent.wt, o.ZR₁)
	  o.twoR₁Zy₁ = 2crossvec(o.ZR₁, o.parent.wt, o.y₁)
	  o.ZR₁ZR₁   = symcross(o.ZR₁, o.parent.wt)
	  o.ZR₁Y₂    = cross(o.ZR₁, o.parent.wt, o.Y₂)
  else
	  o.Y₂y₁par    = o.Y₂y₁
	  o.X₂y₁par    = o.X₂y₁
	  o.X₁y₁par    = o.X₁y₁
	  o.Zy₁par     = o.Zy₁
	  o.y₁pary₁par = o.y₁y₁
	  o.Xy₁par     = [o.X₁y₁ ; o.X₂y₁]
	  o.y₁par      = o.y₁
  end

  o.V =  o.invXX * o.XZ # in 2SLS case, StrEstimator is (V' XZ)^-1 * (V'Xy₁). Also used in k-class and LIML robust VCV by Stata convention
  o.H_2SLS = Symmetric(o.V'o.XZ)  # Hessian
  (o.LIML || o.κ ≠ 1) && (o.H_2SLSmZZ = o.H_2SLS - o.ZZ)

  if o.isDGP
	  !o.LIML && MakeH!(o, !isempty(Rperp)) # DGP is LIML except possibly when getting confidence peak for A-R plot; but LIML=0 when exactly id'd, for then κ=1 always and Hessian doesn't depend on r₁ and can be computed now
  else
	  o.kZ = ncols(o.Rpar)
	  o.Yendog = [true; o.RparY==I ? fill(true, o.kZ) : vec(colsum(o.RparY.≠0)).>0]  # columns of Y = [y₁par Zpar] that are endogenous (normally all)

	  if o.parent.robust && o.parent.bootstrapt  # for WRE replication regression, prepare for CRVE
			o.ScapYX       = Vector{Matrix{T}}(undef, o.kZ+1)
			o.ScapPXYZperp = Vector{Matrix{T}}(undef, o.kZ+1)
			o.XinvXX = X₁₂B(o.X₁, o.X₂, o.invXX)
			
			o.PXZ = X₁₂B(o.X₁, o.X₂, o.invXXXZ)

			o.FillingT₀ = Matrix{Matrix{T}}(undef, o.kZ+1, o.kZ+1)  # fixed component of groupwise term in sandwich filling
			o.parent.NFE>0 &&
				(o.CT_FEcapPY = Vector{Matrix{T}}(undef, o.kZ+1))
			for i ∈ 1:o.kZ
				uwt = vHadw(view(o.PXZ,:,i), o.parent.wt)
				o.parent.NFE>0 &&
				(o.CT_FEcapPY[i+1] = crosstabFEt(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)
				tmp = @panelsum(o.Z, uwt, o.parent.infoCapData)
				for j ∈ 1:o.kZ
					o.FillingT₀[i+1,j+1] = view(tmp,:,j)
				end
			end

			for i ∈ 1:o.kZ  # precompute various clusterwise sums
				o.ScapPXYZperp[i+1] = @panelsum(o.Zperp, vHadw(view(o.PXZ,:,i), o.parent.wt), o.parent.infoCapData)  # Scap(P_(MZperpX) * Z .* Zperp)
				!o.parent.granular &&
					(o.ScapYX[i+1] = @panelsum2(o.X₁, o.X₂, vHadw(view(o.Z,:,i), o.parent.wt), o.parent.infoCapData))  # Scap(M_Zperp[Z or y₁] .* P_(MZperpX)])
	    end
	  end
  end
end


# do most of estimation; for LIML r₁ must be passed now in order to solve eigenvalue problem involving it
# inconsistency: for replication regression of Anderson-Rubin, r₁ refers to the *null*, not the maintained constraints, because that's what affects the endogenous variables
# For OLS, compute β₀ (β when r=0) and ∂β∂r without knowing r₁, for efficiency
# For WRE, should only be called once for the replication regressions, since for them r₁ is the unchanging model constraints
function Estimate!(o::StrEstimator{T,OLS} where T, r₁::AbstractVector)
  o.β = o.β₀ - o.∂β∂r * r₁
end

function Estimate!(o::StrEstimator{T,ARubin} where T, r₁::AbstractVector)
  o.β = o.β₀ - o.∂β∂r * r₁
  o.y₁par = o.parent.y₁ - o.parent.Y₂ * (isone(length(r₁)) ? r₁[1] : r₁)
end

function MakeH!(o::StrEstimator{T,IVGMM} where T, makeXAR::Bool=false)
  H = isone(o.κ) ? o.H_2SLS : o.ZZ + o.κ * o.H_2SLSmZZ
  o.invH = invsym(H)
  if makeXAR  # for replication regression in score bootstrap of IV/GMM
	  o.A = length(o.Rperp)>0 ? o.Rperp * invsym(o.Rperp'H*o.Rperp) * o.Rperp' : o.invH
	  o.AR = o.A * (o.Rpar'o.parent.R')
	  o.XAR = X₁₂B(o.X₁, o.X₂, o.V * o.AR)
  end
end

function Estimate!(o::StrEstimator{T,IVGMM} where T, r₁::AbstractVector)
  if length(o.R₁invR₁R₁)>0
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
  o.YY = [o.y₁pary₁par o.Zy₁par' ; o.Zy₁par o.ZZ]
  o.YPXY = [o.invXXXy₁par'o.Xy₁par  o.ZXinvXXXy₁par' ; o.ZXinvXXXy₁par o.ZXinvXXXZ]

  if o.isDGP
	  if o.LIML
	    o.κ = 1/(1 - eigvals(invsym(o.YY) * o.YPXY)[1])  # like Fast & Wild (81), but more stable, at least in Mata
	    !iszero(o.Fuller) && (o.κ -= o.Fuller / (o.parent._Nobs - o.parent.kX))
	    MakeH!(o)
	  end

	  o.β = o.invH * (isone(o.κ) ? o.ZXinvXXXy₁par : o.κ * (o.ZXinvXXXy₁par - o.Zy₁par) + o.Zy₁par)
	  o.t₁Y = o.R₁invR₁R₁Y * r₁
  elseif o.parent.WREnonARubin  # if not score bootstrap of IV/GMM
	  o.Rt₁ = o.RR₁invR₁R₁ * r₁

	  if o.parent.robust && o.parent.bootstrapt  # prepare WRE replication regressions
	    o.PXy₁ = X₁₂B(o.X₁, o.X₂, o.invXXXy₁par)

	    uwt = vHadw(o.PXy₁, o.parent.wt)
	    tmp = @panelsum(o.Z, uwt, o.parent.infoCapData)
	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[1,i+1] = view(tmp,:,i)
	    end
	    o.ScapPXYZperp[1] = @panelsum(o.Zperp, uwt, o.parent.infoCapData)  # Scap(P_(MZperpX) * y₁ .* Zperp)

	    o.parent.NFE>0 &&
	  		(o.CT_FEcapPY[1] = crosstabFEt(o.parent, uwt, o.parent.infoCapData) .* o.parent.invFEwt)

	    uwt = vHadw(o.y₁par, o.parent.wt)
	    tmp = @panelsum(o.PXZ, uwt, o.parent.infoCapData)
	    for i ∈ 1:o.kZ
	  	  o.FillingT₀[i+1,1] = view(tmp,:,i)
	    end
	    o.FillingT₀[1] = @panelsum(o.PXy₁, uwt, o.parent.infoCapData)
	    !o.parent.granular &&
	  	(o.ScapYX[1] = @panelsum2(o.X₁, o.X₂, uwt, o.parent.infoCapData))  # Scap(M_Zperp*y₁ .* P_(MZperpX)])
	  end
  end
end


@inline function MakeResiduals!(o::StrEstimator{T,<:Union{OLS,ARubin}} where T)
  o.ü₁ = o.y₁par - X₁₂B(o.parent.X₁, o.parent.X₂, o.β)
end

function MakeResiduals!(o::StrEstimator{T,IVGMM} where T)
  o.ü₁ = o.y₁par - o.Z * o.β

  if !o.parent.scorebs
	  _β = [1 ; -o.β]
	  uu = _β'o.YY * _β

	  Xu = o.Xy₁par - o.XZ * o.β  # after DGP regression, compute Y₂ residuals by regressing Y₂ on X while controlling for y₁ residuals, done through FWL
	  negXuinvuu = Xu / -uu
	  o.Ü₂ = o.Y₂ - X₁₂B(o.X₁, o.X₂, invsym(o.XX + negXuinvuu * Xu') * (negXuinvuu * (o.Y₂y₁par - o.ZY₂'o.β)' + o.XY₂))  # large expression is Π̂

	  o.u⃛₁ = o.ü₁ + o.Ü₂ * (o.t₁Y + o.RparY * o.β)
  end
end


# non-WRE stuff that only depends on r in A-R case, for test stat denominators in replication regressions
# since the non-AR OLS code never creates an object for replication regresssions, in that case this is called on the DGP regression object
# depends on results of Estimate() only when doing OLS-style bootstrap on an overidentified IV/GMM regression--score bootstrap or A-R. Then κ from DGP LIML affects Hessian, H.
function InitTestDenoms!(o::StrEstimator)
  if o.parent.bootstrapt && (o.parent.scorebs || o.parent.robust)
	  (o.parent.granular || o.parent.purerobust) && (o.WXAR = vHadw(o.XAR, o.parent.wt))

	  if o.parent.robust && o.parent.NFE>0 && !(o.parent.FEboot || o.parent.scorebs) && o.parent.granular < o.parent.NErrClustCombs  # make first factor of second term of (64) for c=cap (c=1)
	    !isdefined(o, :WXAR) && (o.WXAR = vHadw(o.XAR, o.parent.wt))
	    o.CT_XAR = crosstabFEt(o.parent, o.WXAR, o.parent.infoCapData)
	  end
  end
end
