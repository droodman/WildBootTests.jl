@inline mysortperm(X::AbstractVecOrMat) = size(X,2)==1 ? sortperm(ndims(X)==1 ? X : vec(X), alg=isone(Threads.nthreads()) ? RadixSort : ThreadsX.MergeSort) : 
                                                         sortperm(collect(eachrow(X))     , alg=isone(Threads.nthreads()) ? TimSort   : ThreadsX.MergeSort)

function Init!(o::StrBootTest{T}) where T  # for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
	o.kX = o.kX₁ + o.kX₂
	o.kX₂==0 && (o.X₂ = zeros(T,o.Nobs,0))
  iszero(o.kY₂) && (o.Y₂ = zeros(T,o.Nobs,0))
  o.kZ = o.kX₁ + o.kY₂
  if o.liml && o.kX₂==o.kY₂  # exactly identified liml = 2SLS
  	o.κ = one(T)
  	o.liml = false
  end
  if iszero(length(o.R₁))  # base model contains no restrictions?
    o.R₁ = zeros(T,0,o.kZ)
    o.r₁ = zeros(T,0)
  end
  isnan(o.κ) && (o.κ = o.kX₂>0 ? one(T) : zero(T))  # if κ in κ-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  o.WRE = !(iszero(o.κ) || o.scorebs) || o.arubin

  iszero(o.B) && (o.scorebs = true)

  o.haswt = !iszero(length(o.wt))
	if o.haswt
		o.sumwt = sum(o.wt)
		o._Nobs = o.fweights ? o.sumwt : T(o.Nobs)
	else
		o.sumwt = one(T)
		o._Nobs = T(o.Nobs)
		o.sqrtwt = T[]
	end

	o.subcluster = o.NClustVar - o.NErrClustVar

	if !(iszero(o.NClustVar) || o.issorted)
		if iszero(o.subcluster)  # sort by err cluster vars, then remaining boot cluster vars
			p = mysortperm(o.clustid)
		else
			p = mysortperm(view(o.clustid, :, [collect(o.subcluster+1:o.NClustVar); collect(1:o.subcluster)]))
		end

		o.clustid = ndims(o.clustid)==1 ? o.clustid[p] : o.clustid[p,:]
		o.X₁ = ndims(o.X₁)==1 ? o.X₁[p] : o.X₁[p,:]
		o.X₂ = ndims(o.X₂)==1 ? o.X₂[p] : o.X₂[p,:]
		o.y₁ = o.y₁[p]
		o.Y₂ = ndims(o.Y₂)==1 ? o.Y₂[p] : o.Y₂[p,:]
		isdefined(o, :FEID) && nrows(o.FEID)>0 && (o.FEID = o.FEID[p])
		o.haswt && (o.wt = o.wt[p])
		o.overwrite = true  # data matrices are no longer those passed by caller
	end
	if o.haswt
		o.sqrtwt = sqrt.(o.wt)
		if o.overwrite
			o.y₁ .*= o.sqrtwt
			length(o.Y₂)>0 && (o.Y₂ .*= o.sqrtwt)
			length(o.X₁)>0 && (o.X₁ .*= o.sqrtwt)
			length(o.X₂)>0 && (o.X₂ .*= o.sqrtwt)
		else
			o.sqrtwt = sqrt.(o.wt)
			o.y₁ = o.y₁ .* o.sqrtwt
			length(o.Y₂)>0 && (o.Y₂ = o.Y₂ .* o.sqrtwt)
			length(o.X₁)>0 && (o.X₁ = o.X₁ .* o.sqrtwt)
			length(o.X₂)>0 && (o.X₂ = o.X₂ .* o.sqrtwt)
			o.overwrite = true
		end
	end

  if o.WREnonARubin
	  if iszero(o.NClustVar)
	    o.info⋂ = o.info✻ = UnitRange{Int64}[]  # [i:i for i in 1:o.Nobs]  # no clustering, so no collapsing by cluster
			o.ID⋂ = o.ID✻ = collect(1:o.Nobs)
	  else
	    o.info✻, o.ID✻ = panelsetupID(o.clustid, 1:o.NBootClustVar)
  	end
		o.willfill = o.robust && o.bootstrapt  # will compute sandwich filling for robust denominator?
		o.not2SLS = o.liml || !o.robust || !isone(o.κ)  # sometimes κ ≠ 1?
  elseif iszero(o.NClustVar)
	  o.info✻ = UnitRange{Int64}[]  # [i:i for i ∈ 1:o.Nobs]  # causes no collapsing of data in panelsum() calls
		o.ID✻ = collect(1:o.Nobs)
  else
	  o.info✻, o.ID✻ = panelsetupID(o.clustid, 1:min(o.NClustVar,o.NBootClustVar))  # bootstrap cluster grouping defs rel to original data
  end
	o.N✻ = iszero(nrows(o.info✻)) ? o.Nobs : nrows(o.info✻)

	if o.NClustVar > o.NBootClustVar  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
		o.info✻⋂, ID✻⋂ = panelsetupID(o.clustid, 1:o.NClustVar)
		if o.subcluster>0 && nrows(o.info✻) ≠ nrows(o.info✻⋂)
			throw(ErrorException("\nThis program can only perform the subcluster bootstrap when the bootstrap clusters are nested within the (intersections of the) error clusters.\n"))
		end
	else
		o.info✻⋂, ID✻⋂ = o.info✻, o.ID✻  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
	end
	o.N✻⋂ = iszero(nrows(o.info✻⋂)) ? o.Nobs : nrows(o.info✻⋂)

	if o.NClustVar > o.NErrClustVar  # info for intersections of error clustering wrt data
		o.info⋂, o.ID⋂ = panelsetupID(o.clustid, o.subcluster+1:o.NClustVar)
		clustid⋂_✻⋂ = length(o.info⋂)==o.Nobs ? o.clustid : o.clustid[first.(o.info⋂ ),:]  # version of clustid matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
		clustid✻⋂ = o.N✻⋂==o.Nobs ? o.clustid : o.clustid[first.(o.info✻⋂),:]  # version of clustid matrix with one row for each all-bootstrap && error cluster-var intersection instead of 1 row for each obs
	else
		o.info⋂, o.ID⋂ = o.info✻⋂, ID✻⋂  # info for intersections of error clustering wrt data
		clustid✻⋂ = clustid⋂_✻⋂ = nrows(o.info⋂)==o.Nobs || iszero(nrows(o.info⋂)) ? o.clustid : o.clustid[first.(o.info⋂),:]  # version of clustid matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
	end
	o.WREnonARubin && 
		(o.info✻_✻⋂ = iszero(o.NBootClustVar) ? o.info✻⋂ : panelsetup(clustid✻⋂, 1:o.NBootClustVar))
	o.N⋂ = iszero(nrows(o.info⋂)) ? o.Nobs : nrows(o.info⋂)

	if o.bootstrapt
    if o.NClustVar>0
	    minN = T(Inf)

	    combs = [x & 2^y > 0 for x in 2^o.NErrClustVar-1:-1:1, y in o.NErrClustVar-1:-1:0]  # represent all error clustering combinations. First is intersection of all error clustering vars
	    o.clust = Vector{StrClust{T}}(undef, nrows(combs))  # leave out no-cluster combination
	    o.NErrClustCombs = length(o.clust)

	    o.BootClust = 2^(o.NClustVar - o.NBootClustVar)  # location of bootstrap clustering within list of cluster combinations

	    @inbounds for c ∈ 1:o.NErrClustCombs  # for each error clustering combination
		    ClustCols = o.subcluster .+ findall(@view combs[c,:])
		    even = isodd(length(ClustCols))  # not a typo

		    if isone(c)
		      if iszero(o.subcluster)
		    	  order = Vector{Int64}(undef,0)
		    	  info  = Vector{UnitRange{Int64}}(undef, 0)  # causes no collapsing of data in panelsum() calls
						N = o.N✻⋂
		      else
		    	  order = mysortperm(@view clustid⋂_✻⋂[:,ClustCols])
		    	  info  = panelsetup(view(clustid⋂_✻⋂,order,:), ClustCols)
						N = nrows(info)
		      end
		    else
		      if any(combs[c, minimum(findall(combs[c,:] .≠ combs[c-1,:])):end])  # if this sort ordering same as last to some point and missing thereafter, no need to re-sort
		    	  order = mysortperm(@view clustid⋂_✻⋂[:,ClustCols])
		    	  info = panelsetup(view(clustid⋂_✻⋂,order,:), ClustCols)
		      else
		    	  order = Vector{Int64}(undef,0)
						info = panelsetup(clustid⋂_✻⋂, ClustCols)
		      end
					N = nrows(info)
		    end

		    minN = min(minN,N)
				multiplier = (even*2-1) * (o.clusteradj && !o.clustermin ? T(N / (N-1)) : one(T))
				o.clust[c] = StrClust{T}(N, multiplier, even, order, info)
	    end

	    (o.scorebs || !o.WREnonARubin) &&
		  	(o.ClustShare = o.haswt ? @panelsum(o.wt, o.info⋂)/o.sumwt : T.(length.(o.info⋂)./o.Nobs)) # share of observations by group
    else  # if no clustering, cast "robust" as clustering by observation
      o.clust = [StrClust{T}(o.Nobs, o.small ? o._Nobs / (o._Nobs - one(T)) : one(T), true, Vector{Int64}(undef,0), Vector{UnitRange{Int64}}(undef,0))]
      o.NErrClustCombs = one(Int16)
      (o.scorebs || !o.WREnonARubin) &&
    		(o.ClustShare = o.haswt ? o.wt/o.sumwt : [one(T)/o.Nobs])
    end
  else
	  minN = T(nrows(o.info✻))
	end

	o.purerobust = o.robust && !o.scorebs && iszero(o.subcluster) && o.N✻==o.Nobs  # do we ever error-cluster *and* bootstrap-cluster by individual?
	o.granular   = o.WREnonARubin ? 4. *o.Nobs/1*o.B < o.N✻/1*(o.Nobs+o.N⋂/1*o.B) :  #  "/1" to convert to Float and avoid overflow
																	o.robust && !o.scorebs && (o.purerobust || (o.N⋂/1+o.N✻)*o.kZ/1*o.B + (o.N⋂/1-o.N✻)*o.B + o.kZ/1*o.B < o.N⋂*o.kZ^2. + o.Nobs/1*o.kZ + o.N⋂/1 * o.N✻/1 * o.kZ + o.N⋂/1 * o.N✻)
	if o.jk
		if !o.WREnonARubin
			o.granularjk = o.kZ^3. + o.N✻ * (o.Nobs/o.N✻*o.kZ^2. + (o.Nobs/o.N✻)^2*o.kZ + (o.Nobs/o.N✻)^2 + (o.Nobs/o.N✻)^3) < o.N✻ * (o.kZ^2. *o.Nobs/o.N✻ + o.kZ^3. + 2o.kZ*(o.kZ + o.Nobs/o.N✻))
		end
		(o.WREnonARubin || o.granularjk) && !o.purerobust && (o.maxNg = mapreduce(length, max, o.info✻))
	end
	
	(o.WREnonARubin && o.willfill && !o.granular || (o.subcluster>0 || o.granular || o.purerobust) && !o.WREnonARubin) && 
		(o.info⋂_✻⋂ = panelsetup(clustid✻⋂, o.subcluster+1:o.NClustVar))  # info for error clusters wrt data collapsed to intersections of all bootstrapping && error clusters; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusterings

	InitFEs!(o)
	if o.B>0 && o.robust && o.granular && o.bootstrapt && !o.WREnonARubin
		if o.purerobust
			if !isdefined(o, :ID✻)
				o.ID✻ = Vector{Int64}[]  # panelsum treats 0 length same as max length
			end
		elseif o.NFE>0 && !o.FEboot
			if !isdefined(o, :ID✻)
				_, o.ID✻ = panelsetupID(o.clustid, 1:o.NBootClustVar)
			end
		else
			if iszero(length(o.ID✻_✻⋂))
				_, o.ID✻_✻⋂ = panelsetupID(clustid✻⋂, 1:o.NBootClustVar)
			end
		end
	end

	o.enumerate = o.B>0 && o.auxtwtype == :rademacher && o.N✻*log(2) < log(o.B)+1e-6  # generate full Rademacher set?
	if o.enumerate
		o.maxmatsize = 0
		o.B = 2^o.N✻
		o.Nw = 1
	else
		o.Nw = iszero(o.maxmatsize) ? one(Int64) : ceil(Int64, (o.B+1.) * max(nrows(o.ID✻), length(o.ID✻_✻⋂), o.N✻) * sizeof(T) / o.maxmatsize / 1073741824) # 1073741824 = giga(byte)
	end

	if isone(o.Nw)
		o.v = Matrix{T}(undef, o.N✻, o.B+1)
		MakeWildWeights!(o, o.B, first=true)  # make all wild weights, once
		o.ncolsv = o.B + 1
		o.WeightGrp = [1:o.ncolsv]
	else
    o.seed = rand(o.rng, UInt64)
		o.ncolsv = ceil(Int64, (o.B+1) / o.Nw)
		o.v = Matrix{T}(undef, o.N✻, o.ncolsv)
		o.Nw = ceil(Int64, (o.B+1) / o.ncolsv)
		o.WeightGrp = [(i-1)*o.ncolsv+1:i*o.ncolsv for i ∈ 1:o.Nw]
		o.B = o.Nw * o.ncolsv - 1  # of replications may be slightly increased so each block of v same size
	end

	if o.ml
		o.dof = nrows(o.R)
	else
		if o.arubin
			o.R = hcat(zeros(o.kX₂,o.kX₁), Matrix(I(o.kX₂)))  # attack surface is all endog vars
			o.R₁ = o.kX₁>0 && nrows(o.R₁)>0 ? hcat(o.R₁[:,1:o.kX₁], zeros(T,nrows(o.R₁),o.kX₂)) : zeros(T,0, o.kX)  # and convert model constraints from referring to X₁, Y₂ to X₁, X₂
		end
		o.dof = nrows(o.R)

		if !o.WRE && iszero(o.κ)  # regular OLS
			o.DGP = StrEstimator{T}(true, o.liml, o.fuller, o.κ)
			o.Repl = StrEstimator{T}(true, false, zero(T), zero(T))  # XXX isDGP=1 for Repl? doesn't matter?
			setR!(o.DGP, o, o.null ? [o.R₁ ; o.R] : o.R₁)  # DGP constraints: model constraints + null if imposed
			setR!(o.Repl, o, o.R₁)  # model constraints only
			InitVarsOLS!(o.DGP, o, o.Repl.R₁perp)
			InitTestDenoms!(o.DGP, o)
			o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
		elseif o.arubin
			if o.willplot  # for plotting/CI purposes get original point estimate since not normally generated
				o.DGP = StrEstimator{T}(true, o.liml, o.fuller, o.κ)
				setR!(o.DGP, o, o.R₁, zeros(T,0,o.kZ))  # no-null model
				InitVarsIV!(o.DGP, o)
				EstimateIV!(o.DGP, o, false, o.r₁)
				o.confpeak = view(o.DGP.β̈  ,:,1)  # estimated coordinate of confidence peak
			end

			o.DGP = StrEstimator{T}(true, false, zero(T), zero(T))
			setR!(o.DGP, o, o.R₁)
			InitVarsARubin!(o.DGP, o)
			InitTestDenoms!(o.DGP, o)
			o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
			o.kZ = o.kX

		elseif o.WREnonARubin
			o.DGP = StrEstimator{T}(true, true, zero(T), one(T))
			if o.null
				setR!(o.DGP, o, [o.R₁ ; o.R], zeros(T,0,o.kZ))  #  DGP constraints: model constraints + imposed null
			else
				setR!(o.DGP, o, o.R₁, o.R)  # when null not imposed, keep it in the attack surface, though not used there, so Zperp is same in DGP and Repl
			end
			InitVarsIV!(o.DGP, o)

			o.Repl = StrEstimator{T}(false, o.liml, o.fuller, o.κ)
			setR!(o.Repl, o, o.R₁, o.R)
			InitVarsIV!(o.Repl, o)
			EstimateIV!(o.Repl, o, false, o.r₁)

			InitWRE!(o)

		else  # the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R₁, one with FWL, one not
			o.DGP = StrEstimator{T}(true, o.liml, o.fuller, o.κ)
			setR!(o.DGP, o, o.null ? [o.R₁ ; o.R] : o.R₁, zeros(T,0,o.kZ))  # DGP constraints: model constraints + null if imposed
			InitVarsIV!(o.DGP, o)

			o.Repl = StrEstimator{T}(true, o.liml, o.fuller, o.κ)
			setR!(o.Repl, o, o.R₁, I)  # process replication restraints = model constraints only; I arg short-circuits FWL processing
			InitVarsIV!(o.Repl, o, o.Repl.R₁perp)
			EstimateIV!(o.Repl, o, false, o.r₁)  # bit inefficient to estimate in both objects, but maintains the conformity
			InitTestDenoms!(o.Repl, o)
			o.M = o.Repl  # StrEstimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scorebs for IV/GMM mixes the two

			if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
				EstimateIV!(o.DGP, o, o.jk, o.r₁)
				MakeResidualsIV!(o.DGP, o)
			end
		end
  end

	o.bootstrapt &&
		(o.denom = [Matrix{T}(undef,1,o.ncolsv) for _ in 1:o.dof, _ in 1:o.dof])

  if !o.WREnonARubin && o.bootstrapt
		if o.robust
			o.Kcd =                        Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof)
			o.Jcd = iszero(o.B) ? o.Kcd : [Matrix{T}(undef, o.clust[c].N, o.ncolsv) for c ∈ 1:o.NErrClustCombs, _ ∈ 1:o.dof]  # if B = 0, Kcd will be multiplied by v, which is all 1's, and will constitute Jcd
		end

		if o.robust && o.granular<o.NErrClustCombs && o.B>0
			inds = o.subcluster>0 ?
				        [CartesianIndex(j, i) for (j,v) ∈ enumerate(o.info⋂_✻⋂) for i ∈ v] :  # crosstab ∩,* is wide
								o.NClustVar == o.NBootClustVar ?
										[CartesianIndex(i, i) for i ∈ 1:o.N✻⋂] :  # crosstab ∩,* is square
										[CartesianIndex(i, j) for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab ∩,* is tall
			o.crosstab⋂✻ind = LinearIndices((1:o.N⋂, 1:o.N✻))[inds]
		end
	end

  o.sqrt = isone(o.dof)  # work with t/z stats instead of F/chi2?

  if o.small
		_dof = max(0, o._Nobs - o.kZ + (o.ml ? 0 : nrows(o.R₁)) - o.FEdfadj)
		o.dof_r = o.NClustVar>0 ? minN - one(T) : _dof  # floating-point dof_r allows for fractional fweights, FWIW...
		o.smallsample = _dof / (o._Nobs - o.robust)
	else
		o.smallsample = one(T)
		o.dof_r = zero(T)
	end
	o.clustermin && (o.smallsample *= (minN - 1) / minN)  # ivreg2-style adjustment when multiway clustering
	o.multiplier = o.small ? o.smallsample / o.dof : o.smallsample  # divide by # of constraints because F stat is so defined
  !(o.robust || o.ml) && (o.multiplier *= o._Nobs)  # will turn sum of squared errors in denom of t/z into mean
  o.sqrt && (o.multiplier = √o.multiplier)

	(o.Nw>1 || !o.WREnonARubin && !(o.robust && o.dof≤2)) && (o.dist = Matrix{T}(undef, 1, o.B+1))
	o.Nw>1 && (o.numer = Matrix{T}(undef, o.dof, o.B+1))

  if !o.WREnonARubin
		o.numerw = Matrix{T}(undef, o.dof, o.ncolsv)
		!o.scorebs && (!o.robust || o.granular || o.purerobust) && (o.β̈dev = Matrix{T}(undef, o.kX, o.ncolsv))

		o.poles = o.anchor = zeros(T,0)
		o.interpolable = o.getci && o.bootstrapt && o.null && o.Nw==1 && (iszero(o.κ) || o.arubin)
		o.interpolate_u = !(o.robust || o.ml)
		(o.interpolate_u || (o.B>0 && o.robust && o.bootstrapt && (o.granular || o.purerobust) && (o.purerobust && !o.interpolable || o.NFE>0 && !o.FEboot))) &&
			(o.u✻ = Matrix{T}(undef,o.Nobs,o.ncolsv))
		if o.interpolable
			o.∂numer∂r = [Matrix{T}(undef, o.dof, o.ncolsv) for _ ∈ 1:o.q]
			o.interpolate_u && (o.∂u∂r = Vector{Matrix{T}}(undef, o.q))
			if o.robust
				o.∂denom∂r   = [Matrix{T}(undef,            1, o.ncolsv) for _ ∈ 1:o.q, _ ∈ 1:o.dof, _ ∈ 1:o.dof]
				o.∂²denom∂r² = [Matrix{T}(undef,            1, o.ncolsv) for _ ∈ 1:o.q, _ ∈ 1:o.q, _ ∈ 1:o.dof, _ ∈ 1:o.dof]
				o.∂Jcd∂r     = [Matrix{T}(undef, o.clust[c].N, o.ncolsv) for _ ∈ 1:o.q, c ∈ 1:o.NErrClustCombs, _ ∈ 1:o.dof]
			end
		end
  end
	o.initialized = true
	nothing
end

function InitFEs!(o::StrBootTest{T}) where T
	if isdefined(o, :FEID) && length(o.FEID)>0
		s = BitSet(o.FEID)
		o.NFE = length(s)
		o._FEID = getindex.(Ref(Dict(zip(s, 1:o.NFE))), o.FEID)  # standardize FE clustid to 1, 2, ...

		sumFEwt = zeros(T, o.NFE)
		if o.haswt
			@inbounds for i ∈ 1:o.Nobs
				sumFEwt[o._FEID[i]] += o.wt[i]
			end
			o.FEwt = o.sqrtwt ./ sumFEwt[o._FEID]
		else
			@inbounds for i ∈ 1:o.Nobs
				sumFEwt[o._FEID[i]] += one(T)
			end
		end
		o.invsumFEwt = one(T) ./ sumFEwt

		# is every FE group inside same bootstrapping?
		o.FEboot = o.B>0 && o.NClustVar>0
		if o.FEboot
			first = fill(true, o.NFE)
			clustrows = Matrix{Int64}(undef, o.NFE, o.NBootClustVar)
			_ID = view(o.clustid, :, 1:o.NBootClustVar)
			js = eachindex(axes(_ID))
			@inbounds for i ∈ eachindex(axes(o.clustid))
				if first[i]
					clustrows[i,:] = _ID[i,:]
					first[i] = false
				else
					for j ∈ js
						if clustrows[i,j] ≠ _ID[i,j]
							o.FEboot = false
							@goto afer_loop
						end
					end
				end
			end
		end
		@label afer_loop

		o.FEdfadj==-1 && (o.FEdfadj = o.NFE)

		o.robust && o.B>0 && o.bootstrapt && !o.FEboot && o.granular < o.NErrClustVar &&
			(o.infoBootAll = panelsetup(clustid✻⋂, 1:o.NBootClustVar))  # info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping && error clusters

		!o.FEboot && o.haswt && (o.Mjw = Vector{T}(undef,o.Nobs))

		if o.overwrite
			partialFE!(o, o.X₁)
			partialFE!(o, o.X₂)
			partialFE!(o, o.y₁)
			partialFE!(o, o.Y₂)
		else
			o.X₁ = partialFE(o, o.X₁)
			o.X₂ = partialFE(o, o.X₂)
			o.y₁ = partialFE(o, o.y₁)
			o.Y₂ = partialFE(o, o.Y₂)
			o.overwrite = o.NFE>0  # safe to further modify data matrices
		end
	else
		o.FEdfadj = 0
	end
end

# draw wild weight matrix of width _B. If first=true, insert column of 1s at front.
const ϕ = (1 + √5)/2

function MakeWildWeights!(o::StrBootTest{T}, _B::Integer; first::Bool=true) where T
	m = o.WREnonARubin && o.jk && o.small ? T(sqrt(1 - 1 / o.N✻)) : one(T)  # finite-sample adjuster for jk regressions
	
	if _B>0  # in scoretest or waldtest WRE, still make v a col of 1's
    if o.enumerate
			o.v[:,2:end] = hcat(digits.(0:2^o.N✻-1, base=2, pad=o.N✻)...)
			lmul!(2*m, o.v); o.v .-= m
		elseif o.auxtwtype == :normal
			randn!(o.rng, o.v); !isone(m) && (o.v .*= m)
		elseif o.auxtwtype == :gamma 
			rand!(o.rng, o.v); o.v .= quantile.(Gamma{T}(4,.5), o.v); o.v .-= T(2); !isone(m) && (o.v .*= m)
		elseif o.auxtwtype == :webb
			rand!(o.rng, o.v, m*T[-√1.5, -1, -√.5, √.5, 1, √1.5])
		elseif o.auxtwtype == :mammen
			rand!(o.rng, o.v); o.v .= getindex.(Ref(m * T[1-ϕ; ϕ]), ceil.(Int16, o.v ./ (ϕ/√5)))
		else
			rand!(o.rng, o.v, T[m, -m])  # Rademacher 
		end
		first && fill!(view(o.v,:,1), one(T))  # keep original residuals in first entry to compute base model stat
  end
	nothing
end
