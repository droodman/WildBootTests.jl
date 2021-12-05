@inline _sortperm(X::AbstractVecOrMat) = size(X,2)==1 ? sortperm(ndims(X)==1 ? X : reshape(X,length(X)), alg=RadixSort) : sortperm(collect(eachrow(X)))  # sort a data matrix

function Init!(o::StrBootTest{T}) where T  # for efficiency when varying r repeatedly to make CI, do stuff once that doesn't depend on r
  o.Nobs = nrows(o.X₁)
  o.NClustVar = ncols(o.ID)
  o.kX = (o.kX₁ = ncols(o.X₁)) + (o.kX₂ = ncols(o.X₂))
  o.kX₂==0 && (o.X₂ = zeros(T,o.Nobs,0))
  o.kY₂ = ncols(o.Y₂)
  iszero(o.kY₂) && (o.Y₂ = zeros(T,o.Nobs,0))
  o.kZ = o.kX₁ + o.kY₂
  if o.LIML && o.kX₂==o.kY₂  # exactly identified LIML = 2SLS
  	o.κ = one(T)
  	o.LIML = false
  end
  if !(o.REst = length(o.R₁)>0)  # base model contains no restrictions?
    o.R₁ = zeros(T,0,o.kZ)
    o.r₁ = zeros(T,0)
  end
  isnan(o.κ) && (o.κ = o.kX₂>0 ? one(T) : zero(T))  # if κ in κ-class estimation not specified, it's 0 or 1 for OLS or 2SLS
  o.WRE = !(iszero(o.κ) || o.scorebs) || o.ARubin
  o.WREnonARubin = o.WRE && !o.ARubin

  iszero(o.B) && (o.scorebs = true)

  o.haswt = !iszero(nrows(o.wt))
  o.sumwt = o.haswt ? sum(o.wt) : 1
  o._Nobs = o.haswt && o.fweights ? o.sumwt : o.Nobs

  if !(iszero(o.NClustVar)  || o.issorted)
		o.subcluster = o.NClustVar - o.nerrclustvar
		p = _sortperm(view(o.ID, :, [collect(o.subcluster+1:o.NClustVar); collect(1:o.subcluster)]))  # sort by err cluster vars, then remaining boot cluster vars
		o.ID = ndims(o.ID)==1 ? o.ID[p] : o.ID[p,:]
		o.X₁ = ndims(o.X₁)==1 ? o.X₁[p] : o.X₁[p,:]
		o.X₂ = ndims(o.X₂)==1 ? o.X₂[p] : o.X₂[p,:]
		o.y₁ = o.y₁[p]
		o.Y₂ = ndims(o.Y₂)==1 ? o.Y₂[p] : o.Y₂[p,:]
		o.haswt && (o.wt = o.wt[p])
		isdefined(o, :FEID) && nrows(o.FEID)>0 && (o.FEID = o.FEID[p])
  end

  if o.WREnonARubin
	  if !iszero(o.NClustVar)
	    o.infoBootData, o.IDBootData = panelsetupID(o.ID, 1:o.nbootclustvar)
	  else
	    o.info⋂Data = o.infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # no clustering, so no collapsing by cluster
	  end
  elseif !iszero(o.NClustVar)
	  o.infoBootData = panelsetup(o.ID, 1:min(o.NClustVar,o.nbootclustvar))  # bootstrap cluster grouping defs rel to original data
  else
	  info⋂Data = infoAllData = o.infoBootData = Vector{UnitRange{Int64}}(undef, o.Nobs, 0)  # causes no collapsing of data in panelsum() calls, only multiplying by weights if any
  end
  o.N✻ = nrows(o.infoBootData)

  if o.bootstrapt
    if o.NClustVar>0
	    minN = Inf; sumN = 0

	    combs = [x & 2^y > 0 for x in 2^o.nerrclustvar-1:-1:1, y in o.nerrclustvar-1:-1:0]  # represent all error clustering combinations. First is intersection of all error clustering vars
	    o.clust = Vector{StrClust{T}}(undef, nrows(combs))  # leave out no-cluster combination
	    o.NErrClustCombs = length(o.clust)

	    if o.NClustVar > o.nbootclustvar  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
		    if o.WREnonARubin && !o.granular
		      o.infoAllData, IDAllData = panelsetupID(o.ID, 1:o.NClustVar)
		    else
		      o.infoAllData            = panelsetup(o.ID, 1:o.NClustVar)
		    end
				if o.subcluster>0 && nrows(o.infoBootData) ≠ nrows(o.infoAllData)
					throw(ErrorException("\nThis program can only perform the subcluster bootstrap when the bootstrap clusters are nested within the (intersections of the) error clusters.\n"))
				end
			else
		    o.infoAllData = o.infoBootData  # info for grouping by intersections of all bootstrap && clustering vars wrt data; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusters
		    o.WREnonARubin && !o.granular && (IDAllData = o.IDBootData)
	    end

			Nall = length(o.infoAllData)

	    if o.NClustVar > o.nerrclustvar  # info for intersections of error clustering wrt data
		    if o.WREnonARubin && !o.granular
		      o.info⋂Data, ID⋂Data = panelsetupID(o.ID, o.subcluster+1:o.NClustVar)
		    else
		      o.info⋂Data            = panelsetup(o.ID, o.subcluster+1:o.NClustVar)
		    end
		    ID⋂ = length(o.info⋂Data)==o.Nobs ? o.ID : @views o.ID[first.(o.info⋂Data),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
		    o.IDAll = Nall==o.Nobs ? o.ID : @view o.ID[first.(o.infoAllData),:]  # version of ID matrix with one row for each all-bootstrap && error cluster-var intersection instead of 1 row for each obs
	    else
		    o.info⋂Data = o.infoAllData  # info for intersections of error clustering wrt data
		    o.WREnonARubin && !o.granular && (ID⋂Data = IDAllData)
		    o.IDAll = ID⋂ = nrows(o.info⋂Data)==o.Nobs ? o.ID : @views o.ID[first.(o.info⋂Data),:]  # version of ID matrix with one row for each all-error-cluster-var intersection instead of 1 row for each obs; gets resorted
	    end

	    o.BootClust = 2^(o.NClustVar - o.nbootclustvar)  # location of bootstrap clustering within list of cluster combinations

	    for c ∈ 1:o.NErrClustCombs  # for each error clustering combination
		    ClustCols = o.subcluster .+ findall(@view combs[c,:])
		    even = isodd(length(ClustCols))  # not a typo

		    if isone(c)
		      if iszero(o.subcluster)
		    	  order = Vector{Int64}(undef,0)
		    	  info  = Vector{UnitRange{Int64}}(undef, Nall)  # causes no collapsing of data in panelsum() calls
		      else
		    	  order = _sortperm(@view ID⋂[:,ClustCols])
		    	  ID⋂ = @view ID⋂[order, :]
		    	  info  = panelsetup(ID⋂, ClustCols)
		      end
		    else
		      if any(combs[c, min(findall(view(combs,c,:) .≠ view(combs,c-1,:))...):end])  # if this sort ordering same as last to some point and missing thereafter, no need to re-sort
		    	  order = _sortperm(@view ID⋂[:,ClustCols])
		    	  ID⋂ = @view ID⋂[order,:]
		      else
		    	  order = Vector{Int64}(undef,0)
		      end
		      info = panelsetup(ID⋂, ClustCols)
		    end

		    N = nrows(info)
		    sumN += N

				if o.small
					multiplier = T(N / (N-1))
					N < minN && (minN = N)
				else
					multiplier = one(T)
				end
				o.clust[c] = StrClust{T}(N, multiplier, even, order, info)
	    end

	    (o.scorebs || !o.WREnonARubin) &&
		  	(o.ClustShare = o.haswt ? @panelsum(o.wt, o.info⋂Data)/o.sumwt : length.(o.info⋂Data)./o.Nobs) # share of observations by group

    else  # if no clustering, cast "robust" as clustering by observation
      o.clust = StrClust{T}(Nobs, small ? _Nobs / (_Nobs - 1) : 1, true, Vector{Int64}(undef,0), Vector{UnitRange{Int64}}(undef,0))
      sumN = o.Nobs
      o.NErrClustCombs = 1
      (o.scorebs || !o.WREnonARubin) &&
    		(o.ClustShare = o.haswt ? o.wt/o.sumwt : 1/o._Nobs)
    end

		o.purerobust = o.robust && !o.scorebs && iszero(o.subcluster) && o.N✻==o.Nobs  # do we ever error-cluster *and* bootstrap-cluster by individual?
		o.granular   = o.WREnonARubin ? 2*o.Nobs*o.B*(2*o.N✻+1) < o.N✻*(o.N✻*o.Nobs+o.clust[1].N*o.B*(o.N✻+1)) :
											o.NClustVar>0 && !o.scorebs && (o.purerobust || (o.clust[1].N+o.N✻)*o.kZ*o.B + (o.clust[1].N-o.N✻)*o.B + o.kZ*o.B < o.clust[1].N*o.kZ^2 + o.Nobs*o.kZ + o.clust[1].N * o.N✻ * o.kZ + o.clust[1].N * o.N✻)

		if o.robust && !o.purerobust
			(o.subcluster>0 || o.granular) && !o.WREnonARubin && 
				(o.infoErrAll = panelsetup(o.IDAll, o.subcluster+1:o.NClustVar))  # info for error clusters wrt data collapsed to intersections of all bootstrapping && error clusters; used to speed crosstab UXAR wrt bootstrapping cluster && intersection of all error clusterings
			((o.scorebs && o.B>0) || (o.WREnonARubin && !o.granular && o.bootstrapt)) &&
				(o.JN⋂N✻ = zeros(T, o.clust[1].N, o.N✻))
		end

		if o.WREnonARubin && o.robust && o.bootstrapt && !o.granular
			if iszero(length(IDAllData))
				_, IDAllData = panelsetupID(o.ID,              1:o.NClustVar)
				_, ID⋂Data = panelsetupID(o.ID, o.subcluster+1:o.NClustVar)
			end
			o.IDCT⋂✻   = Vector{Vector{Int64}}(undef, o.N✻)
			o.infoCT⋂✻ = Vector{Vector{UnitRange{Int64}}}(undef, o.N✻)
			for i ∈ 1:o.N✻
				tmp = IDAllData[o.infoBootData[i]]                        # ID numbers w.r.t. intersection of all bootstrap/error clusterings contained in bootstrap cluster i
				o.infoCT⋂✻[i] = o.infoAllData[tmp[1]:tmp[end]]       # for each of those ID's, panel info for the all-bootstrap/error-clusterings data row groupings
				o.IDCT⋂✻[i] = ID⋂Data[first.(o.infoCT⋂✻[i])]  # ID numbers of those groupings w.r.t. the all-error-clusterings grouping
			end
		end
  else
	  minN = nrows(o.infoBootData)
	end

	if isdefined(o, :FEID) && nrows(o.FEID)>0
		p = _sortperm(o.FEID)
		sortID = o.FEID[p]
		i_FE = 1; o.FEboot = o.B>0 && !o.WREnonARubin && o.NClustVar>0; j = o.Nobs; o._FEID = ones(Int64, o.Nobs)
		o.invFEwt = zeros(T, o.NFE>0 ? o.NFE : o.Nobs)
		o.FEs = Vector{StrFE{T}}(undef, o.NFE>0 ? o.NFE : o.Nobs)
		@inbounds for i ∈ o.Nobs-1:-1:1
			if sortID[i] ≠ sortID[i+1]
				is = @view p[i+1:j]
				if o.haswt
					tmp  = o.wt[is]
					wt = tmp / (sumFEwt = sum(tmp))
				else
					sumFEwt = j - i
					wt = fill(1/sumFEwt, j-i)
				end
				o.FEs[i_FE] = StrFE{T}(is, wt)
				if (o.B>0 && o.robust && o.granular < o.nerrclustvar) || (o.WREnonARubin && o.robust && o.granular && o.bootstrapt)
					o.invFEwt[i_FE] = 1 / sumFEwt
				end

				j = i

				if o.FEboot  # are all of this FE's obs in same bootstrapping cluster? (But no need to check if B=0 for then CT_WE in 2nd term of (62) orthogonal to v = col of 1's)
					tmp = o.ID[is, 1:o.nbootclustvar]
					o.FEboot = all(tmp .== view(tmp, 1,:)')
				end
				i_FE += 1
			end
			o._FEID[p[i]] = i_FE
		end
		is = @view p[1:j]
		if o.haswt
			tmp = o.wt[is]
			wt = tmp / (sumFEwt = sum(tmp))
		else
			sumFEwt = j
			wt = fill(1/sumFEwt, j)
		end
		o.FEs[i_FE] = StrFE{T}(is, wt)
		o.robust && ((o.B>0 && o.granular < o.nerrclustvar) || (o.WREnonARubin && o.granular && o.bootstrapt)) &&
			(o.invFEwt[i_FE] = 1 / sumFEwt)
		o.NFE = i_FE
		resize!(o.invFEwt, o.NFE)
		resize!(o.FEs, o.NFE)
		if o.FEboot  # are all of this FE's obs in same bootstrapping cluster?
			tmp = o.ID[is, 1:o.nbootclustvar]
			o.FEboot = all(tmp .== @view tmp[1,:])
		end

		if o.robust && o.B>0 && o.bootstrapt && !o.FEboot && o.granular < o.nerrclustvar
			o.infoBootAll = panelsetup(o.IDAll, 1:o.nbootclustvar)  # info for bootstrapping clusters wrt data collapsed to intersections of all bootstrapping && error clusters
		end

		o.X₁ = partialFE(o, o.X₁)  # don't overwrite caller's data
		o.X₂ = partialFE(o, o.X₂)
		o.y₁ = partialFE(o, o.y₁)
		o.Y₂ = partialFE(o, o.Y₂)
	end

	if o.B>0 && o.robust && o.granular && !o.purerobust && o.bootstrapt && !o.WREnonARubin
		if o.NFE>0 && !o.FEboot
			_, o.IDBootData = panelsetupID(o.ID   , 1:o.nbootclustvar)
		else
			_, o.IDBootAll  = panelsetupID(o.IDAll, 1:o.nbootclustvar)
		end
	end

	o.enumerate = o.B>0 && o.auxtwtype==rademacher && o.N✻*log(2) < log(o.B)+1e-6  # generate full Rademacher set?
	o.enumerate && (o.maxmatsize = 0)

	o.Nw = iszero(o.maxmatsize) ? 1 : ceil((o.B+1) * Float64(max(nrows(o.IDBootData), length(o.IDBootAll), o.N✻) * sizeof(T)) / o.maxmatsize / 1073741824) # 1073741824 = giga(byte)
	if isone(o.Nw)
		MakeWildWeights!(o, o.B, first=true)  # make all wild weights, once
		o.enumerate && (o.B = ncols(o.v) - 1)  # replications reduced to 2^G
		o.WeightGrp = [1:ncols(o.v)]
	else
    o.seed = rand(o.rng,UInt64)
		_B = ceil(Int64, (o.B+1) / o.Nw)
		o.Nw = ceil(Int64, (o.B+1) / _B)
		o.WeightGrp = [(i-1)*_B+1:i*_B for i ∈ 1:o.Nw]
		o.WeightGrp[end] = first(o.WeightGrp[end]):o.B+1
	end

	if o.ML
		o.dof = nrows(o.R)
	else
		if o.ARubin
			o.R = hcat(zeros(o.kX₂,o.kX₁), Matrix(I(o.kX₂)))  # attack surface is all endog vars
			o.R₁ = o.kX₁>0 && nrows(o.R₁)>0 ? hcat(o.R₁[:,1:kX₁], zeros(nrows(o.R₁),o.kX₂)) : zeros(0, o.kX)  # and convert model constraints from referring to X₁, Y₂ to X₁, X₂
		end
		o.dof = nrows(o.R)

		if !o.WRE && iszero(o.κ)  # regular OLS
			o.DGP = StrEstimator{T}(o, true, o.LIML, o.Fuller, o.κ)
			o.Repl = StrEstimator{T}(o, true, false, zero(T), zero(T))  # XXX isDGP=1 for Repl? doesn't matter?
			setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁)  # DGP constraints: model constraints + null if imposed
			setR!(o.Repl, o.R₁)  # model constraints only
			InitVarsOLS!(o.DGP, o.Repl.R₁perp)
			InitTestDenoms!(o.DGP)
			o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
		elseif o.ARubin
			if o.willplot  # for plotting/CI purposes get original point estimate since not normally generated
				o.DGP = StrEstimator{T}(o, true, o.LIML, o.Fuller, o.κ)
				setR!(o.DGP, o.R₁, zeros(T,0,o.kZ))  # no-null model
				InitVarsIVGMM!(o.DGP)
				EstimateIVGMM!(o.DGP, o.r₁)
				o.confpeak = o.DGP.β̂  # estimated coordinate of confidence peak
			end

			o.DGP = StrEstimator{T}(o, true, false, zero(T), zero(T))
			setR!(o.DGP, o.R₁)
			InitVarsARubin!(o.DGP)
			InitTestDenoms!(o.DGP)
			o.M = o.DGP  # StrEstimator object from which to get A, AR, XAR
			o.kZ = o.kX

		elseif o.WREnonARubin

			o.DGP = StrEstimator{T}(o, true, T(o.kX₂ ≠ o.kY₂), zero(T), one(T))
			setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁, zeros(T,0,o.kZ))  # DGP constraints: model constraints + null if imposed
			InitVarsIVGMM!(o.DGP)
			if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
				EstimateIVGMM!(o.DGP, o.r₁)
				MakeResidualsIVGMM!(o.DGP)
			end

			o.Repl = StrEstimator{T}(o, false, o.LIML, o.Fuller, o.κ)
			setR!(o.Repl, o.R₁, o.R)
			InitVarsIVGMM!(o.Repl)
			EstimateIVGMM!(o.Repl, o.r₁)

			o.LIML && o.Repl.kZ==1 && o.Nw==1 && (o.As = o.β̂s = zeros(1, o.B+1))
			o.S✻UZperpinvZperpZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			o.S✻UZperp              = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			o.S✻uY                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			o.S✻UXinvXX             = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			o.S✻UX                  = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
			o.S✻UU                  = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.Repl.kZ+1)
			o.T1L = isone(o.Nw) ? [Matrix{T}(undef, o.Repl.kX, ncols(o.v))] :
									          [Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[1])), Matrix{T}(undef, o.Repl.kX, length(o.WeightGrp[end]))]
			o.T1R = deepcopy(o.T1L)

			if o.bootstrapt
				o.δdenom_b = zeros(o.Repl.kZ, o.Repl.kZ)
				o.S✻UMZperp = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
				o.S✻UPX     = Vector{Matrix{T}}(undef, o.Repl.kZ+1)
				o._J⋂ = zeros(o.clust[1].N, o.Repl.kZ)
				!o.granular && (o.SCT⋂uXinvXX = Matrix{Matrix{T}}(undef, o.Repl.kZ+1, o.N✻))
				if o.LIML || !o.robust
					o.YY✻_b   = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
					o.YPXY✻_b = zeros(o.Repl.kZ+1, o.Repl.kZ+1)
				end
				o.NFE>0 && (o.bootstrapt || !isone(o.κ) || o.LIML) && (o.CTFEU = Vector{Matrix{T}}(undef, o.Repl.kZ+1))
			end

		else  # the score bootstrap for IV/GMM uses a IV/GMM DGP but then masquerades as an OLS test because most factors are fixed during the bootstrap. To conform, need DGP and Repl objects with different R, R₁, one with FWL, one not

			o.DGP = StrEstimator{T}(o, true, o.LIML, o.Fuller, o.κ)
			setR!(o.DGP, o.null ? [o.R₁ ; o.R] : o.R₁, zeros(T,0,o.kZ))  # DGP constraints: model constraints + null if imposed
			InitVarsIVGMM!(o.DGP)
			o.Repl = StrEstimator{T}(o, true, o.LIML, o.Fuller, o.κ)
			setR!(o.Repl, o.R₁, I)  # process replication restraints = model constraints only
			InitVarsIVGMM!(o.Repl, o.Repl.R₁perp)
			EstimateIVGMM!(o.Repl, o.r₁)  # bit inefficient to estimate in both objects, but maintains the conformity
			InitTestDenoms!(o.Repl)
			o.M = o.Repl  # StrEstimator object from which to get A, AR, XAR; DGP follows WRE convention of using FWL, Repl follows OLS convention of not; scorebs for IV/GMM mixes the two
			if !o.null  # if not imposing null, then DGP constraints, κ, Hessian, etc. do not vary with r and can be set now
				EstimateIVGMM!(o.DGP, o.r₁)
				MakeResidualsIVGMM!(o.DGP)
			end
		end
  end

  if !o.WREnonARubin && o.bootstrapt
		o.denom = [Matrix{T}(undef,0,0) for _ in 1:o.dof, _ in 1:o.dof]
		if o.robust
			o.Kcd =                       Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof)
			o.Jcd = iszero(o.B) ? o.Kcd : Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof)  # if B = 0, Kcd will be multiplied by v, which is all 1's, and will constitute Jcd
		end

		if o.robust && o.granular<o.NErrClustCombs && o.B>0
			inds = o.subcluster>0 ?
				        [CartesianIndex(j, d, i) for d ∈ 1:o.dof for (j,v) ∈ enumerate(o.infoErrAll) for i ∈ v] :  # crosstab c,c* is wide
			       o.NClustVar == o.nbootclustvar ?
				        [CartesianIndex(i, d, i) for d ∈ 1:o.dof for i ∈ 1:Nall] :  # crosstab c,c* is square
			          [CartesianIndex(i, d, j) for d ∈ 1:o.dof for (j,v) ∈ enumerate(o.clust[o.BootClust].info) for i ∈ v]  # crosstab c,c* is tall
			o.crosstab⋂✻ind = LinearIndices(FakeArray(Tuple(max(inds...))...))[inds]
		end
	end

  o.small && (o.dof_r = o.NClustVar>0 ? minN - 1 : o._Nobs - o.kZ - o.NFE)

  o.sqrt = isone(o.dof)  # work with t/z stats instead of F/chi2?

  if o.small
		o.multiplier = (o.smallsample = (o._Nobs - o.kZ - o.FEdfadj * o.NFE) / (o._Nobs - o.robust)) / o.dof  # divide by # of constraints because F stat is so defined
  else
		o.multiplier = o.smallsample = 1
  end

  !(o.robust || o.ML) && (o.multiplier *= o._Nobs)  # will turn sum of squared errors in denom of t/z into mean
  o.sqrt && (o.multiplier = √o.multiplier)

  ((!o.bootstrapt && o.dof==1) || o.bootstrapt && (o.WREnonARubin || o.dof>1+o.robust || !isnan(o.maxmatsize))) &&  # unless nonWRE or dof=1 or splitting weight matrix, code will create dist element-by-element, so pre-allocate vector now
		(o.dist = fill(T(NaN), o.B+1))
  (o.Nw>1 || o.WREnonARubin || (!o.null && o.dof≤2)) && (o.numer = fill(T(NaN), o.dof, o.B+1))

  if o.WREnonARubin
		if o.Repl.kZ>1
			o.bootstrapt && o.robust &&
      	(o.Zyg = Vector{Matrix{T}}(undef,o.Repl.kZ))
			o.numer_b = Vector{T}(undef,nrows(o.Repl.RRpar))
		end
		o.bootstrapt && o.robust &&
    	(o.crosstabBootind = o.Nobs==o.N✻ ? diagind(FakeArray(o.N✻,o.N✻)) : 
			                                    LinearIndices(FakeArray(o.Nobs,o.N✻))[CartesianIndex.(1:o.Nobs, o.IDBootData)])
	else
		o.poles = o.anchor = zeros(T,0)
		o.interpolable = o.bootstrapt && o.B>0 && o.null && o.Nw==1 && (iszero(o.κ) || o.ARubin)
		if o.interpolable
			o.∂numer∂r = Vector{Matrix{T}}(undef, o.q)
			o.interpolate_u = !(o.robust || o.ML)
			o.interpolate_u && (o.∂u∂r = Vector{Matrix{T}}(undef, o.q))
			if o.robust
				o.∂denom∂r   = [Matrix{Matrix{T}}(undef, o.dof, o.dof) for _ in 1:o.q]
				o.∂²denom∂r² = [Matrix{Matrix{T}}(undef, o.dof, o.dof) for _ in 1:o.q, _ in 1:o.q]
				o.∂Jcd∂r     = [Matrix{Matrix{T}}(undef, o.NErrClustCombs, o.dof) for _ in 1:o.q]
			end
		end
  end
	o.initialized = true
	nothing
end


# draw wild weight matrix of width _B. If first=true, insert column of 1s at front. For non-Anderson-Rubin WRE, subtract 1 from all
const ϕ = (1 + √5)/2

function MakeWildWeights!(o::StrBootTest{T}, _B::Integer; first::Bool=true) where T
  if _B>0  # in scoretest or waldtest WRE, still make v a col of 1's
    if o.enumerate
			o.v = o.WREnonARubin ? [zeros(o.N✻) count_binary(o.N✻, -2, 0)] :  # complete Rademacher set
								[ones( o.N✻) count_binary(o.N✻, -1, 1)]
		elseif o.auxtwtype==normal
			o.v = randn(o.rng, T, o.N✻, _B+first)
			o.WREnonARubin && (o.v .-= one(T))
		elseif o.auxtwtype==gamma
			tmp = quantile.(Gamma(4,.5), rand(o.rng, o.N✻, _B+first))
			o.v = T==Float64 ? tmp : T.(tmp)
			o.WREnonARubin && (o.v .-= one(T))
		elseif o.auxtwtype==webb
			o.v = rand(o.rng, T.([-√1.5, -1, -√.5, √.5, 1, √1.5] .- o.WREnonARubin), o.N✻, _B+first)
		elseif o.auxtwtype == mammen
			o.v = getindex.(Ref(T.([1-ϕ; ϕ] .- o.WREnonARubin)), ceil.(Int16, rand(o.rng, o.N✻, _B+first) ./ (ϕ/√5)))
		elseif o.WREnonARubin  # Rademacher
			o.v = -2rand(o.rng, Bool, o.N✻, _B+first)
		else
			o.v = rand(o.rng, Bool, o.N✻, _B+first) .- T(.5)  # rand(o.rng, Bool, o.N✻, _B+first) .- T(.5)
			o.v_sd = .5
		end

		first && !(o.enumerate && isone(o.v_sd)) && (o.v[:,1] .= o.WREnonARubin ? zero(T) : o.v_sd)  # keep original residuals in first entry to compute base model stat
  else
		o.v = Matrix{T}(undef,0,1)  # in places, ncols(v) indicates B -- 1 for classical tests
  end
	nothing
end
