# given a pre-configured boottest linear model with one-degree null imposed, compute boostrapped p value associated with r
function r_to_p(o::StrBootTest{T}, r::AbstractVector{T}) where T
  o.r = r
  getp(o)
end


# Chandrupatla 1997, "A new hybrid quadratic-bisection algorithm for finding the zero of a nonlinear function without using derivatives"
# x₁, x₂ must bracket the true value, with f₁=f(x₁) and f₂=f(x₂)
function search(o::StrBootTest{T}, α::T, f₁::T, x₁::T, f₂::T, x₂::T) where T<:Real
  t = half = T(.5)
  while true
		x = x₁ + t * (x₂ - x₁)
		fx = r_to_p(o, [x])
		((fx>f₁) == (fx>f₂)) && return x  # violation of monotonicity because of precision problems? That's as good as it gets.

		if (fx<α) == (f₁<α)
			x₃, x₁, f₃, f₁ = x₁, x, f₁, fx
		else
			x₃, x₂, x₁, f₃, f₂, f₁ = x₂, x₁, x, f₂, f₁, fx
		end

		((o.B>0 && abs(fx - α) < (1+(o.ptype==equaltail)) / o.BFeas * 1.000001) || ≈(x₂, x₁, rtol=o.rtol)) &&
			return abs(f₁ - α) < abs(f₂ - α) ? x₁ : x₂

		ϕ₁ = (f₁ - f₂) / (f₃ - f₂)
		ϕ₁² = ϕ₁^2
		xi1 = (x₁ - x₂) / (x₃ - x₂)
		t = ϕ₁² > xi1 || xi1 > 2 * ϕ₁ - ϕ₁² ?
					half :
					clamp(((f₃ - α) / (f₁ - f₂) + (x₃ - x₁) / ((x₂ - x₁) * (f₃ - f₁)) * (f₂ - α)) * (f₁ - α) / (f₃ - f₂), T(0.000001), T(0.999999))
  end
end


# derive wild bootstrap-based CI, for case of linear model with one-degree null imposed
# and generate plot data
function plot!(o::StrBootTest{T}) where T
  _r = copy(o.r)
  α = one(T) - o.level

	iszero(nrows(o.gridmin   )) && (o.gridmin    = fill(T(NaN), o.q))
	iszero(nrows(o.gridmax   )) && (o.gridmax    = fill(T(NaN), o.q))
	iszero(nrows(o.gridpoints)) && (o.gridpoints = fill(Float32(NaN), o.q))

	o.gridpoints[isnan.(o.gridpoints)] .= 25

  o.boottest!(o)

	Phi = quantile(Normal{T}(zero(T),one(T)), α/2)
  if o.ARubin
		p = o.dist[1] * o.multiplier
		p = ccdf(Chisq{T}(T(o.dof)), o.sqrt ? p^2 : p)
		halfwidth = abs.(o.confpeak) * quantile(Normal{T}(zero(T),one(T)), p/2) / Phi
  else
		halfwidth = T(-1.5) * Phi .* sqrtNaN.(diag(getV(o)))
		o.confpeak = getb(o) + o.r
  end

	if isone(o.q)  # 1D plot
		α≤0 && (α = T(.05))  # if level=100, no CI constructed, but we need a reasonable α to choose graphing bounds
		p_lo = p_hi = T(NaN)
		if isnan(o.gridmin[1]) || isnan(o.gridmax[1])
			if o.B>0  # initial guess based on classical distribution
				lo = isnan(o.gridmin[1]) ? o.confpeak - halfwidth : o.gridmin
				hi = isnan(o.gridmax[1]) ? o.confpeak + halfwidth : o.gridmax
			else
				tmp = vec(sqrtNaN.(o.statDenom)) * (o.small ? cquantile(TDist{T}(o.dof_r), α/2) : cquantile(Normal{T}(zero(T),one(T)), α/2))
				lo = isnan(o.gridmin[1]) ? o.confpeak - tmp : o.gridmin
				hi = isnan(o.gridmax[1]) ? o.confpeak + tmp : o.gridmax
				if o.scorebs && !o.null && !o.willplot  # if doing simple Wald test with no graph, we're done
					o.CI = [lo hi]
					return
				end
			end

			if abs(lo[1] - o.r[1]) > abs(hi[1] - o.r[1])  # brute force way to ensure that first trial bound tested is the farther one from r, for better interpolation
				if isnan(o.gridmin[1]) && o.ptype≠lower  # unless lower-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
					for _ ∈ 1:10
						p_lo = r_to_p(o, lo)
						p_lo < α && break
						diff = hi - lo
						lo .-= diff
						isnan(o.gridmax[1]) && o.twotailed && (hi .+= diff)  # maintain rough symmetry unless user specified upper bound
					end
				end
				if isnan(o.gridmax[1]) && o.ptype≠upper  # ditto for high side
					for _ ∈ 1:10
						p_hi = r_to_p(o, hi)
						p_hi < α && break
						diff = hi - lo
						isnan(o.gridmin[1]) && o.twotailed && (lo .-= diff)
						hi .+= diff
					end
				end
			else
				if isnan(o.gridmax[1]) && o.ptype≠upper  # ditto for high side
					for _ ∈ 1:10
						p_hi = r_to_p(o, hi)
						p_hi < α && break
						diff = hi - lo
						isnan(o.gridmin[1]) && o.twotailed && (lo .-= diff)
						hi .+= diff
					end
				end
				if isnan(o.gridmin[1]) && o.ptype≠lower  # unless upper-tailed p value, try at most 10 times to bracket confidence set by symmetrically widening
					for _ ∈ 1:10
						p_lo = r_to_p(o, lo)
						p_lo < α && break
						diff = hi - lo
						lo .-= diff
						isnan(o.gridmax[1]) && o.twotailed && (hi .+= diff)  # maintain rough symmetry unless user specified upper bound
					end
				end
			end
		else  # both grid bounds pre-specified
			lo = [o.gridmin[1]]
			hi = [o.gridmax[1]]
		end

		_gridpoints = round.(Int32, o.gridpoints)
		o.plotX = [collect(range(lo[1], hi[1], length=_gridpoints[1]))]  # store as 1-vector since q=1
		o.plotY = fill(T(NaN), _gridpoints[1])
		o.plotY[1  ] = p_lo
		o.plotY[end] = p_hi
		p_confpeak = o.WREnonARubin ? T(NaN) : o.twotailed ? one(T) : T(.5)

		c = clamp((floor(Int, (o.confpeak[1] - lo[1]) / (hi[1] - lo[1]) * (_gridpoints[1] - 1)) + 2), 1, _gridpoints[1]+1)  # insert original point estimate into grid
		insert!(o.plotX[1], c, o.confpeak[1])
		insert!(o.plotY   , c, p_confpeak)

		isnan(o.plotY[1]) && (o.plotY[1] = r_to_p(o, [o.plotX[1][1]]))  # do in this order for widest interpolation base range
		for (i,v) ∈ Iterators.reverse(enumerate(o.plotX[1]))
			i>1 && isnan(o.plotY[i]) && (o.plotY[i] = r_to_p(o, [v]))
		end

	else  # 2D plot
		_gridpoints = round.(Int32, o.gridpoints)

		lo = Vector{T}(undef, 2)
		hi = Vector{T}(undef, 2)
		for d ∈ 1:o.dof
			lo[d] = isnan(o.gridmin[d]) ? o.confpeak[d] - halfwidth[d] : o.gridmin[d]
			hi[d] = isnan(o.gridmax[d]) ? o.confpeak[d] + halfwidth[d] : o.gridmax[d]
		end
		o.plotX = [collect(range(lo[1], hi[1], length=_gridpoints[1])), 
		           collect(range(lo[2], hi[2], length=_gridpoints[2]))]
		o.plotY = fill(T(NaN), _gridpoints[1]*_gridpoints[2])

		isnan(o.plotY[1]) && (o.plotY[1] = r_to_p(o, [o.plotX[i][1] for i in 1:o.q]))  # do in this order for widest interpolation base range
		for (i,v) ∈ Iterators.reverse(enumerate(Iterators.product(o.plotX[1],o.plotX[2])))
			i>1 && isnan(o.plotY[i]) && (o.plotY[i] = r_to_p(o, [v...]))
		end
	end

	if any(isnan.(o.plotY))
		o.CI = [T(-Inf) T(Inf)]
	elseif isone(o.q) && o.level<100 # find CI bounds
		_CI = Vector{T}(undef, nrows(o.plotY))
		for i in eachindex(_CI)  # map() version hampers type inference in Julia 1.6.2
			_CI[i] = isnan(o.plotY[i]) ? o.plotY[i] : T(o.plotY[i] > α)
		end
		_CI = _CI[2:end] - _CI[1:end-1]
		lo = T.(findall(x->x== 1, _CI))
		hi = T.(findall(x->x==-1, _CI))
		if iszero(length(lo)) && iszero(length(hi))
			o.CI = [T(-Inf) T(Inf)]
		else
			if iszero(length(lo))
				lo = [T(-Inf)]
			elseif iszero(length(hi))
				hi = [T(Inf)]
			else
				lo[1  ] > hi[1  ] && (lo = [T(-Inf) ; lo    ]) # non-rejection ranges that are not within grid range
				lo[end] > hi[end] && (hi = [hi      ; T(Inf)])
			end
			o.CI = [lo hi]

			for i ∈ 1:length(lo), j ∈ 1:2
				if !isinf(o.CI[i,j])
					t = Int(o.CI[i,j])
					o.CI[i,j] = search(o, α, o.plotY[t], o.plotX[1][t], o.plotY[t+1], o.plotX[1][t+1])
				end
			end
		end
		nothing
	end

  if @isdefined c  # now that it's done helping 1-D graph look good, remove peak point from returned grid for evenness, for Bayesian sampling purposes
		o.peak = (X = [o.plotX[1][c]], p = o.plotY[c])
		deleteat!(o.plotX[1], c)
		deleteat!(o.plotY   , c)
  end

	o.r = _r  # restore backup
	o.notplotted = false
	nothing
end