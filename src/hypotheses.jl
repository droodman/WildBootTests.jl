using StatFiles, StatsModels, DataFrames, DataFramesMeta, GLM
df = DataFrame(load(raw"d:\OneDrive\Documents\Macros\nlsw88.dta"))[:,[:wage; :tenure; :ttl_exp; :collgrad; :industry]]
dropmissing!(df)
desc = describe(df, :eltype)
for i in axes(desc, 1)  # needed only for lm()
  if desc[i,:eltype] == Float32
    sym = desc[i,:variable]
    @transform!(df, @byrow $sym=Float64($sym))
  end
end
f = @formula(wage ~ 1 + tenure + ttl_exp * collgrad)
f = apply_schema(f, schema(f, df))
resp, predexog = modelcols(f, df)

fit = lm(f, df)

isa(fit.model, LinearModel)
_terms = fit.mf.f.rhs.terms
k = length(_terms)
_I = zeros(k,k); _I[diagind(_I)] .= 1.

termcols = Dict(t => _I[:,i] for (i,t) in enumerate(_terms))
newtermcols = Dict{AbstractTerm, Vector{Float64}}()
for pair in termcols
  if isa(pair.first, InteractionTerm)
    newtermcols[InteractionTerm((pair.first.terms[2], pair.first.terms[1]))] = pair.second
  end
end
merge!(termcols, newtermcols)

Hstr = "wage + 2 * tenure ==3"
H = Meta.parse(Hstr)

dump(Meta.parse("1 + ttl_exp & collgrad"))

function processH(H::Expr)
  if H.head==:call
    if H.args[1] == :(==)
      return Expr(:(-), H.args[2], H.args[3])
    end
  end
end
