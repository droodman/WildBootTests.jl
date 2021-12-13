const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Matrix{Float64},Vector{Int64}})   # time: 29.61254
    let fbody = try __lookup_kwbody__(which(wildboottest, (Type,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 28.634747
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.1897423
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),Matrix{Float64},Vector{Int64}})   # time: 0.0375642
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0167941
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0094284
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0043929
end
