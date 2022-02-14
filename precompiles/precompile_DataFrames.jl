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
    Base.precompile(Tuple{typeof(describe),DataFrame,Symbol})   # time: 0.3938833
    Base.precompile(Tuple{Core.kwftype(typeof(fromcolumns)),NamedTuple{(:copycols,), Tuple{Nothing}},typeof(fromcolumns),Tables.CopiedColumns{NamedTuple{(:year, :selfemployed, :hasinsurance, :post, :post_self), Tuple{Vector{Union{Missing, Float32}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}}}},Vector{Symbol}})   # time: 0.2365463
    Base.precompile(Tuple{typeof(dropmissing!),DataFrame})   # time: 0.0704978
    isdefined(DataFrames, Symbol("#DataFrame#152#154")) && Base.precompile(Tuple{getfield(DataFrames, Symbol("#DataFrame#152#154")),Bool,Type{DataFrame},Vector{Any},Index})   # time: 0.0639164
    isdefined(DataFrames, Symbol("#74#82")) && Base.precompile(Tuple{getfield(DataFrames, Symbol("#74#82")),Vector{Float32}})   # time: 0.0562037
    isdefined(DataFrames, Symbol("#74#82")) && Base.precompile(Tuple{getfield(DataFrames, Symbol("#74#82")),Vector{Int8}})   # time: 0.0461868
    Base.precompile(Tuple{typeof(Tables.schema),DataFrameColumns{DataFrame}})   # time: 0.0378148
    isdefined(DataFrames, Symbol("#DataFrame#152#154")) && Base.precompile(Tuple{getfield(DataFrames, Symbol("#DataFrame#152#154")),Bool,Type{DataFrame},Vector{AbstractVector{T} where T},Index})   # time: 0.0124559
    let fbody = try __lookup_kwbody__(which(manipulate, (DataFrame,Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,Bool,Bool,typeof(manipulate),DataFrame,Vector{Int64},))
        end
    end   # time: 0.0098782
    let fbody = try __lookup_kwbody__(which(describe, (DataFrame,Symbol,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Function,typeof(describe),DataFrame,Symbol,))
        end
    end   # time: 0.0077509
    Base.precompile(Tuple{typeof(get_stats),Union{Base.SkipMissing, AbstractVector{T} where T},Vector{Symbol}})   # time: 0.0071073
    Base.precompile(Tuple{Core.kwftype(typeof(fromcolumns)),NamedTuple{(:copycols,), Tuple{Nothing}},typeof(fromcolumns),Tables.CopiedColumns{NamedTuple{(:coll, :merit, :male, :black, :asian, :year, :state, :chst), NTuple{8, Vector{Union{Missing, Float32}}}}},Vector{Symbol}})   # time: 0.0053007
    Base.precompile(Tuple{Type{Matrix{Int8}},DataFrame})   # time: 0.0052499
    Base.precompile(Tuple{Core.kwftype(typeof(fromcolumns)),NamedTuple{(:copycols,), Tuple{Nothing}},typeof(fromcolumns),Tables.CopiedColumns{NamedTuple{(:idcode, :age, :race, :married, :never_married, :grade, :collgrad, :south, :smsa, :c_city, :industry, :occupation, :union, :wage, :hours, :ttl_exp, :tenure), Tuple{Vector{Union{Missing, Int16}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Int8}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}}}},Vector{Symbol}})   # time: 0.0052133
    let fbody = try __lookup_kwbody__(which(make_unique!, (Vector{Symbol},Vector{Symbol},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,typeof(make_unique!),Vector{Symbol},Vector{Symbol},))
        end
    end   # time: 0.0051304
    Base.precompile(Tuple{Core.kwftype(typeof(fromcolumns)),NamedTuple{(:copycols,), Tuple{Nothing}},typeof(fromcolumns),Tables.CopiedColumns{_A} where _A,Any})   # time: 0.005003
    Base.precompile(Tuple{Core.kwftype(typeof(fromcolumns)),NamedTuple{(:copycols,), Tuple{Nothing}},typeof(fromcolumns),Tables.CopiedColumns{NamedTuple{(:name, :id_split, :fid_843_gr, :pixcluster, :pixwbcode, :pixlat, :pixlon, :pixwater, :pixnmbrdia, :pixnmbrpetro, :lnkm, :pixpetro, :pixdia, :pixwaterd, :pixcapdist, :pixmal, :pixsead, :pixsuit, :pixelev, :pixbdist, :pd0, :lnpd0, :l0708, :l0708d, :lnl0708s, :lnwaterkm, :km2split, :lnkm2split, :wbcode, :mean_elev, :mean_suit, :malariasuit, :petroleum, :diamondd, :capdistance1, :seadist1, :borderdist1, :centr_tribe, :gr, :unique, :centr_tribe1), Tuple{Vector{Union{Missing, String}}, Vector{Union{Missing, Int16}}, Vector{Union{Missing, Int32}}, Vector{Union{Missing, String}}, Vector{Union{Missing, String}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, String}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, Float32}}, Vector{Union{Missing, String}}, Vector{Union{Missing, Float32}}}}},Vector{Symbol}})   # time: 0.0047691
    Base.precompile(Tuple{typeof(getindex),DataFrame,Colon,Vector{Symbol}})   # time: 0.0040114
    Base.precompile(Tuple{Type{Matrix{Int16}},DataFrame})   # time: 0.0035682
    Base.precompile(Tuple{typeof(setindex!),DataFrame,Vector{DataType},typeof(!),Symbol})   # time: 0.0023972
    Base.precompile(Tuple{typeof(setproperty!),DataFrame,Symbol,Vector{Int16}})   # time: 0.0017911
    let fbody = try __lookup_kwbody__(which(select, (DataFrame,Any,))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,Bool,typeof(select),DataFrame,Any,))
        end
    end   # time: 0.0017687
    Base.precompile(Tuple{typeof(nrow),DataFrame})   # time: 0.0015653
    let fbody = try __lookup_kwbody__(which(manipulate, (DataFrame,Vector{Symbol},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Bool,Bool,Bool,typeof(manipulate),DataFrame,Vector{Symbol},))
        end
    end   # time: 0.0012022
    Base.precompile(Tuple{typeof(setindex!),DataFrame,Vector{Int8},typeof(!),Int64})   # time: 0.0011949
    Base.precompile(Tuple{typeof(setindex!),DataFrame,Vector{String},typeof(!),Int64})   # time: 0.0011189
    Base.precompile(Tuple{typeof(setindex!),DataFrame,Vector{Float32},typeof(!),Int64})   # time: 0.0010698
end
