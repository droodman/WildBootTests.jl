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
    let fbody = try __lookup_kwbody__(which(sort!, (Matrix{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sort!),Matrix{Float64},))
        end
    end   # time: 0.0503258
    let fbody = try __lookup_kwbody__(which(sort!, (Matrix{Float32},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sort!),Matrix{Float32},))
        end
    end   # time: 0.0485605
    let fbody = try __lookup_kwbody__(which(sort!, (Vector{Int16},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sort!),Vector{Int16},))
        end
    end   # time: 0.0089843
    let fbody = try __lookup_kwbody__(which(sort!, (Vector{Int8},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sort!),Vector{Int8},))
        end
    end   # time: 0.0086583
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{SubArray{Int64, 1, Matrix{Int64}, Tuple{Int64, Vector{Int64}}, false}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sortperm),Vector{SubArray{Int64, 1, Matrix{Int64}, Tuple{Int64, Vector{Int64}}, false}},))
        end
    end   # time: 0.0076983
    let fbody = try __lookup_kwbody__(which(sortperm, (Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (QuickSortAlg,Function,Function,Nothing,ForwardOrdering,typeof(sortperm),Vector{Int64},))
        end
    end   # time: 0.0053763
end
