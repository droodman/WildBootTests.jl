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
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 9.424309
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 7.1170835
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 4.4328876
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 3.3402567
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int32}, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 1.0051409
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float64}, Bool, Float16, Symbol, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float64, Float64, Symbol, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 0.8410579
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float32}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.2917101
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Symbol, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.1963444
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.1256788
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.1190136
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0825599
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0763919
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0755425
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0686343
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0639479
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0631091
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0582189
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0557746
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.054092
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Vector{Int8}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0530406
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float64}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0522802
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0515514
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0507282
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Int64, Vector{Int8}, Bool, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0467907
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0462374
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype, :gridmin, :gridmax, :gridpoints), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Symbol, Bool, Symbol, Vector{Int64}, Vector{Int64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0455705
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0454792
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0454122
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0453752
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0450485
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0444796
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0444228
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0443829
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Float32}, Bool, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0440791
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0439299
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.043891
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0438733
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0437858
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0437764
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.043586
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0435672
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Symbol, Bool, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0435606
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0435263
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0435028
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0434664
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0434421
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0432034
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0430884
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0430132
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0429584
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0429279
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0428602
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0427947
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Symbol, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0427481
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0426704
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0426554
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0422199
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int8}, Int64, Int64, Vector{Int8}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0420223
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.041879
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Symbol, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0400881
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0400193
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0378436
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Bool, Symbol, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.035665
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0349056
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0347015
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0346204
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0345729
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0343215
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Bool, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0341215
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0340818
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0339565
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0338443
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0337248
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0329669
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0325287
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0321077
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :gridmax, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small, :gridmin, :gridpoints), Tuple{Vector{Float32}, Vector{Int64}, Symbol, Bool, Vector{Float32}, Symbol, Matrix{Float64}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Int64}, Vector{Int64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0320616
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0317482
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int16},Int64,Int64,Bool,Bool,Int64,Vector{Int16},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0314408
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0304135
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0302242
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0295772
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.029221
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0287836
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0287798
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int32},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0280387
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.027044
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0269575
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.026761
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Int64,Vector{Int16},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0265651
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Int8},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0263423
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Int8},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0262867
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0259589
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0258994
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0258914
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0258535
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int32},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0257285
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0255012
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0253583
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0253581
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0253155
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0252575
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int64},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0252214
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0250365
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0249334
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0249257
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0248242
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.024802
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0247368
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0247333
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0247094
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0244452
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0244329
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0243585
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Int64},Vector{Int64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.024358
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0241595
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0240106
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.0227377
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0226905
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Int64,Vector{Int64},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0226598
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0220486
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.021369
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.0189111
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.0183914
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int32}})   # time: 0.0171966
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.015809
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float32}, Bool, Float16, Symbol, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float32, Float32, Symbol, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0152942
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0148344
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0105446
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0104286
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0096984
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0095542
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0088032
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0078392
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Symbol, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0078214
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0077752
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Symbol, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0077631
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0076988
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.007581
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0075542
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.007317
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0071282
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0071056
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0070324
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0069828
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068888
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068853
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068156
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0067286
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0064627
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0063863
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0063593
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0063524
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.006324
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0061131
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0060602
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0060118
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Symbol, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0060052
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Symbol, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0060037
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0059884
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0058447
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0044944
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0040738
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.003453
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0033407
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0032591
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.003168
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0031518
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0028458
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0027897
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float64, 3},Array{Float64, 3},Vector{UnitRange{Int64}}})   # time: 0.0027608
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float32, 3},Array{Float32, 3},Vector{UnitRange{Int64}}})   # time: 0.0026849
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0024562
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0024532
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0021926
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0015101
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0014257
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0012034
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int32}})   # time: 0.0011392
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0010934
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010288
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010092
end
