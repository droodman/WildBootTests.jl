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
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 10.142595
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 6.838591
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 3.5100586
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 3.486345
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 1.166926
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int64, Int64, Bool, Bool, Vector{Int64}, Bool, Vector{Float64}, Bool, Float16, PType, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float64, Float64, MAdjType, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float64}, DistStatType, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 0.8081782
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float32}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.2994384
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.2909995
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.2358318
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0542653
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0540688
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.053942
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype, :gridmin, :gridmax, :gridpoints), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, AuxWtType, Bool, PType, Vector{Int64}, Vector{Int64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0494241
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float64}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0458381
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Int64, Vector{Int8}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.044926
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0446042
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0445613
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0444836
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0444415
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0443967
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0443018
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0442717
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0441478
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0440003
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.043695
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.041654
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.041264
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0411461
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0410104
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0407966
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0407937
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Float32}, Bool, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0404564
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0404003
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.040321
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0401969
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0400682
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.039985
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0397801
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.039773
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0397226
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0396545
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0396309
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0396076
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int8}, Int64, Int64, Vector{Int8}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0393941
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.039388
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0393363
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0389003
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0387108
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0385769
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0385503
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0370368
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0354654
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0341098
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Vector{Int8}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0336724
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int64, Int64, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0310178
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0292368
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.027108
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0265656
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0259375
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0254639
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.0232933
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.0225032
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0220379
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0210451
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.020798
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.0202629
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Bool, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0196278
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0194216
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0193154
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0192693
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.0190526
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0176583
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int32}})   # time: 0.0176205
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0174285
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.017264
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0167448
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0160921
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0158118
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.015796
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0157214
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Bool, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0156146
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0150576
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0150508
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0150365
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0150219
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0147477
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0146242
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :gridmax, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small, :gridmin, :gridpoints), Tuple{Vector{Float32}, Vector{Int64}, AuxWtType, Bool, Vector{Float32}, PType, Matrix{Float64}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Int64}, Vector{Int64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0146101
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0145945
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0145188
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0143801
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0127131
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0118659
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0117199
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0112982
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.010724
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0104409
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0100291
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0100164
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0099857
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0098922
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.009886
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0095972
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.009
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0089558
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0086957
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.008658
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0085233
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0085056
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0084637
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0084628
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0083897
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0083311
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0083299
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0082275
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0081576
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.007831
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0077957
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0076506
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0076463
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0075363
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0074696
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0074519
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0074096
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0071953
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0071884
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0071779
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0071444
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0071261
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0070466
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0070388
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0070258
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0070015
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.006985
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0069833
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0069726
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0069305
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0069275
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0068882
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0068473
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0067601
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0067516
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0067141
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0067044
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0066575
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0066539
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0066452
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0065806
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0065091
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0064631
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0064377
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0063104
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.006299
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0062624
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.006135
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0060406
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0059604
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0058467
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0058065
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0057494
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0057213
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0048821
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0048394
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0045202
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.00392
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0034946
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0034847
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0034836
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0033679
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0032519
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float64, 3},Array{Float64, 3},Vector{UnitRange{Int64}}})   # time: 0.0032406
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.002824
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0027506
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float32, 3},Array{Float32, 3},Vector{UnitRange{Int64}}})   # time: 0.0025863
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0018161
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0015376
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0014931
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0013998
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0013271
end
