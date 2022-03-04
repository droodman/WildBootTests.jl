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
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int64, Int64, Bool, Bool, Vector{Int64}, Bool, Vector{Float64}, Bool, Float16, PType, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float64, Float64, MAdjType, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float64}, DistStatType, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 1.1760427
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.9591153
    Base.precompile(Tuple{var"#172#threadsfor_fun#9"{Array{Float32, 3}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.3361888
    Base.precompile(Tuple{var"#172#threadsfor_fun#9"{Array{Float64, 3}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.2977641
    Base.precompile(Tuple{var"#220#threadsfor_fun#12"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0971249
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Vector{Int8}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.079906
    Base.precompile(Tuple{var"#220#threadsfor_fun#12"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0740018
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0726683
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.070277
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype, :gridmin, :gridmax, :gridpoints), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, AuxWtType, Bool, PType, Vector{Int64}, Vector{Int64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0657717
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.064222
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0620719
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0618202
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0617094
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float32}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0613731
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0609149
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0601135
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0592161
    Base.precompile(Tuple{var"#78#threadsfor_fun#4"{Float32, Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0575074
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0542725
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0534279
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int8}, Int64, Int64, Vector{Int8}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0530015
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0519324
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0517563
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0503532
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0500647
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Int64, Vector{Int8}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0496295
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0479905
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Float32}, Bool, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0479145
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0476936
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0472709
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0471536
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0469281
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0465767
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0459559
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0459384
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0454606
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0448377
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0446809
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0445768
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0438531
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0437793
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0436906
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0435365
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0434861
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0434237
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0430992
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0427381
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float64}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0425185
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0424628
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :gridpoints, :reps), Tuple{BitVector, Matrix{Float64}, Vector{Int64}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0423528
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0417798
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0415289
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0412585
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0412343
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0397971
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.038659
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0380971
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.03371
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0323118
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0288406
    Base.precompile(Tuple{var"#78#threadsfor_fun#4"{Float64, Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.02808
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0278399
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0263511
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.025937
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0234222
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Bool, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0230071
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0229733
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.0229164
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0223521
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0221212
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :bootstrapc, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Bool, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0221179
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0209407
    Base.precompile(Tuple{var"#307#threadsfor_fun#17"{StrBootTest{Float32}, Vector{UnitRange{Int64}}, Vector{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0204411
    Base.precompile(Tuple{var"#340#threadsfor_fun#19"{StrBootTest{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0203299
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0190884
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0190212
    Base.precompile(Tuple{var"#340#threadsfor_fun#19"{StrBootTest{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0188126
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0187744
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.018042
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0178043
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0175534
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0174569
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0170534
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.016782
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0166998
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0165059
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0163984
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int64, Int64, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0163561
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0163541
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0161929
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :gridmax, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small, :gridmin, :gridpoints), Tuple{Vector{Float32}, Vector{Int64}, AuxWtType, Bool, Vector{Float32}, PType, Matrix{Float64}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Int64}, Vector{Int64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0159534
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0156369
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0153287
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0152542
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.014722
    Base.precompile(Tuple{var"#307#threadsfor_fun#17"{StrBootTest{Float64}, Vector{UnitRange{Int64}}, Vector{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0142109
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0125035
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0123775
    Base.precompile(Tuple{var"#269#threadsfor_fun#15"{StrBootTest{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0120511
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0120028
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0108486
    Base.precompile(Tuple{var"#236#threadsfor_fun#13"{StrBootTest{Float64}, Vector{UnitRange{Int64}}, Vector{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0107813
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0105288
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0103064
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0102787
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.010278
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0102124
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0101977
    Base.precompile(Tuple{var"#269#threadsfor_fun#15"{StrBootTest{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0101295
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0099126
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0098244
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0097863
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0094458
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0092473
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.009127
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.009043
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0090402
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0089448
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0089181
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0088485
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0088345
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0087651
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.00869
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0086492
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0086427
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0084702
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0084631
    Base.precompile(Tuple{var"#236#threadsfor_fun#13"{StrBootTest{Float32}, Vector{UnitRange{Int64}}, Vector{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0084534
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.008259
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0082311
    Base.precompile(Tuple{var"#188#threadsfor_fun#10"{Array{Float32, 3}, Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0081618
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0080555
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0079749
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0079651
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0079278
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0079121
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0077094
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0076955
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0076734
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0076701
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0076532
    Base.precompile(Tuple{var"#204#threadsfor_fun#11"{Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0076422
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0076321
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.007625
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0076202
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0075869
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0075775
    Base.precompile(Tuple{var"#61#threadsfor_fun#1"{SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0075525
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0075402
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0075306
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0075191
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.007518
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0074532
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0074465
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Vector{Float64}, AuxWtType, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.007411
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0073499
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0072912
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0071896
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.007182
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.007178
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0071752
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0070692
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0070255
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0068882
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
    end
end   # time: 0.0068637
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0068122
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0067273
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0066351
    Base.precompile(Tuple{var"#204#threadsfor_fun#11"{Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0066104
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.006484
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0064721
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0064681
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0064454
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0063684
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0062879
    Base.precompile(Tuple{var"#188#threadsfor_fun#10"{Array{Float64, 3}, Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0059578
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0058685
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0058546
    Base.precompile(Tuple{var"#150#threadsfor_fun#8"{Array{Float32, 3}, Array{Float32, 3}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0057326
    Base.precompile(Tuple{var"#95#threadsfor_fun#5"{Float32, Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0057122
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.005554
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0055059
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0054396
    Base.precompile(Tuple{var"#95#threadsfor_fun#5"{Float64, Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.00532
    Base.precompile(Tuple{var"#61#threadsfor_fun#1"{SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0052754
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.005236
    Base.precompile(Tuple{var"#78#threadsfor_fun#4"{Float64, Adjoint{Float64, Vector{Float64}}, Int64, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0051943
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float32, 3},Array{Float32, 3},Vector{UnitRange{Int64}}})   # time: 0.0051594
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0050898
    Base.precompile(Tuple{var"#114#threadsfor_fun#6"{Vector{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0049825
    Base.precompile(Tuple{var"#78#threadsfor_fun#4"{Float32, Adjoint{Float32, Vector{Float32}}, Int64, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0049499
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0048271
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0047944
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0045993
    Base.precompile(Tuple{var"#150#threadsfor_fun#8"{Array{Float64, 3}, Array{Float64, 3}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0045694
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0040639
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float64, 3},Array{Float64, 3},Vector{UnitRange{Int64}}})   # time: 0.0040484
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0038147
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0037492
    Base.precompile(Tuple{var"#132#threadsfor_fun#7"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0036898
    Base.precompile(Tuple{var"#114#threadsfor_fun#6"{Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0035959
    Base.precompile(Tuple{var"#114#threadsfor_fun#6"{Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0035386
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0035273
    Base.precompile(Tuple{var"#114#threadsfor_fun#6"{Vector{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0034955
    Base.precompile(Tuple{var"#61#threadsfor_fun#1"{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0034938
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0033488
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0031236
    Base.precompile(Tuple{var"#61#threadsfor_fun#1"{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0027655
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0017025
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0015431
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0013642
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0012236
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0011558
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010695
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010404
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0010074
end
