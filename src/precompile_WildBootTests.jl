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
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 7.384286
    Base.precompile(Tuple{typeof(coldotplus_turbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 2.500972
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 1.4456258
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.9750638
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float64}, Bool, Float16, PType, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float64, Float64, MAdjType, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 0.91104
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.8067636
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.666288
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.6174332
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.6012798
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.4452934
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.4330556
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.3836992
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.363491
    Base.precompile(Tuple{typeof(coldotplus_turbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.3479443
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.2722983
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.1743681
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0570123
    isdefined(WildBootTests, Symbol("#5#6")) && Base.precompile(Tuple{getfield(WildBootTests, Symbol("#5#6")),SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0533983
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0520964
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0455373
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.043764
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int32}})   # time: 0.0433843
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0425812
    isdefined(WildBootTests, Symbol("#5#6")) && Base.precompile(Tuple{getfield(WildBootTests, Symbol("#5#6")),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0424144
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0406334
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0389851
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0385406
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.037189
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0354867
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0349768
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0348549
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0347912
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0342838
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0337724
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0323543
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0319488
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0318842
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0317586
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0310479
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0308801
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0307885
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0307542
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0306166
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0299047
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.029723
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0289742
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0288412
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0287667
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0281132
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0276471
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0272731
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0262976
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0262113
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0253069
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.024738
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0241874
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.023958
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0193806
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.017551
    Base.precompile(Tuple{typeof(coldotplus_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0149678
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0116493
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0116131
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0108701
    Base.precompile(Tuple{typeof(coldotplus_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0104375
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0098953
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0094551
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Adjoint{Float32, Vector{Float32}}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0092785
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0086905
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0084532
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0082392
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0081392
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0080218
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0080189
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0079765
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0079355
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0079051
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0077141
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0076637
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.00764
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0075996
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0075994
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0075087
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0074761
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0074607
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0074573
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0074423
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0074143
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0073953
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0073914
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0073854
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0073529
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float32}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.007342
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0072671
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0072626
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0072398
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0072317
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.007188
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0071586
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float64}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.007118
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0070688
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070548
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0069394
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0069274
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.006921
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0069168
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0069036
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068803
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0068752
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0068746
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068692
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068628
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068438
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0067989
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0067943
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0067474
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Adjoint{Float64, Vector{Float64}}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0067424
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0067322
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0067128
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0066089
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0065626
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0065271
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0065238
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0064848
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0064835
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0064675
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0064569
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0064243
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0064011
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0063515
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0063231
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0062897
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0062747
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0062575
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0062296
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.006212
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0062015
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0061511
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0061264
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.006109
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0060649
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.006026
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0060203
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0060203
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0059863
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0059799
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0057983
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0057843
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0057228
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0056936
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0056928
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0056676
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0055455
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0055296
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0054927
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0054812
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.005287
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0051644
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0050048
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0049917
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0047952
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0047049
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0045985
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0045779
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0045418
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.004447
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0043522
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0043072
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0040963
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0040637
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0039611
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0038154
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0037924
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.003771
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0036673
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0036586
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0036565
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0035769
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0034827
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.003455
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.003378
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0033606
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0032939
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0032667
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Vector{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0032371
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0032292
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0031693
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0031429
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0031266
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0030679
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.003033
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.002948
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0029396
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.002681
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0026617
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0026395
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Vector{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0023104
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0013688
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0013143
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0010054
end
