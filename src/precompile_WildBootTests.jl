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
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 7.669295
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 6.8374147
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 5.6122956
    Base.precompile(Tuple{typeof(coldotplus_turbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 2.6222863
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 2.5752935
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 2.0898058
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 1.9083456
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 1.4733161
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float64}, Bool, Float16, PType, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float64, Float64, MAdjType, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 1.1085432
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.851465
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.7183718
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.6621906
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.5566583
    Base.precompile(Tuple{typeof(colquadformminus_turbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.5209184
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.4769747
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.4556499
    Base.precompile(Tuple{typeof(coldotplus_turbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.3960819
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.3618275
    Base.precompile(Tuple{typeof(matmulplus_turbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.3509651
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.272368
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.1947931
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.1682014
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.1348114
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.1127438
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Vector{Int8}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.1078734
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.1059839
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0717236
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float32}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0632265
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0559807
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0536664
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float64}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0508321
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0492981
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype, :gridmin, :gridmax, :gridpoints), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, AuxWtType, Bool, PType, Vector{Int64}, Vector{Int64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0429594
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0411336
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0403518
    Base.precompile(Tuple{typeof(panelsum_turbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0370208
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0367734
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0365918
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0361303
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0359108
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0354628
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :imposenull), Tuple{Bool, BitVector, Matrix{Float64}, Vector{Int64}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0354159
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0343173
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0342058
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0336494
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.033589
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0334645
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0333738
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0332957
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0329933
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0328417
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0326096
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0324992
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.032332
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.032191
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0318666
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, BitVector, Matrix{Float64}, Vector{Int64}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.03141
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0309138
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0307178
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0306998
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int8}, Int64, Int64, Vector{Int8}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0305718
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0305214
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float64}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0302963
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0302612
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0301704
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0300879
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Int64, Vector{Int8}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0300079
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.029972
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0298903
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0298481
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0297988
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0296938
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0296576
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0296558
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0294805
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0293815
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0293169
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Float32}, Bool, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0290433
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0285769
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0277658
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0269433
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0268919
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0268369
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0267735
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0250034
    Base.precompile(Tuple{typeof(panelsum_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0228285
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.0221804
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.022066
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.0212147
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.0207681
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int32}})   # time: 0.0198441
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float32}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0178588
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0168542
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.016527
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0154042
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :gridmax, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small, :gridmin, :gridpoints), Tuple{Bool, Vector{Float32}, Vector{Int64}, AuxWtType, Bool, Vector{Float32}, PType, Matrix{Float64}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Int64}, Vector{Int64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0148782
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0147169
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.014423
    Base.precompile(Tuple{typeof(coldotplus_turbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0142845
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0142721
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0140695
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0140062
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0139724
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0138922
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0138263
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0138249
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0136074
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0136015
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0135892
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0135712
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.013525
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0133081
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0132867
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0132392
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0131409
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.013086
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0130649
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0130122
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0126844
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0123184
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0120081
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0114825
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0114772
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0108422
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.010587
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0105624
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0094603
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0094106
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.009336
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0092028
    Base.precompile(Tuple{typeof(coldotplus_turbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0090629
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0090276
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0090157
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0084102
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0083155
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0083082
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Adjoint{Float32, Vector{Float32}}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.008116
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0080825
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0080601
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.007976
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0079624
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0079463
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0079455
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0078216
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Adjoint{Float64, Vector{Float64}}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0078159
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0076952
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0076564
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0076061
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0075549
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0075437
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0075119
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Vector{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0075086
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0074755
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0074641
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},Vector{UnitRange{Int64}}})   # time: 0.0074396
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0074316
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0074202
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.007409
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0073862
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0073694
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0073303
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0072792
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.007266
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0072617
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0072466
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0072237
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0072202
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0072163
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0072135
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0071774
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0071469
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0071424
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0071292
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070991
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0070687
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070633
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.007059
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0070409
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0070393
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0070257
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070195
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0069938
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0069246
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0069145
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0069019
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0068978
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0068856
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0068848
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068781
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0068767
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0068581
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0068563
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0068151
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0068097
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.006772
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0067649
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.006732
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0067282
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0066758
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0066605
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0066495
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0066309
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0066258
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0066043
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0065928
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0065736
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0065368
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0065189
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0062938
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0062527
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0062146
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0061803
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0061635
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0061518
    Base.precompile(Tuple{typeof(panelsum_turbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.006142
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0061206
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0060723
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.006013
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0058719
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0058315
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0058117
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0057976
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0057849
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0055772
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0055215
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0054493
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0054191
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0052787
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0050819
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0050561
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.004997
    Base.precompile(Tuple{var"#150#threadsfor_fun#2"{Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0048654
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.004852
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0048508
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0048143
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0047492
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0047086
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0046448
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0045657
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0045369
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0045266
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0044609
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0044603
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0044171
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0042283
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0041676
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0041123
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0040317
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0039732
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Vector{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0038176
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Vector{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0037039
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0036245
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0035412
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0034894
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0034398
    Base.precompile(Tuple{var"#303#threadsfor_fun#4"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Vector{Float32}, SubArray{Float32, 1, Vector{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0033734
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0032131
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0031914
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0030588
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0029822
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0029424
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0029084
    Base.precompile(Tuple{var"#103#threadsfor_fun#1"{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0026717
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0024415
    Base.precompile(Tuple{var"#257#threadsfor_fun#3"{Vector{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0024048
    Base.precompile(Tuple{Type{BoottestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0015056
    Base.precompile(Tuple{Type{BoottestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0014248
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0011542
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0011381
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0010797
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010331
end
