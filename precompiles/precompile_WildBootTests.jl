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
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 20.230019
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 9.026225
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 5.3694577
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 4.3672113
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 2.6926303
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float64}, Bool, Float16, PType, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float64, Float64, MAdjType, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 1.0552992
    Base.precompile(Tuple{var"#381#threadsfor_fun#9"{Array{Float64, 3}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.2966534
    Base.precompile(Tuple{var"#381#threadsfor_fun#9"{Array{Float32, 3}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.2678224
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.1964834
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float32}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0713363
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0689124
    Base.precompile(Tuple{var"#429#threadsfor_fun#12"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0673696
    Base.precompile(Tuple{var"#429#threadsfor_fun#12"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0637404
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0599707
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype, :gridmin, :gridmax, :gridpoints), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, AuxWtType, Bool, PType, Vector{Int64}, Vector{Int64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0585602
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int16}, Int64, Vector{Int16}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0555774
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Int32}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0539143
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0476525
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0462442
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0455072
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0444034
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps, :imposenull), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.043615
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float64}, Matrix{Float64}, Vector{Int32}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0435343
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0435259
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Int64, Vector{Int32}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0433607
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0433165
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0425376
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Vector{Float64}, Vector{Int64}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0425086
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0418811
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0415969
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0409442
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0408027
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0405764
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0395512
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0395165
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0394982
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0393891
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :Fuller, :clustid, :small, :reps, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Int64, Vector{Int8}, Bool, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0391586
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :ARubin, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0390069
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0389842
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :imposenull), Tuple{Bool, BitVector, Matrix{Float64}, Vector{Int64}, Int64, Bool}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0387751
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.03876
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0387156
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:R1, :r1, :resp, :predexog, :clustid), Tuple{Matrix{Int64}, Vector{Float64}, Vector{Float32}, Matrix{Float64}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0385285
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0385213
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Float64}, Bool, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0384435
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :LIML, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Float32}, Bool, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0382108
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0382066
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, BitVector, Matrix{Float64}, Matrix{Int64}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0380349
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :ptype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0377716
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0374253
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype, :getCI), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.037407
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :gridpoints, :reps), Tuple{Bool, BitVector, Matrix{Float64}, Vector{Int64}, Vector{Int64}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0370105
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0369444
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :feid, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int16}, Int64, Int64, Vector{Int16}, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0369271
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :reps, :auxwttype), Tuple{Bool, Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, AuxWtType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.03679
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :obswt, :feid), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int8}, Int64, Int64, Vector{Int8}, Vector{Int8}}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.036546
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.0362783
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0352735
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Vector{Int8}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0324258
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.023919
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.0227054
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.0223477
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int32}})   # time: 0.0204977
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0182973
    Base.precompile(Tuple{var"#549#threadsfor_fun#19"{StrBootTest{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0176686
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float64}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0176633
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :gridmax, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small, :gridmin, :gridpoints), Tuple{Bool, Vector{Float32}, Vector{Int64}, AuxWtType, Bool, Vector{Float32}, PType, Matrix{Float64}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Int64}, Vector{Int64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.017648
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0169326
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0169091
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :ptype, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, PType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0167606
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights, :turbo), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int8, Int8, Bool, Bool, Vector{Int64}, Bool, Vector{Float32}, Bool, Float16, PType, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Int64, Bool, AuxWtType, MersenneTwister, Float32, Float32, MAdjType, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, DistStatType, Bool, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.0164997
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :obswt, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Vector{Float64}, Vector{Int64}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0164338
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0161948
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0155525
    Base.precompile(Tuple{var"#516#threadsfor_fun#17"{StrBootTest{Float32}, Vector{UnitRange{Int64}}, Vector{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0151424
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :LIML, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Float32}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0146932
    Base.precompile(Tuple{var"#549#threadsfor_fun#19"{StrBootTest{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0145894
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0145647
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:turbo, :resp, :auxwttype, :bootstrapc, :predendog, :ptype, :predexog, :inst, :clustid, :reps, :small), Tuple{Bool, Vector{Float32}, AuxWtType, Bool, Vector{Float64}, PType, Matrix{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0144714
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0144089
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :feid, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Vector{Int32}, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0144025
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{BitVector, Vector{Int64}, Int64, Bool, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0143695
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float64}, Matrix{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0142229
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0138877
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0137866
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0137621
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int32}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0137452
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0136181
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :gridpoints, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Vector{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0135404
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0135346
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0135306
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :r1, :R1, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Vector{Float64}, Matrix{Int64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0133244
    Base.precompile(Tuple{var"#516#threadsfor_fun#17"{StrBootTest{Float64}, Vector{UnitRange{Int64}}, Vector{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0120408
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0119949
    Base.precompile(Tuple{var"#478#threadsfor_fun#15"{StrBootTest{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0109363
    Base.precompile(Tuple{var"#445#threadsfor_fun#13"{StrBootTest{Float64}, Vector{UnitRange{Int64}}, Vector{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0106926
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0097998
    Base.precompile(Tuple{var"#445#threadsfor_fun#13"{StrBootTest{Float32}, Vector{UnitRange{Int64}}, Vector{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.009773
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0094469
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0093579
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.009223
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.009115
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0086943
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.008654
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0085076
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0084128
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float64}, Matrix{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0083539
    Base.precompile(Tuple{var"#478#threadsfor_fun#15"{StrBootTest{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0083422
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0083344
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0083173
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float64}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0079561
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0078658
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0078409
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0077455
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0077407
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :turbo, :predendog, :auxwttype, :predexog, :Fuller), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Vector{Float32}, AuxWtType, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0077077
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0077067
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0075557
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int32}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0075095
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0074952
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0074316
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0074177
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0073824
    Base.precompile(Tuple{var"#397#threadsfor_fun#10"{Array{Float64, 3}, Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0073504
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0073434
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float32}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float32}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0071793
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0071698
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0071647
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0071259
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0071057
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Matrix{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0071045
    Base.precompile(Tuple{var"#137#threadsfor_fun#4"{Float32, Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0070903
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070291
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0070267
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0070245
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Int64},Vector{Int64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.007021
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070147
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0070091
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0069025
    Base.precompile(Tuple{var"#137#threadsfor_fun#4"{Float64, Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.006878
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0068506
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0068156
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0068047
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :ARubin, :small, :turbo, :predendog, :predexog), Tuple{Vector{Float32}, Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, Bool, Vector{Float64}, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0068015
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0068011
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
        end
    end   # time: 0.0067902
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int16},Int64,Int64,Bool,Bool,Vector{Int16},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0067564
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0066937
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0066911
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int64},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0066867
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int32},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0066355
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0066234
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{BitVector, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0065958
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Vector{Float64},Vector{Float64},Matrix{Float32},Vector{Float32},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0065859
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float64}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0065603
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (BitVector,Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0064959
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Int64},Vector{Float64},Vector{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float32},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0064879
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0064318
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :imposenull, :turbo, :nbootclustvar, :predexog, :nerrclustvar), Tuple{Vector{Float64}, Matrix{Int64}, Int64, Bool, Bool, Int64, Matrix{Float64}, Int64}},typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0063756
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0063282
    Base.precompile(Tuple{var"#137#threadsfor_fun#4"{Float64, Adjoint{Float64, Vector{Float64}}, Int64, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0063024
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Int8},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Int8},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
        end
    end   # time: 0.0062447
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0062229
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{Int32},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0062076
    Base.precompile(Tuple{var"#397#threadsfor_fun#10"{Array{Float32, 3}, Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0061539
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Float64},Vector{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Vector{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Int64},Int64,Int64,Bool,Bool,Vector{Int8},Bool,Vector{Float64},Bool,Int64,PType,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Int64,Bool,AuxWtType,MersenneTwister,Float64,Float64,MAdjType,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int32},DistStatType,Bool,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64},))
        end
    end   # time: 0.0061374
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.006096
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :reps, :turbo, :getCI, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Int64, Bool, Bool, AuxWtType, Matrix{Float32}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0060497
    Base.precompile(Tuple{var"#137#threadsfor_fun#4"{Float32, Adjoint{Float32, Vector{Float32}}, Int64, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0060165
    Base.precompile(Tuple{var"#413#threadsfor_fun#11"{Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0056843
    Base.precompile(Tuple{var"#413#threadsfor_fun#11"{Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0056447
    Base.precompile(Tuple{var"#89#threadsfor_fun#1"{SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.0053679
    Base.precompile(Tuple{var"#359#threadsfor_fun#8"{Array{Float32, 3}, Array{Float32, 3}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0052728
    Base.precompile(Tuple{var"#262#threadsfor_fun#6"{Matrix{Float64}, Matrix{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0046918
    Base.precompile(Tuple{var"#154#threadsfor_fun#5"{Float32, Matrix{Float32}, Int64, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0045997
    Base.precompile(Tuple{var"#89#threadsfor_fun#1"{SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0044912
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0043941
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0042473
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0041455
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0041254
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0039616
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{Matrix{Float64}, Matrix{Float64}, SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0039184
    Base.precompile(Tuple{var"#359#threadsfor_fun#8"{Array{Float64, 3}, Array{Float64, 3}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}}})   # time: 0.0038579
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{Matrix{Float32}, Matrix{Float32}, SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0038253
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0037284
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0036663
    Base.precompile(Tuple{var"#154#threadsfor_fun#5"{Float64, Matrix{Float64}, Int64, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Int64}, UnitRange{Int64}}})   # time: 0.0035629
    Base.precompile(Tuple{var"#341#threadsfor_fun#7"{Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0033649
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0033231
    Base.precompile(Tuple{var"#89#threadsfor_fun#1"{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Base.OneTo{Int64}}})   # time: 0.003247
    Base.precompile(Tuple{var"#262#threadsfor_fun#6"{Vector{Float32}, Vector{Float32}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.003156
    Base.precompile(Tuple{var"#262#threadsfor_fun#6"{Vector{Float64}, Vector{Float64}, Vector{UnitRange{Int64}}, CartesianIndices{0, Tuple{}}, CartesianIndices{0, Tuple{}}, Base.OneTo{Int64}}})   # time: 0.0031457
    Base.precompile(Tuple{var"#262#threadsfor_fun#6"{Matrix{Float32}, Matrix{Float32}, Vector{UnitRange{Int64}}, Base.OneTo{Int64}, CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, Base.OneTo{Int64}}})   # time: 0.0031167
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0030614
    Base.precompile(Tuple{var"#89#threadsfor_fun#1"{Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Base.OneTo{Int64}}})   # time: 0.0026603
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing,StrBootTest{Float64}})   # time: 0.0022327
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0019095
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :clustid, :nbootclustvar, :nerrclustvar, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Matrix{Int32}, Int64, Int64, Int64}},typeof(wildboottest),Type,Matrix{Float64},Vector{Int64}})   # time: 0.0012871
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing,StrBootTest{Float32}})   # time: 0.0012558
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0010817
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int32}})   # time: 0.0010784
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:turbo, :resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :auxwttype, :bootstrapc, :ptype), Tuple{Bool, Vector{Float32}, Matrix{Float32}, Vector{Float64}, Vector{Float64}, Vector{Int32}, Bool, Int64, AuxWtType, Bool, PType}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0010023
end
