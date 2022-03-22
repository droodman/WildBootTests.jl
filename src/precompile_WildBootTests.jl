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
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 8.548078
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 7.01583
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 3.9323723
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 3.7068539
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 1.5354811
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float32}, Bool, Float16, Symbol, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float32, Float32, Symbol, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.9318629
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :LIML, :Fuller, :kappa, :ARubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :NH0, :ML, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getCI, :getplot, :getauxweights), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float64}, Bool, Float16, Symbol, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float64, Float64, Symbol, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 0.8815
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.3031049
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0915435
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0544844
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int16}})   # time: 0.0478876
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :hetrobust, :scorebs), Tuple{Vector{Float32}, Matrix{Float64}, Bool, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0475819
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :hetrobust, :scorebs, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0454296
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :reps, :scorebs, :predexog, :hetrobust), Tuple{Vector{Float32}, Int64, Bool, Matrix{Float64}, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0448079
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0446246
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull, :getplot, :ptype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool, Bool, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0425788
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int8}})   # time: 0.0392885
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0354483
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0328798
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0296697
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Bool,Matrix{Int64},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0293988
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Bool,Matrix{Int64},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0285995
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0285214
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Bool,Vector{Int8},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0276205
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Bool,Vector{Int32},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0271597
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Bool,Vector{Int8},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0258724
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Bool,Vector{Int32},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0252568
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.024167
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.023559
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.0226892
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0225655
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0224888
    Base.precompile(Tuple{typeof(vecconvert),DataType,BitVector})   # time: 0.0224703
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int8}})   # time: 0.0223958
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0222066
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.0207188
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.01979
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.0196539
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Symmetric{Float32, Matrix{Float32}},Matrix{Float32}})   # time: 0.0179556
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Symmetric{Float64, Matrix{Float64}},Matrix{Float64}})   # time: 0.0175957
    Base.precompile(Tuple{typeof(vecconvert),Type{Int64},Vector{Int16}})   # time: 0.0165086
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0103308
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0094425
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.00942
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :scorebs, :predexog, :hetrobust), Tuple{Vector{Float32}, Bool, Matrix{Float64}, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0092454
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),SubArray{Float32, 2, Array{Float32, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0089779
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0087563
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0075981
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :getplot, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0070773
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0061267
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0047216
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0041153
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0040239
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0040137
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0038687
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0036602
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0036449
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0034332
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0030278
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0029808
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float32, 3},Array{Float32, 3},Vector{UnitRange{Int64}}})   # time: 0.0029532
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0028374
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float64, 3},Array{Float64, 3},Vector{UnitRange{Int64}}})   # time: 0.002758
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0024467
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0024202
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int64}})   # time: 0.0019062
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.001591
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0015697
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0014787
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0013633
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Float32}})   # time: 0.0013216
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}})   # time: 0.0012732
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0012667
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}})   # time: 0.0011392
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.001121
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0010394
end
