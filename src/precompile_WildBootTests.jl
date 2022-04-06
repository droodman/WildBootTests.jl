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
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float32}})   # time: 8.2192955
    Base.precompile(Tuple{typeof(boottestOLSARubin!),StrBootTest{Float64}})   # time: 6.797035
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float64}})   # time: 4.739786
    Base.precompile(Tuple{typeof(boottestWRE!),StrBootTest{Float32}})   # time: 4.028054
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :clustid, :auxwttype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Int32}, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 1.4766428
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :liml, :fuller, :kappa, :arubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :nH0, :ml, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getci, :getplot, :getauxweights), Tuple{Vector{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float32}, Bool, Float16, Symbol, Bool, Bool, Float32, Float32, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float32, Float32, Symbol, Int16, Bool, Matrix{Float32}, Vector{Float32}, Symmetric{Float32, Matrix{Float32}}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float32},Vector{Float32}})   # time: 0.9131017
    Base.precompile(Tuple{Core.kwftype(typeof(__wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :R1, :r1, :clustid, :nbootclustvar, :nerrclustvar, :issorted, :hetrobust, :nfe, :feid, :fedfadj, :obswt, :fweights, :maxmatsize, :ptype, :bootstrapc, :liml, :fuller, :kappa, :arubin, :small, :clusteradj, :clustermin, :scorebs, :reps, :imposenull, :auxwttype, :rng, :level, :rtol, :madjtype, :nH0, :ml, :scores, :beta, :A, :gridmin, :gridmax, :gridpoints, :diststat, :getci, :getplot, :getauxweights), Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Matrix{Int64}, Int64, Int64, Bool, Bool, Int64, Vector{Int64}, Int64, Vector{Float64}, Bool, Float16, Symbol, Bool, Bool, Float64, Float64, Bool, Bool, Bool, Bool, Bool, Int64, Bool, Symbol, MersenneTwister, Float64, Float64, Symbol, Int16, Bool, Matrix{Float64}, Vector{Float64}, Symmetric{Float64, Matrix{Float64}}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Symbol, Bool, Bool, Bool}},typeof(__wildboottest),Matrix{Float64},Vector{Float64}})   # time: 0.7967473
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :clustid, :auxwttype, :predexog), Tuple{Vector{Float32}, Vector{Int32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.2241874
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0757203
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple,typeof(_wildboottest),DataType,Matrix{Float64},Vector{Int64}})   # time: 0.0512666
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :reps, :scorebs, :predexog, :hetrobust), Tuple{Vector{Float32}, Int64, Bool, Matrix{Float64}, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0399314
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0370134
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :hetrobust, :scorebs), Tuple{Vector{Float32}, Matrix{Float64}, Bool, Bool}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.036235
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :small, :predendog, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Vector{Float32}, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0359046
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :hetrobust, :scorebs, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Bool, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Float64}})   # time: 0.0341064
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0334703
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0327042
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0324644
    Base.precompile(Tuple{Core.kwftype(typeof(wildboottest)),NamedTuple{(:resp, :predexog, :predendog, :inst, :clustid, :small, :reps, :imposenull, :getplot, :ptype), Tuple{Vector{Float32}, Matrix{Float64}, Vector{Float32}, Matrix{Int8}, Vector{Int8}, Bool, Int64, Bool, Bool, Symbol}},typeof(wildboottest),Type,Matrix{Int64},Vector{Int64}})   # time: 0.0324147
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Bool,Matrix{Int64},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0299475
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float64},Vector{Float64},Bool,Vector{Int8},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.026748
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Vector{Float32},Matrix{Int8},Matrix{Float32},Vector{Float32},Bool,Vector{Int8},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64},))
    end
end   # time: 0.0262361
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Vector{Float64},Bool,Vector{Int32},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float64},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float64},Matrix{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0259131
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0253607
    let fbody = try __lookup_kwbody__(which(_wildboottest, (DataType,Matrix{Int64},Vector{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Vector{Float32},Matrix{Float64},Matrix{Float32},Matrix{Float32},Matrix{Float32},Vector{Float32},Bool,Vector{Int32},Int64,Int64,Bool,Int64,Vector{Int8},Int64,Vector{Float32},Bool,Int64,Symbol,Bool,Bool,Int64,Float64,Bool,Bool,Bool,Bool,Bool,Int64,Bool,Symbol,MersenneTwister,Float64,Float64,Symbol,Int64,Bool,Matrix{Float32},Vector{Float32},Matrix{Float32},Vector{Float32},Vector{Float32},Vector{Int64},Symbol,Bool,Bool,Bool,typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64},))
    end
end   # time: 0.0253141
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0238532
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.022691
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0226257
    Base.precompile(Tuple{typeof(vecconvert),DataType,BitVector})   # time: 0.0223217
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true},Matrix{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Vector{UnitRange{Int64}}})   # time: 0.0221149
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int8}})   # time: 0.0216603
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int16}})   # time: 0.0211393
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int32}})   # time: 0.0209911
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int16}})   # time: 0.0208953
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Symmetric{Float64, Matrix{Float64}},Matrix{Float64}})   # time: 0.0207441
    Base.precompile(Tuple{typeof(vecconvert),DataType,Vector{Int8}})   # time: 0.020406
    Base.precompile(Tuple{typeof(matconvert),DataType,Matrix{Int16}})   # time: 0.0191304
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int8}})   # time: 0.0187837
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Symmetric{Float32, Matrix{Float32}},Matrix{Float32}})   # time: 0.0172341
    Base.precompile(Tuple{typeof(vecconvert),Type{Int64},Vector{Int16}})   # time: 0.0142607
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0126324
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float64, Vector{Float64}},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0094184
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Adjoint{Float32, Vector{Float32}},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0093277
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),SubArray{Float32, 2, Array{Float32, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0087927
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :scorebs, :predexog, :hetrobust), Tuple{Vector{Float32}, Bool, Matrix{Float64}, Bool}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Float64}})   # time: 0.0080437
    Base.precompile(Tuple{Core.kwftype(typeof(_wildboottest)),NamedTuple{(:resp, :inst, :clustid, :reps, :imposenull, :getplot, :small, :predendog, :ptype, :predexog), Tuple{Vector{Float32}, Matrix{Int8}, Vector{Int8}, Int64, Bool, Bool, Bool, Vector{Float32}, Symbol, Matrix{Float64}}},typeof(_wildboottest),DataType,Matrix{Int64},Vector{Int64}})   # time: 0.0079362
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float64},Int64,Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0074214
    Base.precompile(Tuple{typeof(colquadformminus_nonturbo!),Matrix{Float32},Int64,Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0071387
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0066268
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0055284
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float64},Matrix{Float64},Vector{Float64}})   # time: 0.0041138
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float32},Matrix{Float32},Vector{UnitRange{Int64}}})   # time: 0.0040009
    Base.precompile(Tuple{typeof(matmulplus_nonturbo!),Vector{Float32},Matrix{Float32},Vector{Float32}})   # time: 0.0038271
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Matrix{Float64},Matrix{Float64},Vector{UnitRange{Int64}}})   # time: 0.0037496
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}, Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0034816
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float64}}, Vector{Float64}}},NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.003184
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float64},Matrix{Float64}})   # time: 0.0031641
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),SubArray{Float32, 2, Matrix{Float32}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}}, false},Matrix{Float32},Matrix{Float32}})   # time: 0.0028257
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0028149
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float64, 3},Array{Float64, 3},Vector{UnitRange{Int64}}})   # time: 0.002783
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,NamedTuple{(:X, :p), Tuple{Tuple{Vector{Float32}, Vector{Float32}}, Vector{Float32}}},NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0027813
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float32},Vector{Float32},Vector{UnitRange{Int64}}})   # time: 0.0027138
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Vector{Float64},Vector{Float64},Vector{UnitRange{Int64}}})   # time: 0.0026729
    Base.precompile(Tuple{typeof(panelsum_nonturbo!),Array{Float32, 3},Array{Float32, 3},Vector{UnitRange{Int64}}})   # time: 0.0026243
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float64},Matrix{Float64},Matrix{Float64}})   # time: 0.0018618
    Base.precompile(Tuple{typeof(coldotplus_nonturbo!),Matrix{Float32},Matrix{Float32},Matrix{Float32}})   # time: 0.0017882
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Int8}})   # time: 0.0016576
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0014966
    Base.precompile(Tuple{typeof(matconvert),DataType,Vector{Float32}})   # time: 0.0013353
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,NamedTuple{(:X, :p), Tuple{Vector{Float32}, Float32}},Matrix{Float32},Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.0013019
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,Nothing,Nothing,Matrix{Float64},Vector{Float64},Matrix{Float64},SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}})   # time: 0.0012976
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},Nothing})   # time: 0.001277
    Base.precompile(Tuple{Type{BootTestResult{Float32}},Float32,String,Float32,Float32,Int64,Int64,Int64,Int64,Float32,Nothing,Nothing,Nothing,Matrix{Float32},Vector{Float32},Matrix{Float32},SubArray{Float32, 2, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, true}})   # time: 0.00116
    Base.precompile(Tuple{Type{BootTestResult{Float64}},Float64,String,Float64,Float64,Int64,Int64,Int64,Int64,Float64,Nothing,NamedTuple{(:X, :p), Tuple{Vector{Float64}, Float64}},Matrix{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Nothing})   # time: 0.0010974
    Base.precompile(Tuple{typeof(matconvert),Type{Int64},Matrix{Int64}})   # time: 0.0010696
end
