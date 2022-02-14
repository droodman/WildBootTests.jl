function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{(:year, :selfemployed, :hasinsurance, :post, :post_self), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float32}, Vector{Float32}, Vector{Float32}}},Dict{Symbol, UnionAll}})   # time: 0.4119124
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, NTuple{22, Term}}})   # time: 0.0910614
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry, :union), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0809501
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term}}})   # time: 0.0702889
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float32}, ContinuousTerm{Float32}}}},NamedTuple{(:year, :selfemployed, :hasinsurance, :post, :post_self), Tuple{Vector{Float32}, Vector{Int8}, Vector{Float32}, Vector{Float32}, Vector{Float32}}}})   # time: 0.0608134
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float32}, ContinuousTerm{Float32}}})   # time: 0.0503639
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term}},Schema})   # time: 0.0405108
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{names, T} where {N, D, names, T<:Tuple{Vararg{AbstractArray{S, D} where S, N}}},Dict{Symbol, UnionAll}})   # time: 0.0353225
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0334001
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, NTuple{22, Term}},Schema})   # time: 0.0331187
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{NTuple{22, ContinuousTerm{Float32}}}},NamedTuple{(:lnkm, :pixpetro, :pixdia, :pixwaterd, :pixcapdist, :pixmal, :pixsead, :pixsuit, :pixelev, :pixbdist, :lnwaterkm, :lnkm2split, :mean_elev, :mean_suit, :malariasuit, :petroleum, :diamondd, :capdistance1, :seadist1, :borderdist1, :lnl0708s, :centr_tribe, :lnpd0, :pixwbcode, :pixcluster, :ccode, :pixcode), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{String}, Vector{String}, Vector{Int16}, Vector{Int16}}}})   # time: 0.0329905
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:lnkm, :pixpetro, :pixdia, :pixwaterd, :pixcapdist, :pixmal, :pixsead, :pixsuit, :pixelev, :pixbdist, :lnwaterkm, :lnkm2split, :mean_elev, :mean_suit, :malariasuit, :petroleum, :diamondd, :capdistance1, :seadist1, :borderdist1, :lnl0708s, :centr_tribe, :lnpd0, :pixwbcode, :pixcluster, :ccode, :pixcode), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{String}, Vector{String}, Vector{Int16}, Vector{Int16}}},Dict{Symbol, Any}})   # time: 0.0327334
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry, :union), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0320158
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:wage, :ttl_exp, :collgrad, :tenure, :age, :industry, :occupation, :hours), Tuple{Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0317444
    Base.precompile(Tuple{typeof(concrete_term),Term,BitVector,Nothing})   # time: 0.0290505
    Base.precompile(Tuple{typeof(getindex),ContrastsMatrix{DummyCoding, Int16, U} where U,BitVector,Colon})   # time: 0.0274611
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple,Dict{Symbol, UnionAll}})   # time: 0.0248716
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float32}}})   # time: 0.0233771
    Base.precompile(Tuple{typeof(collect_matrix_terms),NTuple{22, ContinuousTerm{Float32}}})   # time: 0.0212824
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0212711
    Base.precompile(Tuple{typeof(getindex),ContrastsMatrix{DummyCoding, Int8, U} where U,BitVector,Colon})   # time: 0.0212044
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}}})   # time: 0.0192244
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term, Term, Term}},Schema})   # time: 0.0182287
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float32}, ContinuousTerm{Float32}}}},NamedTuple})   # time: 0.0178658
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float64}}})   # time: 0.0172021
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, CategoricalTerm{DummyCoding, Int16, 11}, CategoricalTerm{DummyCoding, Int8, 41}}})   # time: 0.0170573
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}})   # time: 0.0168416
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0166512
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0157777
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term}},Schema})   # time: 0.0157046
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry, :union), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0146824
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float32}}}},NamedTuple})   # time: 0.0144165
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0143037
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}}}},NamedTuple})   # time: 0.0139188
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float64}, ContinuousTerm{Float32}}})   # time: 0.0134713
    Base.precompile(Tuple{typeof(fuzzymatch),NTuple{7, Symbol},Symbol})   # time: 0.013466
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float32}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0133498
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}}}},NamedTuple})   # time: 0.0132285
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term}},Schema})   # time: 0.0130178
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}}}},NamedTuple})   # time: 0.0126271
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float64}}}},NamedTuple})   # time: 0.0124935
    Base.precompile(Tuple{typeof(collect_matrix_terms),Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}})   # time: 0.012483
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{Term, Term, Term}},Schema})   # time: 0.012461
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{Term, Term}},Schema})   # time: 0.012455
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float64}, ContinuousTerm{Float32}}}},NamedTuple})   # time: 0.0124502
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{NTuple{22, ContinuousTerm{Float32}}}},NamedTuple})   # time: 0.0123608
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float64}, ContinuousTerm{Float32}}}},NamedTuple{(:wage, :ttl_exp, :collgrad, :tenure, :age, :industry, :occupation, :hours), Tuple{Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0122418
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}}},NamedTuple})   # time: 0.0122003
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float64}}}},NamedTuple})   # time: 0.0121144
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float64}, MatrixTerm{Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, CategoricalTerm{DummyCoding, Int16, 11}, CategoricalTerm{DummyCoding, Int8, 41}}}},NamedTuple})   # time: 0.0120546
    Base.precompile(Tuple{Type{ContrastsMatrix},DummyCoding,Vector{Int16}})   # time: 0.0115995
    Base.precompile(Tuple{Type{ContrastsMatrix},DummyCoding,Vector{Int8}})   # time: 0.0111271
    Base.precompile(Tuple{typeof(modelcols),FormulaTerm{ContinuousTerm{Float32}, MatrixTerm{Tuple{ContinuousTerm{Float64}, ContinuousTerm{Float64}}}},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}}})   # time: 0.0105979
    Base.precompile(Tuple{Type{ContrastsMatrix},Matrix{Float64},Vector{Int16},Vector{Int16},DummyCoding})   # time: 0.0094744
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term}},Schema,Type})   # time: 0.0093121
    Base.precompile(Tuple{Type{ContrastsMatrix},Matrix{Float64},Vector{Int8},Vector{Int8},DummyCoding})   # time: 0.0085177
    Base.precompile(Tuple{typeof(contrasts_matrix),DummyCoding,Int64,Int64})   # time: 0.0072118
    Base.precompile(Tuple{typeof(termnames),DummyCoding,Vector{Int16},Int64})   # time: 0.0057629
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{Term, Term}},Schema,Type})   # time: 0.0056834
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, NTuple{22, Term}},Schema,Type})   # time: 0.0056355
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{Term, Term}}})   # time: 0.0052524
    Base.precompile(Tuple{typeof(termnames),DummyCoding,Vector{Int8},Int64})   # time: 0.0052279
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{names, T} where {N, D, names, T<:Tuple{Vararg{AbstractArray{S, D} where S, N}}},Dict{Symbol, Any}})   # time: 0.005111
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term}},Schema,Type})   # time: 0.0049698
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term, Term, Term}}})   # time: 0.0049159
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term, Term, Term}},Schema,Type})   # time: 0.0047003
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term, Term, Term}}})   # time: 0.0046312
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term}},Schema,Type})   # time: 0.0044502
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Tuple{Term, Term, Term}},Schema,Type})   # time: 0.004168
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{Term, Term, Term}}})   # time: 0.0041097
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0040841
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Tuple{ConstantTerm{Int64}, Term, Term}}})   # time: 0.0039307
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :smsa, :race, :age, :union, :married, :industry), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0037805
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry, :union), Tuple{Vector{Float32}, Vector{Float32}, Vector{Float32}, Vector{Int8}, Vector{Int8}, Vector{Int8}}},Dict{Symbol, Any}})   # time: 0.0036193
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, ConstantTerm{Int64}},Schema})   # time: 0.0035027
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, Term}})   # time: 0.0029041
    Base.precompile(Tuple{typeof(apply_schema),FormulaTerm{Term, Term},Schema})   # time: 0.0027635
    Base.precompile(Tuple{typeof(terms),FormulaTerm{Term, ConstantTerm{Int64}}})   # time: 0.0025727
    Base.precompile(Tuple{typeof(schema),Vector{AbstractTerm},NamedTuple{names, T} where {N, D, names, T<:Tuple{Vararg{AbstractArray{S, D} where S, N}}},Dict{Symbol, Any}})   # time: 0.0024811
    Base.precompile(Tuple{typeof(+),Tuple{ContinuousTerm{Float32}, Vararg{ContinuousTerm{Float32}, N} where N},ContinuousTerm{Float32}})   # time: 0.002198
    Base.precompile(Tuple{typeof(apply_schema),ConstantTerm{Int64},Schema,Type})   # time: 0.00189
    Base.precompile(Tuple{typeof(+),NTuple{20, Term},Term})   # time: 0.0018271
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, CategoricalTerm{DummyCoding, Int16, 11}},CategoricalTerm{DummyCoding, Int8, 41}})   # time: 0.0017845
    Base.precompile(Tuple{typeof(+),Tuple{ConstantTerm{Int64}, Term, Term, Term, Term, Term},Term})   # time: 0.0017826
    Base.precompile(Tuple{typeof(+),NTuple{6, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0017443
    Base.precompile(Tuple{typeof(+),NTuple{11, Term},Term})   # time: 0.0016931
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}},ContinuousTerm{Float64}})   # time: 0.0016676
    Base.precompile(Tuple{typeof(+),NTuple{5, Term},Term})   # time: 0.0016642
    Base.precompile(Tuple{typeof(+),NTuple{7, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.001651
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}},ContinuousTerm{Float64}})   # time: 0.0016416
    Base.precompile(Tuple{typeof(+),NTuple{12, Term},Term})   # time: 0.0016288
    Base.precompile(Tuple{typeof(+),NTuple{21, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0016229
    Base.precompile(Tuple{typeof(+),NTuple{14, Term},Term})   # time: 0.001615
    Base.precompile(Tuple{typeof(+),NTuple{18, Term},Term})   # time: 0.0016114
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float32}, ContinuousTerm{Float32}},ContinuousTerm{Float64}})   # time: 0.0015875
    Base.precompile(Tuple{typeof(+),NTuple{21, Term},Term})   # time: 0.0015869
    Base.precompile(Tuple{typeof(+),NTuple{18, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015836
    Base.precompile(Tuple{typeof(+),Tuple{ConstantTerm{Int64}, Term, Term, Term, Term},Term})   # time: 0.0015823
    Base.precompile(Tuple{typeof(+),NTuple{16, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015704
    Base.precompile(Tuple{typeof(+),NTuple{9, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015592
    Base.precompile(Tuple{typeof(+),NTuple{6, Term},Term})   # time: 0.0015573
    Base.precompile(Tuple{typeof(+),NTuple{19, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015572
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}},ContinuousTerm{Float32}})   # time: 0.001554
    Base.precompile(Tuple{typeof(+),NTuple{8, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015528
    isdefined(StatsModels, Symbol("#34#35")) && Base.precompile(Tuple{getfield(StatsModels, Symbol("#34#35")),ContinuousTerm{Float64}})   # time: 0.0015496
    Base.precompile(Tuple{typeof(+),NTuple{17, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015447
    Base.precompile(Tuple{typeof(+),Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float32}, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015433
    Base.precompile(Tuple{typeof(+),NTuple{12, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015409
    Base.precompile(Tuple{typeof(+),NTuple{10, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015407
    Base.precompile(Tuple{typeof(+),NTuple{17, Term},Term})   # time: 0.0015404
    Base.precompile(Tuple{typeof(+),NTuple{5, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015404
    Base.precompile(Tuple{typeof(+),NTuple{20, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015377
    Base.precompile(Tuple{typeof(+),NTuple{11, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015345
    Base.precompile(Tuple{typeof(+),NTuple{14, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015298
    Base.precompile(Tuple{typeof(+),NTuple{19, Term},Term})   # time: 0.0015243
    Base.precompile(Tuple{typeof(+),NTuple{13, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015234
    Base.precompile(Tuple{typeof(+),NTuple{15, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.001516
    Base.precompile(Tuple{typeof(+),NTuple{4, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0015039
    Base.precompile(Tuple{typeof(+),NTuple{16, Term},Term})   # time: 0.0014947
    Base.precompile(Tuple{typeof(+),NTuple{15, Term},Term})   # time: 0.0014715
    Base.precompile(Tuple{typeof(+),NTuple{13, Term},Term})   # time: 0.0014694
    Base.precompile(Tuple{typeof(+),NTuple{9, Term},Term})   # time: 0.0014666
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.001461
    Base.precompile(Tuple{typeof(+),NTuple{7, Term},Term})   # time: 0.0014602
    Base.precompile(Tuple{typeof(+),NTuple{10, Term},Term})   # time: 0.0014597
    Base.precompile(Tuple{typeof(+),NTuple{8, Term},Term})   # time: 0.0014481
    Base.precompile(Tuple{Type{CategoricalTerm},Symbol,ContrastsMatrix{DummyCoding, Int16, Int16}})   # time: 0.0014392
    Base.precompile(Tuple{typeof(+),Tuple{Term, Term, Term},Term})   # time: 0.0014206
    Base.precompile(Tuple{Type{CategoricalTerm},Symbol,ContrastsMatrix{DummyCoding, Int8, Int8}})   # time: 0.0014149
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0013906
    Base.precompile(Tuple{typeof(+),Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float64}},ContinuousTerm{Float32}})   # time: 0.0013795
    Base.precompile(Tuple{typeof(+),NTuple{4, Term},Term})   # time: 0.0013752
    Base.precompile(Tuple{typeof(+),Tuple{ContinuousTerm{Float32}, ContinuousTerm{Float32}},ContinuousTerm{Float32}})   # time: 0.0012995
    isdefined(StatsModels, Symbol("#34#35")) && Base.precompile(Tuple{getfield(StatsModels, Symbol("#34#35")),ContinuousTerm{Float64}})   # time: 0.0011803
    Base.precompile(Tuple{typeof(+),Tuple{InterceptTerm{true}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}, ContinuousTerm{Float64}},CategoricalTerm{DummyCoding, Int16, 11}})   # time: 0.0010619
    Base.precompile(Tuple{typeof(+),Tuple{Union{Tuple{AbstractTerm, Vararg{AbstractTerm, N} where N}, AbstractTerm}, Vararg{Union{Tuple{AbstractTerm, Vararg{AbstractTerm, N} where N}, AbstractTerm}, N} where N},CategoricalTerm{DummyCoding, Int8, 41}})   # time: 0.001006
end
