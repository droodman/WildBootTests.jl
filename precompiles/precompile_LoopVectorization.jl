function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(add_loops!),LoopSet,Vector{Any},Core.SimpleVector})   # time: 0.086123
    Base.precompile(Tuple{typeof(_turbo_loopset),Any,Any,Any,Any,Core.SimpleVector,Core.SimpleVector,Tuple{Bool, Int8, Int8, Int8, Bool, Int64, Int64, Int64, Int64, Int64, Int64, Int64, UInt64}})   # time: 0.0786073
    Base.precompile(Tuple{typeof(fill_offset_memop_collection!),LoopSet})   # time: 0.0602998
    Base.precompile(Tuple{typeof(lower_and_split_loops),LoopSet,Int64})   # time: 0.0527324
    Base.precompile(Tuple{typeof(setup_preamble!),LoopSet,UnrollSpecification,Int64})   # time: 0.0423409
    Base.precompile(Tuple{typeof(_choose_num_blocks),UInt64,StaticInt{2},UInt64,StaticInt{8}})   # time: 0.0410271
    Base.precompile(Tuple{typeof(evaluate_cost_tile!),Vector{Float64},Array{Bool, 3},LoopSet,Vector{Symbol},UnrollSymbols,Bool,Vector{Vector{Symbol}},Vector{Bool}})   # time: 0.0401029
    Base.precompile(Tuple{typeof(process_metadata!),LoopSet,Vector{Any},Int64})   # time: 0.0367302
    Base.precompile(Tuple{typeof(ifelselast),typeof(vfmadd_fast),EVLMask{8, UInt8},StaticInt{1},StaticInt{0},Vec{8, Float32},Vec{8, Float32},VecUnroll{1, 8, Float32, Vec{8, Float32}}})   # time: 0.0351193
    Base.precompile(Tuple{typeof(stride_penalty),LoopSet,Vector{Symbol}})   # time: 0.0351045
    Base.precompile(Tuple{typeof(lower_compute!),Expr,Operation,LoopSet,UnrollArgs,Bool})   # time: 0.032106
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64, Int64}})   # time: 0.028323
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Ptr{Float32}, Ptr{Float32}, StaticType{Float32}}})   # time: 0.0257205
    Base.precompile(Tuple{typeof(cost),LoopSet,Operation,Tuple{Symbol, Symbol},Symbol,Int64,Int64})   # time: 0.0255747
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64}})   # time: 0.0253814
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Ptr{Float64}, Ptr{Float64}, StaticType{Float64}}})   # time: 0.0207195
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64, Int64}})   # time: 0.0198846
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Ptr{Float64}, StaticType{Float64}}})   # time: 0.0192131
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}})   # time: 0.0191582
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Ptr{Float32}, StaticType{Float32}}})   # time: 0.0191096
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}})   # time: 0.0189187
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}})   # time: 0.0185157
    Base.precompile(Tuple{typeof(partialmap),typeof(vfnmadd_fast),Vec{8, Float32},StaticInt{2},StaticInt{3},Vec{8, Float32},VecUnroll{1, 8, Float32, Vec{8, Float32}},Vec{8, Float32}})   # time: 0.0181857
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64}})   # time: 0.0180185
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64}})   # time: 0.0179397
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}})   # time: 0.0178446
    Base.precompile(Tuple{typeof(setup_turbo_threads!),Ptr{UInt64},Ptr{Nothing},Tuple{Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64}})   # time: 0.0176595
    Base.precompile(Tuple{typeof(which(partialmap,(F,D,Static.StaticInt{M},Static.StaticInt{S},Vararg{Any, K},)).generator.gen),Any,Any,Any,Any,Any,Any,Any,Any,Type,Any,Any})   # time: 0.0175455
    Base.precompile(Tuple{typeof(add_loop_start_stop_manager!),LoopSet})   # time: 0.0169426
    Base.precompile(Tuple{typeof(tupleranks),NTuple{8, Int64}})   # time: 0.0167104
    Base.precompile(Tuple{typeof(which(_choose_num_blocks,(UInt64,Static.StaticInt{U},Any,Static.StaticInt{NTMAX},)).generator.gen),Any,Any,Any,Any,Type,Any,Any})   # time: 0.016018
    Base.precompile(Tuple{typeof(check_valid_reorder_dims!),LoopSet})   # time: 0.0146581
    Base.precompile(Tuple{typeof(thread_loop_summary!),LoopSet,UnrollArgs,Loop,Bool})   # time: 0.0140986
    Base.precompile(Tuple{typeof(iterate),LoopOrders})   # time: 0.0135674
    Base.precompile(Tuple{typeof(calc_Ureduct!),LoopSet,UnrollSpecification})   # time: 0.0134936
    Base.precompile(Tuple{typeof(avx_loopset!),LoopSet,Vector{Instruction},Vector{OperationStruct},Vector{ArrayRefStruct},Vector{Any},Vector{Any},Core.SimpleVector,Core.SimpleVector})   # time: 0.013463
    Base.precompile(Tuple{typeof(cost_vec_buf),LoopSet})   # time: 0.0129979
    Base.precompile(Tuple{typeof(_lower_load!),Expr,LoopSet,Operation,UnrollArgs,Bool,Vector{Bool}})   # time: 0.012487
    Base.precompile(Tuple{typeof(mem_offset_u),Operation,UnrollArgs,Vector{Bool},Bool,Int64,LoopSet,Bool})   # time: 0.0120651
    Base.precompile(Tuple{Type{Loop},LoopSet,Expr,Symbol,Int64,Int64,Nothing,Int64})   # time: 0.0120031
    Base.precompile(Tuple{typeof(lower_store!),Expr,LoopSet,Operation,UnrollArgs,Bool,Symbol,Vector{Bool}})   # time: 0.0112794
    Base.precompile(Tuple{typeof(add_prefetches!),Expr,LoopSet,Operation,UnrollArgs,Int64})   # time: 0.0110603
    Base.precompile(Tuple{typeof(lower_block),LoopSet,UnrollSpecification,Int64,Bool,Int64})   # time: 0.0107019
    Base.precompile(Tuple{typeof(lower_unrolled_dynamic),LoopSet,UnrollSpecification,Int64,Bool})   # time: 0.0094568
    Base.precompile(Tuple{typeof(lower_load_no_optranslation!),Expr,LoopSet,Operation,UnrollArgs,Bool,Vector{Bool}})   # time: 0.0092513
    Base.precompile(Tuple{typeof(unrolled_curly),Operation,Int64,Loop,Loop,Bool,Int64})   # time: 0.0091925
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64},Tuple{Tuple{Int64, Int64}, Tuple{Int64}, Tuple{Int64}}})   # time: 0.0091543
    Base.precompile(Tuple{typeof(_add_mref!),Expr,LoopSet,ArrayReferenceMeta,Symbol,Int64,Int64,Vector{Int64},Symbol})   # time: 0.0089063
    Base.precompile(Tuple{typeof(valid_thread_loops),LoopSet})   # time: 0.0088883
    Base.precompile(Tuple{typeof(thread_one_loops_expr),LoopSet,UnrollArgs,Vector{Bool},UInt64,Float64,Tuple{Bool, Int8, Int8, Int8, Bool, Int64, Int64, Int64, Int64, Int64, Int64, Int64, UInt64},Expr,Expr,Expr,Expr})   # time: 0.0087365
    Base.precompile(Tuple{typeof(parent_op_name!),Expr,LoopSet,Vector{Operation},Int64,Int64,Int64,Vector{Bool},Vector{Bool},Int64,Int64,Bool,Operation,Int64})   # time: 0.0083052
    Base.precompile(Tuple{Type{LoopSet},Symbol})   # time: 0.0080399
    Base.precompile(Tuple{typeof(evaluate_cost_unroll),LoopSet,Vector{Symbol},Symbol,Float64,Vector{Vector{Symbol}}})   # time: 0.0077822
    Base.precompile(Tuple{typeof(reg_pres_buf),LoopSet})   # time: 0.00742
    Base.precompile(Tuple{typeof(_choose_num_blocks),UInt64,StaticInt{1},UInt64,StaticInt{8}})   # time: 0.0073074
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64},Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Tuple{Int64}}})   # time: 0.0072448
    Base.precompile(Tuple{typeof(lower_licm_constants!),LoopSet})   # time: 0.0071505
    Base.precompile(Tuple{Type{ArrayReferenceMeta},LoopSet,UInt128,UInt128,UInt128,UInt128,Symbol,Symbol,Vector{Symbol},Vector{Symbol},Vector{Int64},Vector{Bool}})   # time: 0.0069964
    Base.precompile(Tuple{typeof(thread_two_loops_expr),LoopSet,UnrollArgs,Vector{Bool},UInt64,Float64,Tuple{Bool, Int8, Int8, Int8, Bool, Int64, Int64, Int64, Int64, Int64, Int64, Int64, UInt64},Expr,Expr,Expr,Expr})   # time: 0.0069828
    Base.precompile(Tuple{typeof(addvectoroffset!),Expr,Bool,Int64,MaybeKnown,MaybeKnown,Int64,Symbol,Int64,Bool,Bool})   # time: 0.0067817
    Base.precompile(Tuple{typeof(pointermax_index),LoopSet,ArrayReferenceMeta,Int64,Int64,Bool,Symbol,MaybeKnown})   # time: 0.0066379
    Base.precompile(Tuple{typeof(choose_tile),LoopSet,Vector{Vector{Symbol}},Int64})   # time: 0.0065461
    Base.precompile(Tuple{typeof(use_loop_induct_var!),LoopSet,Expr,ArrayReferenceMeta,Vector{ArrayReferenceMeta},Int64,Bool})   # time: 0.0064208
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),Tuple{Int64, Int64},Tuple{Int64, Int64},Tuple{Tuple{Int64}, Tuple{Int64}}})   # time: 0.0061578
    Base.precompile(Tuple{typeof(mem_offset),Operation,UnrollArgs,Vector{Bool},Bool,LoopSet,Bool})   # time: 0.0059353
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64},Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Tuple{Int64, Int64}}})   # time: 0.0058442
    Base.precompile(Tuple{typeof(rejectcurly),LoopSet,Operation,Symbol,Symbol})   # time: 0.0058014
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),NTuple{4, Int64},NTuple{4, Int64},Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Tuple{Int64, Int64}, Tuple{Int64}}})   # time: 0.0057038
    Base.precompile(Tuple{typeof(lower_zero!),Expr,Operation,LoopSet,UnrollArgs,NumberType})   # time: 0.0054717
    Base.precompile(Tuple{typeof(ifelselastexpr),Bool,Int64,Tuple{DataType, DataType, DataType},Int64,Int64,Bool})   # time: 0.0054586
    Base.precompile(Tuple{typeof(which(ifelselast,(F,VectorizationBase.AbstractMask{W, U} where U<:Union{UInt128, UInt16, UInt32, UInt64, UInt8},Static.StaticInt{M},Static.StaticInt{S},Vararg{Any, K},)).generator.gen),Any,Any,Any,Any,Any,Any,Any,Any,Type,Any,Any})   # time: 0.0054196
    Base.precompile(Tuple{typeof(stride_penalty),LoopSet,Operation,Vector{Symbol},Vector{Int64}})   # time: 0.0052777
    Base.precompile(Tuple{typeof(ifelselastexpr),Bool,Int64,Tuple{DataType, DataType},Int64,Int64,Bool})   # time: 0.0050588
    Base.precompile(Tuple{typeof(add_op!),LoopSet,Instruction,Vector{OperationStruct},Vector{Int64},Vector{Bool},Int64,Vector{ArrayReferenceMeta},Symbol,Vector{Int64}})   # time: 0.0049776
    Base.precompile(Tuple{typeof(determine_unroll_factor),LoopSet,Vector{Symbol},Symbol,Symbol})   # time: 0.0049769
    Base.precompile(Tuple{typeof(unrolledindex),Operation,UnrollArgs,Bool,Vector{Bool},LoopSet})   # time: 0.0048681
    Base.precompile(Tuple{typeof(add_memory_mask!),Expr,Operation,UnrollArgs,Bool,LoopSet,Int64})   # time: 0.0047874
    Base.precompile(Tuple{typeof(add_constant_offset_load_elmination_cost!),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Base.RefValue{Bool},LoopSet,Operation,Float64,UnrollSymbols,Bool,Bool,Int64,Int64,Bool})   # time: 0.0047594
    Base.precompile(Tuple{typeof(pushgespind!),Expr,LoopSet,Symbol,Int64,Int64,Symbol,Bool,Bool,Bool})   # time: 0.0047195
    Base.precompile(Tuple{typeof(iterate),LoopOrders,Tuple{SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}, SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}}})   # time: 0.0046751
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}}}})   # time: 0.0045884
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{OptionallyStaticUnitRange{StaticInt{0}, Int64}}}})   # time: 0.0045172
    Base.precompile(Tuple{typeof(terminatecondition),LoopSet,UnrollSpecification,Int64,Bool,Int64})   # time: 0.004378
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, StaticInt{4}}}})   # time: 0.0043283
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64, Int64}}})   # time: 0.0042894
    Base.precompile(Tuple{typeof(setunrolled!),LoopSet,Operation,Symbol,Symbol,Symbol})   # time: 0.0041489
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}}}})   # time: 0.0040705
    Base.precompile(Tuple{typeof(expandbyoffset!),Vector{Tuple{Int64, Tuple{Int64, Int32, Bool}}},Vector{Any},Vector{Int64},Bool})   # time: 0.0040366
    Base.precompile(Tuple{typeof(which(of_same_size,(Type{T},Type{S},)).generator.gen),Any,Any,Any,Type,Any})   # time: 0.0039963
    Base.precompile(Tuple{typeof(stabilize_grouped_stridedpointer_type),Tuple{Int64},Tuple{Int64},Tuple{Tuple{Int64}}})   # time: 0.0039683
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}}}})   # time: 0.0039661
    Base.precompile(Tuple{typeof(expandbyoffset!),Vector{Tuple{Int64, Float64}},Vector{Any},Vector{Int64},Bool})   # time: 0.0039578
    Base.precompile(Tuple{typeof(expandbyoffset!),Vector{Tuple{Int64, NumberType}},Vector{Any},Vector{Int64},Bool})   # time: 0.0039466
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64}}})   # time: 0.0039161
    Base.precompile(Tuple{typeof(expandbyoffset!),Vector{Int64},Vector{Any},Vector{Int64},Bool})   # time: 0.0038964
    Base.precompile(Tuple{typeof(_create_mrefs!),LoopSet,Vector{ArrayRefStruct},Vector{Symbol},Vector{Symbol},Vector{Int64},Vector{Bool},Core.SimpleVector,Vector{Int64},Vector{Int64},Vector{Tuple{NTuple{8, Int64}, Int64}}})   # time: 0.0038292
    Base.precompile(Tuple{typeof(isoptranslation),LoopSet,Operation,UnrollSymbols})   # time: 0.0037769
    Base.precompile(Tuple{typeof(incrementloopcounter!),Expr,LoopSet,UnrollSpecification,Int64,Int64})   # time: 0.0037603
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64}}})   # time: 0.0037112
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, StaticInt{8}}}})   # time: 0.0036926
    Base.precompile(Tuple{typeof(maxnegativeoffset),LoopSet,Operation,Symbol})   # time: 0.0036215
    Base.precompile(Tuple{typeof(choose_unroll_order),LoopSet,Float64,Vector{Vector{Symbol}},Int64})   # time: 0.0036
    Base.precompile(Tuple{typeof(lower_no_unroll),LoopSet,UnrollSpecification,Int64,Bool})   # time: 0.0035774
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}}}})   # time: 0.0035308
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64, Int64}}})   # time: 0.0034931
    Base.precompile(Tuple{typeof(getu₁forreduct),LoopSet,Operation,Int64})   # time: 0.003492
    Base.precompile(Tuple{typeof(append_pointer_maxes!),Expr,LoopSet,ArrayReferenceMeta,Int64,Int64,Bool,Symbol,MaybeKnown})   # time: 0.0034506
    Base.precompile(Tuple{typeof(loopdepindices),LoopSet,Operation})   # time: 0.0034144
    Base.precompile(Tuple{typeof(prefetchisagoodidea),LoopSet,Operation,UnrollArgs})   # time: 0.003407
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64, Int64}}})   # time: 0.0034013
    Base.precompile(Tuple{typeof(initialize_outer_reductions!),Expr,LoopSet,Operation,Int64,UnrollSpecification,StaticInt{32}})   # time: 0.0033111
    Base.precompile(Tuple{typeof(which(_turbo_!,(Val{var"#UNROLL#"},Val{var"#OPS#"},Val{var"#ARF#"},Val{var"#AM#"},Val{var"#LPSYM#"},Val{Tuple{var"#LB#", var"#V#"}},Vararg{Any, var"#num#vargs#"},)).generator.gen),Any,Any,Any,Any,Any,Any,Any,Any,Any,Type,Type,Type,Type,Type,Any,Any})   # time: 0.003311
    Base.precompile(Tuple{typeof(determine_unroll_factor),LoopSet,Vector{Symbol},Symbol,Int64})   # time: 0.0031851
    Base.precompile(Tuple{typeof(choose_order_cost),LoopSet,Int64})   # time: 0.0031416
    Base.precompile(Tuple{typeof(add_upper_outer_reductions),LoopSet,Expr,Int64,Int64,Loop,Bool})   # time: 0.0031169
    Base.precompile(Tuple{typeof(outer_reduct_combine_expressions),LoopSet,Symbol})   # time: 0.002968
    Base.precompile(Tuple{typeof(solve_unroll_lagrange),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Float64,Float64,Int64,Int64,Bool})   # time: 0.0029427
    Base.precompile(Tuple{typeof(solve_unroll),Symbol,Symbol,SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Symbol,Loop,Loop,Int64,Int64,Bool})   # time: 0.0028719
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, StaticInt{4}}}})   # time: 0.0028025
    Base.precompile(Tuple{typeof(unitstride),LoopSet,Operation,Symbol})   # time: 0.0027631
    Base.precompile(Tuple{typeof(avx_body),LoopSet,Tuple{Bool, Int8, Int8, Int8, Bool, Int64, Int64, Int64, Int64, Int64, Int64, Int64, UInt64}})   # time: 0.0026968
    Base.precompile(Tuple{typeof(fillorder!),LoopSet,Vector{Symbol},Symbol,Symbol,Int64,Symbol})   # time: 0.002668
    Base.precompile(Tuple{typeof(isouterreduction),LoopSet,Operation})   # time: 0.0026568
    Base.precompile(Tuple{typeof(set_upstream_family!),Vector{Bool},Operation,Bool,Vector{Symbol},Int64})   # time: 0.0026366
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}}}})   # time: 0.0026313
    Base.precompile(Tuple{typeof(reduce_expr!),Expr,LoopSet,Int64})   # time: 0.0026012
    Base.precompile(Tuple{typeof(mismatchedstorereductions),LoopSet})   # time: 0.0025742
    Base.precompile(Tuple{typeof(load_short_static_reduction_first!),LoopSet,Symbol,Symbol,Symbol})   # time: 0.002554
    Base.precompile(Tuple{typeof(indices_calculated_by_pointer_offsets),LoopSet,ArrayReferenceMeta})   # time: 0.0025245
    Base.precompile(Tuple{typeof(lower_load!),Expr,Operation,LoopSet,UnrollArgs,Bool})   # time: 0.0025039
    Base.precompile(Tuple{typeof(definemask),Loop})   # time: 0.0024782
    Base.precompile(Tuple{typeof(loopvarremcomparison),Loop,Int64,Bool,Bool})   # time: 0.0024755
    Base.precompile(Tuple{typeof(outer_reduction_zero),Operation,Bool,Int64,Float64,StaticInt{32}})   # time: 0.0024187
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}, Int64}}})   # time: 0.0023735
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{NTuple{4, Ptr{Float32}}}})   # time: 0.0023521
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, StaticInt{8}}}})   # time: 0.0023513
    Base.precompile(Tuple{typeof(offset_ptr),ArrayReferenceMeta,UnrollSpecification,Symbol,Int64,Int64,Vector{Bool},Loop})   # time: 0.0023281
    Base.precompile(Tuple{typeof(extract_external_functions!),LoopSet,Int64,Core.SimpleVector})   # time: 0.0023237
    Base.precompile(Tuple{typeof(add_parents_to_ops!),LoopSet,Vector{OperationStruct},Int64})   # time: 0.0023027
    Base.precompile(Tuple{typeof(advance_state!),SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true}})   # time: 0.0022406
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64}}})   # time: 0.0022014
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}}}})   # time: 0.0021653
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}, Int64, Int64}}})   # time: 0.0021538
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{4}, Int64, StaticInt{4}, Int64, StaticInt{4}}}})   # time: 0.0021521
    Base.precompile(Tuple{typeof(partialmap),typeof(vfmadd_fast),VecUnroll{1, 8, Float32, Vec{8, Float32}},StaticInt{2},StaticInt{0},VecUnroll{1, 8, Float32, Vec{8, Float32}},VecUnroll{1, 8, Float32, Vec{8, Float32}},VecUnroll{1, 8, Float32, Vec{8, Float32}}})   # time: 0.0020751
    Base.precompile(Tuple{typeof(push_loop_length_expr!),Expr,LoopSet})   # time: 0.0020709
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{StaticInt{8}, Int64, StaticInt{8}, Int64, StaticInt{8}}}})   # time: 0.0020457
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{NTuple{4, Ptr{Float32}}}})   # time: 0.0020417
    Base.precompile(Tuple{typeof(offsetindex),Int64,Int64,Int64,Bool,MaybeKnown})   # time: 0.0020255
    Base.precompile(Tuple{typeof(determine_eltype),LoopSet,Bool})   # time: 0.0019972
    Base.precompile(Tuple{typeof(startloop),LoopSet,UnrollSpecification,Int64,Bool})   # time: 0.0019783
    Base.precompile(Tuple{typeof(holdopinregister),LoopSet})   # time: 0.0019548
    Base.precompile(Tuple{typeof(rank_to_sortperm),Tuple{NTuple{8, Int64}, Int64}})   # time: 0.0018683
    Base.precompile(Tuple{typeof(determine_unroll_factor),LoopSet,Vector{Symbol},Symbol})   # time: 0.0018665
    Base.precompile(Tuple{typeof(reduce_expr!),Expr,Symbol,Operation,Int64,Int64,Bool,Bool})   # time: 0.0018258
    Base.precompile(Tuple{typeof(load_constrained),Operation,Symbol,Symbol,Symbol,Bool})   # time: 0.0018232
    Base.precompile(Tuple{typeof(loopindex),LoopSet,UInt128,UInt8})   # time: 0.0017983
    Base.precompile(Tuple{typeof(solve_unroll_iter),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Float64,Float64,StepRange{Int64, Int64},StepRange{Int64, Int64}})   # time: 0.0017763
    Base.precompile(Tuple{typeof(partialmap),typeof(vfmadd_fast),VecUnroll{1, 4, Float64, Vec{4, Float64}},StaticInt{2},StaticInt{0},VecUnroll{1, 4, Float64, Vec{4, Float64}},VecUnroll{1, 4, Float64, Vec{4, Float64}},VecUnroll{1, 4, Float64, Vec{4, Float64}}})   # time: 0.0017736
    Base.precompile(Tuple{typeof(resize!),LoopOrder,Int64})   # time: 0.0017671
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{NTuple{4, Ptr{Float64}}}})   # time: 0.0017404
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{Ptr{Float64}}}})   # time: 0.0017383
    Base.precompile(Tuple{typeof(load_elimination_cost_factor!),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Base.RefValue{Bool},LoopSet,Operation,Float64,UnrollSymbols,Int64,Int64})   # time: 0.0017375
    Base.precompile(Tuple{typeof(isunrolled_sym),Operation,Symbol,Symbol,Symbol,Tuple{Bool, Bool}})   # time: 0.0017233
    Base.precompile(Tuple{typeof(add_upper_comp_check),Loop,Expr})   # time: 0.0017202
    Base.precompile(Tuple{typeof(child_cost_untill_vectorized),Operation})   # time: 0.001709
    Base.precompile(Tuple{typeof(outer_reduct_loopordersplit),LoopSet})   # time: 0.001691
    Base.precompile(Tuple{typeof(solve_unroll),LoopSet,Symbol,Symbol,SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Symbol,Int64})   # time: 0.001685
    Base.precompile(Tuple{typeof(uniquearrayrefs),LoopSet})   # time: 0.0016764
    Base.precompile(Tuple{typeof(add_parents_to_op!),LoopSet,Operation,UInt128,UInt128,Int64,Int64})   # time: 0.0016716
    Base.precompile(Tuple{typeof(parent_unroll_status),Operation,Symbol,UnrollSpecification})   # time: 0.0016556
    Base.precompile(Tuple{typeof(hoist_constant_memory_accesses!),LoopSet})   # time: 0.0016322
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{NTuple{4, Ptr{Float64}}}})   # time: 0.0015907
    Base.precompile(Tuple{typeof(addoptoorder!),LoopSet,Vector{Bool},Vector{Bool},Operation,Symbol,Int64,Symbol,Symbol,Symbol,Int64})   # time: 0.0015862
    Base.precompile(Tuple{typeof(mulexpr),Symbol,StaticInt{6}})   # time: 0.0015837
    Base.precompile(Tuple{typeof(reject_candidate),Operation,Symbol,Symbol})   # time: 0.0015819
    Base.precompile(Tuple{typeof(extract_outerreduct_types!),LoopSet,Int64,Core.SimpleVector})   # time: 0.0015672
    Base.precompile(Tuple{typeof(mulexpr),Symbol,StaticInt{5}})   # time: 0.0014888
    Base.precompile(Tuple{typeof(storeinstr_preprend),Operation,Symbol})   # time: 0.0014679
    Base.precompile(Tuple{typeof(gc_preserve),LoopSet,Expr})   # time: 0.0014603
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{Ptr{Float64}}}})   # time: 0.0014358
    Base.precompile(Tuple{typeof(depchain_cost!),LoopSet,Vector{Bool},Operation,Symbol,Symbol,Int64,Int64,Float64,Int64})   # time: 0.0014336
    Base.precompile(Tuple{typeof(mulexpr),Symbol,StaticInt{3}})   # time: 0.0014265
    Base.precompile(Tuple{typeof(parents),LoopSet,UInt128,UInt128})   # time: 0.0014204
    Base.precompile(Tuple{typeof(init_loop_map!),LoopSet})   # time: 0.0014033
    Base.precompile(Tuple{typeof(ifelselast),typeof(vfmadd_fast),EVLMask{4, UInt8},StaticInt{1},StaticInt{0},Vec{4, Float64},Vec{4, Float64},VecUnroll{1, 4, Float64, Vec{4, Float64}}})   # time: 0.0013832
    Base.precompile(Tuple{typeof(make_partial_map!),Expr,Symbol,Int64,Int64})   # time: 0.0013452
    Base.precompile(Tuple{typeof(init_remblock),Loop,LoopStartStopManager,Int64})   # time: 0.0013351
    Base.precompile(Tuple{typeof(typeof_sym),LoopSet,Operation,NumberType})   # time: 0.0013331
    Base.precompile(Tuple{typeof(store_load_deps!),Vector{Symbol},Operation,ArrayReferenceMeta})   # time: 0.0013329
    Base.precompile(Tuple{typeof(solve_unroll),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64,Int64,Float64,Float64,Int64,Int64,Bool})   # time: 0.0013136
    Base.precompile(Tuple{typeof(avx_threads_expr),LoopSet,Tuple{Bool, Int8, Int8, Int8, Bool, Int64, Int64, Int64, Int64, Int64, Int64, Int64, UInt64},UInt64,Expr,Expr,Expr,Expr})   # time: 0.0013109
    Base.precompile(Tuple{typeof(looprange),Int64,Symbol,Symbol})   # time: 0.0013093
    Base.precompile(Tuple{typeof(reduce_parent!),Expr,LoopSet,Operation,Operation,Symbol})   # time: 0.0012891
    Base.precompile(Tuple{Type{UnrollArgs},LoopSet,Int64,UnrollSymbols,Int64,Int64})   # time: 0.0012795
    Base.precompile(Tuple{typeof(which(reassemble_tuple,(Type{T},Tuple,)).generator.gen),Any,Any,Any,Any})   # time: 0.0012734
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{Ptr{Float64}, Ptr{Float64}}}})   # time: 0.0012617
    Base.precompile(Tuple{typeof(determine_width),LoopSet,Symbol})   # time: 0.0012531
    Base.precompile(Tuple{typeof(ifelselast),typeof(add_fast),EVLMask{8, UInt8},StaticInt{1},StaticInt{0},Vec{8, Float32},VecUnroll{1, 8, Float32, Vec{8, Float32}}})   # time: 0.0012241
    Base.precompile(Tuple{typeof(add_ops!),LoopSet,Vector{Instruction},Vector{OperationStruct},Vector{ArrayReferenceMeta},Vector{Int64},Vector{Symbol},Int64,Vector{Int64},Vector{Bool}})   # time: 0.0012235
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{Ptr{Float32}}}})   # time: 0.0012166
    Base.precompile(Tuple{typeof(parents_symvec),LoopSet,UInt128,Bool,Int64})   # time: 0.0012159
    Base.precompile(Tuple{typeof(solve_unroll_constU),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Int64})   # time: 0.0011979
    Base.precompile(Tuple{typeof(loopordersplit),LoopSet})   # time: 0.0011689
    Base.precompile(Tuple{typeof(otherindexunrolled),LoopSet,Symbol,ArrayReferenceMeta})   # time: 0.001151
    Base.precompile(Tuple{typeof(lower),LoopSet,Vector{Symbol},Symbol,Symbol,Symbol,Int64,Int64,Bool})   # time: 0.0011488
    Base.precompile(Tuple{typeof(getroots!),Vector{Bool},LoopSet})   # time: 0.0011439
    Base.precompile(Tuple{typeof(arithmeticexpr),typeof(*),Int64,Symbol,Symbol,Int64,MaybeKnown})   # time: 0.0011411
    Base.precompile(Tuple{typeof(ifelselast),typeof(add_fast),EVLMask{4, UInt8},StaticInt{1},StaticInt{0},Vec{4, Float64},VecUnroll{1, 4, Float64, Vec{4, Float64}}})   # time: 0.0011314
    Base.precompile(Tuple{typeof(loopset_return_value),LoopSet,Val{false}})   # time: 0.0011087
    Base.precompile(Tuple{typeof(loopset_return_value),LoopSet,Val{true}})   # time: 0.0011073
    Base.precompile(Tuple{typeof(fill_children!),LoopSet})   # time: 0.001084
    Base.precompile(Tuple{typeof(loop_boundary!),Expr,LoopSet,Loop,Bool})   # time: 0.0010813
    Base.precompile(Tuple{typeof(_append_fields!),Expr,Expr,Symbol,Type{Tuple{Ptr{Float32}}}})   # time: 0.0010787
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{Ptr{Float32}, Ptr{Float32}, Ptr{Float32}}}})   # time: 0.0010673
    Base.precompile(Tuple{typeof(rebuild_fields),Int64,Type{Tuple{Ptr{Float32}, Ptr{Float32}}}})   # time: 0.0010581
    Base.precompile(Tuple{typeof(pointerremcomparison),LoopSet,Int64,Int64,Int64,Bool,Bool,Loop})   # time: 0.001054
    Base.precompile(Tuple{typeof(accept_reorder_according_to_tracked_reductions),LoopSet,Symbol})   # time: 0.001044
    Base.precompile(Tuple{typeof(lower_unrollspec),LoopSet})   # time: 0.0010273
    Base.precompile(Tuple{typeof(mulexpr),Symbol,StaticInt{4}})   # time: 0.0010255
    Base.precompile(Tuple{typeof(parent_unroll_status),Operation,Symbol,Symbol,Symbol,Int64,UnrollSpecification})   # time: 0.0010172
    Base.precompile(Tuple{typeof(partialmap),typeof(add_fast),VecUnroll{1, 8, Float32, Vec{8, Float32}},StaticInt{2},StaticInt{0},VecUnroll{1, 8, Float32, Vec{8, Float32}},VecUnroll{1, 8, Float32, Vec{8, Float32}}})   # time: 0.0010077
    Base.precompile(Tuple{typeof(update_cost_vec!),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Float64,Bool,Bool})   # time: 0.0010058
end
