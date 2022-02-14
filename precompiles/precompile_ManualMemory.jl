function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(load_aggregate),Type,Int64})   # time: 0.0142389
    Base.precompile(Tuple{typeof(store_aggregate!),Expr,Symbol,Type{NTuple{8, VecElement{Float32}}},Int64})   # time: 0.0059738
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64}},Int64})   # time: 0.00563
    Base.precompile(Tuple{typeof(which(load,(Ptr{UInt64},Type{T},Int64,)).generator.gen),Any,Any,Any,Any,Any})   # time: 0.0055803
    Base.precompile(Tuple{typeof(store_aggregate!),Expr,Symbol,Type{NTuple{4, VecElement{Float64}}},Int64})   # time: 0.0055523
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64, Int64}},Int64})   # time: 0.0048986
    Base.precompile(Tuple{typeof(store!),Ptr{UInt64},NTuple{8, VecElement{Float32}},Int64})   # time: 0.0045713
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64}},Int64})   # time: 0.0041913
    Base.precompile(Tuple{typeof(load_aggregate),Type{NTuple{4, VecElement{Float64}}},Int64})   # time: 0.0036152
    Base.precompile(Tuple{typeof(load_aggregate),Type{NTuple{8, VecElement{Float32}}},Int64})   # time: 0.0035637
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}},Int64})   # time: 0.0035612
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}},Int64})   # time: 0.0034061
    Base.precompile(Tuple{typeof(which(store!,(Ptr{UInt64},T,Int64,)).generator.gen),Any,Any,Any,Any,Any})   # time: 0.0033045
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64, Int64}},Int64})   # time: 0.0030105
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}},Int64})   # time: 0.0029805
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64}},Int64})   # time: 0.0029741
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64}},Int64})   # time: 0.0029582
    Base.precompile(Tuple{typeof(load_aggregate),Type{Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}},Int64})   # time: 0.0029542
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64, Int64}},Int64})   # time: 0.0018586
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64}},Int64})   # time: 0.0018212
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}},Int64})   # time: 0.0015102
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Int64, Int64, Int64}},Int64})   # time: 0.0014476
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64, Int64}},Int64})   # time: 0.0013898
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{NTuple{8, VecElement{Float32}}},Int64})   # time: 0.00133
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}},Int64})   # time: 0.0012666
    Base.precompile(Tuple{typeof(load),Ptr{UInt64},Type{Tuple{Int64, Int64, Int64, Int64, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64, Int64}},Int64})   # time: 0.0012552
end
