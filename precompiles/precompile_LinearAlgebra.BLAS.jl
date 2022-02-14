function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(gemm!),Char,Char,Float64,Matrix{Float64},Vector{Float64},Float64,Vector{Float64}})   # time: 0.004215
    Base.precompile(Tuple{typeof(gemm!),Char,Char,Float32,Matrix{Float32},Vector{Float32},Float32,Vector{Float32}})   # time: 0.0041723
    Base.precompile(Tuple{typeof(symv!),Char,Float32,Matrix{Float32},Vector{Float32},Float32,Vector{Float32}})   # time: 0.0040382
    Base.precompile(Tuple{typeof(symv!),Char,Float64,Matrix{Float64},Vector{Float64},Float64,Vector{Float64}})   # time: 0.003765
    Base.precompile(Tuple{typeof(dot),Vector{Float64},SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0032586
    Base.precompile(Tuple{typeof(symm!),Char,Char,Float32,Matrix{Float32},Matrix{Float32},Float32,Matrix{Float32}})   # time: 0.0032472
    Base.precompile(Tuple{typeof(symm!),Char,Char,Float64,Matrix{Float64},Matrix{Float64},Float64,Matrix{Float64}})   # time: 0.0032373
    Base.precompile(Tuple{typeof(dot),Vector{Float32},SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0030391
    Base.precompile(Tuple{typeof(dot),Vector{Float64},Vector{Float64}})   # time: 0.002641
    Base.precompile(Tuple{typeof(dot),Vector{Float32},Vector{Float32}})   # time: 0.0022609
end
