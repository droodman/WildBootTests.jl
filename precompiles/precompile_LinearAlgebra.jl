function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},Matrix{Float32},Vararg{Union{UniformScaling, AbstractVecOrMat{T} where T}, N} where N})   # time: 0.1407018
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},Matrix{Float64},Vararg{Union{UniformScaling, AbstractVecOrMat{T} where T}, N} where N})   # time: 0.110474
    Base.precompile(Tuple{typeof(hvcat),Tuple{Int64, Int64},UniformScaling{Bool},Vararg{Union{UniformScaling, AbstractVecOrMat{T} where T}, N} where N})   # time: 0.1083293
    Base.precompile(Tuple{Type{Symmetric},Matrix{Float32}})   # time: 0.0017065
end
