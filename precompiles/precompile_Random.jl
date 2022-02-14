function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(rand),MersenneTwister,Vector{Float64},Int64,Int64})   # time: 0.0071572
    Base.precompile(Tuple{typeof(rand),MersenneTwister,Vector{Float32},Int64,Int64})   # time: 0.0061043
    Base.precompile(Tuple{typeof(rand),MersenneTwister,Type{Bool},Int64,Int64})   # time: 0.0015785
end
