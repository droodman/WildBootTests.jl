function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(sort!),SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Base.Sort.QuickSortAlg,ForwardOrdering})   # time: 0.0129526
    Base.precompile(Tuple{typeof(sort!),SubArray{Float32, 1, Matrix{Float32}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true},Base.Sort.QuickSortAlg,ForwardOrdering})   # time: 0.0125629
end
