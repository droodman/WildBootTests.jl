function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(sort!),Vector{Int64},Int64,Int64,RadixSortAlg,Perm{_A, Base.ReshapedArray{Int64, 1, SubArray{Int64, 2, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Vector{Int64}}, false}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}} where _A<:Ordering})   # time: 0.1644129
    Base.precompile(Tuple{typeof(sort!),Vector{Int64},Int64,Int64,RadixSortAlg,Perm{_A, Vector{Int64}} where _A<:Ordering})   # time: 0.0446849
end
