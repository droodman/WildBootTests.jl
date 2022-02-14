function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(simd_outer_range),SubArray{CartesianIndex{2}, 1, ReshapedArray{CartesianIndex{2}, 1, CartesianIndices{2, Tuple{OneTo{Int64}, OneTo{Int64}}}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}, Tuple{UnitRange{Int64}}, false}})   # time: 0.0068847
end
