function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{Perm},O<:Base.Order.Ordering,Vector{Int64}})   # time: 0.0011304
    Base.precompile(Tuple{Type{Perm},O<:Base.Order.Ordering,Base.ReshapedArray{Int64, 1, SubArray{Int64, 2, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Vector{Int64}}, false}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}})   # time: 0.0010038
end
