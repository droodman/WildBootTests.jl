function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(iterate),PartitionIterator{OneTo{Int64}}})   # time: 0.001183
end
