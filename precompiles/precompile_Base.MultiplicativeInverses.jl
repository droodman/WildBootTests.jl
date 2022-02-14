function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{SignedMultiplicativeInverse{Int64}},Int64})   # time: 0.002934
end
