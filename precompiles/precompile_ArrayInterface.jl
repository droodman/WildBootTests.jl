function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(static_length),OptionallyStaticUnitRange{StaticInt{0}, Int64}})   # time: 0.0022504
end
