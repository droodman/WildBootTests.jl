function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(pick_vector_width),Type{Float32}})   # time: 0.0074239
end
