function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(_erfcinv),Float32})   # time: 0.022767
    Base.precompile(Tuple{typeof(_erfcinv),Float64})   # time: 0.0082021
end
