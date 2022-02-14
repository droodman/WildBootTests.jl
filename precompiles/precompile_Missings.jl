function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(disallowmissing),Vector{Union{Missing, Int8}}})   # time: 0.0135295
    Base.precompile(Tuple{typeof(disallowmissing),Vector{Union{Missing, String}}})   # time: 0.0109559
    Base.precompile(Tuple{typeof(disallowmissing),Vector{Union{Missing, Float32}}})   # time: 0.0106874
end
