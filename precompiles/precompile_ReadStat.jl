function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(readfield!),DataValueVector{Int16},Int64,ReadStatValue})   # time: 0.002808
    Base.precompile(Tuple{typeof(readfield!),DataValueVector{Int8},Int64,ReadStatValue})   # time: 0.002005
    Base.precompile(Tuple{typeof(readfield!),DataValueVector{String},Int64,ReadStatValue})   # time: 0.0015393
    Base.precompile(Tuple{typeof(readfield!),DataValueVector{Int32},Int64,ReadStatValue})   # time: 0.0013232
    Base.precompile(Tuple{typeof(readfield!),DataValueVector{Float32},Int64,ReadStatValue})   # time: 0.0012853
end
