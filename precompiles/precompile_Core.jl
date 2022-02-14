function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{NamedTuple{_A, T} where T<:Tuple} where _A,NTuple{27, Base.AbstractVector{T} where T}})   # time: 0.0073286
    Base.precompile(Tuple{Type{NamedTuple{_A, T} where T<:Tuple} where _A,NTuple{10, Base.AbstractVector{T} where T}})   # time: 0.0061337
    Base.precompile(Tuple{Type{NamedTuple{_A, T} where T<:Tuple} where _A,NTuple{8, Base.AbstractVector{T} where T}})   # time: 0.0051536
    Base.precompile(Tuple{Type{NamedTuple{_A, T} where T<:Tuple} where _A,NTuple{7, Base.AbstractVector{T} where T}})   # time: 0.0040615
    Base.precompile(Tuple{Type{NamedTuple{_A, T} where T<:Tuple} where _A,NTuple{6, Base.AbstractVector{T} where T}})   # time: 0.00406
    Base.precompile(Tuple{Type{NamedTuple{(:year, :selfemployed, :hasinsurance, :post, :post_self), T} where T<:Tuple},NTuple{5, Base.AbstractVector{T} where T}})   # time: 0.0037473
    Base.precompile(Tuple{Type{NamedTuple{(:wage, :tenure, :ttl_exp, :collgrad, :industry), T} where T<:Tuple},NTuple{5, Base.AbstractVector{T} where T}})   # time: 0.0023074
    Base.precompile(Tuple{Core.Type{Base.Vector{T} where T},Base.Vector{Core.Float64}})   # time: 0.0021404
    Base.precompile(Tuple{Core.Type{Base.Matrix{T} where T},Base.Matrix{Core.Float64}})   # time: 0.0013016
end
