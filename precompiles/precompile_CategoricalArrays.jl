function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(categorical),Vector{Int16}})   # time: 0.1703477
    Base.precompile(Tuple{typeof(categorical),Vector{Int8}})   # time: 0.1616339
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:levels, :ordered), Tuple{Nothing, Bool}},Type{CategoricalVector{String, UInt8, V, C, U} where {V, C, U}},Vector{String}})   # time: 0.0728081
    Base.precompile(Tuple{typeof(unique),CategoricalVector{Int16, UInt32, Int16, CategoricalValue{Int16, UInt32}, Union{}}})   # time: 0.0127185
    Base.precompile(Tuple{Core.kwftype(typeof(categorical)),NamedTuple{(:compress,), Tuple{Bool}},typeof(categorical),Vector{String}})   # time: 0.0084864
    Base.precompile(Tuple{typeof(copy),CategoricalVector{Int16, UInt32, Int16, CategoricalValue{Int16, UInt32}, Union{}}})   # time: 0.0077128
    Base.precompile(Tuple{typeof(copy),CategoricalVector{Int8, UInt32, Int8, CategoricalValue{Int8, UInt32}, Union{}}})   # time: 0.007601
    Base.precompile(Tuple{typeof(getindex),CategoricalVector{Int8, UInt32, Int8, CategoricalValue{Int8, UInt32}, Union{}},Vector{Int64}})   # time: 0.00703
    Base.precompile(Tuple{typeof(unique),CategoricalVector{Int8, UInt32, Int8, CategoricalValue{Int8, UInt32}, Union{}}})   # time: 0.0057476
    Base.precompile(Tuple{typeof(Base.Broadcast.broadcasted),typeof(levelcode),CategoricalVector{String, UInt8, String, CategoricalValue{String, UInt8}, Union{}}})   # time: 0.0050804
    Base.precompile(Tuple{typeof(Base.Broadcast.broadcasted),typeof(levelcode),CategoricalVector{Int8, UInt32, Int8, CategoricalValue{Int8, UInt32}, Union{}}})   # time: 0.0036352
    Base.precompile(Tuple{typeof(getindex),CategoricalVector{Int16, UInt32, Int16, CategoricalValue{Int16, UInt32}, Union{}},Vector{Int64}})   # time: 0.0023095
    Base.precompile(Tuple{typeof(getindex),CategoricalVector{Int16, UInt32, Int16, CategoricalValue{Int16, UInt32}, Union{}},Int64})   # time: 0.001259
    Base.precompile(Tuple{typeof(hash),CategoricalValue{Int16, UInt32},UInt64})   # time: 0.0012191
    Base.precompile(Tuple{typeof(==),Int64,CategoricalValue{Int8, UInt32}})   # time: 0.0011278
    Base.precompile(Tuple{typeof(getindex),CategoricalVector{Int8, UInt32, Int8, CategoricalValue{Int8, UInt32}, Union{}},Int64})   # time: 0.0010941
end
