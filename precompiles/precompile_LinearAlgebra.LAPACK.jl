function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(syevr!),Char,Char,Char,Matrix{Float32},Float64,Float64,Int64,Int64,Float64})   # time: 0.0175139
    Base.precompile(Tuple{typeof(syevr!),Char,Char,Char,Matrix{Float64},Float64,Float64,Int64,Int64,Float64})   # time: 0.0154793
end
