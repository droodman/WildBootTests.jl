using SnoopPrecompile, StableRNGs
@precompile_setup begin
  for T in (Float32, Float64)
    predexog = rand(T, 1000, 4)
    resp = rand(T, 1000)
    predendog = rand(T, 1000, 1)
    inst = rand(T, 1000, 2)
    idcoarse =  floor.(Int, collect(0:999)/T(100))
    idgranular =  floor.(Int, collect(0:999)/T(2))
    feid = mod.(collect(0:999), 100)
    
    @precompile_all_calls begin
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, clustid=idcoarse)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, clustid=idgranular)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, rng=StableRNG(1), small=false)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, clustid=idcoarse, R1=T[0 0 1 0], r1=T[.2])
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, feid)
      wildboottest(T, T[0 0 0 0 1], T[.04]; resp, predexog, predendog, inst, clustid=idcoarse)
      wildboottest(T, T[0 0 0 0 1], T[.04]; resp, predexog, predendog, inst, clustid=idcoarse, arubin=true)
      wildboottest(T, T[0 0 0 1 0], T[.04]; resp, predexog, predendog, inst, clustid=idgranular, liml=true)
      nothing
    end
  end
end