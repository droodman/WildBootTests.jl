using SnoopPrecompile, StableRNGs
@precompile_setup begin
  rng = StableRNG(2093488)
  for T in (Float64,)
    predexog = rand(rng, T, 1000, 4)
    resp = rand(rng, T, 1000)
    predendog = rand(rng, T, 1000, 1)
    inst = rand(rng, T, 1000, 2)
    idcoarse =  floor.(Int, collect(0:999)/100)
    idgranular =  floor.(Int, collect(0:999)/2)
    feid = mod.(collect(0:999), 100)
    
    @precompile_all_calls begin
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog,      clustid=idcoarse)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, rng, clustid=idgranular)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, rng, small=false)
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, rng, clustid=idcoarse, R1=T[0 0 1 0], r1=T[.2])
      wildboottest(T, T[0 0 0 1]  , T[.04]; resp, predexog, rng, feid)
      wildboottest(T, T[0 0 0 0 1], T[.04]; resp, predexog, predendog, inst, rng, clustid=idcoarse)
      wildboottest(T, T[0 0 0 0 1], T[.04]; resp, predexog, predendog, inst, rng, clustid=idcoarse, arubin=true)
      wildboottest(T, T[0 0 0 1 0], T[.04]; resp, predexog, predendog, inst, rng, clustid=idgranular, liml=true)
    end
  end
end
