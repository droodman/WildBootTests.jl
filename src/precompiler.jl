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
      test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, clustid=idcoarse)
      test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, clustid=idgranular)
      test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, rng=StableRNG(1), small=false)
      test = wildboottest(T, [0 0 0 1], [.04]; R1=[0 0 1 0], r1=[.2], resp, predexog, clustid=idcoarse)
      test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, feid)
      test = wildboottest(T, [0 0 0 0 1], [.04]; resp, predexog, predendog, inst, clustid=idcoarse)
      test = wildboottest(T, [0 0 0 1 0], [.04]; resp, predexog, predendog, inst, clustid=idgranular, liml=true)
    end
  end
end