# run this in a fresh Julia session after having commented out the precompile lines at the bottom of WildBootTests.jl, both to avoid complications and force re-precompile
push!(LOAD_PATH, ".")
using WildBootTests

using SnoopCompileCore
tinf = @snoopi_deep begin
  for T in (Float32, Float64)
    resp, predexog, clustid = rand(T, 400), rand(T, 400, 4), Int32.(rand(1900:1920, 400))
    test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, clustid, auxwttype=WildBootTests.webb)
    teststat(test); p(test); CI(test); plotpoints(test);

    test = wildboottest(T, [0 0 0 1], [.04]; resp, predexog, clustid, reps=9999999, auxwttype=WildBootTests.webb, getCI=false)
    teststat(test); p(test);

    test = wildboottest(T, [0 0 0 1; 0 0 1 0], [0.05; -0.02]; resp, predexog, clustid, reps=9999, auxwttype=WildBootTests.webb)
    teststat(test); p(test);

    resp, predexog, clustid = rand(T, 2000), rand(T, 2000, 4), Int32.(rand(1:12, 2000))
    test = wildboottest(T, [0 1 0 0], [0]; R1=[0 0 1 0], r1=[.2], resp, predexog, clustid)
    teststat(test); p(test); CI(test);

    resp, predexog, clustid, predendog, inst = rand(T, 2000), rand(T, 2000, 3), Int32.(rand(1:12, 2000)), rand(2000), rand(2000)
    test = wildboottest(T, [0 0 0 1], [0]; resp, predexog, predendog, inst, clustid, small=false, reps=9999, ptype=WildBootTests.equaltail)
    teststat(test); p(test); CI(test);

    test = wildboottest(T, [0 0 0 1], [.0]; resp, predexog, predendog, inst, clustid, small=false, reps=9999, auxwttype=WildBootTests.webb, bootstrapc=true, ptype=WildBootTests.equaltail)
    teststat(test); p(test); CI(test);

    test = wildboottest(T, [0 0 0 1], [0]; resp, predexog, predendog, inst, clustid, small=false, ARubin=true, reps=9999)
    teststat(test); p(test); CI(test);

    test = wildboottest(T, [0 0 0 1], [0]; resp, predexog, predendog, inst, clustid, small=false, ARubin=true, reps=9999, imposenull=false)
    teststat(test); p(test); CI(test);

    test = wildboottest(T, [0 0 0 1], [0]; resp, predexog, predendog, inst, clustid, small=false, reps=0)
    teststat(test); p(test); CI(test);

    test = wildboottest(T, [0 0 0 1], [0]; resp, predexog, predendog, inst, clustid, small=false, reps=0, imposenull=false)
    teststat(test); p(test); CI(test);

    resp, predexog, clustid, predendog, inst = rand(T, 2000), rand(T, 2000, 1), Int32.(rand(1:12, 2000)), rand(2000), rand(2000,2)
    test = wildboottest(T, [0 1], [0]; resp, predexog, predendog, inst, LIML=true, clustid, small=false, reps=999)
    teststat(test); p(test); CI(test);

    resp, predexog, clustid, predendog, inst = rand(T, 2000,), rand(T, 2000, 5), Int32.(rand(1:12, 2000)), rand(2000), rand(2000,2)
    test = wildboottest(T, [0 0 0 0 0 1], [0]; resp, predexog, predendog, inst, Fuller=1, clustid, small=false, reps=9999, auxwttype=WildBootTests.webb)
    teststat(test); p(test); CI(test);

    resp, predexog, clustid, obswt, feid = rand(T, 2000,), rand(T, 2000, 3), Int32.(rand(1:12, 2000, 2)), rand(2000), Int64.(rand(1:12,2000))
    test = wildboottest(T, [0 0 1], [0]; resp, predexog, clustid, nbootclustvar=1, nerrclustvar=2, obswt, feid)
    teststat(test); p(test); CI(test);

    resp, predexog, clustid, feid = rand(T, 20000,), rand(T, 20000, 12), Int32.(rand(1:12, 20000, 2)), Int32.(rand(1:12,20000))
    test = wildboottest(T, [1 zeros(1,size(predexog,2)-1)], [0]; resp, predexog, clustid, nbootclustvar=1, nerrclustvar=2, feid, reps=9999)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [1 zeros(1,size(predexog,2)-1)], [0]; resp, predexog, clustid, nbootclustvar=1, nerrclustvar=2, feid, reps=9999)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [1 zeros(1,size(predexog,2)-1)], [0]; resp, predexog, clustid, nbootclustvar=2, nerrclustvar=2, feid, reps=9999)
    teststat(test); p(test); CI(test);


    resp, predexog, clustid = rand(T, 20000,), rand(T, 20000, 12), Int32.(rand(1:12, 20000, 2))
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid=clustid[:,1], gridpoints=[10], reps=9999)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid=clustid[:,1], reps=9999, imposenull=false)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid=clustid[:,[2,1]], nbootclustvar=2, nerrclustvar=1, reps=9999)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid=clustid[:,[2,1]], nbootclustvar=2, nerrclustvar=1, reps=9999, imposenull=false)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid= [collect(1:size(resp,1)) clustid[:,2]], nbootclustvar=1, nerrclustvar=1, reps=9999)
    teststat(test); p(test); CI(test);
    test = wildboottest(T, [0 1 zeros(1,size(predexog,2)-2)], [0]; resp, predexog, clustid= [collect(1:size(resp,1)) clustid[:,2]], nbootclustvar=1, nerrclustvar=1, reps=9999, imposenull=false)
    teststat(test); p(test); CI(test);
  end
end

using SnoopCompile
ttot, pcs = SnoopCompile.parcel(tinf);
SnoopCompile.write("/src", [pcs[end]])