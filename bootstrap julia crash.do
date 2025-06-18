sysuse auto, clear
ivreg2 trunk (weight turn = length rep78)
boottest {weight} {turn}, nonull noci svmat(numer) seed(1234) julia
mata st_matrix("weightdist", st_matrix("r(dist_1)") :+ `=_b[weight]')
mata st_matrix("turndist", st_matrix("r(dist_2)") :+ `=_b[turn]')

jl GetMatFromMat __000001, source(_boottest_jl.numerdist)
