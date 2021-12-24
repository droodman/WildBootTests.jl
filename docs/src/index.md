WildBootTests.jl performs wild bootstrap-based hypothesis tests at extreme speed. It is intended mainly for linear models: ordinary least squares (OLS) and instrumental variables/two-stage least squares (IV/2SLS). For an introduction to the wild bootstrap and the algorithms deployed here, see [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf).

The package offers and/or supports:
* The wild bootstrap for OLS ([Wu 1986](https://doi.org/10.1214/aos/1176350142)).
* The Wild Restricted Efficient bootstrap (WRE) for IV/2SLS/LIML ([Davidson and MacKinnon 2010](https://doi.org/10.1198/jbes.2009.07221)).
* The subcluster bootstrap ([MacKinnon and Webb 2018](https://doi.org/10.1111/ectj.12107)).
* Non-bootstrapped Wald, Rao, and Anderson-Rubin tests, optionally with multiway clustering.
* Confidence intervals formed by inverting the test and iteratively searching for bounds.
* Multiway clustering.
* Arbitrary and multiple linear hypotheses in the parameters.
* Maintained linear constraints on the model (restricted OLS, IV/2SLS/LIML).
* One-way fixed effects.
* Generation of data for plotting of confidence curves or surfaces after one- or two-dimensional hypothesis tests.

WildBootTests.jl incorporates order-of-magnitude algorithmic speed-ups developed since [Roodman et al. (2019)](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf) for [OLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1586107-boottest-just-as-wild-10x-faster) and [IV/2SLS](https://www.statalist.org/forums/forum/general-stata-discussion/general/1597888-boottest-~100x-faster-after-iv-gmm). And it exploits the efficiency of Julia, for example by offering single-precision (`Float32`) computation.

The interface is low-level: the exported function `wildboottest()` accepts scalars, vectors, and matrices, not [DataFrame](https://github.com/JuliaData/DataFrames.jl)s or results from estimation functions such as [lm()](https://juliastats.org/GLM.jl/v1.5/). This design minimizes the package's dependency footprint while making the core functionality available to multiple programming environments, including Julia, R (through [JuliaConnectoR](https://cran.r-project.org/web/packages/JuliaConnectoR/index.html)), and Python (through [PyJulia](https://github.com/JuliaPy/pyjulia)). A separate package will provide a higher-level Julia interface.

`wildboottest()` accepts many optional arguments. Most correspond to options of the Stata package `boottest`, which are documented in [Roodman et al. (2019), §7](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/qed_wp_1406.pdf#page=28). Julia-specific additions include an optional first argument `T`, which can be `Float32` or `Float64` to specify the precision of computation; and `rng`, which takes a random number generator such as `MersenneTwister(2302394)`.

## On latency
The first time you run `wildboottest()` in a session, Julia's just-in-time compilation will take ~10 seconds. The same will happen the first time you switch between turbo and non-turbo modes or between Float32 and Float64 calculations, or between OLS and IV/2SLS estimation. (Non-turbo and Float32 are defaults.)