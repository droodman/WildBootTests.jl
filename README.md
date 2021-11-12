# WildBootTests.jl

This package performs wild bootstrap-based hypothesis tests in Julia at extreme speed. It is intended mainly for linear models: OLS and IV/2SLS. For an introduction to the wild bootstrap and the algorithms deployed here see [Roodman et al. (2019)](https://journals.sagepub.com/doi/abs/10.1177/1536867X19830877?journalCode=stja). This implementation, like the [current implementation for Stata](https://ideas.repec.org/c/boc/bocode/s458121.html), incorporates order-of-magnitude algorithmic speed-ups developed since that publication. The Julia implementation brings further efficiency gains, such as by offering single- along with double-precision computation.


