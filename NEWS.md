# huge 1.5

* BLAS acceleration for graphical lasso and TIGER solvers (2-3x speedup at d >= 200).
* Added strong screening rule for graphical lasso.
* Removed RcppEigen dependency; the package now uses lightweight C++ headers with direct BLAS calls.
* Added testthat test suite.
* Added TIGER to DESCRIPTION text.
* Fixed Makevars.win BLAS linking.
* Fixed nonparanormal transform (`huge.npn`) rownames bug.
* Fixed documentation typos and parameter reference errors.
