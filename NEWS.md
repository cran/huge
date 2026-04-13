# huge 1.6

* Fixed `huge.roc()` F1 score: precision was computed from rates instead of
  counts, inflating F1 in sparse graphs (TP rate and FP rate were correct;
  only F1 was affected).
* Fixed `huge.select()` missing support for tiger method (no default criterion,
  RIC/StARS refit, or opt.icov extraction).
* Fixed C++ performance regression for glasso (3-4x) and tiger (7x) vs 1.4:
  replaced scalar loops with BLAS ddot calls, precomputed column norms, and
  reverted glasso screening to proven basic method.
* Fixed `plot.sim()` graphics parameter leak (`par()` shadowing and missing
  `on.exit()` restore).
* Fixed `ebic.score` field naming mismatch in documentation.
* Fixed `scr` parameter documentation (incorrectly stated MB not supported).
* Fixed `align` parameter documentation in `plot.huge` (logically inverted).
* Removed dead `fit$rss` matrix from `huge.mb()` and `huge.tiger()`.
* Added CI workflow for R CMD check on Linux and Windows.

# huge 1.5

* BLAS acceleration for graphical lasso and TIGER solvers (2-3x speedup at d >= 200).
* Added strong screening rule for graphical lasso.
* Removed RcppEigen dependency; the package now uses lightweight C++ headers with direct BLAS calls.
* Added testthat test suite.
* Added TIGER to DESCRIPTION text.
* Fixed Makevars.win BLAS linking.
* Fixed nonparanormal transform (`huge.npn`) rownames bug.
* Fixed documentation typos and parameter reference errors.
