# preENM 2.0.0

* The package is renamed from `ENMselect` to **`preENM`** and the repository
  moved to <https://github.com/kkonowalik/preENM>. The package will grow into a
  toolkit for data preparation for ENM/SDM (point filtering, background point
  selection, environmental thinning, correlation analysis); the existing
  evaluation and model-selection functions are unchanged. Replace
  `library(ENMselect)` with `library(preENM)`.
* New hex logo.
* Horizontal error bars use `geom_errorbar(orientation = "y")` instead of
  `geom_errorbarh()`, which is deprecated in ggplot2 4.0.
* `Metrics` is declared in Suggests (used only in a test that checks AUC, MAE
  and bias against it).

# ENMselect 1.2.0

* Choose which metrics to compare: `evaluate_models(criteria = ...)` accepts
  any metrics with a direction, `"max"`, `"min"` or the new `"zero"` (closer
  to 0 is better, e.g. `bias`). The default is still AUC vs MAE. The Pareto
  front, compromise distance, bootstrap intervals and `p_pareto` all follow
  the chosen criteria; `tol` is also available.
* Add your own metrics with `metrics = list(name = function(obs, pred) ...)`
  in `evaluate_models()` and `evaluate_cv()`.
* `plot_pareto()` gains `x_metric` / `y_metric`, defaults to the chosen
  criteria and orients axes so that better models are top-right (the MAE axis
  is now reversed, matching `plot_selection()`).
* `select_models()` and `pareto_rank()` support "closer to 0 is better"
  criteria (`"zero"` / `zero =`).
* `evaluate_models(weights = )` now defaults to equal weights for whichever
  criteria are used.

# ENMselect 1.1.0

New tools generalised from a full biomod2 workflow (spatial cross-validation,
fossil validation and selection of algorithms for an ensemble):

* `spatial_folds()`: spatial block cross-validation folds, balanced for
  presences and background points.
* `evaluate_cv()`: hold-out metrics per cross-validation run, summarised per
  model across folds and background sets, with diagnostics (SD between sets
  vs between folds, AUC per fold) and optional fixed thresholds for TSS.
* `temporal_validation()`: validation of hindcasts against dated records
  (e.g. fossils) ranked against background cells of the same time slice, by
  period, for pooled periods, for analog climates only, and as sensitivity.
* `select_models()`: selection by conventional minimum values or by the
  Pareto front (with tolerance for ties), mean rank, manual override.
* `plot_selection()`: selection figure for any pair of metrics.
* New metrics `tss_score()`, `boyce_index()` (identical to
  `ecospat::ecospat.boyce()`), `threshold_p10()`.
* Helpers `anomaly_climate()` (delta method for past climate),
  `outside_range()` (non-analog climate), `age_to_period()`.
* `pareto_rank()` gains `tol`.
* `evaluate_models()` now also reports TSS and the Boyce index.
* `simulate_enm_data()` also returns environmental layers.

# ENMselect 1.0.0

* First release: `evaluate_models()`, `pareto_rank()`, `plot_pareto()`,
  `extract_predictions()`, AUC/MAE/bias, bootstrap robustness of selection.
