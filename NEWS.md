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
