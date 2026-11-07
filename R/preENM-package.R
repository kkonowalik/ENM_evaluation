#' preENM: data preparation and model selection for ENM/SDM
#'
#' Main functions:
#' * [evaluate_models()]: AUC, MAE and bias for one or many models, Pareto
#'   ranking and bootstrap robustness.
#' * [plot_pareto()]: publication-ready AUC-MAE plot with the Pareto front.
#' * [extract_predictions()]: model predictions at evaluation points.
#' * [pareto_rank()]: generic multi-criteria Pareto ranking.
#' * [simulate_enm_data()]: example data.
#'
#' Please cite Konowalik & Nosol (2021) \doi{10.1038/s41598-020-80062-1}
#' when using this package (see `citation("preENM")`).
#'
#' @keywords internal
"_PACKAGE"
