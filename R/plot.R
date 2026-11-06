# Visualisation ---------------------------------------------------------------

#' Plot models in AUC-MAE space with the Pareto front
#'
#' Each model is a point (x = MAE, y = AUC); the best models sit towards the
#' **top-left**. The Pareto front connects the non-dominated models. With a
#' bootstrapped evaluation, bars show the confidence intervals and the fill
#' shows `p_pareto`, the share of resamples in which the model was on the
#' front.
#'
#' @param x An `enm_evaluation` object from [evaluate_models()].
#' @param labels Which models to label: `"all"`, `"front"` (Pareto models only)
#'   or `"none"`. Default: all when there are 15 models or fewer, otherwise
#'   front only.
#' @param error_bars Show bootstrap confidence intervals (if available).
#' @param zoom If `TRUE`, zoom the axes onto the two best Pareto fronts so
#'   that close competitors can be told apart. Poor models outside that range
#'   are then not shown.
#' @param base_size Base font size passed to [ggplot2::theme_minimal()].
#' @return A `ggplot` object that can be customised further or saved with
#'   [ggplot2::ggsave()].
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' res <- evaluate_models(sim$models, sim$presence, sim$absence,
#'                        n_boot = 100, seed = 1)
#' p <- plot_pareto(res)
#' p
#' # ggplot2::ggsave("pareto.png", p, width = 7, height = 5, dpi = 300)
plot_pareto <- function(x, labels = NULL, error_bars = TRUE, zoom = FALSE,
                        base_size = 12) {
  if (!inherits(x, "enm_evaluation")) {
    stop("`x` must be the output of evaluate_models().", call. = FALSE)
  }
  df <- as.data.frame(x)
  if (is.null(labels)) labels <- if (nrow(df) <= 15) "all" else "front"
  labels <- match.arg(labels, c("all", "front", "none"))
  has_boot <- "p_pareto" %in% names(df)

  df$status <- factor(ifelse(df$pareto_front == 1, "Pareto front",
                             "Dominated"),
                      levels = c("Pareto front", "Dominated"))
  front <- df[df$pareto_front == 1, , drop = FALSE]
  front <- front[order(front$MAE), , drop = FALSE]

  ink <- "#3a3a38"
  muted <- "#8a8a86"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$MAE, y = .data$AUC))

  if (nrow(front) > 1) {
    p <- p + ggplot2::geom_step(data = front, direction = "vh",
                                colour = muted, linewidth = 0.6,
                                linetype = "dashed")
  }
  if (has_boot && error_bars) {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$AUC_lower,
                                          ymax = .data$AUC_upper),
                             width = 0, colour = muted, linewidth = 0.4) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$MAE_lower,
                                           xmax = .data$MAE_upper),
                              height = 0, colour = muted, linewidth = 0.4)
  }
  if (has_boot) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(fill = .data$p_pareto,
                                       shape = .data$status),
                          size = 3.6, colour = "white", stroke = 0.8) +
      ggplot2::scale_fill_gradient(
        low = "#cde2fb", high = "#104281", limits = c(0, 1),
        name = "P(on Pareto front)\nacross bootstraps",
        guide = ggplot2::guide_colourbar(order = 2))
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(fill = .data$status,
                                       shape = .data$status),
                          size = 3.6, colour = "white", stroke = 0.8) +
      ggplot2::scale_fill_manual(values = c("Pareto front" = "#1c5cab",
                                            "Dominated" = "#a3a39e"),
                                 name = NULL, drop = FALSE)
  }
  p <- p + ggplot2::scale_shape_manual(
    values = c("Pareto front" = 21, "Dominated" = 22), name = NULL,
    drop = FALSE, guide = ggplot2::guide_legend(
      order = 1, override.aes = list(fill = "#5c5c58", colour = "#5c5c58")))

  if (labels != "none") {
    lab <- if (labels == "front") front else df
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = .data$model), colour = ink,
        size = base_size / 3.6, box.padding = 0.5, min.segment.length = 0.2,
        segment.colour = muted, segment.size = 0.3, seed = 1,
        max.overlaps = Inf)
    } else {
      p <- p + ggplot2::geom_text(data = lab,
                                  ggplot2::aes(label = .data$model),
                                  colour = ink, size = base_size / 3.6,
                                  vjust = -1.1)
    }
  }

  sub <- if (has_boot) {
    sprintf("Bars: %d%% bootstrap intervals (%d resamples). Better models are top-left.",
            round(100 * attr(x, "conf")), attr(x, "n_boot"))
  } else {
    "Better models are top-left (high AUC, low MAE)."
  }
  if (zoom) {
    z <- df[df$pareto_front <= 2, , drop = FALSE]
    xr <- range(c(z$MAE, z$MAE_lower, z$MAE_upper), na.rm = TRUE)
    yr <- range(c(z$AUC, z$AUC_lower, z$AUC_upper), na.rm = TRUE)
    pad <- function(r) r + c(-1, 1) * max(diff(r), 1e-3) * 0.15
    p <- p + ggplot2::coord_cartesian(xlim = pad(xr), ylim = pad(yr))
  }

  p +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(x = "MAE (lower is better)", y = "AUC (higher is better)",
                  title = "Model selection: AUC vs MAE", subtitle = sub) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e5e1",
                                               linewidth = 0.3),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      text = ggplot2::element_text(colour = ink),
      plot.subtitle = ggplot2::element_text(colour = muted),
      legend.position = "right")
}

#' @importFrom rlang .data
NULL

#' Plot a model selection for any pair of metrics
#'
#' Shows the models from [select_models()] in the space of two metrics.
#' Selected models are filled, Pareto-front models are circles (others are
#' squares), dashed lines mark the minimum values and bars show +/- 1 SD
#' across cross-validation runs when `<metric>_sd` columns are present. Call
#' it several times for several panels, e.g. AUC vs MAE and AUC vs fossil AUC.
#'
#' @param x An `enm_selection` object from [select_models()].
#' @param x_metric,y_metric Columns to plot on the x and y axes.
#' @param labels Label the points with model names.
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' df <- data.frame(model = c("GLM", "GBM", "RF", "SRE"),
#'                  AUC = c(0.86, 0.88, 0.87, 0.71),
#'                  TSS = c(0.62, 0.65, 0.66, 0.38),
#'                  MAE = c(0.30, 0.27, 0.24, 0.41))
#' plot_selection(select_models(df))
#' plot_selection(select_models(df), x_metric = "TSS")
plot_selection <- function(x, x_metric = "MAE", y_metric = "AUC",
                           labels = TRUE, base_size = 12) {
  if (!inherits(x, "enm_selection")) {
    stop("`x` must be the output of select_models().", call. = FALSE)
  }
  df <- as.data.frame(x)
  miss <- setdiff(c(x_metric, y_metric), names(df))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  df$.x <- df[[x_metric]]
  df$.y <- df[[y_metric]]
  df$status <- factor(ifelse(df$selected, "Selected", "Not selected"),
                      levels = c("Selected", "Not selected"))
  df$front <- factor(ifelse(df$pareto, "Pareto front", "Dominated"),
                     levels = c("Pareto front", "Dominated"))
  crit <- attr(x, "criteria")
  mins <- attr(x, "min_metrics")
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  axis_lab <- function(m) {
    dir <- if (m %in% names(crit)) crit[[m]] else NA
    if (is.na(dir)) m else paste0(m, if (dir == "max") " (higher is better)"
                                  else " (lower is better)")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y))
  if (x_metric %in% names(mins)) {
    p <- p + ggplot2::geom_vline(xintercept = mins[[x_metric]],
                                 linetype = "dashed", colour = muted,
                                 linewidth = 0.4)
  }
  if (y_metric %in% names(mins)) {
    p <- p + ggplot2::geom_hline(yintercept = mins[[y_metric]],
                                 linetype = "dashed", colour = muted,
                                 linewidth = 0.4)
  }
  xs <- paste0(x_metric, "_sd")
  ys <- paste0(y_metric, "_sd")
  if (ys %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$.y - .data[[ys]],
                   ymax = .data$.y + .data[[ys]]),
      width = 0, colour = muted, linewidth = 0.4)
  }
  if (xs %in% names(df)) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$.x - .data[[xs]],
                   xmax = .data$.x + .data[[xs]]),
      height = 0, colour = muted, linewidth = 0.4)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(fill = .data$status,
                                     shape = .data$front),
                        size = 3.6, colour = "white", stroke = 0.8) +
    ggplot2::scale_fill_manual(values = c("Selected" = "#1c5cab",
                                          "Not selected" = "#b4b3ae"),
                               name = NULL, drop = FALSE,
                               guide = ggplot2::guide_legend(
                                 override.aes = list(shape = 21))) +
    ggplot2::scale_shape_manual(values = c("Pareto front" = 21,
                                           "Dominated" = 22),
                                name = NULL, drop = FALSE,
                                guide = ggplot2::guide_legend(
                                  override.aes = list(fill = "#5c5c58")))
  if (labels) p <- p + label_layer(df, ink, muted, base_size)
  if (x_metric %in% names(crit) && crit[[x_metric]] == "min") {
    p <- p + ggplot2::scale_x_reverse()
  }
  if (y_metric %in% names(crit) && crit[[y_metric]] == "min") {
    p <- p + ggplot2::scale_y_reverse()
  }
  p +
    ggplot2::labs(x = axis_lab(x_metric), y = axis_lab(y_metric),
                  title = paste("Model selection:", y_metric, "vs", x_metric),
                  subtitle = paste0(
                    "Selection rule: ", attr(x, "rule"),
                    ". Better models are top-right.\n",
                    if (length(mins)) "Dashed lines: minimum values. ",
                    if (ys %in% names(df) || xs %in% names(df))
                      "Bars: \u00b11 SD across runs.")) +
    theme_enm(base_size, ink, muted)
}

label_layer <- function(df, ink, muted, base_size) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = df, ggplot2::aes(label = .data$model),
                             colour = ink, size = base_size / 3.6,
                             box.padding = 0.5, min.segment.length = 0.2,
                             segment.colour = muted, segment.size = 0.3,
                             seed = 1, max.overlaps = Inf)
  } else {
    ggplot2::geom_text(data = df, ggplot2::aes(label = .data$model),
                       colour = ink, size = base_size / 3.6, vjust = -1.1)
  }
}

theme_enm <- function(base_size, ink, muted) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e5e1",
                                               linewidth = 0.3),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      text = ggplot2::element_text(colour = ink),
      plot.subtitle = ggplot2::element_text(colour = muted),
      legend.position = "right")
}
