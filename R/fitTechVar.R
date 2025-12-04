# Overall functions for fitting technical variance.

#' Estimate technical variance and normalization bias
#'
#' Given a \code{\link{proObj}}, this function estimate technical variance and
#' normalization bias for each proteins in each samples. First, it fit global
#' technical elements locally by using sliding windows for each sample. In this
#' way, local technical variances and normalization bias could be obtained. And
#' then, a globally fitted curve for variance and bias would be fitted to
#' capture the variation trends of variance and bias.
#'
#'
#' @param x A \code{\link{proObj}}.
#' @param A.method A character string indicating the method to calculate
#'   A-value.
#' @param filter.low Whether to filter proteins whose raw intensity is under 10?
#'   The default is TRUE.
#' @param offset A scalar used as an offset added to the raw intensity, avoiding
#'   infinity while log-transformation. Default is 1.
#' @param window.size A number which indicates the window size of sliding window
#'   when estimate technical variance and normalization bias locally.
#' @param step.size A number which indicates the step size of sliding window
#'   when estimate technical variance and normalization bias locally.
#' @param middle.percent A number that specifies the percentage of parameters used
#'   for linear model fitting within each sliding window.
#' @param var.fit.method A character string indicating the method to be used for
#'   fitting variance curve globally. Either \code{"natural_cubic_spline"} or
#'   \code{"log_trans"} (default) or \code{"local_regression"}.
#' @param bias.fit.method A character string indicating the method to be used
#'   for fitting bias curve globally. Either \code{"natural_cubic_spline"}
#'   (default) or \code{"local_regression"}.
#' @param used.global.fit A number which indicates which sliding windows should
#'   be used for global curve fitting.
#' @param correct.M Whether to correct M-value by fitted global bias curve?
#' @param cpus The number of CPUs requested for a parallel computation. By
#'   default, no parallel computation is applied.
#' @param PLOT Whether to draw the process plots? If PLOT is a vector whose
#'   length is more than 1, the first would be used as parameter passed to
#'   \code{estTechVarLocal} and the second to \code{estTechVarGlobal}.
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return \code{estTechVar} returns a \code{\link{proObj}} object, it has has
#'   an added (updated) \code{resTech} field, and it is a list consist of
#'   \code{resTechLocal} field indicating results of locally estimated technical
#'   effects and \code{resTechGlobal} field indicating results of globally
#'   estimated technical effects. Besides, it include a \code{MValue} field
#'   recording M-value of all samples, a \code{ref.mean} field indicating mean
#'   of all reference samples, and a \code{adj.MValue} field clarifying adjusted
#'   \code{MValue} if \code{correct.M} is TRUE.
#' @export
#'
estTechVar <- function(x,
                       A.method = "ref",
                       filter.low = T,
                       offset = 1,
                       window.size = 400, step.size = 100,
                       middle.percent = 0.5,
                       var.fit.method = "log_trans",
                       bias.fit.method = "natural_cubic_spline",
                       used.global.fit = c(0.05, 0.95),
                       correct.M = T,
                       cpus = 1, PLOT = FALSE,
                       outdir = NULL) {
  if (!is(x, "proObj")) {
    stop("x must be of the class \"proObj\".")
  }

  if (!is.null(x$resTech)) {
    x$resTech <- NULL
  }
  if (!is.null(x$resBio)) {
    x$resBio <- NULL
  }
  if (!is.null(x$MValue)) {
    x$MValue <- NULL
  }
  if (!is.null(x$ref.mean)) {
    x$ref.mean <- NULL
  }
  if (!is.null(x$adj.MValue)) {
    x$adj.MValue <- NULL
  }
  if (!is.null(x$M.Transfer)) {
    x$M.Transfer <- NULL
  }
  if (is.null(outdir)) {
    outdir <- "./est_tech_res/"
  }

  if(length(PLOT) > 1) {
    PLOT.Local <- PLOT[1]
    PLOT.Global <- PLOT[2]
  } else {
    PLOT.Local <- PLOT[1]
    PLOT.Global <- PLOT[1]
  }

  x <- estTechVarLocal(x,
                       A.method = A.method,
                       filter.low = filter.low,
                       offset = offset,
                       window.size = window.size, step.size = step.size,
                       middle.percent = middle.percent, PLOT = PLOT.Local,
                       cpus = cpus, outdir = outdir)
  x <- estTechVarGlobal(x,
                        var.fit.method = var.fit.method,
                        bias.fit.method = bias.fit.method,
                        used.global.fit = used.global.fit,
                        correct.M = correct.M,
                        PLOT = PLOT.Global,
                        outdir = outdir)

  x
}
