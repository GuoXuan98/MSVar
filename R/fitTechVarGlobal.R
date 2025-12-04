# Estimate technical variance globally.

#' Globally fit the technical variance and bias for each sample
#'
#' Given a \code{\link{proObj}}, this function globally fit the nonlinear model
#' for technical variance and bias in each sample by using the local results
#' stored in \code{resTechLocal} field.
#'
#' @param x A \code{\link{proObj}} with an available \code{resTechLocal} field
#'   in its \code{resTech} field.
#' @param var.fit.method A character string indicating the method to be used for
#'   fitting variance curve globally. One of \code{"natural_cubic_spline"} or
#'   \code{"log_trans"} (default) or \code{"local_regression"}.
#' @param bias.fit.method A character string indicating the method to be used
#'   for fitting bias curve globally. Either \code{"natural_cubic_spline"}
#'   (default) or \code{"local_regression"}.
#' @param used.global.fit A vector of 2 with each element between 0 and 1
#'   indicates which parts of sliding windows could be used for globally curve
#'   fitting. Since the locally estimated parameters in first and end parts of
#'   the sliding window may not be stable and unreliable, \code{used.global.fit}
#'   can specify that only this part of the sliding windows are used. The
#'   default is \eqn{[0.05, 0.95]}.
#' @param correct.M Whether to correct M-value by fitted global bias curve?
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return A \code{\link{proObj}} with globally estimated technical variance and
#'   normalization bias.
#' @return \code{estTechVarGlobal} returns a \code{\link{proObj}} object, it has
#'   has an added (updated) \code{resTechGlobal} in \code{resTech} field
#'   describing globally estimated technical variance and normalization bias. If
#'   \code{correct.M} is TRUE, the \code{\link{proObj}} object would have an
#'   added \code{adj.MValue} in addition. The \code{resTechGlobal} field is
#'   itself a list of all samples, and each one is consist of the following
#'   components:
#'     \describe{
#'         \item{\code{nlc.var}}{A list of all samples, each of them records the
#'          globally fitted technical variance prediction function of current
#'          sample.}
#'         \item{\code{nlc.bias}}{A list of all samples, each of them records the
#'          globally fitted normalization bias prediction function of current
#'          sample.}
#'         \item{\code{est.tech.var}}{A matrix whose number of rows is the same
#'         as number of all proteins and number of columns is the same as sample
#'         size. It indicates the final estimated technical variance.}
#'         \item{\code{est.bias}}{A matrix whose number of rows is the same
#'         as number of all proteins and number of columns is the same as sample
#'         size. It indicates the final estimated normalization bias.}
#'     }
#'   \code{adj.MValue} is a matrix whose number of rows is the same as number of
#'   all proteins and number of columns is the same as sample size.
#'   \code{MValue} is adjusted by normalization bias.
#'
#' @export
#'
estTechVarGlobal <- function(x,
                             var.fit.method = "log_trans",
                             bias.fit.method = "log_trans",
                             used.global.fit = c(0.05, 0.95),
                             correct.M = T,
                             PLOT = F,
                             outdir = "./") {
  if (!is(x, "proObj")) {
    stop("x must be of the class \"proObj\".")
  }
  if (is.null(x$resTech$resTechLocal)) {
    stop(paste("No local fitted technical variance is found in x.",
               "Try calling estTechVarLocal() to initialize one.",
               sep = "\n"))
  }
  if (!is.null(x$resTech$resTechGlobal)) {
    x$resTech$resTechGlobal <- NULL
  }
  if (!is.null(x$adj.MValue)) {
    x$adj.MValue <- NULL
  }
  if (!is.null(x$resBio)) {
    x$resBio <- NULL
  }
  if (!is.null(x$M.Transfer)) {
    x$M.Transfer <- NULL
  }

  resTechLocal <- x$resTech$resTechLocal

  resTechVarGlobal <- lapply(resTechLocal, estTechParaGlobal_Each,
                             var.fit.method,
                             used.global.fit,
                             para = "var", PLOT, outdir)
  est.var <- lapply(resTechVarGlobal, function(x) x$est.var)
  est.var <- as.data.frame(est.var)
  nlc.var <- lapply(resTechVarGlobal, function(x) x$nlc.var)
  nlc.var.R2 <- sapply(nlc.var, function(x) attr(x, "R2"))

  resTechBiasGlobal <- lapply(resTechLocal, estTechParaGlobal_Each,
                              bias.fit.method,
                              used.global.fit,
                              para = "bias", PLOT, outdir)
  est.bias <- lapply(resTechBiasGlobal, function(x) x$est.bias)
  est.bias <- as.data.frame(est.bias)
  nlc.bias <- lapply(resTechBiasGlobal, function(x) x$nlc.bias)
  nlc.bias.R2 <- sapply(nlc.bias, function(x) attr(x, "R2"))

  rownames(est.var) <- rownames(est.bias) <- x$pro.name
  colnames(est.var) <- colnames(est.bias) <- colnames(x$MValue)

  if (PLOT) {
    barplot_NLC_R2(nlc.var.R2,
                   para = "var", para.fit.method = var.fit.method,
                   main = x$name, outdir = outdir)
    barplot_NLC_R2(nlc.bias.R2,
                   para = "bias", para.fit.method = bias.fit.method,
                   main = x$name, outdir = outdir)
  }

  resTechGlobal <- list(nlc.var = nlc.var,
                        nlc.bias = nlc.bias,
                        est.tech.var = est.var,
                        est.bias = est.bias)

  x$resTech$resTechGlobal <- resTechGlobal

  if (correct.M) {
    x$adj.MValue <- x$MValue - est.bias
  }

  x
}

#' Globally fit the technical variance and bias for current sample
#'
#' By using the locally fitted technical parameters \code{resTechLocal}, fit the
#' technical variance and normalization bias globally for current sample. The
#' nonlinear fitted curve is a function of \code{AValue}.
#'
#' @param resTechLocal A list which includes all the locally fitted technical
#'   results.
#' @param para A character string to clarify which parameter would be fitted in
#'   this round. Either \code{"var"} or \code{"bias"}.
#' @param para.fit.method A character string indicating the method to be used
#'   for fitting the nonlinear curve globally. When \code{para} is \code{"var"},
#'   \code{para.fit.method} should be one of \code{"natural_cubic_spline"},
#'   \code{"log_trans"} or \code{"local_regression"}; and when \code{para} is
#'   \code{"bias"}, \code{para.fit.method} should be either
#'   \code{"natural_cubic_spline"} or \code{"local_regression"}
#' @param used.global.fit A vector of 2 with each element between 0 and 1
#'   indicates which parts of sliding windows could be used for globally curve
#'   fitting. Since the locally estimated parameters in first and end parts of
#'   the sliding window may not be stable and unreliable, \code{used.global.fit}
#'   can specify that only this part of the sliding windows are used. The
#'   default is \eqn{[0.05, 0.95]}.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return A list include globally estimated technical curves and the
#'   corresponding predicted parameters.
estTechParaGlobal_Each <- function(resTechLocal, para = c("var", "bias"),
                                   para.fit.method,
                                   used.global.fit,
                                   PLOT = FALSE, outdir) {
  MA.linear.fit <- resTechLocal$MA.linear.fit
  sample <- colnames(resTechLocal$log2.norm.int)[1]
  AValue <- resTechLocal$log2.norm.int$AValue

  if (para == "var") {
    estTechParaGlobal <- Var_NonLinear_Fit(para.fit.method = para.fit.method,
                                           windowMean = MA.linear.fit$window.mean.int,
                                           windowVar = MA.linear.fit$window.var,
                                           used.global.fit = used.global.fit,
                                           AValue = AValue, PLOT = PLOT,
                                           outdir = outdir, sample = sample)
  } else if (para == "bias") {
    estTechParaGlobal <- Bias_NonLinear_Fit(para.fit.method = para.fit.method,
                                            windowMean = MA.linear.fit$window.mean.int,
                                            windowBias = MA.linear.fit$window.bias,
                                            used.global.fit = used.global.fit,
                                            AValue = AValue, PLOT = PLOT,
                                            outdir = outdir, sample = sample)
  }
  estTechParaGlobal
}
