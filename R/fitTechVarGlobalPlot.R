# Plot observed slope variance / intercept bias and fitted curves.

#' Plot mean-parameter in all sliding windows and fitted curve from natural
#' cubic model
#'
#' The final curve is fitted from \code{fit_natural_cubic_model}.
#' \code{NCS_plot} would plot mean-parameter in all sliding windows and fitted
#' natural cubic spline curve.
#'
#' @param windowMean A vector of mean A-value in all sliding windows.
#' @param windowPara A vector of estimated local parameters in all sliding
#'   windows.
#' @param outdir The output directory of the plot.
#' @param sample The name of current sample to be fitted.
#' @param para A character string to clarify which parameter would be fitted in
#'   this round. Either \code{"var"} or \code{"bias"}.
#' @param spc_info A character string containing the specific path to output the
#'   process plots of current nonlinear model.
#' @param outname A character string to describe more about which kind of
#'   nonlinear model was used for fitting.
#' @param fit.region A vector of 2 specifying a subset of sliding windows as the
#'   input of the whole globally curve fitting process.
#' @param fit The prediction function which accepts a vector of means and
#'   returns the predicted parameters.
#' @param R2 The goodness of fit \eqn{R^2} of the nonlinear fitting
#' @param df Degrees of freedom of nonlinear model. See \code{df} in
#'   \code{\link[splines]{ns}} for more details.
#'
#' @return This function returns \code{NULL} invisibly.
NCS_plot <- function(windowMean, windowPara, outdir, sample,
                     para = c("var", "bias"), spc_info = "",
                     outname = "NCS_model",
                     fit.region, fit, R2, df) {
  nonlinear_dir <- paste0(outdir, "/", sample, "/nonlinear_model_", para, "/", spc_info)
  mkdir(nonlinear_dir)

  if (para == "var") {
    ylab <- expression(paste(sigma^2))
    ylim <- c(0, max(windowPara) + 0.2)
  } else if (para == "bias") {
    ylab <- expression(paste(mu))
    ylim <- c(min(windowPara) - 0.1, max(windowPara) + 0.1)
  }

  png(paste0(nonlinear_dir, outname, ".png"),
      width = 1000, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  plot(windowMean, windowPara,
       xlab = expression(paste(log[2], " intensity")),
       ylab = ylab, ylim = ylim,
       main = substitute(atop(sample, R^2 * " = " * R2 * ", df = " * df),
                         list(sample = sample, R2 = round(R2, 3), df = df)),
       pch = 21, col = "red", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowMean[fit.region], windowPara[fit.region],
         ylim = ylim,
         pch = 21, col = "#008000", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  lines(windowMean, fit(windowMean),
        ylim = ylim,
        col = "darkorange", lwd = 2.5, lty = 5)

  legend("topright", legend = c("Fitted nonlinear model",
                                "Not used for model fitting",
                                "Used for model fitting"),
         pch = c(NA, 21, 21), lty = c(5, NA, NA),
         col = c("darkorange", "red", "#008000"), pt.cex = 1.2)

  dev.off()

  invisible()
}

#' Plot mean-parameter in all sliding windows and fitted curve from
#' \eqn{log}-transformed natural cubic spline model
#'
#' The final curve is fitted from \code{fit_log_trans_model}, so the parameters
#' was \eqn{log}-transformed first. \code{log_trans_NCS_plot} would plot 2
#' subplots. Subplot1 would draw mean-\eqn{log}-parameter in all sliding windows
#' and the final curve in \eqn{log} scale, subplot2 draws mean-parameter in all
#' sliding windows and the curve.
#'
#' @inheritParams NCS_plot
#'
#' @return This function returns \code{NULL} invisibly.
log_trans_NCS_plot <- function(windowMean, windowPara, outdir, sample,
                               para = c("var", "bias"), spc_info = "",
                               outname = "log_trans_NCS_model",
                               fit.region, fit, R2, df) {
  nonlinear_dir <- paste0(outdir, "/", sample, "/nonlinear_model_", para, "/", spc_info)
  mkdir(nonlinear_dir)

  if (para == "var") {
    ylab1 <- expression(paste(log[10], " ", sigma^2))
    ylim1 <- c(min(log(windowPara, 10)) - 0.1, max(log(windowPara, 10)) + 0.1)
    ylab2 <- expression(paste(sigma^2))
    ylim2 <- c(0, max(windowPara) + 0.2)
  } else if (para == "bias") {
    ylab <- expression(paste(mu))
    ylim <- c(min(windowPara) - 0.1, max(windowPara) + 0.1)
  }

  png(paste0(nonlinear_dir, outname, ".png"),
      width = 1000 * 2, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3), mfrow = c(1, 2))

  # subplot1
  plot(windowMean, log(windowPara, 10),
       xlab = expression(paste(log[2], " intensity")),
       ylab = ylab1, ylim = ylim1,
       main = substitute(atop(sample * ", " * log[10] ~ "transformed",
                              R^2 * " = " * R2 * ", df = " * df),
                         list(sample = sample, R2 = round(R2, 3), df = df)),
       pch = 21, col = "red", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowMean[fit.region], log(windowPara, 10)[fit.region],
         ylim = ylim1,
         pch = 21, col = "#008000", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  lines(windowMean,log(fit(windowMean), 10),
        ylim = ylim1,
        col = "darkorange", lwd = 2.5, lty = 5)

  legend("topright", legend = c("Fitted nonlinear model",
                                "Not used for model fitting",
                                "Used for model fitting"),
         pch = c(NA, 21, 21), lty = c(5, NA, NA),
         col = c("darkorange", "red", "#008000"), pt.cex = 1.2)

  # subplot2
  plot(windowMean, windowPara,
       xlab = expression(paste(log[2], " intensity")),
       ylab = ylab2, ylim = ylim2,
       main = substitute(atop(sample, R^2 * " = " * R2 * ", df = " * df),
                         list(sample = sample, R2 = round(R2, 3), df = df)),
       pch = 21, col = "red", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowMean[fit.region], windowPara[fit.region],
         ylim = ylim2,
         pch = 21, col = "#008000", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  lines(windowMean, fit(windowMean),
        ylim = ylim2,
        col = "darkorange", lwd = 2.5, lty = 5)

  legend("topright", legend = c("Fitted nonlinear model",
                                "Not used for model fitting",
                                "Used for model fitting"),
         pch = c(NA, 21, 21), lty = c(5, NA, NA),
         col = c("darkorange", "red", "#008000"), pt.cex = 1.2)

  dev.off()

  invisible()
}

#' Plot mean-parameter in all sliding windows and fitted curve from local
#' regression model
#'
#' The final curve is fitted from \code{fit_local_model}. \code{localfit_plot}
#' would plot mean-parameter in all sliding windows and fitted local regression
#' curve.
#'
#' @inheritParams NCS_plot
#' @param nn Nearest neighbor component of the smoothing parameter.
#'   \code{\link[locfit]{lp}} for more details.
#'
#' @return This function returns \code{NULL} invisibly.
localfit_plot <- function(windowMean, windowPara, outdir, sample,
                          para = c("var", "bias"), spc_info = "",
                          outname = "localfit_model",
                          fit.region, fit, R2, nn) {
  nonlinear_dir <- paste0(outdir, "/", sample, "/nonlinear_model_", para, "/", spc_info)
  mkdir(nonlinear_dir)

  if (para == "var") {
    ylab <- expression(paste(sigma^2))
    ylim <- c(0, max(windowPara) + 0.2)
  } else if (para == "bias") {
    ylab <- expression(paste(mu))
    ylim <- c(min(windowPara) - 0.1, max(windowPara) + 0.1)
  }

  png(paste0(nonlinear_dir, outname, ".png"),
      width = 1000, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  plot(windowMean, windowPara,
       xlab = expression(paste(log[2], " intensity")),
       ylab = ylab, ylim = ylim,
       main = substitute(atop(sample, R^2 * " = " * R2 * ", nn = " * nn),
                         list(sample = sample, R2 = round(R2, 3), nn = nn)),
       pch = 21, col = "red", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowMean[fit.region], windowPara[fit.region],
         ylim = ylim,
         pch = 21, col = "#008000", cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  lines(windowMean, fit(windowMean),
        ylim = ylim,
        col = "darkorange", lwd = 2.5, lty = 5)

  legend("topright", legend = c("Fitted nonlinear model",
                                "Not used for model fitting",
                                "Used for model fitting"),
         pch = c(NA, 21, 21), lty = c(5, NA, NA),
         col = c("darkorange", "red", "#008000"), pt.cex = 1.2)

  dev.off()
}

#' Plot the final estimated global variance or bias curve
#'
#' Although nonlinear curves are obtained from \code{fit_natural_cubic_model},
#' \code{fit_log_trans_model} or \code{fit_local_model}, more judgement need to
#' be done to decide whether use estimation from the curve directly. The final
#' estimation curve may not be strictly the same as the fitted curve.
#' \code{est_para_plot} would plot the final estimation curves.
#'
#' @inheritParams NCS_plot
#' @param AValue A vector of A-value. \code{AValue} would be used to predict the
#'   corresponding estimated variance from the curve.
#' @param estPara A numeric vector of final estimated parameters.
#'
#' @return This function returns \code{NULL} invisibly.
est_para_plot <- function(AValue, estPara, R2,
                          outdir, sample, para = c("var", "bias"),
                          spc_info = "", outname = "nonlinear_model_fitted") {
  nonlinear_dir <- paste0(outdir, "/", sample, "/nonlinear_model_", para, "/", spc_info)
  mkdir(nonlinear_dir)

  nonNA <- which(!is.na(AValue) & !is.na(estPara))
  AValue_nonNA <- AValue[nonNA]
  estPara_nonNA <- estPara[nonNA]

  if (para == "var") {
    ylab <- expression(paste("estimated ", sigma^2))
    ylim <- c(0, max(estPara_nonNA) + 0.2)
  } else if (para == "bias") {
    ylab <- expression(paste("estimated ", mu))
    ylim <- c(min(estPara_nonNA) - 0.1, max(estPara_nonNA) + 0.1)
  }

  png(paste0(nonlinear_dir, outname, ".png"),
      width = 1000, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  plot(sort(AValue_nonNA), estPara_nonNA[order(AValue_nonNA)],
       xlab = expression(paste("mean ", log[2], " intensity")),
       ylab = ylab, ylim = ylim,
       main = substitute(atop("estimated " * para, R^2 * " = " * R2),
                         list(para = para, R2 = round(R2, 3))),
       type = "l", lty = 1,
       col = "darkorange", cex.main = 1.6, lwd = 1.3, cex.lab = 1.4, cex.axis = 1.3)

  dev.off()

  invisible()
}

#' Plot \eqn{R^2} of nonlinear model fitting for all samples
#'
#' @param nonlinearR2 A vector of all \eqn{R^2} of nonlinear models.
#' @param para A character indicates which parameter to be estimated in the
#'   nonlinear models.
#' @param para.fit.method A character string indicating the method to be used
#'   for nonlinear model fitting globally. Either \code{"natural_cubic_spline"}
#'   or \code{"log_trans"} (default) or \code{"local_regression"}.
#' @param main The title of the plot, usually it is the name of the
#'   \code{\link{proObj}}.
#' @param outdir The out put directory.
#'
#' @return This function returns \code{NULL} invisibly.
barplot_NLC_R2 <- function(nonlinearR2, para = c("var", "bias"),
                           para.fit.method, main = "",
                           outdir) {
  NLC_R2_dir <- paste0(outdir, "/All_R2/")
  mkdir(NLC_R2_dir)

  if (para.fit.method == "natural_cubic_spline") {
    spc_info <- "NCS"
  } else if (para.fit.method == "log_trans") {
    spc_info <- "log_trans_NCS"
  } else if (para.fit.method == "local_regression") {
    spc_info <- "localfit"
  }

  R2 <- unlist(nonlinearR2)

  if (para == "var") {
    xlab <- expression(paste(R^2, " of non-linear model for ", sigma^2))
  } else if (para == "bias") {
    xlab <- expression(paste(R^2, " of non-linear model for ", mu))
  }

  write.table(R2,
              file = paste0(NLC_R2_dir, "/NLC_R2_for_", para,
                            "_", spc_info, ".txt"),
              quote = F, row.names = T, col.names = F)

  png(paste0(NLC_R2_dir, "/NLC_R2_for_", para, "_", spc_info, ".png"),
    width = 600, height = 300, res = 72)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  q1 <- quantile(R2, 0.25)
  q2 <- quantile(R2, 0.5)
  q3 <- quantile(R2, 0.75)

  hist(R2,
       breaks = 50, col = "#2980b9",
       xlab = xlab,
       ylab = "Sample number", main = main,
       xlim = c(0, 1),
       cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)

  lines(c(q1, q1), c(0, 1.5), col = "red", lwd = 3)
  lines(c(q2, q2), c(0, 1.5), col = "orange", lwd = 3)
  lines(c(q3, q3), c(0, 1.5), col = "green", lwd = 3)
  text(q1 - 0.025, 3, "Q1")
  text(q2 - 0.025, 3, "Q2")
  text(q3 - 0.025, 3, "Q3")

  dev.off()

  invisible()
}
