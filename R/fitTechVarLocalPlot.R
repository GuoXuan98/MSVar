# Draw figures while local linear fitting for technical variance.

#' The expression of linear regression function
#'
#' @param coef The coefficients of the fitted linear model.
#'
#' @return The expression of linear model.
linear_regression_function <- function(coef) {
  if (coef[1] > 0) {
    fit_function <- sprintf("y = %.2f x + %.2f", coef[2], coef[1])
  } else if (coef[1] < 0) {
    fit_function <- sprintf("y = %.2f x - %.2f", coef[2], abs(coef[1]))
  } else {
    fit_function <- sprintf("y = %.2f x", coef[2])
  }

  fit_function
}

# sliding window plot -----------------------------------------------------

#' Draw the MA-plot in each sliding window
#'
#' \code{MA_sliding_plot} draw the MA-plot of all proteins. The X-axis is
#' \eqn{log_2}-intensity of reference sample, named A-value. The Y-axis is
#' \eqn{log_2} fold change, named M-value. To highlight the current fitted
#' window, proteins outside the window is in \code{"lightgrey"}. In the fitting
#' window, the region of proteins used for fitting are in \code{"red"} and
#' others are in \code{"grey"}.
#'
#' @param A A vector of total A-value, the basic intensity of the observations.
#' @param M A vector of total M-value, the \eqn{log_2} fold change.
#' @param windowA A vector of A-value in the given sliding window.
#' @param windowM A vector of M-value in the given sliding window.
#' @param outdir A character which represent the output directory.
#' @param sample The name of current sample to be fitted.
#' @param window.Num The current number of sliding window.
#' @param fit.region The region of proteins to be used for fitting the linear
#'   model.
#' @param nonfit.region The region not used.
#'
#' @return This function returns \code{NULL} invisibly.
MA_sliding_plot <- function(A, M, windowA, windowM, outdir, sample, window.Num,
                            fit.region, nonfit.region) {
  MA.plot.path <- paste0(outdir, "/", sample, "/MA_plot/")

  mkdir(MA.plot.path)

  png(paste0(MA.plot.path, "window_", window.Num, "_MAplot", ".png"),
    width = 1000, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  plot(A, M,
       xlab = expression(paste(log[2], " intensity")),
       ylab = expression(paste(log[2], " ratio")),
       main = paste0("window ", window.Num, " M-A plot"),
       pch = 21, col = "lightgrey", bg = "white",
       cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowA[nonfit.region], windowM[nonfit.region],
         pch = 21, col = "black", bg = alpha("grey", 0.9),
         cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  points(windowA[fit.region], windowM[fit.region],
         pch = 21, col = "black", bg = alpha("red", 0.9),
         cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)

  legend("bottomright",legend = c("Proteins out of sliding window",
                                  "The other proteins in sliding window",
                                  "Proteins in sliding window used \nfor linear fitting"),
         pch = c(21, 21, 21), pt.bg = c("white", "grey", "red"),
         col = c("lightgrey", "black", "black"), pt.cex = 1.2)

  abline(h = 0, col = "black", lty = 2)

  dev.off()

  invisible()
}

#' Draw the linear fitting plot in each sliding window
#'
#' \code{local_linear_fitting_plot} draw the linear fitting plot of all
#' proteins. The X-axis demonstrates theoretical quantiles of standard normal
#' distribution. The Y-axis demonstrates M-values in current window. Proteins
#' used for linear fitting are in \code{"red"} and others are in \code{"grey"}.
#' The fitted linear model is also added in the plot.
#'
#' @param windowM A vector of M-value in current sliding window.
#' @param outdir A character which represent the output directory.
#' @param sample The name of current sample to be fitted.
#' @param window.Num The current number of sliding window.
#' @param fit.region The region of proteins to be used for fitting the linear
#'   model.
#' @param theoretical.quantiles A vector of theoretical quantile correspond to
#'   \code{windowM}.
#' @param R2 The goodness of fit \eqn{R^2} of linear model.
#' @param coef The coefficients of linear model.
#'
#' @return This function returns \code{NULL} invisibly.
local_linear_fitting_plot <- function(windowM, outdir, sample, window.Num,
                                      fit.region, theoretical.quantiles, R2, coef) {
  linear.regression.path <- paste0(outdir, "/", sample, "/local_linear_regression/")

  mkdir(linear.regression.path)

  linear.function <- linear_regression_function(coef)

  png(paste0(linear.regression.path, "window_", window.Num, "_linear_regression", ".png"),
      width = 1000, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  plot(theoretical.quantiles, windowM,
       xlab = "Theoretical quantile of N(0,1)",
       ylab = expression(paste(log[2], " ratio")),
       main = substitute(atop("window " ~ window.Num ~ " linear regression",
                              linear.function ~ ", " ~ R^2 ~ " = " ~ R2),
                         list(window.Num = window.Num, linear.function = linear.function, R2 = R2)),
       col = "black", bg = alpha("grey", 0.7), pch = 21,
       cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3, lwd = 0.2)

  points(theoretical.quantiles[fit.region], windowM[fit.region],
         col = "black", bg = alpha("red", 0.9), pch = 21,
         cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3, lwd = 0.2)
  abline(a = round(coef[1], 2), b = round(coef[2], 2),
         col = "orange", lty = 5, lwd = 2)
  legend("bottomright", legend = c(
    "Fitted linear model", "The other proteins",
    "Proteins used for linear fitting"),
    pch = c(NA, 21, 21), pt.bg = c(NA, "grey", "red"), lty = c(5, NA, NA),
    col = c("orange", "black", "black"), pt.cex = 1.2)

  dev.off()

  invisible()
}

# Boxplot of R2 of linear model for each sample.
#' Draw a \code{boxplot} of \eqn{R^2} of linear model for each sample
#'
#' @param linear_r2 A numeric vector of \eqn{R^2} in all sliding windows.
#' @param outdir A character which represent the output directory.
#' @param sample The name of current fitted sample.
#' @param ylim Range of Y coordinate.
#'
#' @return This function returns \code{NULL} invisibly.
r2_linear_boxplot <- function(linear_r2, outdir, sample, ylim = c(0.95, 1)) {
  write.table(linear_r2,
              file = paste0(outdir, "/", sample,
                            "/r2_linear_model_boxplot", ".txt"),
              quote = F, row.names = T, col.names = F)

  png(paste0(outdir, "/", sample, "/r2_linear_model_boxplot", ".png"),
      width = 600, height = 1000, res = 72 * 2)
  dev.control("enable")

  par(family = "sans", mar = c(6, 6, 5, 3))

  boxplot(linear_r2,
          xlab = sample,
          ylab = expression(paste(R^2, " of linear model")),
          ylim = ylim,
          cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)

  dev.off()

  invisible()
}
