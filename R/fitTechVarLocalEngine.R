# Main functions of linear model fitting based on M-A plot ---------------------

#' Calculate the theoretical quantile of standard normal distribution
#'
#' According to the start number and end number of M-values, calculate their
#' theoretical quantile of standard normal distribution.
#' @param start.pos The start of theoretical plot position.
#' @param end.pos The end of theoretical plot position.
#' @param aa A scalar of scale parameter.
#'
#' @return The theoretical plot position in standard normal distribution.
MS_quantile_position <- function(start.pos, end.pos, aa) {
  get_inner <- function(x) {
    r.digital <- (x - aa) / (x - aa + 1)
    r.digital
  }

  constant <- (end.pos - aa + 1) / (end.pos - 2 * aa + 1)

  quantile.position <- sapply(start.pos:end.pos, function(j) {
    pi.list <- vapply(j:end.pos, get_inner, numeric(1))
    constant * prod(pi.list)
  })

  quantile.position
}

# Calculate R2 ------------------------------------------------------------

#' Calculate goodness of fit \eqn{R^2}
#'
#' @param y.predict A numeric vector of prediction.
#' @param y.data A numeric vector of observation.
#'
#' @return The goodness of fit \eqn{R^2}.
getIndexes <- function(y.predict, y.data) {
  n <- length(y.data)

  SSE <- sum((y.data - y.predict)^2)
  MSE <- SSE / n
  RMSE <- sqrt(MSE)

  u <- mean(y.data)
  SST <- sum((y.data - u)^2)
  SSR <- SST - SSE

  R_square <- SSR / SST

  R_square
}

# Linear model for each sliding window

#' Fit linear model in current sliding window
#'
#' \code{window.Num} together with \code{groups} figures out the M-values used
#' for linear model fitting. \code{\link[stats]{lm}} models the linear
#' relationship between M-values and quantiles of standard normal distribution.
#' The slope of the linear model is local technical standard deviation of this
#' window and the intercept is the normalization bias.
#'
#' @param window.Num The current number of sliding window.
#' @param groups A list of vectors for indexing observations. Each vector
#'   effectively specifies a group of observations. \code{window.Num} together
#'   with \code{groups} determines the index of observations which would be used
#'   for linear model fitting.
#' @param data A matrix which records M-values and A-values, they are sorted by
#'   the order of A-values. Those observations whose M-values and A-values are
#'   \code{NA} will be removed.
#' @param sample The name of current sample to be fitted.
#' @param middle.percent A number which indicates which parameters should be
#'   used for linear model fitting in the sliding window.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return A vector records all of the results about linear model in current
#'   sliding window, which includes fitted variances, bias, mean of A-value in
#'   the window and goodness of fit \eqn{R^2}.
MA_linear_regression <- function(window.Num, groups, data,
                                 sample, middle.percent,
                                 PLOT, outdir) {
  index <- groups[[window.Num]]
  sortM <- data$sortM
  sortA <- data$sortA

  windowM <- data[index, "sortM"]
  windowA <- data[index, "sortA"]

  windowM.sort <- sort(windowM)
  windowA.sort <- windowA[order(windowM)]

  Npro.window <- length(windowM.sort)

  fit.region.start <- floor((1 - middle.percent) / 2 * Npro.window) + 1
  fit.region.end <- floor((0.5 + middle.percent / 2) * Npro.window)
  fit.region.length <- fit.region.end - fit.region.start + 1

  fit.region <- fit.region.start:fit.region.end
  nonfit.region <- c(1:(fit.region.start - 1), (fit.region.end + 1):Npro.window)

  quantile.position <- MS_quantile_position(
    start.pos = 1,
    end.pos = Npro.window,
    aa = 0.375
  )
  theoretical.quantiles <- qnorm(quantile.position, 0, 1)

  fit <- lm(windowM.sort[fit.region] ~ theoretical.quantiles[fit.region])

  R2 <- getIndexes(predict(fit, data.frame(theoretical.quantiles[fit.region])),
                   windowM.sort[fit.region])

  coef <- coefficients(fit)
  res <- c(window.bias = coef[1],
           window.var = coef[2]^2,
           window.mean.int = mean(windowA, na.rm = T),
           window.R2 = R2)

  if (PLOT) {
    MA_sliding_plot(sortA, sortM,
                    windowA = windowA.sort, windowM = windowM.sort,
                    outdir = outdir, sample = sample, window.Num = window.Num,
                    fit.region = fit.region, nonfit.region = nonfit.region)

    local_linear_fitting_plot(windowM = windowM.sort,
                              outdir = outdir, sample = sample,
                              window.Num = window.Num,
                              fit.region = fit.region,
                              theoretical.quantiles = theoretical.quantiles,
                              R2 = R2, coef = coef)

    r2_linear_boxplot(R2, outdir = outdir, sample = sample)
  }
  res
}
