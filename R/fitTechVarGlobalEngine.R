# Fit nonlinear model for global technical variance
# according to local technical variance.


# Main functions of nonlinear model fitting ------------------------------

#' \code{Var_NonLinear_Fit} returns estimated technical variance
#'
#' Technical variance is a function of reference for each sample, nonlinear
#' models would be fitted to describe the relationship between them. By using
#' \code{estTechVarLocal}, technical variance in each sliding window are
#' obtained and stored as \code{windowVar} and \code{windowMean}. Since the
#' locally estimated parameters in first and end parts of the sliding windows
#' may be unstable and unreliable, \code{used.global.fit} can specify that only
#' this part of the sliding windows are used during the curve fitting process.
#'
#' While fitting the nonlinear curve, there are 3 methods to choose.
#' \code{"natural_cubic_spline"} and \code{"log_trans"} actually use natural
#' cubic spline, but \code{"local_regression"} uses local regression. No matter
#' which method is used, 2 or 3 models would be fitted depends on \code{choice}
#' and the one with higher goodness of fit \eqn{R^2} is the final nonlinear
#' model.
#'
#' After fitting the curve, \code{AValue} is used as input of prediction of
#' technical variance and more judgement need to be done to decide whether use
#' estimation from the curve directly. If \code{AValue} is more than or less
#' than the boundary of \code{windowMean} which were used for nonlinear model
#' fitting, the positive predictions at the boundary would be set as the
#' corresponding estimated variance. Besides, if the goodness of fit \eqn{R^2}
#' of this nonlinear curve is less than 0.5, mean of \code{windowVar} would be
#' set as estimated variance for all observations.
#'
#' @param para.fit.method A character string indicating the method to be used
#'   for fitting the variance curve globally. One of
#'   \code{"natural_cubic_spline"}, \code{"log_trans"} (default) or
#'   \code{"local_regression"}
#' @param windowMean A vector of mean A-value in all sliding windows.
#' @param windowVar A vector of estimated local variance in all sliding windows.
#' @param used.global.fit A vector of 2 with each element between 0 and 1
#'   specifying a subset of sliding windows as the input of the whole globally
#'   curve fitting process. Since the locally estimated parameters in first and
#'   end parts of the sliding windows may be unstable and unreliable,
#'   \code{used.global.fit} can specify that only this part of the sliding
#'   windows are used. The default is \code{c(0.05, 0.95)}.
#' @param choice If \code{para.fit.method} is \code{"natural_cubic_spline"},
#'   \code{"log_trans"}, \code{choice} is a numeric to specify how many natural
#'   cubic spline nonlinear models to fit to choose from. By default, 2 non
#'   linear models with \code{df} in \code{\link[splines]{ns}} would be fitted
#'   and the one with higher goodness of fit \eqn{R^2} would be served as the
#'   final nonlinear model. Otherwise, if \code{para.fit.method} is
#'   \code{"local_regression"}, \code{choice} is a numeric to specify how many
#'   local regression nonlinear models to fit to choose from. By default, 2 non
#'   linear models with \code{nn} in \code{\link[locfit]{lp}} would be fitted
#'   and the one with higher goodness of fit \eqn{R^2} would be served as the
#'   final nonlinear model.
#' @param args.locfit A named list of extra arguments to
#'   \code{\link[locfit]{locfit}}. Only used when \code{method} is
#'   \code{"local"}.
#' @param AValue A vector of A-value. \code{AValue} would be used to predict the
#'   technical variance from the nonlinear variance prediction function. If
#'   \code{AValue} is more than or less than the boundary of \code{windowMean}
#'   which were used for nonlinear model fitting, the predictions at the
#'   boundary would be set as the estimated variance at the endpoints.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#' @param sample The name of current sample to be fitted.
#'
#' @return A list include globally fitted technical variance prediction
#'   function and the final estimated technical variance.
Var_NonLinear_Fit <- function(para.fit.method, windowMean, windowVar,
                              used.global.fit, choice = 2, args.locfit = list(),
                              AValue, PLOT = F, outdir = "./", sample) {
  x <- match.arg(para.fit.method,
                 c("natural_cubic_spline",
                   "local_regression", "log_trans"))

  Npro <- length(AValue)
  nonNA <- which(!is.na(AValue))

  est.var <- rep(NA, Npro)

  fit.start <- max(round(used.global.fit[1] * length(windowMean)), 1)
  fit.end <- min(round(used.global.fit[2] * length(windowMean)),
                 length(windowMean))

  if (x == "natural_cubic_spline") {
    spc_info <- "NCS/"
    nlc <- fit_natural_cubic_model(windowMean = windowMean,
                                   windowPara = windowVar,
                                   used.global.fit = used.global.fit,
                                   choice = choice,
                                   PLOT = PLOT, outdir = outdir,
                                   para = "var", sample = sample,
                                   outname = "NCS_model", spc_info = spc_info)
    nonlinearR2 <- attr(nlc, "R2")

    outname <- "NCS_model_fitted"
  } else if (x == "log_trans") {
    spc_info <- "log_trans/"
    nlc <- fit_log_trans_model(windowMean = windowMean,
                               windowPara = windowVar,
                               used.global.fit = used.global.fit,
                               choice = choice,
                               PLOT = PLOT, outdir = outdir,
                               para = "var", sample = sample,
                               outname = "log_trans_NCS_model", spc_info = spc_info)
    nonlinearR2 <- attr(nlc, "R2")

    outname <- "log_trans_NCS_model_fitted"
  } else if (x == "local_regression") {
    spc_info <- "local_regression/"
    nlc <- fit_local_model(windowMean = windowMean,
                           windowPara = windowVar,
                           used.global.fit = used.global.fit,
                           choice = choice, args.locfit = args.locfit,
                           PLOT = PLOT, outdir = outdir,
                           para = "var", sample = sample,
                           outname = "localfit_model", spc_info = spc_info)
    nonlinearR2 <- attr(nlc, "R2")

    outname <- "localfit_model_fitted"
  }

  est.var <- nlc(AValue)
  est.var.start <- nlc(windowMean[fit.start])
  est.var.end <- max(nlc(windowMean[fit.end]),
                     min(est.var[which(est.var > 0)]))

  if (nonlinearR2 > 0.5) {
    A_exceed_left <- which(AValue < windowMean[fit.start])
    est.var[A_exceed_left] <- est.var.start

    A_exceed_right <- which(AValue > min(min(AValue[which(est.var < 0)], windowMean[fit.end])))
    est.var[A_exceed_right] <- est.var.end
  } else {
    est.var[nonNA] <- rep(mean(windowVar), length(nonNA))
  }

  # Plot estimated variance.
  if (PLOT) {
    est_para_plot(AValue, est.var,
                  R2 = nonlinearR2[length(nonlinearR2)],
                  outdir = outdir, sample = sample,
                  para = "var", spc_info = spc_info,
                  outname = outname)
  }

  list(est.var = est.var,
       nlc.var = nlc,
       nlc.var.R2 = nonlinearR2)
}

#' \code{Bias_NonLinear_Fit} returns estimated normalization bias
#'
#' @inheritParams Var_NonLinear_Fit
#' @param para.fit.method A character string indicating the method to be used
#'   for fitting the bias curve globally. One of \code{"natural_cubic_spline"}
#'   (default) or \code{"local_regression"}
#' @param windowBias A vector of estimated local bias in all sliding windows.
#'
#' @return A list include globally fitted normalization bias prediction function
#'   and the final estimated normalization bias.
Bias_NonLinear_Fit <- function(para.fit.method, windowMean, windowBias,
                               used.global.fit, choice = 2, args.locfit = list(),
                               AValue, PLOT = F, outdir, sample) {
  x <- match.arg(para.fit.method,
                 c("natural_cubic_spline", "local_regression"))

  Npro <- length(AValue)

  fit.start <- max(round(used.global.fit[1] * length(windowMean)), 1)
  fit.end <- min(round(used.global.fit[2] * length(windowMean)),
                 length(windowMean))

  if (x == "natural_cubic_spline") {
    spc_info <- "NCS/"
    nlc <- fit_natural_cubic_model(windowMean = windowMean,
                                   windowPara = windowBias,
                                   used.global.fit = used.global.fit,
                                   choice = choice,
                                   PLOT = PLOT, outdir = outdir,
                                   para = "bias", sample = sample,
                                   outname = "NCS_model", spc_info = spc_info)
    nonlinearR2 <- attr(nlc, "R2")

    outname <- "NCS_model_fitted"
  } else if (x == "local_regression") {
    spc_info <- "local_regression/"
    nlc <- fit_local_model(windowMean = windowMean,
                           windowPara = windowBias,
                           used.global.fit = used.global.fit,
                           choice = choice, args.locfit = args.locfit,
                           PLOT = PLOT, outdir = outdir,
                           para = "bias", sample = sample,
                           outname = "localfit_model", spc_info = spc_info)
    nonlinearR2 <- attr(nlc, "R2")

    outname <- "localfit_model_fitted"
  }

  est.bias <- nlc(AValue)
  est.bias.start <- nlc(windowMean[fit.start])
  est.bias.end <- nlc(windowMean[fit.end])

  if (nonlinearR2 > 0.5) {
    A_exceed_left <- which(AValue < windowMean[fit.start])
    est.bias[A_exceed_left] <- est.bias.start

    A_exceed_right <- which(AValue > windowMean[fit.end])
    est.bias[A_exceed_right] <- est.bias.end
  } else {
    est.bias <- rep(mean(windowBias), Npro)
  }

  # Plot estimated bias.
  if (PLOT) {
    est_para_plot(AValue, est.bias,
                  R2 = nonlinearR2[length(nonlinearR2)],
                  outdir = outdir, sample = sample,
                  para = "bias", spc_info = spc_info,
                  outname = outname)
  }

  list(est.bias = est.bias,
       nlc.bias = nlc,
       nlc.bias.R2 = nonlinearR2)
}


# The fundamental functions for model fitting -----------------------------

#' Fit global mean-para curve by natural cubic spline
#'
#' @param windowMean A vector of mean A-value in all sliding windows.
#' @param windowPara A vector of estimated local parameters in all sliding
#'   windows.
#' @param used.global.fit A vector of 2 with each element between 0 and 1
#'   specifying a subset of sliding windows as the input of the whole globally
#'   curve fitting process. Since the locally estimated parameters in first and
#'   end parts of the sliding windows may be unstable and unreliable,
#'   \code{used.global.fit} can specify that only this part of the sliding
#'   windows are used. The default is \code{c(0.05, 0.95)}.
#' @param choice A numeric to specify how many natural cubic spline nonlinear
#'   models to fit to choose from. By default, 2 non linear models with
#'   \code{df} in \code{\link[splines]{ns}} would be fitted and the one with
#'   higher goodness of fit \eqn{R^2} would be served as the final nonlinear
#'   model.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#' @param sample The name of current sample to be fitted.
#' @param para A character string to clarify which parameter would be fitted in
#'   this round. Either \code{"var"} or \code{"bias"}.
#' @param spc_info A character string containing the specific path to output the
#'   process plots of current nonlinear model.
#' @param outname A character string to describe more about which kind of
#'   nonlinear model was used for fitting.
#'
#' @return A prediction function which accepts a vector of means and returns the
#'   predicted parameters.
fit_natural_cubic_model <- function(windowMean, windowPara,
                                    used.global.fit = c(0.05, 0.95),
                                    choice = 2,
                                    PLOT = FALSE, outdir,
                                    sample, para = c("var", "bias"),
                                    spc_info = "", outname = "NCS_model") {
  fit.range <- c(max(floor(used.global.fit[1] * length(windowMean)) + 1, 1),
                 min(floor(used.global.fit[2] * length(windowMean)), length(windowMean)))
  fit.region <- fit.range[1]:fit.range[2]

  x0 <- windowMean[fit.region]
  y0 <- windowPara[fit.region]

  NCS.res <- NCS_df_selection(x0 = x0, y0 = y0,
                              choice = choice)
  fit <- NCS.res$fit

  nlc <- function(x) {
    nD <- data.frame(x0 = x)
    predict(fit, nD)
  }
  attr(nlc, "R2") <- NCS.res$R2
  attr(nlc, "df") <- NCS.res$df

  if (PLOT) {
    NCS_plot(windowMean, windowPara, outdir, sample,
             para = para, spc_info = spc_info,
             outname = outname,
             fit.region, nlc, NCS.res$R2, NCS.res$df)
  }

  nlc
}

#' Fit global mean-para trend by \eqn{log}-transformed natural cubic spline
#'
#' While fitting global mean-var curve, it should be ensured that all the fitted
#' variance curves are positive. \code{fit_log_trans_model} first transform all
#' \code{windowPara} in \eqn{log} scale, and then fit a mean-\eqn{log var} curve
#' by natural cubic spline. Finally it would return a re-scaled prediction
#' function in exponential level. This process can not only make sure that all
#' the data read out of the curve are positive, but also that the curve is
#' smoother because of the transformation.
#'
#' \code{fit_log_trans_model} is more suitable to fit the technical variance
#' curve.
#'
#' @inheritParams fit_natural_cubic_model
#'
#' @return A prediction function which accepts a vector of means and returns the
#'   predicted parameters.
fit_log_trans_model <- function(windowMean, windowPara,
                                used.global.fit = c(0.05, 0.95),
                                choice = 2,
                                PLOT = FALSE, outdir, sample,
                                para = c("var", "bias"),
                                spc_info = "", outname = "log_trans_NCS_model") {
  fit.range <- c(max(floor(used.global.fit[1] * length(windowMean)) + 1, 1),
                 min(floor(used.global.fit[2] * length(windowMean)), length(windowMean)))
  fit.region <- fit.range[1]:fit.range[2]

  log_windowPara <- log(windowPara, base = 10)
  x0 <- windowMean[fit.region]
  y0 <- log_windowPara[fit.region]

  NCS.res <- NCS_df_selection(x0 = x0, y0 = y0,
                              choice = choice)
  fit <- NCS.res$fit

  nlc <- function(x) {
    nD <- data.frame(x0 = x)
    10^(predict(fit, nD))
  }
  attr(nlc, "R2") <- NCS.res$R2
  attr(nlc, "df") <- NCS.res$df

  if (PLOT) {
    log_trans_NCS_plot(windowMean, windowPara, outdir, sample,
                       para = para, spc_info = spc_info,
                       outname = outname,
                       fit.region, nlc, NCS.res$R2, NCS.res$df)
  }

  nlc
}


#' Fit global mean-para trend by local regression
#'
#' @inheritParams fit_natural_cubic_model
#' @param choice A numeric to specify how many local regression nonlinear
#'   models to fit to choose from. By default, 2 non linear models with
#'   \code{nn} in \code{\link[locfit]{lp}} would be fitted and the one with
#'   higher goodness of fit \eqn{R^2} would be served as the final nonlinear
#'   model.
#' @param args.locfit A named list of extra arguments to
#'   \code{\link[locfit]{locfit}}. Only used when \code{method} is
#'   \code{"local"}.
#'
#' @return A prediction function which accepts a vector of means and returns the
#'   predicted parameters.
fit_local_model <- function(windowMean, windowPara,
                            used.global.fit = c(0.05, 0.95), choice = 2,
                            args.locfit = list(),
                            PLOT = FALSE, outdir, sample,
                            para = c("var", "bias"),
                            spc_info = "", outname = "localfit_model") {
  fit.range <- c(max(floor(used.global.fit[1] * length(windowMean)) + 1, 1),
                 min(floor(used.global.fit[2] * length(windowMean)), length(windowMean)))
  fit.region <- fit.range[1]:fit.range[2]

  x0 <- windowMean[fit.region]
  y0 <- windowPara[fit.region]

  if (para == "var") {
    f0 <- "gamma"
  } else if (para == "bias") {
    f0 <- "gaussian"
  }

  lp_res6 <- do.call(lp, c(list(x0), nn = 0.6))
  model6 <- do.call(locfit, c(list(y0 ~ lp_res6, family = f0), args.locfit))
  R2_6 <- getIndexes(predict(model6, x0), y0)

  lp_res7 <- do.call(lp, c(list(x0), nn = 0.7))
  model7 <- do.call(locfit, c(list(y0 ~ lp_res7, family = f0), args.locfit))
  R2_7 <- getIndexes(predict(model7, x0), y0)

  lp_res8 <- do.call(lp, c(list(x0), nn = 0.8))
  model8 <- do.call(locfit, c(list(y0 ~ lp_res8, family = f0), args.locfit))
  R2_8 <- getIndexes(predict(model8, x0), y0)

  if (choice == 3) {
    if (R2_8 > R2_7 & R2_8 > R2_6) {
      R2 <- R2_8
      fit <- model8
      nn <- 0.8
    } else if (R2_7 > R2_8 & R2_7 > R2_6) {
      R2 <- R2_7
      fit <- model7
      nn <- 0.7
    } else if (R2_6 > R2_8 & R2_6 > R2_7) {
      R2 <- R2_6
      fit <- model6
      nn <- 0.6
    }
  } else if (choice == 2) {
    if (R2_8 > R2_7) {
      R2 <- R2_8
      fit <- model8
      nn <- 0.8
    } else {
      R2 <- R2_7
      fit <- model7
      nn <- 0.7
    }
  }

  nlc <- function(x) predict(fit, newdata = x)
  attr(nlc, "R2") <- R2
  attr(nlc, "nn") <- nn

  if (PLOT) {
    localfit_plot(windowMean, windowPara, outdir, sample,
                  para, spc_info, outname,
                  fit.region, nlc, R2, nn)
  }

  nlc
}

#' Select a natural cubic spline model with high goodness of fit \eqn{R^2}

#' @param x0,y0 Two numeric vectors of reference means and local variances in
#'   the sliding window
#' @param choice A numeric to specify how many natural cubic spline nonlinear
#'   models to fit to choose from. By default, 2 non linear models with
#'   \code{df} in \code{\link[splines]{ns}} would be fitted and the one with
#'   higher goodness of fit \eqn{R^2} would be served as the final nonlinear
#'   model.
#'
#' @return A list which contains the final nonlinear model, goodness of fit
#'   \eqn{R^2} and \code{df}.
NCS_df_selection <- function(x0, y0, choice) {
  model5 <- lm(y0 ~ ns(x0, df = 5))
  R2_5 <- getIndexes(predict(model5, data.frame(x0 = x0)), y0)

  model4 <- lm(y0 ~ ns(x0, df = 4))
  R2_4 <- getIndexes(predict(model4, data.frame(x0 = x0)), y0)

  model3 <- lm(y0 ~ ns(x0, df = 3))
  R2_3 <- getIndexes(predict(model3, data.frame(x0 = x0)), y0)

  if (choice == 3) {
    if (R2_5 > R2_4 & R2_5 > R2_3) {
      R2 <- R2_5
      fit <- model5
      df <- 5
    } else if (R2_4 > R2_5 & R2_4 > R2_3) {
      R2 <- R2_4
      fit <- model4
      df <- 4
    } else if (R2_3 > R2_5 & R2_3 > R2_4) {
      R2 <- R2_3
      fit <- model3
      df <- 3
    }
  } else if (choice == 2) {
    if (R2_4 > R2_3) {
      R2 <- R2_4
      fit <- model4
      df <- 4
    } else {
      R2 <- R2_3
      fit <- model3
      df <- 3
    }
  }

  list(R2 = R2, fit = fit, df = df)
}
