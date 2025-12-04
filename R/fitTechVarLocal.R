# Estimate technical variance locally.

#' Locally fit the technical variance and bias for each sample
#'
#' Given a \code{\link{proObj}}, this function locally fit the linear model for
#' technical variance and bias in each sample by linear fitting.
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
#' @param cpus The number of CPUs requested for a parallel computation. By
#'   default, no parallel computation is applied.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return \code{estTechVarLocal} returns a \code{\link{proObj}} object, it has
#'   has an added (updated) \code{resTechLocal} in \code{resTech} field
#'   describing locally estimated technical variance and normalization bias, a
#'   \code{MValue} field recording M-value of all samples and a \code{ref.mean}
#'   field indicating mean of all reference samples. The \code{resTechLocal}
#'   field is itself a list of all samples, and each one is consist of the
#'   following components:
#'     \describe{
#'         \item{\code{log2.norm.int}}{A normalized intensity 3-column matrix
#'         with normalized intensity, M-value and A-value in the \eqn{log_2}
#'         scale.}
#'         \item{\code{MA.linear.fit}}{A 4-column matrix records all of the
#'         results about linear fitting in each sliding window of current
#'         sample, which includes fitted variances, bias, mean of A-value in
#'         each window and goodness of fit \eqn{R^2}.}
#'         \item{\code{name}}{The name of current sample.}
#'     }
#'   \code{MValue} is a matrix whose number of rows is the same as number of all
#'   proteins and number of columns is the same as sample size. It record
#'   M-value (i.e. \eqn{log_2}-transformed fold change) of all samples.
#'
#'   \code{AValue} is a matrix whose number of rows is the same as \code{MValue}
#'   and the number of columns is the same as the number of reference samples.
#'
#'   \code{ref.mean} is a vector whose length is the same as \code{MValue}. It
#'   indicates the mean of normalized \eqn{log_2}-intensity among all reference
#'   samples (i.e. mean of A-value among all reference samples), which is
#'   regarded as the basic intensity of the observations in the whole group of
#'   samples.
#'
#' @export
#'
estTechVarLocal <- function(x, A.method = "ref",
                            filter.low = T, offset = 1,
                            window.size = 400, step.size = 100,
                            middle.percent = 0.5, cpus = 1,
                            PLOT = FALSE, outdir = "./") {
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

  pro.name <- x$pro.name
  sample.intensity <- x$sample.intensity[pro.name, , drop = F]
  ref.intensity <- x$ref.intensity[pro.name, , drop = F]

  sample.ref.cor <- x$sample.ref.cor

  resTechLocal <- lapply(1:length(sample.ref.cor), function(i) {
    sample <- names(sample.ref.cor)[i]
    ref.sample <- sample.ref.cor[i]

    paired.intensity <- cbind(sample.intensity[, sample],
                              ref.intensity[, ref.sample])
    colnames(paired.intensity) <- c(sample, ref.sample)
    rownames(paired.intensity) <- pro.name

    resTechLocal_pair <- estTechVarLocal_pair(sample, ref.sample,
                                              paired.intensity,
                                              A.method, filter.low, offset,
                                              window.size, step.size,
                                              middle.percent, PLOT,
                                              cpus = cpus, outdir)
    resTechLocal_pair
  })
  names(resTechLocal) <- names(sample.ref.cor)

  MValue <- lapply(1:length(resTechLocal), function(i) {
    resTechLocal[[i]]$log2.norm.int[x$pro.name, "MValue"]
  })
  MValue <- do.call(cbind, MValue)

  AValue <- lapply(1:length(resTechLocal), function(i) {
    resTechLocal[[i]]$log2.norm.int[x$pro.name, "AValue"]
  })
  AValue <- do.call(cbind, AValue)

  rownames(MValue) <- rownames(AValue) <- x$pro.name
  colnames(MValue) <- colnames(AValue) <- names(resTechLocal)

  ref.mean <- apply(AValue, 1, mean, na.rm = T)

  x$resTech <- list(resTechLocal = resTechLocal)
  x$MValue <- MValue
  x$ref.mean <- ref.mean

  x
}


#' Normalize raw intensity and fit local parameters for each sample
#'
#' @param sample The name of current sample to be fitted.
#' @param ref.sample The name of reference sample related to current sample.
#' @param paired.int A matrix with 2 columns. The first column is the raw
#'   intensity of current sample, and the second column is the raw intensity of
#'   its corresponding reference sample.
#' @param A.method A character string indicating the method to calculate
#'   A-value. One of \code{"mean_of_2samples"}, \code{"ref"} (default) or
#'   \code{"sample"}. See \code{normalize} for more details.
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
#' @param cpus The number of CPUs requested for a parallel computation. By
#'   default, no parallel computation is applied.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the output directory.
#'
#' @return A \code{list} records the normalized \eqn{log2}-intensity,
estTechVarLocal_pair <- function(sample, ref.sample,
                                 paired.int,
                                 A.method = "ref",
                                 filter.low = T,
                                 offset = 1,
                                 window.size, step.size,
                                 middle.percent, cpus,
                                 PLOT, outdir) {
  # log2 intensity with normalization
  log2.norm.int <- normalize(raw.intensity = paired.int,
                             A.method = A.method,
                             filter.low = filter.low,
                             offset = offset,
                             PLOT = PLOT, outdir = outdir)

  MA.linear.fit <- estTechVarLocal_pairFit(log2.norm.int = log2.norm.int,
                                           sample = sample,
                                           window.size = window.size,
                                           step.size = step.size,
                                           middle.percent = middle.percent,
                                           PLOT = PLOT,
                                           cpus = cpus, outdir = outdir)

  list(log2.norm.int = log2.norm.int,
       MA.linear.fit = MA.linear.fit,
       name = sample)
}


#' Estimate local technical variance and normalization bias
#'
#' Given normalized \eqn{log2}-intensity and the corresponding A-value and
#' M-value, scan the M-A plot by sliding window with \code{window.size} as
#' window size and \code{step.size} as step size. In each sliding window, plot
#' M-value against the corresponding quantile of standard normal distribution
#' and use a linear fitting model to clarify their relationship. The slope of
#' the linear model is local technical standard deviation of this window and the
#' intercept is the normalization bias.
#'
#' @inheritParams estTechVarLocal
#' @param log2.norm.int A matrix with normalized \eqn{log2}-intensity and the
#'   corresponding A-value and M-value.
#' @param sample The name of current sample to be fitted.
#'
#' @return A matrix records all of the results about linear fitting in each
#'   sliding window, which includes fitted variances, bias, mean of A-value in
#'   each window and goodness of fit \eqn{R^2}.
estTechVarLocal_pairFit <- function(log2.norm.int, sample,
                                    window.size = 400, step.size = 100,
                                    middle.percent = 0.5,
                                    PLOT = FALSE, cpus = 1, outdir = "./") {
  AValue <- log2.norm.int$AValue
  MValue <- log2.norm.int$MValue

  Nrow <- length(MValue)

  # Index of unique non-NA A-value.
  UniqueA <- which(duplicated(AValue) == FALSE)
  if (length(which(is.na(AValue))) != 0) {
    UniqueA <- UniqueA[-which(is.na(AValue))[1]]
  }

  # Index of unique non-NA M-value.
  UniqueM <- which(duplicated(MValue) == FALSE)
  if (length(which(is.na(MValue))) != 0) {
    UniqueM <- UniqueM[-which(is.na(MValue))[1]]
  }

  UniqueAM <- intersect(UniqueA, UniqueM)
  A_unique <- AValue[UniqueAM]
  M_unique <- MValue[UniqueAM]

  A_unique_sort <- sort(A_unique)
  M_unique_sort <- M_unique[order(A_unique)]
  index <- 1:length(A_unique_sort)

  data_sort <- data.frame(sortA = A_unique_sort,
                          sortM = M_unique_sort,
                          index = 1:length(A_unique_sort),
                          row.names = rownames(log2.norm.int)[UniqueAM][order(A_unique)])

  Ndot <- length(A_unique_sort)

  window.mu <- numeric(0)
  window.var <- numeric(0)
  window.mean.int <- numeric(0)
  window.R2 <- numeric(0)

  ends <- seq(from = window.size, to = Ndot, by = step.size)
  groups <- lapply(ends, function(i) index[(i - window.size + 1):i])
  last <- ends[length(ends)]
  groups[[length(groups)]] <- index[(last - window.size + 1):length(index)]
  names(groups) <- paste0("window", 1:length(groups))

  if (length(groups) == 0) {
    stop("No valid groups of observations at all.")
  }

  # Apply linear fitting in each sliding window.
  cpus <- round(as.numeric(cpus)[1])
  if (cpus > 1) {
    sfInit(parallel = TRUE, cpus = cpus)
    res <- sfLapply(1:length(groups), MA_linear_regression,
                    data = data_sort, groups = groups, sample = sample,
                    middle.percent = middle.percent,
                    PLOT = PLOT, outdir = outdir)
    sfStop()
  } else {
    res <- lapply(1:length(groups), MA_linear_regression,
                  data = data_sort, groups = groups, sample = sample,
                  middle.percent = middle.percent,
                  PLOT = PLOT, outdir = outdir)
  }

  # Return.
  res <- as.data.frame(res)
  colnames(res) <- names(groups)
  rownames(res) <- c("window.bias", "window.var",
                     "window.mean.int", "window.R2")
  res <- as.data.frame(t(res))

  attr(res, "groups") <- groups

  if (!is.null(window.size)) {
    attr(res, "window.size") <- window.size
    attr(res, "step.size") <- step.size
  }

  res
}
