# Normalize the raw intensity.

# Normalization -----------------------------------------------------------

#' Normalize the intensity for each sample against the reference sample
#'
#' Given the raw intensity for each sample and the corresponding reference
#' sample, conduct 1v1 normalization. Only use proteins which are non-NA in both
#' samples to calculate trimmed total intensity, and then obtain the
#' normalization factor for each sample. For all proteins in the sample, its raw
#' intensity would divide the normalization factor to get the normalized
#' intensity. However, the reference sample remain unchanged.
#'
#' After normalization, intensity in each sample and reference sample would take
#' \eqn{log_2}-transformation and \code{normalize} would calculate A-values and
#' M-values in \eqn{log_2}-scale. A-value is regarded as the basic intensity of
#' the observations, and M-value is \eqn{log_2}-transformed fold change.
#'
#' @param raw.intensity A raw intensity matrix which is composed by a sample and
#'   the corresponding reference sample.
#' @param filter.low Whether to filter proteins whose raw intensity is under 10?
#'   The default is TRUE.
#' @param offset A scalar used as an offset added to the raw intensity, avoiding
#'   infinity while \eqn{log_2}-transformation. Default is 1.
#' @param L A scalar to scale the cutoff of outliers' boundary.
#' @param A.method A character string indicating the method to calculate
#'   A-value. One of \code{"mean_of_2samples"}, \code{"ref"} or \code{"sample"}.
#'   Usually A-value is regarded as the basic intensity of the observations and
#'   it will take the intensity of reference.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is true, the directory of output figures.
#'
#' @return A normalized intensity 3-column matrix with normalized intensity,
#'   M-value and A-value in the \eqn{log_2}-scale.
#' @export
#'
normalize <- function(raw.intensity, filter.low = T, offset = 1,
                      L = 1.5, A.method = c("mean_of_2samples", "ref", "sample"),
                      PLOT = FALSE, outdir) {
  A.method <- match.arg(A.method)
  offset <- as.numeric(offset[1])

  log2.int <- log(raw.intensity + offset, base = 2)

  new.raw.intensity <- raw.intensity
  new.raw.intensity <- new.raw.intensity[which(new.raw.intensity[, 1] != 0 & new.raw.intensity[, 2] != 0), ]
  new.raw.intensity <- new.raw.intensity[which(!is.na(new.raw.intensity[, 1]) & !is.na(new.raw.intensity[, 2])), ]

  # Drop proteins which have intensity lower than 10 among all samples.
  good.index <- numeric(0)
  if (filter.low) {
    for (i in 1:nrow(new.raw.intensity)) {
      index.df <- new.raw.intensity[i, ] > 10

      num <- sum(index.df)
      if (num == ncol(new.raw.intensity)) {
        good.index <- c(good.index, i)
      } else {
        next
      }
    }
  } else {
    good.index <- 1:nrow(new.raw.intensity)
  }

  new.raw.intensity <- new.raw.intensity[good.index, ]

  Ncol <- ncol(new.raw.intensity)
  L <- L

  sample <- colnames(new.raw.intensity)[1]
  ref <- colnames(new.raw.intensity)[2]
  sample_outdir <- paste0(outdir, "/", sample, "/")

  # Begin normalization.
  trimmed.intensity <- rep(0, Ncol)
  NonOutlier <- matrix(0, nrow = nrow(new.raw.intensity), ncol = Ncol)
  for (j in 1:Ncol)
  {
    Q1 <- quantile(new.raw.intensity[, j], 0.25)
    Q3 <- quantile(new.raw.intensity[, j], 0.75)

    cut <- Q3 + L * (Q3 - Q1)

    # Non-outlier in sample j.
    NonOutlier[, j] <- new.raw.intensity[, j] < cut
  }

  # Non-outlier proteins in 2 samples.
  NonOutlier.all <- which(apply(NonOutlier, 1, sum) == Ncol)

  trimmed.intensity <- apply(new.raw.intensity[NonOutlier.all, ], 2, sum)

  norm.factor <- rep(1, 2)
  norm.factor[1] <- trimmed.intensity[1] / mean(trimmed.intensity)

  norm.int <- as.data.frame(lapply(1:Ncol, function(j) {
    raw.intensity[, j] / norm.factor[j]
  }))
  rownames(norm.int) <- rownames(raw.intensity)
  colnames(norm.int) <- colnames(raw.intensity)

  log2.norm.int <- log(norm.int + offset, base = 2)

  if (PLOT) {
    mkdir(sample_outdir)
    MA_before_after_norm_plot(log2.int = log2.int,
                              log2.norm.int = log2.norm.int,
                              sample_outdir = sample_outdir, sample = sample, ref = ref,
                              A.method = A.method)
  }

  log2.norm.int["MValue"] <- log2.norm.int[, sample] - log2.norm.int[, ref]

  log2.norm.int["AValue"] <- AVal(x.sample = log2.norm.int[, sample],
                                  y.ref = log2.norm.int[, ref],
                                  A.method = A.method)

  log2.norm.int <- log2.norm.int[, c(sample, "MValue", "AValue")]

  attr(log2.norm.int, "norm.factor") <- norm.factor
  log2.norm.int
}

#' Calculate A-value
#'
#' @param x.sample A vector of normalized \eqn{log_2} intensity of tissue
#'   sample.
#' @param y.ref A vector of normalized \eqn{log_2} intensity of reference
#'   sample.
#' @param A.method A character string indicating the method to calculate
#'   A-value. One of \code{"mean_of_2samples"}, \code{"ref"} or \code{"sample"}.
#'   If it is \code{"mean_of_2samples"}, A-value is the mean of tissue sample
#'   and reference sample; else if it is \code{"ref"}, A-value take as the
#'   intensity of reference sample; else if it is \code{"sample"}, A-value take
#'   as the intensity of tissue sample. \code{"ref"} is recommended.
#'
#' @return A vector of A-value.
AVal <- function(x.sample, y.ref, A.method) {
  if (A.method == "mean_of_2samples") {
    AVal <- (x.sample + y.ref) / 2
  } else if (A.method == "ref") {
    AVal <- y.ref
  } else if (A.method == "sample") {
    AVal <- x.sample
  }
}

# Plot --------------------------------------------------------------------

#' Draw the M-A plot in \eqn{log_2}-scale before and after normalization
#'
#' \code{MA_before_after_norm_plot} would draw 4 subplots. Subplot1 and subplot3
#' would draw \eqn{log_2}-intensity of the tissue sample against the reference
#' sample before and after normalization. Subplot2 and subplot4 would draw
#' M-value against A-value before and after normalization.
#'
#' @param log2.int A matrix with 2 columns. The first column is the \eqn{log_2}
#'   intensity of current sample, and the second column is the \eqn{log_2}
#'   intensity of its corresponding reference sample. They are not normalized.
#' @param log2.norm.int A matrix with 2 columns. The first column is the
#'   \eqn{log_2}-intensity of current sample, and the second column is the
#'   \eqn{log_2}-intensity of its corresponding reference sample. They are
#'   normalized.
#' @param sample_outdir A character which represent the output directory.
#' @param sample The name of current sample to be fitted.
#' @param ref The name of reference sample related to current sample.
#' @param A.method A character string indicating the method to calculate
#'   A-value. One of \code{"mean_of_2samples"}, \code{"ref"} or \code{"sample"}.
#'   Usually A-value is regarded as the basic intensity of the observations and
#'   it will take the intensity of reference.
#'
#' @return This function returns \code{NULL} invisibly.
MA_before_after_norm_plot <- function(log2.int,
                                      log2.norm.int,
                                      sample_outdir, sample, ref,
                                      A.method) {
  log2.int <- log2.int[which(!is.na(log2.int[, 1]) & !is.na(log2.int[, 2])), ]
  log2.norm.int <- log2.norm.int[which(!is.na(log2.norm.int[, 1]) & !is.na(log2.norm.int[, 2])), ]

  # Plot before normalization.
  png(paste0(sample_outdir, "intensity_MA_scatter.png"),
      width = 540 * 4, height = 4 * 540, res = 72 * 3)
  dev.control("enable")

  par(mfcol = c(2, 2),
    mar = c(6, 6, 5, 3))

  # subplot1
  R_square <- round(getIndexes(log2.int[, sample], log2.int[, ref]), 2)

  min_value <- min(min(log2.int[, sample]), min(log2.int[, ref]))
  max_value <- max(max(log2.int[, sample]), max(log2.int[, ref]))

  plot(log2.int[, sample], log2.int[, ref],
       pch = 20, col = alpha("#1f77b4", 0.2), family = "sans",
       xlab = sample, ylab = "ref", cex.main = 1.5, cex.lab = 1.5,
       xlim = c(min_value - 1, max_value + 1), ylim = c(min_value - 1, max_value + 1),
       main = substitute(atop("Raw intensity", R^2 ~ " = " ~ R_square)))
  abline(a = 0, b = 1, col = "black", lwd = 1.5, lty = 2)

  # subplot2
  plot(AVal(log2.int[, sample], log2.int[, ref], A.method),
       log2.int[, sample] - log2.int[, ref],
       pch = 20, col = alpha("#1f77b4", 0.2), family = "sans",
       main = "M-A plot of raw intensity", cex.main = 1.5, cex.lab = 1.5,
       xlab = substitute(atop("", log[2] ~ "intensity," ~ A.method),
                         list(A.method = A.method)),
       ylab = expression(paste(log[2], " ratio")),
       xlim = c(min_value - 1, max_value + 1))
  abline(h = 0, col = "black", lwd = 1.5, lty = 2)

  # subplot3
  R_square <- round(getIndexes(log2.norm.int[, sample], log2.norm.int[, ref]), 2)

  min_value <- min(min(log2.norm.int[, sample]), min(log2.norm.int[, ref]))
  max_value <- max(max(log2.norm.int[, sample]), max(log2.norm.int[, ref]))

  plot(log2.norm.int[, sample], log2.norm.int[, ref],
       pch = 20, col = alpha("#1f77b4", 0.2), family = "sans",
       xlab = sample, ylab = "ref", cex.main = 1.5, cex.lab = 1.5,
       xlim = c(min_value - 1, max_value + 1), ylim = c(min_value - 1, max_value + 1),
       main = substitute(atop("norm intensity", R^2 ~ " = " ~ R_square)))
  abline(a = 0, b = 1, col = "black", lwd = 1.5, lty = 2)

  # subplot4
  plot(AVal(log2.norm.int[, sample], log2.norm.int[, ref], A.method),
       log2.norm.int[, sample] - log2.norm.int[, ref],
       pch = 20, col = alpha("#1f77b4", 0.2), family = "sans",
       main = "M-A plot of norm intensity", cex.main = 1.5, cex.lab = 1.5,
       xlab = substitute(atop("", log[2] ~ "intensity," ~ A.method),
                         list(A.method = A.method)),
       ylab = expression(paste(log[2], " ratio")),
       xlim = c(min_value - 1, max_value + 1))
  abline(h = 0, col = "black", lwd = 1.5, lty = 2)

  dev.off()

  invisible()
}


# Create directory ------------------------------------------------------

#' Create new folders
#'
#' Create new folders according to the directory provided by users to restore
#' the process plots during the estimation.
#'
#' @param path_dir A character declares the directory to create new folder.
#'
#' @return A newly created folder.
#' @export
#'
mkdir <- function(path_dir) {
  isexists <- dir.exists(path_dir)
  if (!isexists) {
    dir.create(path_dir, recursive = T)
  }
}
