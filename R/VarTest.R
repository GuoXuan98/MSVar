#' Identify hypervariable proteins
#'
#' Given a \code{\link{proObj}}, this function identifies hypervariable proteins
#' (HVPs) by performing a statistical test for each protein, with the resulting
#' \emph{p}-value assessing whether the observed \emph{log}-transformed variance
#' \emph{v} is significantly larger than a pre-specified constant.
#'
#' Proteins which are located in the lower part of \emph{v} are obviously
#' non-HVP, their \emph{p}-values are set as \code{NA} directly.
#'
#' @param x A \code{\link{proObj}} with estimated \emph{v} and fisher
#'   information.
#' @param null.v0 A pre-specified constant \eqn{v_0} used for HVP
#'   identification. The default setting of \eqn{v_0} amounts to an average
#'   fluctuation of 20\% around mean expression, that is, the standard deviation
#'   at the scale of M-values is \eqn{log_2 (1+20\%)} and \eqn{v_0 =
#'   log(log_2(1+20\%)^2)}.
#' @param p.adjust.method The method for adjusting \emph{p}-values for multiple
#'   testing. See \code{\link[stats]{p.adjust}} for all options.
#'
#' @return An object of the \code{\link[base]{class}} \code{c("varTest",
#'   "data.frame")}, which records the hyper variability test results for each
#'   region by each row.
#'
#' @export
#'
VarTest <- function(x, null.v0 = 0.2,
                    p.adjust.method = "BH") {
  if (!is(x, "proObj")) {
    stop("x must be of the class \"proObj\".")
  }
  if (is.null(x$resBio$est.bio.res)) {
    stop(paste("The biological variance have not been estimated yet.",
               "Try applying first the biological-variance-fitting process.",
               sep = "\n"
    ))
  }

  theo.v <- log(log2(1+null.v0)^2)

  # Derive p-values.
  est.bio.res <- x$resBio$est.bio.res
  v <- est.bio.res$v
  fisher.v <- est.bio.res$fisher.v
  Sta <- (v - theo.v) * sqrt(fisher.v)
  names(Sta) <- x$pro.name

  group1 <- attr(est.bio.res, "group1")
  # group2 <- setdiff(1:nrow(est.bio.res), group1)

  highv <- x$pro.name[group1]
  lowv <- setdiff(x$pro.name, highv)

  pval <- pnorm(Sta, 0, 1, lower.tail = F)
  pval[lowv] <- NA
  padj <- p.adjust(pval, method = p.adjust.method)

  # Return.
  res <- data.frame(Sta = Sta, v = v,
                    fisher.v = fisher.v,
                    pval = pval, padj = padj,
                    row.names = x$pro.name)
  attr(res, "null.v0") <- null.v0
  attr(res, "theo.v") <- theo.v
  attr(res, "p.adjust.method") <- p.adjust.method

  class(res) <- c("varTest", "data.frame")
  res
}

#' Plot a \code{varTest} Object
#'
#' Given a \code{varTest} object, this function creates a scatter plot of
#' \code{(ref.mean, v)} pairs from all observations, marking specifically the
#' ones that are associated with statistical significance. Besides, the
#' theoretical \eqn{v_0} associated with this round of test of \emph{v} is added
#' to the plot, serving as theoretical average fluctuation of \emph{v}.
#'
#' @param res A \code{varTest} object.
#' @param ref.mean A numeric vector of mean among all reference samples.
#' @param P.cutoff A scalar between 0 and 1. Proteins with \emph{p}-values below
#'   this cutoff are considered statistically significant under the null
#'   hypothesis and are correspondingly highlighted in a distinct color in the
#'   plot.
#' @param Nsample The sample size in this group.
#' @param sample.type The sample type of this group of samples.
#' @param args.points A named list of graphical parameters for plotting the
#'   point estimates (see \code{\link[graphics]{points}} for commonly used
#'   graphical parameters).
#' @param args.abline A named list of graphical parameters for plotting the null
#'   \eqn{v_0} (see \code{\link[graphics]{abline}} for commonly used graphical
#'   parameters).
#' @param args.legend A named list of graphical parameters for legend (see
#'   \code{\link[graphics]{legend}} for commonly used graphical parameters).
#' @param xlab Label for the X axis.
#' @param ylab Label for the Y axis.
#' @param xlim Range of X coordinate.
#' @param ylim Range of Y coordinate.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}}.
#'
#' @return This function returns \code{NULL} invisibly.
#' @export
#'
#' @seealso \code{\link{proObj}} for creating a \code{proObj} object;
#'   \code{\link{estTechVar}} for fitting technical variance given a
#'   \code{proObj} objects; \code{\link{estBioVar}} for estimating biological
#'   mean \emph{m} and \emph{log}-transformed biological variance \emph{v};
#'   \code{varTest} for a test of hypervariability.
varTestPlot <- function(res, ref.mean,
                        P.cutoff = 0.05,
                        Nsample, sample.type,
                        args.points = list(pch = 21, cex = 1.2),
                        args.abline = list(col = "darkorange", lwd = 2.5, lty = 5),
                        args.legend = list(x = "topright"),
                        xlab = "mean", ylab = "v = log(Var)",
                        xlim = NULL, ylim = NULL, ...) {
  if (!is(res, "varTest")) {
    stop("res must be of the class \"varTest\".")
  }

  null.v0 <- attr(res, "null.v0")

  v <- res$v
  pval <- res$pval
  ref.mean <- ref.mean[rownames(res)]
  theo.v <- log(log2(1+null.v0)^2)

  pval[which(is.na(pval))] <- 1

  main <- substitute(expression(atop("Mean-LogBioVar," ~ "H0: raw FC = " ~ null.v0,
                                     Npro ~ "proteins," ~ Nsample ~ sample.type ~ "samples")),
                     list(null.v0 = 1+null.v0,
                          Nsample = Nsample, sample.type = sample.type,
                          Npro = length(ref.mean)))
  if (is.null(ylim)) {
    ylim <- c(max(min(v), -22) - 0.2, max(v) + 0.2)
  }
  if (is.null(xlim)) {
    xlim <- c(min(ref.mean) - 0.01, max(ref.mean) + 0.01)
  }

  non_sig <- which(pval >= P.cutoff)
  sig <- which(pval < P.cutoff)

  args.plot <- list(main = main, xlab = xlab, ylab = ylab,
                    xlim = xlim, ylim = ylim,
                    cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4, ...)

  # Plot the observed mean-variance pairs.
  do.call(plot, c(list(ref.mean[non_sig], v[non_sig], type = "p"),
                  col = alpha("blue", 0.2),
                  args.points, args.plot))
  do.call(abline, c(list(h = theo.v), args.abline))
  do.call(points, c(list(ref.mean[sig], v[sig]),
                    col = alpha("red", 0.2), args.points))

  # Add the legend.
  temp <- list(legend = c("mean v scatter",
                          paste0("theoratical v0 of the test"),
                          paste0("sig protein, P < ", P.cutoff, ", ", length(sig))))

  args.legend <- c(args.legend, list(pch = c(21, NA, 21), lty = c(NA, 5, NA),
                                     col = c(alpha("blue", 0.2), "darkorange", alpha("red", 0.2)),
                                     pt.cex = 0.8, cex = 0.8))
  do.call(legend, c(temp, args.legend))

  invisible()
}
