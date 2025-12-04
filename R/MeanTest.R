#' Identify differentially expressed proteins
#'
#' Conduct hypothesis test for differentially expressed proteins. All \code{\link{proObj}}
#' objects should be associated with the same scaling of estimated \emph{v}.
#' Thus, they must have the same scaling results of \emph{v} and \emph{m}.
#'
#' @param x,y \code{x} is any R object for which a \code{diffTest} method has
#'   been defined. For the method for class \code{\link{proObj}}, \code{x} and
#'   \code{y} are two \code{proObj} objects to be compared. They must be
#'   associated with estimated \emph{m}, Fisher information and the same scaled
#'   \emph{m}.
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. The default setting of \code{scaleV.raw.FC} is
#'   10\%, which amounts to an average fluctuation of 10\% around mean expression.
#'   Accordingly, the baseline fluctuation for each protein is defined as
#'   \eqn{log(log_2(1+10\%)^2)}. If \code{scaleV.raw.FC} is 0, the differential
#'   expression test would be conducted without scaling of estimated \emph{m}.
#' @param alternative The alternative hypothesis. Must be one of \code{"upper"}
#'   (default), \code{"lower"}, and \code{"two.sided"}. Can be abbreviated.
#' @param p.adjust.method The method for adjusting \emph{p}-values for multiple
#'   testing. See \code{\link[stats]{p.adjust}} for all options.
#'
#' @return An object of the \code{\link[base]{class}} \code{c("meanTest",
#'   "data.frame")}, which records the differential expression test results for
#'   each protein by each row.
#' @export
#'
DiffMean <- function(x, y, scaleV.raw.FC = 0.1,
                     alternative = c("upper", "lower", "two"),
                     p.adjust.method = "BH") {
  if (!is(x, "proObj") | !is(y, "proObj")) {
    stop("x and y must be of the class \"proObj\".")
  }

  if (scaleV.raw.FC == 0) {
    if (is.null(x$resBio$est.bio.res) | is.null(y$resBio$est.bio.res)) {
      stop(paste("The biological variance have not been estimated yet.",
                 "Try applying first the biological-variance-fitting process.",
                 sep = "\n"))
    } else {
      est.bio.res1 <- x$resBio$est.bio.res
      est.bio.res2 <- y$resBio$est.bio.res
    }
  } else {
    if (is.null(x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]) |
        is.null(y$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]])) {
      stop(paste("The biological variance have not been scaled yet.",
                 "Try applying first the biological-variance-scaling process.",
                 sep = "\n"))
    }
    est.bio.res1 <- x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]
    est.bio.res2 <- y$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]
  }

  pro_name <- intersect(rownames(est.bio.res1), rownames(est.bio.res2))
  est.bio.res1 <- est.bio.res1[pro_name, ]
  est.bio.res2 <- est.bio.res2[pro_name, ]

  # Derive p-values.
  m1 <- est.bio.res1$m
  fisher.m1 <- est.bio.res1$fisher.m
  m2 <- est.bio.res2$m
  fisher.m2 <- est.bio.res2$fisher.m
  Sta <- (m1 - m2) / sqrt(1 / fisher.m1 + 1 / fisher.m2)

  alternative <- match.arg(alternative)
  if (alternative == "upper") {
    pval <- pnorm(Sta, 0, 1, lower.tail = F)
  } else if (alternative == "lower") {
    pval <- pnorm(Sta, 0, 1, lower.tail = T)
  } else if (alternative == "two") {
    pval <- 2 * pnorm(-abs(Sta), 0, 1, lower.tail = T)
  }
  padj <- p.adjust(pval, method = p.adjust.method)

  # Return.
  res <- data.frame(Sta = Sta,
                    m1 = m1, fisher.m1 = fisher.m1,
                    m2 = m2, fisher.m2 = fisher.m2,
                    pval = pval, padj = padj,
                    row.names = pro_name)

  sample_type1 <- strsplit(x$name,"[.]")[[1]][2]
  sample_type1 <- sub("(.)", "\\U\\1", sample_type1, perl = TRUE)

  attr(res, "scaleV.raw.FC") <- scaleV.raw.FC
  attr(res, "p.adjust.method") <- p.adjust.method
  if(alternative == "two") {
    attr(res, "alternative") <- "two sided"
  } else if(alternative == "upper") {
    attr(res, "alternative") <- paste0(sample_type1, " High")
  } else if(alternative == "lower") {
    attr(res, "alternative") <- paste0(sample_type1, " Low")
  }

  class(res) <- c("meanTest", "data.frame")
  res
}

#' Plot a \code{meanTest} Object
#'
#' Given a \code{meanTest} object, which records the results of calling
#' HVPs and/or LVPs from a \code{\link{proObj}}, this method creates a scatter
#' plot of \code{(ref.mean, v)} pairs from all observations, marking
#' specifically the ones that are associated with statistical significance.
#' Besides, the theoretical \eqn{v_0} associated with this round of test of
#' \eqn{v} is added to the plot, serving as theoretical mean of \emph{v}.
#'
#' @param res A \code{meanTest} object.
#' @param ref.mean A numeric vector of mean among all reference samples.
#' @param P.cutoff A scalar between 0 and 1. Proteins with \emph{p}-values below
#'   this cutoff are considered statistically significant under the null
#'   hypothesis and are correspondingly highlighted in a distinct color in the
#'   plot.
#' @param p.type A character which indicates the type of \emph{p}-value to show
#'   in the plot. Must be one of \code{"pval"} (default), and \code{"padj"}. Can
#'   be abbreviated.
#' @param mainAdd The added information to describe samples used to do test.
#' @param raw_Diff Whether to use the raw diff of m as y.
#' @param args.points A named list of graphical parameters for plotting the
#'   point estimates (see \code{\link[graphics]{points}} for commonly used
#'   graphical parameters).
#' @param args.abline A named list of graphical parameters for plotting the null
#'   \emph{v_0} (see \code{\link[graphics]{abline}} for commonly used graphical
#'   parameters).
#' @param args.legend A named list of graphical parameters for legend (see
#'   \code{\link[graphics]{legend}} for commonly used graphical parameters).
#' @param xlab Label for the X axis.
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
#'   \code{varTest} for a test of hyper-variability.
meanTestPlot <- function(res, ref.mean,
                         P.cutoff = 0.05,
                         p.type = c("pval", "padj"),
                         mainAdd, raw_Diff = F,
                         args.points = list(pch = 21, cex = 1.2),
                         args.abline = list(col = "darkorange", lwd = 2.5, lty = 5),
                         args.legend = list(x = "topright"),
                         xlab = "mean",
                         xlim = NULL, ylim = NULL, ...) {
  if (!is(res, "meanTest")) {
    stop("res must be of the class \"meanTest\".")
  }

  scaleV.raw.FC <- attr(res, "scaleV.raw.FC")
  alternative <- attr(res, "alternative")

  m <- res$m
  Sta <- res$Sta
  ref.mean <- ref.mean[rownames(res)]

  if(raw_Diff) {
    ylab <- expression(m[1]-m[2])
    Sta <- res$m1 - res$m2
  } else {
    ylab <- "diff m Sta"
    Sta <- res$Sta
  }

  if(p.type == "pval") {
    PVal <- res$pval
    PLeg <- "P val"
  } else if (p.type == "padj") {
    PVal <- res$padj
    PLeg <- "BH P"
  }
  PVal[which(is.na(PVal))] <- 1

  main <- substitute(expression(atop("Mean-DiffM," ~ alternative * ", scale v by" ~ scaleV.raw.FC,
                                     Npro ~ "proteins," ~ mainAdd)),
                     list(alternative = alternative, scaleV.raw.FC = 1+scaleV.raw.FC,
                          Npro = length(ref.mean), mainAdd = mainAdd))

  if (is.null(ylim)) {
    ylim <- c(min(Sta, na.rm = T) - 0.2, max(Sta, na.rm = T) + 0.2)
  }
  if (is.null(xlim)) {
    xlim <- c(min(ref.mean) - 0.01, max(ref.mean) + 0.01)
  }

  non_sig <- which(PVal >= P.cutoff)
  sig <- which(PVal < P.cutoff)

  args.plot <- list(main = main, xlab = xlab, ylab = ylab,
                    xlim = xlim, ylim = ylim,
                    cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4, ...)

  # Plot the observed mean-variance pairs.
  do.call(plot, c(list(ref.mean[non_sig], Sta[non_sig], type = "p"),
                  col = alpha("blue", 0.2),
                  args.points, args.plot))
  do.call(abline, c(list(h = 0), args.abline))
  do.call(points, c(list(ref.mean[sig], Sta[sig]),
                    col = alpha("red", 0.2), args.points))

  # Add the legend.
  temp <- list(legend = c("mean Diff-m scatter",
                          paste0("theoratical 0"),
                          paste0("sig protein, ", PLeg, " < ", P.cutoff, ", ", length(sig))))

  args.legend <- c(args.legend, list(pch = c(21, NA, 21), lty = c(NA, 5, NA),
                                     col = c(alpha("blue", 0.2), "darkorange", alpha("red", 0.2)),
                                     pt.cex = 0.8, cex = 0.8))
  do.call(legend, c(temp, args.legend))

  invisible()
}
