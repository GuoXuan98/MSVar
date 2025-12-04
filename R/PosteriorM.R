#' Derive posterior M-value
#'
#' MSVar derives a posterior estimate of each \eqn{M'_ij} by finding the maximum
#' a posteriori.
#'
#' @param x A \code{\link{proObj}}.
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. The default setting of \code{scaleV.raw.FC} is
#'   10\%, which amounts to an average fluctuation of 10\% around mean expression.
#'   Accordingly, the baseline fluctuation for each protein is defined as
#'   \eqn{log(log_2(1+10\%)^2)}.
#'
#' @return A \code{\link{proObj}} with posterior M-value.
#' @export
#'
PostM <- function(x, scaleV.raw.FC = 0.1) {
  if (!is(x, "proObj")) {
    stop("x must be of the class \"proObj\".")
  }

  if (!is.null(x$adj.MValue)) {
    MValue <- x$adj.MValue
  } else {
    MValue <- x$MValue
  }

  if (scaleV.raw.FC != 0) {
    if (is.null(x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC*100, "%")]])) {
      stop(paste("No scaled biological results with this specified scaleV.raw.FC are found in x.",
                 "Try calling estBioVar() to initialize one.",
                 sep = "\n"))
    } else {
      bio.res <- x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC*100, "%")]]
      bio.var <- bio.res$bio.var.scale
    }
  } else {
    bio.res <- x$resBio$est.bio.res
    bio.var <- bio.res$bio.var
  }

  tech.var <- x$resTech$resTechGlobal$est.tech.var

  # Derive posterior M-value.
  MPost <- (bio.res$m * tech.var + MValue * bio.var) / (tech.var + bio.var)

  # Center and scale each row(protein) of the posterior M-value matrix.
  ZScrPost <- t(scale(t(MPost)))

  x$PostM <- MPost
  x$PostZScr <- ZScrPost

  x
}
