#' Estimate biological variance and mean
#'
#' Given a \code{\link{proObj}}, \code{estBioVar} estimate biological variance
#' and mean for each proteins among all given samples. To estimate these
#' parameters, maximum likelihood estimation method is adopted and Newton
#' iteration to solve the likelihood function.
#'
#' \code{estTechVar} well estimated the technical effects, giving
#' \code{adj.MValue} (if \code{correct.M} is TRUE) by normalization bias and
#' technical variance \code{est.tech.var} in \code{resTechGlobal} field.
#'
#' Then, the technical variances of M-value can be regarded as known parameters
#' in the log-likelihood function. For each protein, its log-likelihood function
#' is with 2 unknown parameters \emph{m} and \emph{v}. \code{estBioVar} would
#' take Newton iteration with the gradient matrix and Hessian matrix by the
#' expression of the function.
#'
#' After the estimation of all proteins, \code{resBioCMB} returns a combined
#' final results matrix with biological elements and their Fisher information.
#'
#' MSVar has set a lower bound for \emph{v}, denoted by \eqn{v_{lb}}. The major
#' reason is that the maximum likelihood estimate of log-biological variance
#' becomes rather unreliable when the underlying true biological variance is
#' close to 0. The large slope of log-transformation near 0 leads to extremely
#' small curvature of likelihood function. Consequently, the amount of Fisher
#' information associated with \emph{v} is extremely small, which leads to a
#' rather large variance of the estimator and can, thus, compromise severely the
#' statistical power for identifying differentially variable proteins. On this
#' account, MSVar replaces \emph{v} with \eqn{v_{lb}} if the former is smaller
#' than the latter, and the associated Fisher information is updated
#' accordingly.
#'
#' Therefore, \code{estBioVar} also include scaled results by a specified
#' \code{scaleV.raw.FC}.
#'
#' @param x A \code{\link{proObj}} with \code{resTech} field.
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. The default setting of \code{scaleV.raw.FC} is
#'   10\%, which amounts to an average fluctuation of 10\% around mean expression.
#'   Accordingly, the baseline fluctuation for each protein is defined as
#'   \eqn{log(log_2(1+10\%)^2)}.
#' @param Nsample The sample number in this group.
#' @param sample.type The sample type of this group of samples. If \code{sample.type}
#'   is NULL, the \code{name} of the \code{\link{proObj}} serves as the default
#'   replacement.
#' @param cpus The number of CPUs requested for a parallel computation. By
#'   default, no parallel computation is applied.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is TRUE, the output directory.
#'
#' @return A \code{\link{proObj}} with estimated biological variance and mean
#'   for each protein. \code{estBioVar} returns a \code{\link{proObj}} object,
#'   it has an added (updated) \code{resBio} field. Again, \code{resBio} field
#'   is a list consist of \code{est.bio.res} field and \code{est.bio.scale.res}
#'   optionally. \code{resBio} field is a matrix with 5 columns. It includes the
#'   estimated biological mean \emph{m}, estimated biological variances
#'   \code{bio.var} and their \eqn{log}-transformation \emph{v}, and the Fisher
#'   information \code{fisher.m} and \code{fisher.v} for \emph{m} and \emph{v}
#'   respectively.
#' @export
#'
estBioVar <- function(x,
                      scaleV.raw.FC = 0.1,
                      Nsample, sample.type = NULL, cpus = 1,
                      PLOT = FALSE,
                      outdir = NULL) {
  if (!is(x, "proObj")) {
    stop("x must be of the class \"proObj\".")
  }
  if (is.null(x$resTech$resTechGlobal$est.tech.var)) {
    stop(paste("No local fitted technical variance is found in x.",
      "Try calling estTechVar() to estimate one.",
      sep = "\n"
    ))
  }
  if (!is.null(x$resBio)) {
    x$resBio <- NULL
  }
  if (is.null(sample.type)) {
    sample.type <- x$name
  }

  if (!is.null(x$adj.MValue)) {
    MValue <- x$adj.MValue
  } else {
    MValue <- x$MValue
  }
  if (is.null(outdir)) {
    outdir <- "./est_bio_res/"
  }

  est.tech.var <- x$resTech$resTechGlobal$est.tech.var
  ref.mean <- x$ref.mean

  # Apply Newton iteration.
  cpus <- round(as.numeric(cpus)[1])
  if (cpus > 1) {
    sfInit(parallel = TRUE, cpus = cpus)
    res <- sfLapply(1:nrow(MValue), function(i) {
      Newton_optimize(unlist(MValue[i, ]), unlist(est.tech.var[i, ]))
    })
    sfStop()
  } else {
    res <- lapply(1:nrow(MValue), function(i) {
      Newton_optimize(unlist(MValue[i, ]), unlist(est.tech.var[i, ]))
    })
  }

  est.fun <- vapply(1:length(res), function(i) {
    res[[i]][["minimum"]]
  }, numeric(1))
  est.v <- vapply(1:length(res), function(i) {
    res[[i]][["estimate"]][2]
  }, numeric(1))
  est.m <- vapply(1:length(res), function(i) {
    res[[i]][["estimate"]][1]
  }, numeric(1))

  # Combine the results of Newton iteration without scaling of v.
  est.res <- resBioCMB(ref.mean = ref.mean, m = est.m, v = est.v,
                       tech.var = est.tech.var,
                       Nsample = Nsample, sample.type = sample.type,
                       pro.name = x$pro.name,
                       PLOT = PLOT,
                       outdir = outdir, outadd = paste0("rawv"))

  # Scale v by scaleV.raw.FC and combine the updated results.
  if (!is.null(scaleV.raw.FC)) {
    est.scale.res <- lapply(scaleV.raw.FC, resBioCMB,
                            ref.mean = ref.mean, m = est.m, v = est.v,
                            tech.var = est.tech.var,
                            Nsample = Nsample, sample.type = sample.type,
                            pro.name = x$pro.name,
                            PLOT = PLOT,
                            outdir = paste0(outdir, "/est_v_rawFCscale_result/"),
                            outadd = paste0("rawFCscale"))

    names(est.scale.res) <- paste0("rawFCscale", scaleV.raw.FC*100, "%")
  }

  x$resBio <- list(est.bio.res = est.res)
  if (!is.null(scaleV.raw.FC)) {
    x$resBio$est.bio.scale.res <- est.scale.res
  }

  x
}
