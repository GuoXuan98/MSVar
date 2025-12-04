# Newton iteration --------------------------------------------------------

# Log-likelihood function
#' Gradient of log-likelihood function
#'
#' @param x A numeric vector contains 2 parameters to be estimated.
#' @param mj A numeric vector of M-value for a protein.
#' @param tj A numeric vector of estimated technical variance for a protein.
#'
#' @return The gradient of log-likelihood function.
#'
logL_gr <- function(x, mj, tj) {
  TVj <- exp(x[2]) + tj
  der <- numeric(2)
  der[1] <- -sum((mj - x[1]) / TVj)
  der[2] <- -exp(x[2]) * (-sum(1 / TVj) + sum(((mj - x[1]) / TVj)^2)) / 2
  der
}

#' Hessian matrix of log-likelihood function
#'
#' @param x A vector contains 2 parameters to be estimated.
#' @param mj A vector of M-value for a protein.
#' @param tj A vector of estimated technical variance for a protein.
#'
#' @return The Hessian matrix of log-likelihood function.
#'
logL_hessian <- function(x, mj, tj) {
  TVj <- exp(x[2]) + tj
  H <- matrix(0, 2, 2)
  H[1, 1] <- sum(1 / TVj)
  H[1, 2] <- H[2, 1] <- exp(x[2]) * sum((mj - x[1]) / TVj^2)
  H[2, 2] <- -(exp(x[2]) * (-sum(1 / TVj) + sum(((mj - x[1]) / TVj)^2)) / 2 + (exp(x[2])^2) * (sum(1 / TVj^2) / 2 - sum((mj - x[1])^2 / TVj^3)))
  H
}

#' Log-likelihood function
#'
#' @param x A vector contains 2 parameters to be estimated.
#' @param mj A vector of M-value for a protein.
#' @param tj A vector of estimated technical variance for a protein.
#'
#' @return The value of log-likelihood function.
logL_fun <- function(x, mj, tj) {
  TVj <- exp(x[2]) + tj
  logL <- -(length(mj) * log(2 * pi) + sum(log(TVj)) + sum((mj - x[1])^2 / TVj)) / 2
  f <- -logL
  attr(f, "gradient") <- logL_gr(x, mj, tj)
  attr(f, "hessian") <- logL_hessian(x, mj, tj)

  f
}

#' Nonlinear minimization of log-likelihood function
#'
#' Use Newton-iteration to minimize the log-likelihood function for the
#' parameter of biological variance and mean.
#'
#' @param MValue A vector of M-value for a protein.
#' @param EstTechVar A vector of estimated technical variance for a protein.
#'
#' @return A list contains the maximum likelihood estimates for
#'   \eqn{log}-transformed biological variance and biological mean, along with
#'   the corresponding minimized value of the log-likelihood function.
#' @export
#'
Newton_optimize <- function(MValue, EstTechVar) {
  nonNA <- which(!is.na(MValue) & !is.na(EstTechVar))
  MValue <- MValue[nonNA]
  EstTechVar <- EstTechVar[nonNA]

  m0 <- mean(MValue)
  v0 <- log(var(MValue))
  if (is.na(v0)) {
    v0 <- 0
  }

  res <- nlm(f = logL_fun,
             p = c(m0, v0),
             mj = MValue,
             tj = EstTechVar,
             hessian = F)

  jac <- logL_gr(res$estimate, mj = MValue, tj = EstTechVar)
  res$jac <- jac
  res
}


#' Calculate Fisher information of \emph{m} and \emph{v}
#'
#' @param bio.var A vector of biological variance for all proteins.
#' @param tech.var A matrix of technical variance for all proteins in all
#'   samples.
#'
#' @return A matrix which records the Fisher information of biological mean
#'   \emph{m}, \eqn{log}-transformed biological variance \emph{v} and
#'   biological variance.
Fisher_Info <- function(bio.var, tech.var) {
  TV <- bio.var + tech.var
  fisher.m <- rowSums(1 / TV, na.rm = T)
  fisher.bio.var <- rowSums(1 / TV^2, na.rm = T) / 2
  fisher.v <- bio.var^2 * fisher.bio.var

  data.frame(fisher.m = fisher.m,
             fisher.bio.var = fisher.bio.var,
             fisher.v = fisher.v)
}

# Exchange K-means cluster result to ensure proteins with high v are in cluster 1.

#' Divide proteins into 2 classes according to \code{v}
#'
#' Carry out K-means clustering on \code{v} and divide them into 2 groups.
#' Besides, check the labels of the observations to make sure the cluster with
#' higher cluster center is class 1 and the lower is class 2.
#'
#' @param v A numeric vector of all proteins' estimated \eqn{log}-biological
#'   variance \emph{v}.
#' @param seed A single value, used to ensure the reproducibility of K-means
#'   clustering result. Passed to \code{\link[base]{set.seed}}
#'
#' @return Labels of two classes of \code{v}.
kmean_clus <- function(v, seed) {
  set.seed(seed)

  kmeans_res <- kmeans(v, centers = 2)

  center <- kmeans_res$centers

  if (center[1] >= center[2]) {
    clus <- kmeans_res$cluster
  } else {
    clus_temp <- kmeans_res$cluster
    clus <- numeric(length = length(clus_temp))
    clus[which(clus_temp == 2)] <- 1
    clus[which(clus_temp == 1)] <- 2
  }
  clus
}

#' Combine the result of biological variance estimation
#'
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. This parameter determines whether to set a lower
#'   bound of \emph{v} and the specific setting of it.
#' @param ref.mean A numeric vector of mean among all reference samples.
#' @param m A numeric vector of biological mean estimates \emph{m} for all
#'   proteins.
#' @param v A numeric vector of \eqn{log}-transformed biological variance
#'   estimates \emph{v} for all proteins.
#' @param tech.var A numeric matrix of technical variance of all proteins in all
#'   samples. The rows are proteins and the columns are samples.
#' @param Nsample The sample size in this group.
#' @param sample.type The sample type of this group of samples.
#' @param pro.name A character vector of all proteins' names.
#' @param PLOT Whether to draw the process plots?
#' @param outdir If \code{PLOT} is TRUE, the output directory.
#' @param outadd If \code{PLOT} is TRUE, the additions to the output documents.
#'
#' @return A combined matrix with estimated \emph{m} & \emph{v} value and their
#'   fisher information matrix.
#' @export
#'
resBioCMB <- function(scaleV.raw.FC = NULL,
                      ref.mean, m, v, tech.var,
                      Nsample, sample.type, pro.name,
                      PLOT, outdir, outadd) {
  if (!is.null(scaleV.raw.FC)) {
    lower.bound <- log(log2(1 + scaleV.raw.FC)^2)

    group1 <- which(v > lower.bound)
    group2 <- which(v < lower.bound)

    scale.v <- v
    scale.v[which(v < lower.bound)] <- lower.bound
    bio.var <- exp(scale.v)

    outdir <- paste0(outdir, "/rawFCscale", scaleV.raw.FC)
    outadd <- paste0(outadd, scaleV.raw.FC)
  } else {
    cluster.v <- kmean_clus(v, 110)
    bio.var <- exp(v)
    scale.v <- v

    group1 <- which(cluster.v == 1)
    group2 <- which(cluster.v == 2)
  }

  fisher.info <- Fisher_Info(bio.var, tech.var)

  res.cmb <- data.frame(m = m,
                        v = scale.v,
                        bio.var = bio.var,
                        fisher.m = fisher.info$fisher.m,
                        fisher.v = fisher.info$fisher.v,
                        row.names = pro.name)

  if (!is.null(scaleV.raw.FC)) {
    colnames(res.cmb)[2:5] <- paste0(colnames(res.cmb)[2:5], ".scale")
    attr(res.cmb, "scaleV.raw.FC") <- scaleV.raw.FC
  }

  if (PLOT) {
    mkdir(outdir)
    FisherInfoPlot(m = m, v = scale.v,
                   fisher.m = fisher.info$fisher.m,
                   fisher.v = fisher.info$fisher.v,
                   Nsample = Nsample, sample.type = sample.type,
                   group1 = group1,
                   group2 = group2,
                   scaleV.raw.FC = scaleV.raw.FC,
                   outdir = outdir,
                   outadd = outadd)

    MeanVarPlot(ref.mean = ref.mean, v = v,
                Nsample = Nsample, sample.type = sample.type,
                group1 = group1,
                scaleV.raw.FC = scaleV.raw.FC,
                outdir = outdir,
                outadd = outadd)
  }

  attr(res.cmb, "group1") <- group1
  attr(res.cmb, "group2") <- group2

  res.cmb
}
