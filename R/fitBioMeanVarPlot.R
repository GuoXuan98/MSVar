#' Plot Fisher information of \emph{m} and \emph{v}
#'
#' \code{FisherInfoPlot} would draw 4 subplots. The first 2 subplots in the
#' first row is the fisher information of \emph{m}, and proteins are divided in
#' 2 groups by \code{group1} and \code{group2}. The other 2 subplots in the
#' second row is the fisher information of \emph{v}, groups as above.
#'
#' @param m A numeric vector of biological mean estimates \emph{m} for all
#'   proteins.
#' @param v A numeric vector of \eqn{log}-transformed biological variance
#'   estimates \emph{v} for all proteins.
#' @param fisher.m A numeric vector of Fisher information of \emph{m}.
#' @param fisher.v A numeric vector of Fisher information of \emph{v}.
#' @param Nsample The sample size in this group.
#' @param sample.type The sample type of this group of samples.
#' @param group1,group2 A vector indicating the first or second group of
#'   proteins according to the classification of \emph{v}. If
#'   \code{scaleV.raw.FC} is not \code{NULL}, the classification is conducted by
#'   the baseline fluctuation \eqn{log(log_2^2(1+scaleV.raw.FC))}.
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. See \code{estBioVar} and \code{resBioCMB} for more
#'   details.
#' @param outdir If \code{PLOT} is TRUE, the output directory.
#' @param outadd If \code{PLOT} is TRUE, the additions to the output documents.
#'
#' @return This function returns \code{NULL} invisibly.
#' @export
FisherInfoPlot <- function(m, v, fisher.m, fisher.v,
                           Nsample, sample.type,
                           group1, group2, scaleV.raw.FC = NULL,
                           outdir = "", outadd) {
  mkdir(paste0(outdir, "/"))
  if (!is.null(scaleV.raw.FC)) {
    main1 <- substitute(atop("Fisher info of m, lower bound" ~ log ~ (log[2]^2 ~ baseline),
                             Nsample ~ sample.type ~ "samples, # v " >= " cut: " ~ Num),
                        list(baseline = 1+scaleV.raw.FC, Nsample = Nsample,
                             sample.type = sample.type, Num = length(group1)))
    main2 <- substitute(atop("Fisher info of m, lower bound" ~ log ~ (log[2]^2 ~ baseline),
                             Nsample ~ sample.type ~ "samples, # v " < " cut: " ~ Num),
                        list(baseline = 1+scaleV.raw.FC, Nsample = Nsample,
                             sample.type = sample.type, Num = length(group2)))
    main3 <- substitute(atop("Fisher info of v, lower bound" ~ log ~ (log[2]^2 ~ baseline),
                             Nsample ~ sample.type ~ "samples, # v " >= " cut: " ~ Num),
                        list(baseline = 1+scaleV.raw.FC, Nsample = Nsample,
                             sample.type = sample.type, Num = length(group1)))
    main4 <- substitute(atop("Fisher info of v, lower bound" ~ log ~ (log[2]^2 ~ baseline),
                             Nsample ~ sample.type ~ "samples, # v " < " cut: " ~ Num),
                        list(baseline = 1+scaleV.raw.FC, Nsample = Nsample,
                             sample.type = sample.type, Num = length(group2)))
  } else {
    main1 <- paste0("Fisher info of m", " \n",
                    Nsample, " ", sample.type, " samples, # class 1: ", length(group1))
    main2 <- paste0("Fisher info of m", " \n",
                    Nsample, " ", sample.type, " samples, # class 2: ", length(group2))
    main3 <- paste0("Fisher info of v", " \n",
                    Nsample, " ", sample.type, " samples, # class 1: ", length(group1))
    main4 <- paste0("Fisher info of v", " \n",
                    Nsample, " ", sample.type, " samples, # class 2: ", length(group2))
  }

  # 4 subplots of Fisher information against m / v
  png(paste0(outdir, "/Fisher_info_", outadd, ".png"),
      height = 2000, width = 2000, res = 150)
  par(mfrow = c(2, 2))
  par(mar = c(6, 6, 5, 3))
  # subplot 1
  plot(m[group1], fisher.m[group1],
       col = alpha("blue", alpha = 0.3),
       xlab = "m", ylab = "Fisher m",
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main1)

  # subplot 2
  plot(m[group2], fisher.m[group2],
       col = alpha("blue", alpha = 0.3),
       xlab = "m", ylab = "Fisher m",
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main2)

  # subplot 3
  plot(v[group1], fisher.v[group1],
       col = alpha("blue", alpha = 0.3),
       xlab = "v", ylab = "Fisher v",
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main3)

  # subplot 4
  plot(v[group2], fisher.v[group2],
       col = alpha("blue", alpha = 0.3),
       xlab = "v", ylab = "Fisher v",
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main4)
  dev.off()

  invisible()
}


#' Plot mean-\emph{v}
#'
#' \code{FisherInfoPlot} would draw 2 pictures. The first one plot mean-\emph{v}
#' of all proteins. If \code{scaleV.raw.FC} is not \code{NULL}, the line of
#' lower bound is also added in the plot. The second one plot mean-\emph{v} of
#' proteins only in group1.
#'
#' @inheritParams resBioCMB
#' @param group1 A vector indicating the proteins in group1 according to the
#'   classification of \emph{v}. If \code{scaleV.raw.FC} is not \code{NULL}, the
#'   classification is by the lower bound \eqn{log(log_2^2(1+scaleV.raw.FC))}.
#'
#' @return This function returns \code{NULL} invisibly.
#' @export
MeanVarPlot <- function(ref.mean, v,
                        Nsample, sample.type,
                        group1, scaleV.raw.FC = NULL,
                        outdir, outadd) {
  mkdir(paste0(outdir, "/"))

  if (!is.null(scaleV.raw.FC)) {
    main1 <- substitute(atop("Mean-LogBioVar, # v " > ~log ~ (log[2]^2 ~ baseline) ~ ":" ~ Num,
                             Npro ~ "proteins," ~ Nsample ~ sample.type ~ "samples"),
                        list(baseline = 1+scaleV.raw.FC, Num = length(group1),
                             Npro = length(ref.mean), Nsample = Nsample, sample.type = sample.type))
    main2 <- substitute(atop("Mean-LogBioVar, # v " > ~log ~ (log[2]^2 ~ baseline) ~ ":" ~ Num,
                             Npro ~ "proteins," ~ Nsample ~ sample.type ~ "samples"),
                        list(baseline = 1+scaleV.raw.FC, Num = length(group1),
                             Npro = length(group1), Nsample = Nsample, sample.type = sample.type))
  } else {
    main1 <- paste0("Mean-LogBioVar, # high v cluster1: ", length(group1), "\n ",
                    length(ref.mean), " proteins, ",
                    Nsample, " ", sample.type, " samples")
    main2 <- paste0("Mean-LogBioVar, # high v cluster1: ", length(group1), "\n ",
                    length(group1), " proteins, ",
                    Nsample, " ", sample.type, " samples")
  }

  png(paste0(outdir, "/meanLogBioRawVar_", outadd, ".png"),
    height = 700, width = 700, res = 100)
  par(mar = c(5, 5, 5, 3))
  plot(ref.mean, v,
       xlab = "ref mean", ylab = "v = log(var)",
       col = alpha("blue", alpha = 0.2),
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main1)
  if (!is.null(scaleV.raw.FC)) {
    abline(h = log(log2(1+scaleV.raw.FC)^2),
           col = "darkorange", lwd = 2.5, lty = 5)
  }
  dev.off()

  png(paste0(outdir, "/meanLogBioHighVar_", outadd, ".png"),
      height = 700, width = 700, res = 100)
  par(mar = c(5, 5, 5, 3))
  plot(ref.mean[group1], v[group1],
       xlab = "ref mean", ylab = "v = log(var)",
       col = alpha("blue", alpha = 0.2),
       cex.main = 1.4, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3,
       main = main2)
  if (!is.null(scaleV.raw.FC)) {
    abline(h = log(log2(1+scaleV.raw.FC)^2),
           col = "darkorange", lwd = 2.5, lty = 5)
  }
  dev.off()

  invisible()
}
