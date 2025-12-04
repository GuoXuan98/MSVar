# Quality control

#' Quality control of raw intensity
#'
#' Given the raw intensity of all tissue samples, plot the number of detected
#' samples for each protein and the number of detected proteins in each sample.
#'
#' @param x A \code{\link{proObj}}.
#' @param outdir The output directory.
#' @param outadd A character string to describe the additional information of
#'   plots.
#' @param cpus The number of CPUs requested for a parallel computation.
#'
#' @return This function returns \code{NULL} invisibly.
#' @export
dataQC <- function(x,
                   outdir = "./", outadd = "", cpus = 5){
  mkdir(outdir)
  sample.intensity <- x$sample.intensity
  Nsample <- ncol(sample.intensity)

  SampleNumber <- function(index, sample.intensity) {
    length(which(!is.na(sample.intensity[index, ]) &
                   sample.intensity[index, ] != 0))
  }

  GeneNumber <- function(index, sample.intensity) {
    length(which(!is.na(sample.intensity[, index]) &
                   sample.intensity[, index] != 0))
  }

  sfInit(parallel = TRUE, cpus = cpus)
  detected.sample.number.each.gene <- sfSapply(1:nrow(sample.intensity),
                                               SampleNumber, sample.intensity)
  sfStop()

  sfInit(parallel = TRUE, cpus = cpus)
  detected.gene.number.each.sample <- sfSapply(1:ncol(sample.intensity),
                                               GeneNumber, sample.intensity)
  sfStop()

  png(paste0(outdir, "/qc_boxplot", outadd, ".png"),
      width = 600, height = 1000, res = 72*2)
  dev.control('enable')

  par(family = "sans", mar = c(6, 6, 5, 3))

  boxplot(detected.gene.number.each.sample,
          xlab = "Sample", ylab = "Detected gene numbers",
          main = paste0("Number of gene detected \n in at least one sample: ",
                        length(which(detected.sample.number.each.gene != 0))),
          cex.main = 1.6, cex = 1.2, cex.lab = 1.4, cex.axis = 1.3)
  dev.off()


  # Plot the table.
  half.detected.num <- sum(detected.sample.number.each.gene > 0.5 * Nsample)
  three.quarters.detected.num <- sum(detected.sample.number.each.gene > 0.75 * Nsample)
  all.detected.num <- sum(detected.sample.number.each.gene > (Nsample - 1))

  describ.table <- data.frame(Description = c("Proteins quantified in at least\nhalf samples at gene level",
                                              "Proteins quantified in at least three\nquarters samples at gene level",
                                              "Proteins quantified in all \nsamples at gene level"),
                              Counts = c(half.detected.num, three.quarters.detected.num, all.detected.num))

  png(paste0(outdir, "/qc_table", outadd, ".png"),
      height = 400, width = 800, res = 72*3)
  dev.control('enable')

  grid.table(describ.table, rows = NULL,
             theme = ttheme_minimal(core = list(bg_params = list(fill = c("#e8f3de", "#d3e8bb"),
                                                                 col = "black"),
                                                fg_params = list(col = "black", fontface = 4L)),
                                    colhead = list(bg_params = list(fill = "#8cc257", col = "black"),
                                                   fg_params = list(col = "white", fontface = 4L))))
  dev.off()

  invisible()
}
