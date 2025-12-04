# For the documentation of MSVar as a whole package, as well as the data sets
# exported by MSVar for demonstrating its usage.



# Display startup message.
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("MSVar 1.0.0 2022-07-14")
# }


#' MSVar for Identifying Hypervariable proteins and Differentially Variable
#' proteins
#'
#' \code{MSVar} is designed for isobaric labeling-based mass spectrometry (ILMS)
#' quantitative proteomic profiling to identify hypervariable proteins for a
#' single group of samples and differentially variable proteins between
#' different sample groups.
#'
#'
#' @section Maintainer: Xuan Guo <\email{guoxuan2020@@sinh.ac.cn}>
#'
#' @docType package
#' @name MSVar
#'
#' @importFrom stats var
#' @importFrom stats predict
#' @importFrom stats lm
#' @importFrom stats nlm
#' @importFrom stats coefficients
#' @importFrom stats quantile
#' @importFrom stats p.adjust
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats kmeans
#'
#' @importFrom graphics plot
#' @importFrom graphics boxplot
#' @importFrom graphics hist
#' @importFrom graphics par
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics text
#'
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.control
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.copy
#' @importFrom grDevices dev.off
#'
#' @importFrom gridExtra grid.table
#' @importFrom gridExtra ttheme_minimal
#'
#' @importFrom methods is
#'
#' @importFrom utils write.table
#'
#' @importFrom locfit lp
#' @importFrom locfit locfit
#'
#' @importFrom scales alpha
#'
#' @importFrom splines ns
#'
#' @importFrom snowfall sfInit
#' @importFrom snowfall sfLibrary
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfSapply
#' @importFrom snowfall sfStop
#'
#'
NULL
