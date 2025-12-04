# For constructing and demonstrating objects of the class "proObj".

#' Create a \code{proObj}
#'
#' This function constructs a new \code{proObj} based on the given raw proteins
#' intensity of a group of samples and reference samples, and the corresponding
#' relationship is indicated by \code{sample.ref.cor}.
#'
#' @param sample.intensity A numeric matrix of a group of samples' raw
#'   intensity.
#' @param ref.intensity A numeric matrix of all reference samples' raw
#'   intensity.
#' @param sample.ref.cor A named character vector which records the
#'   corresponding reference sample name of each tissue sample in
#'   \code{sample.intensity}.
#' @param IDs A character vector of proteins' names.
#' @param name A character scalar specifying the name of this \code{proObj}.
#'   Used only for demonstration.
#'
#' @return A newly constructed object of the \code{\link[base]{class}}
#'   \code{"proObj"}.
#'
#' @seealso \code{\link{proObjFromTMT}} for creating a \code{proObj} from
#'   general result from TMT experiment.
#'
#' @export
#'
proObj <- function(sample.intensity, ref.intensity, sample.ref.cor,
                   IDs = NULL, name = "NA") {
  Npro <- nrow(sample.intensity)
  Nsample <- ncol(sample.intensity)

  if (is.null(names(sample.ref.cor))) {
    stop("sample.ref.cor can not indicate the relationship \nbetween tissue samples and reference samples.")
  }
  if (sum(names(sample.ref.cor) == colnames(sample.intensity)) != Nsample) {
    stop("The length of samples-reference correspondence \nmust match that of sample.intensity.")
  }
  if (any(is.na(sample.ref.cor))) {
    stop("NA values are not allowed in sample.ref.cor.")
  }

  if (nrow(ref.intensity) != Npro) {
    stop("The protein number in the reference samples \nare not the same as that in the group of samples.")
  }

  if (is.null(IDs)) {
    if (!is.null(rownames(sample.intensity))) {
      IDs <- rownames(sample.intensity)
    } else {
      IDs <- paste0("protein", 1:Npro)
      rownames(sample.intensity) <- rownames(ref.intensity) <- IDs
    }
  } else {
    IDs <- as.character(IDs)
    if (length(IDs) != Npro) {
      stop("The length of IDs must match that of sample.intensity.")
    }
    if (any(is.na(IDs))) {
      stop("NA values are not allowed in IDs.")
    }
    if (length(unique(IDs)) != Npro) {
      stop("Duplicates are not allowed in IDs.")
    }
  }

  sample.intensity <- sample.intensity[intersect(IDs, rownames(sample.intensity)), , drop = F]
  ref.intensity <- ref.intensity[intersect(IDs, rownames(ref.intensity)), , drop = F]

  all_detected_pair <- sapply(1:nrow(sample.intensity), function(i) {
    pair <- sapply(1:length(sample.ref.cor), function(k)
      sample.intensity[i, names(sample.ref.cor)[k]] - ref.intensity[i, sample.ref.cor[k]])
    !all(is.na(pair))
  })

  all_detect <- which(all_detected_pair == 1)

  x <- list(sample.intensity = sample.intensity[all_detect, , drop = F],
            ref.intensity = ref.intensity[all_detect, , drop = F],
            sample.ref.cor = sample.ref.cor,
            name = as.character(name)[1],
            pro.name = IDs[all_detect])
  class(x) <- "proObj"
  x
}

#' Create a \code{proObj}
#'
#' This function constructs a new \code{proObj} based on the general TMT
#' experiment result, which include raw intensity of all samples and batch
#' information of the MS experiment.
#'
#' @param raw.intensity A numeric matrix of raw intensity, which includes a
#'   group of samples and reference samples.
#' @param batch.info A character matrix of batch information, which indicates
#'   the corresponding relationship of a group of samples among batches. Each
#'   row includes the reference sample name and sample names in that batch.
#' @param IDs A character vector of proteins' names.
#' @param name A character scalar specifying the name of this \code{proObj}.
#'   Used only for demonstration.
#'
#' @return A newly constructed object of the \code{\link[base]{class}}
#'   \code{"proObj"}.
#'
#' @seealso \code{\link{proObj}} for another constructor of \code{proObj}s.
#'
#' @export
#'
proObjFromTMT <- function(raw.intensity, batch.info,
                          IDs = NULL, name = "NA") {
  batch.col <- grepl("ref", colnames(batch.info))
  sample.ref.cor <- rep(batch.info[, batch.col],
                        each = length(which(batch.col == F)))
  sample.name <- as.vector(t(batch.info[, !batch.col]))

  sample.ref.cor <- sample.ref.cor[which(nchar(sample.name) > 0)]
  names(sample.ref.cor) <- sample.name[which(nchar(sample.name) > 0)]
  sample.ref.cor <- sample.ref.cor[intersect(names(sample.ref.cor), colnames(raw.intensity))]

  sample.intensity <- raw.intensity[, names(sample.ref.cor), drop = F]
  ref.intensity <- raw.intensity[, unique(sample.ref.cor), drop = F]

  proObj(sample.intensity = sample.intensity,
         ref.intensity = ref.intensity,
         sample.ref.cor = sample.ref.cor,
         IDs = IDs, name = name)
}


#' Print a \code{proObj}
#'
#' This function prints its argument, which is a \code{\link{proObj}}, and
#' returns it invisibly (via \code{\link[base]{invisible}(x)}).
#'
#' This function implements the \code{\link[base]{print}} method for the
#' \code{"\link{proObj}"} class.
#'
#' @param x A \code{\link{proObj}}.
#' @param ... Arguments passed from other methods.
#'
#' @return The function returns \code{x} invisibly.
#'
#' @seealso \code{\link{proObj}} for creating a \code{proObj}.
#'
#' @export
#' @export print.proObj
#'
print.proObj <- function(x, ...) {
  temp <- paste0(x$name, ":")
  cat(paste("proObj", temp,
            nrow(x$sample.intensity), "proteins"),
      sep = "\n")
  invisible(x)
}
