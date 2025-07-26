#' Simulate cytometry data
#'
#' Simulate some cytometry data for use in testing or documenting functions
#' which require some cytometry data.
#' Please run \code{set.seed} before running the function if you want to ensure
#' reproducibility.
#'
#' @param nmarkers A numeric value specifying number of markers to simulate.
#' @param ncells A numeric vector specifying the number of cells to simulate 
#' per sample.
#' 1 vector element per sample.
#'
#' @return
#' A \link{data.table} object containing the simulated cytometry data where rows 
#' represent cells and columns represent markers.
#'
#' @examples
#' set.seed(42)
#' cyto_dat <- simCytoData()
#' head(cyto_dat)
#' dim(cyto_dat)
#'
#' @author
#' Givanna Putri
#'
#' @import data.table
#' @importFrom stats rnorm runif
#'
#' @export
#'
simCytoData <- function(nmarkers = 10, ncells = rep(10000, 2)) {

    cyto_data <- lapply(seq(length(ncells)), function(samp_idx) {
        rnorm_mean <- runif(nmarkers, min = 5, max = 20)
        markers <- lapply(rnorm_mean, function(m) {
            return(rnorm(ncells[samp_idx], mean = m))
        })

        out <- data.table(do.call(cbind, markers))
        names(out) <- paste0("Marker_", seq(nmarkers))
        out$Sample <- paste0("Sample_", samp_idx)
        return(out)
    })


    cyto_data <- rbindlist(cyto_data)
    cyto_data$Cell_Id <- paste0("Cell_", seq(nrow(cyto_data)))

    return(cyto_data)
}
