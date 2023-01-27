#' Simulates cytometry data
#' 
#' Simulates some cytometry data for use in testing or documenting functions
#' which require some cytometry data.
#'
#' @param nmarkers Numeric scalar specifying number of markers to simulate.
#' @param ncells Numeric scalar specifying number of cells to simulate per sample.
#' @param nsample Numeric scalar specifying number of samples to simulate.
#'
#' @return
#' A \link{data.table} containing the simulated data. 
#' Rows are cells, columns are markers. 
#'
#' @author 
#' Givanna Putri
#' 
#' @import data.table
#' 
#' @export
simCytoData <- function(nmarkers=10, ncells=10000, nsample=2) {
    
    cyto_data <- lapply(seq_len(nsample), function(samp) {
        rnorm_mean <- runif(nmarkers, min=5, max=20)
        markers <- lapply(rnorm_mean, function(m) {
            return(rnorm(ncells, mean=m))
        })
        
        out <- data.table(do.call(cbind, markers))
        names(out) <- paste0("Marker_", seq_len(nmarkers))
        out$Sample <- paste0("Sample_", samp)
        return(out)
    })
    cyto_data <- rbindlist(cyto_data)
    cyto_data$Cell_Id <- paste0("Cell_", seq_len(nrow(cyto_data)))
    
    return(cyto_data)
    
}