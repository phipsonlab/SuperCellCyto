#' Run SuperCellCyto on Cytometry Data
#'
#' @description
#' This function creates supercells for a cytometry data formatted
#' as a \pkg{data.table} object using the SuperCellCyto algorithm.
#'
#' Please make sure you read additional details below to better understand
#' what the function does and how it works.
#'
#' @param dt A \pkg{data.table} object containing cytometry data where rows
#' represent cells and columns represent markers.
#' If this is a `data.frame` object, the function will try to convert it
#' to a \pkg{data.table} object.
#' A warning message will be displayed when this happens.
#' Otherwise, it will terminate.
#' @param markers A character vector identifying the markers to create
#' supercells with.
#' @param sample_colname A character string identifying the column in
#' \code{dt} that denotes the sample of a cell.
#' @param cell_id_colname A character string identifying the column in
#' \code{dt} representing each cell's unique ID.
#' @param n_pc A numeric value specifying the number of principal components
#' (PCs) to compute from the marker expression matrix.
#' Defaults to 10.
#' If there are less than 10 markers in the \code{markers} parameter,
#' then the number of PCs is set to however many markers there are in the
#' \code{markers} parameter.
#' @param gam A numeric value specifying the gamma value which regulates the
#' number of supercells generated.
#' Defaults to 20.
#' @param k_knn A numeric value specifying the k value (number of neighbours)
#' used to build the kNN network.
#' Defaults to 5.
#' @param BPPARAM A \link[BiocParallel]{BiocParallelParam-class}
#' object specifying the parallel processing settings.
#' Defaults to \link[BiocParallel]{SerialParam-class}, meaning the samples
#' will be processed sequentially one after the other.
#' Refer to additional details section below on parallel processing for
#' more details.
#' @param load_balancing A logical value  indicating whether to use a custom
#' load balancing scheme when processing multiple samples in parallel.
#' Defaults to FALSE.
#' Refer to additional details section below on parallel processing
#' for more details.
#' @param aggregation_method A character string specifying the method to be
#' used for calculating the marker expression of the supercells.
#' Accepted values are "mean" and "median".
#' Based on the choice, the supercells' marker expression are computed by
#' computing either the mean or median of the marker expression of the
#' cells therein.
#' The default value is "mean".
#' If any other value is provided, the function will return an error.
#'
#' @section Parallel Processing:
#' SuperCellCyto can process multiple samples simultaneously in parallel.
#' This can drastically bring down processing time for dataset with a large
#' number of samples.
#' To enable this feature, set the `BPPARAM` parameter to either a
#' \link[BiocParallel]{MulticoreParam-class} object or a
#' \link[BiocParallel]{SnowParam-class} object.
#' Importantly, it is also recommended to set the number of tasks (i.e., the
#' `task` parameter in either \link[BiocParallel]{MulticoreParam-class} or
#' \link[BiocParallel]{SnowParam-class} object) to **the number of samples in
#' the dataset**.
#'
#' Furthermore, we also recommend setting `load_balancing` parameter to TRUE.
#' This ensures optimal distribution of samples across multiple cores, and is
#' particularly important if your samples are of varying sizes
#' (number of cells).
#'
#' @section Cell ID and Sample Definitions:
#' The `cell_id_colname` parameter specifies the column in \code{dt} that
#' denotes the unique identifier for each cell.
#' It is perfectly normal to not have this column in your dataset by default.
#' The good news is that **it is trivial to create one**.
#' You can create a new vector containing a sequence of numbers from 1 to
#' however many cells you have, and append this vector as a new column in
#' your dataset.
#' Refer to our vignette on how to do this.
#'
#' The `sample_colname` parameter specifies the column in \code{dt} that
#' denotes the sample a cell came from.
#' By default, SuperCellCyto creates supercells for each sample independent of
#' other samples.
#' This ensures each supercell to only contain cells from **exactly**
#' one sample.
#'
#' What constitute a sample?
#' For most purposes, a sample represents a biological sample in
#' your experiment.
#' You may be thinking, is it then possible to use this in a different context,
#' say creating supercells for each population or cluster rather than a
#' biological sample?
#' The short answer is yes, and we address this in our vignette.
#'
#' @section Computing PCA:
#' The function will start by computing PCA from all the markers
#' specified in `markers` parameter.
#' By default, the number of PCs calculated is set to 10.
#' If there is less than 10 markers in the \code{markers} parameter, then the
#' number of PCs is set to however many markers there are in the \code{markers}
#' parameter.
#'
#' Notably, \emph{no} scaling or transformation were done on the markers'
#' expressions prior to computing the PCs.
#'
#' The function does not use `irlba` to calculate PCA.
#' There is very little gain to use it for cytometry data because of the
#' relatively tiny number of features (markers) in the data.
#'
#' @section Setting Supercell's Resolution:
#' The `gam` parameter influences the number of supercells created per sample.
#' A lower `gam` value results in more, and thus generally smaller supercells,
#' and vice versa.
#'
#' To estimate how many supercells we will get for our dataset, it is important
#' to understand how the `gam` value is interpreted in the context of number of
#' cells in a sample.
#'
#' `gam=n_cells/n_supercells`
#' where `n_cells` denotes the number of cells and `n_supercells` denotes the
#' number of supercells to be created.
#'
#' By resolving the formula above, we can roughly estimate how many supercells
#' we will get per sample.
#' For example, say we have 2 samples, sample A and B.
#' Sample A has 10,000 cells, while sample B has 5,000 cells:
#'
#' * If `gam` is set to 10, we will end up with 1,000 supercells for sample A
#' and 500 supercells for sample B, a total of 1,500 supercells.
#' * If `gam` is set to 50, we will end up with 200 supercells for sample A
#' and 100 supercells for sample B, a total of 300 supercells.
#'
#' *Importantly*, one cannot expect all the supercells to be of the same size.
#' Some will capture more/less cells than others.
#' It is not trivial to estimate how many cells will be captured in each
#' supercell beforehand.
#'
#' @section Computing kNN network:
#' To create supercells, a kNN (k-Nearest Neighbors) network is constructed
#' based on the `k_knn` parameter which dictates the number of neighbours
#' (for each cell) used to create the network.
#' An actual (not approximate) kNN network is created.
#'
#' A walktrap algorithm then uses this network to group cells into supercells.
#'
#'
#' @return
#' A list with the following components:
#'
#' * `supercell_object`: A list containing the object returned by
#' SCimplify function.
#' One object per sample.
#' This object is critical for recomputing supercells in the future.
#' Hence do not discard it.
#' * `supercell_expression_matrix`:  A \pkg{data.table} object that contains
#' the marker expression for each supercell.
#' These marker expressions are computed by calculating the mean of the marker
#' expressions across all cells
#' within each individual supercell.
#' * `supercell_cell_map`: A \pkg{data.table} that maps each cell to its
#' corresponding supercell.
#' This table is essential for identifying the specific supercell each cell
#' has been allocated to.
#' It proves particularly useful for analyses that require one to expand
#' the supercells to the individual cell level.
#'
#' @author
#' Givanna Putri
#'
#' @examples
#' # Simulate some data
#' set.seed(42)
#' cyto_dat <- simCytoData(nmarkers = 10, ncells = rep(2000,2))
#'
#' # Setup the columns designating the markers, samples, and cell IDs
#' marker_col <- paste0("Marker_", seq_len(10))
#' sample_col <- "Sample"
#' cell_id_col <- "Cell_Id"
#'
#' supercell_dat <- runSuperCellCyto(
#'     cyto_dat, marker_col,
#'     sample_col, cell_id_col
#' )
#'
#' @export
#'
#' @import data.table
#' @import BiocParallel
#' @importFrom Matrix Matrix
#' @importFrom SuperCell SCimplify supercell_GE
#' @importFrom stats median
#'
runSuperCellCyto <- function(
    dt, markers, sample_colname, cell_id_colname,
    n_pc = 10, aggregation_method = c("mean", "median"),
    gam = 20, k_knn = 5, BPPARAM = SerialParam(), load_balancing = FALSE
) {
    stopifnot(is.data.frame(dt) == TRUE)
    if (!is.data.table(dt)) {
        warning(
            "dt is not a data.table object. Converting to a data.table object."
        )
        dt <- as.data.table(dt)
    }

    # Convert data.table to a list of matrices, one per sample
    # If load_balancing is TRUE, sort the samples by number of cells
    matrix_per_samp <- .get_all_sample_marker_matrices(
        dt = dt,
        samples = unique(dt[[sample_colname]]),
        markers = markers,
        sample_colname = sample_colname,
        cell_id_colname = cell_id_colname,
        load_balancing = load_balancing
    )
    # If there are more n_pc than markers, set n_pc to the number of markers
    n_pc <- .adjust_n_pc(n_pc, markers)

    # How to aggregate the cells in supercells to get expression matrix?
    aggregation_method <- match.arg(aggregation_method)

    supercell_res <- .run_supercells_bplapply(
        matrix_per_samp = matrix_per_samp,
        n_pc = n_pc,
        gam = gam,
        k_knn = k_knn,
        aggregation_method = aggregation_method,
        sample_colname = sample_colname,
        cell_id_colname = cell_id_colname,
        markers = markers,
        BPPARAM = BPPARAM
    )

    # Reshape the output
    reshaped_res <- .create_supercell_output(
        supercell_res = supercell_res,
        samples = names(matrix_per_samp)
    )

    return(reshaped_res)
}

#' Adjust n_pc if it exceeds the number of markers
#'
#' Internal helper to ensure n_pc does not exceed the number of markers.
#' 
#' @param n_pc Numeric value specifying the number of principal components.
#' @param markers Character vector of markers used to compute supercells.
#' 
#' @return Adjusted n_pc value.
#' 
#' @keywords internal
#' @noRd
#' 
.adjust_n_pc <- function(n_pc, markers) {
    if (n_pc > length(markers)) {
        warning(
            sprintf(
                "n_pc (%d) > number of markers (%d). Setting n_pc to %d.",
                n_pc, length(markers), length(markers)
            )
        )
        n_pc <- length(markers)
    }
    return(n_pc)
}

#' Convert data.table into list of matrices
#'
#' Internal helper to convert a data.table containing marker expression into
#' a list of matrices, one matrix per sample, and 
#' optionally reorder for load balancing.
#' 
#' @keywords internal
#' @noRd
.get_all_sample_marker_matrices <- function(
    dt, samples, markers, sample_colname, 
    cell_id_colname, load_balancing = FALSE
) {
    matrix_per_samp <- lapply(samples, function(samp) {
        .get_sample_marker_matrix(
            dt = dt,
            sample_name = samp,
            markers = markers,
            sample_colname = sample_colname,
            cell_id_colname = cell_id_colname
        )
    })
    names(matrix_per_samp) <- samples

    if (load_balancing) {
        matrix_per_samp <- .sort_samples_by_cell_count(
            matrix_per_samp, samples
        )
    }
    return(matrix_per_samp)
}


#' Sort samples and their matrices by cell count (descending)
#'
#' Internal helper to sort sample matrices by number of cells
#' for load balancing.
#' Sort the samples based on how many cells they have, starting from
#' the one with the most cells, all the way to the one with the
#' least number of cells.
#' The idea is that after we do this, we can create MulticoreParam
#' where we set the number of tasks to as many sample as we have, such that
#' each sample is sent to each worker. then when the worker returns
#' it will then send the next sample to process to the worker.
#'
#' @param matrix_per_samp List of matrices, one per sample.
#' @param samples Character vector of sample names.
#'
#' @return List with reordered matrix_per_samp and ncells_per_sample.
#' @keywords internal
#' @noRd
#'
.sort_samples_by_cell_count <- function(matrix_per_samp, samples) {

    ncells_per_sample <- vapply(
        matrix_per_samp, ncol,
        FUN.VALUE = integer(1)
    )
    names(ncells_per_sample) <- samples

    ncells_per_sample <- ncells_per_sample[order(
        ncells_per_sample, decreasing = TRUE
    )]
    matrix_per_samp_ordered <- matrix_per_samp[names(ncells_per_sample)]

    return(matrix_per_samp_ordered)
}

#' Run supercells computation for all samples in parallel.
#' @param matrix_per_samp List of marker expression matrices, one per sample.
#' @param n_pc Numeric value specifying the number of principal components
#' to compute.
#' @param gam Numeric value specifying the gamma value for supercell resolution.
#' @param k_knn Numeric value specifying the number of neighbours for kNN.
#' @param aggregation_method Character string specifying the method to compute
#' supercell's centroid.
#' @param sample_colname Name of the sample column.
#' @param cell_id_colname Name of the cell ID column.
#' @param markers Character vector of marker names.
#' @param BPPARAM A \link[BiocParallel]{BiocParallelParam-class}
#' object specifying the parallel processing settings.  
#' 
#' @return A list of supercell related results, one per sample.
#' Each element in the list contains:
#' * `supercell_object`: The supercell object returned by SCimplify.
#' * `supercell_expression_matrix`: A matrix of supercells' centroids.
#' * `supercell_cell_map`: A data.table mapping each cell to its supercell.
#' 
#' @keywords internal
#' @noRd
#' 
.run_supercells_bplapply <- function(
    matrix_per_samp, n_pc, gam, k_knn, aggregation_method,
    sample_colname, cell_id_colname, markers, BPPARAM
) {
    res <- bplapply(names(matrix_per_samp),
        function(sample_name, gam, k_knn, aggregation_method) {
            mt <- matrix_per_samp[[sample_name]]
            res <- .compute_supercells_one_sample(
                sample_name = sample_name,
                mt = mt,
                n_pc = n_pc,
                gam = gam,
                k_knn = k_knn,
                aggregation_method = aggregation_method,
                sample_colname = sample_colname,
                cell_id_colname = cell_id_colname,
                markers = markers
            )
            return(res)
        }, gam = gam, k_knn = k_knn, 
        aggregation_method = aggregation_method, BPPARAM = BPPARAM
    )
    return(res)
}

#' Create supercells for one sample
#'
#' Internal helper to run SCimplify to compute supercell outputs 
#' for one sample.
#' 
#' @param sample_name Character string representing the sample name.
#' @param mt A matrix of marker expression for the sample.
#' @param n_pc Numeric value specifying the number of principal components
#' to compute.
#' @param gam Numeric value specifying the gamma value for supercell resolution.
#' @param k_knn Numeric value specifying the number of neighbours for kNN.
#' @param aggregation_method Character string specifying the method to compute
#' supercell's centroid.
#' @param sample_colname Character string identifying the column in \code{mt}
#' that denotes the sample of a cell.
#' @param cell_id_colname Character string identifying the column in \code{mt}
#' representing each cell's unique ID.
#' @param markers Character vector of markers used to compute supercells.
#' 
#' @return A list containing:
#' * `supercell_object`: The object returned by SCimplify.
#' * `supercell_expression_matrix`: A matrix of supercells' centroids.
#' * `supercell_cell_map`: A data.table mapping each cell to its supercell.
#' 
#' @keywords internal
#' @noRd
#' 
.compute_supercells_one_sample <- function(
    sample_name, mt, n_pc, gam, k_knn, aggregation_method,
    sample_colname, cell_id_colname, markers
) {
    
    # Calculate supercells using SCimplify
    res <- SCimplify(
        X = mt,
        genes.use = rownames(mt),
        do.scale = FALSE,
        do.approx = FALSE,
        gamma = gam,
        k.knn = k_knn,
        fast.pca = FALSE,
        n.pc = n_pc
    )

    supercell_exp_mat <- .compute_supercell_centroid(
        mt = mt,
        sc_membership = res$membership,
        aggregation_method = aggregation_method,
        sample_colname = sample_colname,
        sample_name = sample_name,
        markers = markers
    )

    supercell_cell_map <- .create_supercell_cell_map(
        sc_membership = res$membership,
        cell_ids = colnames(mt),
        sample_name = sample_name
    )

    return(list(
        supercell_object = res,
        supercell_expression_matrix = supercell_exp_mat,
        supercell_cell_map = supercell_cell_map
    ))
}

#' Create supercell output
#'
#' Internal helper to create supercell output.
#' 
#' @param supercell_res List of supercell results for each sample.
#' @param samples Character vector of sample names.
#' 
#' @return A list containing reshaped supercell results:
#' * `supercell_expression_matrix`: A data.table of supercell marker 
#' expressions.
#' * `supercell_cell_map`: A data.table mapping cells to supercells.
#' * `supercell_object`: A list of supercell objects, one per sample.
#' 
#' @keywords internal
#' @noRd
#' 
.create_supercell_output <- function(supercell_res, samples) {
    reshaped_res <- .reshape_supercell_output(
        supercells = supercell_res
    )
    reshaped_res$supercell_object <- lapply(
        supercell_res, function(res_i) res_i$supercell_object
    )
    names(reshaped_res$supercell_object) <- samples
    return(reshaped_res)
}