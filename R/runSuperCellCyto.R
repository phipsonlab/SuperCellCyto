#' Run SuperCellCyto on Cytometry Data
#'
#' @description
#' This function creates supercells for a cytometry data formatted 
#' as a \link{data.table} object using the SuperCellCyto algorithm. 
#' 
#' Please make sure you read additional details below to better understand
#' what the function does and how it works.
#' 
#' @param dt A \link{data.table} object containing cytometry data where rows 
#' represent cells and columns represent markers.
#' If this is a \link{data.frame} object, the function will try to convert it 
#' to a \link{data.table} object.
#' A warning message will be displayed when this happens.
#' Otherwise, it will terminate.
#' @param markers A character vector identifying the markers to create 
#' supercells with.
#' @param sample_colname A character string identifying the column in 
#' \code{dt} that denotes the sample of a cell.
#' @param cell_id_colname A character string identifying the column in 
#' \code{dt} representing each cell's unique ID.
#' @param gam A numeric value specifying the gamma value which regulates the 
#' number of supercells generated.
#' Defaults to 20.
#' @param k_knn A numeric value specifying the k value (number of neighbours)
#' used to build the kNN network.
#' Defaults to 5.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying parallel 
#' processing settings.
#' Defaults to `SerialParam`, meaning the samples will be processed 
#' sequentially one after the other.
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
#' \linkS4class{MulticoreParam} object or a \linkS4class{SnowParam} object.
#' Importantly, it is also recommended to set the number of tasks (i.e., the 
#' `task` parameter in either \linkS4class{MulticoreParam} or 
#' \linkS4class{SnowParam} object) to **the number of samples in the dataset**.
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
#' The `sample_colname` parameter specifies the column in \code{dt} that denotes
#' the sample a cell came from.
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
#' * `supercell_expression_matrix`:  A \link{data.table} object that contains 
#' the marker expression for each supercell.
#' These marker expressions are computed by calculating the mean of the marker 
#' expressions across all cells
#' within each individual supercell.
#' * `supercell_cell_map`: A \link{data.table} that maps each cell to its 
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
runSuperCellCyto <- function(dt,
                             markers,
                             sample_colname,
                             cell_id_colname,
                             aggregation_method = c("mean", "median"),
                             gam = 20,
                             k_knn = 5,
                             BPPARAM = SerialParam(),
                             load_balancing = FALSE) {
    # Check data type first, and error out if dt is not a data.frame
    stopifnot(is.data.frame(dt) == TRUE)

    # Convert dt to data.table if it is not
    if (!is.data.table(dt)) {
        message("dt is not a data.table object. Converting it to a data.table object")
        dt <- as.data.table(dt)
    }

    samples <- unique(dt[[sample_colname]])

    matrix_per_samp <- lapply(samples, function(samp) {
        dt_sub <- dt[dt[[sample_colname]] == samp, ]
        trans_dt_sub <- Matrix(t(dt_sub[, markers, with = FALSE]))
        colnames(trans_dt_sub) <- dt_sub[[cell_id_colname]]
        return(trans_dt_sub)
    })
    names(matrix_per_samp) <- samples
    
    # ---- Load balancing ----
    # Sort the samples based on how many cells they have, starting from
    # the one with the most cells, all the way to the one with the 
    # least number of cells.
    # The idea is that after we do this, we can create MulticoreParam
    # where we set the number of tasks to as many sample as we have, such that
    # each sample is sent to each worker. then when the worker returns 
    # it will then send the next sample to process to the worker. 
    
    if (load_balancing) {
        ncells_per_sample <- vapply(matrix_per_samp, ncol)
        names(ncells_per_sample) <- samples
        
        ncells_per_sample <- ncells_per_sample[order(ncells_per_sample, 
                                                     decreasing = TRUE)]
        matrix_per_samp <- matrix_per_samp[names(ncells_per_sample)]
    }
    
    # Number of PCs are set to 10 by default. We can have panel size 
    # less than 10.
    # If this is the case, we just set PCA to be the number of markers
    if (length(markers) < 10) {
        n_pc <- length(markers)
    } else {
        n_pc <- 10
    }
    
    # How to aggregate the cells in supercells to get expression matrix?
    aggregation_method <- match.arg(aggregation_method)
    
    supercell_res <- bplapply(names(matrix_per_samp), function(
        sample_name, gam, k_knn, aggregation_method) {
        mt <- matrix_per_samp[[sample_name]]
        
        # ---- Run supercell ----
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
        
        # ---- Calculate supercell expression matrix ----
        if (aggregation_method == "mean") {
            supercell_exp_mat <- data.table(
                t(
                    as.matrix(
                        supercell_GE(
                            ge = mt, 
                            groups = res$membership
                        )
                    )
                )
            )
            supercell_exp_mat[[sample_colname]] <- sample_name
            
            # Create a unique supercell id concatenating the sample name
            supercell_exp_mat[["SuperCellId"]] <- paste0(
                "SuperCell_",
                seq(1, nrow(supercell_exp_mat)), "_Sample_", sample_name
            )
        } else if (aggregation_method == "median") {
            supercell_exp_mat <- data.table(t(as.matrix(mt)))
            
            # grab it here now!
            markers_name <- colnames(supercell_exp_mat)
            
            supercell_exp_mat$cell_id <- colnames(mt)
            supercell_membership <- data.table(
                cell_id = names(res$membership),
                SuperCellId = paste0("SuperCell_", res$membership, 
                                     "_Sample_", sample_name)
            )
            supercell_exp_mat <- merge.data.table(
                supercell_exp_mat,
                supercell_membership,
                by = "cell_id"
            )
            
            # this is where you calculate the expression
            supercell_exp_mat <- supercell_exp_mat[, lapply(.SD, median), 
                                                   .SDcols = markers, 
                                                   by='SuperCellId'] 
        }
        
        
        
        # ---- Create supercell and cell mapping ----
        supercell_cell_map <- data.table(
            SuperCellID = paste0(
                "SuperCell_", res$membership,
                "_Sample_", sample_name
            ),
            CellId = colnames(mt),
            Sample = sample_name
        )
        
        # Return a list containing all the objects
        return(list(
            supercell_object = res,
            supercell_expression_matrix = supercell_exp_mat,
            supercell_cell_map = supercell_cell_map
        ))
        
    }, gam = gam, k_knn = k_knn, aggregation_method = aggregation_method, BPPARAM = BPPARAM)

    # Now the messy reshaping so each element is not the output for a sample
    # but either a supercell object, expression matrix or supercell cell map
    reshaped_res <- list(
        supercell_expression_matrix = do.call(rbind, lapply(supercell_res, function(res_i) res_i$supercell_expression_matrix)),
        supercell_cell_map = do.call(rbind, lapply(supercell_res, function(res_i) res_i$supercell_cell_map)),
        supercell_object = lapply(supercell_res, function(res_i) res_i$supercell_object)
    )
    names(reshaped_res$supercell_object) <- names(matrix_per_samp)
    
    return(reshaped_res)
}

