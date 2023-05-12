#' Run SuperCell for cytometry data
#'
#' @description
#' Run SuperCell on cytometry data stored as a \link{data.table} object.
#' This is a wrapper function around the `SCImplify` function in the
#' SuperCell R package by Bilous et.al, 2022.
#' We have enhanced it by adding the capacity to "supercell" multiple samples
#' in parallel through the use of \link{BiocParallel}, and by adding support to
#' "supercell" cytometry data.
#' More explanations are given in various sections below, and are further
#' expanded in our vignette.
#'
#' @param dt \link{data.table} object containing the cytometry data.
#' Rows represent cells, columns represent markers.
#' If this is not a \link{data.table} object, the function will warn you about it,
#' and then try to convert it to a \link{data.table} object.
#' @param markers character vector specifying the markers to run SuperCell on.
#' @param sample_colname character specifying the column in \code{dt} that denotes
#' the sample of a cell.
#' @param cell_id_colname character specifying the column in \code{dt} that denotes
#' the unique ID of a cell.
#' @param gam numeric specifying the gamma value to be used by SuperCell.
#' Default to 20.
#' @param k_knn numeric specifying the k value to be used by SuperCell's knn.
#' Default to 5.
#' @param n_parallel_worker numeric specifying the number of parallel jobs/workers
#' to use to supercell the samples.
#' Default to 1, meaning the samples are processed sequentially. 
#' If you want to process the samples in parallel, set this number to 1, but
#' no more that the maximum number of cores you have in the computer.
#' You can use parallel::detectCores to find out how many cores you have in the computer.
#' It is not recommended to set this to the number returned by parallel::detectCores
#' as it will render your computer unusable for anything else. 
#' It is best to set this to either 1 or 2 cores less than the total you have in the computer.
#'
#' @section What is \code{cell_id_colname}:
#' This is a column in \code{dt} containing a unique identifier for each cell.
#' Commonly, you will have to manually create this column as a FCS file does not
#' typically has a field which uniquely identify each cell.
#' You can create this ID by giving the cells a numeric value of 1 to however many 
#' you have, and store this as a column in \code{dt}.
#' If you don't know how to do this, refer to our vignette.
#'
#' @section Processing one sample independent of the others:
#' This function is designed such that all the samples are processed,
#' independent of each other.
#' Because of this, you can safely assume that each supercell will contain only
#' cells from exactly one sample.
#'
#' The function will work out which cells come from which sample based on what is
#' specified in the \code{sample_colname} column.
#' This is why it is critical to specify what this column is.
#'
#' For most purposes, a sample represents a biological sample in your experiment.
#' You may be thinking, is it then possible to use this in a different context,
#' say creating supercells for each population or cluster rather than a biological
#' sample?
#' The short answer is yes, and we address this in our vignette.
#'
#' @section Computing PCA:
#' By default, the function will start by computing PCA from all the markers
#' specified in \code{markers} parameter, and that 10 PCs are computed.
#' If there are less than 10 markers in the \code{markers} parameter, then the
#' number of PCs are set to however many markers there are in the \code{markers} parameter.
#'
#' Notably, \emph{no} scaling or transformation were done on the markers' expressions
#' prior to computing the PCs.
#'
#' \code{irlba} is not used to calculate PCA as cytometry data tend to only have a
#' handful of features (markers) compared to scRNAseq data.
#' Hence there is very little gain.
#'
#' @section Setting the supercell graining level:
#' How many supercells will I get for my dataset? or to phrase it in another way,
#' can I estimate how many cells will be captured within each supercell?
#' That depends on what you set the \code{gam} parameter to.
#'
#' The \code{gam} parameter is represented by the formula `gamma=n_cells/n_supercells`
#' where `n_cells` denotes the number of cells and `n_supercells` denotes the
#' number of supercells to be created.
#' By resolving this formula, we can roughly estimate how many supercells you will
#' get at the end, and thus, \emph{approximately} how many cells will be captured
#' within each supercell.
#'
#' Generally, the smaller the \code{gam} parameter is, the more supercells
#' you will get.
#' Say for instance you have 10,000 cells.
#' If \code{gam} is set to 10, you will end up with about 1,000 supercells, whereas
#' if \code{gam} is set to 50, you will end up with about 200 supercells.
#'
#' Conversely, as you get more supercells (i.e. smaller \code{gam} value),
#' the smaller their size will be.
#' In other words, each of them will be capturing less cells.
#' Importantly, one cannot expect all the supercells to be of the same size.
#' Some will capture more/less cells than the other, and that is it not trivial
#' to estimate how many will be captured beforehand.
#' We may look into swapping the \code{gam} value to how many cells to be captured
#' within each supercell \emph{in the future}.
#' More thoughts are required into whether this make sense.
#'
#' Lastly, for now, a \code{gam} value for all your samples
#' (read the section above if you are not sure what I mean by samples here).
#' \emph{In the future}, we can perhaps look into setting different \code{gam} values
#' for different samples.
#'
#' @section Computing kNN network:
#' The parameter \code{k_knn} governs and the k value used to compute the
#' single-cell kNN network.
#' Actual (not approximate) kNN network is created, and walktrap algorithm is
#' used to form supercells from the kNN network.
#'
#' @return
#' \code{runSuperCellCyto} will return a list with the following components:
#' \describe{
#' \item{\code{supercell_object}:}{A list containing the object returned by
#' SCimplify function. One object per sample.}
#' \item{\code{supercell_expression_matrix}:}{A \link{data.table} containing
#' the marker expression of all the supercells.
#' These are computed by taking the average marker expression of the cells
#' captured by each supercell.}
#' \item{\code{supercell_cell_map}:}{A \link{data.table} showing which cell is
#' captured by which supercell. This is very useful if you intend to work out
#' which supercell captures which cell.}
#' }
#'
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
runSuperCellCyto <- function(dt,
                             markers,
                             sample_colname,
                             cell_id_colname,
                             gam = 20,
                             k_knn = 5,
                             n_parallel_worker = 1) {
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
    
    if (n_parallel_worker == 1) {
        # Load balancing is only needed if we are running processing samples in 
        # parallel, i.e. n_parallel_worker > 1
        BPPARAM <- SerialParam()
    } else {
        ncells_per_sample <- sapply(matrix_per_samp, ncol)
        names(ncells_per_sample) <- samples
        
        ncells_per_sample <- ncells_per_sample[order(ncells_per_sample, decreasing = TRUE)]
        matrix_per_samp <- matrix_per_samp[names(ncells_per_sample)]
        
        BPPARAM <- MulticoreParam(workers = n_parallel_worker, tasks = length(samples))
    }
    
    # Number of PCs are set to 10 by default. We can have panel size less than 10.
    # If this is the case, we just set PCA to be the number of markers
    if (length(markers) < 10) {
        n_pc <- length(markers)
    } else {
        n_pc <- 10
    }
    
    supercell_res <- bplapply(names(matrix_per_samp), function(sample_name, gam, k_knn) {
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
        
    }, gam = gam, k_knn = k_knn, BPPARAM = BPPARAM)

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

