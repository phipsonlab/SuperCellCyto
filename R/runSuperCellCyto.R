#' Run SuperCell for cytometry data
#'
#' Run SuperCell on cytometry data stored as \link{data.table} object.
#' The vanilla SuperCell algorithm is provided by the SuperCell R package 
#' (Bilous et.al, 2022). 
#' We have enhanced it by adding the capacity to "supercell" multiple samples 
#' in parallel
#' through the use of \link{BiocParallel}, and by adding support to "supercell"
#' cytometry data. 
#' See Details below on how SuperCell was used.
#'
#' @param dt A \link{data.table} object containing the cytometry data.
#' Rows represent cells, columns represent markers.
#' @param markers A character vector specifying the markers to run SuperCell on.
#' @param sample_colname String specifying the column in \code{dt} that denotes
#' the sample of a cell.
#' @param cell_id_colname String specifying the column in \code{dt} that denotes
#' the unique ID of a cell.
#' @param gam Numeric scalar specifying the gamma value to be used by SuperCell.
#' @param k_knn Numeric scalar specifying the k value to be used by SuperCell's knn.
#' @param BPPARAM \linkS4class{BiocParallelParam} object specifying the 
#' configuration parameters for parallel execution.
#' Default to \linkS4class{SerialParam}, i.e., not parallelisation to be used.
#' 
#' #' @section Details about SuperCell:
#' We used SuperCell's SCimplify method to generate supercells for cytometry data.
#' To install SuperCell package, please use \link[remotes]{install_github} as 
#' it is only available on github:
#' \code{remotes::install_github("GfellerLab/SuperCell")}.
#' 
#' By default, all the markers specified in \code{markers} parameter are used 
#' to compute PCA,and that the marker expressions are \emph{not} scaled when 
#' computing PCA.
#' \code{irlba} is not used to calculate PCA as cytometry data only have a 
#' handful of features (markers) in general.
#' Number of PCs are set to 10, default in SCimplify.
#' 
#' \code{gam} and \code{k_knn} are passed on as it is to indicate the graining 
#' level of supercells and the k value used to compute the single-cell 
#' kNN network.
#' Actual (not approximate) kNN network is created, and walktrap algorithm was 
#' used to detect the supercells from the kNN network. 
#' 
#' Cells from each sample are processed independent of those from different samples. 
#' To differentiate the same supercell ID across multiple samples, 
#' the supercell IDs are prepended with the sample ID.
#' For instance, supercell 1 from sample 1 will be named SuperCell_1_Sample_1, 
#' whereas supercell 1 from sample 2 will be named SuperCell_2_Sample_1.
#' 
#' If none of the above make sense to you, please read Bilous et.al, 2022 
#' manuscript to see how SuperCell works.
#' 
#' @return
#' \code{runSuperCellCyto} will return a list with the following components:
#' \describe{
#' \item{\code{supercell_object}:}{A list containing a list returned by 
#' SCimplify function for each sample.}
#' \item{\code{supercell_expression_matrix}:}{A \link{data.table} containing 
#' the marker expression of all the supercells.}
#' \item{\code{supercell_cell_map}:}{A \link{data.table} containing the 
#' supercell ID of all the cells in the data.}
#' }
#' 
#' Each supercell's marker expression in \code{supercell_expression_matrix} is 
#' computed based on the average marker expression of all the cells captured 
#' by the supercell.
#'
#' @author
#' Givanna Putri
#' 
#' @examples 
#' # Simulate some data
#' set.seed(42)
#' cyto_dat <- simCytoData(nmarkers=10, ncells=2000, nsample=2)
#' 
#' # Setup the columns designating the markers, samples, and cell IDs
#' marker_col <- paste0("Marker_", seq_len(10))
#' sample_col <- "Sample"
#' cell_id_col <- "Cell_Id"
#' 
#' supercell_dat <- runSuperCellCyto(cyto_dat, marker_col, 
#' sample_col, cell_id_col)
#'
#'
#' @export
#' 
#' @import data.table
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom Matrix Matrix
#' @importFrom SuperCell SCimplify supercell_GE
runSuperCellCyto <- function(
        dt,
        markers,
        sample_colname,
        cell_id_colname,
        gam=20,
        k_knn=5,
        BPPARAM=SerialParam()
) {
    
    samples <- unique(dt[[sample_colname]])
    
    matrix_per_samp <- lapply(samples, function(samp) {
        dt_sub <- dt[dt[[sample_colname]] == samp, ]
        trans_dt_sub <- Matrix(t(dt_sub[, markers, with=FALSE]))
        colnames(trans_dt_sub) <- dt_sub[[cell_id_colname]]
        return(trans_dt_sub)
    })
    
    # Number of PCs are set to 10 by default. We can have panel size less than 10.
    # If this is the case, we just set PCA to be the number of markers
    n_markers <- length(markers)
    if (n_markers < 10) {
        n_pc <- n_markers
    } else {
        n_pc <- 10
    }
    
    supercell_res <- bplapply(matrix_per_samp, function(mt, seed, gam, k_knn) {
        # Note to self: there is no need to pass on the seed as it is only
        # used to subsample cells when approximate kNN is required, 
        # i.e., when do.approx is set to TRUE.
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
        return(res)
    }, gam=gam, k_knn=k_knn, BPPARAM=BPPARAM)
    
    # Get the expression matrix for each supercell
    # A bit weird how we need to create a list of list so we can pass the data
    # matrix, sample name, and supercell object into bplapply.
    # Can't think of any other less ugly way
    supercell_exp_materials <- lapply(seq_along(matrix_per_samp), function(i) {
        combo <- list(
            mtx=matrix_per_samp[[i]],
            sample=samples[i],
            supercell_obj=supercell_res[[i]]
        )
        return(combo)
    })
    
    # See above on what is in each "material"
    supercell_exp_mat <- bplapply(supercell_exp_materials, function(material) {
        mat <- material$mtx
        samp <- material$sample
        supercell_obj <- material$supercell_obj
        
        supercell_exp <- data.table(
            t(
                as.matrix(
                    supercell_GE(mat, supercell_obj$membership)
                )
            )
        )
        supercell_exp[[sample_colname]] <- samp
        
        # Have to create a unique id concatenating the sample as well
        supercell_exp[['SuperCellId']] <- paste0("SuperCell_", 
                seq_len(nrow(supercell_exp)), "_Sample_", samp)
        
        return(supercell_exp)
    }, BPPARAM=BPPARAM)
    supercell_exp_mat <- rbindlist(supercell_exp_mat)
    
    # Mapping between supercell ID and actual cell
    # Could've combined this with the bplapply above, but it can be messy
    # to detangle the result
    supercell_cell_map <- lapply(supercell_exp_materials, function(material) {
        mat <- material$mtx
        samp <- material$sample
        supercell_obj <- material$supercell_obj
        
        res <- data.table(
            SuperCellID = paste0("SuperCell_", supercell_obj$membership, 
                            "_Sample_", samp),
            CellId = colnames(mat),
            Sample = samp
        )
        return(res)
    })
    supercell_cell_map <- rbindlist(supercell_cell_map)
    
    res <- list(
        supercell_object=supercell_res,
        supercell_expression_matrix=supercell_exp_mat,
        supercell_cell_map=supercell_cell_map
    )
    return(res)
}


