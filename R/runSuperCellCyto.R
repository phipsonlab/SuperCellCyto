#' Run SuperCell for cytometry data
#'
#' Run SuperCell on cytometry data stored as \link{data.table} object.
#' This is a wrapper function around the `SCImplify` function in the 
#' SuperCell R package by Bilous et.al, 2022.
#' We have enhanced it by adding the capacity to "supercell" multiple samples 
#' in parallel through the use of \link{BiocParallel}, and by adding support to 
#' "supercell" cytometry data. 
#' More explanations are given in different sections below, and are somewhat
#' expanded in our vignette.
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
#' @section What is \code{cell_id_colname}:
#' This is a column in \code{dt} containing a unique identifier for each cell.
#' Commonly, you will have to manually create this column as FCS file does not
#' typically contain a field which can uniquely identify each cell.
#' You can do this quite easily by numbering the cells 1 to however many you have,
#' and store this in a column in \code{dt}.
#' If you don't know how to do this, refer to our vignette.
#' 
#' @section Processing one sample independent of the others:
#' This function is designed such that all the samples are processed in parallel,
#' independent of each other. 
#' Because of this, you can safely assume that each supercell will only contain
#' cells from exactly one sample. 
#' 
#' The function will work out which cells come from which sample based on what is
#' specified in the \code{sample_colname} column.
#' This is why it is critical to specify what this column is.
#' 
#' Here, a sample is basically a biological sample in your experiment.
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
#' Generally speaking, the smaller the \code{gam} parameter is, the more supercells
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
#' Lastly, for now, you can only set the same \code{gam} value for all your samples
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


