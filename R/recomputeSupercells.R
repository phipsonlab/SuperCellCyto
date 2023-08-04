#' Recompute supercells
#' 
#' Given a supercell object, recompute the supercell granularity using a different
#' gamma value.
#' Gamma value controls the number of supercells generated.
#' The smaller the value, the more supercells you get, and vice versa.
#'
#' @param dt \link{data.table} object containing the cytometry data.
#' Rows represent cells, columns represent markers.
#' If this is not a \link{data.table} object, the function will warn you about it,
#' and then try to convert it to a \link{data.table} object.
#' @param sc_objects The \code{supercell_object} returned by \link{runSuperCellCyto} function.
#' @param markers character vector specifying the markers in \code{dt}.
#' @param sample_colname character specifying the column in \code{dt} that denotes
#' the sample of a cell.
#' @param cell_id_colname character specifying the column in \code{dt} that denotes
#' the unique ID of a cell.
#' @param gam numeric specifying the gamma value to be used by SuperCell.
#'
#' @return
#' A list with the following components:
#' \describe{
#' \item{\code{supercell_expression_matrix}:}{A \link{data.table} containing
#' the marker expression of all the supercells.
#' These are computed by taking the average marker expression of the cells
#' captured by each supercell.}
#' \item{\code{supercell_cell_map}:}{A \link{data.table} showing which cell is
#' captured by which supercell. This is very useful if you intend to work out
#' which supercell captures which cell.}
#' }
#' @export
#'
#' @examples
#' set.seed(42)
#' cyto_dat <- simCytoData(10, rep(1000, 3))
#' markers <- paste0("Marker_", seq_len(10))
#' out_gam20 <- runSuperCellCyto(
#'     dt = cyto_dat,
#'     markers = markers,
#'     sample_colname = "Sample",
#'     cell_id_colname = "Cell_Id",
#'     gam = 20
#' )
#' recomputed_sc <- recomputeSupercells(
#'     dt = cyto_dat,
#'     sc_objects = out_gam20$supercell_object,
#'     markers = markers,
#'     sample_colname = "Sample",
#'     cell_id_colname = "Cell_Id",
#'     gam = 50
#' )
#' 
#' @author
#' Givanna Putri
#' 
#' @import data.table
#' @importFrom Matrix Matrix
#' @importFrom SuperCell supercell_rescale supercell_GE
recomputeSupercells <- function(dt,
                                sc_objects,
                                markers,
                                sample_colname,
                                cell_id_colname,
                                gam = 20) {
    
    samples <- names(sc_objects)
    
    supercell_rescaled <- lapply(samples, function(sample_name) {
        
        sc_rescaled <- supercell_rescale(sc_objects[[sample_name]], gamma = gam)
        
        # need this to recompute the supercell expression matrix
        mat <- dt[dt[[sample_colname]] == sample_name,]
        cell_id <- mat[[cell_id_colname]]
        mat <- Matrix(t(mat[, markers, with = FALSE]))
        colnames(mat) <- cell_id
        
        
        # ---- Calculate supercell expression matrix ----
        supercell_exp_mat <- data.table(
            t(
                as.matrix(
                    supercell_GE(
                        ge = mat, 
                        groups = sc_rescaled$membership
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
                "SuperCell_", sc_rescaled$membership,
                "_Sample_", sample_name
            ),
            CellId = colnames(mat),
            Sample = sample_name
        )
        
        # Return a list containing all the objects
        return(list(
            supercell_object = sc_rescaled,
            supercell_expression_matrix = supercell_exp_mat,
            supercell_cell_map = supercell_cell_map
        ))
    })
    
    # Now the messy reshaping so each element is not the output for a sample
    # but either a supercell object, expression matrix or supercell cell map
    reshaped_res <- list(
        supercell_expression_matrix = do.call(rbind, lapply(supercell_rescaled, function(res_i) res_i$supercell_expression_matrix)),
        supercell_cell_map = do.call(rbind, lapply(supercell_rescaled, function(res_i) res_i$supercell_cell_map))
    )
    
    return(reshaped_res)
    
}


