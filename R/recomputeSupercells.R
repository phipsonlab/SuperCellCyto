#' Recompute supercells
#' 
#' @description
#' Given a supercell object, recreate the supercells using a different
#' gamma value.
#' 
#' Gamma value controls the number of supercells generated.
#' The smaller the value, the more supercells you get, and vice versa.
#' 
#' For this function to run, you need to have at least run \link{runSuperCellCyto}
#' function **once**!
#' 
#' @param dt A \link{data.table} object containing cytometry data where rows represent 
#' cells and columns represent markers.
#' @param sc_objects The `supercell_object` returned by \link{runSuperCellCyto} function.
#' @param markers A character vector identifying the markers to create supercells with.
#' @param sample_colname A character string identifying the column in \code{dt} that denotes the sample of a cell.
#' @param cell_id_colname A character string identifying the column in \code{dt} representing each cell's unique ID.
#' @param gam A numeric value specifying the gamma value which regulates the number of supercells generated.
#' Defaults to 20.
#'
#' @return
#' A list with the following components:
#' 
#' * `supercell_expression_matrix`:  A \link{data.table} object that contains the marker expression for each supercell.
#' These marker expressions are computed by calculating the mean of the marker expressions across all cells
#' within each individual supercell.
#' * `supercell_cell_map`: A \link{data.table} that maps each cell to its corresponding supercell. 
#' This table is essential for identifying the specific supercell each cell has been allocated to. 
#' It proves particularly useful for analyses that require one to expand the supercells
#' to the individual cell level.
#' 
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
#' 
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


