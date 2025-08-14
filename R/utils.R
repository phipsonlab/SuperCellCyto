#' Create a mapping between supercells and cells for a sample.
#'
#' Internal function to generate a data.table mapping each cell to 
#' its supercell.
#'
#' @param sc_membership Integer vector. Supercell membership for each cell.
#' @param cell_ids Character vector. Cell IDs for each cell in sc_membership.
#' @param sample_name Character. Sample name.
#'
#' @return A data.table with columns: SuperCellID, CellId, Sample.
#' 
#' @import data.table
#' 
#' @keywords internal
#' @noRd 
#' 
.create_supercell_cell_map <- function(sc_membership, cell_ids, sample_name) {
    data.table(
        SuperCellID = paste0(
            "SuperCell_", sc_membership,
            "_Sample_", sample_name
        ),
        CellId = cell_ids,
        Sample = sample_name
    )
}


#' Extracts a marker expression matrix for a sample.
#'
#' Internal function to retrieve the marker expression values of 
#' for a given sample.
#'
#' @param dt A data.table containing cell-level data.
#' @param sample_name Character. The name of the sample to filter.
#' @param markers Character vector. The marker columns to extract.
#' @param sample_colname Character. 
#' The column name in `dt` that identifies the sample.
#' @param cell_id_colname Character. 
#' The column name in `dt` that identifies the cell.
#'
#' @return A matrix with rows as cells and columns as markers for 
#' the specified sample.
#' 
#' @keywords internal
#' @import data.table
#' @import Matrix
#' 
#' @noRd 
#'
.get_sample_marker_matrix <- function(
    dt, sample_name, markers, sample_colname, cell_id_colname
) {
    mat <- dt[dt[[sample_colname]] == sample_name, ]
    cell_id <- mat[[cell_id_colname]]
    mat <- Matrix(t(mat[, markers, with = FALSE]))
    colnames(mat) <- cell_id
    
    return(mat)
}

#' Compute centroid for supercells
#'
#' Internal helper to compute centroids for supercells using mean or median.
#'
#' @param mt Marker expression matrix for a sample.
#' @param sc_membership Integer vector. Supercell membership for each cell.
#' @param aggregation_method "mean" or "median".
#' @param sample_colname Name of the sample column.
#' @param sample_name Name of the sample.
#' @param markers Character vector of marker names.
#'
#' @return data.table with aggregated marker expression per supercell.
#' @keywords internal
#' @noRd
#'
#' @import data.table
#' @importFrom SuperCell supercell_GE
#' @importFrom stats median
#'
#'
.compute_supercell_centroid <- function(
    mt, 
    sc_membership, 
    aggregation_method, 
    sample_colname, 
    sample_name, 
    markers
) {
    if (aggregation_method == "mean") {
        supercell_exp_mat <- data.table(t(
            as.matrix(supercell_GE(ge = mt, groups = sc_membership))
        ))
        supercell_exp_mat[[sample_colname]] <- sample_name

        # Create a unique supercell id concatenating the sample name
        supercell_exp_mat[["SuperCellId"]] <- paste0(
            "SuperCell_",
            seq(1, nrow(supercell_exp_mat)),
            "_Sample_",
            sample_name
        )
    } else if (aggregation_method == "median") {
        supercell_exp_mat <- data.table(t(as.matrix(mt)))
        supercell_exp_mat$cell_id <- colnames(mt)
        supercell_membership <- data.table(
            cell_id = names(sc_membership),
            SuperCellId = paste0(
                "SuperCell_", sc_membership, "_Sample_", sample_name
            )
        )
        supercell_exp_mat <- merge.data.table(
            supercell_exp_mat,
            supercell_membership,
            by = "cell_id"
        )
        supercell_exp_mat <- supercell_exp_mat[
            , lapply(.SD, median),
            .SDcols = markers,
            by = "SuperCellId"
        ]
    }
    return(supercell_exp_mat)
}

#' Combine results from supercell calculations into a single list
#'
#' Internal helper to reshape the output of the per-sample supercell 
#' recalculation into a single list containing the combined
#' supercell expression matrix and supercell cell map.
#'
#' @param supercells List of per-sample recalculated supercells.
#' Each element should contain `supercell_expression_matrix` containing the
#' recalculated supercells centroid, a `supercell_cell_map` which maps 
#' each cell to its corresponding supercell.
#' 
#' @return List with combined supercell_expression_matrix and 
#' supercell_cell_map.
#' 
#' @keywords internal
#' @noRd 
#' 
.reshape_supercell_output <- function(supercells) {
    list(
        supercell_expression_matrix = do.call(
            rbind,
            lapply(supercells, function(res_i) {
                res_i$supercell_expression_matrix
            })
        ),
        supercell_cell_map = do.call(
            rbind,
            lapply(supercells, function(res_i) {
                res_i$supercell_cell_map
            })
        )
    )
}
