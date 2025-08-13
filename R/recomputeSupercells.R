#' Recompute supercells
#'
#' @description
#' Given a supercell object, recreate the supercells using a different
#' gamma value.
#'
#' Gamma value controls the number of supercells generated.
#' The smaller the value, the more supercells you get, and vice versa.
#'
#' For this function to run, you need to have at least run
#' [runSuperCellCyto()] function **once**!
#'
#' @param dt A \pkg{data.table} object containing cytometry data where rows
#' represent cells and columns represent markers.
#' @param sc_objects The `supercell_object` returned by
#' [runSuperCellCyto()] function.
#' @param markers A character vector identifying the markers to create
#' supercells with.
#' @param sample_colname A character string identifying the column in
#' \code{dt} that denotes the sample of a cell.
#' @param cell_id_colname A character string identifying the column in
#' \code{dt} representing each cell's unique ID.
#' @param gam A numeric value specifying the gamma value which regulates
#' the number of supercells generated.
#' Defaults to 20.
#'
#' @return
#' A list with the following components:
#'
#' * `supercell_expression_matrix`:  A \pkg{data.table} object that contains
#' the marker expression for each supercell.
#' These marker expressions are computed by calculating the mean of the marker
#' expressions across all cells
#' within each individual supercell.
#' * `supercell_cell_map`: A \pkg{data.table} that maps each cell to its
#' corresponding supercell.
#' This table is essential for identifying the specific supercell each cell has
#' been allocated to.
#' It proves particularly useful for analyses that require one to expand the
#' supercells to the individual cell level.
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

        exp_mat <- .get_sample_marker_matrix(
            dt = dt,
            sample_name = sample_name,
            markers = markers,
            sample_colname = sample_colname,
            cell_id_colname = cell_id_colname
        )
        cell_ids <- colnames(exp_mat)

        supercell_exp_mat <- .recompute_supercell_centroid(
            mat = exp_mat,
            sample_name = sample_name,
            sc_membership = sc_rescaled$membership,
            sample_colname = sample_colname
        )

        # ---- Create supercell and cell mapping ----
        supercell_cell_map <- .create_supercell_cell_map(
            sc_membership = sc_rescaled$membership,
            cell_ids = cell_ids,
            sample_name = sample_name
        )

        # Return a list containing all the objects
        return(list(
            supercell_expression_matrix = supercell_exp_mat,
            supercell_cell_map = supercell_cell_map
        ))
    })

    # Now the messy reshaping so each element is not the output for a sample
    # but either a supercell object, expression matrix or supercell cell map
    reshaped_res <- .combine_supercell_results(supercell_rescaled)

    return(reshaped_res)
}
#' Recompute the centroid of a supercell
#'
#' Internal function to recalculate the centroid for a given supercell 
#' based on its constituent cells.
#' 
#' @param mat A matrix containing the marker expression values for the cells
#' in a sample.
#' @param sample_name A character string identifying the sample for which
#' the supercell centroid is being recalculated.
#' This should match the sample names in the `dt` data.table.
#' @param sc_membership The membership information for the supercell, i.e.,
#' which cells belong to which supercell.
#' This is typically obtained from after running [supercell_rescale()] function.
#' @param sample_colname A character string identifying the column in
#' \code{dt} that denotes the sample of a cell.
#' 
#' @return A \pkg{data.table} object containing the recalculated centroid for
#' the supercell.
#' 
#' @keywords internal
#' 
#' @import data.table
#' @importFrom SuperCell supercell_GE
#' @importFrom Matrix Matrix
#'
.recompute_supercell_centroid <- function(
    mat, sample_name, sc_membership, sample_colname
) {
    # Calculate supercell expression matrix
    supercell_exp_mat <- data.table(t(
        as.matrix(supercell_GE(ge = mat, groups = sc_membership))
    ))
    supercell_exp_mat[[sample_colname]] <- sample_name

    # Create a unique supercell id concatenating the sample name
    supercell_exp_mat[["SuperCellId"]] <- paste0(
        "SuperCell_",
        seq(1, nrow(supercell_exp_mat)), "_Sample_", sample_name
    )

    return(supercell_exp_mat)
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
#' @keywords internal
#' @import data.table
#' @import Matrix
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
#' @keywords internal
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

#' Combine results from recalculated supercells into a single list
#'
#' Internal helper to reshape the output of the per-sample supercell 
#' recalculation into a single list containing the combined
#' supercell expression matrix and supercell cell map.
#'
#' @param supercell_rescaled List of per-sample recalculated supercells.
#' Each element should contain `supercell_expression_matrix` containing the
#' recalculated supercells centroid, a `supercell_cell_map` which maps 
#' each cell to its corresponding supercell.
#' 
#' @return List with combined supercell_expression_matrix and 
#' supercell_cell_map.
#' 
#' @keywords internal
#' 
.combine_supercell_results <- function(supercell_rescaled) {
    list(
        supercell_expression_matrix = do.call(
            rbind,
            lapply(supercell_rescaled, function(res_i) {
                res_i$supercell_expression_matrix
            })
        ),
        supercell_cell_map = do.call(
            rbind,
            lapply(supercell_rescaled, function(res_i) {
                res_i$supercell_cell_map
            })
        )
    )
}


