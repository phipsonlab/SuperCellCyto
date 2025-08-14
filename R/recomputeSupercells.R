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
#' @param aggregation_method A character string indicating how to aggregate
#' the cells in supercells to get expression matrix.
#' Options are "mean" or "median". Defaults to "mean".
#' If "mean", the mean of the marker expressions across all cells
#' within each individual supercell is computed.
#' If "median", the median of the marker expressions across all cells
#' within each individual supercell is computed.
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
                                aggregation_method = c("mean", "median"),
                                gam = 20) {

    samples <- names(sc_objects)

    # How to aggregate the cells in supercells to get expression matrix?
    aggregation_method <- match.arg(aggregation_method)

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
        supercell_exp_mat <- .compute_supercell_centroid(
            mt = exp_mat,
            sc_membership = sc_rescaled$membership,
            aggregation_method = aggregation_method,
            markers = markers,
            sample_colname = sample_colname,
            sample_name = sample_name
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

    # reshape the output into a single list
    reshaped_res <- .reshape_supercell_output(supercell_rescaled)

    return(reshaped_res)
}
