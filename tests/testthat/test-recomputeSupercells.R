test_that("recompute give the same results", {
    cyto_dat <- simCytoData(10, rep(1000, 3))

    markers <- paste0("Marker_", seq_len(10))

    out_gam20 <- runSuperCellCyto(
        dt = cyto_dat,
        markers = markers,
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        gam = 20
    )

    out_gam50 <- runSuperCellCyto(
        dt = cyto_dat,
        markers = markers,
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        gam = 50
    )

    recomputed_sc <- recomputeSupercells(
        dt = cyto_dat,
        sc_objects = out_gam20$supercell_object,
        markers = markers,
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        gam = 50
    )

    expected_mat <- out_gam50$supercell_expression_matrix
    given_mat <- recomputed_sc$supercell_expression_matrix

    expected_map <- out_gam50$supercell_cell_map
    given_map <- recomputed_sc$supercell_cell_map

    expect_equal(nrow(expected_mat), nrow(given_mat))
    supercell_ids <- expected_mat$SuperCellId

    for (sc in supercell_ids) {
        expected_marker <- expected_mat[SuperCellId == sc, markers, with = FALSE]
        given_marker <- given_mat[SuperCellId == sc, markers, with = FALSE]
        expect_identical(expected_marker, given_marker)

        # test membership
        expected_membership <- expected_map[SuperCellID == sc, ]$CellId
        given_membership <- given_map[SuperCellID == sc, ]$CellId
        expect_identical(expected_membership, given_membership)
    }

})


