# Test the input format for runSuperCellCyto

test_that("data.frame is converted to data.table", {
    cyto_dat <- as.data.frame(simCytoData(10, rep(1000, 2)))

    expect_warning(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(10)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id"
        ),
        "dt is not a data.table object. Converting it to a data.table object"
    )
})

test_that("error is sent if dt is not at least a data.frame", {
    cyto_dat <- matrix(c(1, 2, 3, 11, 12, 13))

    expect_error(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(3)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id"
        )
    )
})
