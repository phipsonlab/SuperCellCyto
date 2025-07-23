# Test runSuperCellCyto can run

test_that("runSuperCellCyto serial works", {
    cyto_dat <- simCytoData(10, c(10000, 200, 30000, 400, 5000))

    expect_no_error(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(10)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id"
        )
    )
})

test_that("runSuperCellCyto parallel without load balancing works", {
    cyto_dat <- simCytoData(10, c(10000, 200, 30000, 400, 5000))

    expect_no_error(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(10)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id",
            BPPARAM = MulticoreParam(tasks = 5)
        )
    )
})

test_that("runSuperCellCyto parallel with load balancing works", {
    cyto_dat <- simCytoData(10, c(10000, 200, 30000, 400, 5000))

    expect_no_error(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(10)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id",
            load_balancing = TRUE,
            BPPARAM = MulticoreParam(tasks = 5)
        )
    )
})
