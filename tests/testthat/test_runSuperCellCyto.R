test_that("Output is a list", {
    cyto_dat <- simCytoData(10, 1000, 2)
    
    out <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(10)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id"
    )
    
    expect_equal(length(out), 3)
    expect_equal(class(out), "list")
    expect_true(all(
        c("supercell_object", 
          "supercell_expression_matrix", 
          "supercell_cell_map") %in% names(out)))
    expect_equal(class(out$supercell_object), "list")
    
})

test_that("Output cell mapping is correct", {
    ncells <- 1000
    nsample <- 3
    cyto_dat <- simCytoData(10, ncells, nsample)
    
    out <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(10)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id"
    )
    
    expect_equal(nrow(out$supercell_cell_map), nrow(cyto_dat))
    expect_true(all(
        paste0("Cell_", seq_len(ncells*nsample)) %in% out$supercell_cell_map$CellId
    ))
    expect_true(all(
        paste0("Sample_", seq_len(nsample)) %in% unique(out$supercell_cell_map$Sample)
    ))
    
})

test_that("Serial and Parallel execution yields the same result", {
    ncells <- 1000
    nsample <- 3
    cyto_dat <- simCytoData(10, ncells, nsample)
    
    out_serial <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(10)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id"
    )
    
    BPPARAM <- BiocParallel::MulticoreParam()
    out_parallel <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(10)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id",
        BPPARAM=BPPARAM
    )
    
    expect_identical(out_serial$SuperCellID, out_parallel$SuperCellID)
    
})

test_that("Data with small number of markers can still be processed", {
    nmarkers <- 7
    cyto_dat <- simCytoData(nmarkers=nmarkers)
    
    expect_error(
        runSuperCellCyto(
            dt=cyto_dat,
            markers=paste0("Marker_", seq_len(nmarkers)),
            sample_colname="Sample",
            cell_id_colname="Cell_Id"
        ),
        NA
    )
})

test_that("Set seed is not required for reproducibility", {
    nmarkers <- 10
    cyto_dat <- simCytoData(nmarkers=nmarkers)
    
    run1 <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(nmarkers)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id"
    )
    
    run2 <- runSuperCellCyto(
        dt=cyto_dat,
        markers=paste0("Marker_", seq_len(nmarkers)),
        sample_colname="Sample",
        cell_id_colname="Cell_Id"
    )
    
    expect_true(
        all.equal(
            run1$supercell_expression_matrix,
            run2$supercell_expression_matrix
        )
    )
    
    expect_true(
        all.equal(
            run1$supercell_cell_map,
            run2$supercell_cell_map
        )
    )
    
    
})