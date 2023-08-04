# Test the correctness of runSuperCellCyto

test_that("Output is a list", {
    cyto_dat <- simCytoData(10, rep(1000, 2))

    out <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(10)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id"
    )

    expect_equal(length(out), 3)
    expect_equal(class(out), "list")
    expect_true(all(
        c(
            "supercell_object",
            "supercell_expression_matrix",
            "supercell_cell_map"
        ) %in% names(out)
    ))
    expect_equal(class(out$supercell_object), "list")
})

test_that("Output cell mapping is correct", {
    ncells <- 1000
    nsample <- 3
    cyto_dat <- simCytoData(10, rep(ncells, nsample))

    out <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(10)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id"
    )

    expect_equal(nrow(out$supercell_cell_map), nrow(cyto_dat))
    expect_true(all(
        paste0("Cell_", seq_len(ncells * nsample)) %in% out$supercell_cell_map$CellId
    ))
    expect_true(all(
        paste0("Sample_", seq_len(nsample)) %in% unique(out$supercell_cell_map$Sample)
    ))
})

test_that("Serial and Parallel execution yields the same result", {
    cyto_dat <- simCytoData(10, c(10000, 200, 30000, 400, 5000))

    out_serial <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(10)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id"
    )

    out_parallel <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(10)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        BPPARAM = MulticoreParam(
            workers = parallel::detectCores() - 1, 
            tasks = 5
        )
    )

    # Check expression matrix
    out_parallel_exp_mat <- out_parallel$supercell_expression_matrix
    out_serial_exp_mat <- out_serial$supercell_expression_matrix
    # Just so we have both matrices in the same order
    out_parallel_exp_mat <- out_parallel_exp_mat[order(match(SuperCellId, out_serial_exp_mat$SuperCellId))]
    
    for (col_name in names(out_parallel_exp_mat)) {
        expect_true(identical(out_parallel_exp_mat[[col_name]], out_serial_exp_mat[[col_name]]))
    }
    
    # Check the cell supercell mapping
    out_parallel_mapping <- out_parallel$supercell_cell_map
    out_serial_mapping <- out_serial$supercell_cell_map
    # Just so we have both tables in the same order
    out_parallel_mapping <- out_parallel_mapping[order(match(CellId, out_serial_mapping$CellId))]
    expect_true(identical(out_parallel_mapping$SuperCellID, out_serial_mapping$SuperCellID))
    
})

test_that("Data with small number of markers can still be processed", {
    nmarkers <- 7
    cyto_dat <- simCytoData(nmarkers = nmarkers)

    expect_error(
        runSuperCellCyto(
            dt = cyto_dat,
            markers = paste0("Marker_", seq_len(nmarkers)),
            sample_colname = "Sample",
            cell_id_colname = "Cell_Id"
        ),
        NA
    )
})

test_that("Set seed is not required for reproducibility", {
    nmarkers <- 10
    cyto_dat <- simCytoData(nmarkers = nmarkers)

    # Serial execution
    run1_serial <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(nmarkers)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id"
    )

    run2_serial <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(nmarkers)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id"
    )

    expect_true(
        all.equal(
            run1_serial$supercell_expression_matrix,
            run2_serial$supercell_expression_matrix
        )
    )

    expect_true(
        all.equal(
            run1_serial$supercell_cell_map,
            run2_serial$supercell_cell_map
        )
    )

    # Parallel execution
    run1_parallel <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(nmarkers)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        BPPARAM = MulticoreParam(
            workers = parallel::detectCores() - 1, 
            tasks = 2
        )
    )

    run2_parallel <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(nmarkers)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        BPPARAM = MulticoreParam(
            workers = parallel::detectCores() - 1, 
            tasks = 2
        )
    )

    expect_true(
        all.equal(
            run1_parallel$supercell_expression_matrix,
            run2_parallel$supercell_expression_matrix
        )
    )

    expect_true(
        all.equal(
            run1_parallel$supercell_cell_map,
            run2_parallel$supercell_cell_map
        )
    )
})

test_that("List containing supercell objects are ordered correctly", {
    cyto_dat <- simCytoData(ncells = c(1000, 30000, 20000, 200))
    
    sc <- runSuperCellCyto(
        dt = cyto_dat,
        markers = paste0("Marker_", seq_len(10)),
        sample_colname = "Sample",
        cell_id_colname = "Cell_Id",
        BPPARAM = MulticoreParam(
            workers = parallel::detectCores() - 1, 
            tasks = 4
        )
    )
    
    samples <- unique(cyto_dat$Sample)
    
    for (s in samples) {
        membership_diff <- union(
            setdiff(names(sc$supercell_object[[s]]$membership), cyto_dat[Sample == s,]$Cell_Id),
            setdiff(cyto_dat[Sample == s,]$Cell_Id, names(sc$supercell_object[[s]]$membership))
        )
        expect_true(length(membership_diff) == 0)
    } 
    
    
})



