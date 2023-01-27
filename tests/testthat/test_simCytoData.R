test_that("Simulated data has correct number of cells and markers", {
    ncells <- 50000
    nmarkers <- 20
    nsamples <- 5
    
    sim_dat <- simCytoData(nmarkers=nmarkers, ncells=ncells, nsample=nsamples)
    
    # +2 because we should have sample and cell id column
    expect_equal(ncol(sim_dat), nmarkers + 2)
    expect_equal(nrow(sim_dat), ncells * nsamples)
    
})

test_that("Simulated data has correct headers", {
    nmarkers <- 20
    sim_dat <- simCytoData(nmarkers=nmarkers)
    
    expect_true("Sample" %in% names(sim_dat))
    expect_true("Cell_Id" %in% names(sim_dat))
    expect_true(all(
        paste0("Marker_", seq_len(nmarkers)) %in% names(sim_dat)
    ))
})

test_that("Simulated data has assigned correct cell ID", {
    ncells <- 100
    sim_dat <- simCytoData(ncells=ncells)
    
    expect_true(all(
        paste0("Cell_", seq_len(ncells)) %in% sim_dat$Cell_Id
    ))
})

test_that("Simulated data has assigned correct sample", {
    nsample <- 3
    sim_dat <- simCytoData(ncells=100, nsample=nsample)
    
    expect_true(all(
        paste0("Sample_", seq_len(nsample)) %in% unique(sim_dat$Sample)
    ))
})