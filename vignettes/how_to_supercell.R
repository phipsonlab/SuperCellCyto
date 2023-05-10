## ---- include = FALSE, echo=FALSE, results="hide", message=FALSE--------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE, message=FALSE-----------------------------------------
library(SuperCellCyto)
library(parallel)

## -----------------------------------------------------------------------------
n_markers <- 15
n_samples <- 3
dat <- SuperCellCyto::simCytoData(nmarkers = n_markers, ncells = rep(10000, n_samples))
head(dat)

## -----------------------------------------------------------------------------
# Specify which columns are the markers to transform
marker_cols <- paste0("Marker_", seq_len(n_markers))
# The co-factor for arc-sinh
cofactor <- 5

# Do the transformation
dat_asinh <- asinh(dat[, marker_cols, with = FALSE] / cofactor)

# Rename the new columns
marker_cols_asinh <- paste0(marker_cols, "_asinh")
names(dat_asinh) <- marker_cols_asinh

# Add them our previously loaded data
dat <- cbind(dat, dat_asinh)

head(dat[, marker_cols_asinh, with = FALSE])

## -----------------------------------------------------------------------------
head(dat$Cell_Id, n = 10)

## -----------------------------------------------------------------------------
dat$Cell_id_dummy <- paste0("Cell_", seq_len(nrow(dat)))
head(dat$Cell_id_dummy, n = 10)

## -----------------------------------------------------------------------------
unique(dat$Sample)

## -----------------------------------------------------------------------------
table(dat$Sample)

## -----------------------------------------------------------------------------
markers_col <- paste0("Marker_", seq_len(n_markers), "_asinh")
sample_col <- "Sample"
cell_id_col <- "Cell_Id_dummy"

## -----------------------------------------------------------------------------
supercells <- runSuperCellCyto(
  dt = dat,
  markers = markers_col,
  sample_colname = sample_col,
  cell_id_colname = cell_id_col
)

## -----------------------------------------------------------------------------
class(supercells)

## -----------------------------------------------------------------------------
names(supercells)

## -----------------------------------------------------------------------------
head(supercells$supercell_expression_matrix)

## -----------------------------------------------------------------------------
names(supercells$supercell_expression_matrix)

## -----------------------------------------------------------------------------
head(unique(supercells$supercell_expression_matrix$SuperCellId))

## -----------------------------------------------------------------------------
supercell_ids <- unique(supercells$supercell_expression_matrix$SuperCellId)
supercell_ids[grep("SuperCell_1_", supercell_ids)]

## -----------------------------------------------------------------------------
head(supercells$supercell_cell_map)

## -----------------------------------------------------------------------------
n_cores <- detectCores()

