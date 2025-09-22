library(HDCytoData)
library(data.table)
library(qs2)


sce <- Levine_32dim_SE()

# isolate data from H1 and H2 patients
sce_H1 <- sce[rowData(sce)$patient_id == "H1", ]
sce_H2 <- sce[rowData(sce)$patient_id == "H2", ]

# remove unassigned
sce_H1 <- sce_H1[rowData(sce_H1)$population_id != "unassigned", ]
sce_H2 <- sce_H2[rowData(sce_H2)$population_id != "unassigned", ]

# randomly select 200 cells from each patient
set.seed(42)
sce_H1_sub <- sce_H1[sample(nrow(sce_H1), 200), ]
sce_H2_sub <- sce_H2[sample(nrow(sce_H2), 200), ]

# convert to data.table
dt_H1_sub <- data.table(assay(sce_H1_sub, "exprs"))
dt_H2_sub <- data.table(assay(sce_H2_sub, "exprs"))

# save them as csv files
fwrite(dt_H1_sub, "inst/extdata/Levine_32dim_H1_sub.csv")
fwrite(dt_H2_sub, "inst/extdata/Levine_32dim_H2_sub.csv")


# create the sce object
dt <- data.table(assay(sce, "exprs"))
dt[, population_id := rowData(sce)$population_id]
dt[, sample := rowData(sce)$patient_id]

# get 50 cells from each population and sample
set.seed(42)
n_cells <- 50
sampled_dt <- dt[, .SD[sample(.N, min(.N, n_cells))], by = c("population_id", "sample")]

# assign cell id
sampled_dt[, cell_id := paste0("cell_", .I)]

# create SCE object
exprs_mat <- t(as.matrix(sampled_dt[, !c("population_id", "sample", "cell_id"), with = FALSE]))
sce_sampled <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = exprs_mat))
SummarizedExperiment::colData(sce_sampled)$population <- sampled_dt$population_id
SummarizedExperiment::colData(sce_sampled)$sample <- sampled_dt$sample
SummarizedExperiment::colData(sce_sampled)$cell_id <- sampled_dt$cell_id
colnames(sce_sampled) <- sampled_dt$cell_id

qs_save(sce_sampled, "inst/extdata/Levine_32dim_sce_sub.qs2")


# create Seurat object
seurat_obj <- Seurat::CreateSeuratObject(
    counts = exprs_mat,
    assay = "originalexp",
)
seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = sampled_dt[, .(population_id, sample, cell_id)])
# rename patient_id to sample


qs_save(seurat_obj, "inst/extdata/Levine_32dim_seurat_sub.qs2")



