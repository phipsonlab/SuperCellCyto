# README

Data in `extdata` folder are obtained from the following sources:

* `Data23_Panel3_base_NR4_Patient9.fcs` and 
`Data23_Panel3_base_R5_Patient15.fcs` are mass cytometry data from anti-PD1 
  immunotherapy study by [Krieg et. al](https://doi.org/10.1038/nm.4466).
  FCS files are available from 
  [FlowRepository](
  http://flowrepository.org/public_experiment_representations/1124).
* `Levine_32dim_H1_sub.csv` and `Levine_32dim_H2_sub.csv` are mass cytometry
  data from AML study by 
  [Levine et. al](https://doi.org/10.1016/j.cell.2015.05.047).
* `Levine_32dim_seurat_sub.qs2` and `Levine_32dim_sce_sub.qs2` are
  Seurat object and SingleCellExperiment object containing subsampled
  mass cytometry data from AML study by 
  [Levine et. al](https://doi.org/10.1016/j.cell.2015.05.047).
  
Scripts to generate the different formats of `Levine_32dim` dataset is stored 
in the `inst/scripts` directory of 
[SuperCellCyto](https://github.com/phipsonlab/SuperCellCyto) repo.
