# SuperCellCyto 0.99.1

* Changed `paste` to `sprintf` for warning messages.
* Moved example data to `inst/extdata` and update vignettes.
* Add chunk labels to vignettes.

# SuperCellCyto 0.99.0

## Major changes
* Added support to use median when recomputing supercells.
* Introduced `n_pc` parameter to allow users to specify the number of 
principal components to use in the PCA step when computing supercells.
* Added `utils.R` file to contain utility functions that are commonly used.

## Minor changes
* Cleaned up functions and added improvements to the codebase, in accordance
with the Bioconductor package development guidelines.
* Updated vignettes to BiocStyle.
* Added more tests.

# SuperCellCyto 0.1.0

This is the first release of SuperCellCyto.

This is also the version of SuperCellCyto that was used in our preprint:

```
@article{
    putri2023supercellcyto,
    title={SuperCellCyto: enabling efficient analysis of large scale cytometry datasets},
    author={Putri, Givanna H and Howitt, George and Marsh-Wakefield, Felix and Ashhurst, Thomas Myles and Phipson, Belinda},
    journal={bioRxiv},
    pages={2023--08},
    year={2023},
    publisher={Cold Spring Harbor Laboratory}
}
```
