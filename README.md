# SuperCellCyto 

<p>
<img src="man/figures/supercellcyto.png" align="right" height="2px" />
</p>

[SuperCellCyto](https://github.com/phipsonlab/SuperCellCyto) is an extension of the [SuperCell R package](https://github.com/GfellerLab/SuperCell). 
Initially developed for scRNAseq data, SuperCell aggregates cells with similar transcriptomic profiles into "supercells" (also known as “metacells” in the scRNAseq literature).

In SuperCellCyto, we've tailored the SuperCell package to specifically cater to cytometry data:

1. Supercells are now aggregating cells that are similar marker expressions.
2. Supercells are now created for each individual sample. 
This adaptation ensures that each supercell encompasses cells from exclusively one sample.
By processing each sample independently, we prevent the intermixing of cells from different samples within supercells.
3. Multiple samples can now be processed in parallel with a custom load balancing strategy.
This enhancement enables simultaneous generation of supercells for multiple samples, significantly reducing the computational time required for processing large datasets.

## Vignettes and Function Documentation

Vignettes and function documentations are available on our website: https://phipsonlab.github.io/SuperCellCyto/.

To understand how SuperCellCyto can be integrated into your workflow, head over to the Articles page.
Click on the dropdown arrow of the `Articles` link in the navbar at the top of the website.
Our vignettes provide step-by-step guides, practical examples, and use-case scenarios that demonstrate the package's application in various research contexts.

Documentation and usage of each function in the SuperCellCyto package 
are available in the Reference page. 
Click on the `Reference` link in the navbar at the top of the website.

## Citation

If you use SuperCellCyto in your study, please kindly cite our manuscript:

Putri, G.H., Howitt, G., Marsh-Wakefield, F. et al. SuperCellCyto: enabling efficient analysis of large scale cytometry datasets. Genome Biol 25, 89 (2024). https://doi.org/10.1186/s13059-024-03229-3

## Installation

The package can be installed using `remotes`:

```
# Install remotes
install.packages("remotes")

# Install SuperCellCyto from this repository
remotes::install_github("phipsonlab/SuperCellCyto")
```

## Code of Conduct

Please note that the SuperCellCyto project is released with a [Contributor Code of Conduct](https://phipsonlab.github.io/SuperCellCyto/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.

## Contribution guide

Please visit [Contributing Guide](https://phipsonlab.github.io/SuperCellCyto/CONTRIBUTING.html)





