# SuperCellCyto
 
SuperCellCyto is an R package to summarise/reduce the size of cytometry data in a more meaningful way using the [SuperCell R package](https://github.com/GfellerLab/SuperCell) by David Gfeller lab from the University of Lausanne. 

This package is still under active development, and we welcome all feedbacks. 
If you have suggestions on how to improve the package or any questions on how to use it, please create a github issue. 

For now, to install the package, please use:

```
devtools::install_github("phipsonlab/SuperCellCyto")
```

SuperCellCyto imports the [SuperCell R package](https://github.com/GfellerLab/SuperCell)
which is only available from Github. 
The installation should automatically take care of this, but in case it does not,
please install it using the following command:

```
devtools::install_github("GfellerLab/SuperCell")
```

# Citations

If you use SuperCellCyto in your publication please cite the original this repository (1st code block) and the original SuperCell paper (2nd code block).

This repository:

[![DOI](https://zenodo.org/badge/592222314.svg)](https://zenodo.org/badge/latestdoi/592222314)

```
@software{supercellcyto,
  author = {Putri, Givanna and Howitt, George and Ashhurst Thomas and Phipson, Belinda},
  doi = {10.5281/zenodo.7601561},
  month = {2},
  title = {{SuperCellCyto}},
  url = {https://github.com/phipsonlab/SuperCellCyto},
  version = {0.99.0},
  year = {2023}
}
```


SuperCell paper:

```
@article{bilous2022metacells,
  title={Metacells untangle large and complex single-cell transcriptome networks},
  author={Bilous, Mariia and Tran, Loc and Cianciaruso, Chiara and Gabriel, Aur{\'e}lie and Michel, Hugo and Carmona, Santiago J and Pittet, Mikael J and Gfeller, David},
  journal={BMC bioinformatics},
  volume={23},
  number={1},
  pages={336},
  year={2022},
  publisher={Springer}
}
```
