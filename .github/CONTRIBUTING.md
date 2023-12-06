# Contributing to SuperCellCyto

This outlines how to propose a change to SuperCellCyto. 

## Making changes to the package

If you want to make a change, please first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer.

*   Install all development dependencies with `devtools::install_dev_deps()`.
*   Create a Git branch for your pull request (PR). You can use `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git.
*   Make sure all the unit tests are passing. You can use `devtools::test()`.
*   Create a PR. Make sure that the title of your PR briefly describe the change, and the body of your PR should contain `Fixes #issue-number`.

## Code of Conduct

Please note that the SuperCellCyto project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
