# destiny

**NOTE (25/1/2024)**

In this fork I am implementing `BiocNeighbors` methods for KNN calculation, which run faster than the original `covertree` method, and offer a more flexible interface to several other exact (`kmknn`, `Vptree`) and approximate (`hsnw`, already implemented in the original `destiny` release, and `annoy`) methods. This fork also adds support for sparse matrices for KNN, including rank correlation distance, through `matrixStats`.

So far build and tests seem to go well, but some polishing can still be done. I will submit a PR in the near future.

An R package for diffusion maps, with additional features for large-scale and single cell data.

-   [Documentation](https://theislab.github.io/destiny/)
-   [GitHub repository](https://github.com/theislab/destiny/)
-   [ICB project](https://www.helmholtz-muenchen.de/icb/destiny)
-   [Bioconductor package](https://bioconductor.org/packages/destiny)
-   [Bioinformatics paper](https://doi.org/10.1093/bioinformatics/btv715)

## builds

-   [![R-CMD-check-bioc](https://github.com/theislab/destiny/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/theislab/destiny/actions/workflows/check-bioc.yml)
-   [![Bioconductor stable](https://bioconductor.org/shields/build/release/bioc/destiny.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/destiny/) -- Bioconductor stable
-   [![Bioconductor devel](https://bioconductor.org/shields/build/devel/bioc/destiny.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/destiny/) -- Bioconductor devel
