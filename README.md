
<!-- README.md is generated from README.Rmd. Please edit that file -->

# annotateR

<!-- badges: start -->
<!-- badges: end -->

The goal of annotateR is to enable annotation of genome wide association
study (GWAS) results (or other tabular data including chromosome,
position, and allele data) with RSIDs from a reference dataset

## Installation

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mglev1n/annotateR")
```

## Reference Datasets

| dataset              |     nrow | ncol | format  | columns                          |
|:---------------------|---------:|-----:|:--------|:---------------------------------|
| 1000 Genomes Phase 3 | 84805772 |    6 | parquet | chr, pos, rsid, ref, alt, marker |
