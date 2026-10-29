# Third-party software

The code in this repository (`R/`, `tests/`, `Dockerfile`, …) is Apache-2.0 (see `LICENSE`).
The image runs that code in R and loads the packages below **in the same R process**.
Cleanet and RcppHNSW are GPL-3.0, and flowCore links cytolib (AGPL-3.0-only), so **the
image as a whole is distributed under the GPL-3.0 / AGPL-3.0 terms of those packages**;
our own source files stay Apache-2.0 (compatible with both). Packages are installed
unmodified from CRAN (Posit Package Manager snapshot 2026-04-23), Bioconductor 3.22 and,
for Cleanet, GitHub at a fixed commit; each package's licence is inside the image under
`/usr/local/lib/R/site-library/<package>/`.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| Cleanet | 1.0.0, commit `d2532096d00f69ab911080e008cc2e27f8cbcb45` (2025-06-24) | The doublet detection (`cleanet()`) | GPL-3 | https://github.com/matei-ionita/Cleanet |
| RcppHNSW | 0.6.0 | Approximate nearest neighbours (hnswlib) used by Cleanet | GPL (>= 3); hnswlib Apache-2.0 | CRAN |
| Rcpp, RcppArmadillo | 1.1.1-1, 15.2.6-1 | Cleanet's C++ code | GPL (>= 2) | CRAN |
| flowCore | 2.22.1 | FCS read/write | Artistic-2.0 | https://bioconductor.org/packages/flowCore |
| cytolib | 2.22.0 | Linked by flowCore | AGPL-3.0-only | Bioconductor |
| RProtoBufLib, Rhdf5lib | 2.22.0, 1.32.0 | Bundled libraries of cytolib | BSD-3-Clause; Artistic-2.0 | Bioconductor |
| Biobase, BiocGenerics, S4Vectors, matrixStats | Bioc 3.22 | Bioconductor infrastructure | Artistic-2.0 | Bioconductor |
| ggplot2, gridExtra, dplyr, tidyr, readr, reshape2, tibble, magrittr | 4.0.3, 2.3, 1.2.1, … | Imported by Cleanet; the diagnostic plots | MIT; gridExtra GPL (>= 2) | CRAN |
| jsonlite | 2.0.0 | provenance.json | MIT | CRAN |
| testthat (+ deps) | 3.3.2 | Unit tests (`/app/tests`) | MIT | CRAN |
| BH | 1.90.0-1 | Boost headers (build of flowCore) | BSL-1.0 | CRAN |
| R | 4.5.3 | Interpreter | GPL-2 / GPL-3 | https://www.r-project.org |
| rocker/r-ver:4.5.3 | sha256:c3f39b36… | Base image (Ubuntu 24.04 + R) | GPL-2+ (rocker scripts); Ubuntu packages under their own free licences | https://hub.docker.com/r/rocker/r-ver |

The complete list of R packages in the image with versions and licences:
`docker run --rm --entrypoint Rscript processor-cleanet -e 'installed.packages(fields="License")[, c("Version","License")]'`.

## Data

The demo data used to test and measure the app (not in the image or the repository):
Koenig Z, Andorf S, Tilburgs T, Eswar S. *CytoBatchFlagR: Spectral flow cytometry
dataset.* Zenodo (2025). doi:10.5281/zenodo.15388817. CC BY 4.0.
