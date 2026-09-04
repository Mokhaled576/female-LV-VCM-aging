# ============================================================
# install_packages.R
# Install packages required by the computational workflow
# ============================================================

cran_packages <- c(
  "tidyverse",
  "janitor",
  "here",
  "ggrepel",
  "uwot",
  "harmony",
  "igraph",
  "irlba",
  "rsvd",
  "msigdbr"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "SingleCellExperiment",
  "scuttle",
  "scran",
  "scater",
  "edgeR",
  "fgsea"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

cat("Package installation check complete.\n")
