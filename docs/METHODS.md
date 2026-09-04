# Methods

## Cohort construction

Donor metadata from Read et al. (2024) were used to identify female donors with left-ventricular tissue. Donors were assigned to a younger group (20–45 years) or an older group (≥55 years); donors aged 46–54 years were excluded from the primary contrast.

The resulting primary cohort comprised 24 female LV donors: 5 younger and 19 older.

## VCM extraction

The harmonized cell metadata were restricted to:

1. donors in the primary cohort,
2. `Anatomical_Site == "LV"`,
3. `Cell_Shared_Label == "Ventricular_Cardiomyocytes"`.

The corresponding columns were extracted from the full sparse Matrix Market count file with a streaming Python helper to avoid loading the approximately 15-GB source matrix into memory.

## Per-nucleus quality control

A `SingleCellExperiment` was constructed from raw counts and nucleus-level metadata. Per-nucleus total counts, detected genes and mitochondrial percentage were calculated with `scuttle::perCellQCMetrics()`.

Low-count and low-feature nuclei were identified with donor-specific `scuttle::isOutlier()` thresholds (`nmads = 3`, lower tail, log scale, `batch = Donor`). Mitochondrial percentage was examined but not included in the final exclusion rule because mitochondrial fractions were extremely low in this snRNA-seq dataset and donor-specific MAD thresholds flagged disproportionate numbers of nuclei.

The final QC object contained 70,364 VCM nuclei while retaining all 24 donors.

## Normalization and dimensionality reduction

Library-size factors were estimated with `scuttle::computeLibraryFactors()` and log-normalized expression was calculated with `scuttle::logNormCounts()`.

Biological variance was modeled with `scran::modelGeneVar()`, and the top 2,000 HVGs were used for PCA. Because the matrix was large and sparse, PCA was computed with `irlba::prcomp_irlba()` rather than dense conversion.

Unintegrated PCA/UMAP diagnostics showed substantial source-study structure.

## Harmony and clustering

Harmony was applied to PCA coordinates with `DataSource` as the correction variable. `age_group` was never used as an integration variable.

An SNN graph was built from the Harmony representation and Leiden clustering was used to define exploratory transcriptional neighborhoods. Because several clusters remained strongly enriched for individual source studies, clusters were not treated as validated biological VCM subtypes and were not used as the primary inferential unit.

## Donor-level pseudobulk

Raw counts from all QC-passing VCM nuclei belonging to the same donor were summed to create one pseudobulk library per donor. This preserves the donor as the independent biological replicate and avoids nucleus-level pseudoreplication.

Low-expression genes were filtered using `edgeR::filterByExpr()` with the study-plus-age design. TMM normalization was applied with `edgeR::calcNormFactors()`.

## Primary differential-expression model

The primary edgeR quasi-likelihood model was:

```r
~ DataSource + age_group
```

The younger group was the reference, so positive log2 fold-change represents higher expression in older females.

## Overlapping-study sensitivity analysis

To reduce study-age confounding, a sensitivity analysis was restricted to source studies containing both younger and older donors:

- Koenig
- Litvinukova
- Read

The restricted cohort contained 15 donors (4 younger and 11 older). Gene filtering, TMM normalization, dispersion estimation, QL fitting and age testing were rerun from the original unfiltered DGEList.

Robustness was assessed using effect-size correlation and directional concordance, not only FDR significance.

## Targeted pathway analysis

A predefined panel was used to inspect genes involved in CYP metabolism, EET/HETE biology, arachidonic-acid/eicosanoid pathways, estrogen/ERR signaling, PPAR/PGC-1 signaling, AHR/xenobiotic pathways and redox biology.

HVG status was not used to determine DE testability; all genes passing donor-level expression filtering could be tested.

## Genome-wide enrichment

All genes tested in the primary edgeR model were ranked with a signed statistic:

```r
sign(logFC) * sqrt(pmax(F, 0))
```

Preranked enrichment was performed with `fgsea` against MSigDB Hallmark, Reactome and GO Biological Process gene sets obtained through `msigdbr`.

Positive NES indicates enrichment toward older females; negative NES indicates enrichment toward younger females.
