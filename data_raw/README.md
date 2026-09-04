# Source data (not tracked by Git)

Place the source files required by the analysis in this directory.

## Read et al. (2024) supplementary donor metadata

Required:

```text
42003_2024_6582_MOESM18_ESM.csv
```

## Harmonized multi-study snRNA-seq files

Required by the scripts:

```text
multi_study_snRNA_Seq_cell_data.csv
multi_study_snRNA_Seq_gene_data.csv
multi_study_snRNA_Seq_counts_mm_file.txt
```

The source publication also provides a processed CDS/RDS object, but this workflow uses the deposited cell metadata, gene metadata and Matrix Market count file.

Source study:

Read DF et al. Single-cell analysis of chromatin and expression reveals age- and sex-associated alterations in the human heart. Communications Biology 7, 1052 (2024). DOI: 10.1038/s42003-024-06582-y

Large source data are intentionally ignored by `.gitignore`.
