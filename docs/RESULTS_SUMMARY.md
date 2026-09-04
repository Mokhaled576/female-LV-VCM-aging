# Results summary

## Primary cohort

- 24 female non-failing LV donors
- 5 younger (20–45 years)
- 19 older (≥55 years)
- 70,364 QC-passing LV VCM nuclei

## Technical structure

Nucleus-level PCA/UMAP and donor-level pseudobulk MDS showed strong source-study structure. Harmony improved mixing for visualization, but residual study-enriched neighborhoods remained.

## Primary donor-level DE

- 15,282 genes tested
- 111 genes at FDR < 0.05
- 108 higher in older
- 3 higher in younger

`CYP4F22` was a prominent targeted CYP-related signal, with an effect of approximately +3.34 log2FC in the primary model.

## Overlapping-study sensitivity

The sensitivity cohort included Koenig, Litvinukova and Read (15 donors total).

No genes reached FDR < 0.05 in the smaller cohort, but the effect estimates remained highly concordant:

- 95 primary-significant genes remained testable
- 95/95 had the same direction
- Pearson correlation of primary vs sensitivity log2FC ≈ 0.999
- median absolute log2FC difference ≈ 0.042

This supports stability of the estimated effects while demonstrating the loss of power caused by the smaller restricted cohort.

`CYP4F22` remained strongly directionally consistent (approximately +3.37 log2FC in the sensitivity analysis).

## Targeted biology

The targeted results do not support a simple global increase of every CYP, estrogen receptor or eicosanoid gene with age.

The strongest targeted gene-level result was `CYP4F22`. Additional nominal signals occurred in genes spanning PPAR/PGC-1, redox, lipid oxidation, nuclear-receptor and eicosanoid biology.

## Genome-wide pathway enrichment

Final corrected Hallmark GSEA identified four FDR-significant pathways, all enriched toward older females:

- Estrogen response late
- KRAS signaling DN
- Coagulation
- Oxidative phosphorylation

Reactome additionally supported younger enrichment of PPARα-regulated lipid metabolism.

Detailed arachidonic-acid/eicosanoid pathways showed some older-directed tendencies but did not support a broad FDR-significant activation claim.

## Interpretation

The analysis supports reproducible age-associated transcriptional remodeling in female LV VCMs despite substantial source-study heterogeneity. The strongest integrated themes involve estrogen-responsive transcription, metabolic remodeling and selected CYP-related signals.

The analysis does **not** establish:

- confirmed menopausal status,
- causal estrogen deficiency,
- altered CYP enzymatic activity,
- altered eicosanoid concentrations,
- or a causal effect of `CYP4F22`.

Those hypotheses require direct experimental validation.
