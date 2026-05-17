# Exercise Mimetics Discovery — MoTrPAC × LINCS × GSFM

Identify drug candidates that reproduce exercise-induced gene expression changes
by integrating MoTrPAC rat training transcriptomics with LINCS L1000 drug
perturbation signatures and Gene Set Foundation Model embeddings.

## Scope decisions

| Dimension | Choice | Rationale |
|-----------|--------|-----------|
| Tissues | Gastrocnemius, Liver, White Adipose (SC) | Strong exercise response; LINCS cell-line matches available |
| Timepoint | 8-week training (primary) | Chronic steady-state; closest to what a mimetic would target |
| Sex | Combined (primary) | Maximises gene count and statistical power |
| Signature size | Top 150 up + 150 down genes | Matches L1000 query tool expectations |
| LINCS access | SigCom LINCS / L2S2 API (cached) | No 50 GB download required |

## Pipeline overview

```
MoTrPAC DE tables
      │
      ▼
01 Rat→human ortholog mapping (HCOP)
      │  ├── Enrichr sanity check
      │  └── data/processed/motrpac_{tissue}_8w_{up,down}.txt
      ▼
02 LINCS L1000 query (L2S2 / SigCom API)
      │  ├── Compound-level aggregation
      │  └── results/ranked_compounds.csv
      ▼
03 GSFM embedding (maayanlab/gsfm)
      │  ├── Cosine similarity ranking
      │  └── UMAP visualization
      ▼
04 Annotation
      │  ├── IDG/Pharos target TDL
      │  ├── GTEx aging concordance
      │  └── MOA clustering
      ▼
05 Composite ranking + figures
      └── results/final_ranked_compounds.csv
```

## Quick start

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate exercise-mimetics

# 2. Get MoTrPAC data (R required) — OR use DEMO mode in notebook 01
Rscript scripts/extract_motrpac.R --outdir data/raw --timepoints 8w

# 3. Run all notebooks sequentially
make all

# OR run individually in Jupyter
jupyter lab
```

## Running individual notebooks

```bash
# Each notebook is self-contained; run in order 01 → 05
jupyter nbconvert --to notebook --execute notebooks/01_motrpac_signatures.ipynb
```

**Demo mode** (no MoTrPAC data required): set `DATA_SOURCE = 'DEMO'` in
notebook 01 cell 1 to generate synthetic DE data and test the full pipeline.

## Project structure

```
├── notebooks/
│   ├── 01_motrpac_signatures.ipynb   # DE loading, ortholog mapping, gene sets
│   ├── 02_lincs_query.ipynb          # LINCS API query, compound ranking
│   ├── 03_gsfm_embedding.ipynb       # GSFM/UMAP compound embedding
│   ├── 04_annotation.ipynb           # IDG/Pharos, aging concordance, MOA
│   └── 05_final_ranking.ipynb        # Composite scoring, final figures
├── src/
│   ├── signatures.py                 # MoTrPAC loading, ortholog mapping
│   ├── connectivity.py               # LINCS query, KS scoring, aggregation
│   └── annotation.py                 # Pharos GraphQL, aging gene sets
├── scripts/
│   └── extract_motrpac.R             # R data extraction from R package
├── data/
│   ├── README.md                     # Data acquisition instructions
│   ├── raw/                          # Downloaded source files
│   ├── processed/                    # Derived gene lists and parquets
│   └── external/                     # Ortholog maps, aging gene sets
├── results/
│   ├── ranked_compounds.csv          # Phase 2 enrichment ranking
│   ├── final_ranked_compounds.csv    # Final composite ranking
│   └── figures/                      # All output figures
├── environment.yml                   # Pinned conda environment
├── requirements.txt                  # Pip-only alternative
└── Makefile                          # make all runs end-to-end pipeline
```

## Key outputs

| File | Description |
|------|-------------|
| `results/final_ranked_compounds.csv` | Final ranked compound list with all annotation columns |
| `results/ranked_compounds.csv` | Phase 2 enrichment-only ranking (minimum deliverable) |
| `results/figures/05_heatmap_top30.png` | Top-30 compounds × tissues heatmap |
| `results/figures/03_umap_embeddings.png` | GSFM compound embedding UMAP |
| `results/figures/05_composite_score_comparison.png` | Ranking sensitivity to weights |

## Composite scoring

```
composite = 0.4 × z(enrichment) + 0.4 × z(gsfm_cosine) + 0.2 × z(anti_aging)
            + 0.5 × [has_novel_IDG_target]
```

Three weight schemes are reported; top hits should be robust across them
(Spearman ρ > 0.9 is the sanity threshold).

## Sanity gates

| Phase | Check | Expected |
|-------|-------|----------|
| 01 | Enrichr Hallmark enrichment of muscle up-set | OXIDATIVE_PHOSPHORYLATION, MYOGENESIS |
| 02 | Known mimetics in top-200 | AICAR, metformin, GW501516, or resveratrol |
| 03 | GSFM: drugs with shared MOA cluster together | PPARδ agonists co-locate in UMAP |

## Known limitations

- LINCS data from cancer-derived cell lines (HA1E, A549, MCF7, PC3) — not rat
  muscle or liver. Cell-line tissue weighting partially mitigates this.
- L1000 covers 978 landmark genes only (~12% of the transcriptome).
- Rat→human ortholog mapping loses ~10–15% of genes (HCOP, 2-database minimum).
- GSFM compound embeddings use overlap genes from API as a proxy when full gene
  signatures from GCTX are not available.

## References

- MoTrPAC Consortium (2024). *Nature* 613, 259–269.
- Subramanian et al. (2017). *Cell* 171, 1437–1452 (LINCS L1000).
- Ma'ayan Lab: SigCom LINCS, L2S2, GSFM.
- Jia et al. (2018). *PNAS* 115, E11425 (GTEx aging signatures).
- Pharos/IDG: pharos.nih.gov
