# Data directory

## Directory layout

```
data/
├── raw/          # Downloaded source files (not committed to git)
├── processed/    # Derived signatures (gene lists, parquet)
└── external/     # Ortholog maps, IDG, aging gene sets
```

## Files required before running the pipeline

### 1. MoTrPAC transcriptomics DE results

**Option A — R package (recommended, easiest)**

```bash
Rscript scripts/extract_motrpac.R --outdir data/raw --timepoints 8w
```

This writes `data/raw/motrpac_{tissue}_{timepoint}_r_export.tsv` for
gastrocnemius, liver, and white_adipose at 8-week timepoint.

Requires the `MotrpacRatTraining6moData` R package:
```r
remotes::install_github("MoTrPAC/MotrpacRatTraining6moData")
```

**Option B — MoTrPAC Data Hub (manual download)**

1. Create a Google account and request access at https://motrpac-data.org/
2. Download transcriptomics DA files from the GCP bucket:
   ```
   gs://motrpac-data-hub/PASS1B-06/T58/TRNSCRPT/{TISSUE}/da-results/
   ```
   Replace `{TISSUE}` with: `SKM-GN`, `LIVER`, `WAT-SC`

3. Place files in `data/raw/` with names matching the URL basenames.
   Notebook 01 auto-downloads from the public GCP HTTP endpoint if the
   bucket is publicly accessible.

**Option C — Demo / synthetic mode**

Notebook 01 can generate synthetic DE data for pipeline testing.
Set `USE_DEMO_DATA = True` in the first cell.

---

### 2. HCOP rat→human ortholog table (auto-downloaded)

Notebook 01 downloads automatically from the EBI FTP:
```
https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_rat_hcop_fifteen_column.txt.gz
```
Saved to: `data/external/human_rat_hcop.txt.gz`

---

### 3. LINCS L1000 signatures (API, no local download required)

Notebooks 02–03 query the SigCom LINCS and L2S2 APIs. Results are
cached as parquet files in `data/processed/` after the first run.
No local LINCS download needed unless APIs are unavailable.

**Fallback — local GCTX (~50 GB):**

If APIs are down, download from GEO:
```
GSE92742 — Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
GSE92742 — Broad_LINCS_gene_info.txt.gz
GSE92742 — Broad_LINCS_sig_info.txt.gz
```
Place in `data/raw/` and set `USE_LOCAL_GCTX = True` in notebook 02.

---

### 4. GTEx aging gene sets (auto-downloaded)

Notebook 04 downloads aging gene sets from Harmonizome via GitHub.
Fallback curated gene sets are included in `src/annotation.py`.

---

### 5. GSFM model (auto-downloaded by HuggingFace)

Notebook 03 pulls `maayanlab/gsfm` from HuggingFace Hub automatically.
Requires internet access and ~2 GB disk space for the model weights.

---

## File sizes (approximate)

| File | Size |
|------|------|
| MoTrPAC DA TSV (per tissue) | ~5 MB |
| HCOP ortholog table (gz) | ~15 MB |
| LINCS Level 5 GCTX (optional) | ~50 GB |
| GSFM model weights | ~2 GB |
| Processed signatures (parquet) | ~5 MB total |

## Notes

- `data/raw/` and `data/external/` large files are excluded from git via `.gitignore`.
- Processed parquet outputs in `data/processed/` are committed (small files).
- All download steps cache to disk — re-running notebooks will use the cache.
