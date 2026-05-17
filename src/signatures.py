"""
MoTrPAC exercise signature extraction and rat→human ortholog mapping.
"""

import re
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

log = logging.getLogger(__name__)

# ── tissue configuration ──────────────────────────────────────────────────────

TISSUES = {
    "gastrocnemius": {
        "motrpac_id": "SKM-GN",
        "display": "Skeletal Muscle (Gastrocnemius)",
        "lincs_cell_lines": ["HA1E", "A549", "PC3"],  # rough tissue relevance
    },
    "liver": {
        "motrpac_id": "LIVER",
        "display": "Liver",
        "lincs_cell_lines": ["HepG2", "HEPG2", "HA1E"],
    },
    "white_adipose": {
        "motrpac_id": "WAT-SC",
        "display": "White Adipose (Subcutaneous)",
        "lincs_cell_lines": ["MCF7", "A549"],
    },
}

TIMEPOINTS = ["8w"]  # primary; add "1w" for acute contrast
PRIMARY_COMPARISON = "trained_vs_control"

# Known exercise-responsive genes for sanity checking (human symbols).
SANITY_GENES = {"PPARGC1A", "NR4A3", "MYC", "FOS", "PDK4", "ANGPTL4", "CPT1B", "ACSL1"}

# ── MoTrPAC data loading ──────────────────────────────────────────────────────

MOTRPAC_GCP_BASE = "https://storage.googleapis.com/motrpac-data-hub"

# Flat-file URLs for the transcriptomics differential analysis results.
# These are the processed outputs from the Nature 2024 MoTrPAC paper.
MOTRPAC_DA_URLS = {
    "gastrocnemius": (
        f"{MOTRPAC_GCP_BASE}/PASS1B-06/T58/TRNSCRPT/SKM-GN/da-results/"
        "PASS1B-06_TRNSCRPT_SKM-GN_DA.txt"
    ),
    "liver": (
        f"{MOTRPAC_GCP_BASE}/PASS1B-06/T58/TRNSCRPT/LIVER/da-results/"
        "PASS1B-06_TRNSCRPT_LIVER_DA.txt"
    ),
    "white_adipose": (
        f"{MOTRPAC_GCP_BASE}/PASS1B-06/T58/TRNSCRPT/WAT-SC/da-results/"
        "PASS1B-06_TRNSCRPT_WAT-SC_DA.txt"
    ),
}


def load_motrpac_da(tissue: str, raw_dir: Path, timepoint: str = "8w") -> pd.DataFrame:
    """Load MoTrPAC differential analysis table for one tissue.

    Tries local cache first, then downloads from GCP public bucket.
    Returns DataFrame with columns: gene_symbol, ensembl_rat, logFC, pvalue, adj_pvalue.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    cache_path = raw_dir / f"motrpac_{tissue}_{timepoint}_da.txt"

    if not cache_path.exists():
        url = MOTRPAC_DA_URLS[tissue]
        log.info("Downloading %s DA results from GCP …", tissue)
        _download(url, cache_path)

    df = pd.read_csv(cache_path, sep="\t")
    df = _standardize_motrpac_columns(df, timepoint)
    log.info("Loaded %d DE records for %s at %s", len(df), tissue, timepoint)
    return df


def _download(url: str, dest: Path, timeout: int = 120) -> None:
    """Stream-download url → dest with progress bar."""
    resp = requests.get(url, stream=True, timeout=timeout)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", 0))
    with open(dest, "wb") as fh, tqdm(
        total=total, unit="B", unit_scale=True, desc=dest.name
    ) as bar:
        for chunk in resp.iter_content(chunk_size=8192):
            fh.write(chunk)
            bar.update(len(chunk))


def _standardize_motrpac_columns(df: pd.DataFrame, timepoint: str) -> pd.DataFrame:
    """Normalise heterogeneous column names from different MoTrPAC data versions."""
    col_map = {}
    for c in df.columns:
        cl = c.lower()
        if re.search(r"gene.?symbol|symbol|gene_name", cl):
            col_map[c] = "gene_symbol"
        elif re.search(r"ensembl|ensg|ensr", cl):
            col_map[c] = "ensembl_rat"
        elif re.search(r"logfc|log2fc|log_fc|estimate", cl):
            col_map[c] = "logFC"
        elif re.search(r"adj.*p|p.*adj|fdr|qval", cl):
            col_map[c] = "adj_pvalue"
        elif re.search(r"^p.?val|pvalue|p_value", cl):
            col_map[c] = "pvalue"
        elif re.search(r"stat|t.?stat|tscore", cl):
            col_map[c] = "t_stat"
    df = df.rename(columns=col_map)

    # Filter to training-vs-control contrast at the requested timepoint if columns exist
    for tp_col in ("comparison", "contrast", "time_point", "timepoint", "week"):
        if tp_col in df.columns:
            mask = df[tp_col].astype(str).str.contains(timepoint, case=False, na=False)
            if mask.sum() > 0:
                df = df[mask].copy()
                break

    required = {"gene_symbol", "logFC", "adj_pvalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Parsed DE table is missing required columns: {missing}. "
            f"Available: {list(df.columns)}"
        )

    return df.dropna(subset=["gene_symbol", "logFC"])


def load_motrpac_from_r_export(tsv_path: Path) -> pd.DataFrame:
    """Load a TSV exported by scripts/extract_motrpac.R."""
    df = pd.read_csv(tsv_path, sep="\t")
    return _standardize_motrpac_columns(df, timepoint="8w")


# ── ortholog mapping ──────────────────────────────────────────────────────────

HCOP_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/"
    "human_rat_hcop_fifteen_column.txt.gz"
)


def load_hcop(external_dir: Path) -> pd.DataFrame:
    """Download and cache the HCOP rat→human ortholog table."""
    external_dir = Path(external_dir)
    external_dir.mkdir(parents=True, exist_ok=True)
    gz_path = external_dir / "human_rat_hcop.txt.gz"

    if not gz_path.exists():
        log.info("Downloading HCOP ortholog table …")
        _download(HCOP_URL, gz_path)

    df = pd.read_csv(
        gz_path,
        sep="\t",
        usecols=[
            "human_symbol",
            "human_entrez_gene",
            "ortholog_species_symbol",
            "ortholog_species_entrez_gene",
            "support",
        ],
        compression="gzip",
        low_memory=False,
    )
    # Rat orthologs only (taxid 10116 rows) already filtered by filename.
    df.columns = ["human_symbol", "human_entrez", "rat_symbol", "rat_entrez", "support"]
    df = df.dropna(subset=["human_symbol", "rat_symbol"])
    # Require at least 2 supporting databases for quality.
    df["n_support"] = df["support"].str.count(",") + 1
    df = df[df["n_support"] >= 2].copy()
    return df


def map_rat_to_human(
    rat_symbols: pd.Series,
    hcop: pd.DataFrame,
    one_to_one_only: bool = True,
) -> pd.DataFrame:
    """Map rat gene symbols → human orthologs via HCOP.

    Returns a DataFrame with columns rat_symbol, human_symbol.
    """
    rat_upper = rat_symbols.str.upper()
    hcop_lookup = hcop.copy()
    hcop_lookup["rat_symbol_upper"] = hcop_lookup["rat_symbol"].str.upper()

    merged = pd.DataFrame({"rat_symbol_upper": rat_upper}).merge(
        hcop_lookup[["rat_symbol_upper", "human_symbol", "rat_symbol"]],
        on="rat_symbol_upper",
        how="left",
    )

    if one_to_one_only:
        # Drop rat genes that map to multiple human genes.
        dup_mask = merged.duplicated(subset=["rat_symbol_upper"], keep=False)
        merged = merged[~dup_mask].copy()

    merged = merged.drop_duplicates(subset=["rat_symbol_upper"])
    n_mapped = merged["human_symbol"].notna().sum()
    log.info(
        "Ortholog mapping: %d/%d rat genes → human symbols (%.0f%%)",
        n_mapped,
        len(rat_symbols),
        100 * n_mapped / max(len(rat_symbols), 1),
    )
    return merged[["rat_symbol", "human_symbol"]].dropna()


# ── signature extraction ──────────────────────────────────────────────────────

def extract_signatures(
    de_df: pd.DataFrame,
    hcop: pd.DataFrame,
    top_n: int = 150,
    fdr_cutoff: float = 0.05,
    logfc_cutoff: float = 0.0,
) -> dict[str, list[str]]:
    """Return {"up": [...human symbols...], "down": [...]} gene sets.

    Filters DE table by FDR and |logFC|, maps to human orthologs,
    then takes top_n up and top_n down by |logFC|.
    """
    df = de_df.copy()
    if "adj_pvalue" in df.columns:
        df = df[df["adj_pvalue"] < fdr_cutoff]
    df = df[df["logFC"].abs() > logfc_cutoff]

    up_df = df[df["logFC"] > 0].nlargest(top_n * 3, "logFC")
    dn_df = df[df["logFC"] < 0].nsmallest(top_n * 3, "logFC")

    orth = map_rat_to_human(
        pd.concat([up_df["gene_symbol"], dn_df["gene_symbol"]]), hcop
    )
    orth_dict = dict(zip(orth["rat_symbol"].str.upper(), orth["human_symbol"]))

    def _map_set(sub_df: pd.DataFrame) -> list[str]:
        mapped = (
            sub_df["gene_symbol"]
            .str.upper()
            .map(orth_dict)
            .dropna()
            .unique()
            .tolist()
        )
        return mapped[:top_n]

    return {"up": _map_set(up_df), "down": _map_set(dn_df)}


def save_gene_set(genes: list[str], path: Path) -> None:
    """Write one gene symbol per line."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(genes) + "\n")


def load_gene_set(path: Path) -> list[str]:
    return [g.strip() for g in Path(path).read_text().splitlines() if g.strip()]


def sanity_check_signature(up_genes: list[str], tissue: str) -> dict:
    """Check overlap with known exercise-responsive genes."""
    found = SANITY_GENES & set(up_genes)
    result = {
        "tissue": tissue,
        "up_n": len(up_genes),
        "sanity_genes_found": sorted(found),
        "sanity_score": len(found) / len(SANITY_GENES),
    }
    if result["sanity_score"] < 0.2:
        log.warning(
            "Sanity check FAILED for %s: only %d/%d expected genes found: %s",
            tissue, len(found), len(SANITY_GENES), found,
        )
    else:
        log.info(
            "Sanity check PASSED for %s: %d/%d expected genes found",
            tissue, len(found), len(SANITY_GENES),
        )
    return result
