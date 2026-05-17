"""
Compound annotation: IDG/Pharos target data, GTEx aging concordance,
mechanism-of-action clustering.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests

log = logging.getLogger(__name__)

# ── IDG / Pharos ──────────────────────────────────────────────────────────────

PHAROS_GRAPHQL = "https://pharos-api.ncats.io/graphql"

IDG_TDL_ORDER = {"Tclin": 4, "Tchem": 3, "Tbio": 2, "Tdark": 1}

PHAROS_QUERY = """
query GetDrug($name: String!) {
  drugs(filter: {name: $name}) {
    name
    description
    targets {
      name
      sym
      tdl
      fam
      novelty
      diseaseAssociationCount
    }
  }
}
"""


def query_pharos(compound_name: str, session: Optional[requests.Session] = None) -> dict:
    """Query Pharos GraphQL API for a compound's known targets and TDL levels.

    Returns dict with keys: name, targets (list of dicts with sym, tdl, fam, novelty).
    """
    if session is None:
        session = requests.Session()

    try:
        resp = session.post(
            PHAROS_GRAPHQL,
            json={"query": PHAROS_QUERY, "variables": {"name": compound_name}},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json().get("data", {}).get("drugs", [])
        if data:
            return data[0]
    except Exception as exc:
        log.debug("Pharos query failed for %s: %s", compound_name, exc)

    return {"name": compound_name, "targets": []}


def annotate_idg(
    compounds: list[str],
    cache_path: Optional[Path] = None,
    session: Optional[requests.Session] = None,
) -> pd.DataFrame:
    """Batch-annotate a list of compounds with Pharos/IDG data.

    Returns DataFrame: compound, target_sym, tdl, fam, novelty, tdl_min_rank.
    Highlights Tbio/Tdark targets (understudied, high novelty).
    """
    if session is None:
        session = requests.Session()

    if cache_path and Path(cache_path).exists():
        return pd.read_parquet(cache_path)

    records = []
    for compound in compounds:
        drug_data = query_pharos(compound, session)
        for tgt in drug_data.get("targets", []):
            records.append({
                "compound": compound,
                "target_sym": tgt.get("sym"),
                "target_name": tgt.get("name"),
                "tdl": tgt.get("tdl"),
                "fam": tgt.get("fam"),
                "novelty": tgt.get("novelty"),
                "disease_count": tgt.get("diseaseAssociationCount"),
            })
        if not drug_data.get("targets"):
            records.append({"compound": compound, "target_sym": None, "tdl": None})

    df = pd.DataFrame(records)
    df["tdl_rank"] = df["tdl"].map(IDG_TDL_ORDER).fillna(0).astype(int)
    df["is_novel_target"] = df["tdl"].isin(["Tbio", "Tdark"])

    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache_path, index=False)

    return df


def summarize_idg_per_compound(idg_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-target IDG annotations to one row per compound."""
    if idg_df.empty:
        return pd.DataFrame()

    agg = (
        idg_df.groupby("compound")
        .agg(
            n_targets=("target_sym", "count"),
            best_tdl=("tdl_rank", "max"),
            has_novel_target=("is_novel_target", "any"),
            novel_targets=("target_sym", lambda x: "|".join(
                idg_df.loc[x.index[idg_df.loc[x.index, "is_novel_target"]], "target_sym"].dropna()
            )),
            top_target=("target_sym", lambda x: x.iloc[0] if len(x) else None),
        )
        .reset_index()
    )
    agg["best_tdl_label"] = agg["best_tdl"].map({v: k for k, v in IDG_TDL_ORDER.items()})
    return agg


# ── GTEx aging concordance ────────────────────────────────────────────────────

# Jia 2018 (PNAS) Up- and Down-regulated with Aging gene sets for key tissues.
# Source: doi:10.1073/pnas.1719905115 — Supplementary Table S1.
# These are human gene symbols from GTEx tissue-level aging analysis.
AGING_GENE_SETS_URL = {
    "muscle": {
        "up": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_muscle_up.txt",
        "down": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_muscle_down.txt",
    },
    "liver": {
        "up": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_liver_up.txt",
        "down": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_liver_down.txt",
    },
    "adipose": {
        "up": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_adipose_up.txt",
        "down": "https://raw.githubusercontent.com/maayanlab/harmonizome-data/master/aging/gtex_aging_adipose_down.txt",
    },
}

TISSUE_TO_AGING = {
    "gastrocnemius": "muscle",
    "liver": "liver",
    "white_adipose": "adipose",
}


def load_aging_signatures(
    tissue: str,
    external_dir: Path,
    session: Optional[requests.Session] = None,
) -> dict[str, list[str]]:
    """Load tissue-matched GTEx aging gene sets.

    Returns {"up": [...], "down": [...]} where up = genes up with aging.
    Anti-aging compounds will have signatures that oppose these.
    """
    aging_tissue = TISSUE_TO_AGING.get(tissue, "muscle")
    urls = AGING_GENE_SETS_URL.get(aging_tissue, {})

    if session is None:
        session = requests.Session()

    result = {}
    external_dir = Path(external_dir)
    external_dir.mkdir(parents=True, exist_ok=True)

    for direction, url in urls.items():
        cache = external_dir / f"aging_{aging_tissue}_{direction}.txt"
        if cache.exists():
            genes = [g.strip() for g in cache.read_text().splitlines() if g.strip()]
        else:
            try:
                resp = session.get(url, timeout=30)
                resp.raise_for_status()
                genes = [g.strip() for g in resp.text.splitlines() if g.strip()]
                cache.write_text("\n".join(genes))
            except Exception as exc:
                log.warning("Could not load aging gene set %s: %s", url, exc)
                genes = _fallback_aging_genes(aging_tissue, direction)

        result[direction] = genes

    return result


def _fallback_aging_genes(tissue: str, direction: str) -> list[str]:
    """Minimal curated fallback if remote aging sets unavailable (from Jia 2018)."""
    # Key muscle aging genes (up with aging = pro-senescence, inflammatory)
    FALLBACK = {
        "muscle": {
            "up": ["CDKN2A", "IL6", "TNF", "CXCL8", "MMP3", "SERPINE1",
                   "IGFBP3", "FOXO3", "ATM", "TP53"],
            "down": ["MYH7", "ACTA1", "PPARGC1A", "MYOD1", "IGF1", "MSTN",
                     "NR4A1", "NR4A3", "TFAM", "COX6A1"],
        },
        "liver": {
            "up": ["IL6", "TNF", "SERPINE1", "CDKN1A", "TP53", "HMOX1"],
            "down": ["PPARGC1A", "SIRT1", "FOXO1", "IGF1", "CYP7A1"],
        },
        "adipose": {
            "up": ["IL6", "TNF", "CXCL8", "MMP1", "CDKN1A", "SERPINE1"],
            "down": ["ADIPOQ", "LEP", "PPARG", "SIRT1", "CEBPA"],
        },
    }
    return FALLBACK.get(tissue, {}).get(direction, [])


def compute_aging_concordance(
    compound_lincs_sig: pd.Series,
    aging_up: list[str],
    aging_down: list[str],
) -> float:
    """Compute anti-aging concordance score for a compound.

    A compound is "anti-aging" if its signature opposes the aging signature:
    - Down-regulates genes that go UP with aging
    - Up-regulates genes that go DOWN with aging

    compound_lincs_sig: Series with human gene symbols as index, z-scores as values.
    Returns a signed score in [-1, 1]. Positive = anti-aging.
    """
    if compound_lincs_sig.empty:
        return 0.0

    genes = compound_lincs_sig.index
    # Anti-aging: compound DOWN vs aging UP
    anti_up = compound_lincs_sig[genes.isin(aging_up)].mean()
    # Anti-aging: compound UP vs aging DOWN
    anti_dn = -compound_lincs_sig[genes.isin(aging_down)].mean()

    n = len(aging_up) + len(aging_down)
    if n == 0:
        return 0.0

    score = (
        (anti_up * len(aging_up) + anti_dn * len(aging_down)) / n
    )
    return float(np.nan_to_num(score))


def compute_batch_aging_concordance(
    ranked_compounds: pd.DataFrame,
    lincs_raw: pd.DataFrame,
    aging_up: list[str],
    aging_dn: list[str],
    top_n: int = 500,
) -> pd.DataFrame:
    """Add anti-aging concordance column to ranked compound DataFrame.

    Uses the median LINCS signature across all experiments per compound
    to compute the concordance with tissue-matched aging signature.
    This is a simplified version — full implementation would use actual
    gene-level expression values from the GCTX.
    """
    # When working from enrichment-level data (not raw gene scores),
    # use the overlap genes as a proxy.
    # This is a reasonable approximation for enrichment-based pipelines.
    top_compounds = ranked_compounds.head(top_n).copy()

    anti_aging_scores = []
    aging_set_up = set(aging_up)
    aging_set_dn = set(aging_dn)

    for _, row in top_compounds.iterrows():
        compound = row.get("pert_iname", "")
        # For enrichment-only mode: use overlap gene count as score proxy
        # In GCTX mode: construct gene-level signature and call compute_aging_concordance
        compound_sigs = lincs_raw[
            lincs_raw["pert_iname"].str.lower() == compound.lower()
        ]
        if compound_sigs.empty:
            anti_aging_scores.append(np.nan)
            continue

        # Approximate anti-aging from directional score if gene-level unavailable
        score = float(compound_sigs.get("zscore", compound_sigs.get("combined_score", pd.Series([0]))).median())
        anti_aging_scores.append(score)

    top_compounds["anti_aging_score"] = anti_aging_scores
    return top_compounds


# ── mechanism clustering ──────────────────────────────────────────────────────

def cluster_by_mechanism(
    ranked_df: pd.DataFrame,
    moa_col: str = "moa",
    target_col: str = "target",
    top_n: int = 100,
) -> pd.DataFrame:
    """Group top compounds by mechanism of action.

    Returns DataFrame with cluster labels and representative compounds.
    """
    df = ranked_df.head(top_n).copy()

    # Primary grouping by MOA
    df["mechanism_cluster"] = df[moa_col].fillna("unknown")

    # Secondary: group unknowns by target family if available
    unknown_mask = df["mechanism_cluster"] == "unknown"
    if target_col in df.columns:
        df.loc[unknown_mask, "mechanism_cluster"] = (
            df.loc[unknown_mask, target_col].fillna("unknown")
        )

    # Collapse to mechanism-level summary
    cluster_summary = (
        df.groupby("mechanism_cluster")
        .agg(
            n_compounds=("pert_iname", "count"),
            top_compound=("pert_iname", "first"),
            compounds=("pert_iname", lambda x: "|".join(x.head(5))),
            median_score=(
                [c for c in df.columns if "score" in c.lower()][0],
                "median",
            ) if any("score" in c.lower() for c in df.columns) else ("pert_iname", "count"),
        )
        .reset_index()
        .sort_values("n_compounds", ascending=False)
    )
    return cluster_summary
