"""
LINCS L1000 connectivity scoring via SigCom LINCS / L2S2 API.
Includes local fallback using cmapPy if APIs are unavailable.
"""

import json
import logging
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

log = logging.getLogger(__name__)

# ── API configuration ─────────────────────────────────────────────────────────

SIGCOM_BASE = "https://maayanlab.cloud/sigcom-lincs/api/v1"
L2S2_BASE = "https://maayanlab.cloud/l2s2"

# Known exercise mimetics for Phase 2 sanity check.
KNOWN_MIMETICS = {
    "AICAR",         # AMPK activator
    "metformin",
    "resveratrol",
    "GW501516",      # PPARδ agonist
    "bezafibrate",
    "fenofibrate",
    "pioglitazone",
    "rosiglitazone",
    "caffeine",
    "epicatechin",
}


# ── SigCom LINCS query ────────────────────────────────────────────────────────

def query_sigcom_lincs(
    up_genes: list[str],
    down_genes: list[str],
    entity_type: str = "signatures",
    limit: int = 5000,
    session: Optional[requests.Session] = None,
) -> pd.DataFrame:
    """Query SigCom LINCS overlap enrichment endpoint.

    Returns a DataFrame of enriched signatures sorted by p-value.
    API docs: https://maayanlab.cloud/sigcom-lincs/api/swagger/
    """
    if session is None:
        session = requests.Session()

    payload = {
        "up_genes": up_genes,
        "down_genes": down_genes,
        "entity_type": entity_type,
    }

    url = f"{SIGCOM_BASE}/enrich/overlap"
    log.info("POSTing to SigCom LINCS (%d up, %d down genes) …", len(up_genes), len(down_genes))

    resp = session.post(url, json=payload, timeout=120)
    resp.raise_for_status()
    results_url = resp.json().get("results")

    # Poll for results if async
    if results_url:
        return _poll_results(session, results_url, limit)

    # Synchronous response
    data = resp.json()
    return _parse_sigcom_response(data, limit)


def _poll_results(
    session: requests.Session,
    results_url: str,
    limit: int,
    max_wait: int = 300,
) -> pd.DataFrame:
    """Poll SigCom async results endpoint until ready."""
    waited = 0
    interval = 5
    while waited < max_wait:
        r = session.get(results_url, params={"limit": limit}, timeout=60)
        if r.status_code == 200:
            return _parse_sigcom_response(r.json(), limit)
        time.sleep(interval)
        waited += interval
    raise TimeoutError(f"SigCom results not ready after {max_wait}s")


def _parse_sigcom_response(data: dict | list, limit: int) -> pd.DataFrame:
    """Normalize SigCom JSON response into a flat DataFrame."""
    if isinstance(data, list):
        records = data
    elif "results" in data:
        records = data["results"]
    else:
        records = [data]

    rows = []
    for r in records[:limit]:
        meta = r.get("metadata", r)
        rows.append({
            "sig_id": r.get("id") or meta.get("sig_id"),
            "pert_iname": meta.get("pert_iname") or meta.get("pert_desc"),
            "cell_id": meta.get("cell_id") or meta.get("cell_line"),
            "pert_type": meta.get("pert_type"),
            "pert_dose": meta.get("pert_dose"),
            "pert_time": meta.get("pert_time"),
            "moa": meta.get("moa"),
            "target": meta.get("target"),
            "pvalue": r.get("pvalue") or r.get("p-value"),
            "qvalue": r.get("qvalue") or r.get("q-value"),
            "zscore": r.get("zscore") or r.get("z-score"),
            "combined_score": r.get("combined_score"),
            "overlap_up": r.get("overlap_up") or r.get("n_up"),
            "overlap_down": r.get("overlap_down") or r.get("n_down"),
        })
    return pd.DataFrame(rows)


# ── L2S2 query ────────────────────────────────────────────────────────────────

def query_l2s2(
    up_genes: list[str],
    down_genes: list[str],
    limit: int = 5000,
    session: Optional[requests.Session] = None,
) -> pd.DataFrame:
    """Query L2S2 (LINCS L2 Signatures v2) API.

    L2S2 searches 1.678M LINCS signatures and returns enrichment statistics.
    Returns sorted by absolute zscore descending.
    """
    if session is None:
        session = requests.Session()

    payload = {"up_genes": up_genes, "down_genes": down_genes}
    url = f"{L2S2_BASE}/enrich/overlap"
    log.info("POSTing to L2S2 …")

    try:
        resp = session.post(url, json=payload, timeout=120)
        resp.raise_for_status()
        data = resp.json()
        df = _parse_sigcom_response(data, limit)
        log.info("L2S2 returned %d signatures", len(df))
        return df
    except Exception as exc:
        log.warning("L2S2 query failed (%s); trying SigCom LINCS fallback", exc)
        return query_sigcom_lincs(up_genes, down_genes, limit=limit, session=session)


# ── result caching ────────────────────────────────────────────────────────────

def cache_query(
    up_genes: list[str],
    down_genes: list[str],
    cache_path: Path,
    force: bool = False,
    api: str = "l2s2",
) -> pd.DataFrame:
    """Run LINCS query once, cache JSON to disk, reload on subsequent calls."""
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists() and not force:
        log.info("Loading cached LINCS results from %s", cache_path)
        return pd.read_parquet(cache_path)

    if api == "l2s2":
        df = query_l2s2(up_genes, down_genes)
    else:
        df = query_sigcom_lincs(up_genes, down_genes)

    df.to_parquet(cache_path, index=False)
    log.info("Cached %d records → %s", len(df), cache_path)
    return df


# ── aggregation ───────────────────────────────────────────────────────────────

def aggregate_by_compound(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-signature scores to one score per compound.

    Strategy: for each compound, take the median z-score across all
    experimental conditions (dose × cell-line × time), then also record
    the max-absolute z-score for sensitivity.
    """
    if df.empty:
        return df

    # Fill missing pert_iname from sig_id
    df = df.copy()
    df["pert_iname"] = df["pert_iname"].fillna(df["sig_id"])

    # Numeric cleanup
    for col in ("zscore", "combined_score", "pvalue", "overlap_up", "overlap_down"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    score_col = "zscore" if "zscore" in df.columns else "combined_score"

    agg = (
        df.groupby("pert_iname")
        .agg(
            n_signatures=(score_col, "count"),
            median_score=(score_col, "median"),
            max_abs_score=(score_col, lambda x: x.abs().max()),
            mean_score=(score_col, "mean"),
            moa=("moa", lambda x: x.dropna().mode().iloc[0] if x.dropna().size else np.nan),
            target=("target", lambda x: x.dropna().mode().iloc[0] if x.dropna().size else np.nan),
            cell_lines=("cell_id", lambda x: "|".join(x.dropna().unique()[:5])),
        )
        .reset_index()
    )

    agg["rank_median"] = agg["median_score"].rank(ascending=False, method="min")
    agg = agg.sort_values("median_score", ascending=False).reset_index(drop=True)
    return agg


def cell_line_weighted_score(
    df: pd.DataFrame,
    tissue: str,
    tissue_weights: Optional[dict[str, float]] = None,
) -> pd.DataFrame:
    """Re-weight signatures by cell-line tissue relevance before aggregation.

    Signatures from cell lines more relevant to the queried tissue
    receive a higher weight in the compound-level aggregation.
    """
    from .signatures import TISSUES
    relevant_lines = set(TISSUES.get(tissue, {}).get("lincs_cell_lines", []))

    if tissue_weights is None:
        tissue_weights = {line: 2.0 for line in relevant_lines}

    df = df.copy()
    df["cell_weight"] = df["cell_id"].map(tissue_weights).fillna(1.0)

    score_col = "zscore" if "zscore" in df.columns else "combined_score"
    df["weighted_score"] = df[score_col] * df["cell_weight"]

    agg = (
        df.groupby("pert_iname")
        .agg(
            n_signatures=(score_col, "count"),
            weighted_median_score=("weighted_score", "median"),
            max_abs_score=(score_col, lambda x: x.abs().max()),
            moa=("moa", lambda x: x.dropna().mode().iloc[0] if x.dropna().size else np.nan),
            target=("target", lambda x: x.dropna().mode().iloc[0] if x.dropna().size else np.nan),
            cell_lines=("cell_id", lambda x: "|".join(x.dropna().unique()[:5])),
        )
        .reset_index()
    )

    agg = agg.sort_values("weighted_median_score", ascending=False).reset_index(drop=True)
    return agg


def sanity_check_ranking(ranked_df: pd.DataFrame, tissue: str) -> dict:
    """Check whether known mimetics appear in the top-200 hits."""
    top200 = set(ranked_df.head(200)["pert_iname"].str.lower())
    found = {m for m in KNOWN_MIMETICS if m.lower() in top200}
    result = {
        "tissue": tissue,
        "n_ranked": len(ranked_df),
        "known_mimetics_in_top200": sorted(found),
        "sanity_score": len(found) / len(KNOWN_MIMETICS),
    }
    if result["sanity_score"] < 0.2:
        log.warning(
            "Phase 2 sanity FAILED for %s: only %d/%d known mimetics in top 200",
            tissue, len(found), len(KNOWN_MIMETICS),
        )
    else:
        log.info(
            "Phase 2 sanity PASSED for %s: %d/%d known mimetics found",
            tissue, len(found), len(KNOWN_MIMETICS),
        )
    return result


# ── local GCTX fallback ───────────────────────────────────────────────────────

def compute_connectivity_local(
    query_up: list[str],
    query_down: list[str],
    gctx_path: Path,
    gene_info_path: Path,
    sig_info_path: Path,
    chunk_size: int = 5000,
) -> pd.DataFrame:
    """Compute weighted KS connectivity scores against LINCS Level 5 GCTX.

    Only use this if the APIs are unavailable. The full GCTX is ~50GB;
    download from GEO GSE92742 or clue.io.

    Parameters
    ----------
    gctx_path:      Path to GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
    gene_info_path: Path to GSE92742_Broad_LINCS_gene_info.txt.gz
    sig_info_path:  Path to GSE92742_Broad_LINCS_sig_info.txt.gz
    """
    try:
        from cmapPy.pandasGEXpress import parse_gctoo
    except ImportError:
        raise ImportError("cmapPy required for local GCTX mode: pip install cmapPy")

    gene_info = pd.read_csv(gene_info_path, sep="\t")
    landmark = gene_info[gene_info["pr_is_lm"] == 1]

    query_rids = list(landmark[landmark["pr_gene_symbol"].isin(query_up + query_down)]["pr_gene_id"].astype(str))

    log.info("Loading %d landmark genes from GCTX …", len(query_rids))
    gct = parse_gctoo.parse(str(gctx_path), rid=query_rids)
    mat = gct.data_df  # genes × signatures

    up_mask = mat.index.isin(
        landmark[landmark["pr_gene_symbol"].isin(query_up)]["pr_gene_id"].astype(str)
    )
    dn_mask = mat.index.isin(
        landmark[landmark["pr_gene_symbol"].isin(query_down)]["pr_gene_id"].astype(str)
    )

    scores = _ks_score_vectorized(mat, up_mask, dn_mask)
    sig_info = pd.read_csv(sig_info_path, sep="\t")
    result = sig_info.copy()
    result["ks_score"] = scores.values
    return result.sort_values("ks_score", ascending=False)


def _ks_score_vectorized(mat: pd.DataFrame, up_mask, dn_mask) -> pd.Series:
    """Vectorised signed KS connectivity score (Lamb 2006)."""
    from scipy.stats import rankdata

    n_genes = mat.shape[0]
    scores = []
    for col in tqdm(mat.columns, desc="KS scoring"):
        ranked = rankdata(mat[col])
        es_up = _es_one_tail(ranked, up_mask, n_genes)
        es_dn = _es_one_tail(ranked, dn_mask, n_genes)
        if np.sign(es_up) == np.sign(es_dn):
            scores.append(0.0)
        else:
            scores.append((es_up - es_dn) / 2)
    return pd.Series(scores, index=mat.columns)


def _es_one_tail(ranked: np.ndarray, mask: np.ndarray, n: int) -> float:
    """One-tail enrichment score for KS connectivity."""
    if mask.sum() == 0:
        return 0.0
    hit_sum = np.cumsum(mask / mask.sum())
    miss_sum = np.cumsum(~mask / (~mask).sum())
    return float(np.max(np.abs(hit_sum - miss_sum)) * np.sign(np.max(hit_sum - miss_sum)))
