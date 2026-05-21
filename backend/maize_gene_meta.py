"""
maize_gene_meta.py — fetches human-readable metadata for maize gene models
(`Zm00001eb*`) so the UI / reports can display labels like
"SDH1_0 — succinate dehydrogenase4" next to bare accession IDs.

Primary source: Gramene Search API (data.gramene.org/search) — same data
MaizeGDB shows on its gene pages, with both gene symbol and full description.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Dict, Iterable, List, Optional

import httpx

logger = logging.getLogger(__name__)

GRAMENE_SEARCH_URL = "https://data.gramene.org/search"

# Module-level cache: gene_id -> {symbol, name, synonyms} | None (negative cache)
_META_CACHE: Dict[str, Optional[dict]] = {}

# Bounded parallelism — Gramene tolerates concurrent queries but no point
# hammering it.
_META_SEMAPHORE = asyncio.Semaphore(8)


def _normalize(gene_id: str) -> str:
    """Strip the UNIPROT: prefix we attach when AlphaFold→gene resolution fails."""
    if gene_id and gene_id.startswith("UNIPROT:"):
        return ""
    return gene_id


async def _fetch_one(client: httpx.AsyncClient, gene_id: str) -> Optional[dict]:
    """Single-gene Gramene lookup. Returns {symbol, name, synonyms} or None."""
    if not gene_id:
        return None
    if gene_id in _META_CACHE:
        return _META_CACHE[gene_id]

    async with _META_SEMAPHORE:
        try:
            r = await client.get(GRAMENE_SEARCH_URL, params={"q": gene_id, "rows": 1}, timeout=15.0)
            if r.status_code != 200:
                _META_CACHE[gene_id] = None
                return None
            docs = (((r.json() or {}).get("response") or {}).get("docs") or [])
            if not docs:
                _META_CACHE[gene_id] = None
                return None
            # Gramene wraps the right doc when the ID matches exactly
            doc = docs[0]
            if doc.get("id") != gene_id:
                # Sometimes the query matches a synonym — accept if so
                if gene_id not in (doc.get("synonyms") or []):
                    _META_CACHE[gene_id] = None
                    return None
            meta = {
                "symbol":      doc.get("name") or "",         # e.g. "SDH1_0"
                "description": doc.get("description") or "",  # e.g. "succinate dehydrogenase4"
                "synonyms":    [s for s in (doc.get("synonyms") or []) if s and s != gene_id][:6],
                "biotype":     doc.get("biotype") or "",
            }
            _META_CACHE[gene_id] = meta
            return meta
        except Exception as e:
            logger.warning(f"Gramene lookup failed for {gene_id}: {e}")
            _META_CACHE[gene_id] = None
            return None


async def fetch_maize_gene_meta_batch(gene_ids: Iterable[str]) -> Dict[str, dict]:
    """
    Look up metadata for a list of maize gene IDs in parallel. Returns
    `{gene_id: {symbol, description, synonyms, biotype}}`. Cached in-process,
    so re-runs in the same server lifetime are free.
    """
    # Dedup + normalise + skip placeholder IDs (UNIPROT:*)
    uniq = sorted({_normalize(g) for g in gene_ids if _normalize(g)})
    if not uniq:
        return {}

    async with httpx.AsyncClient(follow_redirects=True) as client:
        results = await asyncio.gather(*(_fetch_one(client, g) for g in uniq))

    out: Dict[str, dict] = {}
    for g, meta in zip(uniq, results):
        if meta is not None:
            out[g] = meta
    return out


def label_for(gene_id: str, meta: Optional[dict]) -> str:
    """
    Produce a 'sdh4 — succinate dehydrogenase4'-style label.
    Falls back to the bare gene_id when no metadata is available.
    """
    if not meta:
        return gene_id
    sym = meta.get("symbol") or ""
    desc = meta.get("description") or ""
    if sym and desc:
        return f"{sym} — {desc}"
    return sym or desc or gene_id
