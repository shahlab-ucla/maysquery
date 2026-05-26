"""
phytozome_lookup.py — batched BioMart query against JGI Phytozome (v14) for
maize gene annotation.

Phytozome's anonymous BioMart at https://phytozome-next.jgi.doe.gov/biomart
exposes per-gene KEGG KO descriptions, Panther family classifications, Pfam
domains, and gene descriptions — all independent of what we get from
Gramene / InterPro / CornCyc. The Panther family in particular is a
useful corroborating signal: when two of our discovery lanes both surface
the same Panther family, that's another vote of confidence.

One BioMart request per pipeline run (all genes in a single batched query),
~1-3 s on a reasonable connection. Cached in-process so re-runs of the same
gene set across queries are free.

Gracefully degrades — if BioMart is down/slow/changed, we log and return
an empty dict; the pipeline keeps working without Phytozome enrichment.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Dict, Iterable, List, Optional

import httpx

logger = logging.getLogger(__name__)

BIOMART_URL = "https://phytozome-next.jgi.doe.gov/biomart/martservice"

# Maximum genes per BioMart request — JGI tolerates ~500 at a time; bigger
# batches risk URL-length / server-timeout issues.
MAX_BATCH = 250

# Module-level cache: gene_id -> {description, panther_id, panther_desc, ...}
_META_CACHE: Dict[str, dict] = {}
# Negative cache for genes BioMart couldn't find — avoid re-querying.
_NEG_CACHE: set = set()


def _biomart_xml(gene_ids: List[str]) -> str:
    """Build the BioMart XML query string for a batch of maize gene IDs."""
    csv = ",".join(gene_ids)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="zome_mart" formatter="TSV" header="1" datasetConfigVersion="0.6">
  <Dataset name="phytozome" interface="default">
    <Filter name="gene_name_filter" value="{csv}"/>
    <Attribute name="organism_name"/>
    <Attribute name="gene_name1"/>
    <Attribute name="gene_description"/>
    <Attribute name="panther_id"/>
    <Attribute name="panther_desc"/>
  </Dataset>
</Query>"""


async def _fetch_batch(client: httpx.AsyncClient, gene_ids: List[str]) -> Dict[str, dict]:
    """One BioMart POST returning {gene_id: {organism, description, panther_id, panther_desc}}."""
    if not gene_ids:
        return {}
    xml = _biomart_xml(gene_ids)
    out: Dict[str, dict] = {}
    try:
        r = await client.post(BIOMART_URL, data={"query": xml}, timeout=90.0)
        if r.status_code != 200:
            logger.warning(f"Phytozome BioMart returned HTTP {r.status_code}")
            return {}
        # TSV with header row. Each row corresponds to a (gene, panther_id) pair —
        # the same gene can appear multiple times if it has several Panther IDs.
        # We collapse by gene_id, keeping the first row's description + the
        # union of distinct Panther IDs.
        lines = r.text.splitlines()
        if len(lines) < 2:
            return {}
        for line in lines[1:]:
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            organism, gene, desc, panther_id, panther_desc = (cols + ["", "", "", "", ""])[:5]
            if not gene:
                continue
            entry = out.setdefault(gene, {
                "organism":     organism,
                "description":  desc.strip(),
                "panther_ids":  [],
                "panther_descs": [],
            })
            if panther_id and panther_id not in entry["panther_ids"]:
                entry["panther_ids"].append(panther_id)
            if panther_desc and panther_desc not in entry["panther_descs"]:
                entry["panther_descs"].append(panther_desc)
    except Exception as e:
        logger.error(f"Phytozome BioMart batch failed ({type(e).__name__}): {e!r}")
    return out


async def fetch_phytozome_meta_batch(gene_ids: Iterable[str]) -> Dict[str, dict]:
    """
    Look up Phytozome annotation for a list of maize gene IDs in batches.

    Returns `{gene_id: {organism, description, panther_ids, panther_descs}}`.
    Only v5 NAM IDs (`Zm00001eb*`) work in Phytozome v14 — v3/v4 IDs are
    auto-skipped (they're surfaced separately as Gramene synonyms).
    """
    # Phytozome v14 only knows v5 NAM gene names. Filter early.
    uniq = sorted({g for g in gene_ids if g and g.startswith("Zm00001eb")})
    # Pull anything already in cache (positive or negative)
    to_query = [g for g in uniq if g not in _META_CACHE and g not in _NEG_CACHE]
    if not to_query:
        return {g: _META_CACHE[g] for g in uniq if g in _META_CACHE}

    async with httpx.AsyncClient(follow_redirects=True) as client:
        batches = [to_query[i:i + MAX_BATCH] for i in range(0, len(to_query), MAX_BATCH)]
        results_list = await asyncio.gather(*(_fetch_batch(client, b) for b in batches))

    fetched: Dict[str, dict] = {}
    for r in results_list:
        fetched.update(r)
    # Populate cache (positive + negative)
    for g in to_query:
        if g in fetched:
            _META_CACHE[g] = fetched[g]
        else:
            _NEG_CACHE.add(g)
    # Merge cached + fresh and return
    return {g: _META_CACHE[g] for g in uniq if g in _META_CACHE}


def url_phytozome_gene(gene_id: str) -> str:
    """Phytozome gene-page deep link — works for v5 NAM IDs."""
    return f"https://phytozome-next.jgi.doe.gov/genePage/{gene_id}"
