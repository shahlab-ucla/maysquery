"""
chebi_utils.py — small helpers for working with ChEBI identifiers across
adjacent databases (Rhea, UniProt) that use different protonation states.

Background: PubChem normalises compounds to a neutral (uncharged) form, so
a m/z → CMM → PubChem → ChEBI lookup typically lands on the neutral parent
(e.g. CHEBI:18012 = fumaric acid). But Rhea encodes catalytic reactions using
the physiological/ionised forms (CHEBI:29806 = fumarate(2-) via the
intermediate CHEBI:37154 = fumarate(1-)). Querying Rhea with the neutral ID
returns zero hits — silently — and the cascade falls apart downstream.

This module expands a single ChEBI ID into the full conjugate-acid/base
family using the OLS4 ontology graph, so downstream queries (Rhea SPARQL,
UniProt search) can OR all related forms together.
"""

from __future__ import annotations

import asyncio
import logging
import urllib.parse
from typing import List, Set

import httpx

logger = logging.getLogger(__name__)

_OLS_GRAPH_URL = "https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms/{enc}/graph"

# Predicate labels in the OLS graph that indicate a related protonation
# state or tautomer. Verified against the CHEBI ontology, 2026.
_CONJUGATE_LABELS = frozenset({
    "is protonated form of",
    "is deprotonated form of",
    "is tautomer of",
    "is conjugate acid of",
    "is conjugate base of",
})

# In-process cache: chebi_id -> sorted list of related ChEBI IDs (incl. self)
_RELATED_CACHE: dict[str, List[str]] = {}


def _normalize(chebi_id: str) -> str:
    """'18012' or 'chebi:18012' → 'CHEBI:18012'."""
    s = (chebi_id or "").strip()
    if not s:
        return ""
    if s.upper().startswith("CHEBI:"):
        return "CHEBI:" + s.split(":", 1)[1]
    return f"CHEBI:{s}"


def _iri(chebi_id: str) -> str:
    return f"http://purl.obolibrary.org/obo/CHEBI_{chebi_id.split(':')[-1]}"


async def _fetch_direct_neighbors(client: httpx.AsyncClient, chebi_id: str) -> Set[str]:
    """Single-hop OLS lookup — returns ChEBI IDs joined by a conjugate label."""
    enc = urllib.parse.quote(urllib.parse.quote(_iri(chebi_id), safe=""), safe="")
    url = _OLS_GRAPH_URL.format(enc=enc)
    out: Set[str] = set()
    try:
        r = await client.get(url, timeout=15.0)
        if r.status_code != 200:
            return out
        d = r.json()
        for e in d.get("edges", []):
            if e.get("label", "") not in _CONJUGATE_LABELS:
                continue
            src = e.get("source", "").rsplit("/", 1)[-1].replace("CHEBI_", "CHEBI:")
            tgt = e.get("target", "").rsplit("/", 1)[-1].replace("CHEBI_", "CHEBI:")
            # The relation is directed: source -[label]-> target. We want
            # whichever endpoint isn't `chebi_id`.
            other = tgt if src == chebi_id else src
            if other and other != chebi_id and other.startswith("CHEBI:"):
                out.add(other)
    except Exception as e:
        logger.error(f"OLS graph fetch failed for {chebi_id}: {e}")
    return out


async def expand_chebi_conjugate_forms(chebi_id: str, max_depth: int = 2) -> List[str]:
    """
    Return all ChEBI IDs (including the seed) related to `chebi_id` through
    protonation/tautomer edges within `max_depth` hops.

    `max_depth=2` is enough for diprotic acids like fumaric acid where the
    neutral form (18012) connects to the dianion (29806) only through the
    monobasic intermediate (37154). Triprotic acids may need 3.
    """
    cid = _normalize(chebi_id)
    if not cid:
        return []
    if cid in _RELATED_CACHE:
        return _RELATED_CACHE[cid]

    visited: Set[str] = {cid}
    frontier: Set[str] = {cid}

    async with httpx.AsyncClient() as client:
        for _ in range(max(0, max_depth)):
            if not frontier:
                break
            neighbor_sets = await asyncio.gather(
                *(_fetch_direct_neighbors(client, n) for n in frontier)
            )
            new_nodes: Set[str] = set()
            for ns in neighbor_sets:
                new_nodes |= ns
            new_nodes -= visited
            visited |= new_nodes
            frontier = new_nodes

    result = sorted(visited, key=lambda x: int(x.split(":")[1]))
    _RELATED_CACHE[cid] = result
    return result
