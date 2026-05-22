"""
corncyc_lookup.py — high-level query against the loaded CornCyc PGDB.

Given a ChEBI ID (the Phase-1 output of our pipeline), produce a structured
annotation containing:

  * the matching CornCyc compound frame(s)
  * the reactions the compound participates in
  * the maize-specific pathways those reactions belong to
  * the maize gene models annotated as catalysts (Zm00001eb* v5 IDs)

The result is consumed in two places downstream:

  1. Phase 4 — each CornCyc-annotated maize gene becomes an `OrthologMapping`
     entry with `sources=["CornCyc"]`, joining the consensus reducer
     alongside Ensembl Compara, Foldseek-structural, etc.
  2. Reports + UI — the pathway list is shown as a dedicated "Maize
     pathway context (CornCyc)" section, with per-pathway gene lists.

If CornCyc isn't installed, `corncyc_annotation_for_chebi()` returns
`None` and the rest of the pipeline keeps working — CornCyc is opt-in.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

from corncyc_loader import get_index
from models import OrthologMapping

logger = logging.getLogger(__name__)

CORNCYC_SOURCE = "CornCyc"


def corncyc_annotation_for_chebi(chebi_id: str) -> Optional[dict]:
    """
    Build a per-compound CornCyc annotation. Returns `None` if CornCyc
    isn't installed OR if the ChEBI ID isn't in CornCyc at all. The
    pathway list is ordered by the number of reactions in the pathway
    that touch this compound (most-involved pathways first).
    """
    idx = get_index()
    if idx is None:
        return None
    compound_frames = idx.compounds_for_chebi(chebi_id)
    if not compound_frames:
        return None

    # Aggregate across (rare) multi-compound matches
    compounds_info: List[dict] = []
    reactions_info: Dict[str, dict] = {}
    pathway_to_rxns: Dict[str, List[str]] = {}
    pathway_to_genes: Dict[str, set] = {}
    gene_to_evidence: Dict[str, dict] = {}

    for cf in compound_frames:
        cmeta = idx.compounds.get(cf, {})
        compounds_info.append({
            "id": cf,
            "name": cmeta.get("name", cf),
            "chebi": cmeta.get("chebi", ""),
            "pubchem": cmeta.get("pubchem", ""),
            "kegg": cmeta.get("kegg", ""),
            "hmdb": cmeta.get("hmdb", ""),
            "synonyms": cmeta.get("synonyms", []),
        })
        for rid in idx.reactions_for_compound(cf):
            rxn = idx.reactions.get(rid, {})
            reactions_info[rid] = {
                "id": rid,
                "common_name": rxn.get("common_name", ""),
                "ec_numbers": rxn.get("ec", []),
                "pathways": list(rxn.get("in_pathway", [])),
            }
            for pid in idx.pathways_for_reaction(rid):
                pathway_to_rxns.setdefault(pid, [])
                if rid not in pathway_to_rxns[pid]:
                    pathway_to_rxns[pid].append(rid)
            for gene in idx.genes_for_reaction(rid):
                for pid in idx.pathways_for_reaction(rid):
                    pathway_to_genes.setdefault(pid, set()).add(gene)
                ev = gene_to_evidence.setdefault(gene, {
                    "gene": gene,
                    "reactions": set(),
                    "pathways": set(),
                    "ec_numbers": set(),
                })
                ev["reactions"].add(rid)
                for pid in idx.pathways_for_reaction(rid):
                    ev["pathways"].add(pid)
                for ec in rxn.get("ec", []):
                    if ec:
                        ev["ec_numbers"].add(ec)

    # Materialise + sort pathways by involvement strength
    pathways_out: List[dict] = []
    for pid, rxns in pathway_to_rxns.items():
        pmeta = idx.pathways.get(pid, {})
        pathways_out.append({
            "id": pid,
            "common_name": pmeta.get("common_name", pid),
            "synonyms": pmeta.get("synonyms", []),
            "types": pmeta.get("types", []),
            "reactions_touching_compound": rxns,
            "n_reactions_in_pathway": len(pmeta.get("reactions", [])),
            "maize_genes": sorted(pathway_to_genes.get(pid, set())),
        })
    pathways_out.sort(key=lambda p: (-len(p["reactions_touching_compound"]), p["common_name"]))

    # Materialise + sort gene evidence by # reactions touched
    genes_out: List[dict] = []
    for g, ev in gene_to_evidence.items():
        genes_out.append({
            "gene": g,
            "reactions": sorted(ev["reactions"]),
            "pathways": sorted(ev["pathways"]),
            "ec_numbers": sorted(ev["ec_numbers"]),
            "n_reactions": len(ev["reactions"]),
        })
    genes_out.sort(key=lambda g: (-g["n_reactions"], g["gene"]))

    return {
        "version": idx.version,
        "compounds": compounds_info,
        "reactions": list(reactions_info.values()),
        "pathways": pathways_out,
        "maize_genes": genes_out,
        "n_compounds": len(compounds_info),
        "n_reactions": len(reactions_info),
        "n_pathways": len(pathways_out),
        "n_maize_genes": len(genes_out),
    }


def corncyc_orthologs_for_chebi(chebi_id: str) -> List[OrthologMapping]:
    """
    Turn the CornCyc annotation into `OrthologMapping` entries so they
    flow through Phase 4's consensus reducer. Each unique maize gene
    becomes one entry with `sources=["CornCyc"]`.

    `similarity_score` is fixed at 100.0 — CornCyc is curated annotation,
    not a homology search, so there's no continuous score. The
    `plaza_orthogroup` field encodes pathway membership for visibility.
    """
    ann = corncyc_annotation_for_chebi(chebi_id)
    if not ann:
        return []

    # Pathway-name lookup for label encoding
    pathway_names = {p["id"]: p["common_name"] for p in ann["pathways"]}

    out: List[OrthologMapping] = []
    for g in ann["maize_genes"]:
        # Show up to 2 pathway names + count in the orthogroup label
        names = [pathway_names.get(pid, pid) for pid in g["pathways"]]
        names = [n for n in names if n]
        label = "CornCyc " + (
            f"{names[0]}"
            + (f" (+{len(names)-1} more pathway{'s' if len(names) > 2 else ''})"
               if len(names) > 1 else "")
            if names else f"({g['n_reactions']} rxn{'s' if g['n_reactions'] != 1 else ''})"
        )
        out.append(OrthologMapping(
            query_uniprot_id="CornCyc-curated",
            maize_gene_model=g["gene"],
            plaza_orthogroup=label,
            similarity_score=100.0,
            sources=[CORNCYC_SOURCE],
            consensus_score=1,
        ))
    return out


def corncyc_pathway_summary_for_chebi(chebi_id: str) -> Optional[dict]:
    """
    Compact summary for Phase 2's reaction-network log line:
        {n_pathways, top_pathway_name, n_maize_genes}
    """
    ann = corncyc_annotation_for_chebi(chebi_id)
    if not ann:
        return None
    top = ann["pathways"][0] if ann["pathways"] else None
    return {
        "n_pathways": ann["n_pathways"],
        "n_maize_genes": ann["n_maize_genes"],
        "top_pathway_name": (top or {}).get("common_name", ""),
        "top_pathway_id": (top or {}).get("id", ""),
    }
