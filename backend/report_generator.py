import os
import csv
import io
import uuid
from datetime import datetime
from jinja2 import Template
from typing import List, Dict

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports")
os.makedirs(REPORTS_DIR, exist_ok=True)

# ---------- URL builders (kept here so HTML + CSV stay consistent) ----------

def url_maize_gene(gene_id: str) -> str:
    return f"https://www.maizegdb.org/gene_center/gene/{gene_id}"

def url_ensembl_plants_gene(species: str, gene_id: str) -> str:
    sp = (species or "").strip()
    if sp:
        sp_url = sp[:1].upper() + sp[1:]
        return f"https://plants.ensembl.org/{sp_url}/Gene/Summary?g={gene_id}"
    return f"https://plants.ensembl.org/Multi/Search/Results?q={gene_id}"

def url_uniprot(accession: str) -> str:
    return f"https://www.uniprot.org/uniprotkb/{accession}/entry"

def url_chebi(chebi_id: str) -> str:
    return f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"

def url_pubchem(cid) -> str:
    return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

def url_rhea(rhea_id: str) -> str:
    return f"https://www.rhea-db.org/rhea/{str(rhea_id).replace('RHEA:', '')}"

def url_pfam(pfam_id: str) -> str:
    return f"https://www.ebi.ac.uk/interpro/entry/pfam/{pfam_id}/"

def url_kegg_pathway(pid: str) -> str:
    cleaned = (pid or "").split(":")[-1]
    return f"https://www.kegg.jp/entry/{cleaned}"

def url_quickgo(go_id: str) -> str:
    return f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"


def gene_label_html(gene_id: str, meta_map: Dict, with_link: bool = True) -> str:
    """
    Render 'Zm00001eb117970 · sdh4 — succinate dehydrogenase4' for the HTML
    report. Falls back to bare ID + MaizeGDB link when no metadata exists.
    """
    if not gene_id:
        return "—"
    is_placeholder = gene_id.startswith("UNIPROT:")
    href = (
        f"https://www.maizegdb.org/search?query={gene_id.replace('UNIPROT:', '')}"
        if is_placeholder else url_maize_gene(gene_id)
    )
    meta = (meta_map or {}).get(gene_id) or {}
    sym = meta.get("symbol", "")
    desc = meta.get("description", "")
    if sym and desc:
        meta_html = f' <span class="gene-meta">{sym} — {desc}</span>'
    elif sym or desc:
        meta_html = f' <span class="gene-meta">{sym or desc}</span>'
    else:
        meta_html = ""
    if with_link:
        return f'<a href="{href}" target="_blank" class="mono"><b>{gene_id}</b></a>{meta_html}'
    return f'<span class="mono"><b>{gene_id}</b></span>{meta_html}'

# ---------- Discovery-lane classification ----------

STRUCTURAL_SOURCE = "Foldseek-structural"
CURATED_SOURCE    = "CornCyc"
SEQUENCE_SOURCES  = {"Ensembl", "PLAZA", "BioMart", "EnsemblPanHomology", "Local_HMMER"}

# All discovery sources we know how to report per-source columns for. Each pairs
# with a short, human-readable description shown in the CSV's discovery_methods_summary.
KNOWN_SOURCES = [
    ("Ensembl",             "Ensembl Compara plants (sequence homology)"),
    ("EnsemblPanHomology",  "Ensembl pan-homology Compara (cross-kingdom sequence)"),
    ("PLAZA",               "PLAZA Monocots v5 (sequence homology)"),
    ("Local_HMMER",         "Local HMMER (phmmer) vs Zm-NAM-5.0 proteome"),
    ("Foldseek-structural", "Foldseek 3D structural alignment vs maize AlphaFold proteome"),
    ("CornCyc",             "CornCyc (Plant Metabolic Network) curated maize annotation"),
]


def classify_ortholog(sources: List[str]) -> str:
    """
    Classify a maize ortholog by which discovery lanes found it:

      - 'consensus'      : evidence from ≥2 INDEPENDENT lanes
                           (sequence + structure, sequence + curated,
                            structure + curated, or all three)
      - 'sequence_only'  : only sequence-based DBs
      - 'structure_only' : only Foldseek structural discovery
      - 'curated_only'   : only CornCyc curation (no homology hit)
    """
    s = set(sources or [])
    has_struct  = STRUCTURAL_SOURCE in s
    has_curated = CURATED_SOURCE in s
    has_seq     = bool(s - {STRUCTURAL_SOURCE, CURATED_SOURCE})
    lanes_hit   = int(has_seq) + int(has_struct) + int(has_curated)
    if lanes_hit >= 2:
        return "consensus"
    if has_struct:  return "structure_only"
    if has_seq:     return "sequence_only"
    if has_curated: return "curated_only"
    return "sequence_only"  # fallback for unknown source labels


def lane_label(cls: str) -> str:
    return {
        "consensus":      "Consensus (≥2 lanes)",
        "sequence_only":  "Sequence-based",
        "structure_only": "Structure-based (Foldseek)",
        "curated_only":   "CornCyc curated",
    }.get(cls, cls)


def split_orthologs_by_lane(orthologs: List[Dict]) -> Dict[str, List[Dict]]:
    """Bucket orthologs by classify_ortholog. Each bucket pre-sorted by similarity desc."""
    buckets = {"consensus": [], "sequence_only": [], "structure_only": [], "curated_only": []}
    for o in orthologs or []:
        buckets[classify_ortholog(o.get("sources", []))].append(o)
    for k in buckets:
        buckets[k].sort(key=lambda o: o.get("similarity_score", 0), reverse=True)
    return buckets


def join_enrichment(orthologs: List[Dict], res: Dict) -> List[Dict]:
    """
    Build per-ortholog rows that flatten everything we know about each maize
    gene: discovery lane, sources, sequence/structural similarity, enrichment
    (pLDDT, expression, Pfam domain), and the top Phase 7 Compara ortholog.
    """
    targets_by_gene  = {t["maize_gene_model"]: t for t in (res.get("targets") or [])}
    domain_by_gene   = {d["maize_gene_model"]: d for d in (res.get("domain_targets") or [])}
    advanced_by_gene = {a["maize_gene_model"]: a for a in (res.get("advanced_homology_targets") or [])}
    gene_meta_map    = res.get("maize_gene_metadata") or {}
    phytozome_map    = res.get("phytozome_metadata") or {}

    # CornCyc per-gene annotation index (gene → {pathway_ids, pathway_names, reactions, ec})
    cc_ann = res.get("corncyc_annotation") or {}
    cc_gene_index: Dict[str, dict] = {}
    cc_pathway_name = {p["id"]: p["common_name"] for p in (cc_ann.get("pathways") or [])}
    for g in (cc_ann.get("maize_genes") or []):
        cc_gene_index[g["gene"]] = {
            "pathway_ids":   g.get("pathways", []),
            "pathway_names": [cc_pathway_name.get(pid, pid) for pid in g.get("pathways", [])],
            "reactions":     g.get("reactions", []),
            "ec_numbers":    g.get("ec_numbers", []),
        }

    # Reaction context that's shared across all orthologs of the query
    reactions = res.get("reactions") or []
    top_pathway_name = ""
    top_pathway_id = ""
    pathway_summary = ""
    if reactions and reactions[0].get("pathway_names"):
        top_pathway_name = reactions[0]["pathway_names"][0]
        if reactions[0].get("pathway_ids"):
            top_pathway_id = reactions[0]["pathway_ids"][0] if reactions[0]["pathway_ids"] else ""
        pathway_summary = " | ".join(reactions[0]["pathway_names"][:5])
    top_reaction = reactions[0] if reactions else {}

    rows = []
    for o in orthologs:
        sources = o.get("sources", []) or []
        cls = classify_ortholog(sources)
        gene = o.get("maize_gene_model", "")
        sim = float(o.get("similarity_score", 0))

        # similarity_score is TM*100 for structural hits, %identity-ish for sequence hits
        seq_similarity = sim if cls != "structure_only" else None
        structural_tm = sim / 100.0 if cls in ("structure_only", "consensus") else None

        t = targets_by_gene.get(gene)
        d = domain_by_gene.get(gene)
        adv = advanced_by_gene.get(gene) or {}
        top_compara = (adv.get("ensembl_orthologs") or [{}])[0] if adv else {}

        # If enriched (Phase 5 ran), the TM-score lives on the ValidatedTarget;
        # for structural hits this matches the discovery TM, for sequence-only it's a fresh value
        # (or 0.0 if the cheap-path skipped the 1-to-1 Foldseek).
        if t and t.get("tm_score") not in (None, 0, 0.0):
            structural_tm = float(t["tm_score"])

        expr_count = (t or {}).get("n_expression_experiments", 0) or 0
        expr_experiments = (t or {}).get("expression_experiments", []) or []
        tissues = (t or {}).get("tissue_expression_fpkm") or {}

        meta = gene_meta_map.get(gene) or {}
        cc_ev = cc_gene_index.get(gene) or {}
        pz = phytozome_map.get(gene) or {}
        # Split Gramene synonyms by version family so the CSV can be filtered cleanly.
        all_syns = meta.get("synonyms", []) or []
        v3_ids = [s for s in all_syns if s.startswith("GRMZM")]
        v4_ids = [s for s in all_syns if s.startswith("Zm00001d")]
        other_syns = [s for s in all_syns if s not in v3_ids and s not in v4_ids]
        source_ev = o.get("source_evidence") or {}
        # Human-readable summary like: "Ensembl Compara plants (82.5% id, ortholog_one2one);
        #   Foldseek 3D (qTM=0.91 tTM=0.88 prob=0.97); CornCyc curated (3 reactions, 2 pathways)"
        methods_summary_parts = []
        for src_key, src_label in KNOWN_SOURCES:
            ev = source_ev.get(src_key)
            if not ev:
                continue
            detail = ""
            if src_key in ("Ensembl", "EnsemblPanHomology"):
                detail = f"{ev.get('pct_identity', 0):.1f}% id, {ev.get('ortholog_type','ortholog')}"
            elif src_key == "Local_HMMER":
                detail = f"E={ev.get('e_value','?'):.1e}, bit={ev.get('bit_score','?')}"
            elif src_key == "Foldseek-structural":
                bits = []
                if ev.get("qtm") is not None: bits.append(f"qTM={ev['qtm']:.2f}")
                if ev.get("ttm") is not None: bits.append(f"tTM={ev['ttm']:.2f}")
                if ev.get("prob") is not None: bits.append(f"prob={ev['prob']:.2f}")
                detail = " ".join(bits)
            elif src_key == "PLAZA":
                detail = f"orthogroup={ev.get('orthogroup','?')}"
            elif src_key == "CornCyc":
                detail = f"{ev.get('n_reactions','?')} reactions, {len(ev.get('pathway_ids',[]))} pathways"
            methods_summary_parts.append(f"{src_label} ({detail})" if detail else src_label)
        methods_summary = " ; ".join(methods_summary_parts)

        rows.append({
            "gene": gene,
            "gene_symbol":      meta.get("symbol", ""),
            "gene_description": meta.get("description", ""),
            "gene_synonyms":    all_syns,
            "gene_v3_ids":      v3_ids,        # legacy v3 NAM (GRMZM2G*)
            "gene_v4_ids":      v4_ids,        # legacy v4 NAM (Zm00001d*)
            "gene_other_synonyms": other_syns, # community names like 'sudh4'
            "phytozome_description":   pz.get("description", ""),
            "phytozome_panther_ids":   pz.get("panther_ids", []),
            "phytozome_panther_descs": pz.get("panther_descs", []),
            # Full per-source evidence dict (preserved as JSON-like for advanced users)
            "source_evidence":      source_ev,
            "methods_summary":      methods_summary,
            "consensus_class": cls,
            "sources": sources,
            "num_sources": len(sources),
            "consensus_score": o.get("consensus_score", 1),
            "query_uniprot_id": o.get("query_uniprot_id", ""),
            "plaza_orthogroup": o.get("plaza_orthogroup", ""),
            "sequence_similarity": seq_similarity,
            "structural_tm_score": structural_tm,
            "plddt": (t or {}).get("plddt"),
            "tissue_count": len(tissues),
            "tissue_summary": "; ".join(f"{k}={v}" for k, v in tissues.items()),
            "n_expression_experiments": expr_count,
            "expression_experiments":   expr_experiments,
            "enrichment_kind":          (t or {}).get("enrichment_kind", ""),
            "pfam_domain_id":   (d or {}).get("pfam_domain_id", ""),
            "pfam_domain_name": (d or {}).get("pfam_domain_name", ""),
            "pfam_start":       (d or {}).get("domain_start", ""),
            "pfam_end":         (d or {}).get("domain_end", ""),
            "domain_tm_score":  (d or {}).get("domain_tm_score"),
            "top_compara_species":     top_compara.get("species", ""),
            "top_compara_gene":        top_compara.get("gene_id", ""),
            "top_compara_pct_id":      top_compara.get("percent_identity", ""),
            "num_compara_orthologs":   len(adv.get("ensembl_orthologs", []) or []),
            "enriched": t is not None,
            # CornCyc per-gene evidence
            "corncyc_pathway_ids":   cc_ev.get("pathway_ids", []),
            "corncyc_pathway_names": cc_ev.get("pathway_names", []),
            "corncyc_reactions":     cc_ev.get("reactions", []),
            "corncyc_ec_numbers":    cc_ev.get("ec_numbers", []),
            # Query-level reaction context (same for every ortholog of this query)
            "top_reaction_rhea_id": top_reaction.get("rhea_id", ""),
            "top_reaction_equation": top_reaction.get("equation") or top_reaction.get("label", ""),
            "top_reaction_ec": "|".join(top_reaction.get("ec_numbers", []) or []),
            "top_reaction_is_transport": top_reaction.get("is_transport"),
            "top_kegg_pathway_name": top_pathway_name,
            "top_kegg_pathway_id": top_pathway_id,
            "kegg_pathway_summary": pathway_summary,
        })

    # Order: consensus first, then structure_only, then sequence_only;
    # within each group, by best of (structural_tm, sequence_similarity)
    rank = {"consensus": 0, "structure_only": 1, "sequence_only": 2}
    def sort_key(r):
        score = max(
            (r["structural_tm_score"] or 0) * 100,
            (r["sequence_similarity"] or 0),
        )
        return (rank.get(r["consensus_class"], 99), -score)
    rows.sort(key=sort_key)
    return rows


# ---------- Per-result summary helpers ----------

def summarize_result(res: Dict) -> Dict:
    """Build a one-line summary record for the executive header."""
    inp = res.get("input_data", {}) or {}
    qtype = inp.get("query_type", "mz")
    if qtype == "mz":
        query_repr = f"m/z {inp.get('mz', '?')} ({inp.get('mode', '?')})"
    elif qtype == "chemical":
        query_repr = inp.get("chemical_name", "?")
    elif qtype == "ec":
        query_repr = f"EC {inp.get('ec_number', '?')}"
    else:
        query_repr = qtype

    ce = res.get("chemical_entity") or {}
    orthologs = res.get("orthologs") or []
    targets = res.get("targets") or []
    domain_targets = res.get("domain_targets") or []
    adv = res.get("advanced_homology_targets") or []

    lanes = split_orthologs_by_lane(orthologs)
    top_target = targets[0] if targets else None
    # Headline gene: prefer a consensus hit if there is one, else best of any lane
    headline = (lanes["consensus"] or lanes["structure_only"] or lanes["curated_only"] or lanes["sequence_only"] or [None])[0]
    cc = res.get("corncyc_annotation") or {}

    status = (
        "error" if res.get("error")
        else ("targets" if (orthologs or targets) else "no_targets")
    )

    return {
        "query_type": qtype,
        "query_repr": query_repr,
        "chebi_id": ce.get("chebi_id", ""),
        "pubchem_cid": ce.get("pubchem_cid", ""),
        "monoisotopic_mass": ce.get("monoisotopic_mass", ""),
        "num_reactions": len(res.get("reactions") or []),
        "num_proteins": len(res.get("proteins") or []),
        "num_orthologs": len(orthologs),
        "num_consensus": len(lanes["consensus"]),
        "num_sequence_only": len(lanes["sequence_only"]),
        "num_structure_only": len(lanes["structure_only"]),
        "num_curated_only": len(lanes["curated_only"]),
        "corncyc_pathways": cc.get("n_pathways", 0),
        "corncyc_maize_genes": cc.get("n_maize_genes", 0),
        "corncyc_top_pathway": (cc.get("pathways") or [{}])[0].get("common_name", ""),
        "num_targets": len(targets),
        "num_domain_targets": len(domain_targets),
        "num_advanced_targets": len(adv),
        "top_gene": top_target.get("maize_gene_model") if top_target else (headline.get("maize_gene_model") if headline else ""),
        "top_tm": top_target.get("tm_score") if top_target else "",
        "top_plddt": top_target.get("plddt") if top_target else "",
        "headline_class": classify_ortholog(headline.get("sources", [])) if headline else "",
        "status": status,
        "error": res.get("error", ""),
    }

# ---------- HTML template ----------

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Maysquery Pipeline Report</title>
<style>
:root { --bg:#f9fafb; --panel:#fff; --ink:#111827; --muted:#6b7280; --line:#e5e7eb; --accent:#2563eb; --ok:#15803d; --warn:#b45309; --err:#b91c1c; }
* { box-sizing: border-box; }
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 24px; background: var(--bg); color: var(--ink); margin:0; }
.container { max-width: 1100px; margin: 0 auto; }
header.report-head { background: linear-gradient(135deg, #1e3a8a 0%, #6366f1 100%); color: #fff; border-radius: 12px; padding: 28px 32px; margin-bottom: 24px; box-shadow: 0 6px 20px rgba(30,58,138,0.25); }
header.report-head h1 { margin: 0 0 6px; font-size: 1.6rem; }
header.report-head .meta { opacity: 0.85; font-size: 0.9rem; }
.summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin-top: 18px; }
.summary-cell { background: rgba(255,255,255,0.12); border-radius: 8px; padding: 10px 14px; }
.summary-cell .lbl { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.5px; opacity: 0.8; }
.summary-cell .val { font-size: 1.3rem; font-weight: 600; margin-top: 2px; }
.summary-table { width: 100%; border-collapse: collapse; margin-top: 18px; background: rgba(255,255,255,0.08); border-radius: 8px; overflow: hidden; font-size: 0.85rem; }
.summary-table th, .summary-table td { padding: 8px 12px; text-align: left; border-bottom: 1px solid rgba(255,255,255,0.15); }
.summary-table th { background: rgba(0,0,0,0.15); font-weight: 600; }
.summary-table tr:last-child td { border-bottom: none; }
.status-targets { color: #86efac; font-weight: 600; }
.status-no_targets { color: #fcd34d; }
.status-error { color: #fca5a5; }

.query-block { background: var(--panel); border-radius: 10px; padding: 24px 28px; margin-bottom: 24px; box-shadow: 0 2px 6px rgba(0,0,0,0.06); border-left: 4px solid var(--accent); }
.query-block h2 { margin: 0 0 14px; color: var(--accent); font-size: 1.15rem; }
.query-block h3 { margin: 22px 0 10px; color: #374151; font-size: 1rem; border-bottom: 1px solid var(--line); padding-bottom: 6px; }
.kvgrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 8px 18px; margin: 8px 0 0; font-size: 0.9rem; }
.kvgrid .k { color: var(--muted); font-size: 0.78rem; }
.kvgrid .v { font-weight: 500; }

table.data { border-collapse: collapse; width: 100%; font-size: 0.86rem; margin-top: 8px; }
table.data th, table.data td { border: 1px solid var(--line); padding: 8px 10px; text-align: left; vertical-align: top; }
table.data th { background: #f3f4f6; font-weight: 600; color: #4b5563; font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.3px; }
table.data tr:nth-child(even) td { background: #fafafa; }

a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 9999px; font-size: 0.72rem; font-weight: 600; background: #dbeafe; color: #1e40af; margin-right: 4px; }
.badge-mode { background: #ede9fe; color: #5b21b6; }
.badge-transport { background: #e0e7ff; color: #3730a3; }
.badge-metabolic { background: #f1f5f9; color: #475569; }
.badge-enzyme { background: #fee2e2; color: #991b1b; }
.badge-transporter { background: #dcfce7; color: #166534; }
.badge-receptor { background: #fef3c7; color: #92400e; }
.badge-ok { background: #dcfce7; color: #166534; }
.badge-warn { background: #fef3c7; color: #92400e; }
.badge-err { background: #fee2e2; color: #991b1b; }
.badge-info { background: #dbeafe; color: #1e40af; }
.badge-consensus      { background: linear-gradient(90deg,#a78bfa 0%,#34d399 100%); color: #fff; font-weight:700; }
.badge-sequence_only  { background: #fde68a; color: #78350f; }
.badge-structure_only { background: #c7d2fe; color: #312e81; }
.badge-curated_only   { background: #fbcfe8; color: #831843; }
.badge-corncyc        { background: #f0fdf4; color: #166534; border: 1px solid #86efac; font-size: 0.7rem; }
.badge-src { background: #e5e7eb; color: #374151; margin-right:2px; font-size:0.7rem; }
.row-consensus td      { background: rgba(167,139,250,0.06) !important; }
.row-structure_only td { background: rgba(199,210,254,0.10) !important; }
.row-sequence_only td  { background: rgba(253,230,138,0.08) !important; }
.row-curated_only td   { background: rgba(251,207,232,0.10) !important; }
.corncyc-block { background: #f0fdf4; border: 1px solid #bbf7d0; border-radius: 6px; padding: 10px 14px; margin-bottom: 14px; }
.corncyc-block h4 { margin: 0 0 6px; color: #166534; }
.corncyc-pathway-row { padding: 4px 0; border-bottom: 1px dashed #bbf7d0; }
.corncyc-pathway-row:last-child { border-bottom: none; }
.dim { color: var(--muted); font-weight: 400; font-size: 0.85em; }
.gene-meta { color: var(--ink); font-weight: 500; font-size: 0.86em; margin-left: 4px; }

/* Per-query executive card */
.exec-card { background: linear-gradient(135deg, #eef2ff 0%, #faf5ff 100%); border: 1px solid #c7d2fe; border-radius: 8px; padding: 14px 18px; margin: 6px 0 22px; }
.exec-kv { display: flex; gap: 12px; padding: 4px 0; align-items: baseline; border-bottom: 1px dashed #c7d2fe; }
.exec-kv:last-child { border-bottom: none; }
.exec-kv-top { margin-top: 6px; padding-top: 10px; border-top: 2px solid #c7d2fe; border-bottom: none; }
.exec-k { color: #4338ca; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.4px; min-width: 180px; flex-shrink: 0; font-weight: 600; }
.exec-v { color: var(--ink); font-size: 0.95rem; flex: 1; }
.exec-stat { color: var(--ink); }

/* Per-query collapsibles */
.r-drill { border: 1px solid #e5e7eb; border-radius: 6px; margin-bottom: 8px; background: #fff; }
.r-drill > summary { cursor: pointer; padding: 9px 14px; user-select: none; list-style: none; }
.r-drill > summary::-webkit-details-marker { display: none; }
.r-drill > summary::marker { display: none; content: ''; }
.r-drill > summary::before { content: '▸ '; color: var(--muted); display: inline-block; transition: transform 0.15s; }
.r-drill[open] > summary::before { content: '▾ '; }
.r-drill > summary:hover { background: #f3f4f6; }
.r-drill[open] > summary { border-bottom: 1px solid #e5e7eb; background: #fafafa; }
.r-drill > *:not(summary) { padding: 0 14px 14px; }
.r-drill > table.data { margin: 12px; width: calc(100% - 24px); }
.r-drill-main { border-color: var(--accent); }
.lane-table th, .lane-table td { vertical-align: top; }
.lane-table .dim { color: var(--muted); }
.mono { font-family: 'JetBrains Mono', Consolas, monospace; font-size: 0.82rem; }
.empty { color: var(--muted); font-style: italic; padding: 8px 0; }
.tissue-list { display: flex; flex-wrap: wrap; gap: 6px; }
.tissue { background: #eff6ff; color: #1d4ed8; padding: 2px 8px; border-radius: 4px; font-size: 0.78rem; }
footer { text-align: center; color: var(--muted); font-size: 0.78rem; margin-top: 30px; }
.tm-score, .plddt-score { font-family: 'JetBrains Mono', Consolas, monospace; font-weight: 600; }
</style>
</head>
<body>
<div class="container">

<header class="report-head">
    <h1>Maysquery — Metabolite-to-Gene Pipeline Report</h1>
    <div class="meta">
        Generated {{ generated_at }}
        &nbsp;·&nbsp; {{ summaries|length }} {{ 'query' if summaries|length == 1 else 'queries' }}
        &nbsp;·&nbsp; {{ totals.targets }} validated maize target{{ 's' if totals.targets != 1 else '' }}
    </div>

    <div class="summary-grid">
        <div class="summary-cell"><div class="lbl">Queries Run</div><div class="val">{{ summaries|length }}</div></div>
        <div class="summary-cell"><div class="lbl">Reactions</div><div class="val">{{ totals.reactions }}</div></div>
        <div class="summary-cell"><div class="lbl">Proteins</div><div class="val">{{ totals.proteins }}</div></div>
        <div class="summary-cell"><div class="lbl">Total Maize Homologs</div><div class="val">{{ totals.orthologs }}</div></div>
        <div class="summary-cell" style="background: linear-gradient(135deg, rgba(167,139,250,0.4), rgba(52,211,153,0.4));"><div class="lbl">Consensus (Seq + Struct)</div><div class="val">{{ totals.consensus }}</div></div>
        <div class="summary-cell"><div class="lbl">Sequence-only</div><div class="val">{{ totals.sequence_only }}</div></div>
        <div class="summary-cell"><div class="lbl">Structure-only</div><div class="val">{{ totals.structure_only }}</div></div>
        <div class="summary-cell"><div class="lbl">CornCyc curated-only</div><div class="val">{{ totals.curated_only }}</div></div>
        <div class="summary-cell"><div class="lbl">Enriched (pLDDT + expr)</div><div class="val">{{ totals.targets }}</div></div>
    </div>

    <table class="summary-table">
        <thead><tr>
            <th>#</th><th>Query</th><th>Resolved</th><th>Headline Gene</th><th>Lane</th><th>Consensus</th><th>Seq-only</th><th>Struct-only</th><th>Status</th>
        </tr></thead>
        <tbody>
        {% for s in summaries %}
            <tr>
                <td>{{ loop.index }}</td>
                <td class="mono">{{ s.query_repr }}</td>
                <td>{{ s.chebi_id or '—' }}</td>
                <td>{% if s.top_gene %}{{ gene_label(s.top_gene, results[loop.index0].maize_gene_metadata)|safe }}{% else %}—{% endif %}</td>
                <td>{% if s.headline_class %}<span class="badge badge-{{ s.headline_class }}">{{ lane_label(s.headline_class) }}</span>{% else %}—{% endif %}</td>
                <td class="mono">{{ s.num_consensus }}</td>
                <td class="mono">{{ s.num_sequence_only }}</td>
                <td class="mono">{{ s.num_structure_only }}</td>
                <td class="status-{{ s.status }}">{{ s.status.replace('_', ' ') }}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
</header>

{% for res in results %}
{% set s = summaries[loop.index0] %}
{% set rows = enriched_rows(res.orthologs or [], res) %}
{% set lanes = split_lanes(res.orthologs or []) %}
{% set top = rows[0] if rows else None %}
<div class="query-block">
    <h2>Query #{{ loop.index }}: {{ s.query_repr }}</h2>

    {# ---------- Executive summary (always expanded) ---------- #}
    <div class="exec-card">
        <div class="exec-kv">
            <span class="exec-k">Resolved entity</span>
            <span class="exec-v">
                {% if res.chemical_entity %}<a href="{{ url_chebi(res.chemical_entity.chebi_id) }}" target="_blank">{{ res.chemical_entity.chebi_id }}</a>
                    {% if res.chemical_entity.pubchem_cid %} · <a href="{{ url_pubchem(res.chemical_entity.pubchem_cid) }}" target="_blank">CID {{ res.chemical_entity.pubchem_cid }}</a>{% endif %}
                    · mass {{ '%.4f'|format(res.chemical_entity.monoisotopic_mass) }}
                {% else %}—{% endif %}
            </span>
        </div>
        <div class="exec-kv">
            <span class="exec-k">Counts</span>
            <span class="exec-v">
                <span class="exec-stat"><b>{{ res.reactions|length }}</b> reactions</span> ·
                <span class="exec-stat"><b>{{ res.proteins|length }}</b> pan-life proteins</span> ·
                <span class="exec-stat"><b>{{ rows|length }}</b> maize candidates</span>
                ({{ s.num_consensus }} consensus, {{ s.num_structure_only }} structure-only, {{ s.num_sequence_only }} sequence-only)
            </span>
        </div>
        {% if res.reactions and res.reactions[0].pathway_names %}
        <div class="exec-kv">
            <span class="exec-k">Top KEGG pathway</span>
            <span class="exec-v">{{ res.reactions[0].pathway_names[0] }}{% if res.reactions[0].pathway_names|length > 1 %} <span class="dim">+{{ res.reactions[0].pathway_names|length - 1 }} more</span>{% endif %}</span>
        </div>
        {% endif %}
        {% if top %}
        <div class="exec-kv exec-kv-top">
            <span class="exec-k">★ Top maize target</span>
            <span class="exec-v">
                {{ gene_label(top.gene, res.maize_gene_metadata)|safe }}
                <span class="badge badge-{{ top.consensus_class }}">{{ lane_label(top.consensus_class) }}</span>
                {% if top.structural_tm_score is not none %}· TM <b>{{ '%.2f'|format(top.structural_tm_score) }}</b>{% endif %}
                {% if top.sequence_similarity is not none %}· seq <b>{{ '%.1f'|format(top.sequence_similarity) }}%</b>{% endif %}
                {% if top.plddt is not none %}· pLDDT <b>{{ '%.1f'|format(top.plddt) }}</b>{% endif %}
                {% if top.n_expression_experiments %}· expressed in <b>{{ top.n_expression_experiments }}</b> Atlas exp{{ 's' if top.n_expression_experiments != 1 else '' }}{% endif %}
            </span>
        </div>
        {% endif %}
    </div>

    {# ---------- Phase 1 (always shown — small) ---------- #}
    {% if res.chemical_entity or res.input_data.query_type == 'ec' %}
    <details open class="r-drill"><summary><b>Phase 1 — Chemical Entity</b></summary>
    {% if res.chemical_entity %}
    <div class="kvgrid">
        <div><div class="k">ChEBI</div><div class="v"><a href="{{ url_chebi(res.chemical_entity.chebi_id) }}" target="_blank">{{ res.chemical_entity.chebi_id }}</a></div></div>
        {% if res.chemical_entity.pubchem_cid %}<div><div class="k">PubChem CID</div><div class="v"><a href="{{ url_pubchem(res.chemical_entity.pubchem_cid) }}" target="_blank">{{ res.chemical_entity.pubchem_cid }}</a></div></div>{% endif %}
        <div><div class="k">Monoisotopic Mass</div><div class="v mono">{{ res.chemical_entity.monoisotopic_mass }}</div></div>
        {% if res.chemical_entity.smiles %}<div><div class="k">SMILES</div><div class="v mono" style="word-break: break-all;">{{ res.chemical_entity.smiles }}</div></div>{% endif %}
    </div>
    {% else %}
    <div class="kvgrid"><div><div class="k">EC Number</div><div class="v mono">{{ res.input_data.ec_number }}</div></div></div>
    {% endif %}
    </details>
    {% endif %}

    {# ---------- CornCyc maize-specific pathway context ---------- #}
    {% if res.corncyc_annotation %}
    <details open class="r-drill"><summary><b>CornCyc maize pathway context</b> <span class="dim">{{ res.corncyc_annotation.n_pathways }} pathway(s), {{ res.corncyc_annotation.n_maize_genes }} annotated maize gene(s) · PMN CornCyc {{ res.corncyc_annotation.version }}</span></summary>
    <div class="corncyc-block">
        <div style="font-size:0.85rem; color:#166534; margin-bottom: 8px;">
            Curated PlantCyc/CornCyc pathways involving
            {% for c in res.corncyc_annotation.compounds %}<b>{{ c.name }}</b>{% if not loop.last %}, {% endif %}{% endfor %}.
            Maize gene annotations from <a href="https://www.plantcyc.org/" target="_blank">Plant Metabolic Network</a>.
        </div>
        {% for p in res.corncyc_annotation.pathways[:25] %}
            <div class="corncyc-pathway-row">
                <a href="https://pmn.plantcyc.org/pathway?orgid=CORN&id={{ p.id }}" target="_blank"><b>{{ p.common_name|safe }}</b></a>
                <span class="dim mono">{{ p.id }}</span>
                <span class="dim">· {{ p.reactions_touching_compound|length }} matching reaction(s), {{ p.maize_genes|length }} maize gene(s)</span>
                {% if p.maize_genes %}<div class="mono dim" style="font-size:0.78rem; margin-top:2px;">
                    {% for g in p.maize_genes[:6] %}{{ gene_label(g, res.maize_gene_metadata)|safe }}{% if not loop.last %} · {% endif %}{% endfor %}{% if p.maize_genes|length > 6 %} <span class="dim">+{{ p.maize_genes|length - 6 }} more</span>{% endif %}
                </div>{% endif %}
            </div>
        {% endfor %}
        {% if res.corncyc_annotation.pathways|length > 25 %}
            <div class="dim" style="margin-top:6px;">+{{ res.corncyc_annotation.pathways|length - 25 }} more pathway(s) — see CSV for the full list.</div>
        {% endif %}
    </div>
    </details>
    {% endif %}

    {# ---------- Phase 2: reactions (collapsed by default) ---------- #}
    {% if res.reactions %}
    <details class="r-drill"><summary><b>Reaction Networks (Rhea)</b> <span class="dim">{{ res.reactions|length }} reactions</span></summary>
    <table class="data">
        <thead><tr><th>Rhea ID</th><th>Equation</th><th>EC</th><th>Type</th><th>Pathway</th></tr></thead>
        <tbody>
        {% for r in res.reactions %}
            <tr>
                <td><a href="{{ url_rhea(r.rhea_id) }}" target="_blank">{{ r.rhea_id }}</a></td>
                <td class="mono" style="font-size:0.78rem; max-width:380px;">{{ r.equation or r.label or '—' }}</td>
                <td class="mono">{% if r.ec_numbers %}{% for e in r.ec_numbers %}<span class="badge badge-ec">{{ e }}</span>{% endfor %}{% else %}—{% endif %}</td>
                <td>{% if r.is_transport %}<span class="badge badge-transport">Transport</span>{% else %}<span class="badge badge-metabolic">Metabolic</span>{% endif %}</td>
                <td>{% if r.pathway_names %}{% for p in r.pathway_names[:2] %}<a href="{{ url_kegg_pathway(r.pathway_ids[loop.index0] if r.pathway_ids and loop.index0 < r.pathway_ids|length else p) }}" target="_blank">{{ p }}</a>{% if not loop.last %}, {% endif %}{% endfor %}{% if r.pathway_names|length > 2 %} <span class="dim">+{{ r.pathway_names|length - 2 }} more</span>{% endif %}{% else %}—{% endif %}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    </details>
    {% endif %}

    {# ---------- Phase 3: proteins (collapsed) ---------- #}
    {% if res.proteins %}
    <details class="r-drill"><summary><b>Pan-life protein pool (UniProt)</b> <span class="dim">{{ res.proteins|length }} unique</span></summary>
    <table class="data">
        <thead><tr><th>UniProt</th><th>Category</th><th>GO Terms (top 3)</th></tr></thead>
        <tbody>
        {% for p in res.proteins %}
            <tr>
                <td><a href="{{ url_uniprot(p.uniprot_accession) }}" target="_blank">{{ p.uniprot_accession }}</a></td>
                <td><span class="badge badge-{{ p.category.lower() }}">{{ p.category }}</span></td>
                <td>{% if p.go_terms %}{% for go in p.go_terms[:3] %}<a href="{{ url_quickgo(go.id) }}" target="_blank">{{ go.id }}</a> <span style="color:#6b7280;">{{ go.name }}</span>{% if not loop.last %}<br>{% endif %}{% endfor %}{% else %}<span class="empty">none</span>{% endif %}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    </details>
    {% endif %}

    {# ---------- Maize candidates (THE main payload — open by default) ---------- #}
    {% if rows %}
    <details open class="r-drill r-drill-main"><summary><b>Maize candidates</b> <span class="dim">{{ rows|length }} total · {{ lanes.consensus|length }} consensus · {{ lanes.structure_only|length }} structure-only · {{ lanes.sequence_only|length }} sequence-only</span></summary>
    <table class="data lane-table">
        <thead><tr>
            <th>Lane</th><th>Maize Gene</th><th>Query Protein</th>
            <th>TM</th><th>Seq</th><th>pLDDT</th>
            <th>Expression</th><th>Pfam Domain</th><th>Top Compara</th>
        </tr></thead>
        <tbody>
        {% for r in rows %}
            <tr class="row-{{ r.consensus_class }}">
                <td><span class="badge badge-{{ r.consensus_class }}">{{ lane_label(r.consensus_class) }}</span></td>
                <td>{{ gene_label(r.gene, res.maize_gene_metadata)|safe }}<div class="dim mono" style="font-size:0.72rem;">{% for s in r.sources %}{{ s }}{% if not loop.last %} · {% endif %}{% endfor %}</div>{% if r.methods_summary %}<div class="dim" style="font-size:0.7rem; margin-top:2px; line-height:1.3;">{{ r.methods_summary }}</div>{% endif %}</td>
                <td>{% if r.query_uniprot_id %}<a href="{{ url_uniprot(r.query_uniprot_id) }}" target="_blank">{{ r.query_uniprot_id }}</a>{% else %}—{% endif %}</td>
                <td class="mono"><b>{% if r.structural_tm_score is not none %}{{ '%.2f'|format(r.structural_tm_score) }}{% else %}—{% endif %}</b></td>
                <td class="mono">{% if r.sequence_similarity is not none %}{{ '%.1f'|format(r.sequence_similarity) }}{% else %}—{% endif %}</td>
                <td class="mono">{% if r.plddt is not none %}{{ '%.1f'|format(r.plddt) }}{% else %}<span class="empty">—</span>{% endif %}</td>
                <td class="mono" style="font-size:0.78rem;">{% if r.n_expression_experiments %}<b>{{ r.n_expression_experiments }}</b> exps<div class="dim" style="font-size:0.7rem;">{% for e in r.expression_experiments[:4] %}<a href="https://www.ebi.ac.uk/gxa/experiments/{{ e }}" target="_blank">{{ e }}</a>{% if not loop.last %}, {% endif %}{% endfor %}{% if r.n_expression_experiments > 4 %} +{{ r.n_expression_experiments - 4 }}{% endif %}</div>{% else %}<span class="empty">—</span>{% endif %}</td>
                <td>{% if r.pfam_domain_id %}<a href="{{ url_pfam(r.pfam_domain_id) }}" target="_blank">{{ r.pfam_domain_id }}</a> {{ r.pfam_domain_name }}<div class="dim mono" style="font-size:0.7rem;">res {{ r.pfam_start }}–{{ r.pfam_end }} · TM {{ '%.2f'|format(r.domain_tm_score) }}</div>{% else %}<span class="empty">—</span>{% endif %}</td>
                <td>{% if r.top_compara_gene %}<a href="{{ url_ensembl_plants_gene(r.top_compara_species, r.top_compara_gene) }}" target="_blank" class="mono">{{ r.top_compara_species }}</a><div class="dim mono" style="font-size:0.7rem;">{{ r.top_compara_gene }} · {{ '%.1f'|format(r.top_compara_pct_id) }}%{% if r.num_compara_orthologs > 1 %} · +{{ r.num_compara_orthologs - 1 }}{% endif %}</div>{% else %}<span class="empty">—</span>{% endif %}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    </details>
    {% endif %}

    {# ---------- Phase 7: full Compara details (collapsed) ---------- #}
    {% if res.advanced_homology_targets %}
    <details class="r-drill"><summary><b>Pan-Plant Compara orthologs</b> <span class="dim">{{ res.advanced_homology_targets|length }} target{{ 's' if res.advanced_homology_targets|length != 1 else '' }}</span></summary>
    {% for at in res.advanced_homology_targets %}
        {% if at.ensembl_orthologs %}
        <div style="margin-bottom: 14px;">
            <div style="margin-bottom: 6px;">{{ gene_label(at.maize_gene_model, res.maize_gene_metadata)|safe }} <span class="dim">{{ at.ensembl_orthologs|length }} orthologs</span></div>
            <table class="data">
                <thead><tr><th>Species</th><th>Gene ID</th><th>Protein ID</th><th>% Identity</th></tr></thead>
                <tbody>
                {% for o in at.ensembl_orthologs %}
                    <tr>
                        <td>{{ o.species }}</td>
                        <td><a href="{{ url_ensembl_plants_gene(o.species, o.gene_id) }}" target="_blank" class="mono">{{ o.gene_id }}</a></td>
                        <td class="mono">{{ o.protein_id or '—' }}</td>
                        <td class="mono">{{ '%.1f'|format(o.percent_identity) }}</td>
                    </tr>
                {% endfor %}
                </tbody>
            </table>
        </div>
        {% endif %}
    {% endfor %}
    </details>
    {% endif %}

    {# ---------- Execution log (collapsed) ---------- #}
    {% if res.execution_logs %}
    <details class="r-drill"><summary><b>Execution log</b> <span class="dim">{{ res.execution_logs|length }} events</span></summary>
    <table class="data">
        <thead><tr><th>Phase</th><th>Database</th><th>Status</th><th>Hits</th><th>Message</th></tr></thead>
        <tbody>
        {% for log in res.execution_logs %}
            <tr>
                <td>P{{ log.phase }}</td>
                <td>{{ log.database }}</td>
                <td><span class="badge badge-{{ 'ok' if log.status == 'success' else ('err' if log.status == 'error' else ('warn' if log.status == 'warning' else 'info')) }}">{{ log.status }}</span></td>
                <td class="mono">{{ log.hits }}</td>
                <td>{{ log.message }}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    </details>
    {% endif %}

    {% if res.error %}
    <h3>Error</h3>
    <div class="empty" style="color: var(--err);">{{ res.error }}</div>
    {% endif %}
</div>
{% endfor %}

{% if totals.corncyc_pathways %}
<footer class="corncyc-attribution" style="background:#f0fdf4; border:1px solid #bbf7d0; border-radius:8px; padding:14px 18px; margin-top:24px; font-size:0.82rem; color:#166534; line-height:1.5;">
    <b style="font-size:0.9rem;">CornCyc / Plant Metabolic Network attribution</b><br>
    This report contains data derived from <b>CornCyc</b> (Plant Metabolic
    Network, PMN). CornCyc is authored by Charles Hawkins, Bo Xue, and
    Seung Yon Rhee, copyright © PMN / Donald Danforth Plant Science Center.
    Original database: <a href="https://www.plantcyc.org/databases/corncyc/" target="_blank" style="color:#15803d;">plantcyc.org/databases/corncyc</a>.
    Maysquery uses CornCyc as a read-only reference (no in-place edits);
    full list of derived presentations and modifications is documented in
    the <code>CORNCYC_ATTRIBUTION.txt</code> file shipped with the Maysquery
    source. For canonical CornCyc pathway pages and citations, follow the
    PMN PlantCyc links beside each pathway above.
</footer>
{% endif %}
<footer>Generated by Maysquery · Metabolite-to-Gene Pipeline · {{ generated_at }}</footer>
</div>
</body>
</html>
"""

# ---------- Public API ----------

def generate_csv_report(results: List[Dict]) -> str:
    """
    One row per discovered maize ortholog (across both discovery lanes), with
    enrichment columns populated where Phase 5/6/7 also ran.
    Queries with zero orthologs still get one row with empty maize fields, so
    the file is a faithful flat audit.
    """
    output = io.StringIO()

    # If any query in this report has CornCyc-derived data, prepend a hash-
    # comment attribution block. Pandas (`comment='#'`) and R's `read.csv`
    # (`comment.char='#'`) skip these lines automatically.
    has_corncyc = any((r.get("corncyc_annotation") or {}).get("n_pathways", 0) for r in results)
    if has_corncyc:
        output.write(
            "# This CSV contains data derived from CornCyc (Plant Metabolic Network, PMN).\n"
            "# CornCyc authors: Charles Hawkins, Bo Xue, Seung Yon Rhee. License + citation:\n"
            "# https://www.plantcyc.org/about/license-agreement  |  https://www.plantcyc.org/databases/corncyc/\n"
            "# Maysquery uses CornCyc as a read-only reference; see CORNCYC_ATTRIBUTION.txt\n"
            "# in the Maysquery source tree for the full list of modifications applied to\n"
            "# CornCyc-derived fields. For canonical pathway data follow the PMN links in the\n"
            "# corncyc_top_pathway / corncyc_gene_pathway_names columns.\n"
            "#\n"
        )

    writer = csv.writer(output)
    writer.writerow([
        # Query identification
        "query_index", "query_type", "query_input", "ion_mode",
        "chebi_id", "pubchem_cid", "monoisotopic_mass",
        # Per-query rollups
        "num_reactions", "num_proteins",
        "num_orthologs_total", "num_consensus", "num_sequence_only", "num_structure_only",
        "num_curated_only", "num_enriched",
        # CornCyc query-level summary (same for every row of this query)
        "corncyc_compound_ids", "corncyc_n_pathways", "corncyc_top_pathway", "corncyc_pathway_list",
        # Query-level reaction context (same on every row of this query — useful for
        # spreadsheet filtering by pathway, e.g. "show me all hits in the TCA cycle")
        "top_reaction_rhea_id", "top_reaction_equation", "top_reaction_ec",
        "top_reaction_type",                # "Transport" | "Metabolic" | ""
        "top_kegg_pathway_name", "top_kegg_pathway_id", "kegg_pathway_summary",
        # Maize gene + discovery lane
        "rank_overall", "rank_in_lane",
        "maize_gene_model", "maize_gene_url",
        "gene_symbol",        # e.g. "SDH1_0" (Gramene/Ensembl Plants display_name)
        "gene_description",   # e.g. "succinate dehydrogenase4"
        "gene_synonyms",      # pipe-separated union of v3+v4+community names
        "gene_v3_ids",        # legacy NAM v3 IDs (GRMZM2G*), pipe-separated
        "gene_v4_ids",        # legacy NAM v4 IDs (Zm00001d*), pipe-separated
        "gene_other_synonyms",# community names (e.g. 'sudh4'), pipe-separated
        # Phytozome (JGI BioMart) — independent KEGG-KO description + Panther family
        "phytozome_description", "phytozome_panther_ids", "phytozome_panther_descs",
        "phytozome_url",
        "consensus_class",                  # consensus | sequence_only | structure_only | curated_only
        "discovery_sources",                # pipe-separated short source names
        "num_sources", "consensus_score",
        "methods_summary",                  # human-readable, includes per-method evidence
        # Boolean (yes/no) columns per discovery method — one-click spreadsheet filter
        "found_by_ensembl_compara",
        "found_by_ensembl_panhomology",
        "found_by_plaza",
        "found_by_local_hmmer",
        "found_by_foldseek_structural",
        "found_by_corncyc_curated",
        # Per-method evidence detail (only populated for methods that found this gene)
        "ensembl_compara_pct_id", "ensembl_compara_ortholog_type",
        "panhomology_pct_id",     "panhomology_ortholog_type",
        "plaza_orthogroup_id",
        "hmmer_e_value", "hmmer_bit_score",
        "foldseek_qtm", "foldseek_ttm", "foldseek_prob", "foldseek_lddt", "foldseek_target_uniprot",
        "corncyc_n_reactions_for_gene",
        # Query pan-life protein
        "query_uniprot_id", "query_uniprot_url",
        # Sequence-lane evidence
        "sequence_similarity_pct",
        # Structure-lane evidence
        "structural_tm_score",
        # Enrichment (Phase 5)
        "enriched", "enrichment_kind", "plddt",
        "n_expression_experiments", "expression_experiment_ids",
        # Domain (Phase 6)
        "pfam_domain_id", "pfam_domain_name", "pfam_domain_url",
        "domain_residue_start", "domain_residue_end", "domain_tm_score",
        # Compara (Phase 7)
        "num_compara_orthologs",
        "top_compara_species", "top_compara_gene", "top_compara_gene_url", "top_compara_percent_identity",
        # CornCyc per-gene annotation (only populated for rows where CornCyc
        # lists this gene as catalysing a reaction involving the query compound)
        "corncyc_gene_pathway_ids", "corncyc_gene_pathway_names", "corncyc_gene_reactions",
        # PLAZA orthogroup hint / Foldseek score-breakdown tag
        "plaza_orthogroup",
        # Run status
        "status", "error",
    ])

    for q_idx, res in enumerate(results, start=1):
        inp = res.get("input_data") or {}
        qtype = inp.get("query_type", "mz")
        query_input = (
            inp.get("mz", "") if qtype == "mz"
            else inp.get("chemical_name", "") if qtype == "chemical"
            else inp.get("ec_number", "")
        )
        mode = inp.get("mode", "") if qtype == "mz" else ""
        ce = res.get("chemical_entity") or {}
        s = summarize_result(res)
        orthologs = res.get("orthologs") or []
        enriched_rows = join_enrichment(orthologs, res)

        # CornCyc query-level summary for the CSV
        cc_ann_q = res.get("corncyc_annotation") or {}
        cc_compound_ids = "|".join(c.get("id", "") for c in (cc_ann_q.get("compounds") or []))
        cc_n_pathways   = cc_ann_q.get("n_pathways", 0)
        cc_top_pathway  = ((cc_ann_q.get("pathways") or [{}])[0]).get("common_name", "") if cc_ann_q.get("pathways") else ""
        cc_pathway_list = "|".join((p.get("common_name") or p.get("id", "")) for p in (cc_ann_q.get("pathways") or []))

        # Reaction context for empty rows
        top_rxn = (res.get("reactions") or [{}])[0] if res.get("reactions") else {}
        rxn_rhea = top_rxn.get("rhea_id", "")
        rxn_eq   = top_rxn.get("equation") or top_rxn.get("label", "")
        rxn_ec   = "|".join(top_rxn.get("ec_numbers", []) or [])
        rxn_type = ("Transport" if top_rxn.get("is_transport") else "Metabolic") if top_rxn else ""
        pw_name  = (top_rxn.get("pathway_names") or [""])[0] if top_rxn else ""
        pw_id    = (top_rxn.get("pathway_ids") or [""])[0] if top_rxn else ""
        pw_summary = " | ".join((top_rxn.get("pathway_names") or [])[:5]) if top_rxn else ""

        # Empty-result audit row
        if not enriched_rows:
            writer.writerow([
                q_idx, qtype, query_input, mode,
                ce.get("chebi_id", ""), ce.get("pubchem_cid", ""), ce.get("monoisotopic_mass", ""),
                s["num_reactions"], s["num_proteins"],
                0, 0, 0, 0, 0, 0,                                          # +num_curated_only
                cc_compound_ids, cc_n_pathways, cc_top_pathway, cc_pathway_list,  # CornCyc query-level
                rxn_rhea, rxn_eq, rxn_ec, rxn_type, pw_name, pw_id, pw_summary,
                "", "", "", "", "", "", "", "", "", 0, 0,    # rank/maize/symbol/desc/synonyms/lane
                "", "", "",                                    # v3, v4, other synonyms
                "", "", "", "",                                # phytozome description/panther_ids/panther_descs/url
                "",                                            # methods_summary
                "no", "no", "no", "no", "no", "no",            # found_by_* booleans
                "", "", "", "", "",                            # ensembl, panhomology, plaza per-method evidence
                "", "",                                        # hmmer e_value + bit_score
                "", "", "", "", "",                            # foldseek qtm/ttm/prob/lddt/target_uniprot
                "",                                            # corncyc n_reactions
                "", "",
                "", "",
                "", "", "", "", "",       # enriched, enrichment_kind, plddt, n_expression_experiments, expression_experiment_ids
                "", "", "", "", "", "",
                0, "", "", "", "",
                "", "", "",               # corncyc_gene_pathway_ids/names/reactions
                "",
                s["status"], res.get("error", ""),
            ])
            continue

        # Per-lane rank counters
        lane_counter: Dict[str, int] = {}
        for overall_rank, r in enumerate(enriched_rows, start=1):
            cls = r["consensus_class"]
            lane_counter[cls] = lane_counter.get(cls, 0) + 1
            gene = r["gene"]
            quniprot = r["query_uniprot_id"]

            rxn_type_row = ("Transport" if r.get("top_reaction_is_transport") else "Metabolic") if r.get("top_reaction_rhea_id") else ""
            writer.writerow([
                q_idx, qtype, query_input, mode,
                ce.get("chebi_id", ""), ce.get("pubchem_cid", ""), ce.get("monoisotopic_mass", ""),
                s["num_reactions"], s["num_proteins"],
                s["num_orthologs"], s["num_consensus"], s["num_sequence_only"], s["num_structure_only"],
                s.get("num_curated_only", 0), s["num_targets"],
                cc_compound_ids, cc_n_pathways, cc_top_pathway, cc_pathway_list,
                r["top_reaction_rhea_id"], r["top_reaction_equation"], r["top_reaction_ec"], rxn_type_row,
                r["top_kegg_pathway_name"], r["top_kegg_pathway_id"], r["kegg_pathway_summary"],
                overall_rank, lane_counter[cls],
                gene, url_maize_gene(gene) if gene else "",
                r.get("gene_symbol", ""), r.get("gene_description", ""),
                "|".join(r.get("gene_synonyms", []) or []),
                "|".join(r.get("gene_v3_ids", []) or []),
                "|".join(r.get("gene_v4_ids", []) or []),
                "|".join(r.get("gene_other_synonyms", []) or []),
                r.get("phytozome_description", ""),
                "|".join(r.get("phytozome_panther_ids", []) or []),
                "|".join(r.get("phytozome_panther_descs", []) or []),
                f"https://phytozome-next.jgi.doe.gov/genePage/{gene}" if gene.startswith("Zm00001eb") else "",
                cls,
                "|".join(r["sources"]),
                r["num_sources"], r["consensus_score"],
                r.get("methods_summary", ""),
                # Boolean per-method columns
                *[("yes" if src_key in r["sources"] else "no") for src_key, _ in KNOWN_SOURCES],
                # Per-method evidence detail — pull straight from source_evidence dict
                (lambda ev: ev.get("pct_identity", "") if ev else "")(r.get("source_evidence", {}).get("Ensembl")),
                (lambda ev: ev.get("ortholog_type", "") if ev else "")(r.get("source_evidence", {}).get("Ensembl")),
                (lambda ev: ev.get("pct_identity", "") if ev else "")(r.get("source_evidence", {}).get("EnsemblPanHomology")),
                (lambda ev: ev.get("ortholog_type", "") if ev else "")(r.get("source_evidence", {}).get("EnsemblPanHomology")),
                (lambda ev: ev.get("orthogroup", "") if ev else "")(r.get("source_evidence", {}).get("PLAZA")),
                (lambda ev: ev.get("e_value", "") if ev else "")(r.get("source_evidence", {}).get("Local_HMMER")),
                (lambda ev: ev.get("bit_score", "") if ev else "")(r.get("source_evidence", {}).get("Local_HMMER")),
                (lambda ev: ev.get("qtm", "") if ev else "")(r.get("source_evidence", {}).get("Foldseek-structural")),
                (lambda ev: ev.get("ttm", "") if ev else "")(r.get("source_evidence", {}).get("Foldseek-structural")),
                (lambda ev: ev.get("prob", "") if ev else "")(r.get("source_evidence", {}).get("Foldseek-structural")),
                (lambda ev: ev.get("lddt", "") if ev else "")(r.get("source_evidence", {}).get("Foldseek-structural")),
                (lambda ev: ev.get("target_uniprot", "") if ev else "")(r.get("source_evidence", {}).get("Foldseek-structural")),
                (lambda ev: ev.get("n_reactions", "") if ev else "")(r.get("source_evidence", {}).get("CornCyc")),
                quniprot, url_uniprot(quniprot) if quniprot else "",
                ("" if r["sequence_similarity"] is None else round(r["sequence_similarity"], 2)),
                ("" if r["structural_tm_score"] is None else round(r["structural_tm_score"], 3)),
                "yes" if r["enriched"] else "no",
                r.get("enrichment_kind", ""),
                ("" if r["plddt"] is None else round(r["plddt"], 1)),
                r["n_expression_experiments"],
                "|".join(r["expression_experiments"]),
                r["pfam_domain_id"], r["pfam_domain_name"],
                url_pfam(r["pfam_domain_id"]) if r["pfam_domain_id"] else "",
                r["pfam_start"], r["pfam_end"],
                ("" if r["domain_tm_score"] is None else round(r["domain_tm_score"], 3)),
                r["num_compara_orthologs"],
                r["top_compara_species"], r["top_compara_gene"],
                url_ensembl_plants_gene(r["top_compara_species"], r["top_compara_gene"]) if r["top_compara_gene"] else "",
                ("" if r["top_compara_pct_id"] == "" else round(float(r["top_compara_pct_id"]), 1)),
                "|".join(r.get("corncyc_pathway_ids", []) or []),
                "|".join(r.get("corncyc_pathway_names", []) or []),
                "|".join(r.get("corncyc_reactions", []) or []),
                r["plaza_orthogroup"],
                s["status"], res.get("error", ""),
            ])

    return output.getvalue()


def generate_html_report(results: List[Dict]) -> str:
    summaries = [summarize_result(r) for r in results]
    totals = {
        "reactions": sum(s["num_reactions"] for s in summaries),
        "proteins": sum(s["num_proteins"] for s in summaries),
        "orthologs": sum(s["num_orthologs"] for s in summaries),
        "consensus": sum(s["num_consensus"] for s in summaries),
        "sequence_only": sum(s["num_sequence_only"] for s in summaries),
        "structure_only": sum(s["num_structure_only"] for s in summaries),
        "curated_only": sum(s["num_curated_only"] for s in summaries),
        "targets": sum(s["num_targets"] for s in summaries),
        "domain_targets": sum(s["num_domain_targets"] for s in summaries),
        "advanced_targets": sum(s["num_advanced_targets"] for s in summaries),
        "corncyc_pathways": sum(s.get("corncyc_pathways", 0) for s in summaries),
    }
    template = Template(HTML_TEMPLATE)
    return template.render(
        results=results,
        summaries=summaries,
        totals=totals,
        generated_at=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        url_maize_gene=url_maize_gene,
        url_ensembl_plants_gene=url_ensembl_plants_gene,
        url_uniprot=url_uniprot,
        url_chebi=url_chebi,
        url_pubchem=url_pubchem,
        url_rhea=url_rhea,
        url_pfam=url_pfam,
        url_kegg_pathway=url_kegg_pathway,
        url_quickgo=url_quickgo,
        # Lane-aware view helpers
        split_lanes=split_orthologs_by_lane,
        enriched_rows=join_enrichment,
        lane_label=lane_label,
        classify_ortholog=classify_ortholog,
        gene_label=gene_label_html,
    )


def save_report(content: str, ext: str) -> str:
    report_id = str(uuid.uuid4())
    filename = f"{report_id}.{ext}"
    filepath = os.path.join(REPORTS_DIR, filename)
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    return report_id
