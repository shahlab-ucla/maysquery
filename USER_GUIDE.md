# Maysquery — User Guide

This document covers the eight pipeline phases, the data sources each one
queries, the configuration knobs, and the shape of the report output. For
install instructions see [INSTALLATION.md](INSTALLATION.md).

---

## 1. Input modes

The Single Query tab offers three entry points:

| Mode | Bypasses | Starts at |
|---|---|---|
| **m/z value** | none | Phase 1 — CMM mass lookup |
| **Chemical name** | Phase 1 | PubChem name → CID → ChEBI |
| **EC number** | Phases 1 + 2 | UniProt direct EC search |

A **Batch Processing** tab accepts a CSV (template downloadable from the
tab) and runs each row sequentially, emitting an aggregate HTML + CSV
report at the end.

---

## 2. The eight phases

### Phase 1 — Chemical entity identification

Resolves an input m/z into a ChEBI identifier plus the conjugate-acid/base
family it belongs to.

- **CEU Mass Mediator** (`/api/v3/batch`) is the primary lookup. Requires
  `metabolites_type`, `masses_mode`, `deuterium` fields or the server
  throws a NullPointerException — those are populated by default.
- **Metabolomics Workbench** `/rest/moverz/MB/{mz}/{adduct}/{tol}/json`
  is the fallback (note: returns tab-separated text despite the `/json`
  URL suffix; the parser handles this).
- **PubChem PUG-REST** resolves a compound name to a CID → SMILES →
  monoisotopic mass.
- **PubChem → ChEBI cross-reference** picks up the canonical ChEBI ID
  (typically the neutral form).
- **ChEBI conjugate-form expansion** (via OLS4) walks the
  `is_protonated_form_of` / `is_deprotonated_form_of` / `is_tautomer_of`
  graph up to 2 hops, so a query that lands on `CHEBI:18012` (fumaric acid)
  automatically also looks up `CHEBI:29806` (fumarate²⁻) and
  `CHEBI:37154` (fumarate¹⁻) — without this expansion, Rhea returns 0 hits
  because it uses the ionised form.

Tunable: `cmm_tolerance_ppm`, `chebi_expansion_depth`.

### Phase 2 — Reaction networks

For the ChEBI family from Phase 1:

- **Rhea SPARQL endpoint** with a `VALUES ?target_chebi { … }` clause
  covering every form. The query selects reaction ID, equation text, label,
  EC number(s), `isTransport`, `isBalanced`. A separate `COUNT(DISTINCT)`
  query runs when the row cap is hit, so the log can report
  "showing top 100 of 247 available".
- **KEGG** `/conv/compound/chebi:{id}` → KEGG compound, then
  `/link/pathway/{cpd_id}` → pathway IDs, then a batched `/list/{ids}` →
  human-readable pathway names (e.g. "Citrate cycle (TCA cycle)" rather
  than just `path:map00020`).

Tunable: `rhea_fetch_limit` (default 100).

### Phase 3 — Pan-life protein pool

For each reaction (or directly for the EC-number entry point):

- **UniProt REST** queries split into three categories:
  - Enzymes — either from Rhea IDs (`rhea:"RHEA:XXX" OR rhea:"RHEA:YYY"`),
    or `chebi:"…" AND keyword:KW-0255` if Phase 2 returned nothing
  - Transporters — `chebi:"…" AND keyword:KW-0813`
  - Receptors — `chebi:"…" AND keyword:KW-0675`
- All three queries OR together the full ChEBI conjugate family.
- The total-matching count from UniProt's `X-Total-Results` header is
  surfaced in the log, so you can see "top 25 of 1028 matching — bump
  UNIPROT_SIZE_PER_CATEGORY to see more".
- Results are deduplicated by UniProt accession.

Tunable: `uniprot_size_per_category` (default 25).

### Phase 4 — Sequence orthology (lane 1 of discovery)

For each pan-life protein, looks up *Zea mays* orthologs by sequence:

- **Ensembl Compara plants** — UniProt → Gramene/TAIR/EnsemblPlants gene
  ID via UniProt REST → `rest.ensembl.org/homology/id/{species}/{gene_id}?compara=plants;target_species=zea_mays`.
  Returns % identity.
- **Ensembl pan-homology** — same lookup but with `compara=pan_homology` for
  cross-kingdom orthology (e.g. cyanobacterial chloroplast lineage genes
  with maize counterparts).
- **PLAZA Monocots v5** — best-effort; the anonymous endpoint started
  returning 403 in 2026 so this typically logs a warning and returns nothing.
- **HMMER `phmmer` fallback** — only fires when all three live APIs return
  zero. Runs locally against the pre-downloaded Zm-B73-REFERENCE-NAM-5.0
  proteome (~12 MB).

All hits land in the **consensus reducer**, which dedupes by maize gene
model and sums sources. Hits get classified into discovery lanes (see
Phase 4.5 below).

Tunable: `hmmer_e_value` (default 1e-5).

### Phase 4.5 — Structural discovery (lane 2 of discovery)

Runs in parallel with Phase 4. Requires the maize AlphaFold + Foldseek
index built by `install_maize_afdb.py`; degrades gracefully if the index
isn't present (just logs a warning and skips).

- For each pan-life protein, downloads its AlphaFold PDB
  (`https://alphafold.ebi.ac.uk/api/prediction/{acc}` → current `pdbUrl`,
  forward-compatible across AF v4 / v5 / v6+).
- Runs `foldseek easy-search --alignment-type 2` against the indexed
  maize AlphaFold proteome (~40k structures, ~1–2 GB on disk).
- Returns full-protein TM scores: `qtmscore` (normalised by query length)
  and `ttmscore` (normalised by target length). The filter is
  `max(qTM, tTM) >= 0.5`, which catches three useful orthology patterns:
  true full-length ortholog (both high), query is a sub-domain of a
  larger maize protein (qTM high, tTM low), or vice versa.
- Hits are tagged with source `"Foldseek-structural"` and merge into the
  same consensus reducer as Phase 4.

Tunable: `foldseek_tm_threshold` (default 0.5), `foldseek_max_hits_per_query`
(default 10), `foldseek_concurrency` (default 2 parallel processes).

### Discovery lanes

After Phases 4 + 4.5 merge, each maize candidate is classified:

- **Consensus** — found by both a sequence DB *and* Foldseek; highest
  confidence, ranks first in the dashboard.
- **Structure-only** — Foldseek alone (the "hidden orthologs" — divergent
  in sequence but structurally conserved).
- **Sequence-only** — sequence DBs alone.

Each lane gets its own collapsible section in the UI with a one-line
summary per row; expand a row to see sources, query protein, Pfam domain
(if validated), top Compara hit, and tissue expression.

### Phase 5 — Enrichment

For every discovered maize ortholog (not just the top-N), Phase 5 adds:

- **AlphaFold pLDDT** (model quality, 0–100)
- **Gramene expression breadth** — count + IDs of EBI Expression Atlas
  experiments where the gene is reported as expressed. Per-tissue FPKM
  values are no longer available from the EBI REST API (EBI deprecated the
  `/gxa/json/search/baseline` endpoint in 2024); the qualitative breadth
  is what's exposed instead. Each experiment ID links out to its GXA page.
- **TM-score** — for consensus / structure-only hits, reused from Phase 4.5
  (no duplicated Foldseek call). For sequence-only hits in the top-N, a
  fresh 1-to-1 Foldseek alignment runs. Beyond the top-N, sequence-only
  hits leave TM blank (cheap-path enrichment).

A pLDDT > 70 filter trims targets whose AlphaFold model is too uncertain.
Targets with no model but ≥1 expression experiment are retained.

Tunable: `enrichment_top_n` (default 10), `plddt_threshold` (default 70).

### Phase 6 — Domain-level structural validation

For the top-N enriched targets:

- **InterPro Pfam API** → longest Pfam domain on the query protein.
- The query AlphaFold PDB is **sliced** to those residue coordinates.
- `foldseek easy-search` aligns the sliced query against the full maize
  target PDB → domain-level TM-score.

This catches cases where the *catalytic domain* is conserved but the
flanking regions diverge — a common pattern for plant secondary-metabolism
enzymes.

### CornCyc curated lane (opt-in, joins Phase 4)

When the user-licensed CornCyc PGDB is installed locally (see INSTALLATION.md),
a fourth discovery lane joins the consensus reducer:

- The query ChEBI is matched against CornCyc's compound frames via the
  `(CHEBI "…")` cross-references in `compounds.dat`. Conjugate-form
  expansion from Phase 1 ensures we still hit, e.g. fumarate(2-)
  (CHEBI:29806 → CornCyc `FUM`) even when PubChem gave us the neutral
  form (CHEBI:18012).
- For each matched compound, we walk `compounds.dat` → reactions touching
  the compound → `enzrxns.dat` (enzymatic-reaction frames) → `proteins.dat`
  (Polypeptide frames) → maize gene IDs (`GENE - ZM00001EB*` v5).
- Each annotated gene becomes an `OrthologMapping` with
  `sources=["CornCyc"]`, `similarity_score=100` (curation, not homology),
  and the `plaza_orthogroup` field encoding the pathway(s) it appears in.
- The consensus reducer then folds these into the same pool as Ensembl and
  Foldseek hits. **A gene found by ≥2 lanes (sequence, structure, or
  CornCyc) is now classified `consensus`** — the strongest evidence class.
- `curated_only` is a new lane class for genes annotated by CornCyc but
  not surfaced by any homology search.

In parallel, the response carries a `corncyc_annotation` block listing every
PlantCyc pathway the query compound participates in, with per-pathway
maize-gene lists. This drives the **"CornCyc maize pathway context"**
section in the dashboard and the green-tinted block in the HTML report.

CornCyc is opt-in — the pipeline runs fine without it and logs an info
message saying it's not installed.

### Phase 7 — Pan-plant Compara cross-check

For the top-N enriched targets:

- **Ensembl Compara plants** — `/homology/id/zea_mays/{gene}?compara=plants`
  returns the gene's orthologs across all plant species in Compara
  (Arabidopsis, rice, sorghum, brachypodium, …).
- Each ortholog gets a clickable Ensembl Plants gene-summary link.

Tunable: `compara_max_orthologs_per_target` (default 10).

### Final step — maize gene metadata enrichment

After all phases complete, a single batched call to Gramene's search API
resolves every unique maize gene ID surfaced anywhere in the response to a
human-readable label:

- `symbol` — e.g. `SDH1_0`
- `description` — e.g. `succinate dehydrogenase4`
- `synonyms` — legacy IDs + maize-community-style symbols (e.g. `sudh4`,
  `GRMZM2G079888`, `Zm00001d007966`)

These flow into the UI, CSV, and HTML report — so every `Zm00001eb*` ID is
shown alongside its readable label and synonyms.

---

## 3. Configuration tab

Every cutoff, cap, and threshold is exposed in the **Configuration** tab
with inline guidance:

| Field | Default | Phase | What it controls |
|---|---|---|---|
| `cmm_tolerance_ppm` | 5.0 | 1 | m/z mass tolerance for CMM + Metabolomics Workbench |
| `rhea_fetch_limit` | 100 | 2 | Rhea SPARQL LIMIT clause |
| `chebi_expansion_depth` | 2 | 1, 2, 3 | Hops through the ChEBI conjugate-form graph |
| `uniprot_size_per_category` | 25 | 3 | Top-N UniProt hits per Enzyme/Transporter/Receptor query |
| `hmmer_e_value` | 1e-5 | 4 | Fallback HMMER E-value cutoff |
| `foldseek_tm_threshold` | 0.5 | 4.5 | Min `max(qTM, tTM)` to accept a structural ortholog |
| `foldseek_max_hits_per_query` | 10 | 4.5 | Top-N maize structural hits per query protein |
| `foldseek_concurrency` | 2 | 4.5 | Parallel Foldseek processes (multi-threaded each) |
| `enrichment_top_n` | 10 | 5 | Top-N hits that get full enrichment (1-to-1 Foldseek for sequence-only) |
| `plddt_threshold` | 70.0 | 5 | Min AlphaFold pLDDT to keep an enriched target |
| `compara_max_orthologs_per_target` | 10 | 7 | Top-N pan-plant orthologs per maize target |

Changes auto-save in `localStorage` (key `maysquery.pipeline_config`) and
apply to every subsequent run. The Configuration tab also has per-field
reset links and a "Reset all to defaults" button.

The terminal log line at the start of each run summarises any non-default
overrides: `cfg-overrides=3 (foldseek_tm_threshold, plddt_threshold, …)`.

---

## 4. Output

### Live UI

The dashboard renders as:

1. **Executive summary** — high-contrast panel with resolved entity, count
   rollups, top KEGG pathway, lane breakdown, and a ⭐ top-target chip with
   gene + symbol + description + key metrics.
2. **Collapsible drill-down sections**:
   - Resolved chemical entity (open by default)
   - Reactions (collapsed; equation + EC + pathway names visible on expand)
   - Pan-life proteins (collapsed; each protein is itself collapsible to GO terms)
   - **Maize candidates** (open by default) with three lane sub-sections
   - Pan-plant Compara details (collapsed)
   - Execution log (collapsed)
   - Additional Homologs (Pending Enrichment) (collapsed) — for on-demand
     1-to-1 Foldseek on hits that didn't survive the Phase 5 filter
   - Pathway network graph (collapsed) — Cytoscape view, capped to top
     8 reactions × 10 proteins × 8 orthologs to stay readable
3. **Live terminal log** with color-coded per-phase events.

### HTML report

`/api/report/single` and `/api/batch` generate a self-contained HTML report
with the same hierarchy: per-query executive card up top, then `<details>`
sections for each phase. The Maize Candidates table includes columns for
Lane, Gene (linked to MaizeGDB), Query Protein, TM, Seq, pLDDT,
Expression (count + first 4 GXA experiment IDs linked), Pfam Domain (with
residue range + TM), and Top Compara.

### CSV report

One row per discovered maize ortholog (55 columns). Key fields:

- **Query context** — `query_index`, `query_type`, `query_input`, `ion_mode`,
  `chebi_id`, `pubchem_cid`, `monoisotopic_mass`, `num_reactions`,
  `num_proteins`, `num_orthologs_total`, lane breakdown counts,
  `top_reaction_equation`, `top_reaction_ec`, `kegg_pathway_summary`
- **Discovery** — `rank_overall`, `rank_in_lane`, `maize_gene_model`,
  `maize_gene_url`, `gene_symbol`, `gene_description`, `gene_synonyms`,
  `consensus_class` (consensus / sequence_only / structure_only),
  `discovery_sources` (pipe-separated)
- **Evidence** — `sequence_similarity_pct`, `structural_tm_score`,
  `query_uniprot_id` + URL
- **Enrichment** — `enriched`, `enrichment_kind` (full / cheap),
  `plddt`, `n_expression_experiments`, `expression_experiment_ids`
- **Domain (Phase 6)** — `pfam_domain_id`, `pfam_domain_name`,
  `pfam_domain_url`, residue range, `domain_tm_score`
- **Compara (Phase 7)** — `num_compara_orthologs`, top species/gene/URL/%ID

A spreadsheet filter like `consensus_class = "consensus" AND
gene_description CONTAINS "succinate"` is a one-click operation.

---

## 5. Failure modes & gotchas

- **PLAZA 403** — the anonymous PLAZA Monocots 5 endpoint started rejecting
  requests in 2026. Logged once per session; the pipeline keeps moving.
- **EBI Expression Atlas per-tissue values** — that REST path was
  deprecated. We now use Gramene's `expressed_in_gxa_attr_ss` for breadth
  (count + experiment IDs).
- **MaizeGDB / qTeller direct API** — Cloudflare blocks anonymous JSON
  requests. We get gene metadata through Gramene instead.
- **AlphaFold version churn** — the public PDBs migrated to v6. We hit the
  AlphaFold prediction API and use whatever `pdbUrl` it returns, so this
  is forward-compatible.
- **Multi-domain enzymes** — Phase 4.5 uses whole-chain TM-align, so a
  query enzyme whose single catalytic domain is fused to an unrelated
  N-terminal regulatory chunk may under-score against maize. The log line
  in Phase 4.5 explicitly calls this out; Phase 6 handles domain-level
  re-search but only on already-enriched hits.
- **Build cost of the structural index** — ~5 GB download + ~20–40 min
  one-time build. The pipeline still runs without it; Phase 4.5 just stays
  empty and the UI banner offers to build it.

---

## 6. Where things live in the repo

```
backend/
  main.py                  FastAPI app + /api/run_pipeline/stream SSE endpoint
  models.py                Pydantic models (PipelineConfig, OrthologMapping, …)
  phase1.py                CMM + MW + PubChem + ChEBI mapping
  phase2.py                Rhea SPARQL + KEGG pathway names
  phase3.py                UniProt enzyme/transporter/receptor search
  phase4.py                Sequence orthology (Ensembl Compara variants)
  phase4_5.py              Foldseek structural discovery vs maize AFDB
  phase5.py                AlphaFold pLDDT + Gramene expression + 1-to-1 Foldseek
  phase6.py                Pfam-sliced domain Foldseek
  phase7.py                Pan-plant Ensembl Compara cross-check
  chebi_utils.py           OLS-based ChEBI conjugate-form expansion
  hmmer_runner.py          Local phmmer subprocess wrapper (WSL on Windows)
  install_foldseek.py      Foldseek binary installer (CLI)
  install_maize_afdb.py    Maize AlphaFold proteome download + foldseek createdb
  maize_gene_meta.py       Gramene-sourced gene symbol/description lookup
  report_generator.py      HTML + CSV report builders
  static/
    index.html / app.js / style.css   Single-page dashboard

setup.ps1 / setup.sh       One-shot installer (Windows / Unix)
run.ps1 / run.sh           Launcher (Windows / Unix)
INSTALLATION.md            Cross-platform install guide
USER_GUIDE.md              This file
```
