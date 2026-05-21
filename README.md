# Maysquery

Maysquery is an automated bioinformatics pipeline that maps untargeted
metabolomics data (m/z values, chemical names, or EC numbers) to candidate
enzyme-coding genes in ***Zea mays*** (maize). It bridges the gap between
raw mass-spectrometry features and plant genetics through two parallel
discovery lanes — sequence homology and 3D structural homology — and
surfaces the highest-confidence hits as a one-page executive summary with
drill-down details.

## What you get per query

- **Resolved chemical entity** — m/z → CMM/Metabolomics Workbench → PubChem
  → ChEBI (with conjugate-acid/base expansion so reactions involving the
  ionised form aren't silently missed).
- **Reaction networks** — Rhea SPARQL with full reaction equations and EC
  numbers, plus KEGG pathway names (not just IDs).
- **Pan-life protein pool** — UniProt Swiss-Prot grouped by Enzyme /
  Transporter / Receptor, with per-category truncation visible in the log.
- **Maize candidates from two parallel discovery lanes:**
  - **Sequence lane**: Ensembl Compara plants + Ensembl pan-homology
    (cross-kingdom) + local HMMER fallback
  - **Structure lane**: Foldseek `--alignment-type 2` (TM-align mode)
    against a locally-indexed maize AlphaFold proteome
  - **Consensus hits** (found by both lanes) are flagged separately
- **Per-target enrichment** — AlphaFold pLDDT, Pfam-sliced domain TM,
  Gramene expression breadth (number of Expression Atlas experiments where
  the gene is detected), and Ensembl Compara cross-checks across pan-plant
  species.
- **Human-readable labels** everywhere — gene symbol + description from
  Gramene (e.g. `Zm00001eb117970 · SDH1_0 — succinate dehydrogenase4`)
  next to every maize gene ID, linked to MaizeGDB.

## Architecture at a glance

```
m/z  ─▶  Phase 1   CMM + PubChem + ChEBI       ─▶  ChEBI ID (+ conjugate forms)
         Phase 2   Rhea SPARQL + KEGG          ─▶  reactions + pathway names
         Phase 3   UniProtKB                   ─▶  pan-life Enzymes/Transporters/Receptors
         Phase 4   Ensembl + PLAZA  ┐
                                    ├──▶  consensus reducer  ─▶  maize candidates
         Phase 4.5 Foldseek vs AFDB ┘
         Phase 5   AlphaFold pLDDT + Gramene expression breadth + 1-to-1 Foldseek (TM)
         Phase 6   InterPro/Pfam domain-sliced Foldseek
         Phase 7   Ensembl Compara plants (pan-plant ortholog cross-check)
```

## UI overview

The dashboard leads with a high-contrast **executive summary** —
resolved chemical, count rollups, top KEGG pathway, top maize target with
TM/seq/pLDDT/expression metrics — followed by **collapsible drill-down
sections** for reactions, pan-life proteins, the three discovery lanes
(consensus / structure-only / sequence-only), pan-plant Compara details,
and execution logs. A **Configuration tab** exposes every tunable cutoff
(mass tolerance, Rhea/UniProt fetch caps, TM threshold, pLDDT filter,
enrichment top-N, Foldseek concurrency, …) with inline guidance.

Reports come out as both HTML (with collapsible sections mirroring the UI)
and CSV (one row per discovered maize candidate, 55 columns including gene
symbol, description, expression-experiment list, lane classification, and
all relevant URLs).

## Quick start

See **[INSTALLATION.md](INSTALLATION.md)** for the full cross-platform
guide. Short version:

**Windows** (needs WSL — installer walks you through it):
```powershell
git clone https://github.com/shahlab-ucla/maysquery.git
cd maysquery
.\setup.ps1
.\run.ps1
```

**macOS / Linux** (Foldseek + HMMER install natively):
```bash
git clone https://github.com/shahlab-ucla/maysquery.git
cd maysquery
./setup.sh
./run.sh
```

Then open <http://127.0.0.1:8008/static/index.html>.

## Phase 4.5 — structural-discovery index

Phase 4.5 is *opt-in*. It needs a one-time download + Foldseek-index build of
the *Zea mays* AlphaFold proteome (~5 GB download, ~1–2 GB persistent
index, ~20–40 min). The app starts without it and Phase 4.5 simply gets
skipped; a banner in the UI prompts you to build it when you're ready.

## Documentation

- **[INSTALLATION.md](INSTALLATION.md)** — install guide for Windows, macOS, Linux
- **[USER_GUIDE.md](USER_GUIDE.md)** — per-phase reference, configuration knobs, output schema
- **[Metabolite-to-Gene Pipeline Specification.md](Metabolite-to-Gene%20Pipeline%20Specification.md)** — original design spec
