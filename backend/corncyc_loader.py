"""
corncyc_loader.py — parse the CornCyc Pathway/Genome Database (PGDB)
flatfiles distributed by the Plant Metabolic Network at
https://plantcyc.org/database_imported/corncyc-12-0-0/ (or later versions).

The PGDB is license-restricted, so it isn't committed to the repo. Users
who have agreed to the PMN license download it themselves and drop the
extracted folder at one of these locations:

  <repo_root>/corncyc/<version>/data/
  <env $CORNCYC_DIR>/<version>/data/

The whole module is **opt-in**: every public function checks if the data
directory exists and returns an empty result if not, so the pipeline keeps
working when CornCyc isn't installed.

We parse five flatfiles into in-memory indices on first use (~150 MB on
disk → ~50 MB resident, ~5 s parse time on a typical workstation):

  compounds.dat   — compound frames with ChEBI/PubChem/KEGG xrefs
  reactions.dat   — reaction frames with LEFT/RIGHT compounds + pathway membership
  enzrxns.dat     — enzyme-reaction associations (links a protein to a reaction)
  proteins.dat    — protein frames with GENE link (Zm00001eb* v5 IDs)
  pathways.dat    — pathway frames with common names + reaction lists

The Pathway Tools attribute-value flatfile format is documented at
http://bioinformatics.ai.sri.com/ptools/flatfile-format.html — each record
is a sequence of `KEY - VALUE` lines terminated by `//`. Multi-valued
attributes repeat the key.
"""

from __future__ import annotations

import logging
import os
import re
import threading
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ----- Locating CornCyc -----

BACKEND_DIR = Path(__file__).resolve().parent
REPO_ROOT = BACKEND_DIR.parent
DEFAULT_CORNCYC_ROOT = Path(os.environ.get("CORNCYC_DIR", REPO_ROOT / "corncyc"))


def _is_corncyc_data_dir(d: Path) -> bool:
    """Heuristic: a `data/` dir that contains the core PGDB flatfiles."""
    if not d.is_dir():
        return False
    return all((d / fname).is_file() for fname in
               ("compounds.dat", "reactions.dat", "pathways.dat", "proteins.dat"))


def _detect_data_dir() -> Optional[Path]:
    """
    Locate a CornCyc `data/` directory. Tolerates several common layouts so
    users don't have to massage the PMN extraction into one canonical form:

      <root>/<version>/data/         ← native PMN tar layout (preferred)
      <root>/default-version + data  ← uses default-version file
      <root>/data/                   ← user already cd'd into a version dir
      <root>/corncyc-<ver>/data/     ← variant naming from some tarballs
      <root>/corncyc/<ver>/data/     ← double-nested ("corncyc/corncyc/...")

    Honours the version named in `<root>/default-version` when present;
    otherwise picks the lexically-latest versioned subdirectory.
    """
    root = DEFAULT_CORNCYC_ROOT
    if not root.is_dir():
        return None

    # 1. <root>/default-version → <root>/<that-version>/data
    default_file = root / "default-version"
    if default_file.is_file():
        version = default_file.read_text(encoding="utf-8").strip()
        candidate = root / version / "data"
        if _is_corncyc_data_dir(candidate):
            return candidate

    # 2. <root>/data/                 (user extracted into a version dir already)
    if _is_corncyc_data_dir(root / "data"):
        return root / "data"

    # 3. Any direct subdir with a populated data/
    matches = []
    for sub in root.iterdir():
        if not sub.is_dir():
            continue
        if _is_corncyc_data_dir(sub / "data"):
            matches.append(sub / "data")

    # 4. One more hop down (handles `corncyc/corncyc/<ver>/data/` and similar)
    if not matches:
        for sub in root.iterdir():
            if not sub.is_dir():
                continue
            for sub2 in sub.iterdir():
                if sub2.is_dir() and _is_corncyc_data_dir(sub2 / "data"):
                    matches.append(sub2 / "data")

    if not matches:
        return None
    # Pick the lexically-latest by the version-bearing parent dir name
    matches.sort(key=lambda p: p.parent.name, reverse=True)
    return matches[0]


def is_available() -> bool:
    """True iff the CornCyc data dir is on disk."""
    return _detect_data_dir() is not None


# ----- Pathway Tools flatfile iterator -----

_KEY_RE = re.compile(r"^([A-Z][A-Z0-9?\-]*)\s+-\s*(.*)$")


def _iter_records(path: Path) -> Iterable[Dict[str, List[str]]]:
    """
    Yield `{KEY: [value, value, …]}` dicts from a Pathway Tools .dat file.
    Records are terminated by `//`. Multi-valued attributes (which repeat
    on multiple lines) collapse into the same list. Comment / continuation
    lines that start with `/` are ignored (they're prose annotations).
    """
    current: Dict[str, List[str]] = {}
    with path.open("r", encoding="latin-1", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line == "//":
                if current:
                    yield current
                    current = {}
                continue
            if line.startswith("#") or line.startswith("/"):
                continue
            m = _KEY_RE.match(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2).strip()
            current.setdefault(key, []).append(val)
        if current:
            yield current


# ----- Loaded indices -----

class CornCycIndex:
    """In-memory cross-reference indices over a parsed CornCyc PGDB."""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.version = data_dir.parent.name

        # compound_id -> {chebi, pubchem, kegg, hmdb, name, synonyms, inchi_key}
        self.compounds: Dict[str, dict] = {}
        # chebi_numeric -> [compound_id, ...]
        self.chebi_to_compounds: Dict[str, List[str]] = {}
        # pubchem_id -> [compound_id, ...]
        self.pubchem_to_compounds: Dict[str, List[str]] = {}
        # inchi_key -> [compound_id, ...] (for ChEBI-miss fallback)
        self.inchikey_to_compounds: Dict[str, List[str]] = {}

        # reaction_id -> {ec, common_name, left, right, pathways, enzrxns}
        self.reactions: Dict[str, dict] = {}
        # compound_id -> [reaction_id, ...]
        self.compound_to_reactions: Dict[str, Set[str]] = {}

        # enzrxn_id -> protein_id
        self.enzrxn_to_protein: Dict[str, str] = {}

        # protein_id -> gene_id
        self.protein_to_gene: Dict[str, str] = {}

        # gene_id -> protein_ids (back-map, for convenience)
        self.gene_to_proteins: Dict[str, Set[str]] = {}

        # pathway_id -> {common_name, synonyms, reactions, types}
        self.pathways: Dict[str, dict] = {}
        # reaction_id -> [pathway_id, ...]
        self.reaction_to_pathways: Dict[str, Set[str]] = {}

    # ----- Parsing per file -----

    def _parse_compounds(self):
        path = self.data_dir / "compounds.dat"
        if not path.is_file():
            return
        for rec in _iter_records(path):
            cid = (rec.get("UNIQUE-ID") or [None])[0]
            if not cid:
                continue
            entry: dict = {
                "id": cid,
                "name": (rec.get("COMMON-NAME") or [""])[0],
                "synonyms": list(rec.get("SYNONYMS") or []),
                "inchi_key": ((rec.get("INCHI-KEY") or [""])[0] or "").replace("InChIKey=", ""),
                "chebi": "",
                "pubchem": "",
                "kegg": "",
                "hmdb": "",
                "metanetx": "",
            }
            for raw in rec.get("DBLINKS") or []:
                # Lines look like: (CHEBI "29806" NIL |taltman| 3452438104 NIL NIL)
                m = re.match(r'^\(([A-Z\-]+)\s+"([^"]+)"', raw)
                if not m:
                    continue
                db, val = m.group(1), m.group(2)
                if   db == "CHEBI":     entry["chebi"] = val
                elif db == "PUBCHEM":   entry["pubchem"] = val
                elif db == "LIGAND-CPD":entry["kegg"] = val
                elif db == "HMDB":      entry["hmdb"] = val
                elif db == "METANETX":  entry["metanetx"] = val
            self.compounds[cid] = entry
            if entry["chebi"]:
                self.chebi_to_compounds.setdefault(entry["chebi"], []).append(cid)
            if entry["pubchem"]:
                self.pubchem_to_compounds.setdefault(entry["pubchem"], []).append(cid)
            if entry["inchi_key"]:
                self.inchikey_to_compounds.setdefault(entry["inchi_key"], []).append(cid)

    def _parse_reactions(self):
        path = self.data_dir / "reactions.dat"
        if not path.is_file():
            return
        for rec in _iter_records(path):
            rid = (rec.get("UNIQUE-ID") or [None])[0]
            if not rid:
                continue
            left = list(rec.get("LEFT") or [])
            right = list(rec.get("RIGHT") or [])
            entry: dict = {
                "id": rid,
                "common_name": (rec.get("COMMON-NAME") or [""])[0],
                "ec": [v.replace("EC-", "") for v in (rec.get("EC-NUMBER") or [])],
                "left": left,
                "right": right,
                "in_pathway": list(rec.get("IN-PATHWAY") or []),
                "enzrxns": list(rec.get("ENZYMATIC-REACTION") or []),
            }
            self.reactions[rid] = entry
            for c in left + right:
                # LEFT/RIGHT often have coefficient prefixes like "2 NADH"; strip the leading number
                cid = c.split()[-1] if " " in c else c
                self.compound_to_reactions.setdefault(cid, set()).add(rid)
            for p in entry["in_pathway"]:
                self.reaction_to_pathways.setdefault(rid, set()).add(p)

    def _parse_enzrxns(self):
        path = self.data_dir / "enzrxns.dat"
        if not path.is_file():
            return
        for rec in _iter_records(path):
            erid = (rec.get("UNIQUE-ID") or [None])[0]
            if not erid:
                continue
            enzyme = (rec.get("ENZYME") or [""])[0]
            if enzyme:
                self.enzrxn_to_protein[erid] = enzyme

    def _parse_proteins(self):
        path = self.data_dir / "proteins.dat"
        if not path.is_file():
            return
        for rec in _iter_records(path):
            pid = (rec.get("UNIQUE-ID") or [None])[0]
            if not pid:
                continue
            gene = (rec.get("GENE") or [""])[0]
            if gene:
                # The .dat file is uppercase; v5 maize models start with "ZM00001EB".
                # Normalise to the standard lowercase casing the rest of the pipeline uses.
                norm = gene
                if norm.upper().startswith("ZM00001EB"):
                    norm = "Zm" + norm[2:].lower()
                self.protein_to_gene[pid] = norm
                self.gene_to_proteins.setdefault(norm, set()).add(pid)

    def _parse_pathways(self):
        path = self.data_dir / "pathways.dat"
        if not path.is_file():
            return
        for rec in _iter_records(path):
            pid = (rec.get("UNIQUE-ID") or [None])[0]
            if not pid:
                continue
            common = (rec.get("COMMON-NAME") or [""])[0]
            syns = list(rec.get("SYNONYMS") or [])
            reactions = list(rec.get("REACTION-LIST") or [])
            types = list(rec.get("TYPES") or [])
            # Some pathways list reactions only in REACTION-LAYOUT — pull from there too
            for layout in rec.get("REACTION-LAYOUT") or []:
                m = re.match(r"^\(([^\s\)]+)", layout)
                if m and m.group(1) not in reactions:
                    reactions.append(m.group(1))
            self.pathways[pid] = {
                "id": pid,
                "common_name": common,
                "synonyms": syns,
                "reactions": reactions,
                "types": types,
            }

    def load(self):
        logger.info(f"Loading CornCyc {self.version} from {self.data_dir}…")
        self._parse_compounds()
        self._parse_reactions()
        self._parse_enzrxns()
        self._parse_proteins()
        self._parse_pathways()
        logger.info(
            f"CornCyc loaded: {len(self.compounds)} compounds "
            f"({len(self.chebi_to_compounds)} ChEBI-mapped), "
            f"{len(self.reactions)} reactions, "
            f"{len(self.pathways)} pathways, "
            f"{len(self.protein_to_gene)} proteins, "
            f"{len(self.gene_to_proteins)} unique maize genes."
        )

    # ----- Public lookup -----

    def compounds_for_chebi(self, chebi_id: str) -> List[str]:
        """ChEBI ID (with or without 'CHEBI:' prefix) → CornCyc compound frame IDs."""
        num = chebi_id.split(":")[-1] if chebi_id else ""
        return list(self.chebi_to_compounds.get(num, []))

    def reactions_for_compound(self, compound_id: str) -> List[str]:
        return sorted(self.compound_to_reactions.get(compound_id, ()))

    def pathways_for_reaction(self, reaction_id: str) -> List[str]:
        return sorted(self.reaction_to_pathways.get(reaction_id, ()))

    def genes_for_reaction(self, reaction_id: str) -> List[str]:
        """Maize genes whose proteins catalyse this reaction (via ENZYMATIC-REACTION → ENZYME → GENE)."""
        rxn = self.reactions.get(reaction_id) or {}
        genes: Set[str] = set()
        for erid in rxn.get("enzrxns") or ():
            pid = self.enzrxn_to_protein.get(erid)
            if not pid:
                continue
            g = self.protein_to_gene.get(pid)
            if g:
                genes.add(g)
        return sorted(genes)


# ----- Module-level singleton (lazy, thread-safe load) -----

_INDEX: Optional[CornCycIndex] = None
_LOAD_LOCK = threading.Lock()


def get_index(force_reload: bool = False) -> Optional[CornCycIndex]:
    """Return the loaded CornCyc index, or None if CornCyc isn't installed."""
    global _INDEX
    if _INDEX is not None and not force_reload:
        return _INDEX
    data_dir = _detect_data_dir()
    if data_dir is None:
        return None
    with _LOAD_LOCK:
        if _INDEX is None or force_reload:
            idx = CornCycIndex(data_dir)
            idx.load()
            _INDEX = idx
    return _INDEX


def get_status() -> dict:
    """Snapshot for `/api/corncyc/status` and the UI banner."""
    data_dir = _detect_data_dir()
    if data_dir is None:
        return {
            "available": False,
            "loaded": False,
            "expected_path": str(DEFAULT_CORNCYC_ROOT),
        }
    idx = _INDEX
    return {
        "available": True,
        "loaded": idx is not None,
        "version": data_dir.parent.name,
        "data_dir": str(data_dir),
        "compounds": len(idx.compounds) if idx else None,
        "chebi_mapped_compounds": len(idx.chebi_to_compounds) if idx else None,
        "reactions": len(idx.reactions) if idx else None,
        "pathways": len(idx.pathways) if idx else None,
        "maize_genes": len(idx.gene_to_proteins) if idx else None,
    }
