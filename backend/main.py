import sys
import phase3
from fastapi import FastAPI, UploadFile, File
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, PlainTextResponse
import os
import csv
import io
from contextlib import asynccontextmanager
from models import MetaboliteInput, ValidatedTarget, OrthologMapping, PipelineConfig
from typing import List
from report_generator import generate_csv_report, generate_html_report, save_report, REPORTS_DIR

@asynccontextmanager
async def lifespan(app: FastAPI):
    yield

app = FastAPI(title="Metabolite-to-Gene Pipeline API", lifespan=lifespan)

# Setup static files directory
STATIC_DIR = os.path.join(os.path.dirname(__file__), "static")
os.makedirs(STATIC_DIR, exist_ok=True)

app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

@app.get("/")
def read_root():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))

from phase1 import execute_phase1, fetch_chemical_entity_by_name
from phase2 import execute_phase2
from phase3 import execute_phase3, execute_phase3_by_ec
from phase4 import execute_phase4
from phase5 import execute_phase5
from phase6 import execute_phase6
from phase7 import execute_phase7
from models import ExecutionLogEntry

import asyncio
import json
from fastapi.responses import StreamingResponse

class StreamingLogList:
    def __init__(self, queue: asyncio.Queue):
        self.queue = queue
        self.items = []

    def append(self, entry: ExecutionLogEntry):
        self.items.append(entry)
        self.queue.put_nowait({"type": "log", "entry": entry.model_dump()})
        
    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

async def _run_pipeline_core(input_data: MetaboliteInput, execution_logs) -> dict:
    chemical_entity = None
    reactions = []
    proteins = []
    orthologs = []
    targets = []
    unvalidated_targets = []
    domain_targets = []
    advanced_homology_targets = []

    cfg = input_data.pipeline_config or PipelineConfig()
    # Mirror config into the top-level legacy fields so existing phase code paths
    # (and tests that read input_data.hmmer_e_value / .tolerance_ppm directly)
    # see the same values the user set in the Configuration tab.
    input_data.hmmer_e_value = cfg.hmmer_e_value
    input_data.tolerance_ppm = cfg.cmm_tolerance_ppm

    top_n = cfg.enrichment_top_n

    if input_data.query_type == "mz":
        chemical_entity = await execute_phase1(input_data, logs=execution_logs)
        reactions = await execute_phase2(chemical_entity, logs=execution_logs, config=cfg)
        proteins = await execute_phase3(chemical_entity, reactions, logs=execution_logs, config=cfg)
        orthologs = await execute_phase4(
            proteins, input_data, logs=execution_logs, config=cfg,
            chebi_id=(chemical_entity.chebi_id if chemical_entity else None),
        )
    elif input_data.query_type == "chemical":
        chemical_entity = await fetch_chemical_entity_by_name(input_data.chemical_name, logs=execution_logs)
        reactions = await execute_phase2(chemical_entity, logs=execution_logs, config=cfg)
        proteins = await execute_phase3(chemical_entity, reactions, logs=execution_logs, config=cfg)
        orthologs = await execute_phase4(
            proteins, input_data, logs=execution_logs, config=cfg,
            chebi_id=(chemical_entity.chebi_id if chemical_entity else None),
        )
    elif input_data.query_type == "ec":
        proteins = await execute_phase3_by_ec(input_data.ec_number, logs=execution_logs, config=cfg)
        orthologs = await execute_phase4(
            proteins, input_data, logs=execution_logs, config=cfg,
            chebi_id=(chemical_entity.chebi_id if chemical_entity else None),
        )

    # Phase 5 now handles ALL orthologs internally:
    #  - Top `enrichment_top_n` sequence-only hits get the expensive 1-to-1 Foldseek
    #  - Consensus / structure-only hits reuse the Phase 4.5 TM (no duplicated work)
    #  - Everything else still gets pLDDT + Gramene expression breadth (cheap)
    targets = await execute_phase5(orthologs, logs=execution_logs, config=cfg)
    # "Unvalidated" now = orthologs that didn't survive the Phase 5 filter
    enriched_gene_ids = {t.maize_gene_model for t in targets}
    unvalidated_targets = [o for o in orthologs if o.maize_gene_model not in enriched_gene_ids]

    # Phase 6 (Pfam-sliced domain Foldseek) and Phase 7 (Compara) are expensive
    # per-target — keep them on the top-N enriched targets only.
    primary_targets = [t for t in targets if t.enrichment_kind == "full"][:top_n]
    primary_mappings = [o for o in orthologs[:top_n] if o.maize_gene_model in {t.maize_gene_model for t in primary_targets}]
    domain_targets = await execute_phase6(primary_targets, primary_mappings, logs=execution_logs)
    advanced_homology_targets = await execute_phase7(primary_targets, logs=execution_logs, config=cfg)

    # Batched lookup of human-readable gene metadata (symbol + description) for every
    # maize gene ID surfaced anywhere in the pipeline. Powers the
    # "Zm00001eb117970 · sdh4 — succinate dehydrogenase4" labels in the UI/reports.
    from maize_gene_meta import fetch_maize_gene_meta_batch
    from phytozome_lookup import fetch_phytozome_meta_batch
    unique_gene_ids = set()
    for o in orthologs:                  unique_gene_ids.add(o.maize_gene_model)
    for t in targets:                    unique_gene_ids.add(t.maize_gene_model)
    for d in domain_targets:             unique_gene_ids.add(d.maize_gene_model)
    for a in advanced_homology_targets:  unique_gene_ids.add(a.maize_gene_model)
    if unique_gene_ids:
        execution_logs.append(ExecutionLogEntry(
            phase=4, database="Gramene gene metadata", status="info", hits=0,
            message=f"Resolving symbol + description for {len(unique_gene_ids)} unique maize gene IDs",
        ))
    maize_gene_metadata = await fetch_maize_gene_meta_batch(unique_gene_ids)
    if unique_gene_ids:
        n_named = sum(1 for v in maize_gene_metadata.values() if v.get("symbol") or v.get("description"))
        sample = next(iter(maize_gene_metadata.values()), None)
        sample_label = ""
        if sample:
            sym, desc = sample.get("symbol", ""), sample.get("description", "")
            sample_label = f" (e.g. '{sym} — {desc}')" if (sym and desc) else ""
        execution_logs.append(ExecutionLogEntry(
            phase=4, database="Gramene gene metadata",
            status="success" if n_named else "warning",
            hits=n_named,
            message=f"{n_named}/{len(unique_gene_ids)} maize genes resolved{sample_label}",
        ))

    # Phytozome (JGI BioMart) — independent gene-family + KEGG-KO descriptions.
    # One batched POST per pipeline run; degrades silently if Phytozome is down.
    phytozome_metadata = await fetch_phytozome_meta_batch(unique_gene_ids)
    if phytozome_metadata:
        n_panther = sum(1 for v in phytozome_metadata.values() if v.get("panther_descs"))
        execution_logs.append(ExecutionLogEntry(
            phase=4, database="Phytozome (JGI BioMart)",
            status="success", hits=len(phytozome_metadata),
            message=f"Phytozome annotation: {len(phytozome_metadata)} maize genes "
                    f"({n_panther} with Panther family). Source: phytozome-next.jgi.doe.gov",
        ))
    else:
        execution_logs.append(ExecutionLogEntry(
            phase=4, database="Phytozome (JGI BioMart)",
            status="info", hits=0,
            message="No Phytozome annotation returned (BioMart may be unreachable or genes not v5 NAM IDs).",
        ))

    # CornCyc pathway annotation for the resolved compound — flows into the
    # UI's "Maize pathway context (CornCyc)" section and the HTML/CSV reports.
    corncyc_annotation = None
    if chemical_entity:
        from corncyc_lookup import corncyc_annotation_for_chebi
        corncyc_annotation = corncyc_annotation_for_chebi(chemical_entity.chebi_id)
        if corncyc_annotation:
            execution_logs.append(ExecutionLogEntry(
                phase=2, database="CornCyc (PMN)", status="success",
                hits=corncyc_annotation["n_pathways"],
                message=(
                    f"CornCyc maize pathway context: {corncyc_annotation['n_pathways']} pathway(s), "
                    f"{corncyc_annotation['n_maize_genes']} annotated maize gene(s)"
                    + (f". Top: '{corncyc_annotation['pathways'][0]['common_name']}'"
                       if corncyc_annotation['pathways'] else "")
                ),
            ))

    return {
        "input_data": input_data.model_dump(),
        "chemical_entity": chemical_entity.model_dump() if chemical_entity else None,
        "reactions": [r.model_dump() for r in reactions],
        "proteins": [p.model_dump() for p in proteins],
        "orthologs": [o.model_dump() for o in orthologs],
        "unvalidated_targets": [o.model_dump() for o in unvalidated_targets],
        "targets": [t.model_dump() for t in targets],
        "domain_targets": [d.model_dump() for d in domain_targets],
        "advanced_homology_targets": [a.model_dump() for a in advanced_homology_targets],
        "maize_gene_metadata": maize_gene_metadata,
        "phytozome_metadata": phytozome_metadata,
        "corncyc_annotation": corncyc_annotation,
        "execution_logs": [l.model_dump() for l in execution_logs]
    }

async def run_single_pipeline(input_data: MetaboliteInput) -> dict:
    execution_logs: List[ExecutionLogEntry] = []
    return await _run_pipeline_core(input_data, execution_logs)

@app.post("/api/run_pipeline/stream")
async def run_pipeline_stream(input_data: MetaboliteInput):
    """
    Streams pipeline progress as Server-Sent Events.
    Each event is a JSON object on a `data:` line:
      {"type": "log",    "entry": {phase, database, status, hits, message, timestamp}}
      {"type": "result", "data": {...full pipeline result...}}
      {"type": "error",  "message": "...", "exception": "ClassName"}
      {"type": "done"}
    """
    queue: asyncio.Queue = asyncio.Queue()
    logs = StreamingLogList(queue)

    async def run_and_finish():
        import traceback
        try:
            result = await _run_pipeline_core(input_data, logs)
            await queue.put({"type": "result", "data": result})
        except Exception as e:
            tb = traceback.format_exc()
            await queue.put({
                "type": "error",
                "message": str(e) or repr(e),
                "exception": type(e).__name__,
                "traceback": tb,
            })
        finally:
            await queue.put({"type": "done"})

    task = asyncio.create_task(run_and_finish())

    async def gen():
        try:
            while True:
                msg = await queue.get()
                yield f"data: {json.dumps(msg)}\n\n"
                if msg.get("type") == "done":
                    break
        finally:
            if not task.done():
                task.cancel()

    return StreamingResponse(
        gen(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )

@app.get("/api/pipeline_config/schema")
def pipeline_config_schema():
    """
    Defaults + per-field descriptions (drawn from the Pydantic Field metadata)
    so the Configuration tab can render labels, help text, and validation
    bounds without duplicating them in JS.
    """
    fields = []
    for name, info in PipelineConfig.model_fields.items():
        # Extract numeric bounds out of Pydantic's constraint metadata
        ge = le = gt = lt = None
        for meta in getattr(info, "metadata", []) or []:
            for attr in ("ge", "le", "gt", "lt"):
                v = getattr(meta, attr, None)
                if v is not None:
                    locals_ = {"ge": ge, "le": le, "gt": gt, "lt": lt}
                    locals_[attr] = v
                    ge, le, gt, lt = locals_["ge"], locals_["le"], locals_["gt"], locals_["lt"]
        annotation = info.annotation
        type_name = getattr(annotation, "__name__", str(annotation))
        fields.append({
            "name": name,
            "default": info.default,
            "description": info.description or "",
            "type": type_name,
            "ge": ge, "le": le, "gt": gt, "lt": lt,
        })
    return {"fields": fields, "defaults": PipelineConfig().model_dump()}


@app.get("/api/corncyc/status")
def corncyc_status():
    """Report whether the CornCyc PGDB is installed locally + summary stats."""
    from corncyc_loader import get_status
    return get_status()


@app.get("/api/maize_afdb/status")
def maize_afdb_status():
    """Report whether the Foldseek-indexed maize AlphaFold DB is ready."""
    from install_maize_afdb import get_status
    return get_status()


@app.post("/api/maize_afdb/install")
async def maize_afdb_install():
    """
    Trigger the maize AlphaFold download + Foldseek index build, streaming
    progress as Server-Sent Events. Safe to call when already built (will
    short-circuit and emit a single `complete` event).

    Event shapes:
      {"type": "progress", "message": "...", "stage": "download|extract|index", "pct": 0..100}
      {"type": "complete"}
      {"type": "error", "message": "..."}
      {"type": "done"}
    """
    from install_maize_afdb import ensure_db_ready, is_db_ready

    queue: asyncio.Queue = asyncio.Queue()

    async def progress_cb(message, stage=None, pct=None):
        await queue.put({"type": "progress", "message": message, "stage": stage, "pct": pct})

    async def run_install():
        try:
            ok = await ensure_db_ready(interactive=False, progress_callback=progress_cb)
            if ok:
                await queue.put({"type": "complete"})
            else:
                await queue.put({"type": "error", "message": "Install did not complete; check server logs."})
        except Exception as e:
            await queue.put({"type": "error", "message": f"{type(e).__name__}: {e}"})
        finally:
            await queue.put({"type": "done"})

    if is_db_ready():
        async def short_gen():
            yield f"data: {json.dumps({'type':'complete','message':'Already built.'})}\n\n"
            yield f"data: {json.dumps({'type':'done'})}\n\n"
        return StreamingResponse(short_gen(), media_type="text/event-stream")

    task = asyncio.create_task(run_install())

    async def gen():
        try:
            while True:
                msg = await queue.get()
                yield f"data: {json.dumps(msg)}\n\n"
                if msg.get("type") == "done":
                    break
        finally:
            if not task.done():
                task.cancel()

    return StreamingResponse(
        gen(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.post("/api/validate_target")
async def validate_target(mapping: OrthologMapping):
    logs = []
    targets = await execute_phase5([mapping], logs=logs)
    domain_targets = await execute_phase6(targets, [mapping], logs=logs)
    advanced_homology_targets = await execute_phase7(targets, logs=logs)
    return {
        "target": targets[0].model_dump() if targets else None,
        "domain_target": domain_targets[0].model_dump() if domain_targets else None,
        "advanced_homology_target": advanced_homology_targets[0].model_dump() if advanced_homology_targets else None,
        "execution_logs": [l.model_dump() for l in logs]
    }

@app.post("/api/run_pipeline")
async def run_pipeline(input_data: MetaboliteInput):
    return await run_single_pipeline(input_data)

@app.get("/api/download/template")
def download_template():
    content = "mz,mode,adducts,tolerance_ppm\n117.0188,negative,M-H,5.0\n"
    return PlainTextResponse(content, media_type="text/csv", headers={"Content-Disposition": "attachment; filename=batch_template.csv"})

@app.get("/api/debug")
def get_debug():
    import phase3
    import sys
    import json
    return {"phase3_file": phase3.__file__, "sys_path": sys.path}

@app.post("/api/batch")
async def run_batch(file: UploadFile = File(...)):
    content = await file.read()
    reader = csv.DictReader(io.StringIO(content.decode("utf-8")))
    
    results = []
    # Process sequentially to avoid API rate limits on massive batches
    for row in reader:
        try:
            adducts_raw = row.get("adducts", "")
            adducts = adducts_raw.split(";") if adducts_raw else []
            
            inp = MetaboliteInput(
                mz=float(row["mz"]),
                mode=row.get("mode", "negative"),
                adducts=adducts,
                tolerance_ppm=float(row.get("tolerance_ppm", 5.0))
            )
            res = await run_single_pipeline(inp)
            results.append(res)
        except Exception as e:
            # Append failure note but continue
            results.append({
                "input_data": row,
                "error": str(e)
            })
            
    # Generate reports
    csv_content = generate_csv_report(results)
    html_content = generate_html_report(results)
    
    csv_id = save_report(csv_content, "csv")
    html_id = save_report(html_content, "html")
    
    return {
        "status": "completed",
        "processed": len(results),
        "csv_report_id": csv_id,
        "html_report_id": html_id,
        "results": results
    }

@app.post("/api/report/single")
async def create_single_report(result_data: dict):
    # Expects a list containing a single result, or just wraps it
    results = [result_data]
    csv_content = generate_csv_report(results)
    html_content = generate_html_report(results)
    
    csv_id = save_report(csv_content, "csv")
    html_id = save_report(html_content, "html")
    
    return {
        "csv_report_id": csv_id,
        "html_report_id": html_id
    }

@app.get("/api/reports/{report_id}")
def get_report(report_id: str):
    # We don't know the extension just from the ID in the route, so we check both
    for ext in ["csv", "html"]:
        filepath = os.path.join(REPORTS_DIR, f"{report_id}.{ext}")
        if os.path.exists(filepath):
            media_type = "text/csv" if ext == "csv" else "text/html"
            return FileResponse(filepath, media_type=media_type, filename=f"report_{report_id}.{ext}")
    return {"error": "Report not found"}
