from fastapi import FastAPI, UploadFile, File
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, PlainTextResponse
import os
import csv
import io
from contextlib import asynccontextmanager
from models import MetaboliteInput, ValidatedTarget
from typing import List
from install_foldseek import check_and_prompt_foldseek
from report_generator import generate_csv_report, generate_html_report, save_report, REPORTS_DIR

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Prompt for foldseek installation on server startup
    check_and_prompt_foldseek()
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

async def run_single_pipeline(input_data: MetaboliteInput) -> dict:
    chemical_entity = None
    reactions = []
    proteins = []
    orthologs = []
    targets = []
    domain_targets = []
    
    if input_data.query_type == "mz":
        chemical_entity = await execute_phase1(input_data)
        reactions = await execute_phase2(chemical_entity)
        proteins = await execute_phase3(chemical_entity, reactions)
        orthologs = await execute_phase4(proteins)
        targets = await execute_phase5(orthologs)
    elif input_data.query_type == "chemical":
        chemical_entity = await fetch_chemical_entity_by_name(input_data.chemical_name)
        reactions = await execute_phase2(chemical_entity)
        proteins = await execute_phase3(chemical_entity, reactions)
        orthologs = await execute_phase4(proteins)
        targets = await execute_phase5(orthologs)
    elif input_data.query_type == "ec":
        proteins = await execute_phase3_by_ec(input_data.ec_number)
        orthologs = await execute_phase4(proteins)
        targets = await execute_phase5(orthologs)
        
    domain_targets = await execute_phase6(targets, orthologs)
        
    return {
        "input_data": input_data.model_dump(),
        "chemical_entity": chemical_entity.model_dump() if chemical_entity else None,
        "reactions": [r.model_dump() for r in reactions],
        "proteins": [p.model_dump() for p in proteins],
        "orthologs": [o.model_dump() for o in orthologs],
        "targets": [t.model_dump() for t in targets],
        "domain_targets": [d.model_dump() for d in domain_targets]
    }

@app.post("/api/run_pipeline")
async def run_pipeline(input_data: MetaboliteInput):
    return await run_single_pipeline(input_data)

@app.get("/api/download/template")
def download_template():
    content = "mz,mode,adducts,tolerance_ppm\n117.0188,negative,M-H,5.0\n"
    return PlainTextResponse(content, media_type="text/csv", headers={"Content-Disposition": "attachment; filename=batch_template.csv"})

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
