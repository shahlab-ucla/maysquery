import os
import csv
import io
import uuid
from jinja2 import Template
from typing import List, Dict

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports")
os.makedirs(REPORTS_DIR, exist_ok=True)

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
<title>DBsearch Pipeline Report</title>
<style>
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; background-color: #f9fafb; color: #111827; }
.container { max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); }
h1 { color: #2563eb; text-align: center; margin-bottom: 40px; }
h2 { color: #374151; border-bottom: 2px solid #e5e7eb; padding-bottom: 10px; margin-top: 30px; }
table { border-collapse: collapse; width: 100%; margin-bottom: 20px; font-size: 14px; }
th, td { border: 1px solid #e5e7eb; padding: 12px; text-align: left; }
th { background-color: #f3f4f6; font-weight: 600; color: #4b5563; }
tr:nth-child(even) { background-color: #f9fafb; }
.badge { display: inline-block; padding: 4px 8px; border-radius: 9999px; font-size: 12px; font-weight: 500; background-color: #dbeafe; color: #1e40af; }
</style>
</head>
<body>
<div class="container">
    <h1>Metabolite-to-Gene Pipeline Report</h1>
    {% for res in results %}
    <div class="query-block">
        <h2>Query: M/Z {{ res.input_data.mz }} <span class="badge">{{ res.input_data.mode }}</span></h2>
        <p><strong>Resolved Entity:</strong> {{ res.chemical_entity.chebi_id }} (Mass: {{ res.chemical_entity.monoisotopic_mass }})</p>
        
        {% if res.targets %}
        <table>
            <tr>
                <th>Maize Gene Model</th>
                <th>TM-Score (Structure)</th>
                <th>pLDDT (Confidence)</th>
            </tr>
            {% for target in res.targets %}
            <tr>
                <td>{{ target.maize_gene_model }}</td>
                <td>{{ target.tm_score }}</td>
                <td>{{ target.plddt }}</td>
            </tr>
            {% endfor %}
        </table>
        {% else %}
        <p><em>No structural targets validated for this query.</em></p>
        {% endif %}
    </div>
    {% endfor %}
</div>
</body>
</html>
"""

def generate_csv_report(results: List[Dict]) -> str:
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["Query_MZ", "Mode", "ChEBI", "Maize_Gene_Model", "TM_Score", "pLDDT"])
    
    for res in results:
        mz = res.get("input_data", {}).get("mz", "")
        mode = res.get("input_data", {}).get("mode", "")
        chebi = res.get("chemical_entity", {}).get("chebi_id", "")
        targets = res.get("targets", [])
        
        if not targets:
            writer.writerow([mz, mode, chebi, "None", "", ""])
            continue
            
        for t in targets:
            writer.writerow([mz, mode, chebi, t.get("maize_gene_model"), t.get("tm_score"), t.get("plddt")])
            
    return output.getvalue()

def generate_html_report(results: List[Dict]) -> str:
    template = Template(HTML_TEMPLATE)
    return template.render(results=results)

def save_report(content: str, ext: str) -> str:
    report_id = str(uuid.uuid4())
    filename = f"{report_id}.{ext}"
    filepath = os.path.join(REPORTS_DIR, filename)
    with open(filepath, "w") as f:
        f.write(content)
    return report_id
