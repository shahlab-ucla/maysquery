// Tab Switching Logic
document.getElementById('tab-single').addEventListener('click', (e) => {
    e.target.classList.add('active');
    document.getElementById('tab-batch').classList.remove('active');
    document.getElementById('pipeline-form-container').style.display = 'block';
    document.getElementById('batch-form').style.display = 'none';
    document.getElementById('tracker-section').style.display = 'block';
    document.getElementById('results-section').style.display = 'block';
});

document.getElementById('tab-batch').addEventListener('click', (e) => {
    e.target.classList.add('active');
    document.getElementById('tab-single').classList.remove('active');
    document.getElementById('pipeline-form-container').style.display = 'none';
    document.getElementById('batch-form').style.display = 'block';
    document.getElementById('tracker-section').style.display = 'none';
    document.getElementById('results-section').style.display = 'none';
});

// Subtab Switching Logic
const setQueryType = (type, activeBtnId, showInputsId) => {
    document.getElementById('query_type').value = type;
    ['subtab-mz', 'subtab-chem', 'subtab-ec'].forEach(id => document.getElementById(id).classList.remove('active'));
    document.getElementById(activeBtnId).classList.add('active');
    
    ['inputs-mz', 'inputs-chem', 'inputs-ec'].forEach(id => document.getElementById(id).style.display = 'none');
    document.getElementById(showInputsId).style.display = 'block';
    
    // Update tracker steps based on type
    document.getElementById('step-1').style.opacity = type === 'ec' || type === 'chemical' ? '0.3' : '1';
    document.getElementById('step-2').style.opacity = type === 'ec' ? '0.3' : '1';
};

document.getElementById('subtab-mz').addEventListener('click', () => setQueryType('mz', 'subtab-mz', 'inputs-mz'));
document.getElementById('subtab-chem').addEventListener('click', () => setQueryType('chemical', 'subtab-chem', 'inputs-chem'));
document.getElementById('subtab-ec').addEventListener('click', () => setQueryType('ec', 'subtab-ec', 'inputs-ec'));



let currentSingleResult = null;

// Single Query Form
document.getElementById('pipeline-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    document.getElementById('results-section').style.display = 'block';
    document.getElementById('single-report-downloads').style.display = 'none';
    
    const query_type = document.getElementById('query_type').value;
    const mz = parseFloat(document.getElementById('mz').value) || 0;
    const mode = document.getElementById('mode').value;
    const chemical_name = document.getElementById('chemical_name').value;
    const ec_number = document.getElementById('ec_number').value;
    
    const btn = e.target.querySelector('button[type="submit"]');
    btn.textContent = 'Running...';
    btn.disabled = true;

    document.querySelectorAll('.step').forEach(el => {
        el.classList.remove('active', 'completed');
    });

    const steps = ['step-1', 'step-2', 'step-3', 'step-4', 'step-5', 'step-6'];
    
    for (let i = 0; i < steps.length; i++) {
        const currentStep = document.getElementById(steps[i]);
        // Skip visual animation for bypassed steps
        if ((query_type === 'chemical' && i === 0) || (query_type === 'ec' && i < 2)) {
            currentStep.classList.add('completed');
            continue;
        }
        
        if (i > 0) {
            document.getElementById(steps[i-1]).classList.remove('active');
            document.getElementById(steps[i-1]).classList.add('completed');
        }
        currentStep.classList.add('active');
        await new Promise(r => setTimeout(r, 600));
    }
    document.getElementById(steps[steps.length-1]).classList.remove('active');
    document.getElementById(steps[steps.length-1]).classList.add('completed');

    try {
        const payload = { query_type, mz, mode, adducts: [], tolerance_ppm: 5.0 };
        if (query_type === 'chemical') payload.chemical_name = chemical_name;
        if (query_type === 'ec') payload.ec_number = ec_number;
        
        const response = await fetch('/api/run_pipeline', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        
        if (!response.ok) {
            const errData = await response.json().catch(() => ({}));
            throw new Error(`Pipeline failed (Status ${response.status}): ${JSON.stringify(errData)}`);
        }
        const data = await response.json();
        
        if (!data || !data.input_data) {
            throw new Error('Server returned malformed data: ' + JSON.stringify(data));
        }
        
        currentSingleResult = data;
        
        renderResults(data);
        document.getElementById('single-report-downloads').style.display = 'flex';
        
    } catch (err) {
        console.error(err);
        document.getElementById('dashboard-container').textContent = 'Error: ' + err.message;
    } finally {
        btn.textContent = 'Run Pipeline';
        btn.disabled = false;
    }
});

// Single Query Reports
const downloadSingleReport = async (ext) => {
    if (!currentSingleResult) return;
    const res = await fetch('/api/report/single', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(currentSingleResult)
    });
    if (res.ok) {
        const data = await res.json();
        const id = ext === 'html' ? data.html_report_id : data.csv_report_id;
        window.open(`/api/reports/${id}`, '_blank');
    }
};

document.getElementById('dl-single-html').addEventListener('click', () => downloadSingleReport('html'));
document.getElementById('dl-single-csv').addEventListener('click', () => downloadSingleReport('csv'));

// Batch Form
document.getElementById('batch-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const fileInput = document.getElementById('csv-file');
    const btn = e.target.querySelector('button');
    const statusDiv = document.getElementById('batch-status');
    const dlDiv = document.getElementById('batch-downloads');
    
    if (fileInput.files.length === 0) return;
    
    btn.textContent = 'Processing Batch...';
    btn.disabled = true;
    statusDiv.textContent = 'Uploading and processing CSV sequentially... This may take a while.';
    dlDiv.style.display = 'none';
    
    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    
    try {
        const response = await fetch('/api/batch', {
            method: 'POST',
            body: formData
        });
        
        if (!response.ok) throw new Error('Batch processing failed');
        const data = await response.json();
        
        statusDiv.textContent = `Successfully processed ${data.processed} queries!`;
        document.getElementById('batch-dl-html').href = `/api/reports/${data.html_report_id}`;
        document.getElementById('batch-dl-csv').href = `/api/reports/${data.csv_report_id}`;
        dlDiv.style.display = 'flex';
        
    } catch (err) {
        console.error(err);
        statusDiv.textContent = 'Error: ' + err.message;
        statusDiv.style.color = '#f87171';
    } finally {
        btn.textContent = 'Run Batch Pipeline';
        btn.disabled = false;
    }
});

function renderResults(data) {
    const dashboard = document.getElementById('dashboard-container');
    dashboard.innerHTML = ''; 

    const getPlddtColor = (score) => {
        if (score > 90) return '#1e3a8a';
        if (score > 70) return '#3b82f6';
        if (score > 50) return '#facc15';
        return '#ea580c';
    };

    let html = '';

    if (data.chemical_entity) {
        html += `
        <div class="dashboard-section">
            <h3>1. Chemical Entity</h3>
            <div class="data-card">
                <a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${data.chemical_entity.chebi_id}" target="_blank" class="badge badge-chebi">${data.chemical_entity.chebi_id}</a>
                ${data.chemical_entity.pubchem_cid ? `<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${data.chemical_entity.pubchem_cid}" target="_blank" class="badge" style="background:#0f766e; color:#fff; text-decoration:none;">CID: ${data.chemical_entity.pubchem_cid}</a>` : ''}
                <div><strong>Monoisotopic Mass:</strong> ${data.chemical_entity.monoisotopic_mass.toFixed(4)}</div>
            </div>
        </div>`;
    } else if (data.input_data.query_type === 'ec') {
        html += `
        <div class="dashboard-section">
            <h3>1. Input Source</h3>
            <div class="data-card">
                <span class="badge" style="background:#475569; color:#fff;">EC Number: ${data.input_data.ec_number}</span>
            </div>
        </div>`;
    }

    if (data.reactions && data.reactions.length > 0) {
        html += `<div class="dashboard-section"><h3>2. Reaction Networks (Rhea)</h3>`;
        data.reactions.forEach(r => {
            const rheaNum = r.rhea_id.replace('RHEA:', '');
            html += `<div class="data-card">
                <a href="https://www.rhea-db.org/rhea/${rheaNum}" target="_blank" class="badge badge-rhea">${r.rhea_id}</a>
                <span class="badge badge-pathway">${r.pathway_names[0] || 'Unknown Pathway'}</span>
                <span class="badge" style="background: ${r.is_transport ? '#6366f1' : '#475569'}; color:#fff">
                    ${r.is_transport ? 'Transport' : 'Metabolic'}
                </span>
            </div>`;
        });
        html += `</div>`;
    }

    if (data.proteins && data.proteins.length > 0) {
        html += `<div class="dashboard-section"><h3>3. Pan-life Protein Pool (UniProt)</h3>`;
        data.proteins.forEach(p => {
            const goTags = p.go_terms.map(go => `<div class="go-item"><a href="https://www.ebi.ac.uk/QuickGO/term/${go.id}" target="_blank" class="badge badge-go">${go.id}</a> <span class="go-name">${go.name}</span></div>`).join('');
            html += `<div class="data-card">
                <a href="https://www.uniprot.org/uniprotkb/${p.uniprot_accession}/entry" target="_blank" class="badge badge-uniprot">${p.uniprot_accession}</a>
                <span class="badge" style="background: ${p.category==='Enzyme' ? '#b91c1c' : '#047857'}; color:#fff">${p.category}</span>
                <div class="go-list">${goTags}</div>
            </div>`;
        });
        html += `</div>`;
    }

    if (data.targets && data.targets.length > 0) {
        html += `<div class="dashboard-section"><h3>5 & 6. Validated Zea mays Targets</h3>`;
        data.targets.forEach(t => {
            const domainTarget = data.domain_targets ? data.domain_targets.find(d => d.maize_gene_model === t.maize_gene_model) : null;
            html += `<div class="data-card">
                <a href="https://beta.plantplaza.org/gene/${t.maize_gene_model}" target="_blank" class="badge" style="background:#10b981; color:#fff; text-decoration:none; font-size:1.1em;">${t.maize_gene_model}</a>
                <div style="margin-top:10px;">
                    <strong style="color: ${getPlddtColor(t.tm_score * 100)}">Global TM-Score: ${t.tm_score.toFixed(2)}</strong><br>
                    <strong>Global pLDDT:</strong> ${t.plddt.toFixed(1)}
                </div>`;
            
            if (domainTarget) {
                html += `<div style="margin-top: 8px; padding: 8px; background: #f8fafc; border-left: 4px solid #3b82f6; border-radius: 4px;">
                    <strong style="color: ${getPlddtColor(domainTarget.domain_tm_score * 100)}; font-size: 1.1em;">Domain TM-Score: ${domainTarget.domain_tm_score.toFixed(2)}</strong><br>
                    <span style="font-size: 0.9em; color: #64748b;">Matched on <strong>${domainTarget.pfam_domain_name}</strong> (${domainTarget.pfam_domain_id})<br>
                    Residues: ${domainTarget.domain_start} - ${domainTarget.domain_end}</span>
                </div>`;
            }
                
            if (t.tissue_expression_fpkm && Object.keys(t.tissue_expression_fpkm).length > 0) {
                html += `<div style="margin-top:10px;"><strong>Key Tissue Expression (FPKM):</strong><ul>`;
                for (const [tissue, val] of Object.entries(t.tissue_expression_fpkm)) {
                    html += `<li>${tissue}: ${val}</li>`;
                }
                html += `</ul></div>`;
            }
            html += `</div>`;
        });
        html += `</div>`;
    }

    dashboard.innerHTML = html;

    const elements = [];
    const queryType = data.input_data.query_type;
    
    let rootId = null;
    
    if (queryType === 'mz') {
        const mzId = data.input_data.mz.toString();
        elements.push({ data: { id: mzId, label: 'm/z ' + mzId, type: 'mz' } });
        
        if (data.chemical_entity) {
            const chebiId = data.chemical_entity.chebi_id;
            elements.push({ data: { id: chebiId, label: chebiId, type: 'compound' } });
            elements.push({ data: { id: 'e1', source: mzId, target: chebiId } });
            rootId = chebiId;
        }
    } else if (queryType === 'chemical') {
        if (data.chemical_entity) {
            const chebiId = data.chemical_entity.chebi_id;
            elements.push({ data: { id: chebiId, label: chebiId, type: 'compound' } });
            rootId = chebiId;
        }
    } else if (queryType === 'ec') {
        const ecId = data.input_data.ec_number;
        elements.push({ data: { id: ecId, label: 'EC ' + ecId, type: 'reaction' } });
        rootId = ecId;
    }

    if (data.reactions && data.reactions.length > 0 && rootId) {
        data.reactions.forEach((r, i) => {
            const rId = r.rhea_id;
            if (!elements.find(e => e.data.id === rId)) {
                elements.push({ data: { id: rId, label: rId, type: 'reaction' } });
                elements.push({ data: { id: 'e1b_' + i, source: rootId, target: rId } });
            }
            
            if (data.proteins && data.proteins.length > 0) {
                data.proteins.forEach((p, j) => {
                    const pId = p.uniprot_accession;
                    if (!elements.find(e => e.data.id === pId)) {
                        elements.push({ data: { id: pId, label: pId, type: 'protein' } });
                    }
                    elements.push({ data: { id: 'e2_' + i + '_' + j, source: rId, target: pId } });
                    
                    if (data.orthologs && data.orthologs.length > 0) {
                        data.orthologs.forEach((o, k) => {
                            const oId = o.maize_gene_model;
                            if (o.query_uniprot_id === pId) {
                                if (!elements.find(e => e.data.id === oId)) {
                                    elements.push({ data: { id: oId, label: oId, type: 'gene' } });
                                }
                                elements.push({ data: { id: 'e3_' + i + '_' + j + '_' + k, source: pId, target: oId } });
                            }
                        });
                    }
                });
            }
        });
    } else if (queryType === 'ec' && data.proteins && data.proteins.length > 0) {
        // Direct link from EC to proteins
        data.proteins.forEach((p, j) => {
            const pId = p.uniprot_accession;
            if (!elements.find(e => e.data.id === pId)) {
                elements.push({ data: { id: pId, label: pId, type: 'protein' } });
            }
            elements.push({ data: { id: 'e2_ec_' + j, source: rootId, target: pId } });
            
            if (data.orthologs && data.orthologs.length > 0) {
                data.orthologs.forEach((o, k) => {
                    const oId = o.maize_gene_model;
                    if (o.query_uniprot_id === pId) {
                        if (!elements.find(e => e.data.id === oId)) {
                            elements.push({ data: { id: oId, label: oId, type: 'gene' } });
                        }
                        elements.push({ data: { id: 'e3_ec_' + j + '_' + k, source: pId, target: oId } });
                    }
                });
            }
        });
    }

    const cy = cytoscape({
        container: document.getElementById('cy'),
        elements: elements,
        style: [
            {
                selector: 'node',
                style: {
                    'label': 'data(label)',
                    'color': '#fff',
                    'text-outline-color': '#0f172a',
                    'text-outline-width': 2,
                    'text-valign': 'top',
                    'text-halign': 'center',
                    'font-size': '12px',
                    'cursor': 'pointer'
                }
            },
            {
                selector: 'node[type="mz"]',
                style: { 'background-color': '#f87171', 'shape': 'diamond' }
            },
            {
                selector: 'node[type="compound"]',
                style: { 'background-color': '#818cf8', 'shape': 'ellipse' }
            },
            {
                selector: 'node[type="reaction"]',
                style: { 'background-color': '#34d399', 'shape': 'round-rectangle' }
            },
            {
                selector: 'node[type="protein"]',
                style: { 'background-color': '#fbbf24', 'shape': 'hexagon' }
            },
            {
                selector: 'node[type="gene"]',
                style: { 'background-color': '#c084fc', 'shape': 'triangle' }
            },
            {
                selector: 'edge',
                style: {
                    'width': 2,
                    'line-color': '#475569',
                    'target-arrow-color': '#475569',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    'arrow-scale': 0.8
                }
            }
        ],
        layout: {
            name: 'dagre',
            rankDir: 'LR',
            spacingFactor: 1.2,
            animate: true
        }
    });

    cy.on('tap', 'node', function(evt){
        var node = evt.target;
        var type = node.data('type');
        var id = node.id();
        
        if (type === 'compound') {
            window.open(`https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${id}`, '_blank');
        } else if (type === 'reaction') {
            const rheaNum = id.replace('RHEA:', '');
            window.open(`https://www.rhea-db.org/rhea/${rheaNum}`, '_blank');
        } else if (type === 'protein') {
            window.open(`https://www.uniprot.org/uniprotkb/${id}/entry`, '_blank');
        } else if (type === 'gene') {
            window.open(`https://www.maizegdb.org/gene_center/gene/${id}`, '_blank');
        }
    });
}
