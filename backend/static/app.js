// ---------- Maize AlphaFold structural-discovery banner ----------
const afdbBanner = {
    root:   () => document.getElementById('maize-afdb-banner'),
    icon:   () => document.getElementById('afdb-banner-icon'),
    title:  () => document.getElementById('afdb-banner-title'),
    detail: () => document.getElementById('afdb-banner-detail'),
    button: () => document.getElementById('afdb-install-btn'),

    setState(stateClass, icon, title, detail, btnLabel) {
        const r = this.root(); if (!r) return;
        r.style.display = 'flex';
        r.classList.remove(
            'afdb-banner-loading','afdb-banner-ready','afdb-banner-missing',
            'afdb-banner-running','afdb-banner-error',
        );
        r.classList.add(stateClass);
        if (this.icon())   this.icon().textContent = icon;
        if (this.title())  this.title().textContent = title;
        if (this.detail()) this.detail().textContent = detail || '';
        const b = this.button();
        if (b) {
            if (btnLabel) { b.style.display = 'inline-block'; b.textContent = btnLabel; b.disabled = false; }
            else          { b.style.display = 'none'; }
        }
    },

    async refresh() {
        try {
            const resp = await fetch('/api/maize_afdb/status');
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const s = await resp.json();
            if (s.ready) {
                this.setState(
                    'afdb-banner-ready', '✓',
                    'Phase 4.5 enabled — structural discovery active',
                    `${s.pdb_count} maize AlphaFold structures indexed.`,
                    null,
                );
            } else {
                this.setState(
                    'afdb-banner-missing', '⚠',
                    'Phase 4.5 disabled — maize AlphaFold index not built',
                    'Click to download ~5 GB and build the Foldseek index (~20–40 min, one-time).',
                    'Build Index',
                );
            }
        } catch (e) {
            this.setState('afdb-banner-error', '✗', 'Could not check structural-discovery status', String(e));
        }
    },

    async install() {
        const b = this.button();
        if (b) { b.disabled = true; b.textContent = 'Building…'; }
        this.setState(
            'afdb-banner-running', '⏳',
            'Building maize structural-discovery index…',
            'Starting…', null,
        );
        // Mirror progress into the terminal too
        term.show();
        term.info('Starting maize AlphaFold + Foldseek index build…');

        try {
            const resp = await fetch('/api/maize_afdb/install', { method: 'POST', headers: { 'Accept': 'text/event-stream' } });
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);

            const reader = resp.body.getReader();
            const dec = new TextDecoder();
            let buf = '';
            let lastErr = null;
            while (true) {
                const { value, done } = await reader.read();
                if (done) break;
                buf += dec.decode(value, { stream: true });
                let sep;
                while ((sep = buf.indexOf('\n\n')) !== -1) {
                    const raw = buf.slice(0, sep);
                    buf = buf.slice(sep + 2);
                    const dataLines = raw.split('\n').filter(l => l.startsWith('data:')).map(l => l.slice(5).trimStart());
                    if (!dataLines.length) continue;
                    let msg;
                    try { msg = JSON.parse(dataLines.join('\n')); } catch { continue; }

                    if (msg.type === 'progress') {
                        const stage = msg.stage ? `[${msg.stage}] ` : '';
                        const pct   = (msg.pct != null) ? ` (${msg.pct.toFixed(1)}%)` : '';
                        const line  = `${stage}${msg.message}${pct}`;
                        term.info(line);
                        this.setState('afdb-banner-running', '⏳', 'Building structural-discovery index…', line, null);
                    } else if (msg.type === 'complete') {
                        term.success(msg.message || 'Index ready.');
                    } else if (msg.type === 'error') {
                        lastErr = msg.message;
                        term.error(`Install error: ${msg.message}`);
                    }
                }
            }
            if (lastErr) {
                this.setState('afdb-banner-error', '✗', 'Index build failed', lastErr, 'Retry');
            } else {
                await this.refresh();
            }
        } catch (e) {
            term.error(`Install request failed: ${e.message || e}`);
            this.setState('afdb-banner-error', '✗', 'Install request failed', String(e), 'Retry');
        }
    },
};

document.addEventListener('DOMContentLoaded', () => {
    afdbBanner.refresh();
    const b = afdbBanner.button();
    if (b) b.addEventListener('click', () => afdbBanner.install());
    corncycBanner.refresh();
});

// ---------- CornCyc availability banner ----------
const corncycBanner = {
    root:   () => document.getElementById('corncyc-banner'),
    icon:   () => document.getElementById('corncyc-banner-icon'),
    title:  () => document.getElementById('corncyc-banner-title'),
    detail: () => document.getElementById('corncyc-banner-detail'),

    setState(stateClass, icon, title, detail) {
        const r = this.root(); if (!r) return;
        r.style.display = 'flex';
        r.classList.remove('afdb-banner-loading','afdb-banner-ready','afdb-banner-missing','afdb-banner-error');
        r.classList.add(stateClass);
        if (this.icon())   this.icon().textContent = icon;
        if (this.title())  this.title().textContent = title;
        if (this.detail()) this.detail().textContent = detail || '';
    },

    async refresh() {
        try {
            const resp = await fetch('/api/corncyc/status');
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const s = await resp.json();
            if (s.available && s.loaded) {
                this.setState('afdb-banner-ready', '✓',
                    `CornCyc ${s.version} loaded — ${s.maize_genes} maize genes annotated`,
                    `${s.compounds} compounds (${s.chebi_mapped_compounds} ChEBI-mapped) · ${s.reactions} reactions · ${s.pathways} pathways`);
            } else if (s.available) {
                this.setState('afdb-banner-missing', '◷',
                    `CornCyc ${s.version} detected (will load on first query)`,
                    s.data_dir);
            } else {
                this.setState('afdb-banner-missing', '○',
                    'CornCyc curated annotation not installed',
                    `Drop the PMN-licensed PGDB at ${s.expected_path || 'corncyc/<version>/data/'} to enable the curated discovery lane and pathway context.`);
            }
        } catch (e) {
            this.setState('afdb-banner-error', '✗', 'CornCyc status check failed', String(e));
        }
    },
};

// ---------- Pipeline Configuration ----------
// Schema (defaults + bounds + descriptions) is fetched from /api/pipeline_config/schema.
// User overrides are persisted in localStorage under PIPELINE_CONFIG_KEY.
const PIPELINE_CONFIG_KEY = 'maysquery.pipeline_config';

const pipelineConfig = {
    schema: null,         // {fields: [...], defaults: {...}}
    overrides: {},        // user-set values; merged onto defaults

    load() {
        try {
            this.overrides = JSON.parse(localStorage.getItem(PIPELINE_CONFIG_KEY) || '{}');
        } catch (e) {
            this.overrides = {};
        }
    },
    save() {
        try { localStorage.setItem(PIPELINE_CONFIG_KEY, JSON.stringify(this.overrides)); }
        catch (e) { console.warn('Could not persist pipeline config:', e); }
    },
    resetAll() {
        this.overrides = {};
        this.save();
    },
    resetField(name) {
        delete this.overrides[name];
        this.save();
    },
    set(name, value) {
        const def = this.schema && this.schema.defaults ? this.schema.defaults[name] : undefined;
        if (value === def || value === '' || value == null) {
            delete this.overrides[name];
        } else {
            this.overrides[name] = value;
        }
        this.save();
    },
    /** The actual config object to send with each pipeline run. */
    current() {
        if (!this.schema) return null;
        return { ...this.schema.defaults, ...this.overrides };
    },
};

// Field → section mapping (drives the rendering order + grouping)
const CONFIG_SECTIONS = [
    { id: 'phase1', title: 'Phase 1 — Chemical Identification',
      subtitle: 'm/z resolution to ChEBI',
      fields: ['cmm_tolerance_ppm'] },
    { id: 'phase2', title: 'Phase 2 — Reaction Networks',
      subtitle: 'Rhea SPARQL + KEGG pathway lookup',
      fields: ['rhea_fetch_limit', 'chebi_expansion_depth'] },
    { id: 'phase3', title: 'Phase 3 — Protein Pool',
      subtitle: 'UniProtKB Swiss-Prot search',
      fields: ['uniprot_size_per_category'] },
    { id: 'phase4', title: 'Phase 4 — Sequence Orthology',
      subtitle: 'Ensembl + PLAZA + BioMart, with local HMMER fallback',
      fields: ['hmmer_e_value'] },
    { id: 'phase4_5', title: 'Phase 4.5 — Structural Discovery',
      subtitle: 'Foldseek vs indexed maize AlphaFold proteome',
      fields: ['foldseek_tm_threshold', 'foldseek_max_hits_per_query', 'foldseek_concurrency'] },
    { id: 'phase5', title: 'Phase 5 — Enrichment',
      subtitle: 'pLDDT + Expression Atlas + 1-to-1 Foldseek for sequence-only hits',
      fields: ['enrichment_top_n', 'plddt_threshold'] },
    { id: 'phase7', title: 'Phase 7 — Pan-Plant Compara',
      subtitle: 'Ensembl Compara cross-check',
      fields: ['compara_max_orthologs_per_target'] },
];

async function initConfigTab() {
    pipelineConfig.load();
    try {
        const r = await fetch('/api/pipeline_config/schema');
        if (!r.ok) throw new Error(`HTTP ${r.status}`);
        pipelineConfig.schema = await r.json();
        renderConfigForm();
    } catch (e) {
        document.getElementById('config-sections').innerHTML =
            `<div class="config-loading" style="color:#f87171;">Failed to load config schema: ${e.message}</div>`;
    }
    const resetBtn = document.getElementById('config-reset-all');
    if (resetBtn) resetBtn.addEventListener('click', () => {
        if (!confirm('Reset all pipeline parameters to their defaults?')) return;
        pipelineConfig.resetAll();
        renderConfigForm();
    });
}

function renderConfigForm() {
    const container = document.getElementById('config-sections');
    if (!container || !pipelineConfig.schema) return;
    const fieldByName = Object.fromEntries(pipelineConfig.schema.fields.map(f => [f.name, f]));

    const html = CONFIG_SECTIONS.map(sec => {
        const rows = sec.fields.map(fname => {
            const f = fieldByName[fname];
            if (!f) return '';
            return renderConfigField(f);
        }).join('');
        return `<div class="config-section">
            <div class="config-section-header">
                <span class="config-section-title">${sec.title}</span>
                <span class="config-section-subtitle">${sec.subtitle}</span>
            </div>
            ${rows}
        </div>`;
    }).join('');
    container.innerHTML = html;

    // Wire change/input listeners
    container.querySelectorAll('[data-config-field]').forEach(inp => {
        const name = inp.getAttribute('data-config-field');
        const display = container.querySelector(`[data-config-display="${name}"]`);
        const def = pipelineConfig.schema.defaults[name];
        const isLogSlider = inp.getAttribute('data-log-slider') === 'true';

        const onInput = () => {
            let val;
            if (isLogSlider) {
                const exp = parseInt(inp.value, 10);
                val = Math.pow(10, -exp);
                if (display) display.textContent = `1e-${exp}`;
            } else if (inp.type === 'range' || inp.type === 'number') {
                val = parseFloat(inp.value);
                if (Number.isNaN(val)) val = def;
                if (display) display.textContent = formatConfigValue(val, fieldByName[name]);
            } else {
                val = inp.value;
                if (display) display.textContent = val;
            }
            pipelineConfig.set(name, val);
            markDirty(name, fieldByName[name]);
        };
        inp.addEventListener('input', onInput);

        const resetBtn = container.querySelector(`[data-config-reset="${name}"]`);
        if (resetBtn) resetBtn.addEventListener('click', () => {
            pipelineConfig.resetField(name);
            applyValueToField(name, fieldByName[name]);
            markDirty(name, fieldByName[name]);
        });
    });
}

function renderConfigField(f) {
    const name = f.name;
    const current = pipelineConfig.current()[name];
    const def = f.default;
    const isLogSlider = name.endsWith('_e_value');  // hmmer_e_value-style log slider
    const isFloat = (typeof def === 'number') && !Number.isInteger(def);
    const ge = f.ge, le = f.le;
    const bounds = (ge != null && le != null) ? `range ${ge} – ${le}` : '';
    const dirty = (current !== def);

    let control;
    if (isLogSlider) {
        // Map e.g. 1e-5 → exp=5; clamp 0..20
        const curExp = Math.max(0, Math.min(20, Math.round(-Math.log10(current))));
        control = `
            <input type="range" min="0" max="20" step="1" value="${curExp}"
                   class="config-field-input" data-config-field="${name}" data-log-slider="true">
            <span class="config-field-value" data-config-display="${name}">1e-${curExp}</span>
        `;
    } else if (isFloat && ge != null && le != null && (le - ge) <= 1.5) {
        // Sliders for narrow float ranges (TM threshold, etc.)
        const step = 0.01;
        control = `
            <input type="range" min="${ge}" max="${le}" step="${step}" value="${current}"
                   class="config-field-input" data-config-field="${name}">
            <span class="config-field-value" data-config-display="${name}">${formatConfigValue(current, f)}</span>
        `;
    } else if (isFloat) {
        // Float number input (e.g. pLDDT 0..100)
        const step = 0.5;
        control = `
            <input type="range" min="${ge ?? 0}" max="${le ?? 100}" step="${step}" value="${current}"
                   class="config-field-input" data-config-field="${name}">
            <span class="config-field-value" data-config-display="${name}">${formatConfigValue(current, f)}</span>
        `;
    } else {
        // Integer number input
        control = `
            <input type="number" min="${ge ?? ''}" max="${le ?? ''}" step="1" value="${current}"
                   class="config-field-input" data-config-field="${name}">
            <span class="config-field-bounds">${bounds}</span>
        `;
    }

    return `<div class="config-field" data-field-row="${name}">
        <div class="config-field-head">
            <label class="config-field-label" for="cfg-${name}">${humanizeFieldName(name)}</label>
            <span class="config-field-default ${dirty ? 'dirty' : ''}" data-config-default-label="${name}">
                default: ${formatConfigValue(def, f)}
            </span>
            <button type="button" class="config-field-reset ${dirty ? 'visible' : ''}"
                    data-config-reset="${name}">reset</button>
        </div>
        <div class="config-field-control">${control}</div>
        <div class="config-field-help">${f.description}</div>
    </div>`;
}

function applyValueToField(name, fieldDef) {
    const cur = pipelineConfig.current()[name];
    const inp = document.querySelector(`[data-config-field="${name}"]`);
    const display = document.querySelector(`[data-config-display="${name}"]`);
    if (!inp) return;
    if (inp.getAttribute('data-log-slider') === 'true') {
        const exp = Math.max(0, Math.min(20, Math.round(-Math.log10(cur))));
        inp.value = exp;
        if (display) display.textContent = `1e-${exp}`;
    } else {
        inp.value = cur;
        if (display) display.textContent = formatConfigValue(cur, fieldDef);
    }
}

function markDirty(name, fieldDef) {
    const cur = pipelineConfig.current()[name];
    const isDirty = (cur !== fieldDef.default);
    const def = document.querySelector(`[data-config-default-label="${name}"]`);
    const reset = document.querySelector(`[data-config-reset="${name}"]`);
    if (def) def.classList.toggle('dirty', isDirty);
    if (reset) reset.classList.toggle('visible', isDirty);
}

function humanizeFieldName(name) {
    return name
        .replace(/_/g, ' ')
        .replace(/\bcmm\b/g, 'CMM')
        .replace(/\bppm\b/g, 'ppm')
        .replace(/\bppm\b/g, 'ppm')
        .replace(/\btm\b/g, 'TM')
        .replace(/\bplddt\b/g, 'pLDDT')
        .replace(/\brhea\b/g, 'Rhea')
        .replace(/\bchebi\b/g, 'ChEBI')
        .replace(/\buniprot\b/g, 'UniProt')
        .replace(/\bhmmer\b/g, 'HMMER')
        .replace(/\bcompara\b/g, 'Compara')
        .replace(/\bfoldseek\b/g, 'Foldseek')
        .replace(/\be value\b/g, 'E-value')
        .replace(/^./, c => c.toUpperCase());
}

function formatConfigValue(v, fieldDef) {
    if (v == null) return '—';
    if (typeof v === 'number') {
        if (Math.abs(v) < 0.001 && v !== 0) return v.toExponential(1);
        if (Number.isInteger(v)) return String(v);
        return v.toFixed(2);
    }
    return String(v);
}

document.addEventListener('DOMContentLoaded', () => { initConfigTab(); });

// ---------- Discovery-lane classification (mirrors report_generator.classify_ortholog) ----------
const STRUCTURAL_SOURCE = 'Foldseek-structural';
const CURATED_SOURCE    = 'CornCyc';

function classifyOrtholog(sources) {
    const s = new Set(sources || []);
    const hasStruct  = s.has(STRUCTURAL_SOURCE);
    const hasCurated = s.has(CURATED_SOURCE);
    s.delete(STRUCTURAL_SOURCE);
    s.delete(CURATED_SOURCE);
    const hasSeq = s.size > 0;
    const lanesHit = (hasSeq ? 1 : 0) + (hasStruct ? 1 : 0) + (hasCurated ? 1 : 0);
    if (lanesHit >= 2) return 'consensus';
    if (hasStruct)     return 'structure_only';
    if (hasSeq)        return 'sequence_only';
    if (hasCurated)    return 'curated_only';
    return 'sequence_only';
}

const laneMeta = {
    consensus:      { label: 'Consensus (≥2 lanes)',      icon: '◆', cls: 'lane-consensus' },
    sequence_only:  { label: 'Sequence-based',            icon: '◐', cls: 'lane-sequence' },
    structure_only: { label: 'Structure-based',           icon: '◑', cls: 'lane-structure' },
    curated_only:   { label: 'CornCyc curated',           icon: '★', cls: 'lane-curated' },
};

function buildEnrichmentRows(orthologs, data) {
    const targetsByGene = Object.fromEntries((data.targets || []).map(t => [t.maize_gene_model, t]));
    const domainsByGene = Object.fromEntries((data.domain_targets || []).map(d => [d.maize_gene_model, d]));
    const advByGene     = Object.fromEntries((data.advanced_homology_targets || []).map(a => [a.maize_gene_model, a]));

    return (orthologs || []).map(o => {
        const sources = o.sources || [];
        const cls = classifyOrtholog(sources);
        const gene = o.maize_gene_model;
        const sim = Number(o.similarity_score || 0);
        const t = targetsByGene[gene] || null;
        const d = domainsByGene[gene] || null;
        const adv = advByGene[gene] || null;
        const topCompara = (adv && adv.ensembl_orthologs && adv.ensembl_orthologs[0]) || null;

        const seqSim = cls === 'structure_only' ? null : sim;
        let tm = (cls === 'sequence_only') ? null : (sim / 100.0);
        if (t && typeof t.tm_score === 'number') tm = t.tm_score;

        return {
            gene, cls, sources,
            consensus_score: o.consensus_score,
            query_uniprot_id: o.query_uniprot_id,
            plaza_orthogroup: o.plaza_orthogroup,
            seq_sim: seqSim,
            tm,
            plddt: t ? t.plddt : null,
            expr_count:        t ? (t.n_expression_experiments || 0) : 0,
            expr_experiments:  t ? (t.expression_experiments || []) : [],
            enrichment_kind:   t ? (t.enrichment_kind || 'full') : null,
            domain: d,
            top_compara: topCompara,
            num_compara: adv ? (adv.ensembl_orthologs || []).length : 0,
            enriched: !!t,
            raw_ortholog: o,
        };
    });
}

// ---------- Maize gene metadata helpers (Gramene-sourced; populated per pipeline run) ----------
let _maizeGeneMeta = {};  // {gene_id: {symbol, description, synonyms, biotype}}

function setMaizeGeneMeta(meta) { _maizeGeneMeta = meta || {}; }

function geneMeta(gene_id) { return _maizeGeneMeta[gene_id] || null; }

/**
 * Render a maize gene anchor: clickable badge linking to MaizeGDB, with the
 * Gramene symbol + description shown alongside when available.
 * `variant`: 'badge' (default, with chip styling) or 'inline'.
 */
function renderMaizeGeneLabel(gene_id, variant = 'badge', titleSuffix = '') {
    if (!gene_id) return '—';
    const m = geneMeta(gene_id);
    const isPlaceholder = gene_id.startsWith('UNIPROT:');
    const url = isPlaceholder
        ? `https://www.maizegdb.org/search?query=${encodeURIComponent(gene_id.replace('UNIPROT:',''))}`
        : dbUrl.maizeGene(gene_id);
    const synList = m && m.synonyms && m.synonyms.length ? ` · synonyms: ${m.synonyms.join(', ')}` : '';
    const title = `Open ${gene_id} in MaizeGDB${synList}${titleSuffix ? ' · ' + titleSuffix : ''}`;
    const idLink = `<a href="${url}" target="_blank" class="badge ${variant === 'badge' ? 'badge-gene' : ''} gene-link" title="${escapeHtmlAttr(title)}">${gene_id}</a>`;
    if (!m) return idLink;
    const sym = m.symbol || '';
    const desc = m.description || '';
    if (!sym && !desc) return idLink;
    const label = sym && desc ? `${sym} — ${desc}` : (sym || desc);
    return `${idLink}<span class="gene-meta-label" title="${escapeHtmlAttr(synList ? synList.slice(3) : '')}">${escapeHtml(label)}</span>`;
}

function escapeHtmlAttr(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

// ---------- External-DB URL helpers ----------
const dbUrl = {
    maizeGene: (id) => `https://www.maizegdb.org/gene_center/gene/${encodeURIComponent(id)}`,
    uniprot:   (id) => `https://www.uniprot.org/uniprotkb/${encodeURIComponent(id)}/entry`,
    chebi:     (id) => `https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${encodeURIComponent(id)}`,
    pubchem:   (id) => `https://pubchem.ncbi.nlm.nih.gov/compound/${encodeURIComponent(id)}`,
    rhea:      (id) => `https://www.rhea-db.org/rhea/${encodeURIComponent(String(id).replace('RHEA:', ''))}`,
    pfam:      (id) => `https://www.ebi.ac.uk/interpro/entry/pfam/${encodeURIComponent(id)}/`,
    quickgo:   (id) => `https://www.ebi.ac.uk/QuickGO/term/${encodeURIComponent(id)}`,
    gxaExperiment: (id) => `https://www.ebi.ac.uk/gxa/experiments/${encodeURIComponent(id)}`,
    ensemblPlantsGene: (species, geneId) => {
        if (!species) return `https://plants.ensembl.org/Multi/Search/Results?q=${encodeURIComponent(geneId)}`;
        const sp = species.charAt(0).toUpperCase() + species.slice(1);
        return `https://plants.ensembl.org/${encodeURIComponent(sp)}/Gene/Summary?g=${encodeURIComponent(geneId)}`;
    },
    ensemblSpecies: (species) => {
        if (!species) return 'https://plants.ensembl.org/';
        const sp = species.charAt(0).toUpperCase() + species.slice(1);
        return `https://plants.ensembl.org/${encodeURIComponent(sp)}/Info/Index`;
    },
};

// ---------- Tab Switching ----------
const tabs = {
    single: { btn: 'tab-single', body: 'pipeline-form-container', showTracker: true,  showResults: true  },
    config: { btn: 'tab-config', body: 'config-form-container',   showTracker: false, showResults: false },
    batch:  { btn: 'tab-batch',  body: 'batch-form',              showTracker: false, showResults: false },
};
function activateTab(name) {
    Object.entries(tabs).forEach(([k, t]) => {
        const btn  = document.getElementById(t.btn);
        const body = document.getElementById(t.body);
        if (btn)  btn.classList.toggle('active', k === name);
        if (body) body.style.display = (k === name) ? 'block' : 'none';
    });
    document.getElementById('tracker-section').style.display = tabs[name].showTracker ? 'block' : 'none';
    document.getElementById('results-section').style.display = tabs[name].showResults ? 'block' : 'none';
}
document.getElementById('tab-single').addEventListener('click', () => activateTab('single'));
document.getElementById('tab-config').addEventListener('click', () => activateTab('config'));
document.getElementById('tab-batch') .addEventListener('click', () => activateTab('batch'));

// Convenience: the "see Configuration tab" link in the single-query help text
const openConfigLink = document.getElementById('open-config-link');
if (openConfigLink) openConfigLink.addEventListener('click', (e) => { e.preventDefault(); activateTab('config'); });

// ---------- Sub-tab Switching ----------
const setQueryType = (type, activeBtnId, showInputsId) => {
    document.getElementById('query_type').value = type;
    ['subtab-mz', 'subtab-chem', 'subtab-ec'].forEach(id => document.getElementById(id).classList.remove('active'));
    document.getElementById(activeBtnId).classList.add('active');

    ['inputs-mz', 'inputs-chem', 'inputs-ec'].forEach(id => document.getElementById(id).style.display = 'none');
    document.getElementById(showInputsId).style.display = 'block';

    document.getElementById('step-1').style.opacity = type === 'ec' || type === 'chemical' ? '0.3' : '1';
    document.getElementById('step-2').style.opacity = type === 'ec' ? '0.3' : '1';
};

document.getElementById('subtab-mz').addEventListener('click', () => setQueryType('mz', 'subtab-mz', 'inputs-mz'));
document.getElementById('subtab-chem').addEventListener('click', () => setQueryType('chemical', 'subtab-chem', 'inputs-chem'));
document.getElementById('subtab-ec').addEventListener('click', () => setQueryType('ec', 'subtab-ec', 'inputs-ec'));

// ---------- Terminal-style logger ----------
const term = {
    wrap: () => document.getElementById('terminal-wrap'),
    out:  () => document.getElementById('terminal-output'),

    show() {
        const w = this.wrap(); if (w) w.style.display = 'block';
    },
    clear() {
        const o = this.out(); if (o) o.innerHTML = '';
    },
    fmtTime(ts) {
        try {
            const d = ts ? new Date(ts) : new Date();
            const hh = String(d.getUTCHours()).padStart(2, '0');
            const mm = String(d.getUTCMinutes()).padStart(2, '0');
            const ss = String(d.getUTCSeconds()).padStart(2, '0');
            return `${hh}:${mm}:${ss}`;
        } catch (e) { return '--:--:--'; }
    },
    escape(s) {
        return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
    },
    line(html) {
        const o = this.out();
        if (!o) return;
        const div = document.createElement('div');
        div.className = 'terminal-line';
        div.innerHTML = html;
        o.appendChild(div);
        o.scrollTop = o.scrollHeight;
    },

    log(entry) {
        const t = this.fmtTime(entry.timestamp);
        const status = (entry.status || 'info').toLowerCase();
        const statusClass = `tl-${status}`;
        const tag = `[${status.toUpperCase()}]`;
        const hitsTxt = (entry.hits != null && entry.hits !== 0)
            ? ` <span class="tl-hits">(${entry.hits} hit${entry.hits === 1 ? '' : 's'})</span>` : '';
        this.line(
            `<span class="tl-time">${t}</span> ` +
            `<span class="tl-phase">[P${entry.phase}]</span> ` +
            `<span class="tl-db">${this.escape(entry.database)}</span> ` +
            `<span class="${statusClass}">${tag}</span> ` +
            `<span class="${statusClass}">${this.escape(entry.message || '')}</span>` +
            hitsTxt
        );
    },

    info(msg) {
        this.line(`<span class="tl-time">${this.fmtTime()}</span> <span class="tl-meta">▸ ${this.escape(msg)}</span>`);
    },
    error(msg, traceback) {
        this.line(`<span class="tl-time">${this.fmtTime()}</span> <span class="tl-error">✗ ${this.escape(msg)}</span>`);
        if (traceback) {
            this.line(`<span class="tl-meta">${this.escape(traceback)}</span>`);
        }
    },
    success(msg) {
        this.line(`<span class="tl-time">${this.fmtTime()}</span> <span class="tl-success">✓ ${this.escape(msg)}</span>`);
    },
};

const termClearBtn = document.getElementById('terminal-clear');
if (termClearBtn) termClearBtn.addEventListener('click', () => term.clear());

// ---------- Step Tracker (driven by real log events) ----------
const tracker = {
    reset(queryType) {
        for (let i = 1; i <= 7; i++) {
            const el = document.getElementById(`step-${i}`);
            if (el) {
                el.classList.remove('active', 'completed', 'errored');
                el.style.opacity = '1';
            }
        }
        if (queryType === 'chemical') document.getElementById('step-1').style.opacity = '0.3';
        if (queryType === 'ec') {
            document.getElementById('step-1').style.opacity = '0.3';
            document.getElementById('step-2').style.opacity = '0.3';
        }
    },
    currentPhase: 0,
    enterPhase(phase) {
        if (phase < 1 || phase > 7) return;
        // Mark every earlier phase that we've touched as completed
        for (let i = 1; i < phase; i++) {
            const el = document.getElementById(`step-${i}`);
            if (el && !el.classList.contains('errored')) {
                el.classList.remove('active');
                el.classList.add('completed');
            }
        }
        const cur = document.getElementById(`step-${phase}`);
        if (cur && !cur.classList.contains('errored')) {
            cur.classList.add('active');
        }
        this.currentPhase = Math.max(this.currentPhase, phase);
    },
    markError(phase) {
        const el = document.getElementById(`step-${phase}`);
        if (el) {
            el.classList.remove('active', 'completed');
            el.classList.add('errored');
        }
    },
    completeAll() {
        for (let i = 1; i <= 7; i++) {
            const el = document.getElementById(`step-${i}`);
            if (el && !el.classList.contains('errored') && el.style.opacity !== '0.3') {
                el.classList.remove('active');
                el.classList.add('completed');
            }
        }
    },
};

// ---------- Single Query Form ----------
let currentSingleResult = null;

document.getElementById('pipeline-form').addEventListener('submit', async (e) => {
    e.preventDefault();
    document.getElementById('results-section').style.display = 'block';
    const reportRow = document.getElementById('single-report-downloads');
    if (reportRow) reportRow.style.display = 'none';

    const queryType = document.getElementById('query_type').value;
    const cfg = pipelineConfig.current();   // current values from Configuration tab + localStorage

    const payload = {
        query_type: queryType,
        // Mirror the headline values into legacy top-level fields for backward compat.
        // The server reads from `pipeline_config` first.
        hmmer_e_value: cfg ? cfg.hmmer_e_value : 1e-5,
        tolerance_ppm: cfg ? cfg.cmm_tolerance_ppm : 5.0,
        mode: document.getElementById('mode')?.value || 'negative',
        adducts: [],
        chemical_name: queryType === 'chemical' ? document.getElementById('chemical_name').value : '',
        ec_number:     queryType === 'ec'       ? document.getElementById('ec_number').value     : '',
        pipeline_config: cfg,
    };
    if (queryType === 'mz') {
        payload.mz = parseFloat(document.getElementById('mz_value').value);
    }

    const btn = e.target.querySelector('button[type="submit"]');
    btn.textContent = 'Running...';
    btn.disabled = true;

    tracker.reset(queryType);
    term.show();
    term.clear();

    const dirtyEntries = pipelineConfig.schema
        ? Object.keys(pipelineConfig.overrides).filter(k =>
              pipelineConfig.overrides[k] !== pipelineConfig.schema.defaults[k])
        : [];
    const cfgSummary = dirtyEntries.length
        ? ` cfg-overrides=${dirtyEntries.length} (${dirtyEntries.join(', ')})`
        : ' cfg=defaults';

    term.info(`POST /api/run_pipeline/stream  query_type=${queryType}` +
        (queryType === 'mz' ? `  m/z=${payload.mz}` : '') +
        (queryType === 'chemical' ? `  name=${payload.chemical_name}` : '') +
        (queryType === 'ec' ? `  EC=${payload.ec_number}` : '') +
        cfgSummary);

    try {
        const response = await fetch('/api/run_pipeline/stream', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', 'Accept': 'text/event-stream' },
            body: JSON.stringify(payload),
        });

        if (!response.ok) {
            const errData = await response.json().catch(() => ({}));
            throw new Error(`Pipeline failed (Status ${response.status}): ${JSON.stringify(errData)}`);
        }

        // Stream parser — read SSE-style `data: <json>\n\n` chunks from the body.
        const reader = response.body.getReader();
        const decoder = new TextDecoder();
        let buffer = '';
        let finalResult = null;
        let serverError = null;

        while (true) {
            const { value, done } = await reader.read();
            if (done) break;
            buffer += decoder.decode(value, { stream: true });

            let sepIdx;
            while ((sepIdx = buffer.indexOf('\n\n')) !== -1) {
                const rawEvent = buffer.slice(0, sepIdx);
                buffer = buffer.slice(sepIdx + 2);

                // Each event is a sequence of lines; we use just the `data:` lines.
                const dataLines = rawEvent.split('\n')
                    .filter(l => l.startsWith('data:'))
                    .map(l => l.slice(5).trimStart());
                if (dataLines.length === 0) continue;

                let msg;
                try { msg = JSON.parse(dataLines.join('\n')); }
                catch (err) { console.warn('Bad SSE event JSON', dataLines, err); continue; }

                handleStreamMessage(msg);
                if (msg.type === 'result') finalResult = msg.data;
                if (msg.type === 'error')  serverError = msg;
            }
        }

        if (serverError) {
            term.error(`${serverError.exception || 'Error'}: ${serverError.message}`, serverError.traceback);
            throw new Error(`${serverError.exception || 'Error'}: ${serverError.message}`);
        }

        if (!finalResult || !finalResult.input_data) {
            throw new Error('Pipeline finished but no result event was received.');
        }

        tracker.completeAll();
        term.success('Pipeline finished.');

        currentSingleResult = finalResult;
        renderResults(finalResult);

        if (reportRow) reportRow.style.display = 'flex';

    } catch (err) {
        console.error(err);
        term.error(err.message || String(err));
        const dashboard = document.getElementById('dashboard-container');
        if (dashboard) {
            dashboard.classList.remove('hidden');
            dashboard.innerHTML = `<div class="data-card" style="border-color:#7f1d1d;background:rgba(127,29,29,0.15);color:#fecaca;">
                <strong>Pipeline Error</strong><br><span style="font-family:monospace;font-size:0.85rem;">${term.escape(err.message || String(err))}</span>
            </div>`;
        }
    } finally {
        btn.textContent = 'Run Pipeline';
        btn.disabled = false;
    }
});

function handleStreamMessage(msg) {
    if (msg.type === 'log' && msg.entry) {
        tracker.enterPhase(msg.entry.phase);
        if ((msg.entry.status || '').toLowerCase() === 'error') {
            tracker.markError(msg.entry.phase);
        }
        term.log(msg.entry);
    } else if (msg.type === 'result') {
        term.info('← received final result payload');
    } else if (msg.type === 'error') {
        term.error(`${msg.exception || 'Error'}: ${msg.message}`, msg.traceback);
    } else if (msg.type === 'done') {
        // no-op
    }
}

// ---------- Single Query Reports ----------
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
const dlHtmlBtn = document.getElementById('dl-single-html');
const dlCsvBtn  = document.getElementById('dl-single-csv');
if (dlHtmlBtn) dlHtmlBtn.addEventListener('click', () => downloadSingleReport('html'));
if (dlCsvBtn)  dlCsvBtn.addEventListener('click',  () => downloadSingleReport('csv'));

// ---------- Batch Form ----------
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
        const response = await fetch('/api/batch', { method: 'POST', body: formData });
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

// ---------- Results renderer ----------
function renderResults(data) {
    const dashboard = document.getElementById('dashboard-container');
    if (!dashboard) return;
    dashboard.innerHTML = '';

    setMaizeGeneMeta(data.maize_gene_metadata || {});

    // Build the per-ortholog enrichment rows once; the lane sections and summary share them.
    const allRowsPrecomputed = buildEnrichmentRows(data.orthologs || [], data);
    const lanes = {
        consensus:      allRowsPrecomputed.filter(r => r.cls === 'consensus'),
        structure_only: allRowsPrecomputed.filter(r => r.cls === 'structure_only'),
        sequence_only:  allRowsPrecomputed.filter(r => r.cls === 'sequence_only'),
        curated_only:   allRowsPrecomputed.filter(r => r.cls === 'curated_only'),
    };

    // ----- Executive summary -----
    let html = renderExecutiveSummary(data, lanes);

    // ----- Drill-down sections (collapsed by default for long ones) -----
    html += `<div class="drill-stack">`;

    // Chemical entity / input source
    html += renderDetail('Resolved chemical entity', renderChemicalEntityCard(data), true);

    // Reactions
    if (data.reactions && data.reactions.length > 0) {
        html += renderDetail(
            `Reactions <span class="detail-count">${data.reactions.length} from Rhea</span>`,
            renderReactionsBlock(data.reactions),
            false,
        );
    }

    // Pan-life proteins
    if (data.proteins && data.proteins.length > 0) {
        const cats = countBy(data.proteins, p => p.category);
        const catStr = Object.entries(cats).map(([k, v]) => `${v} ${k.toLowerCase()}${v===1?'':'s'}`).join(' · ');
        html += renderDetail(
            `Pan-life proteins <span class="detail-count">${data.proteins.length} unique (${catStr})</span>`,
            renderProteinsBlock(data.proteins),
            false,
        );
    }

    // Maize candidates — each lane in its own collapsible block
    const laneBody = `
        ${renderLaneSection('consensus', lanes.consensus,
            'Detected by ≥2 independent lanes (sequence, structure, or CornCyc curation) — highest confidence.')}
        ${renderLaneSection('structure_only', lanes.structure_only,
            'Recovered only by Foldseek against the maize AlphaFold proteome — "hidden orthologs" with divergent sequence but conserved fold.')}
        ${renderLaneSection('curated_only', lanes.curated_only,
            'Annotated by CornCyc (PlantCyc/PMN) as catalysing a reaction involving this compound, but not surfaced by sequence or structural homology search.')}
        ${renderLaneSection('sequence_only', lanes.sequence_only,
            'Detected only by sequence homology (Ensembl Compara / pan-homology, HMMER fallback).')}
    `;
    const laneCountsLabel = `${allRowsPrecomputed.length} total · ${lanes.consensus.length} consensus · ${lanes.structure_only.length} structure · ${lanes.curated_only.length} curated · ${lanes.sequence_only.length} sequence`;
    html += renderDetail(
        `Maize candidates <span class="detail-count">${laneCountsLabel}</span>`,
        laneBody,
        true,
    );

    // CornCyc maize-specific pathway context (when CornCyc is installed and the compound matched)
    if (data.corncyc_annotation) {
        html += renderDetail(
            `CornCyc maize pathway context <span class="detail-count">${data.corncyc_annotation.n_pathways} pathway(s), ${data.corncyc_annotation.n_maize_genes} annotated maize gene(s) · PMN CornCyc ${data.corncyc_annotation.version}</span>`,
            renderCornCycBlock(data.corncyc_annotation, data.maize_gene_metadata || {}),
            true,
        );
    }

    // Pan-plant Compara details (full per-target ortholog tables)
    if (data.advanced_homology_targets && data.advanced_homology_targets.length > 0) {
        html += renderDetail(
            `Pan-plant Compara <span class="detail-count">${data.advanced_homology_targets.length} target${data.advanced_homology_targets.length===1?'':'s'}</span>`,
            renderComparaBlock(data.advanced_homology_targets),
            false,
        );
    }

    // Execution log
    if (data.execution_logs && data.execution_logs.length > 0) {
        html += renderDetail(
            `Execution log <span class="detail-count">${data.execution_logs.length} events</span>`,
            renderExecutionLogTable(data.execution_logs),
            false,
        );
    }

    html += `</div>`;  // close drill-stack

    dashboard.innerHTML = html;
    dashboard.classList.remove('hidden');

    // Unvalidated targets render outside the drill stack
    renderUnvalidatedSection(data.unvalidated_targets || []);

    // Cytoscape graph (collapsed by default — see index.html `<details>` around #cy)
    renderCytoscape(data);
}

// ----- Helpers used by the new renderResults -----

function countBy(arr, keyFn) {
    const out = {};
    for (const x of arr) {
        const k = keyFn(x);
        out[k] = (out[k] || 0) + 1;
    }
    return out;
}

function renderDetail(summaryHTML, bodyHTML, openByDefault) {
    return `<details class="drill"${openByDefault ? ' open' : ''}>
        <summary class="drill-summary"><span class="drill-caret"></span><span class="drill-title">${summaryHTML}</span></summary>
        <div class="drill-body">${bodyHTML}</div>
    </details>`;
}

function renderExecutiveSummary(data, lanes) {
    const ce = data.chemical_entity;
    const inp = data.input_data || {};
    const rxns = data.reactions || [];
    const proteins = data.proteins || [];
    const allRows = [...lanes.consensus, ...lanes.structure_only, ...lanes.sequence_only];
    const top = allRows[0] || null;

    let queryLabel = '';
    if (inp.query_type === 'mz' && inp.mz != null) queryLabel = `m/z ${inp.mz} (${inp.mode||'?'})`;
    else if (inp.query_type === 'chemical')        queryLabel = `compound name "${inp.chemical_name||'?'}"`;
    else if (inp.query_type === 'ec')              queryLabel = `EC ${inp.ec_number||'?'}`;

    const resolvedHTML = ce
        ? `<a href="${dbUrl.chebi(ce.chebi_id)}" target="_blank" class="exec-link">${ce.chebi_id}</a>${ce.pubchem_cid ? ` · <a href="${dbUrl.pubchem(ce.pubchem_cid)}" target="_blank" class="exec-link">CID ${ce.pubchem_cid}</a>` : ''} · mass ${Number(ce.monoisotopic_mass||0).toFixed(4)}`
        : '<span class="empty-inline">no chemical entity resolved</span>';

    const topPathway = rxns.find(r => (r.pathway_names || []).length > 0);
    const pathwayCount = topPathway ? topPathway.pathway_names.length : 0;
    const pathwaySummary = topPathway
        ? `${topPathway.pathway_names[0]}${pathwayCount > 1 ? ` <span class="exec-dim">+${pathwayCount-1} more</span>` : ''}`
        : '<span class="empty-inline">no KEGG pathway mapped</span>';

    const cats = countBy(proteins, p => p.category);
    const proteinSummary = Object.entries(cats).map(([k, v]) => `${v} ${k.toLowerCase()}${v===1?'':'s'}`).join(' · ') || '0';

    let topRowHTML = '';
    if (top) {
        const tm = (top.tm != null && top.tm > 0) ? top.tm.toFixed(2) : '—';
        const pl = (top.plddt != null && top.plddt > 0) ? top.plddt.toFixed(1) : '—';
        const seq = (top.seq_sim != null) ? top.seq_sim.toFixed(1) + '%' : '—';
        const exprCount = top.expr_count || 0;
        const tissueLabel = exprCount ? `expressed in ${exprCount} Atlas exp${exprCount===1?'':'s'}` : 'no expression evidence';
        const laneIcon = laneMeta[top.cls].icon;
        const laneLbl  = laneMeta[top.cls].label;
        topRowHTML = `
        <div class="exec-top-row">
            <span class="exec-top-icon">⭐</span>
            <span class="exec-top-gene-wrap">${renderMaizeGeneLabel(top.gene, 'inline')}</span>
            <span class="lane-pill lane-${top.cls === 'sequence_only' ? 'sequence' : (top.cls === 'structure_only' ? 'structure' : 'consensus')}">${laneIcon} ${laneLbl}</span>
            <span class="exec-top-metrics">
                Struct TM <strong>${tm}</strong> · Seq sim <strong>${seq}</strong> · pLDDT <strong>${pl}</strong> · ${tissueLabel}
            </span>
        </div>`;
    }

    return `<section class="exec-summary glass-panel">
        <header class="exec-header">
            <div>
                <div class="exec-eyebrow">Query</div>
                <div class="exec-query">${queryLabel}</div>
            </div>
            <div class="exec-headline-counts">
                <div class="exec-count-cell"><div class="ec-val">${rxns.length}</div><div class="ec-lbl">reactions</div></div>
                <div class="exec-count-cell"><div class="ec-val">${proteins.length}</div><div class="ec-lbl">pan-life proteins</div></div>
                <div class="exec-count-cell exec-count-cell-emph"><div class="ec-val">${allRows.length}</div><div class="ec-lbl">maize candidates</div></div>
            </div>
        </header>
        <div class="exec-rows">
            <div class="exec-row"><span class="exec-row-label">Resolved chemical</span><span class="exec-row-val">${resolvedHTML}</span></div>
            <div class="exec-row"><span class="exec-row-label">Top KEGG pathway</span><span class="exec-row-val">${pathwaySummary}</span></div>
            <div class="exec-row"><span class="exec-row-label">Protein breakdown</span><span class="exec-row-val">${proteinSummary}</span></div>
            <div class="exec-row"><span class="exec-row-label">Discovery lanes</span><span class="exec-row-val">
                <span class="lane-pill lane-consensus">${laneMeta.consensus.icon} ${lanes.consensus.length} consensus</span>
                <span class="lane-pill lane-structure">${laneMeta.structure_only.icon} ${lanes.structure_only.length} structure-only</span>
                <span class="lane-pill lane-curated">${laneMeta.curated_only.icon} ${lanes.curated_only.length} CornCyc-only</span>
                <span class="lane-pill lane-sequence">${laneMeta.sequence_only.icon} ${lanes.sequence_only.length} sequence-only</span>
            </span></div>
            ${data.corncyc_annotation ? `
            <div class="exec-row"><span class="exec-row-label">CornCyc context</span><span class="exec-row-val">
                <b>${data.corncyc_annotation.n_pathways}</b> pathway${data.corncyc_annotation.n_pathways===1?'':'s'} ·
                <b>${data.corncyc_annotation.n_maize_genes}</b> curated maize gene${data.corncyc_annotation.n_maize_genes===1?'':'s'}
                ${data.corncyc_annotation.pathways && data.corncyc_annotation.pathways[0] ? `· top: <i>${escapeHtml(data.corncyc_annotation.pathways[0].common_name)}</i>` : ''}
            </span></div>` : ''}
        </div>
        ${topRowHTML}
    </section>`;
}

function renderChemicalEntityCard(data) {
    const ce = data.chemical_entity;
    if (!ce) return `<div class="data-card empty-inline">No chemical entity (EC-number query bypasses Phase 1).</div>`;
    return `<div class="data-card">
        <a href="${dbUrl.chebi(ce.chebi_id)}" target="_blank" class="badge badge-chebi">${ce.chebi_id}</a>
        ${ce.pubchem_cid ? `<a href="${dbUrl.pubchem(ce.pubchem_cid)}" target="_blank" class="badge" style="background:#0f766e;color:#fff;text-decoration:none;">CID: ${ce.pubchem_cid}</a>` : ''}
        <div style="margin-top:6px;font-size:0.88rem;"><strong>Monoisotopic mass:</strong> ${Number(ce.monoisotopic_mass||0).toFixed(4)}</div>
        ${ce.smiles ? `<div style="margin-top:4px;font-size:0.82rem;font-family:'JetBrains Mono',monospace;word-break:break-all;color:#94a3b8;">SMILES: ${ce.smiles}</div>` : ''}
    </div>`;
}

function renderReactionsBlock(reactions) {
    return reactions.map(r => {
        const rheaNum = r.rhea_id.replace('RHEA:', '');
        const eqHTML = r.equation
            ? `<div class="rxn-equation">${escapeHtml(r.equation)}</div>`
            : '';
        const ecHTML = (r.ec_numbers || []).length
            ? (r.ec_numbers || []).map(e => `<span class="badge badge-ec" title="Enzyme Commission #">EC ${e}</span>`).join('')
            : '';
        const pathwayChips = (r.pathway_names || []).slice(0, 3).map(p => `<span class="badge badge-pathway" title="KEGG pathway">${p}</span>`).join('');
        const morePathways = (r.pathway_names || []).length > 3 ? `<span class="exec-dim">+${r.pathway_names.length - 3} more</span>` : '';
        return `<div class="data-card rxn-card">
            <div class="rxn-head">
                <a href="${dbUrl.rhea(r.rhea_id)}" target="_blank" class="badge badge-rhea">${r.rhea_id}</a>
                <span class="badge badge-${r.is_transport ? 'transport' : 'metabolic'}">${r.is_transport ? 'Transport' : 'Metabolic'}</span>
                ${ecHTML}
            </div>
            ${eqHTML}
            ${pathwayChips ? `<div class="rxn-pathways">${pathwayChips} ${morePathways}</div>` : ''}
        </div>`;
    }).join('');
}

function renderProteinsBlock(proteins) {
    return proteins.map(p => {
        const goTags = (p.go_terms || []).slice(0, 6).map(go =>
            `<div class="go-item"><a href="${dbUrl.quickgo(go.id)}" target="_blank" class="badge badge-go">${go.id}</a> <span class="go-name">${escapeHtml(go.name)}</span></div>`
        ).join('');
        const moreGo = (p.go_terms || []).length > 6 ? `<div class="exec-dim" style="font-size:0.75rem;">+${p.go_terms.length - 6} more GO terms</div>` : '';
        const catColor = p.category === 'Enzyme' ? '#b91c1c' : p.category === 'Transporter' ? '#047857' : p.category === 'Receptor' ? '#92400e' : '#475569';
        return `<details class="protein-card data-card">
            <summary class="protein-summary">
                <a href="${dbUrl.uniprot(p.uniprot_accession)}" target="_blank" class="badge badge-uniprot">${p.uniprot_accession}</a>
                <span class="badge" style="background:${catColor};color:#fff;">${p.category}</span>
                <span class="exec-dim">${(p.go_terms || []).length} GO terms</span>
            </summary>
            <div class="go-list">${goTags}${moreGo}</div>
        </details>`;
    }).join('');
}

function renderComparaBlock(advTargets) {
    return advTargets.map(at => {
        if (!at.ensembl_orthologs || at.ensembl_orthologs.length === 0) {
            return `<details class="data-card">
                <summary class="protein-summary">
                    ${renderMaizeGeneLabel(at.maize_gene_model)}
                    <span class="exec-dim">no Compara orthologs</span>
                </summary>
            </details>`;
        }
        const rows = at.ensembl_orthologs.map(o => `
            <tr>
                <td><a href="${dbUrl.ensemblSpecies(o.species)}" target="_blank" class="ortho-link">${o.species}</a></td>
                <td><a href="${dbUrl.ensemblPlantsGene(o.species, o.gene_id)}" target="_blank" class="ortho-link mono">${o.gene_id}</a></td>
                <td class="mono-dim">${o.protein_id || '—'}</td>
                <td class="mono"><strong>${o.percent_identity.toFixed(1)}%</strong></td>
            </tr>`).join('');
        return `<details class="data-card">
            <summary class="protein-summary">
                ${renderMaizeGeneLabel(at.maize_gene_model)}
                <span class="exec-dim">${at.ensembl_orthologs.length} pan-plant ortholog${at.ensembl_orthologs.length===1?'':'s'}</span>
            </summary>
            <table class="compara-table">
                <thead><tr><th>Species</th><th>Gene</th><th>Protein</th><th>%ID</th></tr></thead>
                <tbody>${rows}</tbody>
            </table>
        </details>`;
    }).join('');
}

function renderCornCycBlock(ann, geneMetaMap) {
    if (!ann || !ann.pathways || ann.pathways.length === 0) {
        return `<div class="corncyc-empty">No CornCyc pathway annotation for this compound.</div>`;
    }
    const compoundsLine = (ann.compounds || []).map(c => `<a href="https://pmn.plantcyc.org/CORN/NEW-IMAGE?object=${encodeURIComponent(c.id)}" target="_blank" class="corncyc-link"><b>${escapeHtml(c.name || c.id)}</b></a>`).join(', ');
    const pathwayRows = ann.pathways.slice(0, 25).map(p => {
        const pwUrl = `https://pmn.plantcyc.org/pathway?orgid=CORN&id=${encodeURIComponent(p.id)}`;
        const geneChips = (p.maize_genes || []).slice(0, 6).map(g => {
            const m = (geneMetaMap || {})[g];
            const label = (m && m.symbol) ? `${g} <span class="cc-pw-sym">${escapeHtml(m.symbol)}</span>` : g;
            return `<a href="${dbUrl.maizeGene(g)}" target="_blank" class="cc-pw-gene-chip" title="${escapeHtmlAttr(g + (m && m.description ? ' — ' + m.description : ''))}">${label}</a>`;
        }).join('');
        const moreGenes = p.maize_genes.length > 6 ? `<span class="exec-dim"> +${p.maize_genes.length - 6} more</span>` : '';
        return `<div class="cc-pw-row">
            <div class="cc-pw-head">
                <a href="${pwUrl}" target="_blank" class="corncyc-link"><b>${p.common_name}</b></a>
                <span class="mono-dim">${p.id}</span>
                <span class="exec-dim">· ${p.reactions_touching_compound.length} matching reaction(s) · ${p.maize_genes.length} maize gene(s)</span>
            </div>
            <div class="cc-pw-genes">${geneChips}${moreGenes}</div>
        </div>`;
    }).join('');
    const overflow = ann.pathways.length > 25
        ? `<div class="exec-dim" style="margin-top:8px;">+${ann.pathways.length - 25} more pathway(s) — see CSV for the full list.</div>` : '';
    return `<div class="corncyc-block">
        <div class="corncyc-intro">
            Curated PlantCyc/CornCyc pathways involving ${compoundsLine}. Maize gene annotations from
            <a href="https://www.plantcyc.org/" target="_blank" class="corncyc-link">Plant Metabolic Network</a>
            (authors: Hawkins, Xue, Rhee; CornCyc ${ann.version || ''}).
            <span class="exec-dim">Full attribution: <code>CORNCYC_ATTRIBUTION.txt</code>.</span>
        </div>
        ${pathwayRows}
        ${overflow}
    </div>`;
}

function renderExecutionLogTable(logs) {
    const rows = logs.map(log => {
        const s = (log.status || '').toLowerCase();
        const color = s === 'success' ? '#4ade80' : s === 'error' ? '#f87171' : s === 'warning' ? '#fbbf24' : '#38bdf8';
        return `<tr>
            <td>P${log.phase}</td>
            <td><strong>${escapeHtml(log.database)}</strong></td>
            <td style="color:${color};">[${(log.status||'').toUpperCase()}]</td>
            <td class="mono">${log.hits}</td>
            <td>${escapeHtml(log.message)}</td>
        </tr>`;
    }).join('');
    return `<table class="execlog-table"><thead><tr><th>Phase</th><th>Database</th><th>Status</th><th>Hits</th><th>Message</th></tr></thead><tbody>${rows}</tbody></table>`;
}

function escapeHtml(s) {
    return String(s == null ? '' : s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
}

let _cyInstance = null;
// Cytoscape only lays out correctly when its container is visible. The graph
// lives inside a collapsed <details>; lay out fresh whenever the user opens it.
document.addEventListener('toggle', (e) => {
    if (e.target && e.target.id === 'network-graph-details' && e.target.open && _cyInstance) {
        try {
            _cyInstance.resize();
            _cyInstance.layout({ name: 'dagre', rankDir: 'LR', spacingFactor: 1.2, animate: true }).run();
        } catch (err) { console.warn('cy relayout failed', err); }
    }
}, true);

function renderCytoscape(data) {
    const cyEl = document.getElementById('cy');
    if (!cyEl || typeof cytoscape === 'undefined') return;

    // Trim payload to top-N of each phase so the graph stays readable.
    // Without this, dozens of reactions × dozens of proteins × dozens of orthologs
    // collapses into an unreadable hairball.
    const GRAPH_CAPS = { reactions: 8, proteins: 10, orthologs: 8 };
    const cappedReactions = (data.reactions || []).slice(0, GRAPH_CAPS.reactions);
    const cappedProteins  = (data.proteins  || []).slice(0, GRAPH_CAPS.proteins);
    // For orthologs, prefer consensus first, then the rest, capped.
    const orthRows = buildEnrichmentRows(data.orthologs || [], data);
    const cappedOrthologs = orthRows.slice(0, GRAPH_CAPS.orthologs).map(r => ({
        maize_gene_model: r.gene,
        query_uniprot_id: r.query_uniprot_id,
    }));
    data = { ...data, reactions: cappedReactions, proteins: cappedProteins, orthologs: cappedOrthologs };

    const elements = [];
    const queryType = data.input_data.query_type;
    let rootId = null;

    if (queryType === 'mz' && data.input_data.mz != null) {
        const mzId = data.input_data.mz.toString();
        elements.push({ data: { id: mzId, label: 'm/z ' + mzId, type: 'mz' } });
        if (data.chemical_entity) {
            const chebiId = data.chemical_entity.chebi_id;
            elements.push({ data: { id: chebiId, label: chebiId, type: 'compound' } });
            elements.push({ data: { id: 'e1', source: mzId, target: chebiId } });
            rootId = chebiId;
        }
    } else if (queryType === 'chemical' && data.chemical_entity) {
        const chebiId = data.chemical_entity.chebi_id;
        elements.push({ data: { id: chebiId, label: chebiId, type: 'compound' } });
        rootId = chebiId;
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

    if (_cyInstance) { try { _cyInstance.destroy(); } catch(e){} }
    const cy = cytoscape({
        container: cyEl,
        elements: elements,
        style: [
            { selector: 'node', style: {
                'label': 'data(label)', 'color': '#fff', 'text-outline-color': '#0f172a',
                'text-outline-width': 2, 'text-valign': 'top', 'text-halign': 'center',
                'font-size': '12px', 'cursor': 'pointer'
            }},
            { selector: 'node[type="mz"]',       style: { 'background-color': '#f87171', 'shape': 'diamond' } },
            { selector: 'node[type="compound"]', style: { 'background-color': '#818cf8', 'shape': 'ellipse' } },
            { selector: 'node[type="reaction"]', style: { 'background-color': '#34d399', 'shape': 'round-rectangle' } },
            { selector: 'node[type="protein"]',  style: { 'background-color': '#fbbf24', 'shape': 'hexagon' } },
            { selector: 'node[type="gene"]',     style: { 'background-color': '#c084fc', 'shape': 'triangle' } },
            { selector: 'edge', style: {
                'width': 2, 'line-color': '#475569', 'target-arrow-color': '#475569',
                'target-arrow-shape': 'triangle', 'curve-style': 'bezier', 'arrow-scale': 0.8
            }},
        ],
        layout: { name: 'dagre', rankDir: 'LR', spacingFactor: 1.2, animate: true },
    });

    cy.on('tap', 'node', function(evt) {
        const node = evt.target;
        const type = node.data('type');
        const id = node.id();
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
    _cyInstance = cy;
}

function renderLaneSection(lane, rows, blurb) {
    const meta = laneMeta[lane];
    if (!rows || rows.length === 0) {
        return `<details class="lane-block ${meta.cls}">
            <summary class="lane-header">
                <span class="drill-caret"></span>
                <span class="lane-badge">${meta.icon} ${meta.label}</span>
                <span class="lane-count">0</span>
            </summary>
            <div class="lane-blurb">${blurb}</div>
            <div class="data-card empty-inline" style="text-align:center;">No hits in this lane.</div>
        </details>`;
    }

    // Auto-open consensus lane; keep others collapsed.
    const open = (lane === 'consensus') ? ' open' : '';
    const cards = rows.map(r => renderLaneRow(r, lane)).join('');

    return `<details class="lane-block ${meta.cls}"${open}>
        <summary class="lane-header">
            <span class="drill-caret"></span>
            <span class="lane-badge">${meta.icon} ${meta.label}</span>
            <span class="lane-count">${rows.length}</span>
            ${renderLaneTopHint(rows)}
        </summary>
        <div class="lane-blurb">${blurb}</div>
        ${cards}
    </details>`;
}

function renderLaneTopHint(rows) {
    const top = rows[0];
    if (!top) return '';
    const bits = [];
    if (top.tm != null)      bits.push(`TM ${top.tm.toFixed(2)}`);
    if (top.seq_sim != null) bits.push(`${top.seq_sim.toFixed(1)}% seq`);
    if (top.plddt != null)   bits.push(`pLDDT ${top.plddt.toFixed(1)}`);
    return `<span class="lane-top-hint">top: ${top.gene}${bits.length ? ` (${bits.join(' · ')})` : ''}</span>`;
}

/**
 * For consensus hits, return a short subtype tag like "seq + cur" or
 * "all 3 lanes" so the UI can show *which* lanes agreed.
 */
function consensusSubtype(sources) {
    const s = new Set(sources || []);
    const seq    = [...s].some(x => x !== STRUCTURAL_SOURCE && x !== CURATED_SOURCE);
    const struct = s.has(STRUCTURAL_SOURCE);
    const cur    = s.has(CURATED_SOURCE);
    if (seq && struct && cur) return 'all 3 lanes';
    if (seq && struct)        return 'seq + struct';
    if (seq && cur)           return 'seq + curated';
    if (struct && cur)        return 'struct + curated';
    return '';
}

function renderLaneRow(r, lane) {
    const seq = (r.seq_sim != null) ? r.seq_sim.toFixed(1) + '%' : '—';
    // For sequence-only hits that got the cheap path, TM is unknown (skipped Foldseek)
    const tm  = (r.tm != null && (r.tm > 0 || (r.sources || []).includes(STRUCTURAL_SOURCE))) ? r.tm.toFixed(2) : '—';
    const pl  = (r.plddt  != null && r.plddt > 0) ? r.plddt.toFixed(1) : '—';
    const exprCount = r.expr_count || 0;
    const exprLabel = exprCount ? `${exprCount}` : '0';

    const queryProteinHtml = r.query_uniprot_id
        ? `<a href="${dbUrl.uniprot(r.query_uniprot_id)}" target="_blank" class="mono">${r.query_uniprot_id}</a>`
        : '—';

    const sourcesHtml = (r.sources || []).map(s => {
        const cls = s === STRUCTURAL_SOURCE ? 'badge-src-struct'
                  : s === CURATED_SOURCE    ? 'badge-src-curated'
                  : 'badge-src-seq';
        return `<span class="badge badge-src ${cls}">${s}</span>`;
    }).join('');

    const exprDetail = exprCount > 0
        ? '<div class="tissue-list">' + (r.expr_experiments || []).map(e =>
              `<a href="${dbUrl.gxaExperiment(e)}" target="_blank" class="tissue-pill" title="Open ${e} on Expression Atlas">${e}</a>`
          ).join('') + (exprCount > (r.expr_experiments || []).length ? ` <span class="exec-dim">+${exprCount - r.expr_experiments.length} more</span>` : '') + '</div>'
        : '<span class="empty-inline">no expression evidence</span>';

    const domainHTML = r.domain
        ? `<div class="meta-line"><span class="meta-label">Pfam:</span>
               <a href="${dbUrl.pfam(r.domain.pfam_domain_id)}" target="_blank">${r.domain.pfam_domain_id}</a>
               ${escapeHtml(r.domain.pfam_domain_name)}
               <span class="mono-dim">(${r.domain.domain_start}–${r.domain.domain_end}, TM=${r.domain.domain_tm_score.toFixed(2)})</span>
           </div>`
        : '';

    const comparaHTML = r.top_compara
        ? `<div class="meta-line"><span class="meta-label">Compara:</span>
               <a href="${dbUrl.ensemblPlantsGene(r.top_compara.species, r.top_compara.gene_id)}" target="_blank" class="ortho-link">${r.top_compara.species}:${r.top_compara.gene_id}</a>
               <span class="mono-dim">(${r.top_compara.percent_identity.toFixed(1)}%${r.num_compara > 1 ? `, +${r.num_compara-1} more` : ''})</span>
           </div>`
        : '';

    const enrichedBadge = r.enrichment_kind === 'full'
        ? '<span class="badge badge-enriched" title="Top-N: full Phase 5 enrichment ran (pLDDT + expression + 1-to-1 Foldseek where needed)">full</span>'
        : (r.enrichment_kind === 'cheap'
            ? '<span class="badge badge-cheap" title="Cheap enrichment ran (pLDDT + Gramene expression breadth). 1-to-1 Foldseek not needed because TM came from Phase 4.5, or was skipped for non-top-N sequence-only hits.">light</span>'
            : '<span class="badge badge-unenriched" title="Filtered out by pLDDT/expression threshold — click Run Foldseek below to retry">filtered</span>');

    const subtype = (lane === 'consensus') ? consensusSubtype(r.sources) : '';
    const subtypeBadge = subtype ? `<span class="lane-row-subtype">${subtype}</span>` : '';

    // ONE-LINE summary in <summary>; details expand to full meta
    return `<details class="data-card lane-row lane-row-${lane}">
        <summary class="lane-row-summary">
            <span class="drill-caret-tiny"></span>
            ${renderMaizeGeneLabel(r.gene)}
            ${subtypeBadge}
            <span class="lane-row-metric">TM <strong>${tm}</strong></span>
            <span class="lane-row-metric">seq <strong>${seq}</strong></span>
            <span class="lane-row-metric">pLDDT <strong>${pl}</strong></span>
            <span class="lane-row-metric">expr <strong>${exprLabel}</strong></span>
            ${enrichedBadge}
        </summary>
        <div class="lane-row-body">
            <div class="lane-row-meta">
                <div class="meta-line"><span class="meta-label">Sources:</span> ${sourcesHtml || '—'}</div>
                <div class="meta-line"><span class="meta-label">Query protein:</span> ${queryProteinHtml}</div>
                ${domainHTML}
                ${comparaHTML}
                <div class="meta-line"><span class="meta-label">Expression Atlas:</span> ${exprDetail}</div>
            </div>
        </div>
    </details>`;
}

async function validateTarget(btnElement, targetMapping) {
    btnElement.innerHTML = 'Running Foldseek...';
    btnElement.disabled = true;
    try {
        const data = await validateOne(targetMapping);
        btnElement.parentElement.remove();
        ingestValidatedTarget(targetMapping, data);
        alert("Target Validated! Gene: " + data.target.maize_gene_model + " (TM-Score: " + data.target.tm_score + ")");
    } catch (error) {
        console.error('Error during validation:', error);
        btnElement.innerHTML = 'Error! Try Again';
        btnElement.disabled = false;
    }
}

// ---------- Bulk Validation ----------

// Single validation HTTP call — returns the parsed response from /api/validate_target.
async function validateOne(targetMapping) {
    const response = await fetch('/api/validate_target', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(targetMapping),
    });
    if (!response.ok) {
        const txt = await response.text().catch(() => '');
        throw new Error(`HTTP ${response.status} ${txt.slice(0, 200)}`);
    }
    return response.json();
}

// Merge a validated target back into `currentSingleResult` so downloaded
// reports include freshly validated entries.
function ingestValidatedTarget(originalMapping, validateResponse) {
    if (!currentSingleResult) return;
    currentSingleResult.targets = currentSingleResult.targets || [];
    currentSingleResult.domain_targets = currentSingleResult.domain_targets || [];
    currentSingleResult.advanced_homology_targets = currentSingleResult.advanced_homology_targets || [];

    if (validateResponse.target) {
        const gene = validateResponse.target.maize_gene_model;
        if (!currentSingleResult.targets.some(t => t.maize_gene_model === gene)) {
            currentSingleResult.targets.push(validateResponse.target);
        }
        if (validateResponse.domain_target &&
            !currentSingleResult.domain_targets.some(d => d.maize_gene_model === gene)) {
            currentSingleResult.domain_targets.push(validateResponse.domain_target);
        }
        if (validateResponse.advanced_homology_target &&
            !currentSingleResult.advanced_homology_targets.some(a => a.maize_gene_model === gene)) {
            currentSingleResult.advanced_homology_targets.push(validateResponse.advanced_homology_target);
        }
    }
    // Drop the validated gene from the unvalidated pool
    currentSingleResult.unvalidated_targets = (currentSingleResult.unvalidated_targets || [])
        .filter(t => t.maize_gene_model !== originalMapping.maize_gene_model);

    // Mirror execution logs into the terminal
    (validateResponse.execution_logs || []).forEach(e => term.log(e));
}

function renderUnvalidatedSection(unvalidated) {
    const unvalContainer = document.getElementById('unvalidated-container');
    const unvalList      = document.getElementById('unvalidated-list');
    const countBadge     = document.getElementById('unval-count-badge');
    const progressEl     = document.getElementById('bulk-progress');
    const numInput       = document.getElementById('bulk-validate-n');
    const allBtn         = document.getElementById('bulk-validate-all-btn');
    const topBtn         = document.getElementById('bulk-validate-n-btn');

    if (!unvalContainer || !unvalList) return;

    if (!unvalidated || unvalidated.length === 0) {
        unvalContainer.classList.add('hidden');
        return;
    }

    unvalContainer.classList.remove('hidden');
    if (countBadge) countBadge.textContent = `${unvalidated.length} pending`;
    if (progressEl) progressEl.textContent = '';
    if (numInput) {
        numInput.max = unvalidated.length;
        if (parseInt(numInput.value, 10) > unvalidated.length) {
            numInput.value = unvalidated.length;
        }
    }

    unvalList.innerHTML = '';
    unvalidated.forEach((t, i) => {
        const card = document.createElement('div');
        card.className = 'data-card unval-row queued';
        card.dataset.gene = t.maize_gene_model;
        card.dataset.index = String(i);
        card.innerHTML = `
            <div style="flex:1; min-width:0;">
                <div class="unval-gene-line">${renderMaizeGeneLabel(t.maize_gene_model)}</div>
                <div style="font-size: 0.85rem; color: #94a3b8; margin-top: 4px;">
                    Source: ${t.sources.join(', ')} · Similarity: ${t.similarity_score.toFixed(1)}
                    ${t.query_uniprot_id ? ` · Query: <a href="${dbUrl.uniprot(t.query_uniprot_id)}" target="_blank" class="ortho-link">${t.query_uniprot_id}</a>` : ''}
                    ${t.plaza_orthogroup ? `<br>Orthogroup: <span style="color: #4ade80;">${t.plaza_orthogroup}</span>` : ''}
                </div>
                <div class="unval-status" data-status>queued</div>
            </div>
            <button type="button" class="btn-glow" style="padding: 6px 12px; font-size: 0.85rem; width:auto;"
                onclick='validateTarget(this, ${JSON.stringify(t).replace(/'/g, "&#39;")})'>Run Foldseek</button>
        `;
        unvalList.appendChild(card);
    });

    // Wire bulk buttons (idempotent — replace handlers each time we render)
    if (allBtn) {
        allBtn.onclick = () => runBulkValidation(unvalidated.length);
    }
    if (topBtn) {
        topBtn.onclick = () => {
            const n = parseInt(numInput?.value, 10);
            if (!n || n < 1) return;
            runBulkValidation(Math.min(n, unvalidated.length));
        };
    }
}

async function runBulkValidation(count) {
    if (!currentSingleResult) return;
    const pending = (currentSingleResult.unvalidated_targets || []).slice(0, count);
    if (pending.length === 0) return;

    const allBtn   = document.getElementById('bulk-validate-all-btn');
    const topBtn   = document.getElementById('bulk-validate-n-btn');
    const progress = document.getElementById('bulk-progress');
    [allBtn, topBtn].forEach(b => { if (b) b.disabled = true; });

    term.show();
    term.info(`Bulk validation: running Foldseek on top ${pending.length} of ${currentSingleResult.unvalidated_targets.length} unvalidated targets`);

    // Mark all queued rows visually (only those we'll touch)
    pending.forEach(t => updateRowStatus(t.maize_gene_model, 'queued', 'queued'));

    let success = 0, failed = 0;
    for (let i = 0; i < pending.length; i++) {
        const t = pending[i];
        if (progress) progress.textContent = `${i + 1}/${pending.length} · ${success} ok · ${failed} fail`;
        updateRowStatus(t.maize_gene_model, 'running', 'running Foldseek...');
        try {
            const resp = await validateOne(t);
            ingestValidatedTarget(t, resp);
            const tm = resp.target?.tm_score?.toFixed(2);
            const pl = resp.target?.plddt?.toFixed(1);
            updateRowStatus(t.maize_gene_model, 'done', `done · TM=${tm} pLDDT=${pl}`);
            success++;
        } catch (err) {
            console.error('Bulk validation failed for', t.maize_gene_model, err);
            term.error(`Validation failed for ${t.maize_gene_model}: ${err.message}`);
            updateRowStatus(t.maize_gene_model, 'failed', `failed: ${err.message}`);
            failed++;
        }
    }

    if (progress) progress.textContent = `Done · ${success} ok · ${failed} fail`;
    term.success(`Bulk validation finished: ${success} succeeded, ${failed} failed.`);

    // Re-render the dashboard with freshly merged data so Phase 5/6/7 sections grow
    renderResults(currentSingleResult);

    [allBtn, topBtn].forEach(b => { if (b) b.disabled = false; });
}

function updateRowStatus(geneModel, cssState, label) {
    const row = document.querySelector(`.unval-row[data-gene="${CSS.escape(geneModel)}"]`);
    if (!row) return;
    row.classList.remove('queued', 'running', 'done', 'failed');
    row.classList.add(cssState);
    const statusEl = row.querySelector('[data-status]');
    if (statusEl) {
        statusEl.classList.remove('s-running', 's-done', 's-failed');
        if (cssState === 'running') statusEl.classList.add('s-running');
        else if (cssState === 'done') statusEl.classList.add('s-done');
        else if (cssState === 'failed') statusEl.classList.add('s-failed');
        statusEl.textContent = label;
    }
    // If done, disable that row's own button
    if (cssState === 'done' || cssState === 'running') {
        const btn = row.querySelector('button');
        if (btn) btn.disabled = true;
    }
}
