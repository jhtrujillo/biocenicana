package org.cenicana.bio.io;

import org.cenicana.bio.core.GeneFeature;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.*;

/**
 * Optimized interactive HTML dashboard for exploring gene annotations near markers.
 */
public class AnnotationDashboardGenerator {

    public static class MarkerMatch {
        public String markerId;
        public String chromosome;
        public long position;
        public String refAllele;
        public String altAllele;
        public double alleleFreq;
        public String[] genotypes;
        public String effect; // e.g. "Missense (V34I)"
        public String source; // Primary or Secondary
        public int kaspScore;
        public String flankingSeq;
        public GeneFeature gene;
        public long distance;

        public MarkerMatch(String markerId, String chromosome, long position, GeneFeature gene) {
            this.markerId = markerId;
            this.chromosome = chromosome;
            this.position = position;
            this.gene = gene;
            this.distance = gene.distanceTo(position);
        }
    }

    private List<String> samples = new ArrayList<>();
    public void setSamples(List<String> samples) { this.samples = samples; }

    private Set<String> uniqueFunctions = new TreeSet<>();
    public void setUniqueFunctions(Set<String> functions) { this.uniqueFunctions = functions; }

    public void generate(List<MarkerMatch> matches, String outputPath) throws IOException {
        Map<String, List<MarkerMatch>> grouped = new java.util.LinkedHashMap<>();
        for (MarkerMatch m : matches) {
            grouped.computeIfAbsent(m.markerId, k -> new ArrayList<>()).add(m);
        }

        try (PrintWriter pw = new PrintWriter(new java.io.FileWriter(outputPath))) {
            pw.println("<!DOCTYPE html>");
            pw.println("<html lang='en'>");
            pw.println("<head>");
            pw.println("    <meta charset='UTF-8'>");
            pw.println("    <title>BioJava Genomic Suite - Functional Dashboard</title>");
            pw.println("    <link rel='stylesheet' href='https://cdn.datatables.net/1.13.4/css/jquery.dataTables.css'>");
            pw.println("    <link rel='stylesheet' href='https://cdn.datatables.net/buttons/2.3.6/css/buttons.dataTables.min.css'>");
            pw.println("    <link href='https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&display=swap' rel='stylesheet'>");
            pw.println("    <style>");
            pw.println("        :root { --primary: #2563eb; --accent: #3b82f6; --bg: #f1f5f9; --card: #ffffff; }");
            pw.println("        body { font-family: 'Inter', sans-serif; background: var(--bg); margin: 0; display: flex; height: 100vh; overflow: hidden; }");
            pw.println("        .sidebar { width: 350px; background: var(--card); border-right: 1px solid #e2e8f0; display: flex; flex-direction: column; }");
            pw.println("        .sidebar-header { padding: 20px; border-bottom: 1px solid #e2e8f0; }");
            pw.println("        .search-box { width: 100%; padding: 10px; border: 1px solid #cbd5e1; border-radius: 8px; box-sizing: border-box; }");
            pw.println("        .marker-list { flex: 1; overflow-y: auto; }");
            pw.println("        .marker-item { padding: 15px 20px; border-bottom: 1px solid #f1f5f9; cursor: pointer; transition: 0.2s; }");
            pw.println("        .marker-item:hover { background: #f8fafc; }");
            pw.println("        .marker-item.active { background: #eff6ff; border-left: 4px solid var(--primary); }");
            pw.println("        .main-content { flex: 1; padding: 30px; overflow-y: auto; display: flex; flex-direction: column; }");
            pw.println("        .card { background: var(--card); border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); padding: 25px; margin-bottom: 25px; }");
            pw.println("        .badge { display: inline-block; padding: 4px 10px; border-radius: 6px; font-size: 0.75rem; font-weight: 600; color: white; margin: 2px; border: none; cursor: pointer; }");
            pw.println("        .bg-mis { background: #f59e0b; } .bg-syn { background: #10b981; } .bg-stop { background: #ef4444; } .bg-intron { background: #64748b; }");
            pw.println("        .bg-go { background: #6366f1; } .bg-pfam { background: #ec4899; } .bg-interpro { background: #8b5cf6; }");
            pw.println("        .gene-map-svg { width: 100%; height: 60px; margin-top: 10px; }");
            pw.println("        .modal { display:none; position:fixed; top:0; left:0; width:100%; height:100%; background:rgba(0,0,0,0.5); z-index:1000; align-items:center; justify-content:center; }");
            pw.println("        .modal-content { background:white; padding:30px; border-radius:15px; width:80%; max-height:80%; overflow-y:auto; position:relative; }");
            pw.println("        .close-btn { position:absolute; top:15px; right:20px; cursor:pointer; font-size:1.5rem; color:#94a3b8; }");
            pw.println("        #splash { position:fixed; top:0; left:0; width:100%; height:100%; background:white; z-index:2000; display:flex; align-items:center; justify-content:center; flex-direction:column; }");
            pw.println("    </style>");
            pw.println("</head>");
            pw.println("<body>");
            pw.println("    <div id='splash'><h1>BioJava Genomic Suite</h1><p>Indexing candidates and population data...</p></div>");
            
            pw.println("    <div class='sidebar'>");
            pw.println("        <div class='sidebar-header'>");
            pw.println("            <h2 style='margin:0 0 5px 0; font-size:1.4rem'>Markers</h2>");
            pw.println("            <input type='text' id='markerSearch' list='functionList' class='search-box' placeholder='Type to search function...'>");
            
            pw.println("            <datalist id='functionList'>");
            for (String func : uniqueFunctions) {
                pw.println("                <option value=\"" + func.replace("\"", "&quot;") + "\">");
            }
            pw.println("            </datalist>");

            pw.println("            <div style='margin-top:10px'>");
            pw.println("                <label style='font-size:0.7rem; color:#64748b; font-weight:600; text-transform:uppercase'>Browse Catalog:</label>");
            pw.println("                <select id='functionDropdown' class='search-box' style='margin-top:5px; background-image:none'>");
            pw.println("                    <option value=''>-- Full Function Index --</option>");
            for (String func : uniqueFunctions) {
                pw.println("                    <option value=\"" + func.replace("\"", "&quot;") + "\">" + func + "</option>");
            }
            pw.println("                </select>");
            pw.println("            </div>");

            pw.println("            <div style='margin-top:10px; display:flex; gap:10px; align-items:center'>");
            pw.println("                <label style='font-size:0.7rem; color:#64748b; font-weight:600'>MIN QUALITY:</label>");
            pw.println("                <select id='starFilter' class='search-box' style='padding:5px; width:auto; border-color:var(--primary)' onchange='renderSidebar($(\"#markerSearch\").val())'>");
            pw.println("                    <option value='0'>Any Quality</option>");
            pw.println("                    <option value='50'>⭐⭐+ (Fair)</option>");
            pw.println("                    <option value='80'>⭐⭐⭐ (Best)</option>");
            pw.println("                </select>");
            pw.println("            </div>");

            pw.println("            <div style='margin-top:10px; display:flex; gap:10px'>");
            pw.println("                <button class='badge' style='background:#f59e0b; border:none; cursor:pointer' onclick='filterHighImpact()'>High Impact</button>");
            pw.println("                <button class='badge' style='background:#6366f1; border:none; cursor:pointer' onclick='downloadAllFASTA()'>FASTA (.faa)</button>");
            pw.println("                <button class='badge' style='background:#10b981; border:none; cursor:pointer' onclick='showHaplotypeMatrix()'>Haplotypes</button>");
            pw.println("            </div>");
            pw.println("        </div>");
            pw.println("        <div class='marker-list' id='markerList'></div>");
            pw.println("    </div>");

            pw.println("    <div class='main-content'>");
            pw.println("        <div id='emptyState' style='text-align:center; margin:auto; color:#94a3b8'><h2>Select a marker or view Haplotype Matrix</h2></div>");
            
            // Haplotype Matrix Container
            pw.println("        <div id='haplotypeView' style='display:none'>");
            pw.println("            <div class='card'>");
            pw.println("                <h2>Population Haplotype Matrix</h2>");
            pw.println("                <p style='color:#64748b; font-size:0.8rem'>Visualizing genotypes across all samples for filtered markers.</p>");
            pw.println("                <div id='matrixLegend' style='margin-bottom:15px; font-size:0.7rem'>");
            pw.println("                    <span style='margin-right:15px'><span style='color:#10b981'>■</span> 0/0 (Ref)</span>");
            pw.println("                    <span style='margin-right:15px'><span style='color:#f59e0b'>■</span> 0/1 (Het)</span>");
            pw.println("                    <span style='margin-right:15px'><span style='color:#ef4444'>■</span> 1/1 (Alt)</span>");
            pw.println("                </div>");
            pw.println("                <div id='haploContainer' style='overflow:auto; max-height:70vh; border:1px solid #e2e8f0; border-radius:8px'>");
            pw.println("                    <table id='haploTable' style='border-collapse:collapse; font-size:0.6rem; width:max-content'>");
            pw.println("                        <thead id='haploHead' style='position:sticky; top:0; background:white; z-index:10'></thead>");
            pw.println("                        <tbody id='haploBody'></tbody>");
            pw.println("                    </table>");
            pw.println("                </div>");
            pw.println("            </div>");
            pw.println("        </div>");

            pw.println("        <div id='detailView' style='display:none'>");
            pw.println("            <div class='card' style='padding:15px 30px; margin-bottom:20px; display:flex; gap:20px; align-items:center; background:#f8fafc'>");
            pw.println("                <span style='font-weight:600; font-size:0.8rem; color:#64748b; text-transform:uppercase'>Legend:</span>");
            pw.println("                <span class='badge bg-mis'>Missense</span> <span class='badge bg-syn'>Synonymous</span>");
            pw.println("                <span class='badge bg-stop'>Stop-Gain</span> <span class='badge bg-intron'>Intronic</span>");
            pw.println("            </div>");
            pw.println("            <div class='card'>");
            pw.println("                <div style='display:flex; justify-content:space-between; align-items:center; margin-bottom:20px'>");
            pw.println("                    <h2 id='selectedMarkerTitle' style='margin:0'></h2>");
            pw.println("                    <div id='markerPosBadge' class='badge' style='background:#f1f5f9; color:#64748b; padding:8px 15px'></div>");
            pw.println("                </div>");
            pw.println("                <table id='geneTable' class='display' style='width:100%'>");
            pw.println("                    <thead><tr><th>Gene (Source) & Structure</th><th>Effect & Freq</th><th>Dist</th><th>KASP Rec.</th><th>Annotation</th><th>Action</th></tr></thead>");
            pw.println("                    <tbody id='geneTableBody'></tbody>");
            pw.println("                </table>");
            pw.println("            </div>");
            pw.println("        </div>");
            pw.println("    </div>");

            pw.println("    <div id='seqModal' class='modal' onclick='$(this).hide()'>");
            pw.println("        <div class='modal-content' onclick='event.stopPropagation()'>");
            pw.println("            <span class='close-btn' onclick='$(\"#seqModal\").hide()'>&times;</span>");
            pw.println("            <h3 id='seqTitle'></h3>");
            pw.println("            <pre id='seqText' style='background:#f8fafc; padding:20px; border-radius:10px; word-break:break-all; white-space:pre-wrap'></pre>");
            pw.println("            <button class='badge' style='background:var(--primary); border:none; cursor:pointer' onclick='copySeq()'>Copy Sequence</button>");
            pw.println("        </div>");
            pw.println("    </div>");

            pw.println("    <script src='https://code.jquery.com/jquery-3.6.0.min.js'></script>");
            pw.println("    <script src='https://cdn.datatables.net/1.13.4/js/jquery.dataTables.js'></script>");
            pw.println("    <script src='https://cdn.datatables.net/buttons/2.3.6/js/dataTables.buttons.min.js'></script>");
            pw.println("    <script src='https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js'></script>");
            pw.println("    <script src='https://cdn.datatables.net/buttons/2.3.6/js/buttons.html5.min.js'></script>");
            
            pw.println("    <script>");
            pw.print("        const samples = [");
            for (String s : samples) pw.print("'" + s + "',");
            pw.println("];");

            pw.println("        const markerIds = [" + String.join(",", grouped.keySet().stream().map(s -> "'" + s + "'").toArray(String[]::new)) + "];");
            pw.println("        const data = {");
            for (Map.Entry<String, List<MarkerMatch>> entry : grouped.entrySet()) {
                pw.print("            '" + entry.getKey() + "': [");
                for (MarkerMatch m : entry.getValue()) {
                    pw.print("{ id: '" + m.gene.getId() + "', gStart: " + m.gene.getStart() + ", gEnd: " + m.gene.getEnd() + ", strand: '" + m.gene.getStrand() + "', ");
                    pw.print("pos: " + m.position + ", d: " + m.distance + ", n: '" + (m.gene.getNote() != null ? m.gene.getNote().replace("'", "\\'").replace("\n", " ") : "-") + "', ");
                    pw.print("ef: '" + (m.effect != null ? m.effect : "N/A") + "', af: " + m.alleleFreq + ", src: '" + m.source + "', ");
                    pw.print("ra: '" + m.refAllele + "', aa: '" + m.altAllele + "', gt: [" + (m.genotypes != null ? String.join(",", Arrays.stream(m.genotypes).map(s -> "'" + s + "'").toArray(String[]::new)) : "") + "], ");
                    pw.print("kasp: " + m.kaspScore + ", fl: '" + (m.flankingSeq != null ? m.flankingSeq : "") + "', ");
                    pw.print("go: " + toJsonArray(m.gene.getGoTerms()) + ", pf: " + toJsonArray(m.gene.getPfamDomains()) + ", ip: " + toJsonArray(m.gene.getInterproDomains()) + ", ");
                    pw.print("pr: '" + (m.gene.getProteinSequence() != null ? m.gene.getProteinSequence().replace("\n", "") : "") + "', cd: '" + (m.gene.getCdsSequence() != null ? m.gene.getCdsSequence().replace("\n", "") : "") + "', ");
                    pw.print("sb: [");
                    for (GeneFeature.SubFeature sub : m.gene.getSubFeatures()) pw.print("{t:'" + sub.type + "',s:" + sub.start + ",e:" + sub.end + "},");
                    pw.print("] },");
                }
                pw.println("],");
            }
            pw.println("        };");

            pw.println("        function showKasp(mId, gIdx) {");
            pw.println("            const g = data[mId][gIdx];");
            pw.println("            $('#seqTitle').text('KASP Flanking Sequence - ' + mId);");
            pw.println("            $('#seqText').text(g.fl);");
            pw.println("            $('#seqModal').css('display', 'flex');");
            pw.println("        }");

            pw.println("        let table; let isHighImpact = false;");
            pw.println("        $(document).ready(() => { renderSidebar(); $('#splash').fadeOut(); });");

            pw.println("        function renderSidebar(filter = '') {");
            pw.println("            const container = $('#markerList').empty();");
            pw.println("            let count = 0; const val = filter.toLowerCase();");
            pw.println("            const minKasp = parseInt($('#starFilter').val() || '0');");
            pw.println("            if (val.length > 2) {");
            pw.println("                container.append(`<div style='padding:10px 25px; background:#eef2ff; color:var(--primary); font-size:0.8rem; font-weight:600; cursor:pointer; border-radius:10px; margin:10px' onclick='showGlobalResults(\"${val}\")'>🔍 View all matches &rarr;</div>`);");
            pw.println("            }");
            pw.println("            for (const id of markerIds) {");
            pw.println("                let matchType = ''; let matches = id.toLowerCase().includes(val);");
            pw.println("                if (matches) matchType = 'Marker ID';");
            pw.println("                if (!matches && data[id]) {");
            pw.println("                    const foundGene = data[id].find(g => g.n.toLowerCase().includes(val) || g.go.some(t => t.toLowerCase().includes(val)));");
            pw.println("                    if (foundGene) { matches = true; matchType = foundGene.n.length > 30 ? foundGene.n.substring(0,30) + '...' : foundGene.n; }");
            pw.println("                }");
            pw.println("                if (isHighImpact && data[id]) matches = matches && data[id].some(g => g.ef.includes('Missense') || g.ef.includes('Stop-Gain'));");
            pw.println("                let bestKasp = data[id] ? Math.max(...data[id].map(g => g.kasp)) : 0;");
            pw.println("                if (bestKasp < minKasp) matches = false;");
            pw.println("                if (matches) {");
            pw.println("                    container.append(`<div class='marker-item' onclick='showMarker(\"${id}\")' data-id='${id}'><div style='font-weight:600'>${id}</div><div style='font-size:0.7rem; color:#64748b'>${data[id].length} genes nearby</div>${val.length > 0 ? `<div style='font-size:0.65rem; color:var(--primary); margin-top:3px'>Matches: ${matchType}</div>` : ''}</div>`);");
            pw.println("                    if (++count >= 1000) break;");
            pw.println("                }");
            pw.println("            }");
            pw.println("        }");

            pw.println("        $('#markerSearch').on('input', function() { renderSidebar($(this).val()); });");
            pw.println("        $('#functionDropdown').on('change', function() { const val = $(this).val(); if (val) { $('#markerSearch').val(val); renderSidebar(val); showGlobalResults(val); } });");
            pw.println("        function filterHighImpact() { isHighImpact = !isHighImpact; renderSidebar($('#markerSearch').val()); }");

            pw.println("        function showHaplotypeMatrix() {");
            pw.println("            $('#emptyState, #detailView').hide(); $('#haplotypeView').show();");
            pw.println("            const head = $('#haploHead').empty();");
            pw.println("            const body = $('#haploBody').empty();");
            pw.println("            let headRow = '<tr><th style=\"padding:5px; border:1px solid #e2e8f0\">Marker / Sample</th>';");
            pw.println("            samples.forEach(s => headRow += `<th style=\"padding:5px; border:1px solid #e2e8f0; writing-mode:vertical-rl; text-orientation:mixed\">${s}</th>`);");
            pw.println("            head.append(headRow + '</tr>');");
            pw.println("            const query = $('#markerSearch').val().toLowerCase();");
            pw.println("            const minKasp = parseInt($('#starFilter').val() || '0');");
            pw.println("            markerIds.forEach(id => {");
            pw.println("                const markerData = data[id];");
            pw.println("                if (!markerData) return;");
            pw.println("                let matches = id.toLowerCase().includes(query) || markerData.some(g => g.n.toLowerCase().includes(query) || g.go.some(t => t.toLowerCase().includes(query)));");
            pw.println("                if (isHighImpact) matches = matches && markerData.some(g => g.ef.includes('Missense') || g.ef.includes('Stop-Gain'));");
            pw.println("                let bestKasp = Math.max(...markerData.map(g => g.kasp));");
            pw.println("                if (bestKasp < minKasp) matches = false;");
            pw.println("                if (matches) {");
            pw.println("                    let row = `<tr><td style=\"padding:5px; border:1px solid #e2e8f0; cursor:pointer; font-weight:bold\" onclick=\"showMarker('${id}')\">${id}</td>`;");
            pw.println("                    const g = markerData[0];");
            pw.println("                    samples.forEach((s, i) => {");
            pw.println("                        let gt = g.gt[i] || './.';");
            pw.println("                        let color = gt.includes('1/1') ? '#ef4444' : (gt.includes('0/1') ? '#f59e0b' : (gt.includes('0/0') ? '#10b981' : '#f1f5f9'));");
            pw.println("                        row += `<td title=\"${s}: ${gt}\" style=\"background:${color}; width:15px; border:1px solid rgba(255,255,255,0.2)\"></td>`;");
            pw.println("                    });");
            pw.println("                    body.append(row + '</tr>');");
            pw.println("                }");
            pw.println("            });");
            pw.println("        }");

            pw.println("        function showMarker(id) {");
            pw.println("            $('#haplotypeView').hide();");
            pw.println("            $('.marker-item').removeClass('active'); $(`.marker-item[data-id=\"${id}\"]`).addClass('active');");
            pw.println("            $('#emptyState').hide(); $('#detailView').show(); $('#selectedMarkerTitle').text(id);");
            pw.println("            const genes = data[id]; $('#markerPosBadge').text(genes[0].pos + ' bp');");
            pw.println("            if (table) table.destroy();");
            pw.println("            const tbody = $('#geneTableBody').empty();");
            pw.println("            genes.forEach((g, idx) => {");
            pw.println("                let go = g.go.map(t => `<span class='badge bg-go'>${t}</span>`).join('');");
            pw.println("                let stars = g.kasp >= 80 ? '⭐⭐⭐' : (g.kasp >= 50 ? '⭐⭐' : '⭐');");
            pw.println("                let kaspColor = g.kasp >= 80 ? '#10b981' : (g.kasp >= 50 ? '#f59e0b' : '#ef4444');");
            pw.println("                let btns = `<button class='badge' style='background:#6366f1' onclick='showSeq(\"${id}\",${idx},\"pr\")'>Protein</button>`;");
            pw.println("                if (g.cd) btns += `<button class='badge' style='background:#8b5cf6' onclick='showSeq(\"${id}\",${idx},\"cd\")'>CDS</button>`;");
            pw.println("                btns += `<button class='badge' style='background:#10b981' onclick='showHeatmap(\"${id}\",${idx})'>Genotypes</button>`;");
            pw.println("                if (g.pr) btns += `<button class='badge' style='background:#f59e0b' onclick='blast(\"${g.pr}\")'>BLAST</button>`;");
            pw.println("                tbody.append(`<tr>");
            pw.println("                    <td><b>${g.id}</b> <small>(${g.src})</small><br>${drawGeneMap(g)}</td>");
            pw.println("                    <td><span class='badge ${g.ef.includes('Missense')?'bg-mis':(g.ef.includes('Synonymous')?'bg-syn':(g.ef.includes('Stop-Gain')?'bg-stop':'bg-intron'))}'>${g.ef}</span><br><small>${g.ra} &rarr; ${g.aa}</small></td>");
            pw.println("                    <td>${g.d} bp</td>");
            pw.println("                    <td><div style='color:${kaspColor}; font-weight:600'>${stars} (${g.kasp})</div><button class='badge' style='background:#64748b; margin-top:5px' onclick='showKasp(\"${id}\",${idx})'>KASP Primer</button></td>");
            pw.println("                    <td><small>${g.n}</small><br>${go}</td>");
            pw.println("                    <td>${btns}</td>");
            pw.println("                </tr>`);");
            pw.println("            });");
            pw.println("            table = $('#geneTable').DataTable({ pageLength: 5, dom: 'Bfrtip', buttons: ['copy', 'csv', 'excel'], order: [[2, 'asc']] });");
            pw.println("        }");

            pw.println("        function showGlobalResults(query) {");
            pw.println("            $('#emptyState').hide(); $('#detailView').show(); $('#selectedMarkerTitle').html(`Search: ${query}`);");
            pw.println("            if (table) table.destroy(); const tbody = $('#geneTableBody').empty(); const val = query.toLowerCase();");
            pw.println("            for (const id in data) {");
            pw.println("                data[id].forEach((g, idx) => {");
            pw.println("                    if (g.n.toLowerCase().includes(val) || g.go.some(t => t.toLowerCase().includes(val)) || id.toLowerCase().includes(val)) {");
            pw.println("                        let btns = `<button class='badge' style='background:#6366f1' onclick='showSeq(\"${id}\",${idx},\"pr\")'>Prot</button>`;");
            pw.println("                        if (g.cd) btns += `<button class='badge' style='background:#8b5cf6' onclick='showSeq(\"${id}\",${idx},\"cd\")'>CDS</button>`;");
            pw.println("                        tbody.append(`<tr><td><b>${g.id}</b> <small>(${id})</small><br>${drawGeneMap(g)}</td><td><span class='badge bg-intron'>${g.ef}</span></td><td>${g.d}</td><td><small>${g.n}</small></td><td>${btns}</td></tr>`);");
            pw.println("                    }");
            pw.println("                });");
            pw.println("            }");
            pw.println("            table = $('#geneTable').DataTable({ pageLength: 10, dom: 'Bfrtip', buttons: ['copy', 'csv', 'excel'] });");
            pw.println("        }");

            pw.println("        function showHeatmap(mId, gIdx) {");
            pw.println("            const g = data[mId][gIdx]; $('#seqTitle').text('Genotypes: ' + mId);");
            pw.println("            let html = '<div style=\"display:grid; grid-template-columns: repeat(auto-fill, minmax(100px, 1fr)); gap:5px\">';");
            pw.println("            samples.forEach((s, i) => { let gt = g.gt[i] || './.'; let color = gt.includes('1') ? '#f59e0b' : '#e2e8f0'; html += `<div style=\"background:${color}; padding:5px; border-radius:5px; font-size:0.6rem\">${s}<br>${gt}</div>`; });");
            pw.println("            $('#seqText').html(html + '</div>'); $('#seqModal').css('display', 'flex');");
            pw.println("        }");

            pw.println("        function showSeq(mId, gIdx, type) { const g = data[mId][gIdx]; const seq = type==='pr' ? g.pr : g.cd; $('#seqTitle').text((type==='pr'?'Protein':'CDS') + ' - ' + g.id); $('#seqText').text(seq); $('#seqModal').css('display', 'flex'); }");
            pw.println("        function copySeq() { const text = $('#seqText').text(); navigator.clipboard.writeText(text); alert('Sequence copied!'); }");
            pw.println("        function blast(seq) { const form = document.createElement('form'); form.method='POST'; form.action='https://blast.ncbi.nlm.nih.gov/Blast.cgi'; form.target='_blank'; const cmd=document.createElement('input'); cmd.type='hidden'; cmd.name='CMD'; cmd.value='Web'; const prog=document.createElement('input'); prog.type='hidden'; prog.name='PROGRAM'; prog.value='blastp'; const db=document.createElement('input'); db.type='hidden'; db.name='DATABASE'; db.value='nr'; const query=document.createElement('input'); query.type='hidden'; query.name='QUERY'; query.value=seq; form.appendChild(cmd); form.appendChild(prog); form.appendChild(db); form.appendChild(query); document.body.appendChild(form); form.submit(); document.body.removeChild(form); }");

            pw.println("        function drawGeneMap(g) {");
            pw.println("            const p = 10, w = 300, h = 40; const geneLen = g.gEnd - g.gStart; const scale = (val) => p + ((val - g.gStart) / geneLen) * (w - 2 * p);");
            pw.println("            let svg = `<svg class='gene-map-svg' viewBox='0 0 ${w} ${h}'><line x1='${p}' y1='25' x2='${w-p}' y2='25' stroke='#cbd5e1' stroke-width='2' />`;");
            pw.println("            g.sb.filter(s => s.t === 'exon' || s.t === 'cds').forEach(s => { svg += `<rect x='${scale(s.s)}' y='15' width='${Math.max(2, scale(s.e)-scale(s.s))}' height='20' fill='${s.t==='cds'?'var(--accent)':'#e2e8f0'}' rx='2' />`; });");
            pw.println("            let mx = scale(g.pos); svg += `<line x1='${mx}' y1='5' x2='${mx}' y2='25' stroke='#ef4444' stroke-width='2' stroke-dasharray='2,2' /><circle cx='${mx}' cy='5' r='3' fill='#ef4444' /></svg>`;");
            pw.println("            return svg;");
            pw.println("        }");
            
            pw.println("        function downloadAllFASTA() { let fasta = ''; for (let id in data) { data[id].forEach(g => { if (g.pr) fasta += `>${g.id} | ${id} | ${g.ef}\\n${g.pr}\\n`; }); } const blob = new Blob([fasta], {type: 'text/plain'}); const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'candidates.faa'; a.click(); }");

            pw.println("    </script>");
            pw.println("</body>");
            pw.println("</html>");
        }
    }

    private String toJsonArray(Set<String> set) {
        if (set == null || set.isEmpty()) return "[]";
        StringBuilder sb = new StringBuilder("[");
        for (String s : set) sb.append("'").append(s.replace("'", "\\'")).append("',");
        sb.setLength(sb.length() - 1);
        sb.append("]");
        return sb.toString();
    }
}
