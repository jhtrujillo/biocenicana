package org.cenicana.bio.io;

import org.cenicana.bio.core.GwasEngine;
import org.cenicana.bio.core.GwasEngine.GwasHit;
import org.cenicana.bio.core.GwasEngine.GwasInteraction;
import java.io.*;
import java.util.*;

/**
 * Generates an interactive, premium HTML dashboard for GWAS results with a LIGHT theme.
 * Fixed QQ plot logic to handle thinning without shifting axes.
 */
public class GwasDashboardGenerator {

    public void generate(GwasEngine.GwasResult result, String traitName, String outputPath, int ploidy) throws IOException {
        List<GwasHit> hits = result.hits;
        List<GwasInteraction> interactions = result.interactions;
        // 1. Calculate Chromosome Max Positions and Group Scaffolds
        Map<String, Long> chromMaxPos = new TreeMap<>();
        for (GwasHit h : hits) {
            chromMaxPos.put(h.chromosome, Math.max(chromMaxPos.getOrDefault(h.chromosome, 0L), h.position));
        }

        List<String> sortedChroms = new ArrayList<>(chromMaxPos.keySet());
        sortedChroms.sort((a, b) -> {
            try {
                int ia = Integer.parseInt(a.replaceAll("\\D", ""));
                int ib = Integer.parseInt(b.replaceAll("\\D", ""));
                return Integer.compare(ia, ib);
            } catch (Exception e) {
                return a.compareTo(b);
            }
        });

        Map<String, Long> chromOffsets = new LinkedHashMap<>();
        long currentOffset = 0;
        List<String> mainChroms = new ArrayList<>();
        List<String> scaffolds = new ArrayList<>();
        
        for (String c : sortedChroms) {
            String lower = c.toLowerCase();
            if (lower.startsWith("chr") || lower.startsWith("ch") || lower.matches("\\d+")) {
                mainChroms.add(c);
            } else {
                scaffolds.add(c);
            }
        }

        for (String chrom : mainChroms) {
            chromOffsets.put(chrom, currentOffset);
            currentOffset += chromMaxPos.get(chrom) + 5_000_000;
        }

        if (!scaffolds.isEmpty()) {
            long scaffoldBaseOffset = currentOffset;
            for (String chrom : scaffolds) {
                chromOffsets.put(chrom, scaffoldBaseOffset);
            }
        }

        // Calculate Lambda GC
        double lambdaGC = calculateLambdaGC(hits);

        // Calculate FDR Threshold for plot (the cutoff point)
        double fdrThresholdLogP = Double.MAX_VALUE;
        int fdrCount = 0;
        for (GwasHit h : hits) {
            if (h.qValue < 0.05) {
                fdrCount++;
                fdrThresholdLogP = Math.min(fdrThresholdLogP, -Math.log10(h.pValue));
            }
        }
        if (fdrCount == 0) fdrThresholdLogP = 0;

        try (PrintWriter pw = new PrintWriter(new FileWriter(outputPath))) {
            pw.println("<!DOCTYPE html>");
            pw.println("<html lang='en'>");
            pw.println("<head>");
            pw.println("    <meta charset='UTF-8'>");
            pw.println("    <title>BioJava GWAS | " + traitName + "</title>");
            pw.println("    <script src='https://cdn.plot.ly/plotly-2.24.1.min.js'></script>");
            pw.println("    <link href='https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap' rel='stylesheet'>");
            pw.println("    <style>");
            pw.println("        :root { --primary: #4f46e5; --secondary: #3b82f6; --bg: #f8fafc; --card: #ffffff; --text: #1e293b; --text-dim: #64748b; }");
            pw.println("        body { font-family: 'Inter', sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 20px; }");
            pw.println("        #loader { position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(255,255,255,0.9); z-index: 9999; display: flex; flex-direction: column; align-items: center; justify-content: center; transition: opacity 0.5s; }");
            pw.println("        .spinner { width: 50px; height: 50px; border: 5px solid #f3f3f3; border-top: 5px solid var(--primary); border-radius: 50%; animation: spin 1s linear infinite; }");
            pw.println("        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }");
            pw.println("        .card { background: rgba(255, 255, 255, 0.8); backdrop-filter: blur(10px); padding: 25px; border-radius: 20px; box-shadow: 0 10px 30px rgba(0,0,0,0.05); }");
            pw.println("        .glass { background: var(--card); backdrop-filter: blur(12px); border: 1px solid rgba(255,255,255,0.4); border-radius: 24px; box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.05); }");
            pw.println("        .header { padding: 40px; display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid rgba(0,0,0,0.05); }");
            pw.println("        h1 { margin: 0; font-weight: 700; letter-spacing: -0.02em; background: linear-gradient(to right, #4f46e5, #9333ea); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }");
            pw.println("        .container { padding: 40px; max-width: 1600px; margin: 0 auto; }");
            pw.println("        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 30px; }");
            pw.println("        .stat-card { padding: 24px; text-align: center; transition: transform 0.3s ease; }");
            pw.println("        .stat-card:hover { transform: translateY(-5px); }");
            pw.println("        .stat-val { font-size: 2rem; font-weight: 700; color: var(--primary); }");
            pw.println("        .stat-label { font-size: 0.8rem; color: var(--text-dim); text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px; }");
            pw.println("        .main-grid { display: grid; grid-template-columns: 2.5fr 1fr; gap: 30px; margin-bottom: 30px; }");
            pw.println("        .card-title { font-size: 1.25rem; font-weight: 600; margin-bottom: 20px; display: flex; align-items: center; gap: 10px; color: var(--text); }");
            pw.println("        .toggle-container { background: #eee; border-radius: 30px; padding: 5px; display: flex; gap: 5px; }");
            pw.println("        .toggle-container button { border: none; padding: 8px 20px; border-radius: 25px; cursor: pointer; transition: 0.3s; background: transparent; }");
            pw.println("        .toggle-container button.active { background: var(--primary); color: white; box-shadow: 0 4px 10px rgba(79, 70, 229, 0.3); }");
            pw.println("        select, input { background: rgba(255, 255, 255, 0.9); border: 1px solid rgba(0,0,0,0.1); color: var(--text); padding: 10px 16px; border-radius: 12px; font-family: inherit; outline: none; }");
            pw.println("        table { width: 100%; border-collapse: separate; border-spacing: 0 8px; margin-top: -8px; }");
            pw.println("        th { text-align: left; padding: 16px; color: var(--text-dim); font-weight: 400; font-size: 0.85rem; }");
            pw.println("        td { padding: 16px; background: rgba(255,255,255,0.5); border-top: 1px solid rgba(0,0,0,0.02); border-bottom: 1px solid rgba(0,0,0,0.02); }");
            pw.println("        .badge { padding: 6px 12px; border-radius: 100px; font-size: 0.75rem; font-weight: 600; }");
            pw.println("        .sig { background: #dcfce7; color: #15803d; }");
            pw.println("        .non-sig { background: #f1f5f9; color: var(--text-dim); }");
            pw.println("    </style>");
            pw.println("</head>");
            pw.println("<body>");
            pw.println("    <div id='loader'><div class='spinner'></div><p style='margin-top:20px; font-weight:600; color:var(--primary)'>Procesando Datos Genómicos...</p></div>");
            pw.println("    <div class='header'>");
            pw.println("        <div><h1>" + traitName + " <span style='font-weight:300; font-size:1.2rem; opacity:0.6; color:var(--text)'>| GWAS Insight</span></h1></div>");
            pw.println("        <div style='opacity: 0.8'><img src='https://www.cenicana.org/wp-content/uploads/2021/04/Logo-Cenicana-color.png' height='45'></div>");
            pw.println("    </div>");
            pw.println("    <div class='container'>");
            
            pw.println("        <div class='stats-grid'>");
            pw.println("            <div class='glass stat-card'><div class='stat-val'>" + hits.size() + "</div><div class='stat-label'>Total Variants</div></div>");
            
            // Dynamic Lambda GC interpretation
            String lambdaStatus;
            String lambdaColor;
            String lambdaDesc;
            if (lambdaGC < 0.9) {
                lambdaStatus = "CONSERVATIVE";
                lambdaColor = "#3b82f6";
                lambdaDesc = "Over-correction: Real peaks might be hidden.";
            } else if (lambdaGC > 1.1) {
                lambdaStatus = "INFLATED";
                lambdaColor = "#ef4444";
                lambdaDesc = "Inflation: Possible false positives.";
            } else {
                lambdaStatus = "IDEAL";
                lambdaColor = "#10b981";
                lambdaDesc = "Perfectly controlled population structure.";
            }
            
            pw.println("            <div class='glass stat-card'>");
            pw.println("                <div class='stat-val' style='color:" + lambdaColor + "'>" + String.format("%.3f", lambdaGC) + "</div>");
            pw.println("                <div class='stat-label'>Lambda GC | <span style='font-weight:700'>" + lambdaStatus + "</span></div>");
            pw.println("                <div style='font-size:0.7rem; color:var(--text-dim); margin-top:5px;'>" + lambdaDesc + "</div>");
            pw.println("            </div>");

            pw.println("            <div class='glass stat-card'><div class='stat-val' id='sigCount'>0</div><div class='stat-label'>Significant (Bonf)</div></div>");
            pw.println("            <div class='glass stat-card'><div class='stat-val' style='color:var(--secondary)'>" + fdrCount + "</div><div class='stat-label'>Significant (FDR)</div></div>");
            pw.println("        </div>");

            pw.println("        <div class='main-grid'>");
            pw.println("            <div class='glass card'>");
            pw.println("                <div style='display:flex; justify-content:space-between; align-items:center;'>");
            pw.println("                    <h2 style='margin:0; font-size:1.25rem;'>Manhattan Discovery Plot</h2>");
            pw.println("                    <div class='toggle-container'>");
            pw.println("                        <button id='btnLinear' class='active' onclick='showView(\"linear\")'>Lineal</button>");
            pw.println("                        <button id='btnCircular' onclick='showView(\"circular\")'>Circular</button>");
            pw.println("                    </div>");
            pw.println("                </div>");
            pw.println("                <div id='manhattanPlot' style='height:500px;'></div>");
            pw.println("                <div id='circularPlot' style='height:600px; display:none;'></div>");
            pw.println("            </div>");
            pw.println("            <div class='glass card'>");
            pw.println("                <div class='card-title'>QQ Inflation</div>");
            pw.println("                <div id='qqplot' style='height: 500px;'></div>");
            pw.println("            </div>");
            pw.println("        </div>");
            
            if (interactions != null && !interactions.isEmpty()) {
                pw.println("        <div class='glass card' style='margin-bottom:30px;'>");
                pw.println("            <div class='card-title'>⚡ Sinergias Genéticas (Epistasia)</div>");
                pw.println("            <div style='max-height:300px; overflow-y:auto;'>");
                pw.println("                <table>");
                pw.println("                    <thead><tr><th>Lead Marker</th><th>Interacts with</th><th>P-Value</th><th>Effect</th></tr></thead>");
                pw.println("                    <tbody>");
                for (GwasInteraction gi : interactions) {
                    pw.println("                        <tr><td>" + gi.marker1 + "</td><td>" + gi.marker2 + "</td><td>" + String.format("%.2e", gi.pValue) + "</td><td>" + String.format("%.3f", gi.effect) + "</td></tr>");
                }
                pw.println("                    </tbody>");
                pw.println("                </table>");
                pw.println("            </div>");
                pw.println("        </div>");
            }

            pw.println("        <div class='glass card'>");
            pw.println("            <div style='display: flex; justify-content: space-between; align-items: center; margin-bottom: 25px;'>");
            pw.println("                <div class='card-title'>Associated Markers</div>");
            pw.println("                <div style='display: flex; gap: 12px; align-items: center;'>");
            pw.println("                    <label style='font-size:0.85rem; color:var(--text-dim); display:flex; align-items:center; gap:6px;'>");
            pw.println("                        <input type='checkbox' id='leadsOnly' onchange='renderTable()' style='width:16px; height:16px;'> Solo Leads (Clumping 5Mb)");
            pw.println("                    </label>");
            pw.println("                    <input type='text' id='search' placeholder='Search marker...' onkeyup='renderTable()'>");
            pw.println("                    <select id='pageSize' onchange='renderTable()'>");
            pw.println("                        <option value='20'>20 Results</option>");
            pw.println("                        <option value='50' selected>50 Results</option>");
            pw.println("                        <option value='100'>100 Results</option>");
            pw.println("                    </select>");
            pw.println("                </div>");
            pw.println("            </div>");
            pw.println("            <div style='overflow-x: auto;'>");
            pw.println("                <table><thead><tr><th>MARKER</th><th>CHR</th><th>POSITION</th><th>BEST MODEL</th><th>AIC</th><th>EFFECT</th><th>R²</th><th>Q-VALUE (FDR)</th><th>P-VALUE</th></tr></thead>");
            pw.println("                <tbody id='tableBody'></tbody></table>");
            pw.println("            </div>");
            pw.println("        </div>");
            pw.println("        <div class='glass card' id='selectionAnalysis' style='margin-top: 30px; display: none;'>");
            pw.println("            <div class='card-title'>Selection Analysis & Allele Effect</div>");
            pw.println("            <div style='display: grid; grid-template-columns: 1fr 2.5fr; gap: 30px;'>");
            pw.println("                <div id='selectionInfo' style='padding: 20px; background: rgba(79, 70, 229, 0.05); border-radius: 16px;'>");
            pw.println("                    <h3 id='selMarkerName' style='margin-top: 0; color: var(--primary);'></h3>");
            pw.println("                    <div id='selMarkerStats'></div>");
            pw.println("                </div>");
            pw.println("                <div id='selBoxplot' style='height: 400px;'></div>");
            pw.println("            </div>");
            pw.println("        </div>");
            pw.println("    </div>");

            pw.println("    <script>");
            pw.println("        window.onload = () => document.getElementById('loader').style.opacity = '0';");
            pw.println("        setTimeout(() => document.getElementById('loader').style.display = 'none', 500);");
            pw.println("        const totalM = " + hits.size() + ";");
            pw.println("        const palette = ['#4f46e5', '#7c3aed', '#db2777', '#dc2626', '#d97706', '#059669', '#0891b2', '#2563eb'];");
            pw.println("        const bonferroni = " + (-Math.log10(0.05 / Math.max(1, hits.size()))) + ";");
            
            pw.println("        const all_data = [");
            int sigCount = 0;
            for (int j = 0; j < hits.size(); j++) {
                GwasHit h = hits.get(j);
                double logP = -Math.log10(h.pValue);
                double displayLogP = Math.min(logP, 80.0); // Cap visual at 80
                if (h.pValue < (0.05 / hits.size())) sigCount++;

                if (h.pValue > 0.05 && j % 10 != 0) continue;

                long offset = chromOffsets.getOrDefault(h.chromosome, 0L);
                String chromNum = h.chromosome.replaceAll("\\D", "");
                int cIdx = chromNum.isEmpty() ? h.chromosome.hashCode() : Integer.parseInt(chromNum);
                String color = "palette[" + (Math.abs(cIdx) % 8) + "]";
                String displayChrom = scaffolds.contains(h.chromosome) ? "Scaffold" : h.chromosome;
                String symbol = "circle";
                if (h.model != null) {
                    switch (h.model.toString()) {
                        case "ADDITIVE": symbol = "triangle-up"; break;
                        case "SIMPLEX_DOMINANT": symbol = "circle"; break;
                        case "DUPLEX_DOMINANT": symbol = "square"; break;
                        case "TRIPLEX_DOMINANT": symbol = "cross"; break;
                        case "SIMPLEX_DOMINANT_REF": symbol = "circle-open"; break;
                        case "DUPLEX_DOMINANT_REF": symbol = "square-open"; break;
                        case "TRIPLEX_DOMINANT_REF": symbol = "cross-open"; break;
                        case "GENERAL": symbol = "diamond"; break;
                    }
                }

                pw.print("{x:" + (offset + h.position) + ",y:" + displayLogP + ",id:'" + h.markerId + "',chr:'" + displayChrom + "',c:" + color + ",s:'" + symbol + "'}");
                if (j < hits.size() - 1) pw.print(",");
            }
            pw.println("        ];");
            pw.println("        document.getElementById('sigCount').innerText = '" + sigCount + "';");

            pw.println("        const ticks = { val: [" + String.join(",", mainChroms.stream().map(c -> String.valueOf(chromOffsets.get(c) + chromMaxPos.get(c)/2)).toArray(String[]::new)) + "], text: [" + String.join(",", mainChroms.stream().map(c -> "'" + c + "'").toArray(String[]::new)) + "] };");

            pw.println("        function updateManhattan() {");
            pw.println("            Plotly.react('manhattanPlot', [{");
            pw.println("                x: all_data.map(d => d.x), y: all_data.map(d => d.y), text: all_data.map(d => d.id),");
            pw.println("                mode: 'markers', type: 'scattergl', marker: { color: all_data.map(d => d.c), size: 7, opacity: 0.7, symbol: all_data.map(d => d.s) }");
            pw.println("            }], {");
            pw.println("                paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',");
            pw.println("                xaxis: { gridcolor: 'rgba(0,0,0,0.05)', tickfont: {color: '#64748b'}, tickvals: ticks.val, ticktext: ticks.text },");
            pw.println("                yaxis: { gridcolor: 'rgba(0,0,0,0.05)', tickfont: {color: '#64748b'}, title: '-log10(p)' },");
            pw.println("                margin: { t: 50, b: 40, l: 50, r: 20 }, hovermode: 'closest',");
            pw.println("                shapes: [");
            pw.println("                    { type: 'line', x0: 0, x1: " + currentOffset + ", y0: bonferroni, y1: bonferroni, line: { color: '#ef4444', width: 2, dash: 'dash' } },");
            pw.println("                    { type: 'line', x0: 0, x1: " + currentOffset + ", y0: " + fdrThresholdLogP + ", y1: " + fdrThresholdLogP + ", line: { color: '#3b82f6', width: 2, dash: 'dot' } }");
            pw.println("                ],");
            pw.println("                annotations: [");
            pw.println("                    { x: " + currentOffset + ", y: bonferroni, xref: 'x', yref: 'y', text: 'Bonferroni', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: {color: '#ef4444', size: 10} },");
            pw.println("                    { x: " + currentOffset + ", y: " + fdrThresholdLogP + ", xref: 'x', yref: 'y', text: 'FDR 0.05', showarrow: false, xanchor: 'right', yanchor: 'bottom', font: {color: '#3b82f6', size: 10} },");
            pw.println("                    { x: 0.05, y: 1.05, xref: 'paper', yref: 'paper', text: '▲ Add | ● 1-Dom | ■ 2-Dom | ✖ 3-Dom | ♦ Gen', showarrow: false, font: {size: 11, color: '#4f46e5'}, xanchor: 'left' }");
            pw.println("                ]");
            pw.println("            }, { responsive: true, displayModeBar: false });");
            pw.println("        }");
            pw.println("        updateManhattan();");


            pw.println("        const top_hits = [");
            for (int j = 0; j < Math.min(2000, hits.size()); j++) {
                GwasHit h = hits.get(j);
                if (!h.isBestModel && h.pValue > 1e-4) continue;
                pw.print("{id:'" + h.markerId + "',chr:'" + h.chromosome + "',pos:" + h.position + ",ref:'" + h.refAllele + "',alt:'" + h.altAllele + "',model:'" + h.model + "',eff:" + h.effect + ",r2:" + h.r2 + ",p:" + h.pValue + ",q:" + h.qValue + ",aic:" + h.aic + ",isBest:" + h.isBestModel);
                if (h.phenotypesByDosage != null && j < 150) {
                    pw.print(",dist:[" + formatDosages(h.phenotypesByDosage) + "],samples:[" + formatSamples(h.samplesByDosage) + "]");
                }
                pw.print("}");
                if (j < Math.min(2000, hits.size()) - 1) pw.print(",");
            }
            pw.println("        ];");

            pw.println("        function getLeads(data) {");
            pw.println("            const significant = data.filter(h => h.p < (0.05 / totalM)).sort((a,b) => a.p - b.p);");
            pw.println("            const leads = [];");
            pw.println("            const windowBP = 5000000; // 5MB window");
            pw.println("            significant.forEach(h => {");
            pw.println("                const isProximal = leads.some(l => l.chr === h.chr && Math.abs(l.pos - h.pos) < windowBP);");
            pw.println("                if (!isProximal) leads.push(h);");
            pw.println("            });");
            pw.println("            return leads;");
            pw.println("        }");

            pw.println("        function renderTable() {");
            pw.println("            const loader = document.getElementById('loader');");
            pw.println("            if(loader) loader.style.display = 'flex';");
            pw.println("            setTimeout(() => {");
            pw.println("                const q = document.getElementById('search').value.toLowerCase();");
            pw.println("                const size = parseInt(document.getElementById('pageSize').value);");
            pw.println("                const leadsOnly = document.getElementById('leadsOnly').checked;");
            pw.println("                ");
            pw.println("                let filtered = top_hits.filter(h => h.isBest && h.id.toLowerCase().includes(q));");
            pw.println("                if (leadsOnly) {");
            pw.println("                    filtered = getLeads(filtered);");
            pw.println("                }");
            pw.println("                ");
            pw.println("                const sliced = filtered.slice(0, size);");
            pw.println("                document.getElementById('tableBody').innerHTML = sliced.map(h => `<tr onclick=\"showSelectionAnalysis('${h.id}')\" style='cursor:pointer'>");
            pw.println("                    <td style='font-weight:600; color:var(--text)'>${h.id}</td><td>${h.chr}</td><td>${h.pos.toLocaleString()}</td>");
            pw.println("                    <td><span class='badge' style='background:#e0e7ff; color:#3730a3'>${h.model}</span></td>");
            pw.println("                    <td style='color:var(--text-dim); font-size:0.8rem'>${h.aic.toFixed(1)}</td>");
            pw.println("                    <td style='font-family:monospace'>${h.eff.toFixed(4)}</td>");
            pw.println("                    <td style='font-weight:600; color:var(--secondary)'>${(h.r2 * 100).toFixed(1)}%</td>");
            pw.println("                    <td style='color:${h.q < 0.05 ? \"#15803d\" : \"inherit\"}'>${h.q.toExponential(2)}</td>");
            pw.println("                    <td><span class='badge ${h.p < (0.05/ totalM) ? \"sig\" : \"non-sig\"}'>${h.p.toExponential(2)}</span></td>");
            pw.println("                </tr>`).join('');");
            pw.println("                if(loader) loader.style.display = 'none';");
            pw.println("            }, 50);");
            pw.println("        }");
            pw.println("        renderTable();");

            pw.println("        function showView(view) {");
            pw.println("            if (view === 'linear') {");
            pw.println("                document.getElementById('manhattanPlot').style.display = 'block';");
            pw.println("                document.getElementById('circularPlot').style.display = 'none';");
            pw.println("                document.getElementById('btnLinear').classList.add('active');");
            pw.println("                document.getElementById('btnCircular').classList.remove('active');");
            pw.println("            } else {");
            pw.println("                document.getElementById('manhattanPlot').style.display = 'none';");
            pw.println("                document.getElementById('circularPlot').style.display = 'block';");
            pw.println("                document.getElementById('btnLinear').classList.remove('active');");
            pw.println("                document.getElementById('btnCircular').classList.add('active');");
            pw.println("                renderCircular();");
            pw.println("            }");
            pw.println("        }");
            pw.println("");
            pw.println("        function renderCircular() {");
            pw.println("            const chroms = [...new Set(all_data.map(d => d.chr))].sort();");
            pw.println("            const traces = chroms.map((c, idx) => {");
            pw.println("                const cData = all_data.filter(d => d.chr === c);");
            pw.println("                const maxPos = Math.max(...cData.map(d => d.x));");
            pw.println("                const offset = (idx / chroms.length) * 360;");
            pw.println("                const span = (1 / chroms.length) * 340;");
            pw.println("                return {");
            pw.println("                    type: 'scatterpolar', mode: 'markers',");
            pw.println("                    r: cData.map(d => d.y),");
            pw.println("                    theta: cData.map(d => offset + ((d.x % maxPos) / maxPos) * span),");
            pw.println("                    text: cData.map(d => `${d.id}<br>P: ${d.y.toFixed(2)}`),");
            pw.println("                    name: c, marker: { size: 4 }");
            pw.println("                };");
            pw.println("            });");
            pw.println("            Plotly.newPlot('circularPlot', traces, {");
            pw.println("                polar: { radialaxis: { title: '-log10(p)' }, angularaxis: { showticklabels: false } },");
            pw.println("                margin: { t: 40, b: 40 }");
            pw.println("            }, { responsive: true, displayModeBar: false });");
            pw.println("        }");
            pw.println("");
            pw.println("        function showSelectionAnalysis(markerId) {");
            pw.println("            const h = top_hits.find(x => x.id === markerId);");
            pw.println("            if (!h || !h.dist) return;");
            pw.println("");
            pw.println("            document.getElementById('selectionAnalysis').style.display = 'block';");
            pw.println("            document.getElementById('selMarkerName').innerText = h.id;");
            pw.println("            ");
            pw.println("            const action = h.eff > 0 ? 'AUMENTAR' : 'DISMINUIR';");
            pw.println("            const inverseAction = h.eff > 0 ? 'disminuir' : 'aumentar';");
            pw.println("            const color = h.eff > 0 ? '#059669' : '#dc2626';");
            pw.println("            ");
            pw.println("            document.getElementById('selMarkerStats').innerHTML = `");
            pw.println("                <p><b>Cromosoma:</b> ${h.chr}</p>");
            pw.println("                <p><b>Modelo:</b> ${h.model}</p>");
            pw.println("                <p><b>Efecto por alelo ALT:</b> <span style='color:${color}; font-weight:bold'>${h.eff.toFixed(4)}</span></p>");
            pw.println("                <p><b>Varianza Explicada (R²):</b> ${(h.r2*100).toFixed(2)}%</p>");
            pw.println("                <div style='margin-top:20px; padding:15px; background:white; border-radius:12px; border-left:4px solid ${color}'>");
            pw.println("                    <small style='color:var(--text-dim)'>GUÍA DE SELECCIÓN GENÓMICA</small><br>");
            pw.println("                    <strong>El alelo ALT (Dosis) tiende a ${action.toLowerCase()} el rasgo.</strong><br><br>");
            pw.println("                    <div style='font-size:13px;'>");
            pw.println("                        • Si busca <b>${action}</b> el rasgo: Seleccione individuos con <b>Dosis " + (ploidy) + "</b>.<br>");
            pw.println("                        • Si busca <b>${inverseAction.toUpperCase()}</b> el rasgo: Seleccione individuos con <b>Dosis 0</b>.");
            pw.println("                    </div>");
            pw.println("                </div>");
            pw.println("                <div style='margin-top:20px;'>");
            pw.println("                    <p><b>Candidatos Élite (Genotipo Ideal):</b></p>");
            pw.println("                    <div id='eliteList' style='max-height:250px; overflow-y:auto; font-size:12px; background:white; padding:10px; border-radius:8px; border:1px solid rgba(0,0,0,0.05)'></div>");
            pw.println("                </div>");
            pw.println("            `;");
            pw.println("");
            pw.println("            const getGenotype = (d, ref, alt) => {");
            pw.println("                return alt.repeat(d) + ref.repeat(" + ploidy + " - d);");
            pw.println("            };");
            pw.println("");
            pw.println("            const bestDosage = h.eff > 0 ? " + ploidy + " : 0;");
            pw.println("            const eliteSamples = h.samples[bestDosage] || [];");
            pw.println("            const geno = getGenotype(bestDosage, h.ref, h.alt);");
            pw.println("            if (eliteSamples.length > 0) {");
            pw.println("                let html = '<table style=\"width:100%; border-collapse:collapse; margin-top:10px\">';");
            pw.println("                html += '<tr style=\"background:rgba(0,0,0,0.02); font-weight:bold\"><td>Individuo</td><td>Dosis</td><td>Genotipo</td></tr>';");
            pw.println("                html += eliteSamples.map(s => `<tr><td style=\"padding:4px\">${s}</td><td>${bestDosage}</td><td style=\"font-family:monospace; color:var(--primary)\">${geno}</td></tr>`).join('');");
            pw.println("                html += '</table>';");
            pw.println("                document.getElementById('eliteList').innerHTML = html;");
            pw.println("            } else {");
            pw.println("                document.getElementById('eliteList').innerHTML = '<div style=\"padding:10px; color:var(--text-dim)\">No hay individuos con el genotipo ideal (Dosis ' + bestDosage + ')</div>';");
            pw.println("            }");
            pw.println("");
            pw.println("            const traces = h.dist.map((values, dosage) => ({");
            pw.println("                y: values, type: 'box', name: 'Dosis ' + dosage, ");
            pw.println("                marker: { color: palette[dosage % palette.length] },");
            pw.println("                boxpoints: 'all', jitter: 0.3, pointpos: -1.8");
            pw.println("            })).filter(t => t.y.length > 0);");
            pw.println("");
            pw.println("            Plotly.newPlot('selBoxplot', traces, {");
            pw.println("                paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',");
            pw.println("                title: 'Phenotype Distribution by Allele Dosage',");
            pw.println("                yaxis: { gridcolor: 'rgba(0,0,0,0.05)', title: 'Phenotype Value' },");
            pw.println("                xaxis: { title: 'Allele Dosage (Copies)' },");
            pw.println("                margin: { t: 40, l: 50, r: 10, b: 50 }, showlegend: false");
            pw.println("            }, { responsive: true, displayModeBar: false });");
            pw.println("            ");
            pw.println("            window.scrollTo({ top: document.getElementById('selectionAnalysis').offsetTop - 20, behavior: 'smooth' });");
            pw.println("        }");

            pw.println("        // QQ Plot Data with Correct Ranks");
            pw.println("        const obsData = " + getRankedPValuesJson(hits) + ";");
            pw.println("        const expLogP = obsData.map(d => -Math.log10(d.r / totalM));");
            pw.println("        const obsLogP = obsData.map(d => -Math.log10(d.p));");
            pw.println("        Plotly.newPlot('qqplot', [");
            pw.println("            { x: expLogP, y: obsLogP, mode: 'markers', type: 'scattergl', marker: { color: '#4f46e5', size: 4, opacity: 0.5 } },");
            pw.println("            { x: [0, Math.max(...expLogP)], y: [0, Math.max(...expLogP)], mode: 'lines', line: { color: '#ef4444', dash: 'dash', width: 1 } }");
            pw.println("        ], {");
            pw.println("            paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',");
            pw.println("            xaxis: { title: 'Expected -log10(p)', gridcolor: 'rgba(0,0,0,0.05)', tickfont: {color: '#64748b'} },");
            pw.println("            yaxis: { title: 'Observed -log10(p)', gridcolor: 'rgba(0,0,0,0.05)', tickfont: {color: '#64748b'} },");
            pw.println("            margin: { t: 20, l: 50, r: 10, b: 50 }, showlegend: false");
            pw.println("        }, { responsive: true, displayModeBar: false });");
            pw.println("    </script>");
            pw.println("</body>");
            pw.println("</html>");
        }
    }

    private double calculateLambdaGC(List<GwasHit> hits) {
        if (hits.isEmpty()) return 1.0;
        double[] pValues = hits.stream().mapToDouble(h -> h.pValue).sorted().toArray();
        double medianP = pValues[pValues.length / 2];
        double chiObs = inverseChiSquare(1.0 - medianP);
        double chiExp = 0.4549; 
        return chiObs / chiExp;
    }

    private double inverseChiSquare(double p) {
        double z = inverseNormal(0.5 + p / 2.0);
        return z * z;
    }

    private double inverseNormal(double p) {
        double a[] = {-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00};
        double b[] = {-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
        double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00};
        double d[] = { 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,  3.754408661907416e+00};
        double p_low = 0.02425, p_high = 1 - p_low;
        if (p < p_low) {
            double q = Math.sqrt(-2 * Math.log(p));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
        }
        if (p <= p_high) {
            double q = p - 0.5, r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
        }
        double q = Math.sqrt(-2 * Math.log(1 - p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }

    private String formatDosages(List<Double>[] dist) {
        StringBuilder sb = new StringBuilder();
        for (int d = 0; d < dist.length; d++) {
            sb.append("[");
            for (int k = 0; k < dist[d].size(); k++) {
                sb.append(dist[d].get(k));
                if (k < dist[d].size() - 1) sb.append(",");
            }
            sb.append("]");
            if (d < dist.length - 1) sb.append(",");
        }
        return sb.toString();
    }

    private String formatSamples(List<String>[] samples) {
        StringBuilder sb = new StringBuilder();
        for (int d = 0; d < samples.length; d++) {
            sb.append("[");
            for (int k = 0; k < samples[d].size(); k++) {
                sb.append("'").append(samples[d].get(k)).append("'");
                if (k < samples[d].size() - 1) sb.append(",");
            }
            sb.append("]");
            if (d < samples.length - 1) sb.append(",");
        }
        return sb.toString();
    }

    private String getRankedPValuesJson(List<GwasHit> hits) {
        // Use only best models for QQ plot to avoid redundancy
        double[] pValues = hits.stream().filter(h -> h.isBestModel).mapToDouble(h -> h.pValue).sorted().toArray();
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < pValues.length; i++) {
            if (pValues[i] > 0.05 && i % 40 != 0) continue;
            sb.append("{\"p\":").append(pValues[i]).append(",\"r\":").append(i + 1).append("}");
            if (i < pValues.length - 1) sb.append(",");
        }
        if (sb.length() > 1 && sb.charAt(sb.length()-1) == ',') sb.setLength(sb.length()-1);
        sb.append("]");
        return sb.toString();
    }
}
