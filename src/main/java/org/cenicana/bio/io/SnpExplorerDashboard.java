package org.cenicana.bio.io;

import org.cenicana.bio.core.SnpClusterAnalyzer.SnpResult;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class SnpExplorerDashboard {

    public static void generate(List<SnpResult> results, int ploidy, String outputPath) throws IOException {
        try (PrintWriter w = new PrintWriter(new FileWriter(outputPath))) {
            w.println("<!DOCTYPE html><html lang='en'><head>");
            w.println("<meta charset='UTF-8'><title>BioCenicana - SNP Group Explorer</title>");
            w.println("<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>");
            w.println("<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&display=swap' rel='stylesheet'>");
            w.println("<style>");
            w.println("body{font-family:'Inter',sans-serif;background:#f8fafc;margin:0;display:flex;height:100vh}");
            w.println(".sidebar{width:350px;background:white;border-right:1px solid #e2e8f0;display:flex;flex-direction:column}");
            w.println(".header{padding:1.5rem;background:#1e293b;color:white}");
            w.println(".header h1{font-size:1.2rem;margin:0}");
            w.println(".search-box{padding:1rem;border-bottom:1px solid #e2e8f0}");
            w.println("input{width:100%;padding:.6rem;border:1px solid #cbd5e1;border-radius:6px;font-size:.9rem}");
            w.println(".snp-list{flex:1;overflow-y:auto}");
            w.println(".snp-item{padding:1rem;border-bottom:1px solid #f1f5f9;cursor:pointer;transition:background:.2s}");
            w.println(".snp-item:hover{background:#f1f5f9}");
            w.println(".snp-item.active{background:#eff6ff;border-left:4px solid #3b82f6}");
            w.println(".snp-id{font-weight:700;color:#1e293b;font-size:.9rem}");
            w.println(".snp-meta{font-size:.75rem;color:#64748b;margin-top:.2rem}");
            w.println(".main{flex:1;padding:2rem;display:flex;flex-direction:column;gap:1.5rem;overflow-y:auto}");
            w.println(".chart-card{background:white;padding:2rem;border-radius:12px;box-shadow:0 1px 3px rgba(0,0,0,.1);min-height:500px}");
            w.println(".info-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:1rem}");
            w.println(".info-card{background:white;padding:1.5rem;border-radius:12px;box-shadow:0 1px 3px rgba(0,0,0,.1)}");
            w.println(".info-label{font-size:.7rem;text-transform:uppercase;color:#64748b;font-weight:700}");
            w.println(".info-value{font-size:1.4rem;font-weight:800;color:#3b82f6}");
            w.println("</style></head><body>");

            // Sidebar
            w.println("<div class='sidebar'>");
            w.println("<div class='header'><h1>🧬 SNP explorer</h1><p style='font-size:0.7rem;opacity:0.8'>Group distribution across population</p></div>");
            w.println("<div class='search-box'><input type='text' id='search' placeholder='Search SNP (e.g. Chr1_100)...' onkeyup='filterList()'></div>");
            w.println("<div class='snp-list' id='list'>");
            for (int i = 0; i < results.size(); i++) {
                SnpResult r = results.get(i);
                w.printf("<div class='snp-item' id='item-%d' onclick='selectSnp(%d)'>%n", i, i);
                w.printf("<div class='snp-id'>%s</div>%n", r.id);
                w.printf("<div class='snp-meta'>Chr: %s | Pos: %d | Clusters: %d</div>%n", r.chr, r.pos, r.centroids.length);
                w.println("</div>");
            }
            w.println("</div></div>");

            // Main View
            w.println("<div class='main'>");
            w.println("<div class='info-grid'>");
            w.println("<div class='info-card'><div class='info-label'>Selected SNP</div><div class='info-value' id='v-id'>-</div></div>");
            w.println("<div class='info-card'><div class='info-label'>Detected Groups</div><div class='info-value' id='v-groups'>-</div></div>");
            w.println("<div class='info-card'><div class='info-label'>Theoretical Ploidy</div><div class='info-value'>" + ploidy + "</div></div>");
            w.println("</div>");
            w.println("<div class='chart-card' id='chart'></div>");
            w.println("</div>");

            // Data & JS
            w.println("<script>");
            w.println("const snps = [");
            for (int i = 0; i < results.size(); i++) {
                SnpResult r = results.get(i);
                w.print("{id:'" + r.id + "',bins:" + arrayToString(r.histogramBins) + ",centroids:" + arrayToString(r.centroids) + ",counts:" + arrayToString(r.clusterCounts) + "}");
                if (i < results.size() - 1) w.println(",");
            }
            w.println("];");

            w.println("const ploidy = " + ploidy + ";");
            w.println("var activeIdx = 0;");

            w.println("function selectSnp(idx) {");
            w.println("  activeIdx = idx;");
            w.println("  document.querySelectorAll('.snp-item').forEach(el => el.classList.remove('active'));");
            w.println("  document.getElementById('item-' + idx).classList.add('active');");
            w.println("  updateView();");
            w.println("}");

            w.println("function updateView() {");
            w.println("  const snp = snps[activeIdx];");
            w.println("  document.getElementById('v-id').innerText = snp.id;");
            w.println("  document.getElementById('v-groups').innerText = snp.centroids.length;");
            w.println("  const binLabels = []; for(let i=0; i<25; i++) binLabels.push((i/24).toFixed(2));");
            
            w.println("  const traceHist = { x: binLabels, y: snp.bins, type: 'bar', marker: {color: 'rgba(59, 130, 246, 0.6)'}, name: 'Dosage Dist' };");
            
            // Theo lines
            w.println("  const shapes = [];");
            w.println("  for(let i=0; i<=ploidy; i++) {");
            w.println("    let x = i/ploidy;");
            w.println("    shapes.push({type: 'line', x0: x, y0: 0, x1: x, y1: 1, xref: 'x', yref: 'paper', line: {color: 'red', width: 1, dash: 'dot'}});");
            w.println("  }");
            // Empirical centroids
            w.println("  const empTraces = [];");
            w.println("  snp.centroids.forEach((c, i) => {");
            w.println("    empTraces.push({ x: [c], y: [Math.max(...snp.bins) * 1.05], mode: 'markers+text', text: ['Group ' + (i+1)], textposition: 'top', marker: {size: 12, color: '#1e293b', symbol: 'diamond'}, name: 'Centroid ' + (i+1) });");
            w.println("    shapes.push({type: 'line', x0: c, y0: 0, x1: c, y1: 1, xref: 'x', yref: 'paper', line: {color: '#1e293b', width: 2}});");
            w.println("  });");

            w.println("  const layout = { title: 'Dosage Distribution and Centroids - ' + snp.id, xaxis: {title: 'Dosage Frequency (0.0 - 1.0)', range: [-0.05, 1.05]}, yaxis: {title: 'Individual Count'}, shapes: shapes, showlegend: false, margin: {t: 50} };");
            w.println("  Plotly.react('chart', [traceHist, ...empTraces], layout, {responsive:true});");
            w.println("}");

            w.println("function filterList() {");
            w.println("  const val = document.getElementById('search').value.toLowerCase();");
            w.println("  document.querySelectorAll('.snp-item').forEach(el => {");
            w.println("    const text = el.innerText.toLowerCase();");
            w.println("    el.style.display = text.includes(val) ? 'block' : 'none';");
            w.println("  });");
            w.println("}");

            w.println("if(snps.length > 0) selectSnp(0);");
            w.println("</script></body></html>");
        }
    }

    private static String arrayToString(int[] arr) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < arr.length; i++) {
            sb.append(arr[i]).append(i == arr.length - 1 ? "" : ",");
        }
        return sb.append("]").toString();
    }

    private static String arrayToString(float[] arr) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < arr.length; i++) {
            sb.append(String.format("%.4f", arr[i])).append(i == arr.length - 1 ? "" : ",");
        }
        return sb.append("]").toString();
    }
}
