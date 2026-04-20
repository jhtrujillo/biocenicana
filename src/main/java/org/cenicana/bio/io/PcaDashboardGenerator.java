package org.cenicana.bio.io;

import org.cenicana.bio.core.PopulationStructureAnalyzer.PcaResult;
import java.io.*;
import java.util.Locale;

/**
 * Generates an interactive HTML dashboard for PCA results using Plotly.js.
 */
public class PcaDashboardGenerator {

    public static void generateReport(PcaResult result, String htmlPath) throws IOException {
        try (PrintWriter w = new PrintWriter(new FileWriter(htmlPath))) {
            w.println("<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'>");
            w.println("<title>BioCenicana - PCA Population Structure</title>");
            w.println("<script src='https://cdn.plot.ly/plotly-2.24.1.min.js'></script>");
            w.println("<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap' rel='stylesheet'>");
            w.println("<style>");
            w.println("  body { font-family: 'Inter', sans-serif; background-color: #f8fafc; color: #1e293b; margin: 0; padding: 20px; }");
            w.println("  .container { max-width: 1200px; margin: 0 auto; }");
            w.println("  .header { display: flex; align-items: center; justify-content: space-between; margin-bottom: 30px; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; }");
            w.println("  .header h1 { margin: 0; color: #0f172a; font-size: 24px; }");
            w.println("  .card { background: white; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); padding: 20px; margin-bottom: 20px; }");
            w.println("  .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }");
            w.println("  .full { grid-column: span 2; }");
            w.println("  h3 { margin-top: 0; color: #334155; border-left: 4px solid #6366f1; padding-left: 10px; }");
            w.println("  .stat { font-size: 14px; color: #64748b; }");
            w.println("</style></head><body>");

            w.println("<div class='container'>");
            w.println("<div class='header'>");
            w.println("  <h1>BioCenicana | PCA Population Structure</h1>");
            w.println("  <div style='display:flex; align-items:center; gap:20px;'>");
            w.println("    <div class='stat'>Color by: ");
            w.println("      <select id='colorMode' onchange='updateColoring()' style='padding:5px; border-radius:5px; border:1px solid #cbd5e1;'>");
            w.println("        <option value='kmeans'>K-Means (Optimal K)</option>");
            w.println("        <option value='dbscan'>DBSCAN (Density/Outliers)</option>");
            w.println("      </select>");
            w.println("    </div>");
            w.println("    <div class='stat'>Samples: " + result.sampleNames.length + " | <b>Global Fst: " + String.format("%.4f", result.fstGlobal) + "</b></div>");
            w.println("  </div>");
            w.println("</div>");

            w.println("<div class='grid'>");
            
            // ── Chart 1: PC1 vs PC2
            w.println("<div class='card full'><h3>Principal Component Analysis</h3><div id='pc12'></div></div>");
            
            // ── Chart 2: Scree Plot
            w.println("<div class='card'><h3>Explained Variance (Scree Plot)</h3><div id='scree'></div></div>");

            // ── Chart 4: Elbow Plot
            w.println("<div class='card'><h3>Elbow Method (Optimal K)</h3><div id='elbow'></div></div>");

            // ── Chart 5: Distance Heatmap
            w.println("<div class='card full'><h3>Genetic Distance Matrix (Heatmap)</h3><div id='heatmap'></div></div>");
            
            // ── Chart 3: PC2 vs PC3 (if available)
            if (result.pcMatrix[0].length >= 3) {
                w.println("<div class='card full'><h3>Principal Component Analysis (3D View: PC1, PC2, PC3)</h3><div id='pc3d' style='height: 600px;'></div></div>");
                w.println("<div class='card full'><h3>Principal Component Analysis (PC2 vs PC3)</h3><div id='pc23'></div></div>");
            }
            
            w.println("</div>"); // end grid
            w.println("</div>"); // end container

            w.println("<script>");
            w.println("const cfg = { responsive: true, displaylogo: false };");
            w.println("const colors = ['#6366f1','#f59e0b','#10b981','#ef4444','#8b5cf6','#ec4899','#06b6d4','#f97316','#84cc16','#64748b'];");
            
            // Data Prep
            StringBuilder labels = new StringBuilder("[");
            StringBuilder pc1 = new StringBuilder("[");
            StringBuilder pc2 = new StringBuilder("[");
            StringBuilder pc3 = new StringBuilder("[");
            StringBuilder clusters = new StringBuilder("[");
            for (int i = 0; i < result.sampleNames.length; i++) {
                if (i > 0) { labels.append(","); pc1.append(","); pc2.append(","); pc3.append(","); clusters.append(","); }
                labels.append("'").append(result.sampleNames[i]).append("'");
                pc1.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][0]));
                pc2.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][1]));
                clusters.append(result.clusterAssignments[i]);
                if (result.pcMatrix[0].length >= 3) pc3.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][2]));
            }
            labels.append("]"); pc1.append("]"); pc2.append("]"); pc3.append("]"); clusters.append("]");

            // Distance Matrix Prep
            StringBuilder distData = new StringBuilder("[");
            for (int i = 0; i < result.sampleNames.length; i++) {
                if (i > 0) distData.append(",");
                distData.append("[");
                for (int j = 0; j < result.sampleNames.length; j++) {
                    if (j > 0) distData.append(",");
                    distData.append(String.format(Locale.US, "%.4f", result.distanceMatrix[i][j]));
                }
                distData.append("]");
            }
            distData.append("]");
            
            // Plot Heatmap
            w.println("Plotly.newPlot('heatmap', [{z: " + distData + ", x: labels, y: labels, type: 'heatmap', colorscale: 'Viridis'}], {");
            w.println("  xaxis: { title: 'Sample / Variety', type: 'category', tickangle: -90 },");
            w.println("  yaxis: { title: 'Sample / Variety', type: 'category' },");
            w.println("  margin: { t: 30, l: 150, b: 150 }");
            w.println("}, cfg);");

            // Initial Renders
            w.println("updateColoring();");

            // Plot Scree
            StringBuilder ev = new StringBuilder("[");
            StringBuilder evLabels = new StringBuilder("[");
            for (int j = 0; j < result.explainedVariance.length; j++) {
                if (j > 0) { ev.append(","); evLabels.append(","); }
                ev.append(result.explainedVariance[j] * 100);
                evLabels.append("'PC").append(j+1).append("'");
            }
            ev.append("]"); evLabels.append("]");

            w.println("Plotly.newPlot('scree', [{");
            w.println("  x: " + evLabels + ", y: " + ev + ", type: 'bar', marker: { color: '#0ea5e9' }");
            w.println("}], { xaxis: { title: 'Principal Component' }, yaxis: { title: '% Explained Variance' }, margin: { t: 10 } }, cfg);");

            // Plot Elbow
            StringBuilder wcssVal = new StringBuilder("[");
            StringBuilder wcssK = new StringBuilder("[");
            for (int k = 1; k <= 10; k++) {
                if (k > 1) { wcssVal.append(","); wcssK.append(","); }
                wcssVal.append(result.wcss[k-1]);
                wcssK.append(k);
            }
            wcssVal.append("]"); wcssK.append("]");
            w.println("Plotly.newPlot('elbow', [{x:" + wcssK + ", y:" + wcssVal + ", type:'scatter', mode:'lines+markers', line:{color:'#f43f5e'}}],");
            w.println("{xaxis:{title:'Number of Clusters (K)'}, yaxis:{title:'WCSS (Inertia)'}, margin:{t:10},");
            w.println(" annotations:[{x:" + result.optimalK + ", y:" + result.wcss[result.optimalK-1] + ", text:'Optimal K', showarrow:true, arrowhead:2, ax:0, ay:-40}]");
            w.println("}, cfg);");

            w.println("</script></body></html>");
        }
        System.out.println("[PCA] Dashboard generated: " + htmlPath);
    }
}
