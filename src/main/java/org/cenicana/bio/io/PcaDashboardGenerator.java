package org.cenicana.bio.io;

import org.cenicana.bio.core.PopulationStructureAnalyzer.PcaResult;
import java.io.*;
import java.util.List;
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
            w.println("    <div class='stat'>Visual Mode: ");
            w.println("      <select id='colorMode' onchange='updateColoring()' style='padding:5px; border-radius:5px; border:1px solid #cbd5e1;'>");
            w.println("        <option value='kmeans'>PCA: K-Means</option>");
            w.println("        <option value='dbscan'>PCA: DBSCAN</option>");
            w.println("        <option value='gmm'>PCA: GMM</option>");
            w.println("        <option value='dapc'>DAPC Analysis</option>");
            w.println("      </select>");
            w.println("    </div>");
            w.println("    <div class='stat'>Samples: " + result.sampleNames.length + " | <b>Global Fst: " + String.format("%.4f", result.fstGlobal) + "</b></div>");
            w.println("  </div>");
            w.println("</div>");

            w.println("<div class='grid'>");
            w.println("<div class='card full'><h3>Cluster Visualization (PCA / DAPC)</h3><div id='pc12' style='height: 500px;'></div></div>");
            w.println("<div class='card'><h3>Explained Variance (Scree Plot)</h3><div id='scree' style='height: 400px;'></div></div>");
            w.println("<div class='card'><h3>Elbow Method (Optimal K)</h3><div id='elbow' style='height: 400px;'></div></div>");
            w.println("<div class='card full'><h3>Kinship Analysis (Dendrogram)</h3><div id='dendrogram' style='height: 600px;'></div></div>");
            w.println("<div class='card full'><h3>Genetic Distance Matrix (Heatmap)</h3><div id='heatmap' style='height: 600px;'></div></div>");
            
            if (result.pcMatrix[0].length >= 3) {
                w.println("<div class='card full'><h3>3D Analysis (PC1, PC2, PC3)</h3><div id='pc3d' style='height: 700px;'></div></div>");
            }
            w.println("</div></div>");

            w.println("<script>");
            w.println("const cfg = { responsive: true, displaylogo: false };");
            w.println("const colors = ['#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#46F0F0', '#F032E6', '#BCf60C', '#FABEBE'];");
            
            // Data Prep
            StringBuilder labels = new StringBuilder("[");
            StringBuilder pc1 = new StringBuilder("["); StringBuilder pc2 = new StringBuilder("["); StringBuilder pc3 = new StringBuilder("[");
            StringBuilder ld1 = new StringBuilder("["); StringBuilder ld2 = new StringBuilder("[");
            StringBuilder kmeans = new StringBuilder("["); StringBuilder dbscan = new StringBuilder("["); StringBuilder gmm = new StringBuilder("[");
            
            for (int i = 0; i < result.sampleNames.length; i++) {
                if (i > 0) { labels.append(","); pc1.append(","); pc2.append(","); pc3.append(","); ld1.append(","); ld2.append(","); kmeans.append(","); dbscan.append(","); gmm.append(","); }
                labels.append("'").append(result.sampleNames[i]).append("'");
                pc1.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][0]));
                pc2.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][1]));
                if (result.pcMatrix[0].length >= 3) pc3.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][2])); else pc3.append("0");
                
                ld1.append(String.format(Locale.US, "%.6f", result.dapcMatrix[i][0]));
                if (result.dapcMatrix[i].length > 1) ld2.append(String.format(Locale.US, "%.6f", result.dapcMatrix[i][1])); else ld2.append("0");
                
                kmeans.append(result.clusterAssignments[i]);
                dbscan.append(result.dbscanAssignments[i]);
                gmm.append(result.gmmAssignments[i]);
            }
            labels.append("]"); pc1.append("]"); pc2.append("]"); pc3.append("]"); ld1.append("]"); ld2.append("]"); kmeans.append("]"); dbscan.append("]"); gmm.append("]");

            w.println("const data_labels = " + labels + ";");
            w.println("const data_pc1 = " + pc1 + ";");
            w.println("const data_pc2 = " + pc2 + ";");
            w.println("const data_pc3 = " + pc3 + ";");
            w.println("const data_ld1 = " + ld1 + ";");
            w.println("const data_ld2 = " + ld2 + ";");
            w.println("const data_kmeans = " + kmeans + ";");
            w.println("const data_dbscan = " + dbscan + ";");
            w.println("const data_gmm = " + gmm + ";");
            
            w.println("function updateColoring() {");
            w.println("  const mode = document.getElementById('colorMode').value;");
            w.println("  let xArr = data_pc1, yArr = data_pc2, zArr = data_pc3, assignments;");
            w.println("  if(mode === 'kmeans') assignments = data_kmeans;");
            w.println("  else if(mode === 'dbscan') assignments = data_dbscan;");
            w.println("  else if(mode === 'gmm') assignments = data_gmm;");
            w.println("  else { xArr = data_ld1; yArr = data_ld2; zArr = null; assignments = data_kmeans; }"); // DAPC uses LDs and K-Means groups
            
            w.println("  renderPlot('pc12', xArr, yArr, null, assignments, mode);");
            w.println("  if(document.getElementById('pc3d')) renderPlot('pc3d', xArr, yArr, zArr, assignments, mode);");
            w.println("}");

            w.println("function renderPlot(id, xArr, yArr, zArr, assignments, mode) {");
            w.println("  const unique = [...new Set(assignments)].sort((a,b)=>a-b);");
            w.println("  var traces = unique.map(k => {");
            w.println("    var cx=[], cy=[], cz=[], ct=[];");
            w.println("    for(var i=0; i<assignments.length; i++) {");
            w.println("      if(assignments[i] === k) { cx.push(xArr[i]); cy.push(yArr[i]); if(zArr) cz.push(zArr[i]); ct.push(data_labels[i]); }");
            w.println("    }");
            w.println("    var name = k === -1 ? 'Outliers' : 'Group '+(k+1);");
            w.println("    var color = k === -1 ? '#94a3b8' : colors[Math.abs(k)%10];");
            w.println("    return { x:cx, y:cy, z:zArr?cz:null, text:ct, name:name, type:zArr?'scatter3d':'scatter', mode:'markers', marker:{size:zArr?8:12, color:color, opacity:0.9, line:{width:1.5, color:'white'}}};");
            w.println("  });");
            w.println("  var layout = { margin: {t:30, b:30, l:50, r:50}, hovermode: 'closest', showlegend: true, legend: { orientation: 'h', y: -0.2 } };");
            w.println("  if(id==='pc3d') layout = { margin: {t:0,b:0,l:0,r:0}, scene: { xaxis:{title:mode==='dapc'?'LD1':'PC1'}, yaxis:{title:mode==='dapc'?'LD2':'PC2'}, zaxis:{title:'PC3'} }, legend: { orientation: 'h', y: 0.05 } };");
            w.println("  layout.xaxis = {title: mode==='dapc' ? 'Linear Discriminant 1 (LD1)' : 'Principal Component 1 (PC1)'};");
            w.println("  layout.yaxis = {title: mode==='dapc' ? 'Linear Discriminant 2 (LD2)' : 'Principal Component 2 (PC2)'};");
            w.println("  Plotly.react(id, traces, layout, cfg);");
            w.println("}");

            // Heatmap
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
            w.println("Plotly.newPlot('heatmap', [{z: " + distData + ", x: data_labels, y: data_labels, type: 'heatmap', colorscale: 'Viridis'}], {xaxis:{type:'category'}, yaxis:{type:'category'}, margin:{t:30, l:150, b:150}}, cfg);");

            // Dendrogram
            StringBuilder treeData = new StringBuilder("[");
            for (int i = 0; i < result.treeSegments.size(); i++) {
                double[] s = result.treeSegments.get(i);
                if (i > 0) treeData.append(",");
                treeData.append(String.format(Locale.US, "{x:[%.2f,%.2f,null], y:[%.4f,%.4f,null]}", s[0], s[2], s[1], s[3]));
            }
            treeData.append("]");
            w.println("const treeTraces = " + treeData + ".map(s => ({x:s.x, y:s.y, type:'scatter', mode:'lines', line:{color:'#475569', width:1}, showlegend:false, hoverinfo:'none'}));");
            w.println("treeTraces.push({x:data_labels.map((_,i)=>i), y:data_labels.map(_=>0), mode:'markers', marker:{size:4, color:'#6366f1'}, text:data_labels, hoverinfo:'text', name:'Samples'});");
            w.println("Plotly.newPlot('dendrogram', treeTraces, { xaxis:{title:'Samples (Indices)', showticklabels:false}, yaxis:{title:'Genetic Distance (Height)'}, margin:{t:10} }, cfg);");

            w.println("updateColoring();");

            // Scree
            StringBuilder ev = new StringBuilder("["); StringBuilder evL = new StringBuilder("[");
            for (int j=0; j<result.explainedVariance.length; j++) {
                if(j>0){ ev.append(","); evL.append(","); }
                ev.append(result.explainedVariance[j]*100); evL.append("'PC"+(j+1)+"'");
            }
            ev.append("]"); evL.append("]");
            w.println("Plotly.newPlot('scree', [{x:"+evL+", y:"+ev+", type:'bar', marker:{color:'#0ea5e9'}}], {xaxis:{title:'PC'}, yaxis:{title:'%'}, margin:{t:10}}, cfg);");

            // Elbow
            StringBuilder wVal = new StringBuilder("["); StringBuilder wK = new StringBuilder("[");
            for (int k=1; k<=10; k++) {
                if(k>1){ wVal.append(","); wK.append(","); }
                wVal.append(result.wcss[k-1]); wK.append(k);
            }
            wVal.append("]"); wK.append("]");
            w.println("Plotly.newPlot('elbow', [{x:"+wK+", y:"+wVal+", type:'scatter', mode:'lines+markers', line:{color:'#f43f5e'}}], {xaxis:{title:'K'}, yaxis:{title:'WCSS'}, margin:{t:10}}, cfg);");

            w.println("</script></body></html>");
        }
        System.out.println("[PCA] Dashboard generated: " + htmlPath);
    }
}
