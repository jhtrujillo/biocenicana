package org.cenicana.bio.io;

import org.cenicana.bio.core.PopulationStructureAnalyzer.PcaResult;
import org.cenicana.bio.utils.ResourceUtils;
import java.io.*;
import java.util.List;
import java.util.Locale;

/**
 * Generates an interactive HTML dashboard for PCA results using Plotly.js.
 */
public class PcaDashboardGenerator {

    public static void generateReport(PcaResult result, String htmlPath) throws IOException {
        String plotlyContent = ResourceUtils.loadResource("plotly.min.js");
        
        try (PrintWriter w = new PrintWriter(new FileWriter(htmlPath))) {
            w.println("<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'>");
            w.println("<title>BioJava - PCA Population Structure</title>");
            if (plotlyContent.isEmpty()) {
                w.println("<script src='https://cdn.plot.ly/plotly-2.24.1.min.js'></script>");
            } else {
                w.println("<script>" + plotlyContent + "</script>");
            }
            w.println("<style>");
            w.println("  body { font-family: Inter, system-ui, -apple-system, sans-serif; background-color: #f8fafc; color: #1e293b; margin: 0; padding: 20px; }");
            w.println("  .container { max-width: 1200px; margin: 0 auto; }");
            w.println("  .header { display: flex; align-items: center; justify-content: space-between; margin-bottom: 30px; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; }");
            w.println("  .header h1 { margin: 0; color: #0f172a; font-size: 24px; }");
            w.println("  .card { background: white; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); padding: 20px; margin-bottom: 20px; position: relative; }");
            w.println("  .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }");
            w.println("  .full { grid-column: span 2; }");
            w.println("  h3 { margin-top: 0; color: #334155; border-left: 4px solid #6366f1; padding-left: 10px; display: flex; align-items: center; justify-content: space-between; }");
            w.println("  .info-btn { background: #e2e8f0; color: #64748b; border: none; width: 24px; height: 24px; border-radius: 50%; cursor: pointer; font-size: 14px; font-weight: bold; transition: all 0.2s; }");
            w.println("  .info-btn:hover { background: #6366f1; color: white; }");
            w.println("  .info-text { display: none; background: #f1f5f9; border-left: 3px solid #6366f1; padding: 10px; margin-top: 10px; font-size: 13px; color: #475569; border-radius: 4px; }");
            w.println("  .stat { font-size: 14px; color: #64748b; }");
            w.println("  .summary-bar { display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; margin-bottom: 25px; }");
            w.println("  .summary-item { background: white; padding: 15px; border-radius: 10px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); border-top: 4px solid #6366f1; }");
            w.println("  .summary-label { font-size: 12px; color: #64748b; text-transform: uppercase; letter-spacing: 0.05em; }");
            w.println("  .summary-value { font-size: 20px; font-weight: bold; color: #0f172a; margin-top: 5px; }");
            w.println("  #sampleSearch { padding: 8px 12px; border: 1px solid #cbd5e1; border-radius: 6px; width: 250px; font-size: 14px; outline: none; }");
            w.println("  #sampleSearch:focus { border-color: #6366f1; }");
            w.println("  .btn-export { background: #6366f1; color: white; border: none; padding: 8px 16px; border-radius: 6px; cursor: pointer; font-size: 14px; font-weight: 500; }");
            w.println("  .btn-export:hover { background: #4f46e5; }");
            w.println("</style>");
            w.println("<script>");
            w.println("  function toggleInfo(id) {");
            w.println("    const el = document.getElementById('info-' + id);");
            w.println("    el.style.display = el.style.display === 'block' ? 'none' : 'block';");
            w.println("  }");
            w.println("</script>");
            w.println("</head><body>");

            w.println("<div class='container'>");
            w.println("<div class='header'>");
            w.println("  <h1>BioJava | PCA Population Structure</h1>");
            w.println("  <div class='stat'>Visual Mode: ");
            w.println("    <select id='colorMode' onchange='updateColoring()' style='padding:5px; border-radius:5px; border:1px solid #cbd5e1;'>");
            w.println("      <option value='kmeans'>PCA: K-Means</option>");
            w.println("      <option value='dbscan'>PCA: DBSCAN</option>");
            w.println("      <option value='gmm'>PCA: GMM</option>");
            w.println("      <option value='dapc'>DAPC Analysis</option>");
            w.println("    </select>");
            w.println("  </div>");
            w.println("</div>");
            
            w.println("<div class='summary-bar'>");
            w.println("  <div class='summary-item'><div class='summary-label'>Samples</div><div class='summary-value'>" + result.sampleNames.length + "</div></div>");
            w.println("  <div class='summary-item'><div class='summary-label'>Principal Components</div><div class='summary-value'>" + result.pcMatrix[0].length + "</div></div>");
            w.println("  <div class='summary-item'><div class='summary-label'>Optimal K</div><div class='summary-value'>" + result.optimalK + "</div></div>");
            w.println("  <div class='summary-item'><div class='summary-label'>Global Fst</div><div class='summary-value'>" + String.format("%.4f", result.fstGlobal) + "</div></div>");
            w.println("</div>");

            w.println("<div style='display:flex; justify-content:space-between; align-items:center; margin-bottom:20px;'>");
            w.println("  <input type='text' id='sampleSearch' placeholder='🔍 Search sample name...' onkeyup='handleSearch()'>");
            w.println("  <button class='btn-export' onclick='exportToCSV()'>📥 Export PCA (.csv)</button>");
            w.println("</div>");

            w.println("<div class='grid'>");
            w.println("<div class='card full'><h3>Ancestry Analysis (Admixture Proportions) <button class='info-btn' onclick=\"toggleInfo('ancestry')\">?</button></h3>");
            w.println("<div id='info-ancestry' class='info-text'>Muestra la proporción de ancestría asignada a diferentes grupos (K). Las barras coloreadas indican la probabilidad de pertenencia de cada individuo a los clústeres detectados.</div>");
            w.println("<div id='ancestry' style='height: 400px;'></div></div>");
            w.println("<div class='card full'><h3>Cluster Visualization (PCA / DAPC) <button class='info-btn' onclick=\"toggleInfo('cluster')\">?</button></h3>");
            w.println("<div id='info-cluster' class='info-text'>Visualización en el espacio genético. Los individuos cercanos son genéticamente similares. PCA muestra la variación natural, mientras que DAPC maximiza la separación entre grupos.</div>");
            w.println("<div id='pc12' style='height: 500px;'></div></div>");
            w.println("<div class='card full'><h3>Genomic Relationship Matrix (Kinship - VanRaden) <button class='info-btn' onclick=\"toggleInfo('kinship')\">?</button></h3>");
            w.println("<div id='info-kinship' class='info-text'>Mapa de calor del parentesco genómico. Los colores claros/rojos indican una relación genética estrecha (hermanos o clones), mientras que los oscuros indican individuos no relacionados.</div>");
            w.println("<div id='kinship' style='height: 600px;'></div></div>");
            w.println("<div class='card'><h3>Explained Variance (Scree Plot) <button class='info-btn' onclick=\"toggleInfo('scree')\">?</button></h3>");
            w.println("<div id='info-scree' class='info-text'>Porcentaje de la variabilidad genética total capturado por cada componente principal. Ayuda a determinar cuántos PCs son biológicamente informativos.</div>");
            w.println("<div id='scree' style='height: 400px;'></div></div>");
            w.println("<div class='card'><h3>Elbow Method (Optimal K) <button class='info-btn' onclick=\"toggleInfo('elbow')\">?</button></h3>");
            w.println("<div id='info-elbow' class='info-text'>Gráfico de WCSS vs K. El 'codo' de la curva sugiere el número óptimo de sub-poblaciones reales presentes en el set de datos.</div>");
            w.println("<div id='elbow' style='height: 400px;'></div></div>");
            w.println("<div class='card full'><h3>Kinship Analysis (Dendrogram) <button class='info-btn' onclick=\"toggleInfo('dendro')\">?</button></h3>");
            w.println("<div id='info-dendro' class='info-text'>Representación jerárquica de la estructura familiar basada en la matriz de parentesco. Permite identificar clanes o grupos de familias.</div>");
            w.println("<div id='dendrogram' style='height: 600px;'></div></div>");
            w.println("<div class='card full'><h3>Genetic Distance Matrix (Heatmap) <button class='info-btn' onclick=\"toggleInfo('heatmap')\">?</button></h3>");
            w.println("<div id='info-heatmap' class='info-text'>Mapa de calor de distancias euclidianas. Ideal para detectar muestras muy divergentes (outliers) o duplicados genéticos exactos.</div>");
            w.println("<div id='heatmap' style='height: 600px;'></div></div>");
            
            if (result.pcMatrix[0].length >= 3) {
                w.println("<div class='card full'><h3>3D Analysis (PC1, PC2, PC3) <button class='info-btn' onclick=\"toggleInfo('pc3d')\">?</button></h3>");
                w.println("<div id='info-pc3d' class='info-text'>Exploración tridimensional de los tres primeros componentes principales. Útil para separar grupos complejos que parecen solaparse en 2D.</div>");
                w.println("<div id='pc3d' style='height: 700px;'></div></div>");
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
            StringBuilder ancestry = new StringBuilder("[");
            
            for (int i = 0; i < result.sampleNames.length; i++) {
                if (i > 0) { labels.append(","); pc1.append(","); pc2.append(","); pc3.append(","); ld1.append(","); ld2.append(","); kmeans.append(","); dbscan.append(","); gmm.append(","); ancestry.append(","); }
                labels.append("\"").append(escapeJs(result.sampleNames[i])).append("\"");
                pc1.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][0]));
                pc2.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][1]));
                if (result.pcMatrix[0].length >= 3) pc3.append(String.format(Locale.US, "%.6f", result.pcMatrix[i][2])); else pc3.append("0");
                
                if (result.dapcMatrix[i].length > 0)
                    ld1.append(String.format(Locale.US, "%.6f", result.dapcMatrix[i][0]));
                else 
                    ld1.append("0");
                
                if (result.dapcMatrix[i].length > 1) 
                    ld2.append(String.format(Locale.US, "%.6f", result.dapcMatrix[i][1])); 
                else 
                    ld2.append("0");
                
                kmeans.append(result.clusterAssignments[i]);
                dbscan.append(result.dbscanAssignments[i]);
                gmm.append(result.gmmAssignments[i]);
                
                ancestry.append("[");
                for (int k = 0; k < result.ancestryProportions[i].length; k++) {
                    if (k > 0) ancestry.append(",");
                    ancestry.append(String.format(Locale.US, "%.6f", result.ancestryProportions[i][k]));
                }
                ancestry.append("]");
            }
            labels.append("]"); pc1.append("]"); pc2.append("]"); pc3.append("]"); ld1.append("]"); ld2.append("]"); kmeans.append("]"); dbscan.append("]"); gmm.append("]"); ancestry.append("]");

            w.println("const data_labels = " + labels + ";");
            w.println("const data_pc1 = " + pc1 + ";");
            w.println("const data_pc2 = " + pc2 + ";");
            w.println("const data_pc3 = " + pc3 + ";");
            w.println("const data_ld1 = " + ld1 + ";");
            w.println("const data_ld2 = " + ld2 + ";");
            w.println("const data_kmeans = " + kmeans + ";");
            w.println("const data_dbscan = " + dbscan + ";");
            w.println("const data_gmm = " + gmm + ";");
            w.println("const data_ancestry = " + ancestry + ";");
            
            // Plot Ancestry Barplot
            w.println("function renderAncestry() {");
            w.println("  const K = data_ancestry[0].length;");
            w.println("  const sortedIndices = data_ancestry.map((a,i)=>[a,i]).sort((a,b)=>{");
            w.println("    const maxA = Math.max(...a[0]); const maxB = Math.max(...b[0]);");
            w.println("    const idxA = a[0].indexOf(maxA); const idxB = b[0].indexOf(maxB);");
            w.println("    if(idxA !== idxB) return idxA - idxB; return maxB - maxA;");
            w.println("  }).map(x=>x[1]);");
            w.println("  const traces = [];");
            w.println("  for(let k=0; k<K; k++) {");
            w.println("    traces.push({ x: sortedIndices.map(idx => data_labels[idx]), y: sortedIndices.map(idx => data_ancestry[idx][k]), type:'bar', name:'Ancestry '+(k+1), marker:{color:colors[k%10]} });");
            w.println("  }");
            w.println("  Plotly.newPlot('ancestry', traces, { barmode:'stack', showlegend:false, xaxis:{showticklabels:false}, yaxis:{range:[0,1], title:'Proportion'}, margin:{t:10, b:20} }, cfg);");
            w.println("}");

            w.println("function updateColoring() {");
            w.println("  const mode = document.getElementById('colorMode').value;");
            w.println("  let xArr = data_pc1, yArr = data_pc2, zArr = data_pc3, assignments;");
            w.println("  if(mode === 'kmeans') assignments = data_kmeans;");
            w.println("  else if(mode === 'dbscan') assignments = data_dbscan;");
            w.println("  else if(mode === 'gmm') assignments = data_gmm;");
            w.println("  else { xArr = data_ld1; yArr = data_ld2; zArr = null; assignments = data_kmeans; }");
            
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

            w.println("function handleSearch() {");
            w.println("  const query = document.getElementById('sampleSearch').value.toLowerCase();");
            w.println("  const plots = ['pc12', 'pc3d'].filter(id => document.getElementById(id));");
            w.println("  plots.forEach(id => {");
            w.println("    const plot = document.getElementById(id);");
            w.println("    if (!plot.data) return;");
            w.println("    const update = { 'marker.opacity': [], 'marker.size': [] };");
            w.println("    plot.data.forEach(trace => {");
            w.println("      if (!trace.text) return;");
            w.println("      const opacities = trace.text.map(t => !query || t.toLowerCase().includes(query) ? 0.9 : 0.1);");
            w.println("      const sizes = trace.text.map(t => !query || t.toLowerCase().includes(query) ? (id === 'pc3d' ? 8 : 12) : (id === 'pc3d' ? 3 : 5));");
            w.println("      update['marker.opacity'].push(opacities);");
            w.println("      update['marker.size'].push(sizes);");
            w.println("    });");
            w.println("    Plotly.restyle(id, update);");
            w.println("  });");
            w.println("}");

            w.println("function exportToCSV() {");
            w.println("  let csv = 'Sample,PC1,PC2,PC3,KMeans_Cluster\\n';");
            w.println("  for(let i=0; i<data_labels.length; i++) {");
            w.println("    csv += data_labels[i] + ',' + data_pc1[i] + ',' + data_pc2[i] + ',' + data_pc3[i] + ',' + (data_kmeans[i]+1) + '\\n';");
            w.println("  }");
            w.println("  const blob = new Blob([csv], { type: 'text/csv' });");
            w.println("  const url = window.URL.createObjectURL(blob);");
            w.println("  const a = document.createElement('a');");
            w.println("  a.href = url;");
            w.println("  a.download = 'biojava_pca_results.csv';");
            w.println("  a.click();");
            w.println("}");

            // Heatmap Genetic Distance
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
            w.println("Plotly.newPlot('heatmap', [{z: " + distData + ", x: data_labels, y: data_labels, type: 'heatmap', colorscale: 'Viridis'}], {xaxis:{type:'category'}, yaxis:{type:'category'}, margin:{t:30, l:150, b:150}, title:'Pairwise Genetic Distance (Euclidean)'}, cfg);");

            // Heatmap Kinship
            StringBuilder kinData = new StringBuilder("[");
            for (int i = 0; i < result.sampleNames.length; i++) {
                if (i > 0) kinData.append(",");
                kinData.append("[");
                for (int j = 0; j < result.sampleNames.length; j++) {
                    if (j > 0) kinData.append(",");
                    kinData.append(String.format(Locale.US, "%.4f", result.kinshipMatrix[i][j]));
                }
                kinData.append("]");
            }
            kinData.append("]");
            w.println("Plotly.newPlot('kinship', [{z: " + kinData + ", x: data_labels, y: data_labels, type: 'heatmap', colorscale: 'Hot'}], {xaxis:{type:'category'}, yaxis:{type:'category'}, margin:{t:30, l:150, b:150}, title:'Genomic Relationship Matrix (VanRaden Kinship)'}, cfg);");

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
            w.println("renderAncestry();");

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

    private static String escapeJs(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("'", "\\'")
                .replace("\b", "\\b")
                .replace("\f", "\\f")
                .replace("\n", "\\n")
                .replace("\r", "\\r")
                .replace("\t", "\\t");
    }
}
