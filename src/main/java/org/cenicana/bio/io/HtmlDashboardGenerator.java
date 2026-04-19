package org.cenicana.bio.io;
import org.cenicana.bio.utils.FileUtils;
import org.cenicana.bio.core.JoinMapCpFormat;
import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.VcfStatisticsCalculator;
import org.cenicana.bio.core.VcfFilter;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

public class HtmlDashboardGenerator {

	public static void generateReport(VcfStatisticsCalculator stats, String outputPath) throws IOException {
		try (PrintWriter w = new PrintWriter(new FileWriter(outputPath))) {

			// ── Pre-computed scalars ──────────────────────────────────────────
			int numVariants = stats.numSnps + stats.numIndels;
			int numSnpsTotal = numVariants == 0 ? 1 : numVariants;
			double tsTvRatio = stats.numTransversions == 0 ? 0
				: (double) stats.numTransitions / stats.numTransversions;
			String tsTvStr   = String.format("%.2f", tsTvRatio);
			String meanEhStr = String.format("%.4f", stats.meanEH);
			String tajimaDStr = Double.isNaN(stats.tajimaD) ? "N/A"
				: String.format("%.3f", stats.tajimaD);
			double hweViolPct = stats.numHweTested > 0
				? 100.0 * stats.numHweViolations / stats.numHweTested : 0.0;

			// ── JS data arrays for per-sample charts ─────────────────────────
			StringBuilder names   = new StringBuilder("[");
			StringBuilder depths  = new StringBuilder("[");
			StringBuilder missing = new StringBuilder("[");
			StringBuilder obsHet  = new StringBuilder("[");
			StringBuilder fstat   = new StringBuilder("[");
			for (int i = 0; i < stats.sampleNames.length; i++) {
				double dp  = stats.sampleDepthCount[i] == 0 ? 0
					: (double) stats.sampleTotalDepth[i] / stats.sampleDepthCount[i];
				double mis = ((double) stats.sampleMissingCount[i] / numSnpsTotal) * 100.0;
				names.append("'").append(stats.sampleNames[i]).append("'");
				depths.append(String.format("%.2f", dp));
				missing.append(String.format("%.2f", mis));
				obsHet.append(String.format("%.4f", stats.sampleObsHet[i]));
				fstat.append(String.format("%.4f", stats.sampleFstat[i]));
				if (i < stats.sampleNames.length - 1) {
					names.append(","); depths.append(","); missing.append(",");
					obsHet.append(","); fstat.append(",");
				}
			}
			names.append("]"); depths.append("]"); missing.append("]");
			obsHet.append("]"); fstat.append("]");

			// ── Histogram arrays ─────────────────────────────────────────────
			StringBuilder mafJs      = buildArray(stats.mafHistogram);
			StringBuilder siteMissJs = buildArray(stats.siteMissingnessHistogram);
			StringBuilder ehJs       = buildArray(stats.ehHistogram);
			StringBuilder hweJs      = buildArray(stats.hweChiSqHistogram);

			// ── Chromosome density ────────────────────────────────────────────
			StringBuilder chromNames  = new StringBuilder("[");
			StringBuilder chromCounts = new StringBuilder("[");
			boolean first = true;
			for (Map.Entry<String, Integer> e : stats.variantsPerChromosome.entrySet()) {
				if (!first) { chromNames.append(","); chromCounts.append(","); }
				chromNames.append("'").append(e.getKey()).append("'");
				chromCounts.append(e.getValue());
				first = false;
			}
			chromNames.append("]"); chromCounts.append("]");

			// ── Fst heatmap data ──────────────────────────────────────────────
			StringBuilder fstHeatZ    = new StringBuilder("[");
			StringBuilder fstHeatText = new StringBuilder("[");
			StringBuilder fstPopLabels = new StringBuilder("[");
			boolean hasFst = stats.pairwiseFst != null && stats.populationNames != null;
			if (hasFst) {
				int np = stats.populationNames.length;
				for (int i = 0; i < np; i++) {
					fstPopLabels.append("'").append(stats.populationNames[i]).append("'");
					if (i < np - 1) fstPopLabels.append(",");
					fstHeatZ.append("[");
					fstHeatText.append("[");
					for (int j = 0; j < np; j++) {
						double v = (i == j) ? 0 : stats.pairwiseFst[i][j];
						fstHeatZ.append(String.format("%.4f", v));
						fstHeatText.append("'").append(String.format("%.4f", v)).append("'");
						if (j < np - 1) { fstHeatZ.append(","); fstHeatText.append(","); }
					}
					fstHeatZ.append("]");
					fstHeatText.append("]");
					if (i < np - 1) { fstHeatZ.append(","); fstHeatText.append(","); }
				}
			}
			fstHeatZ.append("]"); fstHeatText.append("]"); fstPopLabels.append("]");

			// ── HTML ──────────────────────────────────────────────────────────
			w.println("<!DOCTYPE html><html lang='en'><head>");
			w.println("<meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1.0'>");
			w.println("<title>BioCenicana – VCF QC Dashboard</title>");
			w.println("<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>");
			w.println("<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&display=swap' rel='stylesheet'>");
			w.println("<style>");
			w.println("*{box-sizing:border-box}");
			w.println("body{font-family:'Inter',sans-serif;background:#f0f4f8;color:#0f172a;margin:0;padding:2rem}");
			w.println(".header{text-align:center;margin-bottom:2.5rem}");
			w.println(".header h1{font-weight:800;font-size:2.2rem;color:#1e293b;margin:0 0 .4rem}");
			w.println(".header p{color:#64748b;font-size:1rem;margin:0}");
			w.println(".kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;margin-bottom:2.5rem}");
			w.println(".kpi{background:white;padding:1.2rem 1rem;border-radius:14px;box-shadow:0 2px 8px rgba(0,0,0,.07);text-align:center}");
			w.println(".kpi-t{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:#64748b;font-weight:700;margin-bottom:.4rem}");
			w.println(".kpi-v{font-size:1.8rem;font-weight:800;color:#3b82f6;line-height:1}");
			w.println(".kpi-v.warn{color:#f59e0b}.kpi-v.good{color:#10b981}.kpi-v.info{color:#6366f1}");
			w.println(".section{font-size:1.3rem;font-weight:700;color:#1e293b;margin:2.5rem 0 1rem;border-left:4px solid #3b82f6;padding-left:.8rem}");
			w.println(".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(500px,1fr));gap:1.5rem}");
			w.println(".card{background:white;padding:1.4rem;border-radius:14px;box-shadow:0 2px 8px rgba(0,0,0,.07)}");
			w.println(".card h3{margin:0 0 .8rem;color:#1e293b;font-size:1rem;border-bottom:2px solid #e2e8f0;padding-bottom:.5rem}");
			w.println(".full{grid-column:1/-1}");
			w.println(".badge{display:inline-block;padding:.2rem .7rem;border-radius:999px;font-size:.8rem;font-weight:600;margin:.2rem}");
			w.println(".badge-ok{background:#dcfce7;color:#166534}.badge-warn{background:#fef3c7;color:#92400e}.badge-bad{background:#fee2e2;color:#991b1b}");
			w.println("</style></head><body>");

			// Header
			w.println("<div class='header'><h1>🧬 BioCenicana VCF Quality Control Dashboard</h1>");
			w.println("<p>Level 1 + Level 2 NGSEP-compatible statistics &bull; GATK &bull; FreeBayes &bull; Interactive & publication-ready</p></div>");

			// ── KPI Cards ─────────────────────────────────────────────────────
			w.println("<div class='kpi-grid'>");
			w.printf("<div class='kpi'><div class='kpi-t'>Total Variants</div><div class='kpi-v'>%d</div></div>%n", numVariants);
			w.printf("<div class='kpi'><div class='kpi-t'>SNPs</div><div class='kpi-v'>%d</div></div>%n", stats.numSnps);
			w.printf("<div class='kpi'><div class='kpi-t'>InDels</div><div class='kpi-v'>%d</div></div>%n", stats.numIndels);
			w.printf("<div class='kpi'><div class='kpi-t'>Ts/Tv Ratio</div><div class='kpi-v good'>%s</div></div>%n", tsTvStr);
			w.printf("<div class='kpi'><div class='kpi-t'>Mean EH</div><div class='kpi-v info'>%s</div></div>%n", meanEhStr);
			w.printf("<div class='kpi'><div class='kpi-t'>Samples</div><div class='kpi-v'>%d</div></div>%n", stats.sampleNames.length);
			w.printf("<div class='kpi'><div class='kpi-t'>Tajima's D</div><div class='kpi-v %s'>%s</div></div>%n",
				Double.isNaN(stats.tajimaD) ? "" : (stats.tajimaD < -2 || stats.tajimaD > 2 ? "warn" : "good"),
				tajimaDStr);
			w.printf("<div class='kpi'><div class='kpi-t'>HWE Violations</div><div class='kpi-v %s'>%.1f%%</div></div>%n",
				hweViolPct > 20 ? "warn" : "good", hweViolPct);
			w.printf("<div class='kpi'><div class='kpi-t'>Biallelic SNPs</div><div class='kpi-v'>%d</div></div>%n", stats.numBiallelic);
			w.printf("<div class='kpi'><div class='kpi-t'>Multiallelic</div><div class='kpi-v warn'>%d</div></div>%n", stats.numMultiallelic);
			w.println("</div>");

			// ── Section 1: Per-Sample ─────────────────────────────────────────
			w.println("<div class='section'>📊 Per-Sample Quality Metrics</div>");
			w.println("<div class='grid'>");
			w.println("<div class='card'><h3>Average Sequencing Depth (DP) per Sample</h3><div id='c-depth'></div></div>");
			w.println("<div class='card'><h3>Sample Missingness (%)</h3><div id='c-miss'></div></div>");
			w.println("<div class='card'><h3>Observed Heterozygosity (OH) per Sample</h3><div id='c-oh'></div></div>");
			w.println("<div class='card'><h3>Inbreeding Coefficient F = 1 − OH/EH</h3><div id='c-fstat'></div></div>");
			w.println("</div>");

			// ── Section 2: Variant Level ──────────────────────────────────────
			w.println("<div class='section'>🔬 Variant-Level Statistics (Level 1)</div>");
			w.println("<div class='grid'>");
			w.println("<div class='card'><h3>Minor Allele Frequency (MAF) Spectrum</h3><div id='c-maf'></div></div>");
			w.println("<div class='card'><h3>Expected Heterozygosity (EH) Distribution per Site</h3><div id='c-eh'></div></div>");
			w.println("<div class='card'><h3>Variant Types (Ts / Tv / InDel)</h3><div id='c-types'></div></div>");
			w.println("<div class='card'><h3>Variant Density per Chromosome</h3><div id='c-chrom'></div></div>");
			w.println("<div class='card full'><h3>Site Missingness Distribution</h3><div id='c-sitemiss'></div></div>");
			w.println("</div>");

			// ── Section 3: Level 2 ────────────────────────────────────────────
			w.println("<div class='section'>🧪 Advanced Population Genetics (Level 2)</div>");
			w.println("<div class='grid'>");
			w.println("<div class='card'><h3>HWE Chi-square Distribution & Exact Test (Fisher p&lt;0.05)</h3><div id='c-hwe'></div></div>");
			w.println("<div class='card'><h3>AN – Allele Number per Site</h3><div id='c-an'></div></div>");
			if (hasFst) {
				w.println("<div class='card full'><h3>Pairwise Fst Between Populations</h3><div id='c-fst'></div></div>");
			}
			w.println("</div>");

			// ── JavaScript ────────────────────────────────────────────────────
			w.println("<script>");
			// Utility
			w.println("function S(a){if(!a.length)return{mean:0,sd:0,med:0,max:0};");
			w.println("let s=a.reduce((x,y)=>x+y,0),m=s/a.length;");
			w.println("let sd=Math.sqrt(a.map(x=>Math.pow(x-m,2)).reduce((x,y)=>x+y,0)/a.length);");
			w.println("let so=[...a].sort((a,b)=>a-b);");
			w.println("return{mean:m.toFixed(3),sd:sd.toFixed(3),med:so[Math.floor(so.length/2)].toFixed(3),max:so[so.length-1].toFixed(3)}}");
			w.println("function ann(t){return[{xref:'paper',yref:'paper',x:.98,y:.98,xanchor:'right',yanchor:'top',text:t,showarrow:false,font:{size:12},bordercolor:'#cbd5e1',borderwidth:1,borderpad:5,bgcolor:'rgba(255,255,255,0.9)'}]}");

			// Data
			w.println("var sN=" + names   + ";");
			w.println("var sD=" + depths  + ";");
			w.println("var sM=" + missing + ";");
			w.println("var sOH="+ obsHet  + ";");
			w.println("var sF=" + fstat   + ";");
			w.println("var cfg={responsive:true,displayModeBar:true,toImageButtonOptions:{format:'png',filename:'biocenicana_qc',height:1080,width:1920,scale:2}};");

			// ── Chart 1: Depth
			w.println("(()=>{var ds=S(sD);");
			w.println("Plotly.newPlot('c-depth',[{x:sN,y:sD,type:'bar',marker:{color:'#3b82f6'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Avg. Depth (reads)'},annotations:ann('Mean: '+ds.mean+'<br>SD: '+ds.sd+'<br>Median: '+ds.med)},cfg)})();");

			// ── Chart 2: Missingness scatter
			w.println("(()=>{var ms=S(sM);");
			w.println("Plotly.newPlot('c-miss',[{x:sN,y:sM,type:'scatter',mode:'markers',marker:{size:9,color:'#ef4444'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Missing Genotypes (%)',rangemode:'tozero'},annotations:ann('Mean: '+ms.mean+'%<br>SD: '+ms.sd+'%<br>Max: '+ms.max+'%')},cfg)})();");

			// ── Chart 3: OH
			w.println("(()=>{var ohs=S(sOH);");
			w.println("Plotly.newPlot('c-oh',[{x:sN,y:sOH,type:'bar',marker:{color:'#8b5cf6'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Observed Heterozygosity (OH)'},annotations:ann('Mean: '+ohs.mean+'<br>SD: '+ohs.sd+'<br>Pop EH: " + meanEhStr + "')},cfg)})();");

			// ── Chart 4: F-statistic
			w.println("(()=>{var fs=S(sF);");
			w.println("Plotly.newPlot('c-fstat',[{x:sN,y:sF,type:'bar',marker:{color:'#f59e0b'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Inbreeding Coefficient F'},");
			w.println("shapes:[{type:'line',x0:0,x1:1,xref:'paper',y0:0,y1:0,line:{color:'red',width:1,dash:'dot'}}],");
			w.println("annotations:ann('Mean F: '+fs.mean+'<br>SD: '+fs.sd)},cfg)})();");

			// ── Chart 5: MAF histogram
			w.println("(()=>{var mafL=[]; for(var i=0; i<50; i++) mafL.push((i*0.01).toFixed(2)+'-'+((i+1)*0.01).toFixed(2));");
			w.println("Plotly.newPlot('c-maf',[{x:mafL,y:" + mafJs + ",type:'bar',marker:{color:'#10b981'}}],");
			w.println("{xaxis:{title:'Minor Allele Frequency (MAF)',tickangle:-45,nticks:10},yaxis:{title:'Number of Variants'}},cfg)})();");

			// ── Chart 6: EH histogram
			w.println("(()=>{var ehL=['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0'];");
			w.println("Plotly.newPlot('c-eh',[{x:ehL,y:" + ehJs + ",type:'bar',marker:{color:'#6366f1'}}],");
			w.println("{xaxis:{title:'Expected Heterozygosity (EH) per Site',tickangle:-45},yaxis:{title:'Number of SNPs'}},cfg)})();");

			// ── Chart 7: Types pie
			w.println("Plotly.newPlot('c-types',[{labels:['Transitions (Ts)','Transversions (Tv)','InDels'],values:[" + stats.numTransitions + "," + stats.numTransversions + "," + stats.numIndels + "],type:'pie',hole:.4,marker:{colors:['#6366f1','#f59e0b','#8b5cf6']}}],");
			w.println("{annotations:[{font:{size:14},showarrow:false,text:'Variants',x:.5,y:.5}]},cfg);");

			// ── Chart 8: Chromosome density
			w.println("Plotly.newPlot('c-chrom',[{x:" + chromNames + ",y:" + chromCounts + ",type:'bar',marker:{color:'#0ea5e9'}}],");
			w.println("{xaxis:{title:'Chromosome / Scaffold'},yaxis:{title:'Number of Variants'}},cfg);");

			// ── Chart 9: Site missingness histogram
			w.println("(()=>{var smL=['0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%'];");
			w.println("Plotly.newPlot('c-sitemiss',[{x:smL,y:" + siteMissJs + ",type:'bar',marker:{color:'#f43f5e'}}],");
			w.println("{xaxis:{title:'% Missing Genotypes per SNP',tickangle:-45},yaxis:{title:'Number of SNPs'}},cfg)})();");

			// ── Level 2 Chart 10: HWE chi-square histogram
			w.println("(()=>{var hweL=['χ²: 0-2 (OK)','χ²: 2-3.84 (OK)','χ²: 3.84-7 (p<0.05)','χ²: 7-10 (p<0.01)','χ²: >10 (p<0.001)'];");
			w.println("var hweColors=['#10b981','#84cc16','#f59e0b','#ef4444','#991b1b'];");
			w.println("Plotly.newPlot('c-hwe',[{x:hweL,y:" + hweJs + ",type:'bar',marker:{color:hweColors}}],");
			w.printf("{xaxis:{title:'Chi-square Value (1 df)'},yaxis:{title:'Number of SNPs'},");
			w.printf("title:{text:'HWE Violations: %d / %d sites (%.1f%%)'}},cfg)})();%n",
				stats.numHweViolations, stats.numHweTested, hweViolPct);

			// ── Level 2 Chart 11: AN distribution (monomorphic / bi / multi)
			w.println("Plotly.newPlot('c-an',[{");
			w.println("labels:['Monomorphic','Biallelic','Multiallelic'],");
			w.printf("values:[%d,%d,%d],%n", stats.numMonomorphic, stats.numBiallelic, stats.numMultiallelic);
			w.println("type:'pie',hole:.45,marker:{colors:['#cbd5e1','#3b82f6','#f59e0b']},");
			w.println("textinfo:'label+percent'}],");
			w.println("{annotations:[{font:{size:13},showarrow:false,text:'Allele<br>Types',x:.5,y:.5}]},cfg);");

			// ── Level 2 Chart 12: Fst heatmap (optional)
			if (hasFst) {
				w.println("(()=>{var popL=" + fstPopLabels + ";");
				w.println("Plotly.newPlot('c-fst',[{z:" + fstHeatZ + ",text:" + fstHeatText + ",x:popL,y:popL,type:'heatmap',colorscale:'RdBu',reversescale:true,zmin:0,zmax:0.5,texttemplate:'%{text}'}],");
				w.println("{xaxis:{title:'Population'},yaxis:{title:'Population'},margin:{l:120,b:100}},cfg)})();");
			}

			w.println("</script></body></html>");
		}
	}

	private static StringBuilder buildArray(int[] arr) {
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < arr.length; i++) {
			sb.append(arr[i]);
			if (i < arr.length - 1) sb.append(",");
		}
		return sb.append("]");
	}

    public static void generateLdDecayDashboard(String outputPath, double[] sumR2, long[] countR2, int binSizeBp, int halfDecayDistanceBp, double thresholdR2) throws java.io.IOException {
        StringBuilder xData = new StringBuilder("[");
        StringBuilder yData = new StringBuilder("[");
        
        for (int i = 0; i < sumR2.length; i++) {
            if (countR2[i] > 0) {
                double avgR2 = sumR2[i] / countR2[i];
                int distance = i * binSizeBp;
                xData.append(distance).append(",");
                yData.append(String.format(java.util.Locale.US, "%.5f", avgR2)).append(",");
            }
        }
        
        if (xData.length() > 1) {
            xData.setLength(xData.length() - 1);
            yData.setLength(yData.length() - 1);
        }
        xData.append("]");
        yData.append("]");

        String htmlTemplate = "<!DOCTYPE html>\n" +
                "<html>\n" +
                "<head>\n" +
                "    <title>Linkage Disequilibrium (LD) Decay</title>\n" +
                "    <script src=\"https://cdn.plot.ly/plotly-2.24.1.min.js\"></script>\n" +
                "    <style>\n" +
                "        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: #f8f9fa; }\n" +
                "        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }\n" +
                "        h1 { color: #2c3e50; text-align: center; }\n" +
                "        .kpi-container { display: flex; justify-content: center; gap: 20px; margin-bottom: 10px; }\n" +
                "        .kpi { background: #ecf0f1; padding: 15px 25px; border-radius: 8px; text-align: center; }\n" +
                "        .kpi span { display: block; font-size: 24px; font-weight: bold; color: #2980b9; }\n" +
                "        .chart-container { width: 100%; height: 600px; margin-top: 20px; }\n" +
                "    </style>\n" +
                "</head>\n" +
                "<body>\n" +
                "    <div class=\"container\">\n" +
                "        <h1>LD Decay Dashboard</h1>\n" +
                "        <p style=\"text-align: center; color: #7f8c8d;\">Average squared Pearson correlation (R&sup2;) vs Physical Distance (bp)</p>\n" +
                "        <div class=\"kpi-container\">\n" +
                "            <div class=\"kpi\">Estimated Half-Decay Distance: <span>" + (halfDecayDistanceBp == -1 ? "Not reached" : halfDecayDistanceBp + " bp") + "</span></div>\n" +
                "            <div class=\"kpi\">Half-Decay R&sup2; Threshold: <span>" + String.format(java.util.Locale.US, "%.4f", thresholdR2) + "</span></div>\n" +
                "        </div>\n" +
                "        <div id=\"ldPlot\" class=\"chart-container\"></div>\n" +
                "    </div>\n" +
                "\n" +
                "    <script>\n" +
                "        var trace1 = {\n" +
                "            x: " + xData.toString() + ",\n" +
                "            y: " + yData.toString() + ",\n" +
                "            mode: 'lines+markers',\n" +
                "            name: 'LD Decay',\n" +
                "            line: {color: 'rgb(55, 128, 191)', width: 3, shape: 'spline'},\n" +
                "            marker: {color: 'rgb(55, 128, 191)', size: 6}\n" +
                "        };\n" +
                "\n" +
                "        var layout = {\n" +
                "            title: 'Linkage Disequilibrium Decay (Bin Size = " + binSizeBp + " bp)',\n" +
                "            xaxis: { title: 'Physical Distance (bp)', gridcolor: '#ecf0f1' },\n" +
                "            yaxis: { title: 'Average R&sup2;', range: [0, 1.05], gridcolor: '#ecf0f1' },\n" +
                "            plot_bgcolor: '#ffffff',\n" +
                "            paper_bgcolor: '#ffffff',\n" +
                "            hovermode: 'closest',\n";
                
        if (halfDecayDistanceBp != -1) {
            htmlTemplate += "            shapes: [{\n" +
                    "                type: 'line',\n" +
                    "                x0: " + halfDecayDistanceBp + ",\n" +
                    "                y0: 0,\n" +
                    "                x1: " + halfDecayDistanceBp + ",\n" +
                    "                y1: " + thresholdR2 + ",\n" +
                    "                line: { color: 'rgba(231, 76, 60, 0.8)', width: 2, dash: 'dash' }\n" +
                    "            }, {\n" +
                    "                type: 'line',\n" +
                    "                x0: 0,\n" +
                    "                y0: " + thresholdR2 + ",\n" +
                    "                x1: " + halfDecayDistanceBp + ",\n" +
                    "                y1: " + thresholdR2 + ",\n" +
                    "                line: { color: 'rgba(231, 76, 60, 0.8)', width: 2, dash: 'dash' }\n" +
                    "            }]\n";
        }
        
        htmlTemplate += "        };\n" +
                "\n" +
                "        Plotly.newPlot('ldPlot', [trace1], layout);\n" +
                "    </script>\n" +
                "</body>\n" +
                "</html>";

        try (java.io.PrintWriter pw = new java.io.PrintWriter(new java.io.FileWriter(outputPath))) {
            pw.write(htmlTemplate);
        }
    }
}
