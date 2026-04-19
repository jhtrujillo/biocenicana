package org.cenicana.bio.io;

import org.cenicana.bio.VcfStatisticsCalculator;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

public class HtmlDashboardGenerator {

	public static void generateReport(VcfStatisticsCalculator stats, String outputPath) throws IOException {
		try (PrintWriter w = new PrintWriter(new FileWriter(outputPath))) {

			// ── Computed values used in KPI cards ─────────────────────────────
			double tsTvRatio = stats.numTransversions == 0 ? 0
				: (double) stats.numTransitions / stats.numTransversions;
			String tsTvStr = String.format("%.2f", tsTvRatio);
			String meanEhStr = String.format("%.4f", stats.meanEH);
			int numVariants = stats.numSnps + stats.numIndels;
			int numSnpsTotal = numVariants == 0 ? 1 : numVariants;

			// ── Build JS data arrays ──────────────────────────────────────────
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

			// Build histogram arrays
			StringBuilder mafJs      = new StringBuilder("[");
			StringBuilder siteMissJs = new StringBuilder("[");
			StringBuilder ehJs       = new StringBuilder("[");
			for (int i = 0; i < 10; i++) {
				mafJs.append(stats.mafHistogram[i]);
				siteMissJs.append(stats.siteMissingnessHistogram[i]);
				ehJs.append(stats.ehHistogram[i]);
				if (i < 9) { mafJs.append(","); siteMissJs.append(","); ehJs.append(","); }
			}
			mafJs.append("]"); siteMissJs.append("]"); ehJs.append("]");

			// Build chromosome density arrays
			StringBuilder chromNames  = new StringBuilder("[");
			StringBuilder chromCounts = new StringBuilder("[");
			boolean firstChrom = true;
			for (Map.Entry<String, Integer> e : stats.variantsPerChromosome.entrySet()) {
				if (!firstChrom) { chromNames.append(","); chromCounts.append(","); }
				chromNames.append("'").append(e.getKey()).append("'");
				chromCounts.append(e.getValue());
				firstChrom = false;
			}
			chromNames.append("]"); chromCounts.append("]");

			// ── HTML output ───────────────────────────────────────────────────
			w.println("<!DOCTYPE html><html lang='en'><head>");
			w.println("<meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1.0'>");
			w.println("<title>BioCenicana – VCF QC Dashboard</title>");
			w.println("<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>");
			w.println("<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&display=swap' rel='stylesheet'>");
			w.println("<style>");
			w.println("body{font-family:'Inter',sans-serif;background:#f8fafc;color:#0f172a;margin:0;padding:2rem}");
			w.println(".header{text-align:center;margin-bottom:3rem}");
			w.println(".header h1{font-weight:800;font-size:2.5rem;color:#1e293b;margin-bottom:.5rem}");
			w.println(".header p{color:#64748b;font-size:1.1rem}");
			w.println(".kpi-container{display:flex;flex-wrap:wrap;gap:1rem;margin-bottom:3rem;justify-content:space-around}");
			w.println(".kpi-card{background:white;padding:1.5rem;border-radius:12px;box-shadow:0 4px 6px -1px rgb(0 0 0/.1);text-align:center;flex:1;min-width:180px}");
			w.println(".kpi-title{font-size:.8rem;text-transform:uppercase;letter-spacing:.05em;color:#64748b;margin-bottom:.5rem;font-weight:600}");
			w.println(".kpi-value{font-size:2rem;font-weight:800;color:#3b82f6}");
			w.println(".section-title{font-size:1.4rem;font-weight:700;color:#1e293b;margin:2.5rem 0 1rem;border-left:4px solid #3b82f6;padding-left:.8rem}");
			w.println(".charts-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(500px,1fr));gap:2rem}");
			w.println(".chart-card{background:white;padding:1.5rem;border-radius:12px;box-shadow:0 4px 6px -1px rgb(0 0 0/.1)}");
			w.println(".chart-card h3{margin-top:0;color:#1e293b;font-size:1.1rem;margin-bottom:.8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.5rem}");
			w.println("</style></head><body>");

			// ── Header ────────────────────────────────────────────────────────
			w.println("<div class='header'><h1>BioCenicana VCF Quality Control Dashboard</h1>");
			w.println("<p>Interactive statistics • NGSEP / GATK / FreeBayes compatible • High-resolution export</p></div>");

			// ── KPI Cards ─────────────────────────────────────────────────────
			w.println("<div class='kpi-container'>");
			w.println("<div class='kpi-card'><div class='kpi-title'>Total Variants</div><div class='kpi-value'>" + numVariants + "</div></div>");
			w.println("<div class='kpi-card'><div class='kpi-title'>SNPs</div><div class='kpi-value'>" + stats.numSnps + "</div></div>");
			w.println("<div class='kpi-card'><div class='kpi-title'>InDels</div><div class='kpi-value'>" + stats.numIndels + "</div></div>");
			w.println("<div class='kpi-card'><div class='kpi-title'>Ts/Tv Ratio</div><div class='kpi-value'>" + tsTvStr + "</div></div>");
			w.println("<div class='kpi-card'><div class='kpi-title'>Mean Exp. Het. (EH)</div><div class='kpi-value'>" + meanEhStr + "</div></div>");
			w.println("<div class='kpi-card'><div class='kpi-title'>Samples</div><div class='kpi-value'>" + stats.sampleNames.length + "</div></div>");
			w.println("</div>");

			// ── Section 1: Per Sample Quality ─────────────────────────────────
			w.println("<div class='section-title'>Per-Sample Quality Metrics</div>");
			w.println("<div class='charts-grid'>");
			w.println("<div class='chart-card'><h3>Average Sequencing Depth (DP) per Sample</h3><div id='c-depth'></div></div>");
			w.println("<div class='chart-card'><h3>Sample Missingness (%)</h3><div id='c-miss'></div></div>");
			w.println("<div class='chart-card'><h3>Observed Heterozygosity (OH) per Sample</h3><div id='c-oh'></div></div>");
			w.println("<div class='chart-card'><h3>Inbreeding Coefficient F = 1 - OH/EH per Sample</h3><div id='c-fstat'></div></div>");
			w.println("</div>");

			// ── Section 2: Variant-Level Statistics ───────────────────────────
			w.println("<div class='section-title'>Variant-Level Statistics</div>");
			w.println("<div class='charts-grid'>");
			w.println("<div class='chart-card'><h3>Minor Allele Frequency (MAF) Spectrum</h3><div id='c-maf'></div></div>");
			w.println("<div class='chart-card'><h3>Expected Heterozygosity (EH) Distribution per Site</h3><div id='c-eh'></div></div>");
			w.println("<div class='chart-card'><h3>Variant Types Distribution (Ts/Tv/InDel)</h3><div id='c-types'></div></div>");
			w.println("<div class='chart-card'><h3>Variant Density per Chromosome</h3><div id='c-chrom'></div></div>");
			w.println("<div class='chart-card' style='grid-column:1/-1'><h3>Site Missingness Distribution</h3><div id='c-sitemiss'></div></div>");
			w.println("</div>");

			// ── JavaScript ────────────────────────────────────────────────────
			w.println("<script>");
			// Utility function
			w.println("function S(a){if(!a.length)return{mean:0,sd:0,med:0,max:0};");
			w.println("let s=a.reduce((x,y)=>x+y,0),m=s/a.length;");
			w.println("let sd=Math.sqrt(a.map(x=>Math.pow(x-m,2)).reduce((x,y)=>x+y,0)/a.length);");
			w.println("let so=[...a].sort((a,b)=>a-b);");
			w.println("return{mean:m.toFixed(3),sd:sd.toFixed(3),med:so[Math.floor(so.length/2)].toFixed(3),max:so[so.length-1].toFixed(3)}}");

			// Data
			w.println("var sN=" + names   + ";");
			w.println("var sD=" + depths  + ";");
			w.println("var sM=" + missing + ";");
			w.println("var sOH="+ obsHet  + ";");
			w.println("var sF=" + fstat   + ";");
			w.println("var cfg={responsive:true,displayModeBar:true,toImageButtonOptions:{format:'png',filename:'biocenicana_qc',height:1080,width:1920,scale:2}};");

			// Annotation helper
			w.println("function ann(t){return[{xref:'paper',yref:'paper',x:.98,y:.98,xanchor:'right',yanchor:'top',text:t,showarrow:false,font:{size:13},bordercolor:'#cbd5e1',borderwidth:1,borderpad:5,bgcolor:'rgba(255,255,255,0.92)'}]}");

			// Chart 1 – Depth
			w.println("var ds=S(sD);");
			w.println("Plotly.newPlot('c-depth',[{x:sN,y:sD,type:'bar',marker:{color:'#3b82f6'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Avg. Depth (reads)'},annotations:ann('Mean: '+ds.mean+'<br>SD: '+ds.sd+'<br>Median: '+ds.med)},cfg);");

			// Chart 2 – Missingness scatter
			w.println("var ms=S(sM);");
			w.println("Plotly.newPlot('c-miss',[{x:sN,y:sM,type:'scatter',mode:'markers',marker:{size:9,color:'#ef4444'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Missing Genotypes (%)',rangemode:'tozero'},annotations:ann('Mean: '+ms.mean+'%<br>SD: '+ms.sd+'%<br>Max: '+ms.max+'%')},cfg);");

			// Chart 3 – Observed Heterozygosity
			w.println("var ohs=S(sOH);");
			w.println("Plotly.newPlot('c-oh',[{x:sN,y:sOH,type:'bar',marker:{color:'#8b5cf6'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Observed Heterozygosity (OH)'},annotations:ann('Mean: '+ohs.mean+'<br>SD: '+ohs.sd+'<br>EH pop: " + meanEhStr + "')},cfg);");

			// Chart 4 – F-statistic
			w.println("var fs=S(sF);");
			w.println("Plotly.newPlot('c-fstat',[{x:sN,y:sF,type:'bar',marker:{color:'#f59e0b'}}],");
			w.println("{margin:{b:130},xaxis:{title:'Sample ID',tickangle:-45},yaxis:{title:'Inbreeding Coefficient F'},shapes:[{type:'line',x0:0,x1:1,xref:'paper',y0:0,y1:0,line:{color:'red',width:1,dash:'dot'}}],annotations:ann('Mean F: '+fs.mean+'<br>SD: '+fs.sd)},cfg);");

			// Chart 5 – MAF histogram
			w.println("var mafL=['0-0.05','0.05-0.1','0.1-0.15','0.15-0.2','0.2-0.25','0.25-0.3','0.3-0.35','0.35-0.4','0.4-0.45','0.45-0.5'];");
			w.println("Plotly.newPlot('c-maf',[{x:mafL,y:" + mafJs + ",type:'bar',marker:{color:'#10b981'}}],");
			w.println("{xaxis:{title:'Minor Allele Frequency (MAF) Range',tickangle:-45},yaxis:{title:'Number of Variants'}},cfg);");

			// Chart 6 – EH distribution
			w.println("var ehL=['0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0'];");
			w.println("Plotly.newPlot('c-eh',[{x:ehL,y:" + ehJs + ",type:'bar',marker:{color:'#6366f1'}}],");
			w.println("{xaxis:{title:'Expected Heterozygosity (EH) per Site',tickangle:-45},yaxis:{title:'Number of SNPs'}},cfg);");

			// Chart 7 – Variant types pie
			w.println("Plotly.newPlot('c-types',[{labels:['Transitions (Ts)','Transversions (Tv)','InDels'],values:[" + stats.numTransitions + "," + stats.numTransversions + "," + stats.numIndels + "],type:'pie',hole:.4,marker:{colors:['#6366f1','#f59e0b','#8b5cf6']}}],");
			w.println("{annotations:[{font:{size:14},showarrow:false,text:'Variants',x:.5,y:.5}]},cfg);");

			// Chart 8 – Density per chromosome
			w.println("Plotly.newPlot('c-chrom',[{x:" + chromNames + ",y:" + chromCounts + ",type:'bar',marker:{color:'#0ea5e9'}}],");
			w.println("{xaxis:{title:'Chromosome / Scaffold'},yaxis:{title:'Number of Variants'}},cfg);");

			// Chart 9 – Site missingness histogram
			w.println("var smL=['0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%'];");
			w.println("Plotly.newPlot('c-sitemiss',[{x:smL,y:" + siteMissJs + ",type:'bar',marker:{color:'#f43f5e'}}],");
			w.println("{xaxis:{title:'% Missing Genotypes per SNP',tickangle:-45},yaxis:{title:'Number of SNPs'}},cfg);");

			w.println("</script></body></html>");
		}
	}
}
