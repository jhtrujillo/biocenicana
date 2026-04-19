package org.cenicana.bio.io;

import org.cenicana.bio.VcfStatisticsCalculator;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class HtmlDashboardGenerator {

	public static void generateReport(VcfStatisticsCalculator stats, String outputPath) throws IOException {
		try (PrintWriter writer = new PrintWriter(new FileWriter(outputPath))) {
			
			writer.println("<!DOCTYPE html>");
			writer.println("<html lang=\"en\">");
			writer.println("<head>");
			writer.println("    <meta charset=\"UTF-8\">");
			writer.println("    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">");
			writer.println("    <title>BioCenicana - VCF Quality Control Dashboard</title>");
			writer.println("    <script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\"></script>");
			writer.println("    <link href=\"https://fonts.googleapis.com/css2?family=Inter:wght@400;600;800&display=swap\" rel=\"stylesheet\">");
			writer.println("    <style>");
			writer.println("        body { font-family: 'Inter', sans-serif; background-color: #f8fafc; color: #0f172a; margin: 0; padding: 2rem; }");
			writer.println("        .header { text-align: center; margin-bottom: 3rem; }");
			writer.println("        .header h1 { font-weight: 800; font-size: 2.5rem; color: #1e293b; margin-bottom: 0.5rem; }");
			writer.println("        .header p { color: #64748b; font-size: 1.1rem; }");
			writer.println("        .kpi-container { display: flex; justify-content: space-around; flex-wrap: wrap; margin-bottom: 3rem; gap: 1rem; }");
			writer.println("        .kpi-card { background: white; padding: 1.5rem; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); text-align: center; flex: 1; min-width: 200px; }");
			writer.println("        .kpi-title { font-size: 0.875rem; text-transform: uppercase; letter-spacing: 0.05em; color: #64748b; margin-bottom: 0.5rem; font-weight: 600; }");
			writer.println("        .kpi-value { font-size: 2.25rem; font-weight: 800; color: #3b82f6; }");
			writer.println("        .charts-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(500px, 1fr)); gap: 2rem; }");
			writer.println("        .chart-card { background: white; padding: 1.5rem; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); }");
			writer.println("        .chart-card h3 { margin-top: 0; color: #1e293b; font-size: 1.25rem; margin-bottom: 1rem; border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; }");
			writer.println("    </style>");
			writer.println("</head>");
			writer.println("<body>");

			writer.println("    <div class=\"header\">");
			writer.println("        <h1>BioCenicana VCF Quality Control Dashboard</h1>");
			writer.println("        <p>Interactive high-resolution statistics ready for scientific publication.</p>");
			writer.println("    </div>");

			// KPIs
			writer.println("    <div class=\"kpi-container\">");
			writer.println("        <div class=\"kpi-card\"><div class=\"kpi-title\">Total Variants</div><div class=\"kpi-value\">" + (stats.numSnps + stats.numIndels) + "</div></div>");
			writer.println("        <div class=\"kpi-card\"><div class=\"kpi-title\">Total SNPs</div><div class=\"kpi-value\">" + stats.numSnps + "</div></div>");
			writer.println("        <div class=\"kpi-card\"><div class=\"kpi-title\">Total InDels</div><div class=\"kpi-value\">" + stats.numIndels + "</div></div>");
			
			double tsTvRatio = stats.numTransversions == 0 ? 0 : (double) stats.numTransitions / stats.numTransversions;
			String tsTvStr = String.format("%.2f", tsTvRatio);
			writer.println("        <div class=\"kpi-card\"><div class=\"kpi-title\">Ts/Tv Ratio</div><div class=\"kpi-value\">" + tsTvStr + "</div></div>");
			writer.println("    </div>");

			writer.println("    <div class=\"charts-grid\">");
			writer.println("        <div class=\"chart-card\"><h3>Average Depth (DP) per Sample</h3><div id=\"chart-depth\"></div></div>");
			writer.println("        <div class=\"chart-card\"><h3>Sample Missingness (%)</h3><div id=\"chart-missing\"></div></div>");
			writer.println("        <div class=\"chart-card\"><h3>Minor Allele Frequency (MAF) Spectrum</h3><div id=\"chart-maf\"></div></div>");
			writer.println("        <div class=\"chart-card\"><h3>Variant Types Distribution</h3><div id=\"chart-types\"></div></div>");
			writer.println("        <div class=\"chart-card\" style=\"grid-column: 1 / -1;\"><h3>Site Missingness Distribution</h3><div id=\"chart-site-missing\"></div></div>");
			writer.println("    </div>");

			// Generate JS Data Arrays
			StringBuilder sampleNamesJs = new StringBuilder("[");
			StringBuilder sampleDepthsJs = new StringBuilder("[");
			StringBuilder sampleMissingJs = new StringBuilder("[");

			int numSnpsTotal = stats.numSnps + stats.numIndels;
			if (numSnpsTotal == 0) numSnpsTotal = 1;

			for (int i = 0; i < stats.sampleNames.length; i++) {
				sampleNamesJs.append("'").append(stats.sampleNames[i]).append("'");
				
				double avgDepth = stats.sampleDepthCount[i] == 0 ? 0 : (double) stats.sampleTotalDepth[i] / stats.sampleDepthCount[i];
				sampleDepthsJs.append(avgDepth);

				double pctMissing = ((double) stats.sampleMissingCount[i] / numSnpsTotal) * 100.0;
				sampleMissingJs.append(pctMissing);

				if (i < stats.sampleNames.length - 1) {
					sampleNamesJs.append(", ");
					sampleDepthsJs.append(", ");
					sampleMissingJs.append(", ");
				}
			}
			sampleNamesJs.append("]");
			sampleDepthsJs.append("]");
			sampleMissingJs.append("]");

			StringBuilder mafJs = new StringBuilder("[");
			StringBuilder siteMissJs = new StringBuilder("[");
			for (int i=0; i<10; i++) {
				mafJs.append(stats.mafHistogram[i]);
				siteMissJs.append(stats.siteMissingnessHistogram[i]);
				if (i < 9) { mafJs.append(", "); siteMissJs.append(", "); }
			}
			mafJs.append("]");
			siteMissJs.append("]");

			writer.println("    <script>");
			writer.println("        var sampleNames = " + sampleNamesJs.toString() + ";");
			writer.println("        var sampleDepths = " + sampleDepthsJs.toString() + ";");
			writer.println("        var sampleMissing = " + sampleMissingJs.toString() + ";");
			
			// Configuration for High Quality Export
			writer.println("        var config = { responsive: true, displayModeBar: true, toImageButtonOptions: { format: 'png', filename: 'biocenicana_chart', height: 1080, width: 1920, scale: 2 } };");

			// Chart 1: Depth
			writer.println("        var trDepth = { x: sampleNames, y: sampleDepths, type: 'bar', marker: { color: '#3b82f6' } };");
			writer.println("        var layDepth = { margin: { b: 120 }, xaxis: { tickangle: -45 } };");
			writer.println("        Plotly.newPlot('chart-depth', [trDepth], layDepth, config);");

			// Chart 2: Sample Missingness
			writer.println("        var trMiss = { x: sampleNames, y: sampleMissing, type: 'scatter', mode: 'markers', marker: { size: 10, color: '#ef4444' } };");
			writer.println("        var layMiss = { margin: { b: 120 }, xaxis: { tickangle: -45 }, yaxis: { title: '% Missing', rangemode: 'tozero' } };");
			writer.println("        Plotly.newPlot('chart-missing', [trMiss], layMiss, config);");

			// Chart 3: MAF Spectrum
			writer.println("        var mafLabels = ['0-0.05', '0.05-0.1', '0.1-0.15', '0.15-0.2', '0.2-0.25', '0.25-0.3', '0.3-0.35', '0.35-0.4', '0.4-0.45', '0.45-0.5'];");
			writer.println("        var trMaf = { x: mafLabels, y: " + mafJs.toString() + ", type: 'bar', marker: { color: '#10b981' } };");
			writer.println("        var layMaf = { xaxis: { title: 'Minor Allele Frequency' }, yaxis: { title: 'Variant Count' } };");
			writer.println("        Plotly.newPlot('chart-maf', [trMaf], layMaf, config);");

			// Chart 4: Variant Types Pie
			writer.println("        var trPie = { labels: ['Transitions (Ts)', 'Transversions (Tv)', 'InDels'], values: [" + stats.numTransitions + ", " + stats.numTransversions + ", " + stats.numIndels + "], type: 'pie', hole: 0.4, marker: { colors: ['#6366f1', '#f59e0b', '#8b5cf6'] } };");
			writer.println("        Plotly.newPlot('chart-types', [trPie], {}, config);");

			// Chart 5: Site Missingness Histogram
			writer.println("        var smLabels = ['0-10%', '10-20%', '20-30%', '30-40%', '40-50%', '50-60%', '60-70%', '70-80%', '80-90%', '90-100%'];");
			writer.println("        var trSiteMiss = { x: smLabels, y: " + siteMissJs.toString() + ", type: 'bar', marker: { color: '#f43f5e' } };");
			writer.println("        var laySiteMiss = { xaxis: { title: 'Percentage of missing genotypes per SNP' }, yaxis: { title: 'Number of SNPs' } };");
			writer.println("        Plotly.newPlot('chart-site-missing', [trSiteMiss], laySiteMiss, config);");

			writer.println("    </script>");
			writer.println("</body>");
			writer.println("</html>");
		}
	}
}
