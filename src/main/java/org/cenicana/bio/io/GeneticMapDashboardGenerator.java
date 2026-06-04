package org.cenicana.bio.io;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Generates a standalone interactive HTML dashboard for a genetic linkage map.
 * Includes: linkage group view (with chromosome coloring), physical vs genetic
 * correlation, LG-chromosome synteny heatmap, and statistics table.
 */
public class GeneticMapDashboardGenerator {

    public static class MapMarker {
        public String id;
        public String lg;
        public double posCm;
        public String chrPhys;
        public long posPhys;

        public MapMarker(String id, String lg, double posCm, String chrPhys, long posPhys) {
            this.id = id; this.lg = lg; this.posCm = posCm;
            this.chrPhys = chrPhys; this.posPhys = posPhys;
        }
    }

    public static List<MapMarker> readMap(String mapPath) throws IOException {
        List<MapMarker> markers = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(mapPath))) {
            String line; boolean header = true;
            while ((line = br.readLine()) != null) {
                if (header) { header = false; continue; }
                String[] cols = line.split("\t");
                if (cols.length < 5) continue;
                try {
                    markers.add(new MapMarker(cols[0], cols[1],
                        Double.parseDouble(cols[2]), cols[3], Long.parseLong(cols[4])));
                } catch (NumberFormatException ignored) {}
            }
        }
        return markers;
    }

    /** Gene entry loaded from genes_en_mapa.tsv */
    public static class GeneEntry {
        public String id, lg, chr, funcion, categoria;
        public double posCm;
        public GeneEntry(String id, String lg, double posCm, String chr, String funcion, String categoria) {
            this.id = id; this.lg = lg; this.posCm = posCm;
            this.chr = chr; this.funcion = funcion; this.categoria = categoria;
        }
    }

    public static List<GeneEntry> readGenes(String genesPath) throws IOException {
        List<GeneEntry> genes = new ArrayList<>();
        if (genesPath == null || genesPath.isEmpty()) return genes;
        Path p = Paths.get(genesPath);
        if (!Files.exists(p)) { System.err.println("[MapViz] Genes file not found: " + genesPath); return genes; }
        try (BufferedReader br = Files.newBufferedReader(p)) {
            String line; boolean header = true;
            while ((line = br.readLine()) != null) {
                if (header) { header = false; continue; }
                String[] c = line.split("\t");
                if (c.length < 8) continue;
                try {
                    String id       = c[0];
                    String chr      = c[1];
                    String lg       = c.length > 5 ? c[5] : "N/A";
                    double cm       = c.length > 6 && !c[6].equals("N/A") ? Double.parseDouble(c[6]) : -1;
                    String funcion  = c.length > 14 ? c[14] : "";
                    String categoria = c.length > 15 ? c[15] : "";
                    if (!lg.equals("N/A") && cm >= 0)
                        genes.add(new GeneEntry(id, lg, cm, chr, funcion, categoria));
                } catch (NumberFormatException ignored) {}
            }
        }
        System.out.println("[MapViz] " + genes.size() + " genes loaded for overlay.");
        return genes;
    }

    public static void generate(String mapPath, String outputHtmlPath, String genesMapFile) throws IOException {
        List<MapMarker> markers = readMap(mapPath);
        if (markers.isEmpty()) { System.err.println("[MapViz] No markers found."); return; }
        List<GeneEntry> genes = readGenes(genesMapFile);
        Files.writeString(Paths.get(outputHtmlPath), buildHtml(markers, genes, mapPath));
        System.out.println("[MapViz] Interactive map dashboard written to: " + outputHtmlPath);
    }

    public static void generate(String mapPath, String outputHtmlPath) throws IOException {
        generate(mapPath, outputHtmlPath, null);
    }

    private static String buildHtml(List<MapMarker> markers, List<GeneEntry> genes, String mapPath) {
        Map<String, List<MapMarker>> byLg = new LinkedHashMap<>();
        for (MapMarker m : markers) byLg.computeIfAbsent(m.lg, k -> new ArrayList<>()).add(m);

        List<String> lgsSorted = byLg.keySet().stream()
            .sorted(Comparator.comparingInt(lg -> {
                try { return Integer.parseInt(lg.replaceAll("[^0-9]", "")); }
                catch (NumberFormatException e) { return Integer.MAX_VALUE; }
            })).collect(Collectors.toList());

        // Chromosome composition per LG: dominant chr, percentage, all chrs
        Set<String> allChroms = new TreeSet<>(Comparator.comparingInt(c -> {
            try { return Integer.parseInt(c); } catch (Exception e) { return 99; }
        }));
        for (MapMarker m : markers) allChroms.add(m.chrPhys);
        List<String> chromList = new ArrayList<>(allChroms);

        // Summary stats
        int totalMarkers = markers.size();
        int totalLgs = lgsSorted.size();
        double totalCm = 0;
        for (String lg : lgsSorted) {
            double max = byLg.get(lg).stream().mapToDouble(m -> m.posCm).max().orElse(0);
            totalCm += max;
        }
        double avgLgLength = totalLgs > 0 ? totalCm / totalLgs : 0;

        String markersJson   = buildMarkersJson(markers);
        String lgStatsJson   = buildLgStatsJson(lgsSorted, byLg);
        String heatmapJson   = buildHeatmapJson(lgsSorted, chromList, byLg);
        String chromListJson = chromList.stream().map(c -> "\"" + c + "\"")
                                .collect(Collectors.joining(",", "[", "]"));
        String genesJson     = buildGenesJson(genes);

        return "<!DOCTYPE html>\n<html lang='es'>\n<head>\n" +
            "<meta charset='UTF-8'>\n<meta name='viewport' content='width=device-width,initial-scale=1'>\n" +
            "<title>BioJava — Mapa Genético</title>\n" +
            "<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>\n" +
            "<style>" + buildCss() + "</style>\n</head>\n<body>\n" +
            "<div class='header'>\n" +
            "  <div class='header-title'>🧬 Mapa Genético Interactivo</div>\n" +
            "  <div class='header-sub'>BioJava · " + Paths.get(mapPath).getFileName() + "</div>\n" +
            "</div>\n" +
            "<div class='stats-bar'>\n" +
            stat(totalMarkers + "", "Marcadores") +
            stat(totalLgs + "", "Grupos de Ligamiento") +
            stat(chromList.size() + "", "Cromosomas") +
            stat(String.format("%.0f", totalCm), "cM totales") +
            stat(String.format("%.1f", avgLgLength), "cM promedio/LG") +
            "</div>\n" +
            // Chromosome legend
            "<div class='chr-legend' id='chrLegend'></div>\n" +
            (genes.isEmpty() ? "" :
            "<div class='gene-legend'>" +
            "<span>Genes candidatos:</span>" +
            "<div class='gene-legend-item'><div class='gene-dot' style='background:#e74c3c'></div>Selección positiva</div>" +
            "<div class='gene-legend-item'><div class='gene-dot' style='background:#f7c948'></div>Neutral</div>" +
            "<div class='gene-legend-item'><div class='gene-dot' style='background:#58d68d'></div>Purificadora</div>" +
            "</div>\n") +
            "<div class='tabs'>\n" +
            "  <button class='tab active' onclick='showTab(\"lg-view\",this)'>📊 Grupos de Ligamiento</button>\n" +
            "  <button class='tab' onclick='showTab(\"genes-view\",this)'>🧬 Genes Candidatos</button>\n" +
            "  <button class='tab' onclick='showTab(\"corr-view\",this)'>📈 Físico vs Genético</button>\n" +
            "  <button class='tab' onclick='showTab(\"heat-view\",this)'>🗺️ Sintenia LG–Cromosoma</button>\n" +
            "  <button class='tab' onclick='showTab(\"stats-view\",this)'>📋 Estadísticas</button>\n" +
            "</div>\n" +
            "<div id='lg-view' class='tab-content active'>\n" +
            "  <div class='controls'>\n" +
            "    <label>Filtrar por cromosoma: <select id='chrFilter' onchange='renderLgMap()'>\n" +
            "      <option value='all'>Todos</option>" +
            chromList.stream().map(c -> "<option value='" + c + "'>chr" + c + "</option>")
                     .collect(Collectors.joining()) + "\n" +
            "    </select></label>\n" +
            "    <label style='margin-left:16px'>Escala: <input type='range' id='cmScale' min='1' max='15' value='5' " +
            "oninput='document.getElementById(\"cmScaleVal\").textContent=this.value;renderLgMap()'>" +
            " <span id='cmScaleVal'>5</span> px/cM</label>\n" +
            "    <label style='margin-left:16px'><input type='checkbox' id='showChrLabel' checked onchange='renderLgMap()'>" +
            " Mostrar cromosoma dominante</label>\n" +
            (genes.isEmpty() ? "" :
            "    <label style='margin-left:16px'><input type='checkbox' id='showGeneNames' checked onchange='renderLgMap()'>" +
            " Mostrar nombres de genes</label>\n") +
            "  </div>\n" +
            "  <div id='lg-canvas-wrap'><div id='lg-canvas'></div></div>\n" +
            "</div>\n" +
            "<div id='genes-view' class='tab-content'><div id='genes-plot' style='height:500px'></div><div id='genes-table-wrap'></div></div>\n" +
            "<div id='corr-view' class='tab-content'><div id='corr-plot' style='height:600px'></div></div>\n" +
            "<div id='heat-view' class='tab-content'><div id='heat-plot' style='height:600px'></div></div>\n" +
            "<div id='stats-view' class='tab-content'>" +
            "  <div id='stats-plot' style='height:460px'></div>" +
            "  <div id='stats-table-wrap'></div>" +
            "</div>\n" +
            "<div id='tooltip' class='tooltip'></div>\n" +
            "<script>\n" +
            "const MARKERS=" + markersJson + ";\n" +
            "const LG_STATS=" + lgStatsJson + ";\n" +
            "const HEATMAP=" + heatmapJson + ";\n" +
            "const CHROMS=" + chromListJson + ";\n" +
            "const GENES=" + genesJson + ";\n" +
            buildJs() +
            "\n</script>\n</body>\n</html>";
    }

    private static String stat(String val, String lbl) {
        return "<div class='stat-card'><div class='stat-val'>" + val +
               "</div><div class='stat-lbl'>" + lbl + "</div></div>\n";
    }

    // ── JSON builders ──────────────────────────────────────────────────────

    private static String buildMarkersJson(List<MapMarker> markers) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < markers.size(); i++) {
            MapMarker m = markers.get(i);
            if (i > 0) sb.append(",");
            sb.append("{\"id\":\"").append(esc(m.id)).append("\",")
              .append("\"lg\":\"").append(esc(m.lg)).append("\",")
              .append("\"cm\":").append(m.posCm).append(",")
              .append("\"chr\":\"").append(esc(m.chrPhys)).append("\",")
              .append("\"pos\":").append(m.posPhys).append("}");
        }
        return sb.append("]").toString();
    }

    private static String buildLgStatsJson(List<String> lgs, Map<String, List<MapMarker>> byLg) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < lgs.size(); i++) {
            String lg = lgs.get(i);
            List<MapMarker> ms = byLg.get(lg);
            if (i > 0) sb.append(",");
            double length = ms.stream().mapToDouble(m -> m.posCm).max().orElse(0);
            // chromosome counts
            Map<String, Long> chrCounts = ms.stream()
                .collect(Collectors.groupingBy(m -> m.chrPhys, Collectors.counting()));
            String dominant = chrCounts.entrySet().stream()
                .max(Map.Entry.comparingByValue()).map(Map.Entry::getKey).orElse("?");
            long domCount = chrCounts.getOrDefault(dominant, 0L);
            double domPct = ms.isEmpty() ? 0 : 100.0 * domCount / ms.size();
            // build chrMap JSON
            String chrMapJson = chrCounts.entrySet().stream()
                .map(e -> "\"" + e.getKey() + "\":" + e.getValue())
                .collect(Collectors.joining(",", "{", "}"));
            sb.append("{\"lg\":\"").append(lg).append("\",")
              .append("\"length\":").append(String.format("%.2f", length)).append(",")
              .append("\"count\":").append(ms.size()).append(",")
              .append("\"dominant\":\"").append(dominant).append("\",")
              .append("\"domPct\":").append(String.format("%.0f", domPct)).append(",")
              .append("\"chrMap\":").append(chrMapJson).append("}");
        }
        return sb.append("]").toString();
    }

    private static String buildGenesJson(List<GeneEntry> genes) {
        if (genes.isEmpty()) return "[]";
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < genes.size(); i++) {
            GeneEntry g = genes.get(i);
            if (i > 0) sb.append(",");
            sb.append("{\"id\":\"").append(esc(g.id)).append("\",")
              .append("\"lg\":\"").append(esc(g.lg)).append("\",")
              .append("\"cm\":").append(g.posCm).append(",")
              .append("\"chr\":\"").append(esc(g.chr)).append("\",")
              .append("\"funcion\":\"").append(esc(g.funcion.replaceAll("%2C", ","))).append("\",")
              .append("\"categoria\":\"").append(esc(g.categoria)).append("\"}");
        }
        return sb.append("]").toString();
    }

    private static String buildHeatmapJson(List<String> lgs, List<String> chroms,
                                            Map<String, List<MapMarker>> byLg) {
        StringBuilder sb = new StringBuilder("{\"lgs\":[");
        sb.append(lgs.stream().map(l -> "\"" + l + "\"").collect(Collectors.joining(",")));
        sb.append("],\"chroms\":[");
        sb.append(chroms.stream().map(c -> "\"chr" + esc(c) + "\"").collect(Collectors.joining(",")));
        sb.append("],\"matrix\":[");
        for (int i = 0; i < lgs.size(); i++) {
            if (i > 0) sb.append(",");
            sb.append("[");
            List<MapMarker> ms = byLg.get(lgs.get(i));
            for (int j = 0; j < chroms.size(); j++) {
                if (j > 0) sb.append(",");
                final String chr = chroms.get(j);
                long count = ms.stream().filter(m -> m.chrPhys.equals(chr)).count();
                sb.append(count);
            }
            sb.append("]");
        }
        sb.append("]}");
        return sb.toString();
    }

    private static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }

    // ── CSS ────────────────────────────────────────────────────────────────

    private static String buildCss() {
        return "*{box-sizing:border-box;margin:0;padding:0}" +
            "body{font-family:'Segoe UI',Arial,sans-serif;background:#0f1117;color:#e0e0e0;min-height:100vh}" +
            ".header{background:linear-gradient(135deg,#1a1d2e,#16213e);padding:20px 30px;border-bottom:2px solid #2a3f6f}" +
            ".header-title{font-size:24px;font-weight:700;color:#7eb8f7}" +
            ".header-sub{font-size:13px;color:#8899aa;margin-top:4px}" +
            ".stats-bar{display:flex;gap:16px;padding:16px 30px;background:#13162a;flex-wrap:wrap}" +
            ".stat-card{background:#1e2235;border:1px solid #2a3f6f;border-radius:10px;padding:14px 20px;min-width:120px;text-align:center}" +
            ".stat-val{font-size:26px;font-weight:700;color:#7eb8f7}" +
            ".stat-lbl{font-size:11px;color:#8899aa;margin-top:4px;text-transform:uppercase;letter-spacing:.5px}" +
            ".chr-legend{display:flex;gap:8px;padding:10px 30px;background:#13162a;flex-wrap:wrap;border-bottom:1px solid #1e2235}" +
            ".chr-chip{display:inline-flex;align-items:center;gap:6px;padding:4px 12px;border-radius:20px;font-size:12px;cursor:pointer;border:1px solid transparent;transition:all .2s}" +
            ".chr-chip.active{border-color:#fff4}" +
            ".chr-dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}" +
            ".tabs{display:flex;gap:4px;padding:12px 30px 0;background:#13162a;border-bottom:2px solid #2a3f6f}" +
            ".tab{background:#1e2235;border:1px solid #2a3f6f;border-bottom:none;border-radius:8px 8px 0 0;padding:9px 20px;cursor:pointer;color:#8899aa;font-size:13px;transition:all .2s}" +
            ".tab:hover{background:#253050;color:#c0d8f8}" +
            ".tab.active{background:#0f1117;color:#7eb8f7;font-weight:600}" +
            ".tab-content{display:none;padding:24px 30px}" +
            ".tab-content.active{display:block}" +
            ".controls{display:flex;align-items:center;gap:16px;margin-bottom:16px;font-size:13px;color:#8899aa;flex-wrap:wrap}" +
            ".controls select,.controls input[type=number]{background:#1e2235;border:1px solid #2a3f6f;border-radius:6px;color:#e0e0e0;padding:4px 8px}" +
            ".controls input[type=range]{accent-color:#7eb8f7}" +
            "#lg-canvas-wrap{overflow:auto;max-height:72vh}" +
            "#lg-canvas{display:flex;gap:12px;align-items:flex-start;padding:10px;min-width:max-content}" +
            ".lg-col{display:flex;flex-direction:column;align-items:center;min-width:28px}" +
            ".lg-label{font-size:10px;font-weight:600;color:#7eb8f7;margin-bottom:4px;white-space:nowrap;text-align:center}" +
            ".lg-chr-badge{font-size:9px;color:#aabbcc;margin-bottom:3px;white-space:nowrap}" +
            ".lg-bar-wrap{position:relative;width:18px;border-radius:4px;background:#1e2235;border:1px solid #2a3f6f}" +
            ".lg-marker{position:absolute;left:-1px;right:-1px;height:3px;border-radius:1px;cursor:pointer;transition:opacity .15s}" +
            ".lg-marker:hover{opacity:1!important;outline:1px solid #fff4;z-index:10}" +
            ".lg-len{font-size:9px;color:#556677;margin-top:3px}" +
            ".tooltip{position:fixed;background:#1e2235;border:1px solid #2a3f6f;border-radius:8px;padding:10px 14px;font-size:12px;pointer-events:none;display:none;z-index:1000;max-width:270px;line-height:1.6}" +
            ".tooltip b{color:#7eb8f7}" +
            "table{width:100%;border-collapse:collapse;font-size:13px;margin-top:16px}" +
            "th{background:#1e2235;color:#7eb8f7;padding:10px 14px;text-align:left;border-bottom:2px solid #2a3f6f}" +
            "td{padding:8px 14px;border-bottom:1px solid #1e2235;color:#c0d0e0}" +
            "tr:hover td{background:#1a1f35}" +
            ".bar-inline{display:inline-block;height:10px;border-radius:3px;vertical-align:middle;margin-right:6px}" +
            ".gene-legend{display:flex;gap:16px;padding:8px 30px;background:#13162a;font-size:12px;color:#8899aa;border-bottom:1px solid #1e2235}" +
            ".gene-legend-item{display:flex;align-items:center;gap:6px}" +
            ".gene-dot{width:12px;height:4px;border-radius:2px}";
    }

    // ── JavaScript ─────────────────────────────────────────────────────────

    private static String buildJs() {
        return "\n// Palette: 10 chromosomes\n" +
            "const CHR_COLORS={\n" +
            "  '1':'#7eb8f7','2':'#f7c948','3':'#58d68d','4':'#e74c3c','5':'#a569bd',\n" +
            "  '6':'#45b39d','7':'#f39c12','8':'#ec407a','9':'#26c6da','10':'#ffb74d'\n" +
            "};\n" +
            "function chrColor(c){return CHR_COLORS[c]||'#8899aa';}\n\n" +

            "// Build legend\n" +
            "function buildLegend(){\n" +
            "  const div=document.getElementById('chrLegend');\n" +
            "  CHROMS.forEach(c=>{\n" +
            "    const chip=document.createElement('div'); chip.className='chr-chip active';\n" +
            "    chip.innerHTML='<div class=\"chr-dot\" style=\"background:'+chrColor(c)+'\"></div>chr'+c;\n" +
            "    chip.dataset.chr=c;\n" +
            "    chip.onclick=()=>{chip.classList.toggle('active');renderLgMap();};\n" +
            "    div.appendChild(chip);\n" +
            "  });\n" +
            "}\n\n" +

            "// Tab switching\n" +
            "function showTab(id,btn){\n" +
            "  document.querySelectorAll('.tab-content').forEach(e=>e.classList.remove('active'));\n" +
            "  document.querySelectorAll('.tab').forEach(e=>e.classList.remove('active'));\n" +
            "  document.getElementById(id).classList.add('active'); btn.classList.add('active');\n" +
            "  if(id==='corr-view') renderCorr();\n" +
            "  if(id==='heat-view') renderHeatmap();\n" +
            "  if(id==='stats-view') renderStats();\n" +
            "}\n\n" +

            "// ── LG Map ──────────────────────────────────────────────────\n" +
            "function renderLgMap(){\n" +
            "  const scale=parseInt(document.getElementById('cmScale').value)||5;\n" +
            "  const showBadge=document.getElementById('showChrLabel').checked;\n" +
            "  const chrFilter=document.getElementById('chrFilter').value;\n" +
            "  // Active chroms from legend chips\n" +
            "  const activeChips=new Set([...document.querySelectorAll('.chr-chip.active')].map(c=>c.dataset.chr));\n" +
            "  const canvas=document.getElementById('lg-canvas'); canvas.innerHTML='';\n" +
            "  const byLg={};\n" +
            "  MARKERS.forEach(m=>{if(!byLg[m.lg])byLg[m.lg]=[]; byLg[m.lg].push(m);});\n" +
            "  LG_STATS.forEach(stat=>{\n" +
            "    const ms=byLg[stat.lg]||[];\n" +
            "    // filter by chr dropdown\n" +
            "    if(chrFilter!=='all' && stat.dominant!==chrFilter) return;\n" +
            "    // filter by legend chips\n" +
            "    if(!activeChips.has(stat.dominant)) return;\n" +
            "    const maxCm=stat.length; const barH=Math.max(20,maxCm*scale);\n" +
            "    const col=document.createElement('div'); col.className='lg-col';\n" +
            "    const lbl=document.createElement('div'); lbl.className='lg-label';\n" +
            "    lbl.textContent=stat.lg; col.appendChild(lbl);\n" +
            "    if(showBadge){\n" +
            "      const badge=document.createElement('div'); badge.className='lg-chr-badge';\n" +
            "      badge.innerHTML='<span style=\"color:'+chrColor(stat.dominant)+'\">chr'+stat.dominant+'</span> '+stat.domPct+'%';\n" +
            "      col.appendChild(badge);\n" +
            "    }\n" +
            "    const wrap=document.createElement('div'); wrap.className='lg-bar-wrap';\n" +
            "    wrap.style.height=barH+'px';\n" +
            "    // color bar background with dominant chr color\n" +
            "    wrap.style.borderColor=chrColor(stat.dominant)+'66';\n" +
            "    ms.forEach(m=>{\n" +
            "      const pct=maxCm>0?(m.cm/maxCm)*100:0;\n" +
            "      const tick=document.createElement('div'); tick.className='lg-marker';\n" +
            "      tick.style.top=pct+'%';\n" +
            "      tick.style.background=chrColor(m.chr);\n" +
            "      tick.style.opacity='0.85';\n" +
            "      tick.addEventListener('mousemove',e=>showTip(e,m,stat));\n" +
            "      tick.addEventListener('mouseleave',hideTip);\n" +
            "      wrap.appendChild(tick);\n" +
            "    });\n" +
            "    col.appendChild(wrap);\n" +
            "    const len=document.createElement('div'); len.className='lg-len';\n" +
            "    len.textContent=maxCm.toFixed(1)+' cM'; col.appendChild(len);\n" +
            "    canvas.appendChild(col);\n" +
            "  });\n" +
            "  if(!canvas.children.length) canvas.innerHTML='<p style=\"color:#556677;padding:20px\">Sin grupos con ese filtro.</p>';\n" +
            "}\n\n" +

            "// Tooltip\n" +
            "const tip=document.getElementById('tooltip');\n" +
            "function showTip(e,m,stat){\n" +
            "  const chrStr=Object.entries(stat.chrMap).sort((a,b)=>b[1]-a[1])\n" +
            "    .map(([c,n])=>'<span style=\"color:'+chrColor(c)+'\">chr'+c+'</span>: '+n).join(', ');\n" +
            "  tip.innerHTML='<b>'+m.id+'</b><br>'+stat.lg+' · <span style=\"color:'+chrColor(m.chr)+'\">chr'+m.chr+'</span><br>'\n" +
            "    +'Pos: '+m.cm.toFixed(2)+' cM<br>Phys: '+m.pos.toLocaleString()+' bp<br>'\n" +
            "    +'<small>Cromosomas del LG: '+chrStr+'</small>';\n" +
            "  tip.style.display='block';\n" +
            "  tip.style.left=(e.clientX+14)+'px'; tip.style.top=(e.clientY-10)+'px';\n" +
            "}\n" +
            "function hideTip(){tip.style.display='none';}\n\n" +

            "// ── Correlation plot ────────────────────────────────────────\n" +
            "let corrDone=false;\n" +
            "function renderCorr(){\n" +
            "  if(corrDone) return; corrDone=true;\n" +
            "  const traces=CHROMS.map(chr=>{\n" +
            "    const ms=MARKERS.filter(m=>m.chr===chr);\n" +
            "    return{x:ms.map(m=>m.pos/1e6),y:ms.map(m=>m.cm),mode:'markers',type:'scatter',\n" +
            "      name:'chr'+chr,marker:{color:chrColor(chr),size:5,opacity:.75},\n" +
            "      text:ms.map(m=>m.id+'<br>'+m.lg+'<br>'+m.cm.toFixed(2)+' cM'),\n" +
            "      hovertemplate:'%{text}<extra></extra>'};\n" +
            "  });\n" +
            "  Plotly.newPlot('corr-plot',traces,{\n" +
            "    paper_bgcolor:'#0f1117',plot_bgcolor:'#13162a',\n" +
            "    font:{color:'#c0d0e0',size:12},\n" +
            "    xaxis:{title:'Posición física (Mb)',gridcolor:'#1e2235',zerolinecolor:'#2a3f6f'},\n" +
            "    yaxis:{title:'Posición genética (cM)',gridcolor:'#1e2235',zerolinecolor:'#2a3f6f'},\n" +
            "    legend:{bgcolor:'#1e2235',bordercolor:'#2a3f6f',borderwidth:1},\n" +
            "    title:{text:'Correlación posición física vs genética por cromosoma',font:{color:'#7eb8f7'}},\n" +
            "    margin:{t:50,r:30,b:60,l:70}\n" +
            "  },{responsive:true});\n" +
            "}\n\n" +

            "// ── Heatmap ─────────────────────────────────────────────────\n" +
            "let heatDone=false;\n" +
            "function renderHeatmap(){\n" +
            "  if(heatDone) return; heatDone=true;\n" +
            "  Plotly.newPlot('heat-plot',[{z:HEATMAP.matrix,x:HEATMAP.chroms,y:HEATMAP.lgs,\n" +
            "    type:'heatmap',colorscale:[[0,'#13162a'],[.01,'#1a2744'],[.3,'#2a5298'],[.6,'#7eb8f7'],[1,'#f7c948']],\n" +
            "    hovertemplate:'LG: %{y}<br>%{x}<br>Marcadores: %{z}<extra></extra>'}],{\n" +
            "    paper_bgcolor:'#0f1117',plot_bgcolor:'#13162a',\n" +
            "    font:{color:'#c0d0e0',size:11},\n" +
            "    xaxis:{title:'Cromosoma físico',tickangle:-45},\n" +
            "    yaxis:{title:'Grupo de Ligamiento',autorange:'reversed'},\n" +
            "    title:{text:'Distribución de marcadores: LG vs Cromosoma físico',font:{color:'#7eb8f7'}},\n" +
            "    margin:{t:50,r:30,b:100,l:80}\n" +
            "  },{responsive:true});\n" +
            "}\n\n" +

            "// ── Stats ────────────────────────────────────────────────────\n" +
            "let statsDone=false;\n" +
            "function renderStats(){\n" +
            "  if(statsDone) return; statsDone=true;\n" +
            "  const sorted=[...LG_STATS].sort((a,b)=>b.count-a.count);\n" +
            "  const top=sorted.slice(0,30);\n" +
            "  Plotly.newPlot('stats-plot',[\n" +
            "    {x:top.map(s=>s.lg),y:top.map(s=>s.count),type:'bar',name:'Marcadores',\n" +
            "     marker:{color:top.map(s=>chrColor(s.dominant))},\n" +
            "     hovertemplate:'%{x}: %{y} marcadores<extra></extra>'},\n" +
            "    {x:top.map(s=>s.lg),y:top.map(s=>s.length),type:'bar',name:'cM',\n" +
            "     marker:{color:top.map(s=>chrColor(s.dominant)),opacity:.5},yaxis:'y2',\n" +
            "     hovertemplate:'%{x}: %{y:.1f} cM<extra></extra>'}\n" +
            "  ],{\n" +
            "    paper_bgcolor:'#0f1117',plot_bgcolor:'#13162a',barmode:'group',\n" +
            "    font:{color:'#c0d0e0',size:11},\n" +
            "    xaxis:{title:'Grupo de Ligamiento',tickangle:-45},\n" +
            "    yaxis:{title:'N° Marcadores',gridcolor:'#1e2235'},\n" +
            "    yaxis2:{title:'Longitud (cM)',overlaying:'y',side:'right',gridcolor:'#1e2235'},\n" +
            "    legend:{bgcolor:'#1e2235',bordercolor:'#2a3f6f',borderwidth:1},\n" +
            "    title:{text:'Top 30 LGs — marcadores coloreados por cromosoma dominante',font:{color:'#7eb8f7'}},\n" +
            "    margin:{t:50,r:80,b:90,l:60}\n" +
            "  },{responsive:true});\n" +
            "  // Table\n" +
            "  let html='<table><thead><tr><th>LG</th><th>Chr dominante</th><th>% dom.</th>" +
            "<th>Composición cromosómica</th><th>N° Marcadores</th><th>Longitud (cM)</th><th>Densidad (mks/cM)</th></tr></thead><tbody>';\n" +
            "  sorted.forEach(s=>{\n" +
            "    const dens=s.length>0?(s.count/s.length).toFixed(2):'—';\n" +
            "    const bar='<div class=\"bar-inline\" style=\"width:'+(s.domPct*0.8).toFixed(0)+'px;background:'+chrColor(s.dominant)+'\"></div>';\n" +
            "    const chrComp=Object.entries(s.chrMap).sort((a,b)=>b[1]-a[1])\n" +
            "      .map(([c,n])=>'<span style=\"color:'+chrColor(c)+'\">chr'+c+'('+n+')</span>').join(' ');\n" +
            "    html+='<tr><td>'+s.lg+'</td><td>'+bar+'chr'+s.dominant+'</td><td>'+s.domPct+'%</td><td>'+chrComp+'</td><td>'+s.count+'</td><td>'+s.length.toFixed(1)+'</td><td>'+dens+'</td></tr>';\n" +
            "  });\n" +
            "  html+='</tbody></table>';\n" +
            "  document.getElementById('stats-table-wrap').innerHTML=html;\n" +
            "}\n\n" +

            "// ── Gene category colors ────────────────────────────────────\n" +
            "function geneColor(cat){\n" +
            "  if(cat&&cat.toLowerCase().includes('positiva')) return '#e74c3c';\n" +
            "  if(cat&&cat.toLowerCase().includes('neutral'))  return '#f7c948';\n" +
            "  if(cat&&cat.toLowerCase().includes('purificadora')) return '#58d68d';\n" +
            "  return '#c0c0c0';\n" +
            "}\n\n" +

            "// ── Overlay genes on LG bars ─────────────────────────────────\n" +
            "function overlayGenes(wrap, lg, maxCm, scale, showNames){\n" +
            "  if(!GENES||GENES.length===0) return;\n" +
            "  GENES.filter(g=>g.lg===lg&&g.cm>=0&&g.cm<=maxCm).forEach(g=>{\n" +
            "    const pct=maxCm>0?(g.cm/maxCm)*100:0;\n" +
            "    const color=geneColor(g.categoria);\n" +
            "    const pin=document.createElement('div');\n" +
            "    pin.style.cssText='position:absolute;left:-2px;right:-2px;height:4px;border-radius:2px;z-index:5;cursor:pointer;';\n" +
            "    pin.style.top=pct+'%';\n" +
            "    pin.style.background=color;\n" +
            "    pin.style.boxShadow='0 0 5px '+color;\n" +
            "    if(showNames){\n" +
            "      const lbl=document.createElement('div');\n" +
            "      lbl.style.cssText='position:absolute;left:22px;white-space:nowrap;font-size:9px;font-weight:600;pointer-events:none;z-index:6;transform:translateY(-50%);';\n" +
            "      lbl.style.top=pct+'%';\n" +
            "      lbl.style.color=color;\n" +
            "      lbl.textContent=g.id.replace(/\\.\\d+$/,'');\n" +
            "      wrap.parentElement.style.position='relative';\n" +
            "      wrap.parentElement.appendChild(lbl);\n" +
            "    }\n" +
            "    pin.addEventListener('mousemove',e=>{\n" +
            "      tip.innerHTML='<b style=\"color:'+color+'\">'+g.id+'</b><br>'\n" +
            "        +g.lg+' · '+g.cm.toFixed(2)+' cM<br>'\n" +
            "        +'<small>'+g.funcion+'</small><br>'\n" +
            "        +'<span style=\"color:'+color+'\">'+g.categoria+'</span>';\n" +
            "      tip.style.display='block';\n" +
            "      tip.style.left=(e.clientX+14)+'px'; tip.style.top=(e.clientY-10)+'px';\n" +
            "    });\n" +
            "    pin.addEventListener('mouseleave',hideTip);\n" +
            "    wrap.appendChild(pin);\n" +
            "  });\n" +
            "}\n\n" +

            "// Patch renderLgMap to call overlayGenes\n" +
            "const _origRenderLgMap=renderLgMap;\n" +
            "renderLgMap=function(){\n" +
            "  _origRenderLgMap();\n" +
            "  if(!GENES||GENES.length===0) return;\n" +
            "  const showNames=document.getElementById('showGeneNames')&&document.getElementById('showGeneNames').checked;\n" +
            "  document.querySelectorAll('.lg-col').forEach(col=>{\n" +
            "    const lbl=col.querySelector('.lg-label');\n" +
            "    if(!lbl) return;\n" +
            "    const lg=lbl.textContent.trim();\n" +
            "    const wrap=col.querySelector('.lg-bar-wrap');\n" +
            "    if(!wrap) return;\n" +
            "    const stat=LG_STATS.find(s=>s.lg===lg);\n" +
            "    if(!stat) return;\n" +
            "    overlayGenes(wrap,lg,stat.length,parseInt(document.getElementById('cmScale').value)||5,showNames);\n" +
            "  });\n" +
            "};\n\n" +

            "// ── Genes tab ────────────────────────────────────────────────\n" +
            "let genesDone=false;\n" +
            "function renderGenes(){\n" +
            "  if(genesDone) return; genesDone=true;\n" +
            "  if(!GENES||GENES.length===0){\n" +
            "    document.getElementById('genes-table-wrap').innerHTML='<p style=\"color:#556677;padding:20px\">No se cargaron genes candidatos. Usa --genes-map al generar el visor.</p>';\n" +
            "    return;\n" +
            "  }\n" +
            "  // Bar chart: genes per LG colored by category\n" +
            "  const cats=['Positiva','Neutral','Purificadora'];\n" +
            "  const colors={'Positiva':'#e74c3c','Neutral':'#f7c948','Purificadora':'#58d68d'};\n" +
            "  const lgSet=[...new Set(GENES.map(g=>g.lg))].sort((a,b)=>{\n" +
            "    const na=parseInt(a.replace(/\\D/g,'')),nb=parseInt(b.replace(/\\D/g,''));\n" +
            "    return na-nb;\n" +
            "  });\n" +
            "  const traces=cats.map(cat=>{\n" +
            "    return{x:lgSet,\n" +
            "      y:lgSet.map(lg=>GENES.filter(g=>g.lg===lg&&g.categoria&&g.categoria.includes(cat.toLowerCase())).length),\n" +
            "      name:cat,type:'bar',marker:{color:colors[cat]},\n" +
            "      hovertemplate:'%{x}: %{y} genes '+cat+'<extra></extra>'};\n" +
            "  });\n" +
            "  Plotly.newPlot('genes-plot',traces,{\n" +
            "    barmode:'stack',\n" +
            "    paper_bgcolor:'#0f1117',plot_bgcolor:'#13162a',\n" +
            "    font:{color:'#c0d0e0',size:12},\n" +
            "    xaxis:{title:'Grupo de Ligamiento',tickangle:-45},\n" +
            "    yaxis:{title:'N° Genes Candidatos',gridcolor:'#1e2235'},\n" +
            "    legend:{bgcolor:'#1e2235',bordercolor:'#2a3f6f',borderwidth:1},\n" +
            "    title:{text:'Distribución de genes de sacarosa por LG y presión selectiva',font:{color:'#7eb8f7'}},\n" +
            "    margin:{t:50,r:30,b:100,l:60}\n" +
            "  },{responsive:true});\n" +
            "  // Table\n" +
            "  const sorted=[...GENES].sort((a,b)=>{\n" +
            "    const order={'positiva':0,'neutral':1,'purificadora':2};\n" +
            "    const ca=(a.categoria||'').toLowerCase(); const cb=(b.categoria||'').toLowerCase();\n" +
            "    const oa=Object.keys(order).find(k=>ca.includes(k));\n" +
            "    const ob=Object.keys(order).find(k=>cb.includes(k));\n" +
            "    return (order[oa]||9)-(order[ob]||9)||a.lg.localeCompare(b.lg,undefined,{numeric:true});\n" +
            "  });\n" +
            "  let html='<table><thead><tr><th>Gen</th><th>LG</th><th>Pos (cM)</th><th>Presión Ka/Ks</th><th>Función</th></tr></thead><tbody>';\n" +
            "  sorted.forEach(g=>{\n" +
            "    const col=geneColor(g.categoria);\n" +
            "    const dot='<div style=\"display:inline-block;width:10px;height:10px;border-radius:50%;background:'+col+';margin-right:6px\"></div>';\n" +
            "    html+='<tr><td><b>'+g.id+'</b></td><td>'+g.lg+'</td><td>'+g.cm.toFixed(2)+'</td><td>'+dot+g.categoria+'</td><td><small>'+g.funcion+'</small></td></tr>';\n" +
            "  });\n" +
            "  html+='</tbody></table>';\n" +
            "  document.getElementById('genes-table-wrap').innerHTML=html;\n" +
            "}\n\n" +

            "// Patch showTab to handle genes-view\n" +
            "const _origShowTab=showTab;\n" +
            "showTab=function(id,btn){\n" +
            "  _origShowTab(id,btn);\n" +
            "  if(id==='genes-view') renderGenes();\n" +
            "};\n\n" +

            "// Init\n" +
            "buildLegend();\n" +
            "renderLgMap();\n";
    }
}
