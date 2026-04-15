package huellamolecular;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Representa la cabecera de un archivo VCF.
 * Contiene las líneas de metadatos (##), los nombres de las muestras
 * y la definición de los campos FORMAT e INFO.
 */
public class VcfHeader {

    /** Líneas de metadatos (##INFO, ##FORMAT, ##contig, etc.) */
    private final List<String> metaLines = new ArrayList<>();

    /** Nombres de las muestras en el orden del archivo (columnas 9+) */
    private final List<String> sampleIds = new ArrayList<>();

    /** Definición de cada campo FORMAT: clave → descripción */
    private final Map<String, String> formatFields = new LinkedHashMap<>();

    /** Línea completa de la cabecera de columnas (#CHROM POS ID ...) */
    private String headerLine = "";

    // -------------------------------------------------------------------------
    // Métodos de construcción (usados por VcfReader)
    // -------------------------------------------------------------------------

    /**
     * Agrega una línea de metadatos (##).
     * Si la línea es ##FORMAT, extrae la clave y la descripción.
     */
    public void addMetaLine(String line) {
        metaLines.add(line);
        if (line.startsWith("##FORMAT=")) {
            parseFormatLine(line);
        }
    }

    /**
     * Procesa la línea #CHROM y extrae los IDs de las muestras (columnas >= 9).
     */
    public void setHeaderLine(String line) {
        this.headerLine = line;
        String[] cols = line.split("\t");
        sampleIds.clear();
        for (int i = 9; i < cols.length; i++) {
            sampleIds.add(cols[i].trim());
        }
    }

    // -------------------------------------------------------------------------
    // Accesores
    // -------------------------------------------------------------------------

    /** Retorna los nombres de las muestras en el orden del archivo. */
    public List<String> getSampleIds() {
        return sampleIds;
    }

    /** Retorna todas las líneas de metadatos (##). */
    public List<String> getMetaLines() {
        return metaLines;
    }

    /** Retorna el mapa de campos FORMAT definidos en la cabecera. */
    public Map<String, String> getFormatFields() {
        return formatFields;
    }

    /** Retorna la línea de cabecera de columnas (#CHROM POS …). */
    public String getHeaderLine() {
        return headerLine;
    }

    /** Retorna el número de muestras en el archivo. */
    public int getSampleCount() {
        return sampleIds.size();
    }

    /**
     * Imprime la cabecera completa (metadatos + línea de columnas) al stream dado.
     *
     * @param out PrintStream de destino (ej: System.out o un FileOutputStream).
     */
    public void print(PrintStream out) {
        for (String meta : metaLines) {
            out.println(meta);
        }
        if (!headerLine.isEmpty()) {
            out.println(headerLine);
        }
    }

    // -------------------------------------------------------------------------
    // Métodos internos
    // -------------------------------------------------------------------------

    /**
     * Extrae la clave ID y la descripción de una línea ##FORMAT.
     * Ejemplo: ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
     */
    private void parseFormatLine(String line) {
        try {
            // Extrae contenido entre < y >
            String content = line.substring(line.indexOf('<') + 1, line.lastIndexOf('>'));
            String id = null;
            String description = "";
            for (String part : content.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)")) {
                if (part.startsWith("ID=")) {
                    id = part.substring(3).trim();
                } else if (part.startsWith("Description=")) {
                    description = part.substring(12).replace("\"", "").trim();
                }
            }
            if (id != null) {
                formatFields.put(id, description);
            }
        } catch (Exception e) {
            // Línea FORMAT malformada — se ignora silenciosamente
        }
    }
}
