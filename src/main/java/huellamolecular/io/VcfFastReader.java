package huellamolecular.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.stream.Stream;

/**
 * Lector VCF de Alto Rendimiento (Fast Reader).
 * Provee herramientas enfocadas en velocidad de procesamiento I/O y baja huella de recolección de basura (GC).
 */
public class VcfFastReader {

	/**
	 * Extrae los IDs de las muestras (Samples) directamente desde la cabecera del VCF.
	 * Asume que los identificadores de muestras vienen tras la columna FORMAT.
	 *
	 * @param vcfFilePath Ruta al archivo VCF.
	 * @return Array con los nombres de las muestras.
	 * @throws IOException Si el archivo no se puede abrir.
	 */
	public static String[] getSampleIds(String vcfFilePath) throws IOException {
		try (BufferedReader reader = Files.newBufferedReader(Paths.get(vcfFilePath))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("#CHROM")) {
					String[] columns = fastSplitByTab(line);
					if (columns.length > 9) {
						String[] samples = new String[columns.length - 9];
						System.arraycopy(columns, 9, samples, 0, samples.length);
						return samples;
					}
					return new String[0];
				}
			}
		}
		throw new RuntimeException("El archivo VCF no contiene una línea de cabecera válida (#CHROM...)");
	}

	/**
	 * Devuelve un Stream nativo de NIO para leer el VCF en frío.
	 * Automáticamente expulsa todas las líneas de cabecera (las que inician con '#').
	 * Ideal para operaciones funcionales: VcfFastReader.getDataLinesStream("archivo").forEach(...)
	 * 
	 * @param vcfFilePath Ruta absoluta o relativa al archivo VCF
	 * @return Stream de líneas (String) puramente de datos.
	 * @throws IOException Si el archivo no se puede abrir.
	 */
	public static Stream<String> getDataLinesStream(String vcfFilePath) throws IOException {
		Path path = Paths.get(vcfFilePath);
		// Aprovechamos Files.lines que carga de a bloques (Buffers optimizados por el SO)
		return Files.lines(path).filter(line -> !line.startsWith("#"));
	}

	/**
	 * Crea un Iterador que permite avanzar línea por línea y bloque a bloque,
	 * devolviendo ya la fila dividida (String[]) lista para ser procesada.
	 * Protege la RAM al no guardar toda la matriz en la memoria.
	 * 
	 * @param vcfFilePath Ruta al archivo VCF
	 * @return Una colección iterable sobre todos los SNPs
	 * @throws IOException Si hay error de apertura
	 */
	public static Iterable<String[]> iterateDataBlocks(String vcfFilePath) throws IOException {
		final BufferedReader reader = Files.newBufferedReader(Paths.get(vcfFilePath));

		return new Iterable<String[]>() {
			@Override
			public Iterator<String[]> iterator() {
				return new Iterator<String[]>() {
					private String nextLine = fetchNextDataLine();

					private String fetchNextDataLine() {
						try {
							String line;
							// Ignorar las cabeceras inmediatamente
							while ((line = reader.readLine()) != null) {
								if (!line.startsWith("#")) {
									return line;
								}
							}
							reader.close();
							return null;
						} catch (IOException e) {
							throw new RuntimeException("Error en lectura iterativa I/O", e);
						}
					}

					@Override
					public boolean hasNext() {
						return nextLine != null;
					}

					@Override
					public String[] next() {
						String[] columns = fastSplitByTab(nextLine);
						nextLine = fetchNextDataLine();
						return columns;
					}
				};
			}
		};
	}

	/**
	 * Separador optimizado. Evita la sobrecarga de Regex al aprovechar un solo char.
	 */
	public static String[] fastSplitByTab(String line) {
		// JVM fast-path: cuando el delimiter es 1 único caracter, NO se instancia Matcher Regex.
		return line.split("\t");
	}
}
