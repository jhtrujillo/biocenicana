package org.cenicana.bio.io;
import org.cenicana.bio.utils.FileUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.stream.Stream;

/**
 * High-Performance VCF Reader (Fast Reader).
 * Provides tools focused on I/O speed and low garbage collection (GC) footprint.
 */
public class VcfFastReader {

	/**
	 * Extracts sample IDs directly from the VCF header.
	 * Assumes sample identifiers follow the FORMAT column.
	 *
	 * @param vcfFilePath Path to the VCF file.
	 * @return Array of sample names.
	 * @throws IOException If the file cannot be opened.
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
		throw new RuntimeException("The VCF file does not contain a valid header line (#CHROM...)");
	}

	/**
	 * Returns a native NIO Stream to read the VCF lazily.
	 * Automatically skips all header lines (starting with '#').
	 * Ideal for functional operations: VcfFastReader.getDataLinesStream("file").forEach(...)
	 * 
	 * @param vcfFilePath Absolute or relative path to the VCF file.
	 * @return Stream of data lines.
	 * @throws IOException If the file cannot be opened.
	 */
	public static Stream<String> getDataLinesStream(String vcfFilePath) throws IOException {
		Path path = Paths.get(vcfFilePath);
		// Leverage Files.lines which loads by blocks (buffered by OS)
		return Files.lines(path).filter(line -> !line.startsWith("#"));
	}

	/**
	 * Creates an Iterator that advances line by line and block by block,
	 * returning the split row (String[]) ready for processing.
	 * Protects RAM by not keeping the entire matrix in memory.
	 * 
	 * @param vcfFilePath Path to the VCF file.
	 * @return An iterable collection over all SNPs.
	 * @throws IOException If there is an I/O error.
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
							// Skip headers immediately
							while ((line = reader.readLine()) != null) {
								if (!line.startsWith("#")) {
									return line;
								}
							}
							reader.close();
							return null;
						} catch (IOException e) {
							try { reader.close(); } catch (IOException ignored) {}
							throw new RuntimeException("Error during iterative I/O reading", e);
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
	 * Optimized separator. Avoids regex overhead by leveraging a single character.
	 */
	public static String[] fastSplitByTab(String line) {
		// JVM fast-path: when the delimiter is a single character, Matcher Regex is NOT instantiated.
		return line.split("\t");
	}
}
