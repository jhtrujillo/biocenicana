package huellamolecular.cli;

import huellamolecular.VCFgetfilterprint;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import java.util.concurrent.Callable;

@Command(name = "similitudGeneitcaCCdist", aliases = {"7"}, mixinStandardHelpOptions = true,
		description = "Calcula las estadísticas de similitud genética.")
public class SimilitudGeneticaCommand implements Callable<Integer> {

	@Parameters(index = "0", description = "Ruta al archivo VCF")
	private String vcfFile;

	@Override
	public Integer call() throws Exception {
		try {
			VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
			vcfmatrix.VCFload(vcfFile);
			System.out.print("p:10,GD:3 ");
			vcfmatrix.getSimilitudeStats(System.out);
			return 0;
		} catch (Exception e) {
			System.err.println("Error ejecutando similitudGeneitcaCCdist: " + e.getMessage());
			e.printStackTrace();
			return 1;
		}
	}
}
