package org.cenicana.bio.cli;

import org.cenicana.bio.TargetedAlleleDosage;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import java.util.concurrent.Callable;

@Command(name = "genDosisTargeted", aliases = {"9"}, mixinStandardHelpOptions = true,
		description = "Genera dosis targeted a partir de un VCF.")
public class TargetedAlleleDosageCommand implements Callable<Integer> {

	@Parameters(index = "0", description = "Ruta al archivo VCF")
	private String vcfFile;

	@Parameters(index = "1", description = "Ploidía (ej. 2, 4, 10)")
	private int ploidy;

	@Override
	public Integer call() throws Exception {
		try {
			TargetedAlleleDosage genDosis = new TargetedAlleleDosage();
			genDosis.generarDosis(vcfFile, ploidy);
			return 0;
		} catch (Exception e) {
			System.err.println("Error ejecutando genDosisTargeted: " + e.getMessage());
			e.printStackTrace();
			return 1;
		}
	}
}
