package org.cenicana.bio;

import java.io.IOException;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import org.cenicana.bio.cli.AlleleDosageCommand;
import org.cenicana.bio.cli.GeneticDistanceCommand;

@Command(name = "biocenicana", mixinStandardHelpOptions = true, version = "1.0", description = "Bioinformatics tools for Cenicana", subcommands = {
		AlleleDosageCommand.class,
		GeneticDistanceCommand.class
})
public class VcfToolkit {

	public static void main(String[] args) {
		int exitCode = new CommandLine(new VcfToolkit()).execute(args);
		if (exitCode != 0) {
			System.exit(exitCode);
		}
	}

	@Command(name = "seleccionarDosisAbanico", aliases = { "2" }, description = "Selecciona dosis abanico")
	public void seleccionarDosisAbanico(@Parameters(index = "0", description = "snps_dosis.txt") String dosisFile)
			throws Exception {
		SelectDosageRange abanicodosis = new SelectDosageRange();
		abanicodosis.loadDosis(dosisFile);
	}

	@Command(name = "generarAlelosVCF", aliases = { "3" }, description = "Generar alelos VCF")
	public void generarAlelosVCF(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		AlleleGenerator alelos = new AlleleGenerator();
		alelos.getalelos(vcfFile, "dosis");
	}



	@Command(name = "frecuenciaAlelos", aliases = { "8" }, description = "Frecuencia Alelos")
	public void frecuenciaAlelos(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		VcfHaplotypeExtractor vcfcounter = new VcfHaplotypeExtractor();
		vcfcounter.CounterHaplotipes(vcfFile);
	}





	@Command(name = "computeDosageTransposed", aliases = { "15" }, description = "Generar Dosis Transpuesta")
	public void computeDosageTransposed(@Parameters(index = "0", description = "path_vcf") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "imputationMethod") String metodo) throws Exception {
		AlleleDosageCalculator dosiscgene = new AlleleDosageCalculator();
		dosiscgene.computeAlleleDosage(vcfFile, ploidy, metodo, true);
		dosiscgene.TransposeDosisMatrix();
		dosiscgene.printTransposeDosisMatrix();
	}




	@Command(name = "joinmap", aliases = { "20" }, description = "Join map format")
	public void joinmap(@Parameters(index = "0", description = "joinmap-file-ngsep") String joinmapFile)
			throws Exception {
		JoinMapCpFormat joinmap = new JoinMapCpFormat();
		joinmap.fixformat(joinmapFile);
	}






}
