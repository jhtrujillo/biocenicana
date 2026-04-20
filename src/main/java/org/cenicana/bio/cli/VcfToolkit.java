package org.cenicana.bio.cli;
import org.cenicana.bio.utils.FileUtils;
import org.cenicana.bio.core.JoinMapCpFormat;
import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.VcfStatisticsCalculator;
import org.cenicana.bio.core.VcfFilter;

import java.io.IOException;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import org.cenicana.bio.cli.AlleleDosageCommand;
import org.cenicana.bio.cli.GeneticDistanceCommand;
import org.cenicana.bio.cli.LinkageDisequilibriumCommand;
import org.cenicana.bio.cli.VcfFilterCommand;
import org.cenicana.bio.cli.VcfStatsCommand;

@Command(name = "biocenicana", mixinStandardHelpOptions = true, version = "1.0", description = "Bioinformatics tools for Cenicana", subcommands = {
		AlleleDosageCommand.class,
		GeneticDistanceCommand.class,
		VcfStatsCommand.class,
		VcfFilterCommand.class,
		LinkageDisequilibriumCommand.class,
		GwasPolyExportCommand.class
})
public class VcfToolkit {

	public static void main(String[] args) {
		int exitCode = new CommandLine(new VcfToolkit()).execute(args);
		if (exitCode != 0) {
			System.exit(exitCode);
		}
	}













	@Command(name = "joinmap", aliases = { "20" }, description = "Join map format")
	public void joinmap(@Parameters(index = "0", description = "joinmap-file-ngsep") String joinmapFile)
			throws Exception {
		JoinMapCpFormat joinmap = new JoinMapCpFormat();
		joinmap.fixformat(joinmapFile);
	}






}
