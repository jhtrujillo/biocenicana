package org.cenicana.bio;

import java.io.IOException;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import org.cenicana.bio.cli.AlleleDosageCommand;
import org.cenicana.bio.cli.SimilitudGeneticaCommand;
import org.cenicana.bio.cli.GenDosisTargetedCommand;

@Command(name = "biocenicana", mixinStandardHelpOptions = true, version = "1.0", description = "Bioinformatics tools for Cenicana", subcommands = {
		AlleleDosageCommand.class,
		SimilitudGeneticaCommand.class,
		GenDosisTargetedCommand.class
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
		SeleccionarAbanicoDosis abanicodosis = new SeleccionarAbanicoDosis();
		abanicodosis.loadDosis(dosisFile);
	}

	@Command(name = "generarAlelosVCF", aliases = { "3" }, description = "Generar alelos VCF")
	public void generarAlelosVCF(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		GenerarAlelos alelos = new GenerarAlelos();
		alelos.getalelos(vcfFile, "dosis");
	}

	@Command(name = "seleccionarDosisAUS", aliases = { "4" }, description = "Filtros Alelológicos Australianos")
	public void seleccionarDosisAUS(@Parameters(index = "0", description = "snps_alelos.txt") String snpsAlelosTxt,
			@Parameters(index = "1", description = "Num Individuos") int numIndividuos) throws Exception {
		FiltrosAlelologicosAUS fa = new FiltrosAlelologicosAUS();
		fa.filtrarporClases(snpsAlelosTxt, numIndividuos);
	}

	@Command(name = "FiltrarVCF", aliases = { "5" }, description = "Filtrar VCF")
	public void FiltrarVCF(
			@Parameters(index = "0", description = "snps_dosis_aus.txt o snps_dosis_abanico.txt") String filterFile,
			@Parameters(index = "1", description = "path_vcf") String vcfFile) throws Exception {
		FiltrarVCF al = new FiltrarVCF();
		al.filtrarvcf(filterFile, vcfFile);
	}

	@Command(name = "ReducirHuellaVCF", aliases = { "6" }, description = "Reducir Huella VCF")
	public void ReducirHuellaVCF(@Parameters(index = "0", description = "path_vcf_original") String vcfOriginal,
			@Parameters(index = "1", description = "path_vcf_filtrado") String vcfFiltrado,
			@Parameters(index = "2", description = "numSNP") int numSnp) throws Exception {
		VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
		vcfmatrix.VCFfingerprint(vcfOriginal, vcfFiltrado, 0.0, 0.0, numSnp);
	}

	@Command(name = "frecuenciaAlelos", aliases = { "8" }, description = "Frecuencia Alelos")
	public void frecuenciaAlelos(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		VCFgetHaplotipes vcfcounter = new VCFgetHaplotipes();
		vcfcounter.CounterHaplotipes(vcfFile);
	}

	@Command(name = "CompararDosisHuellaVsTargeted", aliases = {
			"10" }, description = "Comparar Dosis Huella Vs Targeted")
	public void CompararDosisHuellaVsTargeted(
			@Parameters(index = "0", description = "dosisSecuenciacion") String dosisSec,
			@Parameters(index = "1", description = "dosisTargeted") String dosisTar) throws Exception {
		CompararDosisHuellaVsTargeted cdht = new CompararDosisHuellaVsTargeted();
		cdht.compararIndividuos(dosisSec, dosisTar);
	}

	@Command(name = "vcf-to-tab-targeted", aliases = { "11" }, description = "VCF to tab targeted")
	public void vcfToTabTargeted(@Parameters(index = "0", description = "vcf_targeted") String vcfTargeted)
			throws Exception {
		VcfToTabTargeted vcttotab = new VcfToTabTargeted();
		vcttotab.vcfToTab(vcfTargeted);
	}

	@Command(name = "vcftargetedTovcfNGSEP", aliases = { "12" }, description = "VCF targeted To VCF NGSEP")
	public void vcftargetedTovcfNGSEP(@Parameters(index = "0", description = "vcf_targeted") String vcfTargeted)
			throws Exception {
		SimilitudGeneticaCCdistTargeted smgt = new SimilitudGeneticaCCdistTargeted();
		smgt.formattoVCF(vcfTargeted);
	}

	@Command(name = "printDistanceMatrix", aliases = { "13" }, description = "Print Distance Matrix")
	public void printDistanceMatrix(@Parameters(index = "0", description = "path_vcf") String vcfFile)
			throws Exception {
		VCFgetfilterprint vcfmatrix = new VCFgetfilterprint();
		vcfmatrix.VCFload(vcfFile);
		System.out.print("p:10,GD:3 ");
		vcfmatrix.ImprimirMatrix();
	}

	@Command(name = "VCF-to-tab", aliases = { "14" }, description = "VCF to tab (geno)")
	public void vcfToTab(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		GenerarAlelos alelos = new GenerarAlelos();
		alelos.getalelos(vcfFile, "geno");
	}

	@Command(name = "generarDosisTranspuesta", aliases = { "15" }, description = "Generar Dosis Transpuesta")
	public void generarDosisTranspuesta(@Parameters(index = "0", description = "path_vcf") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "metodoImputar") String metodo) throws Exception {
		GeneDosis dosiscgene = new GeneDosis();
		dosiscgene.genDosisAlelicas(vcfFile, ploidy, metodo, true);
		dosiscgene.TransposeDosisMatrix();
		dosiscgene.printTransposeDosisMatrix();
	}

	@Command(name = "ordenarVCFxListadoIndividuos", aliases = {
			"16" }, description = "Ordenar VCF x Listado Individuos")
	public void ordenarVCFxListadoIndividuos(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ListadoIndividuos") String listado) throws Exception {
		OrdenarVCFporIndividuos sortVCF = new OrdenarVCFporIndividuos();
		sortVCF.ordenarVCFxListadoIndividuos(vcfFile, listado);
	}

	@Command(name = "vcfToStructure", aliases = { "17" }, description = "VCF To Structure")
	public void vcfToStructure(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "option") String option) throws Exception {
		VcfToStructure2 vcftosctructure = new VcfToStructure2();
		vcftosctructure.vcfconverTostructure(vcfFile, ploidy, option);
		vcftosctructure.printMatrix();
	}

	@Command(name = "vcfToStructureAllele", aliases = { "18" }, description = "VCF To Structure Allele")
	public void vcfToStructureAllele(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "option") String option,
			@Parameters(index = "3", description = "impute") String impute) throws Exception {
		VcfToStructure2 vcftosctructure = new VcfToStructure2();
		vcftosctructure.vcfconverTostructureAlleles(vcfFile, ploidy, option, impute);
		vcftosctructure.printMatrixTranspuesta();
	}

	@Command(name = "addfuntionstogff", aliases = { "19" }, description = "Add functions to gff")
	public void addfuntionstogff(@Parameters(index = "0", description = "gffFile") String gffFile) throws Exception {
		AddFunctionsGff gff = new AddFunctionsGff();
		gff.loadgff(gffFile);
	}

	@Command(name = "joinmap", aliases = { "20" }, description = "Join map format")
	public void joinmap(@Parameters(index = "0", description = "joinmap-file-ngsep") String joinmapFile)
			throws Exception {
		JoinMapCpFormat joinmap = new JoinMapCpFormat();
		joinmap.fixformat(joinmapFile);
	}

	@Command(name = "protocolnotes", aliases = { "21" }, description = "Protocol notes")
	public void protocolnotes(
			@Parameters(index = "0", description = "opcionPN (kasp-gbsrad|kasp-tr)") String opcionPN,
			@Parameters(index = "1", description = "file_dosis_kasp") String fileDosisKasp,
			@Parameters(index = "2", description = "file_dosis_gbs_rad_OR_tr") String fileDosisGbsRad)
			throws Exception {
		if (opcionPN.compareTo("kasp-gbsrad") == 0) {
			Protocol_notes pn = new Protocol_notes(fileDosisKasp, fileDosisGbsRad, "", "");
			pn.get_all_dosis_kasp_gbs_rad();
		} else if (opcionPN.compareTo("kasp-tr") == 0) {
			Protocol_notes pn = new Protocol_notes(fileDosisKasp, "", "", fileDosisGbsRad);
			pn.get_all_dosis_kasp_tr();
		}
	}

	@Command(name = "get_depth_seq", aliases = { "22" }, description = "Get depth seq")
	public void get_depth_seq(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		GeneDosis dosiscgene = new GeneDosis();
		dosiscgene.get_depth_sequ(vcfFile);
		dosiscgene.printDosisMatrix();
	}

	@Command(name = "get_dosis_freebayes", aliases = { "23" }, description = "Get dosis freebayes")
	public void get_dosis_freebayes(@Parameters(index = "0", description = "path_vcf") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy) throws Exception {
		GeneDosisRapidGenomicVCF dosiscgene = new GeneDosisRapidGenomicVCF();
		dosiscgene.genDosisAlelicas(vcfFile, ploidy);
	}

	@Command(name = "rapidgenomic", aliases = { "24" }, description = "Rapid genomic")
	public void rapidgenomic(@Parameters(index = "0", description = "fastafile") String fastafile,
			@Parameters(index = "1", description = "patron") String patron) throws Exception {
		RapidGenomic rg = new RapidGenomic();
		rg.leerarchivo2(fastafile, patron);
	}

	@Command(name = "generarSNPsformatoRG", aliases = { "25" }, description = "Generar SNPs formato RG")
	public void generarSNPsformatoRG(@Parameters(index = "0", description = "posicionesSNPsgenoCompleto") String pos,
			@Parameters(index = "1", description = "VCF") String vcf) throws Exception {
		RapidGenomic rg = new RapidGenomic();
		rg.generarSNPsformatoRG(pos, vcf);
	}

	@Command(name = "generarSNPsformatoRG2", aliases = { "26" }, description = "Generar SNPs formato RG2")
	public void generarSNPsformatoRG2(@Parameters(index = "0", description = "posiciones") String pos,
			@Parameters(index = "1", description = "VCF") String vcf,
			@Parameters(index = "2", description = "Chr") String chr,
			@Parameters(index = "3", description = "SNP") String snp) throws Exception {
		RapidGenomic rg = new RapidGenomic();
		rg.generarSNPsformatoRG(pos, vcf);
	}

	@Command(name = "generarSNPsformatoRG3", aliases = { "27" }, description = "Generar SNPs formato RG3")
	public void generarSNPsformatoRG3(@Parameters(index = "0", description = "posiciones") String pos,
			@Parameters(index = "1", description = "VCF") String vcf) throws Exception {
		RapidGenomic rg = new RapidGenomic();
		rg.generarSNPsformatoRG3(pos, vcf);
	}

	@Command(name = "fixgffformat", aliases = { "28" }, description = "Fix gff format")
	public void fixgffformat(@Parameters(index = "0", description = "gffFile") String gffFile) throws Exception {
		AddFunctionsGff gff = new AddFunctionsGff();
		gff.fixgffformat(gffFile);
	}

	@Command(name = "filtrargffporTamano", aliases = { "29" }, description = "Filtrar gff por tamaño")
	public void filtrargffporTamano(@Parameters(index = "0", description = "gffFile") String gffFile,
			@Parameters(index = "1", description = "minSize") String minSize,
			@Parameters(index = "2", description = "maxSize") String maxSize) throws Exception {
		AddFunctionsGff gff = new AddFunctionsGff();
		gff.filtrargffporTamano(gffFile, minSize, maxSize);
	}

	@Command(name = "ajustargfffunciones", description = "Ajustar gff funciones")
	public void ajustargfffunciones(@Parameters(index = "0", description = "gffFile") String gffFile) throws Exception {
		AddFunctionsGff gff = new AddFunctionsGff();
		gff.ajustargfffunciones(gffFile);
	}

	@Command(name = "vcfToACGT", aliases = { "30" }, description = "VCF To ACGT")
	public void vcfToACGT(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "option") String option,
			@Parameters(index = "3", description = "impute") String impute) throws Exception {
		VcfToStructure2 vcftosctructure = new VcfToStructure2();
		vcftosctructure.vcfconverTostructureAlleles(vcfFile, ploidy, option, impute);
		vcftosctructure.printMatrixTranspuesta();
	}
}
