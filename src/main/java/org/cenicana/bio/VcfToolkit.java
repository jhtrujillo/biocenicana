package org.cenicana.bio;

import java.io.IOException;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import org.cenicana.bio.cli.AlleleDosageCommand;
import org.cenicana.bio.cli.GeneticDistanceCommand;
import org.cenicana.bio.cli.TargetedAlleleDosageCommand;

@Command(name = "biocenicana", mixinStandardHelpOptions = true, version = "1.0", description = "Bioinformatics tools for Cenicana", subcommands = {
		AlleleDosageCommand.class,
		TargetedAlleleDosageCommand.class,
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

	@Command(name = "seleccionarDosisAUS", aliases = { "4" }, description = "Filtros Alelológicos Australianos")
	public void seleccionarDosisAUS(@Parameters(index = "0", description = "snps_alelos.txt") String snpsAlelosTxt,
			@Parameters(index = "1", description = "Num Individuos") int numIndividuos) throws Exception {
		AllelicFiltersAus fa = new AllelicFiltersAus();
		fa.filtrarporClases(snpsAlelosTxt, numIndividuos);
	}

	@Command(name = "VcfFilter", aliases = { "5" }, description = "Filtrar VCF")
	public void VcfFilter(
			@Parameters(index = "0", description = "snps_dosis_aus.txt o snps_dosis_abanico.txt") String filterFile,
			@Parameters(index = "1", description = "path_vcf") String vcfFile) throws Exception {
		VcfFilter al = new VcfFilter();
		al.filtrarvcf(filterFile, vcfFile);
	}

	@Command(name = "ReducirHuellaVCF", aliases = { "6" }, description = "Reducir Huella VCF")
	public void ReducirHuellaVCF(@Parameters(index = "0", description = "path_vcf_original") String vcfOriginal,
			@Parameters(index = "1", description = "path_vcf_filtrado") String vcfFiltrado,
			@Parameters(index = "2", description = "numSNP") int numSnp) throws Exception {
		VcfFilterPrint vcfmatrix = new VcfFilterPrint();
		vcfmatrix.VCFfingerprint(vcfOriginal, vcfFiltrado, 0.0, 0.0, numSnp);
	}

	@Command(name = "frecuenciaAlelos", aliases = { "8" }, description = "Frecuencia Alelos")
	public void frecuenciaAlelos(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		VcfHaplotypeExtractor vcfcounter = new VcfHaplotypeExtractor();
		vcfcounter.CounterHaplotipes(vcfFile);
	}

	@Command(name = "CompareDosageFingerprintVsTargeted", aliases = {
			"10" }, description = "Comparar Dosis Huella Vs Targeted")
	public void CompareDosageFingerprintVsTargeted(
			@Parameters(index = "0", description = "dosisSecuenciacion") String dosisSec,
			@Parameters(index = "1", description = "dosisTargeted") String dosisTar) throws Exception {
		CompareDosageFingerprintVsTargeted cdht = new CompareDosageFingerprintVsTargeted();
		cdht.compararIndividuos(dosisSec, dosisTar);
	}

	@Command(name = "vcf-to-tab-targeted", aliases = { "11" }, description = "VCF to tab targeted")
	public void vcfToTabTargeted(@Parameters(index = "0", description = "vcf_targeted") String vcfTargeted)
			throws Exception {
		TargetedVcfToTab vcttotab = new TargetedVcfToTab();
		vcttotab.vcfToTab(vcfTargeted);
	}

	@Command(name = "vcftargetedTovcfNGSEP", aliases = { "12" }, description = "VCF targeted To VCF NGSEP")
	public void vcftargetedTovcfNGSEP(@Parameters(index = "0", description = "vcf_targeted") String vcfTargeted)
			throws Exception {
		TargetedGeneticSimilarity smgt = new TargetedGeneticSimilarity();
		smgt.formattoVCF(vcfTargeted);
	}

	@Command(name = "printDistanceMatrix", aliases = { "13" }, description = "Print Distance Matrix")
	public void printDistanceMatrix(@Parameters(index = "0", description = "path_vcf") String vcfFile)
			throws Exception {
		VcfFilterPrint vcfmatrix = new VcfFilterPrint();
		vcfmatrix.VCFload(vcfFile);
		System.out.print("p:10,GD:3 ");
		vcfmatrix.ImprimirMatrix();
	}

	@Command(name = "VCF-to-tab", aliases = { "14" }, description = "VCF to tab (geno)")
	public void vcfToTab(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		AlleleGenerator alelos = new AlleleGenerator();
		alelos.getalelos(vcfFile, "geno");
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

	@Command(name = "ordenarVCFxListadoIndividuos", aliases = {
			"16" }, description = "Ordenar VCF x Listado Individuos")
	public void ordenarVCFxListadoIndividuos(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ListadoIndividuos") String listado) throws Exception {
		SortVcfBySample sortVCF = new SortVcfBySample();
		sortVCF.ordenarVCFxListadoIndividuos(vcfFile, listado);
	}

	@Command(name = "vcfToStructure", aliases = { "17" }, description = "VCF To Structure")
	public void vcfToStructure(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "option") String option) throws Exception {
		VcfToStructure vcftosctructure = new VcfToStructure();
		vcftosctructure.vcfconverTostructure(vcfFile, ploidy, option);
		vcftosctructure.printMatrix();
	}

	@Command(name = "vcfToStructureAllele", aliases = { "18" }, description = "VCF To Structure Allele")
	public void vcfToStructureAllele(@Parameters(index = "0", description = "VCFfile") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy,
			@Parameters(index = "2", description = "option") String option,
			@Parameters(index = "3", description = "impute") String impute) throws Exception {
		VcfToStructure vcftosctructure = new VcfToStructure();
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
			ProtocolNotes pn = new ProtocolNotes(fileDosisKasp, fileDosisGbsRad, "", "");
			pn.get_all_dosis_kasp_gbs_rad();
		} else if (opcionPN.compareTo("kasp-tr") == 0) {
			ProtocolNotes pn = new ProtocolNotes(fileDosisKasp, "", "", fileDosisGbsRad);
			pn.get_all_dosis_kasp_tr();
		}
	}

	@Command(name = "get_depth_seq", aliases = { "22" }, description = "Get depth seq")
	public void get_depth_seq(@Parameters(index = "0", description = "path_vcf") String vcfFile) throws Exception {
		AlleleDosageCalculator dosiscgene = new AlleleDosageCalculator();
		dosiscgene.get_depth_sequ(vcfFile);
		dosiscgene.printDosisMatrix();
	}

	@Command(name = "get_dosis_freebayes", aliases = { "23" }, description = "Get dosis freebayes")
	public void get_dosis_freebayes(@Parameters(index = "0", description = "path_vcf") String vcfFile,
			@Parameters(index = "1", description = "ploidy") int ploidy) throws Exception {
		RapidGenomicAlleleDosage dosiscgene = new RapidGenomicAlleleDosage();
		dosiscgene.computeAlleleDosage(vcfFile, ploidy);
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
		VcfToStructure vcftosctructure = new VcfToStructure();
		vcftosctructure.vcfconverTostructureAlleles(vcfFile, ploidy, option, impute);
		vcftosctructure.printMatrixTranspuesta();
	}
}
