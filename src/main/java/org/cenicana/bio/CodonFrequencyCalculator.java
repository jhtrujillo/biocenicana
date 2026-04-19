package org.cenicana.bio;

import java.io.*;
import java.util.*;

public class CodonFrequencyCalculator {

	// Mapa de codones a aminoácidos
	private static final Map<String, String> CODON_TABLE = new HashMap<>();
	static {
		CODON_TABLE.put("TTT", "Phe");
		CODON_TABLE.put("TTC", "Phe");
		CODON_TABLE.put("TTA", "Leu");
		CODON_TABLE.put("TTG", "Leu");
		CODON_TABLE.put("CTT", "Leu");
		CODON_TABLE.put("CTC", "Leu");
		CODON_TABLE.put("CTA", "Leu");
		CODON_TABLE.put("CTG", "Leu");
		CODON_TABLE.put("ATT", "Ile");
		CODON_TABLE.put("ATC", "Ile");
		CODON_TABLE.put("ATA", "Met");
		CODON_TABLE.put("ATG", "Met");
		CODON_TABLE.put("GTT", "Val");
		CODON_TABLE.put("GTC", "Val");
		CODON_TABLE.put("GTA", "Val");
		CODON_TABLE.put("GTG", "Val");
		CODON_TABLE.put("TCT", "Ser");
		CODON_TABLE.put("TCC", "Ser");
		CODON_TABLE.put("TCA", "Ser");
		CODON_TABLE.put("TCG", "Ser");
		CODON_TABLE.put("CCT", "Pro");
		CODON_TABLE.put("CCC", "Pro");
		CODON_TABLE.put("CCA", "Pro");
		CODON_TABLE.put("CCG", "Pro");
		CODON_TABLE.put("ACT", "Thr");
		CODON_TABLE.put("ACC", "Thr");
		CODON_TABLE.put("ACA", "Thr");
		CODON_TABLE.put("ACG", "Thr");
		CODON_TABLE.put("GCT", "Ala");
		CODON_TABLE.put("GCC", "Ala");
		CODON_TABLE.put("GCA", "Ala");
		CODON_TABLE.put("GCG", "Ala");
		CODON_TABLE.put("TAT", "Tyr");
		CODON_TABLE.put("TAC", "Tyr");
		CODON_TABLE.put("TAA", "STOP");
		CODON_TABLE.put("TAG", "STOP");
		CODON_TABLE.put("CAT", "His");
		CODON_TABLE.put("CAC", "His");
		CODON_TABLE.put("CAA", "Gln");
		CODON_TABLE.put("CAG", "Gln");
		CODON_TABLE.put("AAT", "Asn");
		CODON_TABLE.put("AAC", "Asn");
		CODON_TABLE.put("AAA", "Lys");
		CODON_TABLE.put("AAG", "Lys");
		CODON_TABLE.put("GAT", "Asp");
		CODON_TABLE.put("GAC", "Asp");
		CODON_TABLE.put("GAA", "Glu");
		CODON_TABLE.put("GAG", "Glu");
		CODON_TABLE.put("TGT", "Cys");
		CODON_TABLE.put("TGC", "Cys");
		CODON_TABLE.put("TGA", "Trp");
		CODON_TABLE.put("TGG", "Trp");
		CODON_TABLE.put("CGT", "Arg");
		CODON_TABLE.put("CGC", "Arg");
		CODON_TABLE.put("CGA", "Arg");
		CODON_TABLE.put("CGG", "Arg");
		CODON_TABLE.put("AGT", "Ser");
		CODON_TABLE.put("AGC", "Ser");
		CODON_TABLE.put("AGA", "Arg");
		CODON_TABLE.put("AGG", "Arg");
		CODON_TABLE.put("GGT", "Gly");
		CODON_TABLE.put("GGC", "Gly");
		CODON_TABLE.put("GGA", "Gly");
		CODON_TABLE.put("GGG", "Gly");
	}

	// Leer el archivo FASTA y obtener la secuencia
	private static String readFasta(String fastaFile) throws IOException {
		StringBuilder sequence = new StringBuilder();
		BufferedReader br = new BufferedReader(new FileReader(fastaFile));
		String line;
		while ((line = br.readLine()) != null) {
			if (!line.startsWith(">")) {
				sequence.append(line.trim());
			}
		}
		br.close();
		return sequence.toString();
	}

	
	  // Contar el número de secuencias en el archivo FASTA (líneas que comienzan con '>')
    private static int countSequencesInFasta(String fastaFile) throws IOException {
        int sequenceCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(fastaFile));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith(">")) {
                sequenceCount++;
            }
        }
        br.close();
        return sequenceCount;
    }

	 
	
	 // Calcular la frecuencia de los codones
    private static Map<String, Integer> calculateCodonCounts(String sequence) {
        Map<String, Integer> codonCounts = new HashMap<>();
        for (int i = 0; i < sequence.length() - 2; i += 3) {
            String codon = sequence.substring(i, i + 3);
            if (CODON_TABLE.containsKey(codon)) {
                codonCounts.put(codon, codonCounts.getOrDefault(codon, 0) + 1);
            }
        }
        return codonCounts;
    }

    // Calcular la frecuencia de los aminoácidos a partir de los codones
    private static Map<String, Integer> calculateAminoAcidCounts(Map<String, Integer> codonCounts) {
        Map<String, Integer> aminoAcidCounts = new HashMap<>();
        for (Map.Entry<String, Integer> entry : codonCounts.entrySet()) {
            String codon = entry.getKey();
            String aminoAcid = CODON_TABLE.get(codon);
            //if (!aminoAcid.equals("STOP")) {
                aminoAcidCounts.put(aminoAcid, aminoAcidCounts.getOrDefault(aminoAcid, 0) + entry.getValue());
            //}
        }
        return aminoAcidCounts;
    }
    

    // Imprimir las frecuencias de codones y aminoácidos
    private static void printCodonFrequencies(Map<String, Integer> codonCounts, Map<String, Integer> aminoAcidCounts) {
        System.out.println("Codón\tAminoácido\tConteo_Codón\tConteo_Aminoácido\tFrecuencia_Relativa\tFrecuencia_porMIles");
        for (Map.Entry<String, String> entry : CODON_TABLE.entrySet()) {
            String codon = entry.getKey();
            String aminoAcid = entry.getValue();
            int codonCount = codonCounts.getOrDefault(codon, 0);
            int aminoAcidCount = aminoAcidCounts.getOrDefault(aminoAcid, 0);
            double frequency = aminoAcidCount > 0 ? (double) codonCount / aminoAcidCount : 0;
            double frequencyPerThousand = aminoAcidCount > 0 ? (double) codonCount / aminoAcidCount * 1000 : 0;
            
            //if (!aminoAcid.equals("STOP")) {
            	 System.out.printf("%s\t%s\t%d\t%d\t%.4f\t%.4f%n", codon, aminoAcid, codonCount, aminoAcidCount, frequency, frequencyPerThousand);
            //}
           

        }
    }
    
    public static void main(String[] args) throws Exception {
//        if (args.length < 1) {
//            System.out.println("Por favor, proporciona un archivo FASTA como argumento.");
//            System.exit(1);
//        }

        String fastaFile = "//Users/estuvar4/Downloads/CC-01-1940.cds.fna";
        //String fastaFile = "//Users/estuvar4/Downloads/Saccharum_hybrid_cultivar_R570.cds.fna";
        //String fastaFile = "//Users/estuvar4/Downloads/Saccharum_officinarum_LA-Purple.cds.fna";
        
    	//String fastaFile = "//Users/estuvar4/Downloads/Erianthus_rufipilus.cds.fna";
    	
      
   
        
        
        String sequence = readFasta(fastaFile);
        
        

        // Calcula las frecuencias de codones y aminoácidos
        Map<String, Integer> codonCounts = calculateCodonCounts(sequence);
        Map<String, Integer> aminoAcidCounts = calculateAminoAcidCounts(codonCounts);
        
        
        // Contar el número de secuencias en el archivo FASTA
        int numberOfSequences = countSequencesInFasta(fastaFile);
        System.out.println("Número de secuencias en el archivo FASTA: " + numberOfSequences);

        
        
        
        // Imprimir la cantidad total de codones
        int totalCodons = 0;
        for (int count : codonCounts.values()) {
            totalCodons += count;
        }
        System.out.println("Número total de codones utilizados: " + totalCodons);
        
        
        // Imprimir los resultados
        printCodonFrequencies(codonCounts, aminoAcidCounts);

//     // Imprimir las entradas del mapa
//        
//        for (Map.Entry<String, Integer> entry : codonCounts.entrySet()) {
//            String codon = entry.getKey();
//            Integer count = entry.getValue();
//            System.out.println("Codón " + codon + ": Frecuencia: " + count);
//        }
//        
//        
//        
//        for (Map.Entry<String, Integer> entry : aminoAcidCounts.entrySet()) {
//            String codon = entry.getKey();
//            Integer count = entry.getValue();
//            System.out.println("AminoAcid " + codon + ": Frecuencia: " + count);
//        }
        
       
        
        
        
    }
    
    
    

}
