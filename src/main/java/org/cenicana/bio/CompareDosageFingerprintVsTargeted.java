package org.cenicana.bio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class CompareDosageFingerprintVsTargeted {

	public void compararIndividuos(String dosisSecuenciacion, String dosisTargeted) {
		try {
			// Preload Sequenced data into memory Map explicitly instead of scanning multiple times
			Map<String, String[]> mapSec = new HashMap<>();
			Map<String, Integer> headerSec = new HashMap<>();
			
			BufferedReader brSec = new BufferedReader(new FileReader(dosisSecuenciacion));
			String lineSec = brSec.readLine();
			if (lineSec != null) {
				String[] colsSec = lineSec.split("\t");
				for(int i = 0; i < colsSec.length; i++) {
				    headerSec.put(colsSec[i], i);
				}
			}
			
			while((lineSec = brSec.readLine()) != null) {
				String[] parts = lineSec.split("\t");
				if (parts.length > 2) {
				    mapSec.put(parts[0] + ":" + parts[1], parts);
				}
			}
			brSec.close();
			
			// Stream targeted data
			BufferedReader brTar = new BufferedReader(new FileReader(dosisTargeted));
			String lineTar = brTar.readLine();
			if (lineTar == null) {
			    brTar.close();
			    return;
			}
			
			String[] colsTar = lineTar.replace("CONTROL_", "").replace("PV_", "").replace("_R", "").split("\t");
			Map<String, Integer> headerTar = new HashMap<>();
			for(int i = 0; i < colsTar.length; i++) {
			    headerTar.put(colsTar[i], i);
			}
			
			while((lineTar = brTar.readLine()) != null) {
				String[] partsTar = lineTar.split("\t");
				if (partsTar.length < 2) continue;
				String key = partsTar[0] + ":" + partsTar[1];
				
				if (mapSec.containsKey(key)) {
				    String[] partsSec = mapSec.get(key);
					for (int j = 4; j < colsTar.length; j++) {
						String individuo = colsTar[j].replace("CONTROL_", "");
						if (headerSec.containsKey(individuo) && headerTar.containsKey(individuo)) {
							int posIndiSec = headerSec.get(individuo);
							int posIndiTar = headerTar.get(individuo);
							
							String dosisSecVal = (posIndiSec < partsSec.length) ? partsSec[posIndiSec] : "";
							String dosisTarVal = (posIndiTar < partsTar.length) ? partsTar[posIndiTar] : "";
							
							if (!dosisTarVal.equals(partsTar[0])) {
								System.out.println(partsTar[0] + "\t" + partsTar[1] + "\t" + individuo + "\t" + dosisSecVal + "\t" + dosisTarVal);
							}
						}
					}
				}
			}
			brTar.close();
			
		} catch (IOException e) {
			System.err.println("Error procesando los archivos: " + e.getMessage());
			e.printStackTrace();
		}
	}
}
