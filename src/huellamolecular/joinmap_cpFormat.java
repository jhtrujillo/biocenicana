package huellamolecular;
public class joinmap_cpFormat {

	public void fixformat(String archivoNGSEP) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(archivoNGSEP);
		
		for (int i=0; i<ar.numerolineas; i++) {
			
			String[] fila = datos[i].split("	");
			
			for (int j = 0; j< fila.length; j++) {
				if (j+1==fila.length) {
					System.out.print(fila[j]); 
				}
				else if (j==2) {
					System.out.print("\t");
					System.out.print(fila[j]+"\t");
				}
				else {
					System.out.print(fila[j]+"\t");
				}
			}
			System.out.println(); 
			
			
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		joinmap_cpFormat joinmap = new joinmap_cpFormat();
		joinmap.fixformat("/home/estuvar4/Downloads/tmp.txt");
	}

}
