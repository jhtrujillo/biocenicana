package org.cenicana.bio;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class FileUtils2 {
	public File archivo;
	public FileReader fr;
	public BufferedReader br;
	public int numerolineas=0;
	public String[] fichero = new String[999999999];
	public double sizefile = 0;
 
	
	public double getsize(){
		return this.sizefile;  
	}
	
	public String[] leerfichero(String path){
	
		//String[] fichero=new String[900000];
	
		 try {
	         // Apertura del fichero y creacion de BufferedReader para poder
	         // hacer una lectura comoda (disponer del metodo readLine()).
	         archivo = new File (path);
	         fr = new FileReader (archivo);
	         br = new BufferedReader(fr);
	         
	         this.sizefile = archivo.getUsableSpace();
	
	         // Lectura del fichero
	         String linea;
	         int cont=0;
	         while((linea=br.readLine())!=null){ 
	        	 	this.fichero[cont]=linea;
		            numerolineas++;
		            cont++;
	         }
	
	      }
	      catch(Exception e){
	         e.printStackTrace();
	      }finally{
	         // En el finally cerramos el fichero, para asegurarnos
	         // que se cierra tanto si todo va bien como si salta
	         // una excepcion.
	         try{
	            if( null != fr ){
	               fr.close();
	            }
	         }catch (Exception e2){
	            e2.printStackTrace();
	         }
	      }
		return this.fichero;
	   }

	
	public String[]  leerfichero2(String path){
		try {
	         @SuppressWarnings("resource")
			BufferedReader in = new BufferedReader
	         (new FileReader(path));
	         String str;
	         while ((str = in.readLine()) != null) {
	        	this.fichero[numerolineas]=str;
	            numerolineas++;
	            
	         }
	         }
	         catch (IOException e) {
	         }
			return this.fichero;
	 }
	
	public void imprimirDatosFichero(){
		for(int i = 0; i < this.numerolineas ; i++){
			System.out.println(this.fichero[i]);
		}
	}
	
	  public int numerolineas(){
	        return this.numerolineas;
	  }
	
	public void crearfichero(String path){
	
		File f;
	    f=new File(path);
	    if(!f.exists()){
	    	try {
				f.createNewFile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	    	//System.out.println("has been createdto the current directory");
	    }
	}
	
	public int  creardirectorio(String path){
		int r=0;
		File Carpeta = new File(path);
		if (Carpeta.exists()){
			System.out.println("El directorio ya existe");
	
		}
		else{
			Carpeta.mkdir();
			r=1;
		}
		return r;
	
	}
	
	public void addlineafichereo(String path, String linea){
	
		try
	    {
	            BufferedWriter out = new BufferedWriter(new FileWriter(path, true));
	            out.write(linea);
	            out.close();
	    } catch (IOException e)
	    {
	    }
	
	}
	
	public int compararArchivos(String fileini, String filefin){
		File file1 = new File(fileini);
	    File file2 = new File(filefin);
	    if(file1.compareTo(file2) == 0) {
	         System.out.println("Both paths are same!");
	      } else {
	         System.out.println("Paths are not same!");
	      }
	    return file1.compareTo(file2); 
	}
	
	public void DuplicarFichero(String fileOriginal, String fileDestino)throws IOException {
		
		File source = new File(fileOriginal);
		File dest = new File(fileDestino);
		
		InputStream is = null;
	    OutputStream os = null;
	    try {
	        is = new FileInputStream(source);
	        os = new FileOutputStream(dest);
	        byte[] buffer = new byte[1024];
	        int length;
	        while ((length = is.read(buffer)) > 0) {
	            os.write(buffer, 0, length);
	        }
	    } finally {
	        is.close();
	        os.close();
	    }
	}
	
	public int existefichero(String path){
		int salida=0;
		File fichero = new File(path);
		if (fichero.exists())
			salida=1;
	
		return salida;
	}
	
	public void eliminarfichero(String path){
	    File fichero = new File(path);
	    fichero.delete();
	
	}
	
	public void eliminarfichero2(String path){
	    File fichero = new File(path);
	    fichero.delete();
	 
	}
	
	public void renameFile(String path, String name){
		File fichero = new File(path);
		File newname = new File(name);
		fichero.renameTo(newname);
	}
	
	
}