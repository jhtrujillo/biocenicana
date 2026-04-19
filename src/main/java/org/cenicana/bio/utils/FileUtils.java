package org.cenicana.bio.utils;
import org.cenicana.bio.core.JoinMapCpFormat;
import org.cenicana.bio.core.AlleleDosageCalculator;
import org.cenicana.bio.core.VcfStatisticsCalculator;
import org.cenicana.bio.core.VcfFilter;

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
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class FileUtils {
    public File archivo;
    public FileReader fr;
    public BufferedReader br;
    public int numerolineas = 0;
    // Se elimina la inicialización con tamaño fijo
    public String[] fichero;
    public double sizefile = 0;

    public double getsize(){
        return this.sizefile;  
    }
    
    /**
     * @deprecated Usar getLineStream() o getBufferedReader() para evitar OutOfMemoryError con archivos grandes.
     * Lee el archivo de texto ubicado en 'path' y retorna un arreglo con todas sus líneas.
     * Se actualizan la variable numerolineas y fichero.
     */
    @Deprecated
    public String[] leerfichero(String path) {
        List<String> lineas = new ArrayList<>();
        try {
            archivo = new File(path);
            fr = new FileReader(archivo);
            br = new BufferedReader(fr);
            
            this.sizefile = archivo.getUsableSpace();
            
            String linea;
            while((linea = br.readLine()) != null){ 
                lineas.add(linea);
            }
        } catch(Exception e) {
            e.printStackTrace();
        } finally {
            try {
                if (fr != null) {
                    fr.close();
                }
            } catch (Exception e2) {
                e2.printStackTrace();
            }
        }
        // Se actualiza el número de líneas y se convierte la lista en arreglo
        this.numerolineas = lineas.size();
        this.fichero = lineas.toArray(new String[0]);
        return this.fichero;
    }

    /**
     * @deprecated Usar getLineStream() o getBufferedReader() para evitar OutOfMemoryError con archivos grandes.
     * Variante para leer el archivo y retornar un arreglo con todas sus líneas.
     */
    @Deprecated
    public String[] leerfichero2(String path) {
        List<String> lineas = new ArrayList<>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(path));
            String str;
            while ((str = in.readLine()) != null) {
                lineas.add(str);
            }
            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Se actualiza el número de líneas y se convierte la lista en arreglo
        this.numerolineas = lineas.size();
        this.fichero = lineas.toArray(new String[0]);
        return this.fichero;
    }
    
    /**
     * Devuelve un Stream de líneas del archivo. El recurso debe cerrarse usándose en un bloque try-with-resources.
     */
    public Stream<String> getLineStream(String path) throws IOException {
        return Files.lines(Paths.get(path));
    }
    
    /**
     * Devuelve un BufferedReader del archivo asignado. El recurso debe cerrarse.
     */
    public BufferedReader getBufferedReader(String path) throws IOException {
        return Files.newBufferedReader(Paths.get(path));
    }
    
    public void imprimirDatosFichero() {
        for (int i = 0; i < this.numerolineas; i++) {
            System.out.println(this.fichero[i]);
        }
    }
    
    public int numerolineas() {
        return this.numerolineas;
    }
    
    public void crearfichero(String path) {
        File f = new File(path);
        if (!f.exists()) {
            try {
                f.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    public int creardirectorio(String path) {
        int r = 0;
        File Carpeta = new File(path);
        if (Carpeta.exists()) {
            System.out.println("El directorio ya existe");
        } else {
            Carpeta.mkdir();
            r = 1;
        }
        return r;
    }
    
    public void addlineafichereo(String path, String linea) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path, true));
            out.write(linea);
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public int compararArchivos(String fileini, String filefin) {
        File file1 = new File(fileini);
        File file2 = new File(filefin);
        if(file1.compareTo(file2) == 0) {
            System.out.println("Both paths are same!");
        } else {
            System.out.println("Paths are not same!");
        }
        return file1.compareTo(file2);
    }
    
    public void DuplicarFichero(String fileOriginal, String fileDestino) throws IOException {
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
            if (is != null) is.close();
            if (os != null) os.close();
        }
    }
    
    public int existefichero(String path) {
        int salida = 0;
        File fichero = new File(path);
        if (fichero.exists())
            salida = 1;
        return salida;
    }
    
    public void eliminarfichero(String path) {
        File fichero = new File(path);
        fichero.delete();
    }
    
    public void eliminarfichero2(String path) {
        File fichero = new File(path);
        fichero.delete();
    }
    
    public void renameFile(String path, String name) {
        File fichero = new File(path);
        File newname = new File(name);
        fichero.renameTo(newname);
    }
}
