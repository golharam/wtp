package com.lifetechnologies.solid.wt;

import java.io.*;
import java.util.zip.*;
import java.util.Enumeration;
import java.util.Date;

/**
 * User: tuchbb
 * Date: Oct 13, 2008
 * Time: 3:16:46 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class CompressionUtilities {

    public static int SIZE_OF_BUFFER = 1024;

    public static void compressFiles(File files[], File fileCompressed) throws IOException {
        
        byte[] buffer = new byte[SIZE_OF_BUFFER];
        ZipOutputStream zipOutputStream = new ZipOutputStream(new FileOutputStream(fileCompressed));
        for (int i = 0; i < files.length; i++) {
            FileInputStream fileInputStream = new FileInputStream(files[i]);
            zipOutputStream.putNextEntry(new ZipEntry(files[i].getPath()));

            int size;
            while ((size = fileInputStream.read(buffer)) > 0)
                zipOutputStream.write(buffer, 0, size);

            zipOutputStream.closeEntry();
            fileInputStream.close();
        }

        zipOutputStream.close();
    }

    public static void compressFile(File file, File fileCompressed) throws IOException {

        byte[] buffer = new byte[SIZE_OF_BUFFER];
        ZipOutputStream zipOutputStream = new ZipOutputStream(new FileOutputStream(fileCompressed));
        FileInputStream fileInputStream = new FileInputStream(file);
        zipOutputStream.putNextEntry(new ZipEntry(file.getPath()));

        int size;
        while ((size = fileInputStream.read(buffer)) > 0)
            zipOutputStream.write(buffer, 0, size);

        zipOutputStream.closeEntry();
        fileInputStream.close();
        zipOutputStream.close();
    }

    public static void decompressFiles(File fileCompressed, File folderDestination) throws IOException {

        ZipFile zipFile = new ZipFile(fileCompressed);
        Enumeration<? extends ZipEntry> enumZipFileEntries = zipFile.entries();
        while (enumZipFileEntries.hasMoreElements()) {
            ZipEntry zipEntry = enumZipFileEntries.nextElement();
            BufferedInputStream bufferedInputStream = new BufferedInputStream(zipFile.getInputStream(zipEntry));
            int size;
            byte[] buffer = new byte[SIZE_OF_BUFFER];
            File fileListedInZip = new File(zipEntry.getName());
            FileOutputStream fileOutputStream = new FileOutputStream(folderDestination.getPath() + "/" + fileListedInZip.getName());
            BufferedOutputStream bufferedOutputStream = new BufferedOutputStream(fileOutputStream, buffer.length);
            while ((size = bufferedInputStream.read(buffer, 0, buffer.length)) != -1)
                bufferedOutputStream.write(buffer, 0, size);

            bufferedOutputStream.flush();
            bufferedOutputStream.close();
            fileOutputStream.close();
            bufferedInputStream.close();
        }
    }   

    public static void main(String args[]) throws IOException {
        File files[] = new File[]{ new File("/home/tuchbb/projects/MAQC/mapping_results/brain/match_genome/tmp/reference.0.fasta"),
                                   new File("/home/tuchbb/projects/MAQC/mapping_results/brain/match_genome/tmp/reference.1.fasta")};
        File fileZip = new File("/home/tuchbb/projects/MAQC/mapping_results/reference.zip");

        System.out.println("\nStarting at " + new Date(System.currentTimeMillis()));
        CompressionUtilities.compressFile(files[0], fileZip);
        System.out.println("\nFinished at " + new Date(System.currentTimeMillis()));
        CompressionUtilities.decompressFiles(fileZip, new File("/home/tuchbb/projects/MAQC/mapping_results/"));
    }
}
