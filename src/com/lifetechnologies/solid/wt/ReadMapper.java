package com.lifetechnologies.solid.wt;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;

/**
 * User: tuchbb
 * Date: Aug 14, 2008
 * Time: 9:20:47 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ReadMapper {

    public final File fileMapReadsExe;

    private File fileReferenceFasta = null;
    private File fileReadsCSFasta = null;

    private int lengthOfRead;
    private String maskPattern = "";

    private static Logger logger = Logger.getLogger(ReadMapper.class.toString());
    
    public ReadMapper(File fileMapReadsExe,
                      File fileReferenceFasta,
                      File filesReadsCSFasta,
                      int lengthOfRead,
                      String maskPattern) throws Exception {

        if (maskPattern.length() != 0 && maskPattern.length() != lengthOfRead)
            throw new Exception("Mask pattern must be of length equal to tag length: '"+maskPattern+"'("+maskPattern.length()+") != "+lengthOfRead);

        this.fileMapReadsExe = fileMapReadsExe;
        this.fileReferenceFasta = fileReferenceFasta;
        this.fileReadsCSFasta = filesReadsCSFasta;
        this.lengthOfRead = lengthOfRead;
        this.maskPattern = maskPattern;
    }

    public ReadMapper(File fileMapReadsExe,
                      File fileReferenceFasta,
                      File filesReadsCSFasta,
                      int lengthOfRead,
                      int countMaskLeft,
                      int countMaskRight) throws Exception {


        if ((countMaskLeft + countMaskRight) > lengthOfRead)
            throw new Exception("Mask pattern would be longer than tag length: [countMaskLeft + countMaskRight > lengthOfRead]");

        this.fileMapReadsExe = fileMapReadsExe;
        this.fileReferenceFasta = fileReferenceFasta;
        this.fileReadsCSFasta = filesReadsCSFasta;
        this.lengthOfRead = lengthOfRead;

        for (int i = 0; i < countMaskLeft; i++)
            this.maskPattern += "0";

        for (int i = 0; i < lengthOfRead - countMaskLeft - countMaskRight; i++)
            this.maskPattern += "1";

        for (int i = 0; i < countMaskRight; i++)
            this.maskPattern += "0";

    }

    public void mapReads(int numberOfJobs, int countMaxErrorsAllowedWhenMatching, boolean validAdjacentMismatchesCountAsOne, boolean iupacMatchesValid, int countMaxMatchingLocationsToReport,
                         ClusterInterface clusterInterface, String prefixForJobNames,
                         File folderForTempFilesOnComputeNodes, File folderContainingSchemaFiles, File folderForScriptIO, File folderForReadMappingOutput, boolean verbose) throws Exception {

        FastaDatabase referenceFastaDB = new FastaDatabase(fileReferenceFasta);
        int numberOfRefrenceSequences = referenceFastaDB.getNumberOfSequencesInDatabase();

        FastaDatabase readsFastaDB = new FastaDatabase(fileReadsCSFasta);
        int numberOfReadSequences = readsFastaDB.getNumberOfSequencesInDatabase();

        startMapReadsJob((numberOfRefrenceSequences == 1), numberOfReadSequences, numberOfJobs,
                        countMaxErrorsAllowedWhenMatching,  validAdjacentMismatchesCountAsOne, iupacMatchesValid, countMaxMatchingLocationsToReport, clusterInterface, prefixForJobNames,
                        folderForTempFilesOnComputeNodes, folderContainingSchemaFiles, folderForScriptIO, folderForReadMappingOutput, verbose);
    }

    public File startMapReadsJob(boolean isMultiFastaReference,
                                     int indexOfFirstRead,
                                     int indexOfLastRead,
                                     int countMaxErrorsAllowedWhenMatching,
                                     boolean validAdjacentMismatchesCountAsOne,
                                     boolean iupacMatchesValid,
                                     int countMaxMatchingLocationsToReport,
                                     ClusterInterface clusterInterface,
                                     String prefixForJobNames,
                                     File folderForTempFilesOnComputeNodes,
                                     File folderContainingSchemaFiles,
                                     File folderForScriptIO,
                                     File fileOutput,
                                     boolean verbose) throws Exception {

        if (verbose) {
        	StringBuffer msg = new StringBuffer();
        	msg.append("\n");
            msg.append("---------------------------------------------------------------------------------------------\n");
            msg.append("reference fasta file:\t" + fileReferenceFasta.getPath()+"\n");
            msg.append("reads csfasta file:\t" + fileReadsCSFasta.getPath()+"\n");
            msg.append("index of first read:\t" + indexOfFirstRead+"\n");
            msg.append("index of last read:\t" + indexOfLastRead+"\n");
            msg.append("output file:\t" + fileOutput.getPath()+"\n");
            msg.append("read length:\t" + lengthOfRead+"\n");
            msg.append("max errors allowed when matching:\t" + countMaxErrorsAllowedWhenMatching+"\n");
            msg.append("max matches returned:\t" + countMaxMatchingLocationsToReport+"\n");
            msg.append("valid adjacent mismatches count as one:\t" + validAdjacentMismatchesCountAsOne+"\n");
            msg.append("iupac matches valid:\t" + iupacMatchesValid+"\n");
            msg.append("---------------------------------------------------------------------------------------------\n");
          	logger.info(msg.toString()); 
        }

        ArrayList<String> commandStrings = new ArrayList<String>();
        File fileSchema = getSchemaFile(maskPattern, countMaxErrorsAllowedWhenMatching, validAdjacentMismatchesCountAsOne, folderContainingSchemaFiles);
        if (!fileSchema.exists())
            throw new Exception("Schema file " + fileSchema.getPath() + " could not be found.");
        File fileSchemaMasked = new File(folderForScriptIO.getParent() + "/" + fileSchema.getName() + ".masked");
        generateMaskedSchemaFile(fileSchema, this.maskPattern, fileSchemaMasked);

        String command = this.fileMapReadsExe +
                            " " + this.fileReadsCSFasta.getAbsolutePath() +
                            " " + this.fileReferenceFasta.getAbsolutePath() +
                            " T=" + fileSchemaMasked.getAbsolutePath() +
                            " S=0 " +   // align in color space
                            " O=0 " +   // offset?
                            " u=2 " +   // what is this?  carried over from fasta2match
                            " R=0 " +   // what is this?  carried over from fasta2match
                            " L=" + this.lengthOfRead +
                            " M=" + countMaxErrorsAllowedWhenMatching +
                            " Z=" + countMaxMatchingLocationsToReport +
                            " b=" + indexOfFirstRead +   // first read to map
                            " e=" + indexOfLastRead; // last read to map

        if (validAdjacentMismatchesCountAsOne)
            command += " A=1 ";

        if (iupacMatchesValid)
            command += " H=1 ";

        if (this.maskPattern.length() > 0)
            command += " X=" + this.maskPattern;

        if (isMultiFastaReference)
            command += " I=1 ";

        command += " > " + fileOutput.getAbsolutePath();

        commandStrings.add(command);

        File fileScript = new File(folderForScriptIO.getPath() + "/" + prefixForJobNames + fileOutput.getName() + ".sh");
        
        //String workingFolderName = prefixForJobNames.concat(ProcessId.getProcessId())+"."+System.currentTimeMillis();
        //File workingFolder = new File(folderForTempFilesOnComputeNodes, workingFolderName);
        File workingFolder = Utilities.getUniqueFileName(folderForTempFilesOnComputeNodes, prefixForJobNames);
        clusterInterface.setLocalTemporaryWorkingPath(workingFolder.getAbsolutePath());

        return clusterInterface.executeJob(fileScript, commandStrings);
    }



    /**
     * 
     * @param fileSchema
     * @param maskPattern
     * @param fileSchemaMasked
     * @throws IOException
     */
    private static void generateMaskedSchemaFile(File fileSchema, 
                                                  String maskPattern,
                                                  File fileSchemaMasked) throws IOException {

        char charsMaskPattern[] = maskPattern.substring(1).toCharArray();
        BufferedReader readerSchemaFile = new BufferedReader(new FileReader(fileSchema));
        BufferedWriter writerSchemaFileMasked = new BufferedWriter(new FileWriter(fileSchemaMasked));
        String line;
        while ((line = readerSchemaFile.readLine()) != null) {
            if (line.startsWith("#"))
                writerSchemaFileMasked.write(line + " # MASKED " +  maskPattern);
            else {
                int indexInSchema = 0;
                char charsSchema[] = line.trim().toCharArray();
                for (int i = 0; i < charsMaskPattern.length; i++) {
                    if (charsMaskPattern[i] == '1')
                        writerSchemaFileMasked.write("" + charsSchema[indexInSchema++]);
                    else
                        writerSchemaFileMasked.write("" + charsMaskPattern[i]);
                }
            }
            writerSchemaFileMasked.newLine();
        }
        readerSchemaFile.close();
        writerSchemaFileMasked.close();
    }
    public static void main(String[] args ) throws Exception {
    	File inFile = new File("etc/schemas/schema_26_2_adj");
    	File outFile = new File("schema");
    	generateMaskedSchemaFile(inFile, "0111111111111111111111111100000", outFile);
    	BufferedReader reader = new BufferedReader(new FileReader(outFile));
    	for (String line=reader.readLine(); line != null; line = reader.readLine()) {
    		System.out.println(line);
    	}
    	outFile.delete();
    }

    static long getMinimumMemoryRequirementForMapreads(long longestSequenceLength) {
    	return (long)(1.3e9 + 7 * longestSequenceLength);
    }
    
    public static File getSchemaFile(String maskPattern, int countMaxErrorsAllowedWhenMatching,
    		boolean validAdjacentMismatchesCountAsOne,
    		File folderContainingSchemaFiles ) {
    	
        int lengthOfUnmaskedSchema = maskPattern.length() - maskPattern.substring(1).replaceAll("1", "").length();
        String schemaFileName = "schema_" + lengthOfUnmaskedSchema + "_" + countMaxErrorsAllowedWhenMatching;
        if (validAdjacentMismatchesCountAsOne && countMaxErrorsAllowedWhenMatching > 1) schemaFileName += "_adj";
        return new File(folderContainingSchemaFiles, schemaFileName );
    }
    
    public static long calculateMaxReferenceSizeForMapReads(long memoryAllocatedInBytes) {
        return (long)(Math.min(1E9, Math.floor((memoryAllocatedInBytes - 1.3E9) / 7.0)));
    }
    
    public static long calculateMemoryRequiredByMapReadsInBytes(long lengthOfReferenceSequence) throws Exception {
        if (lengthOfReferenceSequence > 1E9)
            throw new Exception("Reference sequence is too long for mapReads (" + lengthOfReferenceSequence + " > " + 1E9 + ")");
        //return (1.1E9 + 6 * lengthOfReferenceSequence);
        return (long)((1.3E9 + 7.0 * lengthOfReferenceSequence));
    }
}
