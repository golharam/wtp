package com.lifetechnologies.solid.wt;

import java.io.File;
import java.util.ArrayList;
import java.util.logging.Logger;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;

/**
 * User: tuchbb
 * Date: Sep 29, 2008
 * Time: 2:37:22 PM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ReadMappingExtender {

    public final File fileExtendMappedReadsExe;

    private File fileReferenceFasta = null;
    private File fileMA = null;
    private File fileFullLengthReadsCSFasta = null;

    private static Logger logger = Logger.getLogger(ReadMappingExtender.class.getName());
    
    private int lengthOfRead;


    public ReadMappingExtender(File fileExtendMappedReadsExe,
                               File fileReferenceFasta,
                               File fileMA,
                               File fileFullLengthReadsCSFasta,
                               int lengthOfRead) throws Exception {


        this.fileExtendMappedReadsExe = fileExtendMappedReadsExe.getAbsoluteFile();
        this.fileReferenceFasta = fileReferenceFasta.getAbsoluteFile();
        this.fileMA = fileMA.getAbsoluteFile();
        this.fileFullLengthReadsCSFasta = fileFullLengthReadsCSFasta.getAbsoluteFile();
        this.lengthOfRead = lengthOfRead;

    }

    /**
     *
     *
     * @param validAdjacentMismatchesCountAsOne
     * @param iupacMatchesValid
     * @param clusterInterface
     * @param folderForScriptIO
     * @param fileOutput
     * @param verbose
     * @return
     * @throws Exception
     */
    public File startExtendMappedReadsJob(int scoreOfMatch,
                                              int scoreOfMismatch,
                                              boolean validAdjacentMismatchesCountAsOne,
                                              boolean iupacMatchesValid,
                                              int positionOfMapStartInRead,
                                              ClusterInterface clusterInterface,
                                              String prefixForJobNames,
                                              File folderForTempFilesOnComputeNodes,
                                              File folderForScriptIO,
                                              File fileOutput,
                                              File fileOutputIndex,
                                              boolean verbose,
                                              String maskOfReads) throws Exception {

        if (verbose) {
        	StringBuffer msg = new StringBuffer();
            msg.append("\n---------------------------------------------------------------------------------------------\n");
            msg.append("reference fasta file:\t" + this.fileReferenceFasta.getPath()+"\n");
            msg.append("full length reads csfasta file:\t" + this.fileFullLengthReadsCSFasta.getPath()+"\n");
            msg.append("MA file:\t" + this.fileMA.getPath()+"\n");
            msg.append("output file:\t" + fileOutput.getPath()+"\n");
            msg.append("read length:\t" + this.lengthOfRead+"\n");
            msg.append("match score:\t" + scoreOfMatch+"\n");
            msg.append("mismatch score:\t" + scoreOfMismatch+"\n");
            msg.append("valid adjacent mismatches count as one:\t" + validAdjacentMismatchesCountAsOne+"\n");
            msg.append("iupac matches valid:\t" + iupacMatchesValid+"\n");
            msg.append("position of map start in read:\t" + positionOfMapStartInRead+"\n");
            msg.append("---------------------------------------------------------------------------------------------");
            logger.info(msg.toString());
        }

        ArrayList<String> commandStrings = new ArrayList<String>();

        String command = this.fileExtendMappedReadsExe +
                            " " + this.fileReferenceFasta.getAbsolutePath() +
                            " " + this.fileMA.getAbsolutePath() +
                            " " + this.fileFullLengthReadsCSFasta.getAbsolutePath() +
                            " " + fileOutputIndex.getAbsolutePath() +
                            " " + this.lengthOfRead +
                            " " + positionOfMapStartInRead +
                            " " + scoreOfMatch +
                            " " + scoreOfMismatch;
        
        //Not relying on the #! path anymore.  User must have correct python in $PATH.
        if (this.fileExtendMappedReadsExe.getName().endsWith(".py"))
        	command = "python " + command;

        if (iupacMatchesValid)
            command += " 1 ";
        else
            command += " 0 ";

        if (validAdjacentMismatchesCountAsOne)
            command += " 1 ";
        else
            command += " 0 ";
        
        command += " " + maskOfReads;

        command += " > " + fileOutput.getAbsolutePath();

        commandStrings.add(command);

        File fileScript = new File(folderForScriptIO, prefixForJobNames + fileOutput.getName() + ".sh");

        //String workingFolderName = prefixForJobNames.concat(ProcessId.getProcessId())+"."+System.currentTimeMillis();
        //File workingFolder = new File(folderForTempFilesOnComputeNodes, workingFolderName);
        File workingFolder = Utilities.getUniqueFileName(folderForTempFilesOnComputeNodes, prefixForJobNames);
        clusterInterface.setLocalTemporaryWorkingPath(workingFolder.getAbsolutePath());
        return clusterInterface.executeJob(fileScript, commandStrings);
    }

    public static long calculateMemoryRequiredByExtendMapReadsInBytes(long lengthOfReferenceSequence) throws Exception {
        return (long)(1.2 * lengthOfReferenceSequence);
    }


}
