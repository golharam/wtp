package com.lifetechnologies.solid.wt;

import java.util.TreeSet;
import java.util.Iterator;
import java.util.Comparator;

/**
 * User: tuchbb
 * Date: Oct 1, 2008
 * Time: 11:26:58 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ExtendedReadMappingSetForRead {

    private String idOfRead;
    private TreeSet<ExtendedReadMapping> setOfSortedReadMappings;

    /**
     *
     * @param headerWithoutGTSign a fasta header produced by the mapped read extension module
     */
    public ExtendedReadMappingSetForRead(String headerWithoutGTSign, Comparator<ExtendedReadMapping> comparatorOfExtendingReadMappings) {

        this.setOfSortedReadMappings = new TreeSet<ExtendedReadMapping>(comparatorOfExtendingReadMappings);

        String tokensForHeader[] = headerWithoutGTSign.split(",");
        this.idOfRead = tokensForHeader[0];
        for (int i = 1; i < tokensForHeader.length; i++) {
            ExtendedReadMapping extendedReadMapping = parseSingleMappingBlockOfHeader(this.idOfRead, tokensForHeader[i]);
            setOfSortedReadMappings.add(extendedReadMapping);
        }

    }

    public void addMoreMappingsFromHeader(String headerWithoutGTSign) throws Exception {

        String tokensForHeader[] = headerWithoutGTSign.split(",");
        if (tokensForHeader[0] != this.idOfRead)
            throw new Exception("Read Ids do not match: " + tokensForHeader[0] + " != " + this.idOfRead);

        for (int i = 1; i < tokensForHeader.length; i++) {
            ExtendedReadMapping extendedReadMapping = parseSingleMappingBlockOfHeader(this.idOfRead, tokensForHeader[i]);
            setOfSortedReadMappings.add(extendedReadMapping);
        }
    }

    public void addMapping(ExtendedReadMapping extendedReadMapping) {
        this.setOfSortedReadMappings.add(extendedReadMapping);
    }

    public int size() {
        return setOfSortedReadMappings.size();
    }

    public String getIdOfRead() {
        return idOfRead;
    }

    public TreeSet<ExtendedReadMapping> getSetOfSortedReadMappings() {
        return setOfSortedReadMappings;
    }

    public String toString() {

        StringBuffer header = new StringBuffer(this.idOfRead);
        Iterator<ExtendedReadMapping> iteratorOverMappings = this.setOfSortedReadMappings.iterator();
        while (iteratorOverMappings.hasNext()) {
            ExtendedReadMapping extendedReadMapping = iteratorOverMappings.next();
            header.append(extendedReadMapping.toString());            
        }

        return header.toString();
    }

    public static ExtendedReadMapping[] parseExtendedReadMappingsFromHeader(String headerWithoutGTSign) {
        String tokensForHeader[] = headerWithoutGTSign.split(",");
        ExtendedReadMapping extendedReadMappings[] = new ExtendedReadMapping[tokensForHeader.length -1];
        for (int i = 0; i < extendedReadMappings.length; i++)
            extendedReadMappings[i] = parseSingleMappingBlockOfHeader(tokensForHeader[0], tokensForHeader[i +1]);
        return extendedReadMappings;
    }

    private static ExtendedReadMapping parseSingleMappingBlockOfHeader(String idOfRead, String blockOfHeader) {

        String tokensForMapping[] = blockOfHeader.split("#");
        String tokensForMappingBasicInfo[] = tokensForMapping[0].split("\\.");

        ExtendedReadMapping extendedReadMapping = new ExtendedReadMapping(idOfRead,
                                                                          tokensForMappingBasicInfo[0],
                                                                          Integer.parseInt(tokensForMappingBasicInfo[0]),
                                                                          Integer.parseInt(tokensForMappingBasicInfo[1]),
                                                                          Integer.parseInt(tokensForMappingBasicInfo[5]),
                                                                          Integer.parseInt(tokensForMappingBasicInfo[4]),
                                                                          Integer.parseInt(tokensForMappingBasicInfo[2]),
                                                                          Integer.parseInt(tokensForMappingBasicInfo[3])) ;

        if (tokensForMapping.length >= 2) {
            String tokensForMappingValidAdjacents[] = tokensForMapping[1].split("\\.");
            for (int j = 0; j < tokensForMappingValidAdjacents.length; j++)
                extendedReadMapping.addValidAdjacenStartPosition(Integer.parseInt(tokensForMappingValidAdjacents[j]));
        }
        if (tokensForMapping.length >= 3) {
            String tokensForMappingInsertionsInRead[] = tokensForMapping[2].split(";");
            for (int j = 0; j < tokensForMappingInsertionsInRead.length; j++) {
                String tokensForInsertion[] = tokensForMappingInsertionsInRead[j].split("\\.");
                extendedReadMapping.addInsertionInRead(new Insertion(Integer.parseInt(tokensForInsertion[0]), Integer.parseInt(tokensForInsertion[1])));
            }
        }
        if (tokensForMapping.length >= 4) {
            String tokensForMappingInsertionsInReference[] = tokensForMapping[3].split(";");
            for (int j = 0; j < tokensForMappingInsertionsInReference.length; j++) {
                String tokensForInsertion[] = tokensForMappingInsertionsInReference[j].split("\\.");
                extendedReadMapping.addInsertionInReference(new Insertion(Integer.parseInt(tokensForInsertion[0]), Integer.parseInt(tokensForInsertion[1])));
            }
        }
        return extendedReadMapping;
    }


}
