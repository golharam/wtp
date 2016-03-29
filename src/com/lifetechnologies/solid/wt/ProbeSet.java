package com.lifetechnologies.solid.wt;

import java.util.HashMap;

/**
 * User: tuchbb
 * Date: Dec 9, 2008
 * Time: 11:17:27 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class ProbeSet {

    private String id;
    private String idsRefSeq;
    private String idRepresentativePublic;
    private String nameOfGene;
    private String idChromosomeLocation;

    private HashMap<String, Double> mapSampleIdToValue;
    private HashMap<String, String> mapSampleIdToCall;

    public ProbeSet(String id, String idsRefSeq, String idRepresentativePublic, String nameOfGene, String idChromosomeLocation) {
        this.id = id;
        this.idsRefSeq = idsRefSeq;
        this.idRepresentativePublic = idRepresentativePublic;
        this.nameOfGene = nameOfGene;
        this.idChromosomeLocation = idChromosomeLocation;

        this.mapSampleIdToValue = new HashMap<String, Double>();
        this.mapSampleIdToCall = new HashMap<String, String>();
    }


    public void setValue(String idSample, Double value) {
        this.mapSampleIdToValue.put(idSample, value);
    }

    public Double getValue(String idSample) {
        return this.mapSampleIdToValue.get(idSample);
    }

    public void setCall(String idSample, String stringCall) {
        this.mapSampleIdToCall.put(idSample, stringCall);
    }

    public String getCall(String idSample) {
        return this.mapSampleIdToCall.get(idSample);
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getIdsRefSeq() {
        return idsRefSeq;
    }

    public void setIdsRefSeq(String idsRefSeq) {
        this.idsRefSeq = idsRefSeq;
    }

    public String getIdRepresentativePublic() {
        return idRepresentativePublic;
    }

    public void setIdRepresentativePublic(String idRepresentativePublic) {
        this.idRepresentativePublic = idRepresentativePublic;
    }

    public String getNameOfGene() {
        return nameOfGene;
    }

    public void setNameOfGene(String nameOfGene) {
        this.nameOfGene = nameOfGene;
    }

    public String getIdChromosomeLocation() {
        return idChromosomeLocation;
    }

    public void setIdChromosomeLocation(String idChromosomeLocation) {
        this.idChromosomeLocation = idChromosomeLocation;
    }
}

