package com.lifetechnologies.solid.wt.cmd;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Properties;

import org.apache.commons.cli.Options;

import com.lifetechnologies.solid.wt.config.ComplexOption;
import com.lifetechnologies.solid.wt.config.ConfigKey;
import com.lifetechnologies.solid.wt.config.validator.CreateableFileValidator;
import com.lifetechnologies.solid.wt.config.validator.MaxValuesValidator;
import com.lifetechnologies.solid.wt.config.validator.ReadableFileValidator;
import com.lifetechnologies.solid.wt.config.validator.RequiredArgumentValidator;
import com.lifetechnologies.util.RefGene2Gtf;

public class RefGene2GffCmd extends AbstractCmd {

	public RefGene2GffCmd() {
		super("refgene2gff");
	}
	
	@Override
	public void addOptions(Options options) {
		ComplexOption option = new ComplexOption(ConfigKey.WT_REFGENE_2_GFF_REFGENE_FILE, true);
		option.addValidator(new RequiredArgumentValidator());
		option.addValidator(new MaxValuesValidator(1));
		option.addValidator(new ReadableFileValidator(option));
		options.addOption(option);
		
		option = new ComplexOption(ConfigKey.WT_REFGENE_2_GFF_FILE, true);
		option.addValidator(new MaxValuesValidator(1));
		option.addValidator(new CreateableFileValidator(option));
		options.addOption(option);
	}

	@Override
	public void runCmd() throws Exception {
		File refGeneFile = new File(config.get(ConfigKey.WT_REFGENE_2_GFF_REFGENE_FILE).iterator().next()).getAbsoluteFile();
		InputStream refGene = new FileInputStream(refGeneFile);
		PrintStream gff;
		if (config.containsKey(ConfigKey.WT_REFGENE_2_GFF_FILE)) {
			gff = new PrintStream(new File(config.get(ConfigKey.WT_REFGENE_2_GFF_FILE).iterator().next()));
		} else {
			gff = System.out;
		}
		Properties header = new Properties();
		header.put(RefGene2Gtf.SOURCE_FILE_HEADER_KEY, refGeneFile.toString());
		RefGene2Gtf refGene2Gff = new RefGene2Gtf(refGene, gff, header);
		refGene2Gff.run();
	}

}
