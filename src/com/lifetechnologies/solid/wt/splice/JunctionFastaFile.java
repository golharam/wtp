package com.lifetechnologies.solid.wt.splice;

import java.io.File;
import java.io.IOException;

import com.lifetechnologies.solid.wt.BufferedRandomAccessFile;

/**
 * A fasta file that contains entries corresponding to regions flanking splice
 * junctions.  Fasta entries should be parseable as FlankJunction objects.
 * @author mullermw
 *
 */
public class JunctionFastaFile {
	
	private BufferedRandomAccessFile braf = null;
	private int numberOfJunctions;

	/**
	 * 
	 * @param file the file containing the junctions.
	 * @throws IOException
	 */
	public JunctionFastaFile(File file) throws IOException {
		this.braf = new BufferedRandomAccessFile(file, "r");
		if (braf.length() > 1000) braf.seek(braf.length() - 1000);
		for (String line = braf.readLine(); line != null; line=braf.readLine()) {
			if (!line.startsWith(">")) continue;
			line = line.substring(1);
			String[] fields = line.split("\\s+");
			numberOfJunctions = parseIntFromEndOfString(fields[0]);
		}
	}
	
	/** 
	 * Retrieves the flanked junction in the file by junction number.
	 * Uses binary search for faster random access.
	 * @param num
	 * @return
	 * @throws IOException
	 */
	public FlankedJunction getFlankedJunction(int num) throws IOException {
		if (num < 1) throw new IllegalArgumentException("num cannot be less than 1.");
		if (num > numberOfJunctions) throw new IllegalArgumentException("Attempted to get Junction #"+num+" when the number of junctions is "+numberOfJunctions+".");
		//Binary search for the junction by number.
		int nextJunctionNumber = Integer.MAX_VALUE;
		long left = 0;
		long right = braf.length() - 1;
		String descLine = null;
		while (descLine == null ) {
			long i = (left + right) / 2;
			braf.seek(i);
			for (String line = braf.readNextLine(); line!= null; line = braf.readNextLine()) {
				if (!line.startsWith(">")) continue;
				nextJunctionNumber = parseIntFromEndOfString(line.split("\\s+")[0]);
				if (nextJunctionNumber == num) {
					descLine = line;
				} else if (nextJunctionNumber > num){
					right = i;
				} else {
					left = i;
				}
				break;
			}
			i = (left + right) / 2;
		}
		FlankedJunction junction = FlankedJunction.parseFlankedJunction(descLine.split("\\s+")[1]);
		String sequence = braf.readNextLine();
		junction.getLeftSequence().append(sequence.substring(0,junction.getLeftEnd()-junction.getLeftStart()+1));
		junction.getRightSequence().append(sequence.substring(junction.getLeftEnd()-junction.getLeftStart()+1));
		return junction;
	}
	
	public void close() throws IOException {
		braf.close();
	}
	
	@Override
	protected void finalize() throws Throwable {
		braf.close();
	}
	
	public static void main(String[] args) throws IOException {
		JunctionFastaFile jff = new JunctionFastaFile(new File("//siena/mullermw/test/1.2_testing/output/tmp/junctions.fa"));
		for (int i=1; i<=9166; i++) {
			System.out.println(i);
			System.out.println(jff.getFlankedJunction(i));
		}
	}
	
	private static int parseIntFromEndOfString(String s) {
		int i = s.length()-1;
		for (;i>0 && Character.isDigit(s.charAt(i)); i--);
		if (s.charAt(i) == '-') i--;
		return Integer.parseInt(s.substring(i+1));
	}
}
