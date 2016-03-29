/*import org.biojava.bio.symbol.*;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;

import java.util.Set;
import java.util.HashSet;
import java.util.ArrayList;
    */
/**
 * User: tuchbb
 * Date: Aug 13, 2008
 * Time: 3:30:21 PM 
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
/*public class SOLiDAlphabet extends SimpleAlphabet {

    CharacterTokenization solidTokenization = null;

    public SOLiDAlphabet() throws BioException {
        super("SOLiD");

        this.solidTokenization = new CharacterTokenization(this, false);
        //this.putTokenization(this.getName(), this.solidTokenization);
        this.putTokenization("token", this.solidTokenization);

        AtomicSymbol newSymbol = AlphabetManager.createSymbol("A");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, 'A');

        newSymbol = AlphabetManager.createSymbol("C");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, 'C');

        newSymbol = AlphabetManager.createSymbol("G");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, 'G');

        newSymbol = AlphabetManager.createSymbol("T");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, 'T');

        newSymbol = AlphabetManager.createSymbol("0");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, '0');

        newSymbol = AlphabetManager.createSymbol("1");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, '1');

        newSymbol = AlphabetManager.createSymbol("2");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, '2');

        newSymbol = AlphabetManager.createSymbol("3");
        this.addSymbol(newSymbol);
        solidTokenization.bindSymbol(newSymbol, '3');

        
        AlphabetManager.registerAlphabet(this.getName(), this);
    }
}
  */