/**
 * 
 */
package japsa.tools.seq;

import java.io.IOException;
import java.util.ArrayList;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * @author minhduc
 * 
 */
@Deployable(
	scriptName = "jsa.seq.extract", 
	scriptDesc = "Extract subsequences"
)
public class SequenceExtract extends CommandLine {
	public SequenceExtract(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options] <chr:start-end> <chr:start-end> ...");
		setDesc(annotation.scriptDesc()); 
		///////////////////////////////////////////////////////////////
		
		addStdInputFile();
		addStdOutputFile();		
		
		addStdAlphabet();
		addBoolean("reverse", false , "Reverse complement the subsequence");
		
		addString("format", "fasta",
				"format of the output file (jsa and fasta)");
		
		addStdHelp();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new SequenceExtract();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String inputFile = cmdLine.getStringVal("input");
		String outputFile = cmdLine.getStringVal("output");
		String format = cmdLine.getStringVal("format").toLowerCase();
		Alphabet alphabet = Alphabet.getAlphabet(cmdLine.getStringVal("alphabet"));
		boolean rev = cmdLine.getBooleanVal("reverse");
		
		/**********************************************************************/		
		ArrayList<Sequence> seqs = SequenceReader.readAll(inputFile, alphabet);		
		SequenceOutputStream ps = SequenceOutputStream
				.makeOutputStream(outputFile);		

		for (int i = 0; i < args.length; i++) {
			String [] toks = args[i].split(":");
			String chr = toks[0];
			toks = toks[1].split("-");
			int start = Integer.parseInt(toks[0]);
			int end = Integer.parseInt(toks[1]);
			
			Sequence seq = findSequence(seqs, chr);
			if (seq == null){
				Logging.error("Sequence " + chr + " not found");
			}else{
				Sequence newSequence = seq.subSequence(start - 1, end);			
				if (rev)
					newSequence = Alphabet.DNA.complement(newSequence);
				
				newSequence.setName(chr+"_"+start+"_"+end);				
				if (format.startsWith("fa"))
					newSequence.writeFasta(ps);
				else
					newSequence.print(ps);
			}
		}// for
		ps.close();
	}
	/**
	 * Find a sequence with ID from a list
	 * @param seqHash
	 * @param id
	 * @return
	 */
	static Sequence findSequence(ArrayList<Sequence> seqs, String id){		
		for (Sequence seq:seqs){
			if (id.equals(seq.getName()))
				return seq;
		}
		return null;
	}
}
