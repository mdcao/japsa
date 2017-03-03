package japsadev.bio.np.phage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.bio.np.ErrorCorrection;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;



public class ConsensusGenerator {
	HashMap<String,Sequence> sequences = new HashMap<String,Sequence>();
	public ConsensusGenerator(String sequenceFile) throws IOException{
		SequenceReader reader = SequenceReader.getReader(sequenceFile);
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
				sequences.put(seq.getName(), seq);
		}
		reader.close();
	}
	public Sequence getSeq(String name){
		return sequences.get(name);
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		ConsensusGenerator cg = new ConsensusGenerator("/home/s.hoangnguyen/Projects/Phage/delta/insert_1.fasta");
		
		BufferedReader groupReader = new BufferedReader(new FileReader("/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_nnp1"));
		SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream("/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_trimmed_nnp1.consensus");
		String s;
		int count=1;
		
		while((s=groupReader.readLine()) != null){
			ArrayList<Sequence> aGroup = new ArrayList<Sequence>();
			String[] toks = s.split(" ");
			for(int i=0; i < toks.length; i++){
				aGroup.add(cg.getSeq(toks[i]));
			}
			Sequence consensus = ErrorCorrection.consensusSequence(aGroup, "insert_nnp1", "kalign");

			//trimming the flanks
			consensus=consensus.subSequence(SequenceExtractor.FLANKING,consensus.length()-SequenceExtractor.FLANKING);
			
			consensus.setName("N1_"+count);
			consensus.writeFasta(outFile);
			count++;
			
		}
		groupReader.close();
		outFile.close();
	}

}
