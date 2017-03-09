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
	HashMap<String,Sequence> sequences;
	public ConsensusGenerator(){
		sequences = new HashMap<String,Sequence>();
	}
	
	public Sequence getSeq(String name){
		return sequences.get(name);
	}
	
	public void generate(String seqFile, String listFile, String outFile, String aligner, boolean trim) throws IOException, InterruptedException{
		
		SequenceReader reader = SequenceReader.getReader(seqFile);
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
				sequences.put(seq.getName(), seq);
		}
		reader.close();
		
		BufferedReader groupReader = new BufferedReader(new FileReader(listFile));
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(outFile);
		String s;
		int count=1;
		
		while((s=groupReader.readLine()) != null){
			ArrayList<Sequence> aGroup = new ArrayList<Sequence>();
			String[] toks = s.split(" ");
			for(int i=0; i < toks.length; i++){
				aGroup.add(getSeq(toks[i]));
			}
			Sequence consensus = ErrorCorrection.consensusSequence(aGroup, "grouping", aligner);

			//trimming the flanks
			if(trim)
				consensus=consensus.subSequence(SequenceExtractor.FLANKING,consensus.length()-SequenceExtractor.FLANKING);
			
				
			consensus.setName("S_"+count);
			consensus.writeFasta(out);
			count++;
			
		}
		groupReader.close();
		out.close();
		
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
//		ConsensusGenerator cg = new ConsensusGenerator("/home/s.hoangnguyen/Projects/Phage/delta/insert_1.fasta");
//		
//		BufferedReader groupReader = new BufferedReader(new FileReader("/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_nnp1"));
//		SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream("/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_trimmed_nnp1.consensus");
		

		ConsensusGenerator cg = new ConsensusGenerator();
//		cg.generate("/home/s.hoangnguyen/Projects/Phage/delta/sangerInserts.fasta", 
//					"/home/s.hoangnguyen/Projects/Phage/delta/blastclust/sanger-blastclust", 
//					"/home/s.hoangnguyen/Projects/Phage/delta/blastclust/sanger.consensus",
//					"kalign",false);
		cg.generate("/home/s.hoangnguyen/Projects/Phage/delta/insert_2.fasta", 
				"/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_nnp2_88", 
				"/home/s.hoangnguyen/Projects/Phage/delta/blastclust/insert_nnp2_88.consensus",
				"kalign",false);
	}

}
