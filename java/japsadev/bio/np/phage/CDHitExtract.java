package japsadev.bio.np.phage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.bio.np.ErrorCorrection;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

public class CDHitExtract {

	public static void main(String[] args) throws IOException, InterruptedException {
		// TODO Auto-generated method stub
		ArrayList<Sequence> representSeq = SequenceReader.readAll("/home/sonhoanghguyen/Projects/Phage/paper/poa-consensus/allNanopore.fasta", Alphabet.DNA());
		HashMap<String,Sequence> map = new HashMap<String,Sequence>();
		for(Sequence seq:representSeq){
			map.put(seq.getName(), seq);
		}
		String aligner = "clustal";
		BufferedReader pathReader = new BufferedReader(new FileReader("/home/sonhoanghguyen/Projects/Phage/paper/poa-consensus/nanopore.clstr"));
		SequenceOutputStream out =  SequenceOutputStream.makeOutputStream("/home/sonhoanghguyen/Projects/Phage/paper/poa-consensus/nanopore.fasta");
		String s;
		//Read contigs from contigs.paths and refer themselves to contigs.fasta
		Sequence consensus = new Sequence( Alphabet.DNA(), 10000);
		int count = 0;
		String seq = "";
		ArrayList<Sequence> aGroup=new ArrayList<Sequence>();;
		while((s=pathReader.readLine()) != null){
			if(s.startsWith(">")){
				if(count>0){
					//group=map.get(seq);
					if(count > 1){ 
						System.out.println("Consensusing group with " + aGroup.size() + " members");						
						consensus = ErrorCorrection.consensusSequence(aGroup, "grouping", aligner);
					}
					else
						consensus = map.get(seq);
					consensus.setName(seq); // name of the CDHit representative sequence, but content is the consensus
					consensus.setDesc(aligner+"="+count);
					//System.out.println(group.getName() + " : " + group.getDesc());
					consensus.writeFasta(out);
				}
				aGroup = new ArrayList<Sequence>();
				count=0;
			}else{
				count++;
				aGroup.add(map.get(s.substring(s.indexOf(">")+1, s.indexOf("..."))));
				if(s.contains("*")){
					seq = s.substring(s.indexOf(">")+1, s.indexOf("..."));
					
				}
			}

		}
		//last round
		if(count>0){
			//group=map.get(seq);
			if(count > 1){ 
				System.out.println("Consensusing group with " + aGroup.size() + " members");						
				consensus = ErrorCorrection.consensusSequence(aGroup, "grouping", aligner);
			}
			else
				consensus = map.get(seq);
			consensus.setName(seq); // name of the CDHit representative sequence, but content is the consensus
			consensus.setDesc(aligner+"="+count);
			//System.out.println(group.getName() + " : " + group.getDesc());
			consensus.writeFasta(out);
		}
		
		pathReader.close();
		out.close();
	}

}
