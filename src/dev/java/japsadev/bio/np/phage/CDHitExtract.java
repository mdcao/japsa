package japsadev.bio.np.phage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

public class CDHitExtract {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		ArrayList<Sequence> representSeq = SequenceReader.readAll("/home/sonhoanghguyen/Projects/Phage/paper/sanger/sanger", Alphabet.DNA());
		HashMap<String,Sequence> map = new HashMap<String,Sequence>();
		for(Sequence seq:representSeq){
			map.put(seq.getName(), seq);
		}
		
		BufferedReader pathReader = new BufferedReader(new FileReader("/home/sonhoanghguyen/Projects/Phage/paper/sanger/sanger.clstr"));
		SequenceOutputStream out =  SequenceOutputStream.makeOutputStream("/home/sonhoanghguyen/Projects/Phage/paper/sanger/sanger.fasta");
		String s;
		//Read contigs from contigs.paths and refer themselves to contigs.fasta
		Sequence group = new Sequence( Alphabet.DNA(), 10000);
		int count = 0;
		String seq = "";
		
		while((s=pathReader.readLine()) != null){
			if(s.startsWith(">")){
				if(count>0){
					group=map.get(seq);
					//group.setName(s.substring(1));
					group.setDesc("count="+count);
					//System.out.println(group.getName() + " : " + group.getDesc());
					group.writeFasta(out);
				}
				
				count=0;
			}else{
				count++;
				if(s.contains("*")){
					seq = s.substring(s.indexOf(">")+1, s.indexOf("..."));
					
				}
			}

		}
		pathReader.close();
	}

}
