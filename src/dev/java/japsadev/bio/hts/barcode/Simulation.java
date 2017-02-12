package japsadev.bio.hts.barcode;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

public class Simulation {
	public static void main(String[] args) throws IOException{
//		ArrayList<Sequence> barCodes = SequenceReader.readAll("/home/s.hoangnguyen/Workspace/barcode/my_barcode.fasta", Alphabet.DNA());
//		int nSamples = barCodes.size();
//		SequenceReader reader = SequenceReader.getReader("/home/s.hoangnguyen/Workspace/barcode/LejlaControl.2D.min500bp.fasta");
//		SequenceOutputStream out = SequenceOutputStream.makeOutputStream("/home/s.hoangnguyen/Workspace/barcode/test.fasta");
//		Sequence seq;
//		Random random = new Random();
//		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
//			int index = random.nextInt(nSamples);
//			Sequence 	barcode = barCodes.get(index),
//						barcode_rc = Alphabet.DNA.complement(barcode);
//			SequenceBuilder seqBuild = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  barcode.getName()+ "_" +seq.getName());
//			boolean reverse = random.nextBoolean();
//			seqBuild.append(reverse?barcode:barcode_rc);
//			seqBuild.append(seq);
//			seqBuild.append(reverse?barcode:barcode_rc);
//			
//			seqBuild.writeFasta(out);
//		}
//		out.close();
		String currentDirectory;
		currentDirectory = System.getProperty("user.dir");
		System.out.println("Current working directory : "+currentDirectory);
	}
}
