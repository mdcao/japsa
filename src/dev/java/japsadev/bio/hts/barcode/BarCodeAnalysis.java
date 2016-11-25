package japsadev.bio.hts.barcode;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;
import japsadev.bio.BarcodeAlignment;
public class BarCodeAnalysis {
	static final int SCAN_WINDOW=60, SCORE_THRES=58; 
	ArrayList<Sequence> barCodes;
	Process[] processes;
	int nSamples;
	SequenceOutputStream[] streamToScaffolder, streamToFile;
	
	public BarCodeAnalysis(String barcodeFile, String scriptFile) throws IOException{
			barCodes = SequenceReader.readAll(barcodeFile, Alphabet.DNA());
			nSamples = barCodes.size();
			processes = new Process[nSamples];
			streamToScaffolder = new SequenceOutputStream[nSamples];
			//streamToFile = new SequenceOutputStream[nSamples];
			
			String id;
			for(int i=0;i<nSamples;i++){
				id = barCodes.get(i).getName();
				System.out.println(i + " >" + id + ":" + barCodes.get(i));
				
				ProcessBuilder pb = new ProcessBuilder(scriptFile, id);
				processes[i]  = pb.start();
				Logging.info("Job for " + id  + " started");
				streamToScaffolder[i] = new SequenceOutputStream(processes[i].getOutputStream());
				//streamToFile[i] = SequenceOutputStream.makeOutputStream(id+"_clustered.fasta");
			}
	}
	/*
	 * Trying to clustering MinION read data into different samples based on the barcode
	 */
	public void clustering(String dataFile) throws IOException, InterruptedException{
		SequenceReader reader;
		if(dataFile.equals("-"))
			reader = SequenceReader.getReader(System.in);
		else
			reader = SequenceReader.getReader(dataFile);
		Sequence seq;

		Sequence t5, t3, c5, c3;
		final double[] 	tf = new double[nSamples],
						tr = new double[nSamples],
						cr = new double[nSamples],
						cf = new double[nSamples];
		Integer[] 	tRank = new Integer[nSamples],
					cRank = new Integer[nSamples];
//		jaligner.Alignment[] 	alignmentsTF = new jaligner.Alignment[pop],
//								alignmentsTR = new jaligner.Alignment[pop],
//								alignmentsCF = new jaligner.Alignment[pop],
//								alignmentsCR = new jaligner.Alignment[pop];
				
		Sequence barcodeSeq = new Sequence(Alphabet.DNA4(),21,"barcode");
		Sequence tipSeq = new Sequence(Alphabet.DNA4(),SCAN_WINDOW,"tip");

		BarcodeAlignment barcodeAlignment = new BarcodeAlignment(barcodeSeq, tipSeq);

		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length() < 200){
				System.err.println("Ignore short sequence " + seq.getName());
				continue;
			}
			//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used

			t5 = seq.subSequence(0, SCAN_WINDOW);
			t3 = seq.subSequence(seq.length()-SCAN_WINDOW,seq.length());
			c5 = Alphabet.DNA.complement(seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()));
			c3 = Alphabet.DNA.complement(seq.subSequence(0, SCAN_WINDOW));

	
			for(int i=0;i<nSamples;i++){
				Sequence 	barcode = barCodes.get(i);

				barcodeAlignment.setBarcodeSequence(barcode);
				
				barcodeAlignment.setReadSequence(t5);
				tf[i]=barcodeAlignment.align();
				barcodeAlignment.setReadSequence(c5);
				cf[i]=barcodeAlignment.align();
						
				barcodeAlignment.setReadSequence(t3);
				tr[i]=barcodeAlignment.align();
				barcodeAlignment.setReadSequence(c3);
				cr[i]=barcodeAlignment.align();

			}
			for(int i=0;i<nSamples;i++)
				tRank[i]=cRank[i]=i;

			//sort the sum of alignment scores between template sequence and all barcode pairs
			Arrays.sort(tRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(tf[o1]+tr[o1], tf[o2]+tr[o2]);
				}			
			}));
			//sort the sum of alignment scores between complement sequence and all barcode pairs
			Arrays.sort(cRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(cf[o1]+cr[o1], cf[o2]+cr[o2]);
				}			
			}));
			
			int index=-1;
			if(Math.max(tf[tRank[0]]+tr[tRank[0]], cf[cRank[0]]+cr[cRank[0]]) <= SCORE_THRES){
				//Logging.info("Unknown sequence " + seq.getName());
				continue;
			}
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			else if(tf[tRank[0]]+tr[tRank[0]] > cf[cRank[0]]+cr[cRank[0]]){
				index = tRank[0];
				Logging.info("Template sequence " + seq.getName() + " might belongs to sample " + barCodes.get(index).getName());

			} else{
				index = cRank[0];
				Logging.info("Complement sequence " + seq.getName() + " might belongs to sample " + barCodes.get(index).getName());

			}
			if(index>=0 && index<nSamples){
				if(processes[index]!=null && processes[index].isAlive())
					seq.writeFasta(streamToScaffolder[index]);
				//seq.writeFasta(streamToFile[index]);
			}
		}
		
		System.out.println("Done all input");
		for (int i = 0; i < nSamples;i++){
			if(processes[i]!=null && processes[i].isAlive()){
				streamToScaffolder[i].close();
				//streamToFile[i].close();
				processes[i].waitFor();
			}
		}
		System.out.println("Done every thing");
		reader.close();
	}

//	//display jaligner.Alignment. TODO: convert to ours
//	public void printAlignment(jaligner.Alignment alignment){
//		String 	origSeq1 = alignment.getOriginalSequence1().getSequence(),
//				origSeq2 = alignment.getOriginalSequence2().getSequence(),
//				alnSeq1 = new String(alignment.getSequence1()),
//				alnSeq2 = new String(alignment.getSequence2());
//		int 	start1 = alignment.getStart1(),
//				start2 = alignment.getStart2(),
//				gap1 = alignment.getGaps1(),
//				gap2 = alignment.getGaps2();
//		
//		String seq1, seq2, mark;
//		if(start1>=start2){
//			seq1=origSeq1.substring(0, start1) + alnSeq1 + origSeq1.substring(start1+alnSeq1.length()-gap1);
//			String 	seq2Filler = start1==start2?"":String.format("%"+(start1-start2)+"s", ""),
//					markFiller = start1==0?"":String.format("%"+start1+"s", "");
//			seq2= seq2Filler + origSeq2.substring(0, start2) + alnSeq2 + origSeq2.substring(start2+alnSeq2.length()-gap2);
//			mark= markFiller+String.valueOf(alignment.getMarkupLine());
//		}else{
//			seq2=origSeq2.substring(0, start2) + alnSeq2 + origSeq2.substring(start2+alnSeq2.length()-gap2);
//			String 	markFiller = start2==0?"":String.format("%"+start2+"s", "");
//			seq1=String.format("%"+(start2-start1)+"s", "") + origSeq1.substring(0, start1) + alnSeq1 + origSeq1.substring(start1+alnSeq1.length()-gap1);
//			mark=markFiller+String.valueOf(alignment.getMarkupLine());
//		}
//		//System.out.println(alignment.getSummary());
//		System.out.println(seq1);
//		System.out.println(mark);
//		System.out.println(seq2);
//	}

}
