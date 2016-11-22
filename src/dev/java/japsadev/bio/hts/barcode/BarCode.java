package japsadev.bio.hts.barcode;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import jaligner.matrix.MatrixLoaderException;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
public class BarCode {
	static final int SCAN_WINDOW=60; 
	HashMap<String, SampleData> samplesMap;
	
	String[] stt = {"First","Second","Third"};
	
	public BarCode(String barcodeFile) throws IOException{
		samplesMap = new HashMap<String, SampleData>();
		SequenceReader reader = new FastaReader(barcodeFile);
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			//header must follow the format [FR]_id
			String[] header = seq.getName().split("_");
			SampleData sample = samplesMap.get(header[1]);
			if(sample==null){
				sample = new SampleData();
				samplesMap.put(header[1], sample);
			}
			sample.setId(header[1]);
			
			if(header[0].equals("F"))
				sample.setFBarcode(seq);
			else if (header[0].equals("R"))
				sample.setRBarcode(seq);
			//TODO: how to set corresponding SPAdes contigs file???
			
		}
		reader.close();
	}
	/*
	 * Trying to clustering MinION read data into different samples based on the barcode
	 */
	public void clustering(String dataFile) throws IOException, MatrixLoaderException{
		int pop = samplesMap.size();
		SequenceReader reader = new FastaReader(dataFile);
		Sequence seq;
		jaligner.Sequence t5, t3, c5, c3;
		String[] samples = new String[pop];
		final float[] tf = new float[pop],
				tr = new float[pop],
				cr = new float[pop],
				cf = new float[pop];
		jaligner.Alignment[] 	alignmentsTF = new jaligner.Alignment[pop],
								alignmentsTR = new jaligner.Alignment[pop],
								alignmentsCF = new jaligner.Alignment[pop],
								alignmentsCR = new jaligner.Alignment[pop];
								
		Integer[] 	tfRank = new Integer[pop],
					trRank = new Integer[pop],
					cfRank = new Integer[pop],
					crRank = new Integer[pop],
					tRank = new Integer[pop],
					cRank = new Integer[pop];
		int readNum=0;								
		
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length() < 500){
				System.out.println("Ignore short sequence " + seq.getName());
				continue;
			}
			//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used
			t5 = new jaligner.Sequence("t5_"+seq.getName(), seq.subSequence(0, SCAN_WINDOW).toString());
			t3 = new jaligner.Sequence("t3_"+seq.getName(), seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()).toString());
			c5 = new jaligner.Sequence("c5_"+seq.getName(), (Alphabet.DNA.complement(seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()))).toString());
			c3 = new jaligner.Sequence("c3_"+seq.getName(), (Alphabet.DNA.complement(seq.subSequence(0, SCAN_WINDOW))).toString());

			int count=0;
			for(String id:samplesMap.keySet()){
				SampleData sample = samplesMap.get(id);
				jaligner.Sequence 	fBarcode = new jaligner.Sequence("F_"+id, sample.getFBarcode().toString()),
									rBarcode = new jaligner.Sequence("R_"+id, sample.getRBarcode().toString());
				alignmentsTF[count] = jaligner.SmithWatermanGotoh.align(t5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentsTR[count] = jaligner.SmithWatermanGotoh.align(t3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentsCF[count] = jaligner.SmithWatermanGotoh.align(c5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentsCR[count] = jaligner.SmithWatermanGotoh.align(c3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				
				samples[count]=id;
				tf[count]=alignmentsTF[count].getScore();
				tr[count]=alignmentsTR[count].getScore();
				cf[count]=alignmentsCF[count].getScore();
				cr[count]=alignmentsCR[count].getScore();
				count++;
			}
			for(int i=0;i<pop;i++)
				tfRank[i]=trRank[i]=cfRank[i]=crRank[i]=tRank[i]=cRank[i]=i;
			//sort the alignment scores between template sequence and all forward barcode
			Arrays.sort(tfRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tf[o1], tf[o2]);
				}			
			}));
			//sort the alignment scores between template sequence and all reverse barcode
			Arrays.sort(trRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tr[o1], tr[o2]);
				}			
			}));
			//sort the alignment scores between complement sequence and all forward barcode
			Arrays.sort(cfRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cf[o1], cf[o2]);
				}			
			}));
			//sort the alignment scores between complement sequence and all reverse barcode
			Arrays.sort(crRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cr[o1], cr[o2]);
				}			
			}));
			//sort the sum of alignment scores between template sequence and all barcode pairs
			Arrays.sort(tRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tf[o1]+tr[o1], tf[o2]+tr[o2]);
				}			
			}));
			//sort the sum of alignment scores between complement sequence and all barcode pairs
			Arrays.sort(cRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cf[o1]+cr[o1], cf[o2]+cr[o2]);
				}			
			}));
			System.out.println("=============================================================================================");
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			if(tf[tRank[0]]+tr[tRank[0]] > cf[cRank[0]]+cr[cRank[0]]){
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[tfRank[0]].equals(samples[trRank[0]])){
					System.out.println("\nTemplate sequence " + seq.getName() + " 100% belongs to sample " + samples[tfRank[0]]);
					//do smt

				} else{
					System.out.print("\nTemplate sequence " + seq.getName() + " might belongs to sample " + samples[tRank[0]]);
					System.out.println(": tfRank=" + indexOf(tfRank,tRank[0]) + " trRank=" + indexOf(trRank,tRank[0]));
					//do smt
				}
				System.out.println();
				for(int i=0;i<3;i++){
					System.out.println(stt[i] + " best:");
					printAlignment(alignmentsTF[tRank[i]]);
					printAlignment(alignmentsTR[tRank[i]]);
				}
				System.out.println();

			} else{
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[cfRank[0]].equals(samples[crRank[0]])){
					System.out.println("\nComplement sequence " + seq.getName() + " 100% belongs to sample " + samples[cfRank[0]]);
					//do smt

				} else{
					System.out.print("\nComplement sequence " + seq.getName() + " might belongs to sample " + samples[cRank[0]]);
					System.out.println(": cfRank=" + indexOf(cfRank,cRank[0]) + " crRank=" + indexOf(crRank,cRank[0]));
					//do smt
				}
				System.out.println();
				for(int i=0;i<3;i++){
					System.out.println(stt[i] + " best:");
					printAlignment(alignmentsCF[cRank[i]]);
					printAlignment(alignmentsCR[cRank[i]]);
				}
				System.out.println();
			}
			//invoke garbage collection explicitly every 50 reads. TODO: implement self-stand Smith-Waterman alignment
			readNum++;
			if(readNum%50==0)
				System.gc();
		}
		reader.close();
	}
	
	//helper
	int indexOf(Integer[] arr, int value){
		int retVal=-1;
		for(int i=0;i<arr.length;i++)
			if(value==arr[i])
				return i;
		return retVal;
	}
	
	public void printAlignment(jaligner.Alignment alignment){
		String 	origSeq1 = alignment.getOriginalSequence1().getSequence(),
				origSeq2 = alignment.getOriginalSequence2().getSequence(),
				alnSeq1 = new String(alignment.getSequence1()),
				alnSeq2 = new String(alignment.getSequence2());
		int 	start1 = alignment.getStart1(),
				start2 = alignment.getStart2(),
				gap1 = alignment.getGaps1(),
				gap2 = alignment.getGaps2();
		
		String seq1, seq2, mark;
		if(start1>=start2){
			seq1=origSeq1.substring(0, start1) + alnSeq1 + origSeq1.substring(start1+alnSeq1.length()-gap1);
			String 	seq2Filler = start1==start2?"":String.format("%"+(start1-start2)+"s", ""),
					markFiller = start1==0?"":String.format("%"+start1+"s", "");
			seq2= seq2Filler + origSeq2.substring(0, start2) + alnSeq2 + origSeq2.substring(start2+alnSeq2.length()-gap2);
			mark= markFiller+String.valueOf(alignment.getMarkupLine());
		}else{
			seq2=origSeq2.substring(0, start2) + alnSeq2 + origSeq2.substring(start2+alnSeq2.length()-gap2);
			String 	markFiller = start2==0?"":String.format("%"+start2+"s", "");
			seq1=String.format("%"+(start2-start1)+"s", "") + origSeq1.substring(0, start1) + alnSeq1 + origSeq1.substring(start1+alnSeq1.length()-gap1);
			mark=markFiller+String.valueOf(alignment.getMarkupLine());
		}
		//System.out.println(alignment.getSummary());
		System.out.println(seq1);
		System.out.println(mark);
		System.out.println(seq2);
	}
	
}
