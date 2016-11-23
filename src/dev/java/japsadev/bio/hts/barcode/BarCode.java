package japsadev.bio.hts.barcode;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.Logging;
import japsadev.bio.BarcodeAlignment;
public class BarCode {
	static final int SCAN_WINDOW=60; 
	HashMap<String, SampleData> samplesMap;
	
	public BarCode(String barcodeFile) throws IOException{
		samplesMap = new HashMap<String, SampleData>();
		SequenceReader reader = new FastaReader(barcodeFile);
		Sequence seq;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			//header must follow the format [F/R]_id
			String 	ori = seq.getName().substring(0, 1),
					id = seq.getName().substring(2);
			SampleData sample = samplesMap.get(id);
			if(sample==null){
				sample = new SampleData(id);
				samplesMap.put(id, sample);
			}

			if(ori.equals("F"))
				sample.setFBarcode(seq);
			else if (ori.equals("R"))
				sample.setRBarcode(seq);
			//TODO: how to set corresponding SPAdes contigs file???

		}
		reader.close();
	}
	/*
	 * Trying to clustering MinION read data into different samples based on the barcode
	 */
	public void clustering(String dataFile) throws IOException{
		int pop = samplesMap.size();
		SequenceReader reader;
		if(dataFile.equals("-"))
			reader = new FastaReader(System.in);
		else
			reader = new FastaReader(dataFile);
		Sequence seq;

		Sequence t5, t3, c5, c3;
		String[] samples = new String[pop];
		final double[] 	tf = new double[pop],
						tr = new double[pop],
						cr = new double[pop],
						cf = new double[pop];
//		jaligner.Alignment[] 	alignmentsTF = new jaligner.Alignment[pop],
//								alignmentsTR = new jaligner.Alignment[pop],
//								alignmentsCF = new jaligner.Alignment[pop],
//								alignmentsCR = new jaligner.Alignment[pop];
								
		Integer[] 	tfRank = new Integer[pop],
					trRank = new Integer[pop],
					cfRank = new Integer[pop],
					crRank = new Integer[pop],
					tRank = new Integer[pop],
					cRank = new Integer[pop];
//		int readNum=0;								


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

			int count=0;
		
			for(String id:samplesMap.keySet()){
				SampleData sample = samplesMap.get(id);
				Sequence 	fBarcode = sample.getFBarcode(),
							rBarcode = sample.getRBarcode();

//				alignmentsTF[count] = jaligner.SmithWatermanGotoh.align(t5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
//				alignmentsTR[count] = jaligner.SmithWatermanGotoh.align(t3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
//				alignmentsCF[count] = jaligner.SmithWatermanGotoh.align(c5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
//				alignmentsCR[count] = jaligner.SmithWatermanGotoh.align(c3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);

				samples[count]=id;
				//TODO: copy the content, not object
				barcodeAlignment.setBarcodeSequence(fBarcode);
				
				barcodeAlignment.setReadSequence(t5);
				tf[count]=barcodeAlignment.align();
				barcodeAlignment.setReadSequence(c5);
				cf[count]=barcodeAlignment.align();
				
				barcodeAlignment.setBarcodeSequence(rBarcode);
				
				barcodeAlignment.setReadSequence(t3);
				tr[count]=barcodeAlignment.align();
				barcodeAlignment.setReadSequence(c3);
				cr[count]=barcodeAlignment.align();
				count++;


			}
			for(int i=0;i<pop;i++)
				tfRank[i]=trRank[i]=cfRank[i]=crRank[i]=tRank[i]=cRank[i]=i;

			//sort the alignment scores between template sequence and all forward barcode
			Arrays.sort(tfRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(tf[o1], tf[o2]);
				}			
			}));
			//sort the alignment scores between template sequence and all reverse barcode
			Arrays.sort(trRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(tr[o1], tr[o2]);
				}			
			}));
			//sort the alignment scores between complement sequence and all forward barcode
			Arrays.sort(cfRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(cf[o1], cf[o2]);
				}			
			}));
			//sort the alignment scores between complement sequence and all reverse barcode
			Arrays.sort(crRank, Collections.reverseOrder(new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Double.compare(cr[o1], cr[o2]);
				}			
			}));
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
			
			if(Math.max(tf[tRank[0]]+tr[tRank[0]], cf[cRank[0]]+cr[cRank[0]]) <= 58){
				Logging.info("Unknown sequence " + seq.getName());
				continue;
			}
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			else if(tf[tRank[0]]+tr[tRank[0]] > cf[cRank[0]]+cr[cRank[0]]){
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[tfRank[0]].equals(samples[trRank[0]])){
					Logging.info("Template sequence " + seq.getName() + " 100% belongs to sample " + samples[tfRank[0]]);
					//do smt

				} else{
					Logging.info("Template sequence " + seq.getName() + " might belongs to sample " + samples[tRank[0]]+": tfRank=" + indexOf(tfRank,tRank[0]) + " trRank=" + indexOf(trRank,tRank[0]));
					//do smt
				}
				samplesMap.get(samples[tRank[0]]).passRead(seq);
//				for(int i=0;i<pop;i++)
//					System.out.printf("%dT:%.2f+%.2f=%.2f ", i,tr[tRank[i]], tf[tRank[i]], tr[tRank[i]] + tf[tRank[i]]);
//				System.out.println();


			} else{
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[cfRank[0]].equals(samples[crRank[0]])){
					Logging.info("Complement sequence " + seq.getName() + " 100% belongs to sample " + samples[cfRank[0]]);
					//do smt

				} else{
					Logging.info("Complement sequence " + seq.getName() + " might belongs to sample " + samples[cRank[0]] + ": cfRank=" + indexOf(cfRank,cRank[0]) + " crRank=" + indexOf(crRank,cRank[0]));
					//do smt
				}
				samplesMap.get(samples[cRank[0]]).passRead(seq);
//				for(int i=0;i<pop;i++)
//					System.out.printf("%dC:%.2f+%.2f=%.2f ", i,cr[cRank[i]], cf[cRank[i]], cr[cRank[i]] + cf[cRank[i]]);
//				System.out.println();

			}
		}
		
		for(SampleData sample:samplesMap.values()){
			if(sample.terminate())
				Logging.info("All done successfully!");
			else
				Logging.error("Cannot finish properly!");
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
	
	//display jaligner.Alignment. TODO: convert to ours
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
