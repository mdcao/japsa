package japsadev.bio.hts.barcode;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import jaligner.matrix.MatrixLoaderException;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
public class BarCode {
	static final int SCAN_WINDOW=120; 
	HashMap<String, SampleData> samplesMap;
	
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
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length() < 500)
				continue;
			//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used
			jaligner.Sequence 	t5 = new jaligner.Sequence("t5_"+seq.getName(), seq.subSequence(0, SCAN_WINDOW).toString()),
								t3 = new jaligner.Sequence("t3_"+seq.getName(), seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()).toString()),
								c5 = new jaligner.Sequence("c5_"+seq.getName(), (Alphabet.DNA.complement(seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()))).toString()),
								c3 = new jaligner.Sequence("c3_"+seq.getName(), (Alphabet.DNA.complement(seq.subSequence(0, SCAN_WINDOW))).toString());
			String[] samples = new String[pop];
			float[] tf = new float[pop],
					tr = new float[pop],
					cr = new float[pop],
					cf = new float[pop];
			int count=0;
			for(String id:samplesMap.keySet()){
				SampleData sample = samplesMap.get(id);
				jaligner.Sequence 	fBarcode = new jaligner.Sequence("F_"+id, sample.getFBarcode().toString()),
									rBarcode = new jaligner.Sequence("R_"+id, sample.getRBarcode().toString());
				jaligner.Alignment 	alignmentsTF = jaligner.SmithWatermanGotoh.align(t5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f),
									alignmentsTR = jaligner.SmithWatermanGotoh.align(t3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f),
									alignmentsCF = jaligner.SmithWatermanGotoh.align(c5, fBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f),
									alignmentsCR = jaligner.SmithWatermanGotoh.align(c3, rBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				samples[count]=id;
				tf[count]=alignmentsTF.getScore();
				tr[count]=alignmentsTR.getScore();
				cf[count]=alignmentsCF.getScore();
				cr[count++]=alignmentsCR.getScore();
				
			}
			Integer[] 	tfRank = new Integer[pop];
			for(int i=0;i<pop;i++)
				tfRank[i]=i;
			Integer[] 	trRank = tfRank.clone(),
						cfRank = tfRank.clone(),
						crRank = tfRank.clone(),
						tRank = tfRank.clone(),
						cRank = tfRank.clone();
			//sort the alignment scores between template sequence and all forward barcode
			Arrays.sort(tfRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tf[o1], tf[o2]);
				}			
			});
			//sort the alignment scores between template sequence and all reverse barcode
			Arrays.sort(trRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tr[o1], tf[o2]);
				}			
			});
			//sort the alignment scores between complement sequence and all forward barcode
			Arrays.sort(cfRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cf[o1], tf[o2]);
				}			
			});
			//sort the alignment scores between complement sequence and all reverse barcode
			Arrays.sort(crRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cr[o1], cr[o2]);
				}			
			});
			//sort the sum of alignment scores between template sequence and all barcode pairs
			Arrays.sort(tRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(tf[o1]+tr[o1], tf[o2]+tr[o1]);
				}			
			});
			//sort the sum of alignment scores between complement sequence and all barcode pairs
			Arrays.sort(cRank, new Comparator<Integer>() {
				@Override 
				public int compare(Integer o1, Integer o2){
					return Float.compare(cf[o1]+cr[o1], cf[o2]+cr[o1]);
				}			
			});
			System.out.println(pop);
			for(int i=0;i<pop;i++)
				System.out.printf("%d - %d -%d -%d -%d -%d\n",tfRank[i], trRank[i], tRank[i], cfRank[i], crRank[i], cRank[i]);
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			if(tf[tRank[0]]+tr[tRank[0]] > cf[cRank[0]]+cr[cRank[0]]){
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[tfRank[0]].equals(samples[trRank[0]])){
					System.out.println("Template sequence " + seq.getName() + " 100% belongs to sample " + samples[tfRank[0]]);
					//do smt
				} else{
					System.out.print("Template sequence " + seq.getName() + " might belongs to sample " + samples[tRank[0]]);
					System.out.println(": tfRank=" + indexOf(tfRank,tRank[0]) + " trRank=" + indexOf(trRank,tRank[0]));
					//do smt
				}
			} else{
				//if both ends of the same sequence report the best alignment with the barcodes
				if(samples[cfRank[0]].equals(samples[crRank[0]])){
					System.out.println("Complement sequence " + seq.getName() + " 100% belongs to sample " + samples[cfRank[0]]);
					//do smt
				} else{
					System.out.print("Complement sequence " + seq.getName() + " might belongs to sample " + samples[cRank[0]]);
					System.out.println(": cfRank=" + indexOf(cfRank,cRank[0]) + " crRank=" + indexOf(crRank,cRank[0]));
					//do smt
				}
			}
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
	
	public static void main(String[] args){
		BarCode bc;
		try {
			bc = new BarCode("/home/s.hoangnguyen/Workspace/poreFUME/inputData/pb_39.fasta");
			bc.clustering("/home/s.hoangnguyen/Workspace/poreFUME/inputData/n.fasta.protein.homolog.fasta");
		} catch (IOException | MatrixLoaderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
	}
	
}
