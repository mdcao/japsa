package japsa.bio.np.barcode;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import jaligner.matrix.MatrixLoaderException;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;
public class BarCodeAnalysisVerbose {
	static final int SCAN_WINDOW=120, SCORE_THRES=70; 
	public static boolean toPrint=false;
	ArrayList<Sequence> barCodes;
	ArrayList<Sequence> barCodeComps;
	Process[] processes;
	int nSamples;
	SequenceOutputStream[] streamToScaffolder, streamToFile;

	public BarCodeAnalysisVerbose(String barcodeFile, String scriptFile) throws IOException{
		barCodes = SequenceReader.readAll(barcodeFile, Alphabet.DNA());
		nSamples = barCodes.size();



		processes = new Process[nSamples];
		streamToScaffolder = new SequenceOutputStream[nSamples];
		if(toPrint)
			streamToFile = new SequenceOutputStream[nSamples];

		barCodeComps = new ArrayList<Sequence> (barCodes.size());
		String id;
		for(int i=0;i<nSamples;i++){		
			Sequence barCode = barCodes.get(i);
			barCodeComps.add(Alphabet.DNA.complement(barCode));

			id = barCode.getName();
			ProcessBuilder pb = new ProcessBuilder(scriptFile, id)
					.redirectError(new File("log_" + id + ".err"))
					.redirectOutput(new File("log_" + id + ".out"));
			pb.directory(new File(System.getProperty("user.dir")));

			processes[i]  = pb.start();

			Logging.info("Job for " + id  + " started");
			streamToScaffolder[i] = new SequenceOutputStream(processes[i].getOutputStream());
			if(toPrint)
				streamToFile[i] = SequenceOutputStream.makeOutputStream(id+"_clustered.fasta");
		}
	}
	/*
	 * Trying to clustering MinION read data into different samples based on the barcode
	 */
	public void clustering(String dataFile) throws IOException, InterruptedException, MatrixLoaderException{
		SequenceReader reader;
		if(dataFile.equals("-"))
			reader = SequenceReader.getReader(System.in);
		else
			reader = SequenceReader.getReader(dataFile);
		Sequence seq;

		Sequence s5, s3;
		final double[] 	tf = new double[nSamples],
						tr = new double[nSamples],
						cr = new double[nSamples],
						cf = new double[nSamples];

		jaligner.Alignment 		alignmentTF = new jaligner.Alignment(),
								alignmentTR = new jaligner.Alignment(),
								alignmentCF = new jaligner.Alignment(),
								alignmentCR = new jaligner.Alignment();
		jaligner.Sequence js5,js3, jBarcode, jBarcodeComp;
		
//		float[][] scores = {
//				{  2.7f, -4.5f, -4.5f, -4.5f},
//				{ -4.5f,  2.7f, -4.5f, -4.5f},
//				{ -4.5f, -4.5f,  2.7f, -4.5f},
//				{ -4.5f, -4.5f, -4.5f,  2.7f}			
//		};
//		jaligner.matrix.Matrix 	matrix = new jaligner.matrix.Matrix("poreFUME",scores);
		
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length() < 200){
				System.out.println("Ignore short sequence " + seq.getName());
				continue;
			}
			//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used

			s5 = seq.subSequence(0, SCAN_WINDOW);
			s3 = seq.subSequence(seq.length()-SCAN_WINDOW,seq.length());


			double bestScore = 0.0;
			int bestIndex = nSamples;
			jaligner.Alignment bestAlignment = new jaligner.Alignment();

			for(int i=0;i<nSamples; i++){
				Sequence barcode = barCodes.get(i);
				Sequence barcodeComp = barCodeComps.get(i);

//				barcodeAlignment.setBarcodeSequence(barcode);				
//				barcodeAlignment.setReadSequence(s5);				
//				tf[i]=barcodeAlignment.align();				
//
//				barcodeAlignment.setReadSequence(s3);
//				tr[i]=barcodeAlignment.align();
//
//				barcodeAlignment.setBarcodeSequence(barcodeComp);
//				barcodeAlignment.setReadSequence(s3);
//				cr[i]=barcodeAlignment.align();
//				barcodeAlignment.setReadSequence(s5);
//				cf[i]=barcodeAlignment.align();
				js5 = new jaligner.Sequence(s5.toString());
				js3 = new jaligner.Sequence(s3.toString());
				jBarcode = new jaligner.Sequence(barcode.toString());
				jBarcodeComp = new jaligner.Sequence(barcodeComp.toString());

				alignmentTF = jaligner.SmithWatermanGotoh.align(js5, jBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentTR = jaligner.SmithWatermanGotoh.align(js3, jBarcode, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentCF = jaligner.SmithWatermanGotoh.align(js5, jBarcodeComp, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				alignmentCR = jaligner.SmithWatermanGotoh.align(js3, jBarcodeComp, jaligner.matrix.MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
				
				tf[i] = alignmentTF.getScore();
				tr[i] = alignmentTR.getScore();
				cf[i] = alignmentCF.getScore();
				cr[i] = alignmentCR.getScore();

				//This is for both end
				//double myScore = Math.max(tf[i], tr[i]) + Math.max(cf[i], cr[i]);

				//but may be we can use onely one end
				double myScore = Math.max(Math.max(tf[i], tr[i]), Math.max(cf[i], cr[i]));
				if (myScore > bestScore){
					//Logging.info("Better score=" + myScore);
					bestScore = myScore;
					bestIndex = i;
					if(Math.max(tf[i], tr[i]) > Math.max(cf[i], cr[i])){
						if(tf[i]>tr[i])
							bestAlignment = alignmentTF;
						else 
							bestAlignment = alignmentTR;
					}else{
						if(cf[i]>cr[i])
							bestAlignment = alignmentCF;
						else 
							bestAlignment = alignmentCR;
					}
				}
			}

			if(bestScore <= SCORE_THRES){
				//Logging.info("Unknown sequence " + seq.getName());
				continue;
			}
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			else {
				System.out.println("Sequence " + seq.getName() + " might belongs to sample " + barCodes.get(bestIndex).getName() + " with score=" + bestScore);
				printAlignment(bestAlignment);
				if(bestIndex<nSamples && processes[bestIndex]!=null && processes[bestIndex].isAlive()){
					if(bestIndex==0){
						System.out.println("...skip to write to stream " + barCodes.get(0).getName());
					}else{
						System.out.println("...writing to stream " + bestIndex);
						seq.writeFasta(streamToScaffolder[bestIndex]);
						if(toPrint)
							seq.writeFasta(streamToFile[bestIndex]);
					}
				}
			}

		}

		System.out.println("Done all input");
		for (int i = 0; i < nSamples;i++){
			if(processes[i]!=null && processes[i].isAlive()){
				streamToScaffolder[i].close();
				if(toPrint)
					streamToFile[i].close();
				processes[i].waitFor();
			}
		}
		System.out.println("Done every thing");
		reader.close();
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

