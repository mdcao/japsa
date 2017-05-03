package japsadev.bio.hts.barcode;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.Logging;

public class BarCodeAnalysis {
	int 	SCAN_WINDOW; 
	double	DIST_THRES,
			SCORE_THRES;
	public static boolean 	print=false, //whether to print to files
							script=false, //whether invoke downstream analysis
							twoends=false; // both-ends-matching 
	ArrayList<Sequence> barCodesLeft = new ArrayList<Sequence>(); //barcode sequences from left end
	ArrayList<Sequence> barCodesRight = new ArrayList<Sequence>(); //barcode from right end
	Process[] processes;
	int nSamples;
	int barcodeLen;
	SequenceOutputStream[] streamToScript, streamToFile;

	public BarCodeAnalysis(String barcodeFile, String scriptFile) throws IOException{
		if(scriptFile!=null)
			script=true;
		
		ArrayList<Sequence> allSeq = SequenceReader.readAll(barcodeFile, Alphabet.DNA());
		if(twoends){
			allSeq.sort(Comparator.comparing(Sequence::getName));
			int i = 0;
			while(i < allSeq.size()){
				barCodesLeft.add(allSeq.get(i++));
				barCodesRight.add(allSeq.get(i++));
			}
		}else{
			barCodesLeft = allSeq;
			for(Sequence seq:barCodesLeft)
				barCodesRight.add(Alphabet.DNA.complement(seq));

		}
		
		nSamples = barCodesLeft.size();

		processes = new Process[nSamples];
		
		if(script)
			streamToScript = new SequenceOutputStream[nSamples];
		if(print)
			streamToFile = new SequenceOutputStream[nSamples+1]; //unknown sequences included

		String id;
		for(int i=0;i<nSamples;i++){		
			Sequence barCode = barCodesLeft.get(i);

			id = barCode.getName();
			//System.out.println(i + " >" + id + ":" + barCode);
			if(script){
				ProcessBuilder pb = new ProcessBuilder(scriptFile, id)
						.redirectError(new File("/dev/null"))
						.redirectOutput(new File("/dev/null"));
				//pb.directory(new File(System.getProperty("user.dir")));
				
				processes[i]  = pb.start();
	
				Logging.info("Job for " + id  + " started");
				streamToScript[i] = new SequenceOutputStream(processes[i].getOutputStream());
			}
			barcodeLen += barCodesLeft.get(i).length();
		}
		
		barcodeLen /= nSamples;
		
		if(print){
			streamToFile = new SequenceOutputStream[nSamples+1]; // plus unknown
			for(int i=0;i<nSamples;i++){		
				streamToFile[i] = SequenceOutputStream.makeOutputStream(barCodesLeft.get(i).getName()+".fastq");
			}
			streamToFile[nSamples] = SequenceOutputStream.makeOutputStream("unknown.fastq");

		}
		
		
		SCAN_WINDOW = barcodeLen * 3;
		SCORE_THRES = .7; //70%
		DIST_THRES = .1; //10%
	}
	public void setScanWindow(int window){
		SCAN_WINDOW = window;
	}
	
	public void setThreshold(double ident){
		SCORE_THRES = ident/100;
	}
	
	public void setDistance(double dist){
		DIST_THRES= dist/100;
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

		Sequence s5, s3;
		final double[] 	lf = new double[nSamples], //left-forward
				lr = new double[nSamples],	//left-reversed
				rr = new double[nSamples], //right-reversed
				rf = new double[nSamples]; //right-forward
		
		SWGAlignment 	alignmentLF = new SWGAlignment(),
						alignmentLR = new SWGAlignment(),
						alignmentRF = new SWGAlignment(),
						alignmentRR = new SWGAlignment();

		SWGAlignment 	bestLeftAlignment = new SWGAlignment(),
						bestRightAlignment = new SWGAlignment();



		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length() < barcodeLen*4){
				System.err.println("Ignore short sequence " + seq.getName());
				continue;
			}
			//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used

			s5 = seq.subSequence(0, SCAN_WINDOW);
			s3 = Alphabet.DNA.complement(seq.subSequence(seq.length()-SCAN_WINDOW,seq.length()));



			double bestScore = 0.0;
			double distance = 0.0; //distance between bestscore and the runner-up
			
			int bestIndex = nSamples;

			for(int i=0;i<nSamples; i++){
				Sequence barcodeLeft = barCodesLeft.get(i);
				Sequence barcodeRight = barCodesRight.get(i); //rc of right barcode sequence

			
				alignmentLF = SWGAlignment.align(s5, barcodeLeft);
				alignmentLR = SWGAlignment.align(s3, barcodeLeft);
				alignmentRF = SWGAlignment.align(s5, barcodeRight);
				alignmentRR = SWGAlignment.align(s3, barcodeRight);
				
				lf[i] = alignmentLF.getIdentity()/(float)Math.max(barcodeLeft.length(),alignmentLF.getLength());
				lr[i] = alignmentLR.getIdentity()/(float)Math.max(barcodeLeft.length(),alignmentLR.getLength());
				rf[i] = alignmentRF.getIdentity()/(float)Math.max(barcodeRight.length(),alignmentRF.getLength());
				rr[i] = alignmentRR.getIdentity()/(float)Math.max(barcodeRight.length(),alignmentRR.getLength());
				

				double myScore = 0.0;
				if(twoends){
					myScore = Math.max(lf[i] + rr[i] , lr[i] + rf[i])/2;
				}
				else{	
					myScore = Math.max(Math.max(lf[i], lr[i]), Math.max(rf[i], rr[i]));
				}
				
				if (myScore > bestScore){
					//Logging.info("Better score=" + myScore);
					distance = myScore-bestScore;
					bestScore = myScore;		
					bestIndex = i;
					if(twoends){
						if(lf[i] + rr[i] > lr[i] + rf[i]){
							bestLeftAlignment = alignmentLF;
							bestRightAlignment = alignmentRR;
						}else{
							bestLeftAlignment = alignmentLR;
							bestRightAlignment = alignmentRF;
						}
						
					}else{
						if(myScore==lf[i] || myScore==rr[i]){
							bestLeftAlignment = alignmentLF;
							bestRightAlignment = alignmentRR;
						}else if(myScore==lr[i] || myScore==rf[i]){
							bestLeftAlignment = alignmentLR;
							bestRightAlignment = alignmentRF;
						}
					}
					
				} else if((bestScore-myScore) < distance){
					distance=bestScore-myScore;
				}
					
			}
			
			
			String retval="";
			DecimalFormat twoDForm =  new DecimalFormat("#.##");
			if(bestScore < SCORE_THRES || distance < DIST_THRES ){
				//Logging.info("Unknown sequence " + seq.getName());
				retval = "unknown:"+Double.valueOf(twoDForm.format(bestScore))+":"+Double.valueOf(twoDForm.format(distance))+"|0-0:0-0|";
				seq.setName(retval + seq.getName());

				if(print)
					seq.print(streamToFile[nSamples]);
			}
			//if the best (sum of both ends) alignment in template sequence is greater than in complement
			else {
//				Logging.info("Sequence " + seq.getName() + " might belongs to sample " + barCodesLeft.get(bestIndex).getName() + " with score=" + bestScore);
				if(bestIndex<nSamples){
					retval = barCodesLeft.get(bestIndex).getName()+":"+Double.valueOf(twoDForm.format(bestScore))+":"+Double.valueOf(twoDForm.format(distance))+"|";
					int s1 = bestLeftAlignment.getStart1(),
						e1 = bestLeftAlignment.getStart1()+bestLeftAlignment.getSequence1().length-bestLeftAlignment.getGaps1(),
						e2 = seq.length()-1-bestRightAlignment.getStart1(),
						s2 = seq.length()-1-(bestRightAlignment.getStart1()+bestRightAlignment.getSequence1().length-bestRightAlignment.getGaps1());
					retval += s1+"-"+e1+":"+s2+"-"+e2+"|";
					seq.setName(retval + seq.getName());
					
					if(script && processes[bestIndex]!=null && processes[bestIndex].isAlive())
						seq.print(streamToScript[bestIndex]);
					if(print)
						seq.print(streamToFile[bestIndex]);
				}
			}

			System.out.println(seq.getName());
			printAlignment(bestLeftAlignment);
			System.out.println();
			printAlignment(bestRightAlignment);
			System.out.println("\n==================================================================================\n");
		}

		System.out.println("Done all input");
		for (int i = 0; i < nSamples;i++){
			if(script && processes[i]!=null && processes[i].isAlive()){
				streamToScript[i].close();
				processes[i].waitFor();
			}
			if(print)
				streamToFile[i].close();
		}
		if(print)
			streamToFile[nSamples].close();
		System.out.println("Done every thing");
		reader.close();
	}

		//display jaligner.Alignment. TODO: convert to ours
		public void printAlignment(SWGAlignment alignment){
			String 	origSeq1 = alignment.getOriginalSequence1().toString(),
					origSeq2 = alignment.getOriginalSequence2().toString(),
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
			System.out.println(alignment.getSummary());
			System.out.println(seq1);
			System.out.println(mark);
			System.out.println(seq2);
		}

}
