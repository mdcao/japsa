package japsadev.seq.nanopore;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.Logging;

public class Demultiplexer {
	int 	SCAN_WINDOW, 
			DIST_THRES,
			SCORE_THRES; 
	
	ArrayList<Sequence> barCodes;
	ArrayList<Sequence> barCodeComps;
	int nSamples;
	int barcodeLen;
	
	int[] readCount;
	
	public Demultiplexer(String barcodeFile) throws IOException{
		barCodes = SequenceReader.readAll(barcodeFile, Alphabet.DNA());
		nSamples = barCodes.size();
		barcodeLen = barCodes.get(0).length();
		
		barCodeComps = new ArrayList<Sequence> (barCodes.size());
		Sequence barCode = null;
		for(int i=0;i<nSamples;i++){		
			barCode = barCodes.get(i);
			barCodeComps.add(Alphabet.DNA.complement(barCode));
		}
		// Default setting for searching parameters
		SCAN_WINDOW = barcodeLen * 3;
		SCORE_THRES = barcodeLen;
		DIST_THRES = barcodeLen / 3;
		
		readCount = new int[nSamples];
	}
	
	public void setThreshold(int score){
		SCORE_THRES=score;
	}
	/*
	 * Trying to clustering MinION read data into different samples based on the barcode
	 */
	public void clustering(Sequence seq) throws IOException, InterruptedException{

		Sequence s5, s3;
		final double[] 	tf = new double[nSamples],
				tr = new double[nSamples],
				cr = new double[nSamples],
				cf = new double[nSamples];

		Sequence barcodeSeq = new Sequence(Alphabet.DNA4(), barcodeLen, "barcode");
		Sequence tipSeq = new Sequence(Alphabet.DNA4(), SCAN_WINDOW, "tip");

		BarcodeAlignment barcodeAlignment = new BarcodeAlignment(barcodeSeq, tipSeq);

		if(seq.length() < barcodeLen * 2 + 200){
//			Logging.info("Ignoring short sequence " + seq.getName());
			seq.setName("Barcode:unknown:0.0:0.0|" + seq.getName());
			return;
		}
		//alignment algorithm is applied here. For the beginning, Smith-Waterman local pairwise alignment is used

		s5 = seq.subSequence(0, SCAN_WINDOW);
		s3 = seq.subSequence(seq.length()-SCAN_WINDOW,seq.length());


		double bestScore = 0.0, secondBestScore = 0.0,
				otherEndOfBest = 0.0, otherEndOfSecondBest = 0.0;
		//double distance = 0.0; //distance between bestscore and the runner-up
		int bestIndex = nSamples, secondBestIndex = nSamples;

		for(int i=0;i<nSamples; i++){
			Sequence barcode = barCodes.get(i);
			Sequence barcodeComp = barCodeComps.get(i);

			barcodeAlignment.setBarcodeSequence(barcode);				
			barcodeAlignment.setReadSequence(s5);				
			tf[i]=barcodeAlignment.align();				

			barcodeAlignment.setReadSequence(s3);
			tr[i]=barcodeAlignment.align();

			barcodeAlignment.setBarcodeSequence(barcodeComp);
			barcodeAlignment.setReadSequence(s3);
			cr[i]=barcodeAlignment.align();
			barcodeAlignment.setReadSequence(s5);
			cf[i]=barcodeAlignment.align();


			//This is for both end
			//double myScore = Math.max(tf[i], tr[i]) + Math.max(cf[i], cr[i]);

			double myScore = Math.max(Math.max(tf[i], tr[i]), Math.max(cf[i], cr[i]));
			if (myScore > bestScore){
				//Logging.info("Better score=" + myScore);
				secondBestScore = bestScore;
				bestScore = myScore;		
				secondBestIndex = bestIndex;
				bestIndex = i;
			} else if(myScore > secondBestScore){
				secondBestScore = bestScore;
				secondBestIndex = bestIndex;
			}
				
		}
		otherEndOfBest = Math.max(tf[bestIndex], tr[bestIndex]) > Math.max(cf[bestIndex], cr[bestIndex]) ?
				(tf[bestIndex]+tr[bestIndex]-bestScore):(cf[bestIndex]+cr[bestIndex]-bestScore);
		otherEndOfSecondBest = Math.max(tf[secondBestIndex], tr[secondBestIndex]) > Math.max(cf[secondBestIndex], cr[secondBestIndex]) ?
				(tf[secondBestIndex]+tr[secondBestIndex]-secondBestScore):(cf[secondBestIndex]+cr[secondBestIndex]-secondBestScore);
		
		

		String retval="";
		DecimalFormat twoDForm =  new DecimalFormat("#.##");
//		if(bestScore < SCORE_THRES || distance < DIST_THRES){
//			//Logging.info("Confounding sequence " + seq.getName() + " with low grouping score " + bestScore);
//			retval = "Barcode:unknown:"+Double.valueOf(twoDForm.format(bestScore))+":"+Double.valueOf(twoDForm.format(distance))+"|";
//
//		}
//		else {
//			//Logging.info("Sequence " + seq.getName() + " might belongs to sample " + barCodes.get(bestIndex).getName() + " with score=" + bestScore);
//			retval = "Barcode:"+barCodes.get(bestIndex).getName()+":"+Double.valueOf(twoDForm.format(bestScore))+":"+Double.valueOf(twoDForm.format(distance))+"|";
//			readCount[bestIndex]++;
//		}
				
		retval = barCodes.get(bestIndex).getName()+":"+Double.valueOf(twoDForm.format(bestScore))+","+Double.valueOf(twoDForm.format(otherEndOfBest))+"|"
				+barCodes.get(secondBestIndex).getName()+":"+Double.valueOf(twoDForm.format(secondBestScore))+","+Double.valueOf(twoDForm.format(otherEndOfSecondBest))
				+"|";
		
		
		readCount[bestIndex]++;
		seq.setName(retval + seq.getName());

	}
	
	public final class BarcodeAlignment {

		/**
		 * Traceback direction stop
		 */
		public static final byte STOP = 0;
		/**
		 * Traceback direction left
		 */
		public static final byte LEFT = 1;
		/**
		 * Traceback direction diagonal
		 */
		public static final byte DIAGONAL = 2;
		/**
		 * Traceback direction up
		 */
		public static final byte UP = 3;

		public BarcodeAlignment(Sequence s1, Sequence s2) {
			super();
			this.barcodeSequence = s1;
			this.readSequence = s2;	

			m = s1.length() + 1;
			n = s2.length() + 1;

			//Initialize the arrays
			pointers = new byte[m * n];
			sizesOfVerticalGaps = new short[m * n];
			sizesOfHorizontalGaps = new short[m * n];
		}	
		
		Sequence barcodeSequence;
		Sequence readSequence;
		int m,n;
		byte[] pointers;
		short[] sizesOfVerticalGaps;
		short[] sizesOfHorizontalGaps;
		//BLOSSOM62
		//double [][] scores = 
		//	{{4.0,0.0,0.0,0.0},
		//			{0.0,9.0,-3.0,-1.0}, 
		//			{0.0,-3.0,6.0,-2.0},
		//			{0.0,-1.0,-2.0,5.0}
		//	};
		
		//poreFUME's scores
		double openPenalty = 4.7;
		double extendPenalty = 1.6;

		double [][] scores = {
				{  2.7, -4.5, -4.5, -4.5},
				{ -4.5,  2.7, -4.5, -4.5},
				{ -4.5, -4.5,  2.7, -4.5},
				{ -4.5, -4.5, -4.5,  2.7}			
		};
			

		/**
		 * Alignment score at this cell
		 */
		private double cellScore;

		
		public void setBarcodeSequence(Sequence seq){
			barcodeSequence = seq;
		}
		
		public void setReadSequence(Sequence seq){
			readSequence = seq;
		}
		
		

		public double align() {		
			// Initializes the boundaries of the traceback matrix to STOP.
			for (int i = 0, k = 0; i < m; i++, k += n) {
				pointers[k] = STOP;
			}
			for (int j = 1; j < n; j++) {
				pointers[j] = STOP;
			}

			for (int i = 0, k = 0; i < m; i++, k += n) {
				for (int j = 0; j < n; j++) {
					sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
				}
			}
			return construct();
		}

		/**
		 * Constructs directions matrix for the traceback
		 * 
		 * @param barcodeSequence
		 *            sequence #1
		 * @param readSequence
		 *            sequence #2
		 * @param scores
		 *            scoring matrix
		 * @param openPenalty
		 *            open gap penalty
		 * @param extendPenalty
		 *            extend gap penalty
		 * @return The cell where the traceback starts.
		 */
		private double construct() {
			//logger.info("Started...");
			//long start = System.currentTimeMillis();

			double f; // score of alignment x1...xi to y1...yi if xi aligns to yi
			double[] g = new double[n]; // score if xi aligns to a gap after yi
			double h; // score if yi aligns to a gap after xi
			double[] v = new double[n]; // best score of alignment x1...xi to
			// y1...yi
			double vDiagonal;

			g[0] = Float.NEGATIVE_INFINITY;
			h = Float.NEGATIVE_INFINITY;
			v[0] = 0;

			for (int j = 1; j < n; j++) {
				g[j] = Float.NEGATIVE_INFINITY;
				v[j] = 0;
			}

			double similarityScore, g1, g2, h1, h2;

			cellScore = Float.NEGATIVE_INFINITY;
			//Cell cell = new Cell();

			for (int i = 1, k = n; i < m; i++, k += n) {
				h = Float.NEGATIVE_INFINITY;
				vDiagonal = v[0];
				for (int j = 1, l = k + 1; j < n; j++, l++) {
					similarityScore = scores[barcodeSequence.getBase(i-1)][readSequence.getBase(j-1)];

					// Fill the matrices
					f = vDiagonal + similarityScore;

					g1 = g[j] - extendPenalty;
					g2 = v[j] - openPenalty;
					if (g1 > g2) {
						g[j] = g1;
						sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
					} else {
						g[j] = g2;
					}

					h1 = h - extendPenalty;
					h2 = v[j - 1] - openPenalty;
					if (h1 > h2) {
						h = h1;
						sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
					} else {
						h = h2;
					}

					vDiagonal = v[j];
					v[j] = maximum(f, g[j], h, 0);

					// Determine the traceback direction
					if (v[j] == 0) {
						pointers[l] = STOP;
					} else if (v[j] == f) {
						pointers[l] = DIAGONAL;
					} else if (v[j] == g[j]) {
						pointers[l] = UP;
					} else {
						pointers[l] = LEFT;
					}

					// Set the traceback start at the current cell i, j and score
					if (v[j] > cellScore) {
						cellScore = v[j];
						//cell.set(i, j, v[j]);
					}
				}
			}		
			return cellScore;
		}


		/**
		 * Returns the maximum of 4 float numbers.
		 * 
		 * @param a
		 *            float #1
		 * @param b
		 *            float #2
		 * @param c
		 *            float #3
		 * @param d
		 *            float #4
		 * @return The maximum of a, b, c and d.
		 */
		private double maximum(double a, double b, double c, double d) {
			if (a > b) {
				if (a > c) {
					return a > d ? a : d;
				} else {
					return c > d ? c : d;
				}
			} else if (b > c) {
				return b > d ? b : d;
			} else {
				return c > d ? c : d;
			}
		}

	}
}

