/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.np;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
//import japsa.util.BetaBinomialModel;
import japsa.util.CommandLine;
import japsa.util.HTSUtilities;
import japsa.util.IntArray;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.JFrame;
import javax.swing.JLabel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.data.xy.YIntervalSeries;
import org.jfree.data.xy.YIntervalSeriesCollection;



/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.np.geneStrainTyping", scriptDesc = "Strain typing using present/absence of gene")
public class GeneStrainTyping {

	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = GeneStrainTyping.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		


		cmdLine.addString("output", "output.dat",  "Output file");
		cmdLine.addString("profile", null,  "Output file containing gene profile of all strains");
		cmdLine.addString("bamFile", null,  "The bam file");
		cmdLine.addString("geneFile", null,  "The gene file");

		cmdLine.addInt("top", 10,  "The number of top strains");
		cmdLine.addInt("scoreThreshold", 0,  "The alignment score threshold");
		cmdLine.addString("tmp", "tmp/t",  "Temporary folder");
		cmdLine.addString("hours", null,  "The file containging hours against yields, if set will output acording to tiime");

		cmdLine.addInt("timestamp", 0,  "Timestamp to check, if <=0 then use read number instead");
		cmdLine.addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");

		cmdLine.addBoolean("twodonly", false,  "Use only two dimentional reads");
		cmdLine.addInt("sim", 0,  "Scale for simulation");
		cmdLine.addBoolean("GUI", false,  "Run on GUI");

		args = cmdLine.stdParseLine(args);		
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String profile = cmdLine.getStringVal("profile");
		String bamFile = cmdLine.getStringVal("bam");
		String geneFile = cmdLine.getStringVal("geneFile");		
		String tmp = cmdLine.getStringVal("tmp");
		String hours = cmdLine.getStringVal("hours");
		int top = cmdLine.getIntVal("top");		
		int read = cmdLine.getIntVal("read");
		boolean GUI = cmdLine.getBooleanVal("GUI");
		int timestamp = cmdLine.getIntVal("timestamp");

		{
			GeneStrainTyping paTyping = new GeneStrainTyping(GUI);	
			paTyping.simulation = cmdLine.getIntVal("sim");
			paTyping.prefix = tmp;
			paTyping.readNumber = read;
			if (hours !=null){
				BufferedReader bf = SequenceReader.openFile(hours);
				String line = bf.readLine();//first line
				paTyping.hoursArray = new IntArray();
				paTyping.readCountArray = new IntArray();

				while ((line = bf.readLine())!= null){
					String [] tokens = line.split("\\s");
					int hrs = Integer.parseInt(tokens[0]);
					int readCount = Integer.parseInt(tokens[2]);

					paTyping.hoursArray.add(hrs);
					paTyping.readCountArray.add(readCount);	
				}
			}


			if (paTyping.readNumber < 1)
				paTyping.readNumber = 1;

			paTyping.datOS = SequenceOutputStream.makeOutputStream(output);
			paTyping.datOS.print("step\treads\tbases\tstrain\tprob\tlow\thigh\tgenes\n");
			paTyping.readGenes(geneFile);
			paTyping.readKnowProfiles(profile);
			Logging.info("Read in " + paTyping.profileList.size() + " gene profiles");

			paTyping.timestamp = timestamp;

			if (GUI)
				paTyping.startGUI();

			paTyping.typing(bamFile,  top);
			paTyping.datOS.close();
		}
	}


	/////////////////////////////////////////////////////////////////////////////


	HashSet<String> addedGenes = new HashSet<String>(); 
	PresenceAbsence lcTyping;

	ArrayList<GeneProfile> profileList;	

	ArrayList<Sequence> geneList;
	HashMap<String, Sequence> geneMap;

	HashMap<String, ArrayList<Sequence>> alignmentMap;

	HashSet<String> targetGenes;


	String prefix = "tmp";	
	int simulation = 0;

	int readNumber = 100;
	SequenceOutputStream datOS = null;


	int currentReadCount = 0;
	long currentBaseCount = 0;
	int currentReadAligned = 0;

	IntArray hoursArray = null;
	IntArray readCountArray = null;
	int arrayIndex = 0;

	YIntervalSeriesCollection dataset = new YIntervalSeriesCollection();
	long startTime;
	long firstReadTime = 0;
	JLabel timeLabel;
	boolean withGUI = false;
	int timestamp = 5000;


	public GeneStrainTyping(boolean withGUI){
		this.withGUI = withGUI;
		startTime = System.currentTimeMillis();		
	}

	public void startGUI(){
		System.setProperty("java.awt.headless", "false");

		GUIStrainTyping myGen = new GUIStrainTyping(this, timestamp);
		new Thread(myGen).start();
		//dataset.addSeries(s1);
		//dataset.addSeries(s2);
		JFreeChart chart = ChartFactory.createTimeSeriesChart(
				"Strain Typing",
				"Time",
				"Value",
				dataset,
				true,
				true,
				false
				);
		final XYPlot plot = chart.getXYPlot();
		plot.setRenderer(new XYErrorRenderer());

		//System.out.println(chart.getXYPlot().getRenderer().getClass().getCanonicalName());

		ValueAxis axis = plot.getDomainAxis();
		axis.setAutoRange(true);
		axis.setAutoRangeMinimumSize(6000);

		ValueAxis yAxis = plot.getRangeAxis();
		yAxis.setRange(0.0, 1.0);
		//axis.set

		//axis.setFixedAutoRange(6000.0);

		JFrame frame = new JFrame("Strain Typing");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		ChartPanel label = new ChartPanel(chart);
		frame.getContentPane().add(label,BorderLayout.CENTER);
		//Suppose I add combo boxes and buttons here later

		timeLabel = new JLabel("Time lapsed : waiting  Total reads: " + currentReadCount + "  Aligned reads: " + currentReadAligned);
		frame.getContentPane().add(timeLabel, BorderLayout.SOUTH);
		frame.pack();
		frame.setVisible(true);  
	}


	/**
	 * Analysis of drug resistance profile
	 * @param mcoordFile
	 * @throws IOException
	 */


	/**
	 * Read genes from gene file to a list + map: for random access
	 * @param geneFile
	 * @throws IOException
	 */
	private void readGenes(String geneFile) throws IOException{
		geneList = SequenceReader.readAll(geneFile, Alphabet.DNA());		 
		geneMap = new HashMap<String, Sequence>();

		for (Sequence gene:geneList){
			geneMap.put(gene.getName(), gene);		
		}
	}


	/**
	 * Read the gene profile of all known strains (for strain typing) 
	 * @param profileFile
	 * @throws IOException
	 */
	private void readKnowProfiles(String profileFile) throws IOException{
		BufferedReader reader = new BufferedReader (new FileReader(profileFile));		
		profileList = new ArrayList<GeneProfile>(); 
		String line;
		String currentStrainID = "";
		GeneProfile profile = null;

		while ((line = reader.readLine()) != null){
			if (line.startsWith("#"))
				continue;
			String [] toks = line.trim().split("\t");
			String strainID = toks[0];
			String geneFamID = toks[1];


			if (strainID.equals(currentStrainID)){
				profile.addGene(geneFamID);
			}else{
				profile = new GeneProfile(strainID);
				currentStrainID = strainID;
				profile.addGene(geneFamID);
				profileList.add(profile);				
			}

		}		
		reader.close();

		lcTyping = new PresenceAbsence(profileList);

	}




	double[] posterior = new double[0];// = lcTyping.calcPosterior();
	double[][] samp = null;// =lcTyping.calcPosterior(1000);
	double[][] ranges = null;// = lcTyping.getRanges(samp, 0.99);


	/**
	 * Compute the alignment score between a gene and a list of (errornous) reads
	 * that were aligned to the gene. The algorithm is:
	 *  - Get all read sequences (which were previously trimmed)
	 *  - Call a MSA method to alignment them
	 *  - Make a consensus sequence of those reads
	 *  - Align the consensus sequence to the gene sequence using needle
	 *  - Return the needle alignment score
	 * @param gene : the gene sequence
	 * @param readList: an array list of read sequences aligned to this gene

	 * @throws IOException 
	 * @throws InterruptedException 
	 */


	private double alignmentScore2(Sequence gene, ArrayList<Sequence> readList){
		double score = 0;
		if (readList != null){
			score = readList.size();
		}
		return score;
	}


	HashSet<String> mentionedStrain = new HashSet<String>(); 
	double threshold = 0;

	private ArrayList<LCTypingResult> makePresenceTyping(int top) throws IOException, InterruptedException{

		int step = currentReadCount;
		if (hoursArray != null) 
			step = hoursArray.get(arrayIndex);

		//HashSet<String> myGenes = new HashSet<String>();	
		boolean compute = false;
		for (Sequence gene:geneList){			
			ArrayList<Sequence> alignmentList =  alignmentMap.get(gene.getName());
			//This method use the simple scoreing as it involes ten thousands of genes
			if (alignmentScore2(gene, alignmentList) > threshold){
				//myGenes.add(gene.getName());

				if (!addedGenes.contains(gene.getName())){
					lcTyping.likelihood(100,gene.getName());
					addedGenes.add(gene.getName());
					compute = true;//only need to compute if new evidence is observed
				}				
			}
		}
		Logging.info(step + ": Found " + addedGenes.size() + "  " + compute);

		if (compute){
			posterior = lcTyping.calcPosterior();
			samp =lcTyping.calcPosterior(1000);
			ranges = lcTyping.getRanges(samp, 0.99);
		}
		ArrayList<LCTypingResult> lcT = new ArrayList<LCTypingResult>(); 
		for(int i=0; i<posterior.length; i++){
			LCTypingResult lts = new LCTypingResult();
			lts.strainID = lcTyping.spl[i].species;
			lts.postProb = posterior[i];
			lts.l = ranges[i][0];
			lts.h = ranges[i][1];
			lcT.add(lts);
		}
		Collections.sort(lcT);



		for (GeneProfile profile:profileList){
			int TP = 0, FN = 0;
			for (String geneID:profile.genes){
				if (addedGenes.contains(geneID))
					TP ++;
				else
					FN ++;
			}//for geneID
			int FP = addedGenes.size() - TP;
			//estimate p:


			double precision = (TP * 1.0) / (TP + FP);
			double recall = (TP * 1.0) / (TP + FN);			
			double f1 = 2 * precision * recall/ (precision + recall);


			profile.precision = precision;
			profile.recall = recall;//
			profile.f1 = profile.score = f1; 
		}//for profile

		Collections.sort(profileList);

		//delay
		if (simulation > 0){
			long delay = step * 60 * 1000 / simulation - (System.currentTimeMillis() - startTime);
			Logging.info("Step " + step + " delay " + delay/1000);
			if(delay > 0){
				Logging.info("Step " + step + " delay " + delay/1000);			
				try {
					Thread.sleep(delay);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}


		for (int i = 0; i < top && i < lcT.size();i++){			
			GeneProfile profile = profileList.get(i);			
			Logging.info("F\t" + step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + profile.strainID + "\t" + profile.f1 + "\t" + profile.precision +"\t" + profile.recall + "\t"+addedGenes.size());

			LCTypingResult lr  = lcT.get(i);

			//if (lr.postProb < 0.10)
			//	break;
			datOS.print(step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + lr.strainID + "\t" + lr.postProb +"\t" + (lr.postProb - lr.l) + "\t" + (lr.h -lr.postProb)  +"\t"+addedGenes.size());
			datOS.println();			
		}
		datOS.flush();

		return lcT;
	}


	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void typing(String bamFile, int top) throws IOException, InterruptedException{		
		alignmentMap = new HashMap<String, ArrayList<Sequence>> ();

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();				

		String readName = "";
		//A dummy sequence
		Sequence readSequence = new Sequence(Alphabet.DNA(),1,"");
		while (samIter.hasNext()){
			SAMRecord record = samIter.next();
			if (firstReadTime <=0)
				firstReadTime = System.currentTimeMillis();


			//if (this.twoDOnly && !record.getReadName().contains("twodim")){
			//	continue;
			//}

			if (!record.getReadName().equals(readName)){
				readName = record.getReadName();

				currentReadCount ++;	
				currentBaseCount += record.getReadLength();

				if (!withGUI){
					if (hoursArray != null){
						if (arrayIndex < hoursArray.size() && currentReadCount >= this.readCountArray.get(arrayIndex)){						
							makePresenceTyping(top);
							arrayIndex ++;
						}
					}else{				
						if (currentReadCount % readNumber == 0){
							makePresenceTyping(top);
						}
					}
				}
				//Get the read
				if (!record.getReadUnmappedFlag()){
					readSequence = new Sequence(Alphabet.DNA(), record.getReadString(), readName);
					if (record.getReadNegativeStrandFlag()){
						readSequence = Alphabet.DNA.complement(readSequence);
						readSequence.setName(readName);
					}
				}
			}

			if (record.getReadUnmappedFlag())
				continue;			
			//assert: the read sequence is stored in readSequence with the right direction

			currentReadAligned ++;
			String	geneID = record.getReferenceName();
			if (!geneMap.containsKey(geneID))
				continue;

			int refLength =  geneMap.get(geneID).length();


			ArrayList<Sequence> alignmentList = alignmentMap.get(geneID);
			if (alignmentList == null){
				alignmentList = new ArrayList<Sequence>();
				alignmentMap.put(geneID, alignmentList);
			}
			//put the sequence into alignment list

			Sequence readSeq = HTSUtilities.spanningSequence(record, readSequence, refLength,0);

			if (readSeq == null){
				Logging.warn("Read sequence is NULL sequence ");
			}else{
				alignmentList.add(readSeq);
			}


		}//while	
		samIter.close();
		samReader.close();


		makePresenceTyping(top);

	}


	public static class LCTypingResult implements Comparable<LCTypingResult>{
		String strainID;
		double postProb, l, h;
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(LCTypingResult o) {
			double comp = postProb - o.postProb;

			if (comp < 0)
				return 1;
			else if (comp > 0)
				return -1;
			else 
				return 0;			

		}
	}


	/**
	 * Set up: -- read gff files, extract gene sequences, generate profile for each strain
	 * 
	 * @param file
	 * @param out
	 * @param profile
	 * @throws IOException
	 */

	public static class GeneProfile implements Comparable<GeneProfile>{
		String strainID;
		double score = 0;
		double f1 = 0, precision, recall;
		HashSet<String> genes;

		public GeneProfile(String id){
			strainID = id;
			genes = new HashSet<String>();
		}

		public void addGene(String geneID){
			genes.add(geneID);
		}

		public HashSet<String>  getGeneList(){
			return genes;
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(GeneProfile o) {
			double comp = score - o.score;

			if (comp < 0)
				return 1;
			else if (comp > 0)
				return -1;
			else 
				return 0;
		}

		public String strainID(){
			return strainID;
		}
	}


	public static String compare(HashSet<String> s1,HashSet<String> s2){
		String ret = "";
		int count1 = 0,count2=0,count=0;

		for (String st:s1){
			if (s2.contains(st))
				count ++;
			else{ count1++;
			ret = ret + ";" +st;
			}
		}
		ret = ret + "#";

		for (String st:s2){
			if (!s1.contains(st)){			
				count2++;			
				ret = ret + ";" +st;
			}
		}		
		return "Common" + count + "#"+count1+"#"+count2+"#"+ret;
	}


	static class GUIStrainTyping implements Runnable {
		GeneStrainTyping typing;
		int timestamp = 5000;
		public GUIStrainTyping (GeneStrainTyping typ, int timeInteval){
			this.typing = typ;
			this.timestamp = timeInteval;
		}	
		HashMap<String, YIntervalSeries> speciesSeries = new HashMap<String, YIntervalSeries>();


		public void run() {
			long lastRun = 0;
			while(true) {	
				long delay = System.currentTimeMillis() - lastRun;
				if (delay < timestamp){
					try {
						Thread.sleep((timestamp - delay));
					} catch (InterruptedException ex) {
						System.out.println(ex);
					}					
				}
				System.out.println("TICK " + delay);
				lastRun = System.currentTimeMillis();
				synchronized(this.typing) {//avoid concurrent update					
					try {						
						if (typing.firstReadTime > 0){							
							long lapsedTime = (lastRun - typing.firstReadTime) / 1000;
							long hours = lapsedTime / 3600;
							lapsedTime = lapsedTime % 3600;
							long mins  = lapsedTime / 60;
							lapsedTime = lapsedTime % 60;						
							typing.timeLabel.setText("Time lapsed : " + hours + ":" + (mins < 10?"0":"") + mins + ":" + (lapsedTime < 10?"0":"") + lapsedTime 
									+ "   Total reads: " + typing.currentReadCount + "  Aligned reads: " + typing.currentReadAligned);


							ArrayList<LCTypingResult> lcT = typing.makePresenceTyping(10);

							for (LCTypingResult r:lcT){
								double mid = r.postProb;
								if (mid < 0.1)
									continue;
								YIntervalSeries series = speciesSeries.get(r.strainID);


								if (series == null){
									series = new YIntervalSeries(r.strainID);
									speciesSeries.put(r.strainID, series);
									typing.dataset.addSeries(series);
								}
								series.add(lastRun, mid, r.l, r.h);
							}

						}
					} catch (IOException e) {						
						e.printStackTrace();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}

				}  

			}
		}
	}
}