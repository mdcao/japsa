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

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.DoubleArray;
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
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;

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
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.Rengine;


/**
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.np.speciesTyping", scriptDesc = "Species typing using Nanopore Sequencing")
public class SpeciesMixtureTyping {
	/**
	 * @param args
	 * @throws InterruptedException 
	 * @throws Exception 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws IOException, InterruptedException{
		/*********************** Setting up script ****************************/
		Deployable annotation = SpeciesMixtureTyping.class.getAnnotation(Deployable.class);		 		
		CommandLine cmdLine = new CommandLine("\nUsage: " + annotation.scriptName() + " [options]", annotation.scriptDesc());		
		/**********************************************************************/		

		cmdLine.addString("output", "output.dat",  "Output file");		
		cmdLine.addString("bamFile", null,  "The bam file");		
		cmdLine.addString("indexFile", null,  "indexFile ");
		cmdLine.addString("hours", null,  "The file containging hours against yields, if set will output acording to tiime");
		cmdLine.addBoolean("GUI", false,  "Run on GUI");
		cmdLine.addInt("number", 50,  "Number of reads");

		cmdLine.addInt("timestamp", 0,  "Timestamp to check, if <=0 then use read number instead");
		//cmdLine.addInt("read", 500,  "Number of reads before a typing, NA if timestamp is set");

		cmdLine.addInt("sim", 0,  "Scale for simulation");
		cmdLine.addDouble("qual", 0,  "Minimum alignment quality");
		

		args = cmdLine.stdParseLine(args);			
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String bamFile = cmdLine.getStringVal("bamFile");			
		String indexFile = cmdLine.getStringVal("indexFile");
		String hours = cmdLine.getStringVal("hours");
		boolean GUI = cmdLine.getBooleanVal("GUI");
		int number  = cmdLine.getIntVal("number");
		double qual  = cmdLine.getDoubleVal("qual");
		
		SpeciesMixtureTyping paTyping = new SpeciesMixtureTyping(GUI);

		paTyping.simulation = cmdLine.getIntVal("sim");
		paTyping.qual = qual;

		if (hours !=null){
			BufferedReader bf = SequenceReader.openFile(hours);
			String line = bf.readLine();//first line -> ignore
			paTyping.hoursArray = new IntArray();
			paTyping.readCountArray = new IntArray();

			while ((line = bf.readLine())!= null){
				String [] tokens = line.split("\\s+");
				int hrs = Integer.parseInt(tokens[0]);
				int readCount = Integer.parseInt(tokens[2]);

				paTyping.hoursArray.add(hrs);
				paTyping.readCountArray.add(readCount);	
			}
			bf.close();
		}

		//	paTyping.prefix = prefix;
		paTyping.countsOS = SequenceOutputStream.makeOutputStream(output);
		paTyping.preTyping(indexFile);
		paTyping.typing(bamFile, number);
		paTyping.countsOS.close();
		paTyping.close();
	}

	boolean withGUI = false;
	double qual = 0;

	Rengine rengine;

	IntArray hoursArray = null;
	IntArray readCountArray = null;

	int currentReadCount = 0;
	int currentReadAligned = 0;
	long currentBaseCount = 0;

	int arrayIndex = 0;
	String prefix;
	SequenceOutputStream countsOS;
	YIntervalSeriesCollection dataset = new YIntervalSeriesCollection();
	
	long firstReadTime = 0;
	JLabel timeLabel;

	/////////////////////////////////////////////////////////////////////////////

	long startTime;
	public SpeciesMixtureTyping(boolean withGUI){
		this.withGUI = withGUI;
		rengine = new Rengine (new String [] {"--no-save"}, false, null);
		if (!rengine.waitForR()){
			Logging.exit("Cannot load R",1);            
		}    
		rengine.eval("library(MultinomialCI)");
		rengine.eval("alpha<-0.05");

		Logging.info("REngine ready");

		if (withGUI){ 
			System.setProperty("java.awt.headless", "false");

			GUITyping myGen = new GUITyping(this);
			new Thread(myGen).start();
			//dataset.addSeries(s1);
			//dataset.addSeries(s2);
			JFreeChart chart = ChartFactory.createTimeSeriesChart(
					"Species Typing",
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

			JFrame frame = new JFrame("Species Typing");
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			ChartPanel label = new ChartPanel(chart);
			frame.getContentPane().add(label,BorderLayout.CENTER);
			//Suppose I add combo boxes and buttons here later

			timeLabel = new JLabel("Time lapsed : waiting  Total reads: " + currentReadCount + "  Aligned reads: " + currentReadAligned);
			frame.getContentPane().add(timeLabel, BorderLayout.SOUTH);
			frame.pack();
			frame.setVisible(true);  	
		}

		startTime = System.currentTimeMillis();
	}

	public void close(){
		rengine.end();
	}
	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	static class SpeciesCount implements Comparable<SpeciesCount>{
		String species;
		int count = 0;

		SpeciesCount (String s){
			species = s;
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(SpeciesCount o) {		
			return o.count - count;
		}

	}

	
	HashMap<String, String> seq2Species = new HashMap<String, String>();
	HashMap<String, SpeciesCount> species2Count = new HashMap<String, SpeciesCount>();
	ArrayList<String> speciesList = new ArrayList<String>(); 


	private void preTyping(String indexFile)throws IOException{
		BufferedReader bf = SequenceReader.openFile(indexFile);
		String line = "";
		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.split(" ");
			String sp =  toks[0];
			String seq = toks[1].substring(1);

			if (seq2Species.put(seq, sp) != null)
				throw new RuntimeException("sequence " + seq +" presents multiple time");

			if (species2Count.get(sp) == null){
				species2Count.put(sp,new SpeciesCount(sp));
			}			
		}//while
		bf.close();
		Logging.info(seq2Species.size() + "   " + species2Count.size());
		speciesList.addAll(species2Count.keySet());

		//Write header
		countsOS.print("step\treads\tbases\tspecies\tprob\terr\n");
		//for (String species:speciesList){
		//	countsOS.print("\t" + species);
		//}	
		//countsOS.println();
	}

	private void simpleAnalysisCurrent(int currentRead) throws IOException{		
		int step = currentRead;
		if (hoursArray != null) 
			step = hoursArray.get(arrayIndex);

		int sum = 0;
		double [] count = new double[speciesList.size()];
		for (int i = 0; i < count.length;i++){			
			count[i] = species2Count.get(speciesList.get(i)).count;			
			sum += count[i];
		}
		DoubleArray countArray = new DoubleArray();
		ArrayList<String> speciesArray = new ArrayList<String> ();

		int minCount = Math.max(1,sum/100);
		for (int i = 0; i < count.length;i++){			
			if (count[i] >= minCount){
				countArray.add(count[i]);
				speciesArray.add(speciesList.get(i));
				Logging.info(step+" : " + speciesList.get(i) + " == " + count[i]);
			}
		}		
		//if (countArray.size() > 10) return;
		countArray.add(1);
		speciesArray.add("others");		

		rengine.assign("count", countArray.toArray());
		rengine.eval("tab = multinomialCI(count,alpha)");        
		REXP tab  = rengine.eval("tab",true);  
		double [][] results = tab.asDoubleMatrix();

		//countsOS.print(step);
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

		for (int i = 0; i < results.length;i++){
			if (results[i][0] <= 0.00001)
				continue;

			double mid = (results[i][0] + results[i][1])/2;
			double err = mid - results[i][0];

			countsOS.print(step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + speciesArray.get(i).replaceAll("_"," ") + "\t" + mid +"\t" + err);
			countsOS.println();
		}

		countsOS.flush();
		Logging.info(step+"  " + countArray.size());
	}

	int simulation = 0;
	private void guiAnalysisCurrent( HashMap<String, YIntervalSeries> speciesSeries){		

		int sum = 0;
		double [] count = new double[speciesList.size()];
		for (int i = 0; i < count.length;i++){			
			count[i] = species2Count.get(speciesList.get(i)).count;			
			sum += count[i];
		}
		DoubleArray countArray = new DoubleArray();
		ArrayList<String> speciesArray = new ArrayList<String> ();

		for (int i = 0; i < count.length;i++){			
			if (count[i] >= sum/50){
				countArray.add(count[i]);
				speciesArray.add(speciesList.get(i));
			}
		}		
		
		Calendar cal = Calendar.getInstance();    	
    	SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");    	
		Logging.info(sdf.format(cal.getTime()) + " Found " + countArray.size() + " counts");
		if (countArray.size() > 20) return;

		countArray.add(1);
		speciesArray.add("others");

		rengine.assign("count", countArray.toArray());
		rengine.eval("tab = multinomialCI(count,alpha)");        
		REXP tab  = rengine.eval("tab",true);  
		double [][] results = tab.asDoubleMatrix();

		//countsOS.print(step);
		long x = System.currentTimeMillis(); 
		for (int i = 0; i < results.length;i++){
			String species = speciesArray.get(i);
			if (species.equals("others"))
				continue;
			
			double mid = (results[i][0] + results[i][1])/2;
			YIntervalSeries series = speciesSeries.get(species);
			if (series == null){
				series = new YIntervalSeries(species);
				speciesSeries.put(species, series);
				dataset.addSeries(series);
			}
			series.add(x, mid, results[i][0], results[i][1]);        	
			//countsOS.print("\t" + mid +"\t" + err);
		}
		//countsOS.println();
		//countsOS.flush();
		//Logging.info(step+"  " + countArray.size());
	}

	private void typing(String bamFile, int readNumber) throws IOException, InterruptedException{
		if (readNumber <= 0)
			readNumber = 1;


		String readName = "";
		//Read the bam file		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();			

		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			if (firstReadTime <=0)
				firstReadTime = System.currentTimeMillis();
			
			if (!sam.getReadName().equals(readName)){
				readName = sam.getReadName();

				currentReadCount ++;
				currentBaseCount += sam.getReadLength();

				if (hoursArray != null){
					if (arrayIndex < hoursArray.size() && currentReadCount >= this.readCountArray.get(arrayIndex)){
						simpleAnalysisCurrent(currentReadCount);
						arrayIndex ++;
					}
				}else{				
					if (currentReadCount % readNumber == 0){
						simpleAnalysisCurrent(currentReadCount);
					}
				}
			}

			if (sam.getReadUnmappedFlag()){				
				continue;			
			}
			
			if (sam.getMappingQuality() < this.qual)
				continue;
			
			currentReadAligned ++;

			String refSequence = sam.getReferenceName();
			String species = seq2Species.get(refSequence);
			if (species == null){
				throw new RuntimeException(" Can find species with ref " + refSequence + " line " + currentReadCount );
			}
			SpeciesCount sCount = species2Count.get(species);
			if (sCount == null){
				throw new RuntimeException(" Can find record with species " + species + " line " + currentReadCount );
			}

			synchronized(this) {
				sCount.count ++;
			}			

		}//while
		
		//final run
		simpleAnalysisCurrent(currentReadCount);

		samIter.close();
		samReader.close();
	}

	static class GUITyping implements Runnable {
		SpeciesMixtureTyping typing;
		public GUITyping (SpeciesMixtureTyping typ){
			this.typing = typ;
		}	
		HashMap<String, YIntervalSeries> speciesSeries = new HashMap<String, YIntervalSeries>();


		public void run() {
			long lastRun = 0;
			while(true) {	
				long delay = System.currentTimeMillis() - lastRun;
				if (delay < 5000){
					try {
						Thread.sleep((5000 - delay));
					} catch (InterruptedException ex) {
						System.out.println(ex);
					}					
				}
				synchronized(this.typing) {//avoid concurrent update
					lastRun = System.currentTimeMillis();
					typing.guiAnalysisCurrent(speciesSeries);
					if (typing.firstReadTime > 0){
						long lapsedTime = (System.currentTimeMillis() - typing.firstReadTime) / 1000;
						long hours = lapsedTime / 3600;
						lapsedTime = lapsedTime % 3600;
						long mins  = lapsedTime / 60;
						lapsedTime = lapsedTime % 60;						
						typing.timeLabel.setText("Time lapsed : " + hours + ":" + (mins < 10?"0":"") + mins + ":" + (lapsedTime < 10?"0":"") + lapsedTime 
								+ "   Total reads: " + typing.currentReadCount + "  Aligned reads: " + typing.currentReadAligned);
					}
				}  

			}
		}
	}
}