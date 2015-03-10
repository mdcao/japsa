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
 * 21/07/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq.nanopore;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.Socket;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.HashSet;

import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.data.time.Second;
import org.jfree.data.time.TimePeriod;
import org.jfree.data.time.TimeTableXYDataset;

import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format. 
 * @author minhduc
 *
 */
@Deployable(scriptName = "jsa.np.f5reader", scriptDesc = "Extract nanopore data (fastq/fasta and native data) from h5 files")
public class NanoporeReaderStream
{
	public static void main(String[] args) throws OutOfMemoryError, Exception {
		/*********************** Setting up script ****************************/
		Deployable annotation = NanoporeReaderStream.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options]",
				annotation.scriptDesc());

		cmdLine.addString("output", "-",
				"Name of the output file, -  for stdout");
		cmdLine.addInt("minLength", 0, 
				"Minimum sequence length");
		cmdLine.addBoolean("stats", false, "Compute statistics of reads");
		cmdLine.addBoolean("number", false, "Add a unique number to read name");
		cmdLine.addString("f5list",null, "File containing list of fast5 files, one file per line");
		cmdLine.addString("folder",null, "The download folder");
		cmdLine.addBoolean("fail",false, "Include fail reads");		
		cmdLine.addString("pFolderName",null, "Folder to move processed files to");
		cmdLine.addBoolean("GUI",false, "Whether run with GUI");

		cmdLine.addString("netAddress",null, "Network Address");		
		cmdLine.addInt("netPort", 3645,  "Network port");

		cmdLine.addInt("interval", 30,  "Interval between check in seconds");
		cmdLine.addInt("age", 30,  "The file has to be this old in seconds");

		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");
		String f5list = cmdLine.getStringVal("f5list");
		String folder = cmdLine.getStringVal("folder");
		String pFolderName = cmdLine.getStringVal("pFolderName");
		int minLength  = cmdLine.getIntVal("minLength");
		boolean stats  = cmdLine.getBooleanVal("stats");
		boolean number  = cmdLine.getBooleanVal("number");
		boolean GUI  = cmdLine.getBooleanVal("GUI");
		boolean fail  = cmdLine.getBooleanVal("fail");
		String netAddress = cmdLine.getStringVal("netAddress");
		int netPort  = cmdLine.getIntVal("netPort");

		int interval = cmdLine.getIntVal("interval") * 1000;//in second
		int age = cmdLine.getIntVal("age") * 1000;//in second


		if (folder == null && f5list == null){
			Logging.exit("One of folder and f5list has to be set", 1);
		}

		NanoporeReaderStream reader = new NanoporeReaderStream();
		if (GUI){ 
			System.setProperty("java.awt.headless", "false");
			TimeTableXYDataset dataset = new TimeTableXYDataset();

			GUIReader mGUI = new GUIReader(reader,dataset);			
			new Thread(mGUI).start();

			JFreeChart chart = ChartFactory.createStackedXYAreaChart(
					"Nanopore Fast5 Reader",      // chart title
					"Time",             // domain axis label
					"read number",                   // range axis label
					dataset   
					);			

			final StackedXYAreaRenderer render = new StackedXYAreaRenderer();

			DateAxis domainAxis = new DateAxis();
			domainAxis.setAutoRange(true);
			domainAxis.setDateFormatOverride(new SimpleDateFormat("HH:mm:ss"));
			//domainAxis.setTickUnit(new DateTickUnit(DateTickUnitType.SECOND, 1));

			XYPlot plot = (XYPlot) chart.getPlot();
			plot.setRenderer(render);
			plot.setDomainAxis(domainAxis);
			plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
			plot.setForegroundAlpha(0.5f);

			NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setNumberFormatOverride(new DecimalFormat("#,###.#"));
			rangeAxis.setAutoRange(true);

			JFrame frame = new JFrame("NP Reader");
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			JPanel mainPanel = new JPanel(new BorderLayout());
			ChartPanel chartLabel = new ChartPanel(chart);

			frame.getContentPane().add(mainPanel);
			mainPanel.add(chartLabel);			

			JPanel controlPanel = new JPanel(new GridLayout(0,2));
			mainPanel.add(controlPanel, BorderLayout.PAGE_END);

			controlPanel.add(new JLabel("Total File Number"));			
			controlPanel.add(mGUI.fileNumberField);

			controlPanel.add(new JLabel("Pass File Number"));
			controlPanel.add(mGUI.pFileNumberField);

			controlPanel.add(new JLabel("Fail File Number"));
			controlPanel.add(mGUI.fFileNumberField);				

			controlPanel.add(new JLabel("2D Reads"));
			controlPanel.add(mGUI.twoDReadField);

			controlPanel.add(new JLabel("Complement Reads"));			
			controlPanel.add(mGUI.compReadField);

			controlPanel.add(new JLabel("Template Reads"));			
			controlPanel.add(mGUI.tempReadField);

			/****************************************************
			JPanel iPanel = new JPanel();
			JPanel oPanel = new JPanel();

			iPanel.setBorder(BorderFactory.createTitledBorder("Input"));
			iPanel.setLayout(new GridLayout(0,3));

			iPanel.add(new JCheckBox("Pass Folder"));
			iPanel.add(new JFileChooser());
			iPanel.add(new JFileChooser());

			iPanel.add(new JCheckBox("Fail Folder"));		

			controlPanel.add(iPanel);
			controlPanel.add(oPanel);


			oPanel.setBorder(BorderFactory.createTitledBorder("Output"));			

			/****************************************************/

			frame.pack();
			frame.setVisible(true);  	
		}

		reader.stats = stats;
		reader.number = number;
		reader.minLength = minLength;

		reader.interval = interval;
		reader.age = age;

		reader.f5List = f5list;
		//if (folder != null)
		reader.folder = folder;		
		reader.doFail = fail;

		if (netAddress != null){
			try{
				reader.networkOS = new SequenceOutputStream((new Socket(netAddress, netPort)).getOutputStream());
			}catch(Exception e){
				Logging.error("Sending over network fail : " + e.getMessage());
			}
		}
		reader.sos = SequenceOutputStream.makeOutputStream(output);
		reader.readFastq(pFolderName);
		reader.sos.close();

	}//main


	double tempLength = 0, compLength = 0, twoDLength = 0;
	int tempCount = 0, compCount = 0, twoDCount = 0;
	IntArray lengths = new IntArray();
	int fileNumber = 0;
	int passNumber = 0, failNumber = 0;
	SequenceOutputStream sos;
	SequenceOutputStream networkOS = null;
	boolean stats, number;
	String f5List = null, folder = null;
	int minLength = 0;
	boolean wait = true;
	int interval = 1000, age = 1000;
	boolean doFail = false;

	public boolean readFastq2(String fileName) throws IOException{
		Logging.info("Open " + fileName);
		try{					
			NanoporeReader npReader = new NanoporeReader(fileName);
			npReader.readFastq();
			npReader.close();

			//Get time & date
			String log = npReader.getLog();					
			if (log != null){
				String [] toks = log.split("\n");
				if (toks.length > 0)
					toks = toks[toks.length - 1].split(",");

				log = toks[0];
			}else
				log = "";

			FastqSequence fq;

			fq = npReader.getSeq2D();
			if (fq != null && fq.length() >= minLength){
				fq.setName((number?(fileNumber *3) + "_":"") + fq.getName() + " " + log);
				fq.print(sos);
				if (stats){						
					lengths.add(fq.length());	
					twoDCount ++;
				}
			}

			fq = npReader.getSeqTemplate();
			if (fq != null && fq.length() >= minLength){
				fq.setName((number?(fileNumber *3 + 1) + "_":"") + fq.getName() + " " + log);
				fq.print(sos);
				if (stats){						
					lengths.add(fq.length());	
					tempCount ++;
				}
			}

			fq = npReader.getSeqComplement();
			if (fq != null && fq.length() >= minLength){						
				fq.setName((number?(fileNumber *3 + 2) + "_":"") + fq.getName() + " " + log);						
				fq.print(sos);
				if (stats){						
					lengths.add(fq.length());	
					compCount ++;
				}
			}

			fileNumber ++;			
		}catch (Exception e){
			Logging.error("Problem with reading " + fileName + ":" + e.getMessage());
			return false;
		}		
		return true;
	}

	public boolean moveFile(File f, String pFolder){
		String fName = f.getName();
		if (f.renameTo(new File(pFolder + fName))){
			Logging.info("Move " + fName + " to " + pFolder);
			return true;
		}
		else
			return false;
	}


	/**
	 * Read read sequence from a list of fast5 files.
	 * @param fileList
	 * @param sos : output stream
	 * @param stats: print out statistics
	 * @throws IOException 
	 */
	public void readFastq(String pFolder) throws IOException{
		if (pFolder != null ){
			pFolder = pFolder + File.separatorChar;
			Logging.info("Copy to " + pFolder);
		}

		if (f5List != null){
			Logging.info("Reading in file " + f5List);
			BufferedReader bf = SequenceReader.openFile(f5List);
			String fileName;
			while ((fileName = bf.readLine())!=null){
				readFastq2(fileName);

				//Move to done folder
				if (pFolder != null){
					moveFile(new File(fileName),  pFolder);
				}					
			}//while
			bf.close();
		}else{//folder
			HashSet<String> filesDone = new HashSet<String>();

			File mainFolder = new File(folder);
			File passFolder = new File(folder + File.separatorChar + "pass");
			File failFolder = new File(folder + File.separatorChar + "fail");			

			while (wait){
				try {
					Thread.sleep(interval);
				} catch (InterruptedException e) {					
					e.printStackTrace();
				}

				//Do main
				long now = System.currentTimeMillis();
				File [] fileList = mainFolder.listFiles();
				Logging.info("Reading in folder " + mainFolder.getAbsolutePath());
				if (fileList!=null){
					for (File f:fileList){
						//directory
						if (!f.isFile())
							continue;						

						if (!f.getName().endsWith("fast5"))
							continue;						

						//File too new
						if (now - f.lastModified() < age)
							continue;

						//if processed already
						String sPath = f.getAbsolutePath();					
						if (filesDone.contains(sPath))
							continue;

						readFastq2(sPath);						
						filesDone.add(sPath);	
						if (pFolder != null){
							moveFile(f,  pFolder);
						}//if		
					}//for			
				}//if

				//Pass folder
				now = System.currentTimeMillis();
				Logging.info("Reading in folder " + passFolder.getAbsolutePath());
				fileList = passFolder.listFiles();
				if (fileList!=null){
					for (File f:fileList){
						//directory
						if (!f.isFile())
							continue;

						if (!f.getName().endsWith("fast5"))
							continue;

						//File too new
						if (now - f.lastModified() < age)
							continue;

						//if processed already
						String sPath = f.getAbsolutePath();					
						if (filesDone.contains(sPath))
							continue;

						readFastq2(sPath);
						passNumber ++;
						filesDone.add(sPath);	
						if (pFolder != null){
							moveFile(f,  pFolder);
						}//if			
					}//for			
				}//if

				//Fail folder
				if (!doFail)
					continue;					
				now = System.currentTimeMillis();
				Logging.info("Reading in folder " + failFolder.getAbsolutePath());
				fileList = failFolder.listFiles();
				if (fileList!=null){
					for (File f:fileList){
						//directory
						if (!f.isFile())
							continue;

						if (!f.getName().endsWith("fast5"))
							continue;

						//File too new
						if (now - f.lastModified() < age)
							continue;

						//if processed already
						String sPath = f.getAbsolutePath();					
						if (filesDone.contains(sPath))
							continue;

						readFastq2(sPath);	
						failNumber ++;
						filesDone.add(sPath);	
						if (pFolder != null){
							moveFile(f,  pFolder);
						}//if								
					}//for			
				}//if
			}//while			
		}

		if (stats){
			Logging.info("Getting stats ... ");
			int [] ls = lengths.toArray();
			Arrays.sort(ls);

			long baseCount = 0;						
			for (int i = 0; i < ls.length; i++)
				baseCount += ls[i];

			double mean = baseCount / ls.length;
			double median = ls[ls.length/2];
			long sum = 0;
			int quantile1st = 0, quantile2nd = 0, quantile3rd = 0;
			for (int i = 0; i < ls.length; i++){
				sum += ls[i];
				if (quantile1st == 0 && sum >= baseCount / 4)
					quantile1st = i;

				if (quantile2nd == 0 && sum >= baseCount / 2)
					quantile2nd = i;

				if (quantile3rd == 0 && sum >= baseCount * 3/ 4)
					quantile3rd = i;
			}

			Logging.info("Open " + fileNumber + " files");
			Logging.info("Read count = " + ls.length + "(" + tempCount + " temppate, " + compCount + " complements and " + twoDCount +"  2D)");
			Logging.info("Base count = " + baseCount);		
			Logging.info("Longest read = " + ls[ls.length - 1] + ", shortest read = " + ls[0]);
			Logging.info("Average read length = " + mean);
			Logging.info("Median read length = " + median);
			Logging.info("Quantile first = " + ls[quantile1st] + " second = " + ls[quantile2nd] + " third = " + ls[quantile3rd]);
		}
	}

	static class GUIReader implements Runnable {
		NanoporeReaderStream reader;
		TimeTableXYDataset dataSet;		

		JTextField fileNumberField = new JTextField("0");
		JTextField pFileNumberField = new JTextField("0");
		JTextField fFileNumberField = new JTextField("0");
		JTextField twoDReadField = new JTextField("0");
		JTextField compReadField = new JTextField("0");
		JTextField tempReadField = new JTextField("0");		


		public GUIReader (NanoporeReaderStream r, TimeTableXYDataset dataset){
			reader = r;
			this.dataSet = dataset;

			fileNumberField.setEnabled(false);
			pFileNumberField.setEnabled(false);
			fFileNumberField.setEnabled(false);
			twoDReadField.setEnabled(false);
			compReadField.setEnabled(false);
			tempReadField.setEnabled(false);
		}		

		public void run() {
			while(true) {				                
				//synchronized(reader) {//avoid concurrent update					
				TimePeriod period = new Second();
				dataSet.add(period, reader.twoDCount,"2D");
				dataSet.add(period, reader.compCount,"complement");
				dataSet.add(period, reader.tempCount,"template");


				fileNumberField.setText(reader.fileNumber+"");	                
				pFileNumberField.setText(reader.passNumber+"");
				fFileNumberField.setText(reader.failNumber+"");
				twoDReadField.setText(reader.twoDCount+"");
				compReadField.setText(reader.compCount+"");
				tempReadField.setText(reader.tempCount+"");	

				//}  
				try {
					Thread.sleep(1000);
				} catch (InterruptedException ex) {
					Logging.error(ex.getMessage());
				}
			}
		}
	}
}
