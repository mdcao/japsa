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

package japsa.tools.bio.np;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5CompoundDS;
import ncsa.hdf.object.h5.H5ScalarDS;
import japsa.seq.Alphabet.DNA;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.IntArray;
import japsa.util.JapsaException;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format. 
 * @author minhduc
 *
 */
//@Deployable(scriptName = "jsa.np.f5reader2", scriptDesc = "Extract nanopore data (fastq/fasta and native data) from h5 files")
public class NanoporeReader// implements Closeable
{
	public static void main(String[] args) throws OutOfMemoryError, Exception {
		/*********************** Setting up script ****************************/
		Deployable annotation = NanoporeReader.class.getAnnotation(Deployable.class);
		CommandLine cmdLine = new CommandLine("\nUsage: "
				+ annotation.scriptName() + " [options] f1.fast5 f2.fast5 ...",
				annotation.scriptDesc());

		cmdLine.addString("output", "-",
				"Name of the output file, -  for stdout");
		cmdLine.addString("type", "fastq", 
				"Type of data to be extracted:" 
						+ "\nfastq: sequence read in fastq format"
						+ "\nevents: get events"
						+ "\nmodels: get models"
						+ "\nkeys:   list all keys"
				);

		cmdLine.addInt("minLength", 0, 
				"Minimum sequence length");

		cmdLine.addBoolean("stats", false, "Compute statistics of reads");
		cmdLine.addBoolean("number", false, "Add a unique number to read name");
		cmdLine.addString("f5list",null, "File containing list of fast5 files, one file per line");			

		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String type   = cmdLine.getStringVal("type");
		String output = cmdLine.getStringVal("output");
		String f5list = cmdLine.getStringVal("f5list");
		int minLength  = cmdLine.getIntVal("minLength");
		boolean stats  = cmdLine.getBooleanVal("stats");
		boolean number  = cmdLine.getBooleanVal("number");

		ArrayList<String> fileList = new ArrayList<String>();
		if (f5list != null){
			BufferedReader bf = SequenceReader.openFile(f5list);
			String line;
			while ((line = bf.readLine()) != null){
				fileList.add(line.trim());
			}
			bf.close();
		}

		for (int i = 0; i < args.length; i++){
			fileList.add(args[i]);
		}

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);
		if (type.equals("fastq"))
			readFastq(fileList, minLength, sos, stats,number);
		else if (type.equals("events"))
			readEvents(fileList, sos, stats);
		else if (type.equals("models"))
			readModels(fileList, sos, stats);
		else if (type.equals("keys"))
			readKeys(fileList, sos, stats);

		sos.close();
		//int maxLength = 0, minLength = Integer.MAX_VALUE;
	}//main


	public static void readKeys(ArrayList<String> fileList, SequenceOutputStream sos, boolean stats){
		for (String fileName:fileList){
			Logging.info("Open " + fileName);
			try{				
				NanoporeReader reader = new NanoporeReader(fileName);
				reader.readKeys();
				reader.close();
				sos.print("Keys in " + fileName +"\n");
				sos.print(reader.keyList);
			}catch (Exception e){
				Logging.error("Problem with reading " + fileName + ":" + e.getMessage());					
			}
		}//for		

	}

	public static void readEvents(ArrayList<String> fileList, SequenceOutputStream sos, boolean stats){
		for (String fileName:fileList){
			Logging.info("Open " + fileName);
			try{				
				NanoporeReader reader = new NanoporeReader(fileName);
				reader.readData();
				reader.close();
				if (reader.events != null){
					int maxIndx = 0, minIndx = 0;
					sos.print("Detected Events:\n");
					for (int i = 0; i < reader.events.mean.length;i++){
						sos.print(i);
						sos.print('\t');
						sos.print(reader.events.mean[i]);
						sos.print('\t');
						sos.print(reader.events.stdv[i]);
						sos.print('\t');
						sos.print(reader.events.length[i]);
						sos.print('\t');
						sos.print(reader.events.start[i]);
						sos.print('\n');
						if (stats){
							if (reader.events.mean[i] < reader.events.mean[minIndx])
								minIndx = i;

							if (reader.events.mean[i] > reader.events.mean[maxIndx])
								maxIndx = i;								
						}
					}
					if (stats){
						Logging.info("Min Event = " + reader.events.mean[minIndx] + " at " + minIndx);
						Logging.info("Max Event = " + reader.events.mean[maxIndx] + " at " + maxIndx);
					}
				}
				if (reader.bcTempEvents != null){
					int maxIndx = 0, minIndx = 0;
					sos.print("Template Events:\n");
					for (int i = 0; i < reader.bcTempEvents.mean.length;i++){
						sos.print(i);
						sos.print('\t');
						sos.print(reader.bcTempEvents.mean[i]);
						sos.print('\t');
						sos.print(reader.bcTempEvents.stdv[i]);
						sos.print('\t');
						sos.print(reader.bcTempEvents.length[i]);
						sos.print('\t');
						sos.print(reader.bcTempEvents.start[i]);
						sos.print('\n');

						if (stats){
							if (reader.bcTempEvents.mean[i] < reader.bcTempEvents.mean[minIndx])
								minIndx = i;

							if (reader.bcTempEvents.mean[i] > reader.bcTempEvents.mean[maxIndx])
								maxIndx = i;								
						}
					}
					if (stats){
						Logging.info("Min Temp = " + reader.bcTempEvents.mean[minIndx] + " at " + minIndx);
						Logging.info("Max Temp = " + reader.bcTempEvents.mean[maxIndx] + " at " + maxIndx);
					}
				}
				if (reader.bcCompEvents != null){
					int maxIndx = 0, minIndx = 0;
					sos.print("Complement Events:\n");
					for (int i = 0; i < reader.bcCompEvents.mean.length;i++){
						sos.print(i);
						sos.print('\t');
						sos.print(reader.bcCompEvents.mean[i]);
						sos.print('\t');
						sos.print(reader.bcCompEvents.stdv[i]);
						sos.print('\t');
						sos.print(reader.bcCompEvents.length[i]);
						sos.print('\t');
						sos.print(reader.bcCompEvents.start[i]);
						sos.print('\n');

						if (stats){
							if (reader.bcCompEvents.mean[i] < reader.bcCompEvents.mean[minIndx])
								minIndx = i;

							if (reader.bcCompEvents.mean[i] > reader.bcCompEvents.mean[maxIndx])
								maxIndx = i;								
						}				
					}
					if (stats){
						Logging.info("Min Comp = " + reader.bcCompEvents.mean[minIndx] + " at " + minIndx);
						Logging.info("Max Comp = " + reader.bcCompEvents.mean[maxIndx] + " at " + maxIndx);
					}
				}



			}catch (Exception e){
				Logging.error("Problem with reading " + fileName + ":" + e.getMessage());					
			}
		}//for		

	}


	public static void readModels(ArrayList<String> fileList, SequenceOutputStream sos, boolean stats){
		for (String fileName:fileList){
			Logging.info("Open " + fileName);
			try{				
				NanoporeReader reader = new NanoporeReader(fileName);
				reader.readData();
				reader.close();

				if (reader.bcTempModel != null){
					int maxIndx = 0, minIndx = 0;
					sos.print("Template model:" + fileName +"\n");
					for (int i = 0; i < reader.bcTempModel.levelMean.length;i++){
						sos.print(reader.bcTempModel.kmer[i]);
						sos.print('\t');
						sos.print(reader.bcTempModel.levelMean[i]);
						sos.print('\t');
						sos.print(reader.bcTempModel.levelStdv[i]);
						sos.print('\t');
						sos.print(reader.bcTempModel.sdMean[i]);
						sos.print('\t');
						sos.print(reader.bcTempModel.sdStdv[i]);
						sos.print('\t');	
						//sos.print(reader.bcTempModel.weigth[i]);
						//sos.print('\n');					
						if (stats){
							if (reader.bcTempModel.levelMean[i] < reader.bcTempModel.levelMean[minIndx])
								minIndx = i;

							if (reader.bcTempModel.levelMean[i] > reader.bcTempModel.levelMean[maxIndx])
								maxIndx = i;								
						}
					}
					if (stats){
						Logging.info("Min Event = " + reader.bcTempModel.levelMean[minIndx] + " at " + minIndx + "(" + reader.bcTempModel.kmer[minIndx] + ")");
						Logging.info("Max Event = " + reader.bcTempModel.levelMean[maxIndx] + " at " + maxIndx + "(" + reader.bcTempModel.kmer[maxIndx] + ")");
					}
				}

				if (reader.bcCompModel != null){
					int maxIndx = 0, minIndx = 0;
					sos.print("Complement model:\n");
					for (int i = 0; i < reader.bcCompModel.levelMean.length;i++){
						sos.print(reader.bcCompModel.kmer[i]);
						sos.print('\t');
						sos.print(reader.bcCompModel.levelMean[i]);
						sos.print('\t');
						sos.print(reader.bcCompModel.levelStdv[i]);
						sos.print('\t');
						sos.print(reader.bcCompModel.sdMean[i]);
						sos.print('\t');
						sos.print(reader.bcCompModel.sdStdv[i]);
						sos.print('\t');	
						//sos.print(reader.bcCompModel.weigth[i]);
						//sos.print('\n');
						if (stats){
							if (reader.bcCompModel.levelMean[i] < reader.bcCompModel.levelMean[minIndx])
								minIndx = i;

							if (reader.bcCompModel.levelMean[i] > reader.bcCompModel.levelMean[maxIndx])
								maxIndx = i;								
						}
					}
					if (stats){
						Logging.info("Min Event = " + reader.bcCompModel.levelMean[minIndx] + " at " + minIndx + "(" + reader.bcCompModel.kmer[minIndx] + ")");
						Logging.info("Max Event = " + reader.bcCompModel.levelMean[maxIndx] + " at " + maxIndx + "(" + reader.bcCompModel.kmer[maxIndx] + ")");
					}
				}
			}catch (Exception e){
				Logging.error("Problem with reading " + fileName + ":" + e.getMessage());					
			}
		}//for		

	}



	/**
	 * Read read sequence from a list of fast5 files.
	 * @param fileList
	 * @param sos : output stream
	 * @param stats: print out statistics
	 */
	public static void readFastq(ArrayList<String> fileList, int minLength, SequenceOutputStream sos, boolean stats,boolean number){
		int tempCount = 0, compCount = 0, twoDCount = 0;
		int fileNumber = 0;
		IntArray lengths = new IntArray();
		{
			for (String fileName:fileList){
				Logging.info("Open " + fileName);
				try{					
					NanoporeReader reader = new NanoporeReader(fileName);
					reader.readFastq();
					reader.close();


					//Get time & date
					String log = reader.getLog();					
					if (log != null){
						String [] toks = log.split("\n");
						if (toks.length > 0)
							toks = toks[toks.length - 1].split(",");

						log = toks[0];
					}else
						log = "";

					FastqSequence fq;

					fq = reader.getSeq2D();
					if (fq != null && fq.length() >= minLength){
						fq.setName((number?(fileNumber *3) + "_":"") + fq.getName() + " " + log);
						fq.print(sos);
						if (stats){						
							lengths.add(fq.length());	
							twoDCount ++;
						}
					}

					fq = reader.getSeqTemplate();
					if (fq != null && fq.length() >= minLength){
						fq.setName((number?(fileNumber *3 + 1) + "_":"") + fq.getName() + " " + log);
						fq.print(sos);
						if (stats){						
							lengths.add(fq.length());	
							tempCount ++;
						}
					}

					fq = reader.getSeqComplement();
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
				}

			}//for			
		}//if - else


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

			Logging.info("Open " + fileNumber + " files from " + fileList.size());
			Logging.info("Read count = " + ls.length + "(" + tempCount + " temppate, " + compCount + " complements and " + twoDCount +"  2D)");
			Logging.info("Base count = " + baseCount);		
			Logging.info("Longest read = " + ls[ls.length - 1] + ", shortest read = " + ls[0]);
			Logging.info("Average read length = " + mean);
			Logging.info("Median read length = " + median);
			Logging.info("Quantile first = " + ls[quantile1st] + " second = " + ls[quantile2nd] + " third = " + ls[quantile3rd]);
		}
	}


	String log = null;
	String keyList = "";//List of key 

	BaseCallAlignment2D bcAlignment2D = null;
	BaseCallAlignmentHairpin bcAlignmentHairpin = null;
	BaseCallModel bcCompModel = null, bcTempModel = null;

	DetectedEvents events;
	BaseCallEvents bcCompEvents = null, bcTempEvents = null;

	FastqSequence seqTemplate = null, seqComplement = null, seq2D = null;


	private FileFormat f5File;

	/**
	 * Open a fast5 file before reading anything from it.
	 * 
	 * The file should be closed before gabbage collected.
	 * 
	 * @param fileName
	 * @throws OutOfMemoryError
	 * @throws Exception
	 */
	public NanoporeReader (String fileName) throws JapsaException, OutOfMemoryError, Exception{		
		FileFormat fileFormat = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);
		
		if (fileFormat == null){
			throw new JapsaException("Cannot read HDF5 file, possily because JHI5 is not installed or configured. Please refer to npReader installation guide or contact the deverlopers.");
		}

		//Logging.info("Open " + fileName);
		f5File = fileFormat.createInstance(fileName, FileFormat.READ);
		if (f5File == null) 
			throw new RuntimeException("Unable to open file " + fileName);

		f5File.open();				
	}	


	public void close() throws Exception{		
		f5File.close();
	}

	public void readFastq() throws OutOfMemoryError, Exception{
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();	
		readData(root, false);
	}

	public void readData() throws OutOfMemoryError, Exception{
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readData(root, true);
	}
	public void readKeys() throws OutOfMemoryError, Exception{
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readMembers(root);
	}


	/**
	 * Get base call events for complement strand
	 * @return the bcCompEvents
	 */
	public BaseCallEvents getBcCompEvents() {
		return bcCompEvents;
	}


	/**
	 * Get base call events for template strand
	 * @return the bcTempEvents
	 */
	public BaseCallEvents getBcTempEvents() {
		return bcTempEvents;
	}

	/**
	 * Get 2D alignment
	 * @return the bcAlignment2D
	 */
	public BaseCallAlignment2D getBcAlignment2D() {
		return bcAlignment2D;
	}

	/**
	 * Get hairpin alignment
	 * @return the bcAlignmentHairpin
	 */
	public BaseCallAlignmentHairpin getBcAlignmentHairpin() {
		return bcAlignmentHairpin;
	}

	/**
	 * Get the model for base call of the complement
	 * @return the bcCompModel
	 */
	public BaseCallModel getBcCompModel() {
		return bcCompModel;
	}

	/**
	 * Get the model for base call of the template
	 * @return the bcTempModel
	 */
	public BaseCallModel getBcTempModel() {
		return bcTempModel;
	}

	/**
	 * Get the events from the pore
	 * @return the events
	 */
	public DetectedEvents getEvents() {
		return events;
	}

	/**
	 * @return the seqTemplate
	 */
	public FastqSequence getSeqTemplate() {
		return seqTemplate;
	}

	public String getLog() {
		return log;
	}

	/**
	 * @return the seqComplement
	 */
	public FastqSequence getSeqComplement() {
		return seqComplement;
	}

	/**
	 * @return the seq2D
	 */
	public FastqSequence getSeq2D() {
		return seq2D;
	}

	/**
	 * Recursively print its member names and types.
	 * 
	 * @throws OutOfMemoryError
	 * @throws Exception
	 */
	private void readMembers(Group g) throws OutOfMemoryError, Exception{

		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		


		for (HObject member:members) {
			String fullName = member.getFullName();
			if (member instanceof Group) {
				this.keyList += "Group        : " + fullName + "\n"; 
				readMembers((Group) member);
			}else if (member instanceof H5CompoundDS){

				//Logging.info(member.getClass() +" ");				
				Object dat =   ((H5CompoundDS) member).getData();
				if (dat != null){
					this.keyList += "H5CompoundDS : " + fullName + "=" + dat.getClass() +"\n";
				}else
					this.keyList += "H5CompoundDS : " + fullName + "=null\n";
			}else if (member instanceof H5ScalarDS){
				Object  dat = ((H5ScalarDS) member).getData();
				if (dat != null){
					this.keyList += "H5ScalarDS   : " + fullName + "=" + dat.getClass() +"\n";
				}else
					this.keyList += "H5ScalarDS   : " + fullName + "=null\n";
			}
		}
	}

	String expStart = null;
	/**
	 * Recursively print a group and its members. Fastq data are read.If all 
	 * flag is turned on, this method will also reads all events and model data.
	 * @throws OutOfMemoryError 
	 * 
	 * @throws Exception
	 */
	private void readData(Group g, boolean all) throws OutOfMemoryError, Exception{

		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		

		for (HObject member:members) {
			//System.out.println(indent + member + " " + member.getPath() + " " + member.getClass());
			//System.out.println(indent + member + " " + member.getPath() + " " + member.getClass());
			String f = member.getFullName();
			if (f.contains("tracking_id")){				
				@SuppressWarnings("unchecked")
				List<ncsa.hdf.object.Attribute> aL =  (List<ncsa.hdf.object.Attribute>) member.getMetadata();
				for (ncsa.hdf.object.Attribute att:aL){
					if (att.getName().equals("exp_start_time")){
						expStart = ((String[])att.getValue())[0];		
					}
				}
			}			

			if (member instanceof Group) {
				readData((Group) member, all);
			}else if (all && member instanceof H5CompoundDS){ 
				String fullName = member.getFullName();

				//Logging.info(member.getClass() +" ");				
				@SuppressWarnings("unchecked")
				List<Object> dat = (List<Object>)  (((H5CompoundDS) member).getData());
				if (dat != null){
					/********************************************************/
					if (fullName.endsWith("BaseCalled_2D/Alignment")){
						Logging.info("Read " + fullName);
						bcAlignment2D = new BaseCallAlignment2D();
						bcAlignment2D.template = (long[]) dat.get(0);
						bcAlignment2D.complement = (long[]) dat.get(1);
						bcAlignment2D.kmer = (String[]) dat.get(2);
					}else if (fullName.endsWith("BaseCalled_complement/Events")){
						Logging.info("Read " + fullName);
						bcCompEvents = new BaseCallEvents();
						bcCompEvents.mean  =  (double[]) dat.get(0);
						bcCompEvents.start =  (double[]) dat.get(1);
						bcCompEvents.stdv  =  (double[]) dat.get(2);
						bcCompEvents.length =  (double[]) dat.get(3);
						bcCompEvents.modelState =  (String[]) dat.get(4);						
						bcCompEvents.modelLevel = (double[]) dat.get(5);
						bcCompEvents.move = (long[]) dat.get(6);						
						bcCompEvents.pModelState = (double[]) dat.get(7);						
						bcCompEvents.mpState = (String[]) dat.get(8);						
						bcCompEvents.pMpState = (double[]) dat.get(9);

						bcCompEvents.pA = (double[]) dat.get(10);
						bcCompEvents.pC = (double[]) dat.get(11);
						bcCompEvents.pG = (double[]) dat.get(12);
						bcCompEvents.pT = (double[]) dat.get(13);
						//bcCompEvents.rawIndex = (long[]) dat.get(14);
					}else if (fullName.endsWith("BaseCalled_template/Events")){
						Logging.info("Read " + fullName);
						bcTempEvents = new BaseCallEvents();
						bcTempEvents.mean  =  (double[]) dat.get(0);
						bcTempEvents.start =  (double[]) dat.get(1);
						bcTempEvents.stdv  =  (double[]) dat.get(2);
						bcTempEvents.length =  (double[]) dat.get(3);
						bcTempEvents.modelState =  (String[]) dat.get(4);						
						bcTempEvents.modelLevel = (double[]) dat.get(5);
						bcTempEvents.move = (long[]) dat.get(6);						
						bcTempEvents.pModelState = (double[]) dat.get(7);						
						bcTempEvents.mpState = (String[]) dat.get(8);						
						bcTempEvents.pMpState = (double[]) dat.get(9);

						bcTempEvents.pA = (double[]) dat.get(10);
						bcTempEvents.pC = (double[]) dat.get(11);
						bcTempEvents.pG = (double[]) dat.get(12);
						bcTempEvents.pT = (double[]) dat.get(13);
						//bcTempEvents.rawIndex = (long[]) dat.get(14);
					}else if (fullName.endsWith("BaseCalled_complement/Model")){
						Logging.info("Read " + fullName);
						bcCompModel = new BaseCallModel();
						bcCompModel.kmer = (String[]) dat.get(0);
						//bcCompModel.variant = (double[]) dat.get(1);
						bcCompModel.levelMean = (double[]) dat.get(2);
						bcCompModel.levelStdv = (double[]) dat.get(3);
						bcCompModel.sdMean = (double[]) dat.get(4);
						bcCompModel.sdStdv = (double[]) dat.get(5);						
						//bcCompModel.weigth = (double[]) dat.get(6);
					}else if (fullName.endsWith("BaseCalled_template/Model")){
						Logging.info("Read " + fullName);
						bcTempModel = new BaseCallModel();
						bcTempModel.kmer = (String[]) dat.get(0);
						//bcTempModel.variant = (double[]) dat.get(1);
						bcTempModel.levelMean = (double[]) dat.get(2);
						bcTempModel.levelStdv = (double[]) dat.get(3);
						bcTempModel.sdMean = (double[]) dat.get(4);
						bcTempModel.sdStdv = (double[]) dat.get(5);						
						//bcTempModel.weigth = (double[]) dat.get(6);
					}else if (fullName.startsWith("/Analyses/EventDetection_000/Reads/") && fullName.endsWith("Events") ){
						Logging.info("Read " + fullName);
						events = new DetectedEvents();						
						events.mean =  (double[]) dat.get(0);						
						events.stdv =  (double[]) dat.get(1);
						events.start =  (long[]) dat.get(2);
						events.length =  (long[]) dat.get(3);
					}else if (fullName.endsWith("HairpinAlign/Alignment")){
						Logging.info("Read " + fullName);
						bcAlignmentHairpin = new BaseCallAlignmentHairpin();
						bcAlignmentHairpin.template =  (long[]) dat.get(0);
						bcAlignmentHairpin.complement =  (long[]) dat.get(1);						
					}
					/********************************************************/
				}
			}else if (member instanceof H5ScalarDS){
				String fullName = member.getFullName(); 
				if (fullName.endsWith("Fastq")){
					Object  data = ((H5ScalarDS) member).getData();
					if (data != null){
						Logging.info("Read " + fullName);
						String [] toks = ((String[]) data)[0].split("\n");						
						if  (fullName.contains("BaseCalled_2D")){
							//toks[0] = toks[0].substring(1) + "_twodimentional#" + f5File.getName().replace("imb13_010577_lt", "imb13-010577-lt") + " length=" + toks[1].length() ;							 
							toks[0] = toks[0].substring(1) + "_twodimentional" + " length=" + toks[1].length() ;
							this.seq2D =  new FastqSequence(DNA.DNA16(), toks);                		
						}else if (fullName.contains("BaseCalled_complement")){
							//toks[0] = toks[0].substring(1) + "_complement#" + f5File.getName().replace("imb13_010577_lt", "imb13-010577-lt") + " length=" + toks[1].length() ;							
							toks[0] = toks[0].substring(1) + "_complement" + " length=" + toks[1].length() ;
							this.seqComplement =  new FastqSequence(DNA.DNA16(), toks);							
						}else if (fullName.contains("BaseCalled_template")){
							//toks[0] = toks[0].substring(1) + "_template#" + f5File.getName().replace("imb13_010577_lt", "imb13-010577-lt") + " length=" + toks[1].length() ;
							toks[0] = toks[0].substring(1) + "_template" + " length=" + toks[1].length() ;
							this.seqTemplate =  new FastqSequence(DNA.DNA16(), toks);
						}
					}
				}else if (fullName.endsWith("Log")){
					Logging.info("Read " + fullName);
					Object  data = ((H5ScalarDS) member).getData();
					if (data != null){
						log =  ((String[]) data)[0];
						//System.out.println("\n\n" + log + "\n\n");
					}					
				}

			}
		}
	}

	public static class BaseCallModel{		
		String [] kmer;
		double[] variant;
		double[] levelMean, levelStdv, sdMean, sdStdv;//, weigth;

	}

	public static class BaseCallEvents{
		int dim;
		double [] mean, start, stdv, length, modelLevel, pModelState, pMpState, pA, pC, pG, pT;
		long []move;//, rawIndex;
		String [] modelState, mpState;
	}

	public static class BaseCallAlignment2D{
		int dim;
		long [] template, complement;
		String [] kmer;
	}

	public static class BaseCallAlignmentHairpin{
		int dim;
		long [] template, complement;
	}

	public static class DetectedEvents{
		int dim;
		double [] mean, stdv;
		long [] start;
		long [] length;

		public double [] getMean(){
			return mean;
		}

		/**
		 * @return the stdv
		 */
		public double[] getStdv() {
			return stdv;
		}

		/**
		 * @return the start
		 */
		public long[] getStart() {
			return start;
		}

		/**
		 * @return the length
		 */
		public long[] getLength() {
			return length;
		}

	}

}
