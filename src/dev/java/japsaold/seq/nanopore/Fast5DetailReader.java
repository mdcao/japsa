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

package japsaold.seq.nanopore;

import java.util.List;

import ncsa.hdf.object.Group;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5CompoundDS;
import ncsa.hdf.object.h5.H5ScalarDS;
import japsa.util.JapsaException;
import japsa.util.Logging;

/**
 * Read detail nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format.
 * Re-implemented from the previous to aim for faster, which static key
 *  
 * @author minhduc
 */
public class Fast5DetailReader extends Fast5NPReader{	
	static String RAW_PREFIX = "/Raw/Reads";
	static String CHANNEL_ID = "/UniqueGlobalKey/channel_id";
	static String TRACKING_ID = "/UniqueGlobalKey/tracking_id";		

	private double samplingRate = 0;//Default
	private  int channelNumber = 0;
	private  long startTime  = 0;	

	public Fast5DetailReader (String fileName) throws JapsaException, OutOfMemoryError, Exception{	
		super(fileName);
		readMetaData();
	}


	public double getSamplingRate() {
		return samplingRate;
	}


	public int getChannelNumber() {
		return channelNumber;
	}

	public long getStartTime() {
		return startTime;
	}


	/**
	 * Extract metadata about the read, including:
	 * - sampling rate
	 * - channel number
	 * @throws Exception
	 */

	private void readMetaData() throws Exception{
		HObject data = f5File.get(CHANNEL_ID);		
		@SuppressWarnings("unchecked")
		List<ncsa.hdf.object.Attribute> aL =  (List<ncsa.hdf.object.Attribute>) data.getMetadata();
		for (ncsa.hdf.object.Attribute att:aL){
			if (att.getName().equals("sampling_rate")){
				samplingRate = ((double[]) att.getValue())[0];				
			}else if (att.getName().equals("channel_number")){
				channelNumber = Integer.parseInt(((String[]) att.getValue())[0]);				
			}			
		}		
	}

	/**
	 * Extract raw event from the read. Start time is also extracted in the process
	 * @return
	 * @throws Exception
	 */

	public RawSignal getRawEvent() throws Exception{
		HObject data = f5File.get(RAW_PREFIX);
		//Logging.info("Read 1 " + (data == null));
		if (data !=null){
			Group group = (Group) data;
			group = (Group) group.getMemberList().get(0);
			@SuppressWarnings("unchecked")
			List<ncsa.hdf.object.Attribute> aL =  (List<ncsa.hdf.object.Attribute>) group.getMetadata();
			for (ncsa.hdf.object.Attribute att:aL){
				if (att.getName().equals("start_time")){
					startTime = ((long[]) att.getValue())[0];					
				}
			}
			H5ScalarDS  myDat =((H5ScalarDS) group.getMemberList().get(0));			
			//Logging.info("Read 2 " + (myDat == null));
			if (myDat != null){
				short [] rawEvent = (short[])myDat.getData();
				//Logging.info("Read 3 " + (rawEvent == null));
				RawSignal rawSignal  = new RawSignal(rawEvent);
				return rawSignal;
			}
		}		
		return null;
	}
	
	public void readData() throws OutOfMemoryError, Exception{
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readData(root);
	}
	
	
	
	/**
	 * Recursively print a group and its members. Fastq data are read.If all 
	 * flag is turned on, this method will also reads all events and model data.
	 * @throws OutOfMemoryError 
	 * 
	 * @throws Exception
	 */
	private void readData(Group g) throws OutOfMemoryError, Exception{

		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		

		for (HObject member:members) {
			//String f = member.getFullName();
			if (member instanceof Group) {
				readData((Group) member);
			}else if (member instanceof H5CompoundDS){ 
				String fullName = member.getFullName();

				//Logging.info(member.getClass() +" ");				
				@SuppressWarnings("unchecked")
				List<Object> dat = (List<Object>)  (((H5CompoundDS) member).getData());
				if (dat != null){
					/********************************************************/
					if (fullName.startsWith("/Analyses/EventDetection_000/Reads/") && fullName.endsWith("Events") ){
						Logging.info("Read " + fullName);
						detectedEvents = new DetectionEvents();						
						detectedEvents.start =  (long[]) dat.get(0);
						detectedEvents.length =  (long[]) dat.get(1);
						detectedEvents.mean =  (double[]) dat.get(2);						
						detectedEvents.stdv =  (double[]) dat.get(3);							
					}else if (fullName.endsWith("BaseCalled_template/Events")){
						Logging.info("Read " + fullName);
						bcTempEvents = new BaseCallEvents();
						bcTempEvents.mean  =  (double[]) dat.get(0);
						bcTempEvents.start =  (double[]) dat.get(1);
						bcTempEvents.stdv  =  (double[]) dat.get(2);
						bcTempEvents.length =  (double[]) dat.get(3);
						bcTempEvents.modelState =  (String[]) dat.get(4);
						bcTempEvents.move = (long[]) dat.get(5);
						bcTempEvents.weight = (float[]) dat.get(6);
						bcTempEvents.pModelState = (float[]) dat.get(7);						
						bcTempEvents.mpState = (String[]) dat.get(8);						
						bcTempEvents.pMpState = (float[]) dat.get(9);

						bcTempEvents.pA = (float[]) dat.get(10);
						bcTempEvents.pC = (float[]) dat.get(11);
						bcTempEvents.pG = (float[]) dat.get(12);
						bcTempEvents.pT = (float[]) dat.get(13);
					}else if (fullName.endsWith("BaseCalled_complement/Events")){
						Logging.info("Read " + fullName);
						bcCompEvents = new BaseCallEvents();
						bcCompEvents.mean  =  (double[]) dat.get(0);
						bcCompEvents.start =  (double[]) dat.get(1);
						bcCompEvents.stdv  =  (double[]) dat.get(2);
						bcCompEvents.length =  (double[]) dat.get(3);
						bcCompEvents.modelState =  (String[]) dat.get(4);
						bcCompEvents.move = (long[]) dat.get(5);
						bcCompEvents.weight = (float[]) dat.get(6);
						bcCompEvents.pModelState = (float[]) dat.get(7);						
						bcCompEvents.mpState = (String[]) dat.get(8);						
						bcCompEvents.pMpState = (float[]) dat.get(9);

						bcCompEvents.pA = (float[]) dat.get(10);
						bcCompEvents.pC = (float[]) dat.get(11);
						bcCompEvents.pG = (float[]) dat.get(12);
						bcCompEvents.pT = (float[]) dat.get(13);						
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
					}else if (fullName.endsWith("HairpinAlign/Alignment")){
						Logging.info("Read " + fullName);
						bcAlignmentHairpin = new BaseCallAlignmentHairpin();
						bcAlignmentHairpin.template =  (long[]) dat.get(0);
						bcAlignmentHairpin.complement =  (long[]) dat.get(1);						
					}else 	
						if (fullName.endsWith("BaseCalled_2D/Alignment")){
							Logging.info("Read " + fullName);
							bcAlignment2D = new BaseCallAlignment2D();
							bcAlignment2D.template = (long[]) dat.get(0);
							bcAlignment2D.complement = (long[]) dat.get(1);
							bcAlignment2D.kmer = (String[]) dat.get(2);
						}
					/********************************************************/
				}
			}
		}
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
	 * Get the events from the pore
	 * @return the events
	 */
	public DetectionEvents getEvents() {
		return detectedEvents;
	}

	

	/**
	 * Class represent RawSignal, which is an array of short
	 * @author minhduc
	 *
	 */
	public static class RawSignal{
		short [] signal;

		private RawSignal(short [] data){
			signal = data;
		}
		public short [] getSignal(){
			return signal;
		}
	}

	/********************************************************************************************************
H5CompoundDS : /Analyses/EventDetection_000/Reads/Read_12/Events=class java.util.Vector


			H5ScalarDS   : /Analyses/Calibration_Strand_000/Log=class [Ljava.lang.String;
			H5ScalarDS   : /Analyses/EventDetection_000/Log=class [Ljava.lang.String;
			H5CompoundDS : /Analyses/EventDetection_000/Reads/Read_5543/Events=class java.util.Vector
			H5ScalarDS   : /Raw/Reads/Read_5543/Signal=class [S

			                                                  /*******************************************************/
	BaseCallAlignment2D bcAlignment2D = null;
	BaseCallAlignmentHairpin bcAlignmentHairpin = null;
	BaseCallModel bcCompModel = null, bcTempModel = null;

	DetectionEvents detectedEvents;
	BaseCallEvents bcCompEvents = null, bcTempEvents = null;
	RawSignal rawSignal = null;
	
	double seqTime = 0;



	
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
	//String expStart = "";

	public static class BaseCallModel{		
		String [] kmer;
		double[] variant;
		double[] levelMean, levelStdv, sdMean, sdStdv;//, weigth;

	}

	public static class BaseCallEvents{
		int dim;
		double [] mean, start, stdv, length;
		float [] pA, pC, pG, pT;
		long []move;//, rawIndex;
		float [] pModelState;
		float [] pMpState;
		float [] weight;
		//long [] modelLevel;
		String [] modelState, mpState;

		public long [] getMove(){
			return move;
		}

		public double [] length(){
			return length;
		}

		public double [] mean(){
			return mean;
		}

		public double [] stdv(){
			return stdv;
		}		

		public float [] weight(){
			return weight;
		}

		public String [] modelState(){
			return modelState;
		}
	}

	public static class BaseCallAlignment2D{
		int dim;
		long [] template, complement;
		String [] kmer;

		public String[] getKmer(){
			return kmer;
		}
		public long[] getComplementKmer(){
			return complement;
		}
		public long[] getTemplateKmer(){
			return template;
		}
	}

	public static class BaseCallAlignmentHairpin{
		int dim;
		long [] template, complement;
	}

	public static class DetectionEvents{
		int dim;
		double [] stdv;
		double [] mean;
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
