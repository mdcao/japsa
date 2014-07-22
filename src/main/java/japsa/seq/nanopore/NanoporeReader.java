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

import java.util.Vector;

import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5CompoundDS;
import ncsa.hdf.object.h5.H5ScalarDS;
import japsa.seq.Alphabet.DNA;
import japsa.seq.FastqSequence;
import japsa.util.Logging;

/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format. 
 * @author minhduc
 *
 */
public class NanoporeReader {
	BaseCallEvents bcCompEvents = null, bcTempEvents = null;
	BaseCallAlignment2D bcAlignment2D = null;
	BaseCallAlignmentHairpin bcAlignmentHairpin = null;
	BaseCallModel bcCompModel = null, bcTempModel = null;
	DetectedEvents events;
	FastqSequence seqTemplate = null, seqComplement = null, seq2D = null;
	
	
	
	
	public NanoporeReader (String fileName) throws OutOfMemoryError, Exception{
		FileFormat fileFormat = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);

		FileFormat f5File = fileFormat.createInstance(fileName, FileFormat.READ);
		if (f5File == null) 
			throw new RuntimeException("Unable to open file " + fileName);

		f5File.open();

		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();

		readGroup(root);

		// close file resource
		f5File.close();		
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
	 * Recursively print a group and its members.
	 * @throws OutOfMemoryError 
	 * 
	 * @throws Exception
	 */
	private void readGroup(Group g) throws OutOfMemoryError, Exception{

		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		

		for (HObject member:members) {
			//System.out.println(indent + member + " " + member.getPath() + " " + member.getClass());
			if (member instanceof Group) {
				readGroup((Group) member);
			}else if (member instanceof H5CompoundDS){ 
				String fullName = member.getFullName();
				
				Vector<Object> dat = (Vector<Object>)  (((H5CompoundDS) member).getData());				
				
				if (dat != null){
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
						bcCompEvents.rawIndex = (long[]) dat.get(14);
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
						bcTempEvents.rawIndex = (long[]) dat.get(14);
					}else if (fullName.endsWith("BaseCalled_complement/Model")){
						Logging.info("Read " + fullName);
						bcCompModel = new BaseCallModel();
						bcCompModel.kmer = (String[]) dat.get(0);
						bcCompModel.variant = (long[]) dat.get(1);
						bcCompModel.levelMean = (double[]) dat.get(2);
						bcCompModel.levelStdv = (double[]) dat.get(3);
						bcCompModel.sdMean = (double[]) dat.get(4);
						bcCompModel.sdStdv = (double[]) dat.get(5);						
						bcCompModel.weigth = (double[]) dat.get(6);
					}else if (fullName.endsWith("BaseCalled_template/Model")){
						Logging.info("Read " + fullName);
						bcTempModel = new BaseCallModel();
						bcTempModel.kmer = (String[]) dat.get(0);
						bcTempModel.variant = (long[]) dat.get(1);
						bcTempModel.levelMean = (double[]) dat.get(2);
						bcTempModel.levelStdv = (double[]) dat.get(3);
						bcTempModel.sdMean = (double[]) dat.get(4);
						bcTempModel.sdStdv = (double[]) dat.get(5);						
						bcTempModel.weigth = (double[]) dat.get(6);
					}else if (fullName.startsWith("/Analyses/EventDetection_000/Reads/") && fullName.endsWith("Events") ){
						Logging.info("Read " + fullName);
						events = new DetectedEvents();						
						events.mean =  (double[]) dat.get(0);						
						events.stdv =  (double[]) dat.get(1);
						events.start =  (double[]) dat.get(2);
						events.length =  (double[]) dat.get(3);
					}else if (fullName.endsWith("HairpinAlign/Alignment")){
						Logging.info("Read " + fullName);
						bcAlignmentHairpin = new BaseCallAlignmentHairpin();
						bcAlignmentHairpin.template =  (long[]) dat.get(0);
						bcAlignmentHairpin.complement =  (long[]) dat.get(1);						
					}
				}

			}else if (member instanceof H5ScalarDS){
				String fullName = member.getFullName(); 
				if (fullName.endsWith("Fastq")){
					Object  data = ((H5ScalarDS) member).getData();
					if (data != null){
						Logging.info("Read " + fullName);
						String [] toks = ((String[]) data)[0].split("\n");						
						if  (fullName.contains("BaseCalled_2D")){
							toks[0] = toks[0].substring(1) + "_twodimentional length=" + toks[1].length() ;							 
							this.seq2D =  new FastqSequence(DNA.DNA16(), toks);;                		
						}else if (fullName.contains("BaseCalled_complement")){
							toks[0] = toks[0].substring(1) + "_complement length=" + toks[1].length() ;							
							this.seqComplement =  new FastqSequence(DNA.DNA16(), toks);							
						}else if (fullName.contains("BaseCalled_template")){
							toks[0] = toks[0].substring(1) + "_template length=" + toks[1].length() ;
							this.seqTemplate =  new FastqSequence(DNA.DNA16(), toks);
						}
					}
				}

			}
		}
	}	

	public static class BaseCallModel{		
		String [] kmer;
		long[] variant;
		double[] levelMean, levelStdv, sdMean, sdStdv, weigth;

	}

	public static class BaseCallEvents{
		int dim;
		double [] mean, start, stdv, length, modelLevel, pModelState, pMpState, pA, pC, pG, pT;
		long []move, rawIndex;
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
		double [] mean, stdv, start, length;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
