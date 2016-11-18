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
 * 23/09/20116 - Minh Duc Cao: Resigned of the reader class                                        
 *  
 ****************************************************************************/

package japsa.seq.nanopore;

import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5ScalarDS;
import japsa.seq.Alphabet.DNA;

import java.util.ArrayList;

import japsa.seq.Alphabet;
import japsa.seq.FastqSequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.JapsaException;


/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format.
 *  
 *  
 * @author minhduc
 */
public class Fast5NPReader{	
	protected FileFormat f5File;
	//	FastqSequence seqTemplate = null, seqComplement = null, seq2D = null;

	ArrayList<BaseCalledFastq> seqList = null;

	/**
	 * Open a fast5 file before reading anything from it.
	 * 
	 * The file should be closed before gabbage collected.
	 * 
	 * @param fileName
	 * @throws OutOfMemoryError
	 * @throws Exception
	 */
	public Fast5NPReader (String fileName) throws JapsaException, OutOfMemoryError, Exception{		
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
		if (seqList !=null) return;
		seqList = new ArrayList<BaseCalledFastq>();
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readFastq(root);
	}

	public ArrayList<BaseCalledFastq> getFastqList(){
		return seqList;
	}

	public void readAllFastq(SequenceOutputStream sos) throws OutOfMemoryError, Exception{
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readAllFastq(root, sos);
	}


	///**
	// * @return the seqTemplate
	// */
	//public FastqSequence getSeqTemplate() {
	//	return seqTemplate;
	//}


	///**
	// * @return the seqComplement
	// */
	//public FastqSequence getSeqComplement() {
	//	return seqComplement;
	//}

	///**
	// * @return the seq2D
	// */
	//public FastqSequence getSeq2D() {
	//	return seq2D;
	//}

	private void readAllFastq(Group g, SequenceOutputStream out) throws OutOfMemoryError, Exception{
		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		

		for (HObject member:members) {
			if (member instanceof Group) {
				readAllFastq((Group) member, out);			
			}else if (member instanceof H5ScalarDS){
				String fullName = member.getFullName(); 
				if (fullName.endsWith("Fastq")){
					Object  data = ((H5ScalarDS) member).getData();
					if (data != null){
						//Logging.info(fullName);
						//out.print(((String[]) data)[0]);
						//out.println();						
						//Logging.info("Read " + fullName);

						String [] toks = ((String[]) data)[0].split("\n",2);						
						if  (fullName.contains("BaseCalled_2D")){							
							out.print(toks[0] + "_twodimentional path="  + fullName);							                		
						}else if (fullName.contains("BaseCalled_complement")){
							out.print(toks[0] + "_complement path="  + fullName);							
						}else if (fullName.contains("BaseCalled_template")){
							out.print(toks[0] + "_template path="  + fullName);
						}else
							out.print(toks[0] + "_unknown path="  + fullName);

						out.print('\n');
						out.print(toks[1]);
						out.print('\n');						
					}
				}
			}
		}	
	}

	/**
	 * Recursively print a group and its members. Fastq data are read.If all 
	 * flag is turned on, this method will also reads all events and model data.
	 * @throws OutOfMemoryError 
	 * 
	 * @throws Exception
	 */
	private void readFastq(Group g) throws OutOfMemoryError, Exception{
		if (g == null) return;
		java.util.List<HObject> members = g.getMemberList();		

		for (HObject member:members) {
			if (member instanceof Group) {
				readFastq((Group) member);			
			}else if (member instanceof H5ScalarDS){
				String fullName = member.getFullName(); 
				if (fullName.endsWith("Fastq")){
					Object  data = ((H5ScalarDS) member).getData();
					String group = fullName.split("/")[2];					
					if (data != null){
						//Logging.info("Read " + fullName);
						String [] toks = ((String[]) data)[0].split("\n");						
						if  (fullName.contains("BaseCalled_2D")){
							toks[0] = toks[0].substring(1) + "_twodimentional" + " length=" + toks[1].length() + " group=" + group;
							seqList.add(new BaseCalledFastq(DNA.DNA16(), toks, BaseCalledFastq.TWODIM));
						}else if (fullName.contains("BaseCalled_complement")){
							toks[0] = toks[0].substring(1) + "_complement" + " length=" + toks[1].length() + " group=" + group ;
							seqList.add(new BaseCalledFastq(DNA.DNA16(), toks, BaseCalledFastq.COMPLEMENT));
						}else if (fullName.contains("BaseCalled_template")){
							toks[0] = toks[0].substring(1) + "_template" + " length=" + toks[1].length() + " group=" + group;
							seqList.add(new BaseCalledFastq(DNA.DNA16(), toks, BaseCalledFastq.TEMPLATE));
						}
					}
				}
			}
		}
	}

	public static class BaseCalledFastq extends FastqSequence{
		public static final int UNKNOWN = 4;
		public static final int TWODIM = 0;
		public static final int TEMPLATE = 1;
		public static final int COMPLEMENT = 2;		
		int myType = 4;

		public BaseCalledFastq(Alphabet alphabet, String [] toks, int type) {	
			super(alphabet, toks);		
			myType = type; 
		}

		public int type(){
			return myType;
		}

		public boolean isTwoDim(){
			return myType == TWODIM;
		}
		public boolean isTemplate(){
			return myType == TEMPLATE;
		}
		public boolean isComplement(){
			return myType == COMPLEMENT;
		}
	}
}
