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
import japsa.seq.FastqSequence;
import japsa.util.JapsaException;
import japsa.util.Logging;


/**
 * Read nanopore data (read sequence, events, alignment, models etc) from a raw
 * (fast5) format.
 * Re-implemented from the previous to aim for faster, which static key
 *  
 * @author minhduc
 */
public class Fast5NPReader{	
	protected FileFormat f5File;
	FastqSequence seqTemplate = null, seqComplement = null, seq2D = null;

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
		Group root = (Group) ((javax.swing.tree.DefaultMutableTreeNode) f5File.getRootNode()).getUserObject();
		readFastq(root);
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
					if (data != null){
						Logging.info("Read " + fullName);
						String [] toks = ((String[]) data)[0].split("\n");						
						if  (fullName.contains("BaseCalled_2D")){
							toks[0] = toks[0].substring(1) + "_twodimentional" + " length=" + toks[1].length() ;
							this.seq2D =  new FastqSequence(DNA.DNA16(), toks);                		
						}else if (fullName.contains("BaseCalled_complement")){
							toks[0] = toks[0].substring(1) + "_complement" + " length=" + toks[1].length() ;
							this.seqComplement =  new FastqSequence(DNA.DNA16(), toks);							
						}else if (fullName.contains("BaseCalled_template")){
							toks[0] = toks[0].substring(1) + "_template" + " length=" + toks[1].length() ;
							this.seqTemplate =  new FastqSequence(DNA.DNA16(), toks);
						}
					}
				}
			}
		}
	}
/**************************************************************************************
	
	static String [] TEMPLATE_FQ_PATHS = {
			"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq",
			"/Analyses/Basecall_2D_000/BaseCalled_template/Fastq",
	};

	static String [] COMPLEMENT_FQ_PATHS = {
			"/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq",
			"/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq"
	};

	static String [] TWODIM_FQ_PATHS = {			
			"/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
	};
	

	public FastqSequence readTemplate() throws Exception{		
		for (String path:TEMPLATE_FQ_PATHS){
			Object object = f5File.get(path);
			if (object != null){
				Object  data = ((H5ScalarDS) object).getData();						
				String [] toks = ((String[]) data)[0].split("\n");		
				toks[0] = toks[0].substring(1) + "_template length=" + toks[1].length() ;			
				FastqSequence seq =  new FastqSequence(DNA.DNA16(), toks);
				return seq;
			}
		}
		return null;		
	}

	public FastqSequence readComplement() throws Exception{
		for (String path:COMPLEMENT_FQ_PATHS){Object object = f5File.get(path);
		if (object != null){
			Object  data = ((H5ScalarDS) object).getData();			
				String [] toks = ((String[]) data)[0].split("\n");		
				toks[0] = toks[0].substring(1) + "_complement length=" + toks[1].length() ;			
				FastqSequence seq =  new FastqSequence(DNA.DNA16(), toks);
				return seq;
			}
		}
		return null;		
	}

	public FastqSequence readTwoDim() throws Exception{
		for (String path:TWODIM_FQ_PATHS){Object object = f5File.get(path);
		if (object != null){
			Object  data = ((H5ScalarDS) object).getData();			
				String [] toks = ((String[]) data)[0].split("\n");		
				toks[0] = toks[0].substring(1) + "_twodimentional length=" + toks[1].length() ;			
				FastqSequence seq =  new FastqSequence(DNA.DNA16(), toks);
				return seq;
			}
		}
		return null;		
	}
	
/**************************************************************************************/	
	
}
