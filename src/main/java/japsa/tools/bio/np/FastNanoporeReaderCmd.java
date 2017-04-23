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

/*****************************************************************************
 *                           Revision History                                
 * 7 Aug 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.tools.bio.np;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import japsa.seq.SequenceOutputStream;
import japsa.seq.nanopore.Fast5NPReader;
import japsa.util.CommandLine;
import japsa.util.JapsaException;
import japsa.util.deploy.Deployable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author minhduc
 *
 */
@Deployable(	
		scriptName = "jsa.np.fastnpreader", 
		scriptDesc = "Fast Extraction of Oxford Nanopore sequencing data in real-time"
		)
public class FastNanoporeReaderCmd extends CommandLine{
    private static final Logger LOG = LoggerFactory.getLogger(FastNanoporeReaderCmd.class);
	public FastNanoporeReaderCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());


		addString("folder", null,"The folder containing base-called reads");		
		addString("output", "-","Name of the output file, - for stdout");		

		addStdHelp();		
	}

	public static void main(String[] args) throws OutOfMemoryError, Exception {		
		CommandLine cmdLine = new FastNanoporeReaderCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/

		String output = cmdLine.getStringVal("output");		
		String folder = cmdLine.getStringVal("folder");		

		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(output);
		readFastq(folder, out);
		out.close();

	}//main


	public static void readFastq (String folderPath, SequenceOutputStream out) throws JapsaException, IOException{		
		File mainFolder = new File(folderPath);
		File passFolder = new File(folderPath + File.separatorChar + "pass");
		File failFolder = new File(folderPath + File.separatorChar + "fail");

		ArrayList<File> folders = new ArrayList<File>();
		folders.add(mainFolder);
		folders.add(passFolder);
		folders.add(failFolder);

		for (File folder:folders){
			File [] fileList = folder.listFiles();
			if (fileList!=null){
				for (File f:fileList){
					//directory
					if (!f.isFile())
						continue;//for						

					if (!f.getName().endsWith("fast5"))
						continue;//for
					String sPath = f.getAbsolutePath();
					try{
						Fast5NPReader npReader = new Fast5NPReader(sPath);
						npReader.readAllFastq(out);						
						npReader.close();
					}catch (JapsaException e){
						throw e;
					}catch (Exception e){
						LOG.error("Problem with reading " + sPath + ":" + e.getMessage());
						e.printStackTrace();					
					}
				}//for f
			}//if
		}//for folder	

	}
}

