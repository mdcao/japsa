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
 * 14 Aug 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.tools.bio.np;

import java.io.IOException;
import java.util.Date;

import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;


/**
 * @author minhduc
 *
 */
@Deployable(
	scriptName = "jsa.np.timeEmulate", 
	scriptDesc = "Regulate time"
	)
public class RegulateTimeCmd extends CommandLine {
	public RegulateTimeCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdInputFile();
		addStdOutputFile();
		
		addDouble("scale",1.0, "Scale");
		
		addStdHelp();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		RegulateTimeCmd cmdTool = new RegulateTimeCmd();
		args = cmdTool.stdParseLine(args);

		String input = cmdTool.getStringVal("input");
		String output = cmdTool.getStringVal("output");		
		double scale = cmdTool.getDoubleVal("scale");

		SequenceOutputStream sos = SequenceOutputStream.makeOutputStream(output);

		SequenceReader reader = SequenceReader.getReader(input);
		
		boolean isFastq = (reader instanceof FastqReader); 

		Sequence seq;

		long timeStart = System.currentTimeMillis();
		long reportTime = timeStart;
		
		String sortKeyOptionPrefix = "cTime=";
		int sortKeyOptionIndex = sortKeyOptionPrefix.length(); 
		long firstReadTime = 0; 
		int numRead = 0;
		long numBase = 0;
		
		while ((seq = reader.nextSequence(Alphabet.DNA()))!= null){			
			double cTime = 0;
			String [] toks = seq.getName().split(" ");
			try{
				for (int i = 0; i < toks.length;i++){
					if (toks[i].startsWith(sortKeyOptionPrefix)){
						cTime = Double.parseDouble(toks[i].substring(sortKeyOptionIndex));
						break;
					}	
				}
			}catch (Exception e){
				Logging.error(e.getMessage());
			}			
			if (cTime == 0){
				Logging.info("Not found timing for sequence " + seq.getName());
				continue;
			}
			if (firstReadTime == 0){
				firstReadTime = (long) cTime;
			}
			
			long reportTimeNow = System.currentTimeMillis();
			if (reportTimeNow - reportTime >= 60000){
				reportTime = reportTimeNow; 
				Logging.info(new Date(reportTime) + " : " + numRead + " reads " + numBase + " bases");
			}			
			
			cTime = 1000* (cTime - firstReadTime) / scale;//scale and convert to milisecond
						
			long timeNow = ((long) cTime)  - (System.currentTimeMillis() - timeStart);			
			if (timeNow > 0){
				sos.flush();
				try {
					Thread.sleep(timeNow);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			//time is up
			
			if (isFastq){
				FastqSequence fq = ((FastqSequence) seq);
				fq.print(sos);				
			}else				
				seq.writeFasta(sos);
			
			numRead ++;
			numBase += seq.length();			
			
		}//while
		sos.close();
	}
}
