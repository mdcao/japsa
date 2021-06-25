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

/*                           Revision History                                
 * 11/01/2012 - Minh Duc Cao: Revised 
 * 01/01/2013 - Minh Duc Cao, revised                                       
 ****************************************************************************/

package japsadev.tools;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import japsa.bio.tr.TandemRepeat;
import japsa.seq.XAFReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;


/**
 * @author Minh Duc Cao
 * 
 */
@Deployable(
		scriptName = "jsa.dev.vntrselect",
		scriptDesc = "Select vntr d"
		)
public class VNTRSelectCmd extends CommandLine{	
	public VNTRSelectCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addStdHelp();		
	} 

	public static void main(String[] args) throws IOException {
		/*********************** Setting up script ****************************/		
		CommandLine cmdLine = new VNTRSelectCmd();		
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		if (args.length != 3){
			System.err.println(" args = 3 ");
		}
		ArrayList<TandemRepeat> myStr = new ArrayList<TandemRepeat>(); 

		String strFile = args[0];
		String hgFile = args[1];
		String trdbFile = args[2];

		XAFReader xafReader = new XAFReader(strFile);

		while (xafReader.next() != null){	
			//_tIndex ++;
			TandemRepeat str = TandemRepeat.read(xafReader);
			myStr.add(str);
			//System.out.println(str.getID() + "  " + str.getChr() + " " + str.getStart() + " " + str.getEnd());			
		}
		xafReader.close();
		
		HashMap<String,String> map = new HashMap<String,String>();
		
				
		BufferedReader bf = new BufferedReader (new FileReader(hgFile));
		String line = "";

		while ( (line = bf.readLine())!= null){
			line = line.trim();
			String [] toks = line.split("\t");
			map.put("trid" + toks[0], line);
		}
		bf.close();

			
		

		bf = new BufferedReader (new FileReader(trdbFile));
		line = "";

		while ( (line = bf.readLine())!= null){
			if (line.startsWith(">")){
				line = line.trim().substring(1);
				String [] toks = line.split(":");
				int    myStart = Integer.parseInt(toks[2]);
				String myChrom = "chr"+toks[1];

				for (TandemRepeat tr:myStr){
					int start = tr.getStart();
					String chrom = tr.getChr();
					if (Math.abs(start - myStart) < 400 && chrom.equals(myChrom)){
						
						System.out.println(tr.getID() + "\t" + 
								tr.getChr() + "\t" + 
								tr.getStart() +"\t" + 
								tr.getEnd() +"\t" +
								tr.getPeriod() +"\t" +
								tr.getUnitNo() +"==" +								
								line +"==" +
								map.get(toks[0]));
					}
				}

			}
		}//while
		bf.close();

	}	
}