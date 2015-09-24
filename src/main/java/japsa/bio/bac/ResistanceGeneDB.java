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
 * 7 Sep 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsa.bio.bac;

import japsa.seq.SequenceReader;
import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Database for resistance gene identification
 * @author minhduc
 *
 */
public class ResistanceGeneDB {

	public static final String SEPARATOR = "\t"; 
	public static final String COMMENT = "#";

	String dbPath;//Act line the ID of the database	
	String fasPath;

	HashMap<String, String> gene2Res;

	public ResistanceGeneDB(String path) throws IOException {
		dbPath = path;
		fasPath = path + ".fas";
		
		//Read 		
		gene2Res = new HashMap<String, String>();
		BufferedReader br = SequenceReader.openFile(dbPath + ".map");
		String line;
		//int lineNo = 0;
		while ((line =  br.readLine()) != null){
			//lineNo ++;
			if (line.startsWith(COMMENT))
				continue;
			String [] toks = line.split(SEPARATOR);
			if (toks.length >= 2){
				gene2Res.put(toks[0], toks[1]);	
			}			

		}
		Logging.info("Read in " + gene2Res.size() + " genes");
		br.close();		
	}	

	public String getRes(String geneID){
		return gene2Res.get(geneID);
	}

	public String getSequenceFile(){
		return fasPath;
	}
}
