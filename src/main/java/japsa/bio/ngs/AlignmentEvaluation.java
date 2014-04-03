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
 * 31/01/2013 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.bio.ngs;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import java.io.File;
import java.util.HashMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 *
 */
@Deployable(scriptName = "jsa.ngs.evalAlignment",
            scriptDesc = "Evaluate the accuracy of read aligmment")

public class AlignmentEvaluation {
	public static void main(String[] args) throws Exception {		
		/*********************** Setting up script ****************************/
		String scriptName = "jsa.ngs.evalAlignment";
		String desc = "Evaluate the accuracy of read aligmment\n";		
		CommandLine cmdLine = new CommandLine("\nUsage: " + scriptName + " [params]");
		/**********************************************************************/	

		cmdLine.addString("ansBam", null, "Name of the answer sam/bam file",true);
		cmdLine.addString("evaBam", null, "Name of the  sam/bam file to evaluate",true);		
		cmdLine.addInt("relax", 20, "Acceptable relax"); 
		cmdLine.addStdHelp();

		/**********************************************************************/
		args = cmdLine.parseLine(args);
		if (cmdLine.getBooleanVal("help")){
			System.out.println(desc + cmdLine.usage());			
			System.exit(0);
		}
		if (cmdLine.errors() != null) {
			System.err.println(cmdLine.errors() + cmdLine.usage());
			System.exit(-1);
		}	
		/**********************************************************************/
		String ansBam = cmdLine.getStringVal("ansBam");
		String evaBam = cmdLine.getStringVal("evaBam");
		int relax = cmdLine.getIntVal("relax");

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);		
		SAMFileReader ansReader = new  SAMFileReader(new File(ansBam));		
		SAMFileReader evaReader = new  SAMFileReader(new File(evaBam));

		SAMRecordIterator ansIter = ansReader.iterator(),
				evaIter = evaReader.iterator();

		HashMap<String, Boolean> ansHash = new HashMap<String, Boolean>();		

		//Timer timer = new Timer();
		int countIn = 0;
		while (ansIter.hasNext()){			
			SAMRecord ansRecord = ansIter.next();
			countIn ++;

			if (ansRecord.getReadPairedFlag()){
				if (ansRecord.getFirstOfPairFlag())
					ansHash.put(ansRecord.getReadName()+"/1", true);
				else 
					ansHash.put(ansRecord.getReadName()+"/2", true);
			}else
				ansHash.put(ansRecord.getReadName(), true);

			if (countIn % 10000000 == 0){
				//timer.mark("# Read in  up to " + countIn);
				Runtime.getRuntime().gc();
			}
		}		
		ansReader.close();
		evaReader.close();
		
		//timer.mark("# Read done");
		//System.out.println(ansHash.size() + " records ");

		int TP = 0, FP = 0, FN = ansHash.size(), dup = 0;
		int count = 0, countAns = FN;

		while (evaIter.hasNext()){
			count ++;
			SAMRecord evaRecord = evaIter.next();
			//check if the record is the right position:
			String readName = evaRecord.getReadName();
			String[] toks = readName.split("_");

			int pos = 0;

			//FIXME:comment out the below only for lobstr
			if (evaRecord.getReadPairedFlag()){
				if (evaRecord.getFirstOfPairFlag())
					pos = Integer.parseInt(toks[1]);
				else
					pos = Integer.parseInt(toks[2]);			
				//if (readName.endsWith("/2"))
				//	pos = Integer.parseInt(toks[2]);			
			}else
				pos = Integer.parseInt(toks[1]);


			if ( evaRecord.getReferenceName().equals(toks[0])
					&& Math.abs(evaRecord.getAlignmentStart() - pos) < relax){//mapped correctly
				TP ++;

				String hashKey = evaRecord.getReadName();

				//FIXME:comment out the below only for lobstr
				if (evaRecord.getReadPairedFlag()){
					if (evaRecord.getFirstOfPairFlag())
						hashKey = evaRecord.getReadName()+"/1";
					else 
						hashKey = evaRecord.getReadName()+"/2";
				}				

				if (ansHash.remove(hashKey) != null){
					FN --;
				}else{
					dup ++;
				}
			}else{
				FP ++;				
			}

			if (count % 10000000 == 0){
				//timer.mark("# Eva up to " + count);
				double precision = 1.0 * TP / (TP +FP), recall = 1.0 * TP /(TP + FN),
						Fscore = 2 *precision * recall /(precision + recall);		

				System.out.println("#Precision " + precision + " recall " + recall + " F2 "  + Fscore + " total " + count + " " + countAns+ " " + TP + " " + FP + " " + FN + " " + dup);
			}
		}			
		//System.out.println(" FN = " + FN + " Left over = " + ansHash.size());

		double precision = 1.0 * TP / (TP +FP), recall = 1.0 * TP /(TP + FN),
				Fscore = 2 *precision * recall /(precision + recall);		

		System.out.println("Precision " + precision + " recall " + recall + " F2 "  + Fscore + " total " + count + " " + countAns+ " " + TP + " " + FP + " " + FN + " " + dup);

	}

	public static int compareSamRecord(SAMRecord record1, SAMRecord record2){
		if (record1.getReferenceIndex() < record2.getReferenceIndex())
			return Integer.MIN_VALUE;
		if (record1.getReferenceIndex() > record2.getReferenceIndex())
			return Integer.MAX_VALUE;
		//they have to be on the same reference
		return record1.getAlignmentStart() - record2.getAlignmentStart();		
	}
}

class LinkRecord<T>{
	LinkRecord<T> mNext, mPrev;
	T mRecord;
	LinkRecord(T record){
		this(record,null,null);		
	}
	LinkRecord (T record, LinkRecord<T> next, LinkRecord<T> prev){
		mRecord = record;
		mNext = next;
		mPrev = prev;
	}

	T getRecord(){
		return mRecord;
	}	
}



