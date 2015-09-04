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
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.bio.np;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.DoubleArray;
import japsa.util.Logging;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.Rengine;


/**
 * @author minhduc
 *
 */
public class RealtimeSpeciesTyping {

	public double qual = 0;
	Rengine rengine;

	int currentReadCount = 0;
	int currentReadAligned = 0;
	long currentBaseCount = 0;

	int arrayIndex = 0;
	String prefix;
	public SequenceOutputStream countsOS;	
	long firstReadTime = 0;


	/////////////////////////////////////////////////////////////////////////////

	long startTime;
	public RealtimeSpeciesTyping(){	
		rengine = new Rengine (new String [] {"--no-save"}, false, null);
		if (!rengine.waitForR()){
			Logging.exit("Cannot load R",1);            
		}    
		rengine.eval("library(MultinomialCI)");
		rengine.eval("alpha<-0.05");

		Logging.info("REngine ready");
		startTime = System.currentTimeMillis();
	}

	public void close(){
		rengine.end();
	}
	/**
	 * @param bamFile
	 * @param geneFile
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	static class SpeciesCount implements Comparable<SpeciesCount>{
		String species;
		int count = 0;

		SpeciesCount (String s){
			species = s;
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(SpeciesCount o) {		
			return o.count - count;
		}

	}


	HashMap<String, String> seq2Species = new HashMap<String, String>();
	HashMap<String, SpeciesCount> species2Count = new HashMap<String, SpeciesCount>();
	ArrayList<String> speciesList = new ArrayList<String>(); 


	public void preTyping(String indexFile)throws IOException{
		BufferedReader bf = SequenceReader.openFile(indexFile);
		String line = "";
		while ( (line = bf.readLine())!=null){
			if (line.startsWith("#"))
				continue;

			String [] toks = line.split(" ");
			String sp =  toks[0];
			String seq = toks[1].substring(1);

			if (seq2Species.put(seq, sp) != null)
				throw new RuntimeException("sequence " + seq +" presents multiple time");

			if (species2Count.get(sp) == null){
				species2Count.put(sp,new SpeciesCount(sp));
			}			
		}//while
		bf.close();
		Logging.info(seq2Species.size() + "   " + species2Count.size());
		speciesList.addAll(species2Count.keySet());

		//Write header
		countsOS.print("step\treads\tbases\tspecies\tprob\terr\ttAligned\tsAligned\n");		
	}

	private void simpleAnalysisCurrent(int currentRead) throws IOException{		
		int step = currentRead;

		int sum = 0;
		double [] count = new double[speciesList.size()];
		for (int i = 0; i < count.length;i++){			
			count[i] = species2Count.get(speciesList.get(i)).count;			
			sum += count[i];
		}
		DoubleArray countArray = new DoubleArray();
		ArrayList<String> speciesArray = new ArrayList<String> ();

		int minCount = Math.max(1,sum/50);

		for (int i = 0; i < count.length;i++){			
			if (count[i] >= minCount){
				countArray.add(count[i]);
				speciesArray.add(speciesList.get(i));
				Logging.info(step+" : " + speciesList.get(i) + " == " + count[i]);
			}
		}		
		//if (countArray.size() > 10) return;
		countArray.add(1);
		speciesArray.add("others");		

		rengine.assign("count", countArray.toArray());
		rengine.eval("tab = multinomialCI(count,alpha)");        
		REXP tab  = rengine.eval("tab",true);  
		double [][] results = tab.asDoubleMatrix();

		for (int i = 0; i < results.length;i++){
			if (results[i][0] <= 0.00001)
				continue;

			double mid = (results[i][0] + results[i][1])/2;
			double err = mid - results[i][0];

			countsOS.print(step + "\t" + currentReadCount + "\t" + currentBaseCount + "\t" + speciesArray.get(i).replaceAll("_"," ") + "\t" + mid +"\t" + err + "\t" + this.currentReadAligned + "\t" + count[i]);
			countsOS.println();
		}

		countsOS.flush();
		Logging.info(step+"  " + countArray.size());
	}



	public void typing(String bamFile, int readNumber) throws IOException, InterruptedException{
		if (readNumber <= 0)
			readNumber = 1;


		String readName = "";
		//Read the bam file		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if ("-".equals(bamFile))
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			samReader = SamReaderFactory.makeDefault().open(new File(bamFile));

		SAMRecordIterator samIter = samReader.iterator();			

		while (samIter.hasNext()){
			SAMRecord sam = samIter.next();
			if (firstReadTime <=0)
				firstReadTime = System.currentTimeMillis();

			if (!sam.getReadName().equals(readName)){
				readName = sam.getReadName();

				currentReadCount ++;
				currentBaseCount += sam.getReadLength();


				if (currentReadCount % readNumber == 0){
					simpleAnalysisCurrent(currentReadCount);
				}

			}

			if (sam.getReadUnmappedFlag()){				
				continue;			
			}

			if (sam.getMappingQuality() < this.qual)
				continue;

			currentReadAligned ++;

			String refSequence = sam.getReferenceName();
			String species = seq2Species.get(refSequence);
			if (species == null){
				throw new RuntimeException(" Can find species with ref " + refSequence + " line " + currentReadCount );
			}
			SpeciesCount sCount = species2Count.get(species);
			if (sCount == null){
				throw new RuntimeException(" Can find record with species " + species + " line " + currentReadCount );
			}

			synchronized(this) {
				sCount.count ++;
			}			

		}//while

		//final run
		simpleAnalysisCurrent(currentReadCount);

		samIter.close();
		samReader.close();
	}	

}