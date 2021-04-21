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

package japsa.bio.hts;

import japsa.bio.np.MultinomialCI;
import japsa.seq.SequenceReader;
import japsa.util.DoubleArray;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author minhduc
 *
 */
public class SpeciesCoverageIdenfication {
	private static final Logger LOG = LoggerFactory.getLogger(SpeciesCoverageIdenfication.class);

	private double qual = 0;
	private MultinomialCI rengine;
	
	private int currentReadCount = 0;
	private int currentReadAligned = 0;
	//private long currentBaseCount = 0;
	//private long currentBaseAligned= 0;

	//int arrayIndex = 0;	
	//String prefix;
	
	//public SequenceOutputStream countsOS;
	private PrintStream outOS;

	public SpeciesCoverageIdenfication(String outputFile, double minQual) throws IOException{		
		
		rengine = new MultinomialCI(0.05);

		//countsOS = SequenceOutputStream.makeOutputStream(outputFile);
		if (outputFile.equals("-"))
			outOS = System.out;
		else
			outOS = new PrintStream (new FileOutputStream(outputFile));
		
		this.qual = minQual;
	}

	public void close() throws IOException{
		outOS.close();
		//rengine.end();
	}


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

			//if (
					seq2Species.put(seq, sp);
					//!= null)
			//	throw new RuntimeException("sequence " + seq +" presents multiple time");

			if (species2Count.get(sp) == null){
				species2Count.put(sp,new SpeciesCount(sp));
			}			
		}//while
		bf.close();
		LOG.info(seq2Species.size() + "   " + species2Count.size());
		speciesList.addAll(species2Count.keySet());

		//Write header
		outOS.println("Species\tProportion\tError\tRead Count\tTotal Reads Aligned\tTotal Reads");
		//for (String species:speciesList){
		//	countsOS.print("\t" + species);
		//}	
		//countsOS.println();
	}

	private void simpleAnalysisCurrent() throws IOException{		
		
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
			}
		}
		
		countArray.add(1);
		speciesArray.add("others");		

		rengine.assignCount(countArray.toArray());
		rengine.eval();        
	//	REXP tab  = rengine.eval("tab",true);  
		double [][] results = rengine.tab();
		
		for (int i = 0; i < results.length;i++){
			if (results[i][0] <= 0.00001)
				continue;

			double mid = (results[i][0] + results[i][1])/2;
			double err = mid - results[i][0];

			//Species 
			outOS.printf("%s\t%.4f\t%.4f\t%d\t%d\t%d\n", speciesArray.get(i).replaceAll("_"," "), mid, err, (int) countArray.get(i),currentReadAligned, currentReadCount);
			
		}
		outOS.flush();
		//LOG.info(step+"  " + countArray.size());
	}

	
	public void typing(String bamFile) throws IOException, InterruptedException{
		
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
			currentReadCount ++;
			
			//int readLength = sam.getReadLength();
			//currentBaseCount += readLength;
			
			if (sam.getReadUnmappedFlag()){				
				continue;			
			}

			if (sam.getMappingQuality() < this.qual)
				continue;

			currentReadAligned ++;
			//currentBaseAligned  += readLength;

			String refSequence = sam.getReferenceName();
			String species = seq2Species.get(refSequence);
			if (species == null){
				throw new RuntimeException(" Can find species with ref " + refSequence + " line " + currentReadCount );
			}
			SpeciesCount sCount = species2Count.get(species);
			if (sCount == null){
				throw new RuntimeException(" Can find record with species " + species + " line " + currentReadCount );
			}			
			sCount.count ++;						

		}
		
		//final run
		simpleAnalysisCurrent();

		samIter.close();
		samReader.close();
	}	
}