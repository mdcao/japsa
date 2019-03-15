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

import com.google.common.collect.ImmutableMap;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import htsjdk.samtools.*;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.DoubleArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

/**
 * @author Minh Duc Cao, Son Hoang Nguyen
 *
 */
public class RealtimeSpeciesTyping {
	private static final Logger LOG = LoggerFactory.getLogger(RealtimeSpeciesTyping.class);
	public static boolean JSON = false;
	public static double ALPHA=0.05;
	public static int MIN_READS_COUNT=0;

	private RealtimeSpeciesTyper typer;
	private OutputStream outputStream;
	private BufferedReader indexBufferedReader;
	private HashSet<String> filterSet = new HashSet<String>();
	/**
	 * Minimum quality of alignment
	 */
	private double minQual = 1;
	private boolean twoDOnly = false;


	Integer currentReadCount = 0;
	Integer currentReadAligned = 0;
	Long currentBaseCount = 0L;


	//seq ID to species name (from index ref file)
	HashMap<String, String> seq2Species = new HashMap<String, String>();
	//HashMap<String, SpeciesCount> species2Count = new HashMap<String, SpeciesCount>();
	ArrayList<String> speciesList = new ArrayList<String>(); 
	
	//to output binned sequences
	public static boolean OUTSEQ=false;
	HashMap<String, ArrayList<String>> species2ReadList = new HashMap<String, ArrayList<String>>();

	public RealtimeSpeciesTyping(String indexFile, String outputFile) throws IOException{
		LOG.debug("string string");
		this.indexBufferedReader = SequenceReader.openFile(indexFile);
		this.outputStream = SequenceOutputStream.makeOutputStream(outputFile);
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}

	public RealtimeSpeciesTyping(String indexFile, OutputStream outputStream) throws IOException {
		LOG.debug("string outputstream");
		this.indexBufferedReader = SequenceReader.openFile(indexFile);
		this.outputStream = outputStream;
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}

	public RealtimeSpeciesTyping(BufferedReader indexBufferedReader, String outputFile) throws IOException {
		LOG.debug("bufferedreader string");
		this.indexBufferedReader = indexBufferedReader;
		this.outputStream = SequenceOutputStream.makeOutputStream(outputFile);
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}

	public RealtimeSpeciesTyping(BufferedReader indexBufferedReader, OutputStream outputStream) throws IOException {
		LOG.debug("bufferedreader outputstream");
		this.indexBufferedReader = indexBufferedReader;
		this.outputStream = outputStream;
		typer = new RealtimeSpeciesTyper(this, outputStream);
		preTyping();
	}

//	static class SpeciesCount implements Comparable<SpeciesCount> {
//		String species;
//		int count = 0;
//
//		SpeciesCount (String s){
//			species = s;
//		}
//
//		/* (non-Javadoc)
//		 * @see java.lang.Comparable#compareTo(java.lang.Object)
//		 */
//		@Override
//		public int compareTo(SpeciesCount o) {		
//			return o.count - count;
//		}
//
//	}

	private void preTyping() throws IOException{
		String line = "";
		while ( (line = indexBufferedReader.readLine())!=null){
			if (line.startsWith("#"))
				continue;


			String sp=null,seq=null;
				
			String [] toks = line.split(">");
			if(toks.length < 2){
				LOG.info("Illegal speciesIndex file!");
				System.exit(1);
			}
				
			sp=toks[0].trim();
			seq=toks[1].split("\\s+")[0];

			if (seq2Species.put(seq, sp) != null)
				throw new RuntimeException("sequence " + seq +" presents multiple time");
//			else
//				LOG.info("==>adding " + seq + " to " + sp);
			
			if (species2ReadList.get(sp) == null){
//				LOG.info("add species: "+sp);
				species2ReadList.put(sp,new ArrayList<String>());
			}			
		}//while
		indexBufferedReader.close();
		LOG.info(seq2Species.size() + "   " + species2ReadList.size());
		speciesList.addAll(species2ReadList.keySet());

		//Write header				
	}

	/**
	 * @param minQual the minQual to set
	 */
	public void setMinQual(double minQual) {
		this.minQual = minQual;
	}

	/**
	 * @param twoOnly the twoOnly to set
	 */
	public void setTwoOnly(boolean twoOnly) {
		this.twoDOnly = twoOnly;
	}

	public void typing(String bamFile, int readNumber, int timeNumber) throws IOException, InterruptedException {
		InputStream bamInputStream;

		if ("-".equals(bamFile))
			bamInputStream = System.in;
		else
			bamInputStream = new FileInputStream(bamFile);

		typing(bamInputStream, readNumber, timeNumber);
	}

	/**
	 * @param filter the species keywords list (separated by comma) to excluded
	 */
	public void setFilter(String filter) {
		if(!filter.isEmpty()){
			String[] toks = filter.split(";");
			for(String tok:toks)
				if(!tok.isEmpty()){
					filterSet.add(tok);
				}
		}
		
	}
	public void typing(InputStream bamInputStream, int readNumber, int timeNumber) throws IOException, InterruptedException{
		//if (readNumber <= 0)
		//	readNumber = 1;			

		typer.setReadPeriod(readNumber);
		typer.setTimePeriod(timeNumber * 1000);

		LOG.info("Species typing ready at " + new Date());

		String readName = "", refName = "";

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamInputStream));
		SAMRecordIterator samIter = samReader.iterator();

		Thread thread = new Thread(typer);
		LOG.info("starting RealtimeSpeciesTyper thread");
		thread.start();
		LOG.info("started  RealtimeSpeciesTyper thread");
		HashSet<String> skipList = new HashSet<>();
		while (samIter.hasNext()){
//			try{
			SAMRecord sam = samIter.next();
			//LOG.info("sam read name = "+sam.getReadName());
			//if (firstReadTime <=0)
			//	firstReadTime = System.currentTimeMillis();

			if (this.twoDOnly && !sam.getReadName().contains("twodim")){
				continue;
			}

			if (!sam.getReadName().equals(readName)){
				readName = sam.getReadName();
				synchronized(this){
					currentReadCount ++;
					currentBaseCount += sam.getReadLength();
				}
			} 
			
			if (sam.getReadUnmappedFlag()){
				LOG.debug("failed unmapped check");
				continue;			
			}

			if (sam.getMappingQuality() < this.minQual) {
				LOG.debug("failed minQual check");
				continue;
			}

			if(skipList.contains(readName)){
				LOG.debug("filter {}", readName);
				continue;
			}
			
			refName = sam.getReferenceName();
			if(filterSet.contains(seq2Species.get(refName))){
				if(!sam.isSecondaryOrSupplementary())
					skipList.add(readName);
				continue;
			}
			
			String species = seq2Species.get(refName);
			if (species == null){
				throw new RuntimeException(" Can't find species with ref " + refName + " line " + currentReadCount );
			}

			//SpeciesCount sCount = species2Count.get(species);
			ArrayList<String> readList = species2ReadList.get(species);
			if (readList == null){
				throw new RuntimeException(" Can't find record with species " + species + " line " + currentReadCount );
			}
			if(readList.size()==0 || !readList.contains(readName))
				synchronized(this) {
					currentReadAligned ++;			

					readList.add(readName);
	
				}
//			}catch(Exception exc){
//				exc.printStackTrace();
//			}
		}//while

		//final run
		//typer.simpleAnalysisCurrent();

		typer.stopWaiting();//Tell typer to stop
		samIter.close();
		samReader.close();
	}	

	public static class RealtimeSpeciesTyper extends RealtimeAnalysis {
		MultinomialCI rengine;
		RealtimeSpeciesTyping typing;
		public SequenceOutputStream countsOS;


		public RealtimeSpeciesTyper(RealtimeSpeciesTyping t, OutputStream outputStream) throws IOException {
			typing = t;
			rengine = new MultinomialCI(ALPHA);

			countsOS = new SequenceOutputStream(outputStream);
			if(!JSON)
				countsOS.print("time\tstep\treads\tbases\tspecies\tprob\terr\ttAligned\tsAligned\n");
		}

		private void simpleAnalysisCurrent() throws IOException {
			//long step = lastTime;

			//Date date = new Date(lastTime);
			Long step = (lastTime - startTime)/1000;//convert to second

			int sum = 0;
			double [] count = new double[typing.speciesList.size()];
			for (int i = 0; i < count.length;i++){			
				count[i] = typing.species2ReadList.get(typing.speciesList.get(i)).size();			
				sum += count[i];
			}
			DoubleArray countArray = new DoubleArray();
			ArrayList<String> speciesArray = new ArrayList<String> ();

			int minCount = MIN_READS_COUNT>0?MIN_READS_COUNT:Math.max(1,sum/50);

			for (int i = 0; i < count.length;i++){			
				if (count[i] >= minCount){
					countArray.add(count[i]);
					speciesArray.add(typing.speciesList.get(i));
					LOG.info(step+" : " + typing.speciesList.get(i) + " == " + count[i]);
				}
			}		
			//if (countArray.size() > 10) return;
			countArray.add(1);
			speciesArray.add("others");		

			rengine.assignCount(countArray.toArray());
			rengine.eval();        
			//REXP tab  = rengine.eval("tab",true);  
			double [][] results =rengine.tab();

			Gson gson = new GsonBuilder().serializeNulls().create();
			List<JsonObject> data = new ArrayList<JsonObject>();

			for (int i = 0; i < results.length;i++){
				if (results[i][0] <= 0.00001)
					continue;

				Double mid = (results[i][0] + results[i][1])/2;
				Double err = mid - results[i][0];
				if(!JSON) {
					countsOS.print(timeNow + "\t" + step + "\t" + lastReadNumber + "\t" + typing.currentBaseCount + "\t" + speciesArray.get(i).replaceAll("_", " ") + "\t" + mid + "\t" + err + "\t" + typing.currentReadAligned + "\t" + countArray.get(i));
					countsOS.println();
				}
				else {
					JsonObject jo = new JsonObject();
					jo.addProperty("species", speciesArray.get(i).replaceAll("_", " "));
					jo.addProperty("step", step.toString());
					jo.addProperty("reads", lastReadNumber.toString());
					jo.addProperty("bases", typing.currentBaseCount.toString());
					jo.addProperty("prob", mid.toString());
					jo.addProperty("err", err.toString());
					jo.addProperty("tAligned", typing.currentReadAligned.toString());
					jo.addProperty("sAligned", Double.valueOf(countArray.get(i)).toString());
					data.add(jo);

				}
			}

			if(JSON) {
				countsOS.print(gson.toJson(ImmutableMap.of(
						"timestamp", timeNow.toString(),
						"data", data
				)));
				countsOS.println();

			}
			countsOS.flush();
			LOG.info(step+"  " + countArray.size());
		}

		protected void close(){
			try{
				//rengine.end();
				countsOS.close();
			}catch (Exception e){
				e.printStackTrace();
			}
			
			//print out
			if(OUTSEQ){
				try (BufferedWriter bw = new BufferedWriter(new FileWriter("species2reads.map"))) {
					for(String sp:typing.species2ReadList.keySet()){
						ArrayList<String> readList = typing.species2ReadList.get(sp);
						if(readList.size()==0)
							continue;
						bw.write(">"+sp+"\n");
						for(String read:readList)
							bw.write(read+"\n");
					}			

				} catch (IOException e) {

					e.printStackTrace();

				}
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#analysis()
		 */
		@Override
		protected void analysis(){
			try{
				simpleAnalysisCurrent();
			}catch (IOException e){
				LOG.warn(e.getMessage());
			}
		}

		/* (non-Javadoc)
		 * @see japsa.bio.np.RealtimeAnalysis#getCurrentRead()
		 */
		@Override
		protected int getCurrentRead() {
			return typing.currentReadCount;
		}
	}

}
