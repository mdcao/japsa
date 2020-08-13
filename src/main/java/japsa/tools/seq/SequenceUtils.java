/*****************************************************************************
 * Copyright (c) Lachlan Coin, University of Meblourne, All rights reserved.         *
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


package japsa.tools.seq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SequenceUtils {
	
	public static Iterator<SAMRecord> getSAMIteratorFromFastq(File inFile, String mm2Index, String mm2_path, 
			int mm2_threads, String mm2Preset, String mm2_mem, String mm2_splicing) throws IOException{
		return new FastqToSAMRecord(inFile, mm2Index, mm2_path,mm2_threads, mm2Preset, mm2_mem , mm2_splicing);
	}
	
	public static Iterator<SAMRecord> getSAMIteratorFromFastq(File inFile, String mm2Index, String mm2_path, 
			int mm2_threads, String mm2Preset, String mm2_mem) throws IOException{
		return new FastqToSAMRecord(inFile, mm2Index, mm2_path,mm2_threads, mm2Preset, mm2_mem , null);
	}
	
	
	private static class FastqToSAMRecord implements Iterator<SAMRecord> {
		// ProcessBuilder pb;
		 SamReader reader;
		 SAMRecordIterator iterator;
		public FastqToSAMRecord(File inFile, String mm2Index, String mm2_path, int mm2_threads, String mm2Preset, String mm2_mem, String splicing) throws IOException{
			ProcessBuilder pb;
			
			if(splicing==null) pb = new ProcessBuilder(mm2_path, 
					"-t",
					"" + mm2_threads,
					"-ax",
					mm2Preset,
				//	"--for-only",
					"-I",
					mm2_mem,
//					"-K",
//					"200M",
					mm2Index,
					"-"
				
					);
			else pb = new ProcessBuilder(mm2_path, 
					"-t",
					"" + mm2_threads,
					"-ax",
					mm2Preset,
					splicing,
				//	"--for-only",
					"-I",
					mm2_mem,
//					"-K",
//					"200M",
					mm2Index,
					"-"
				
					);
			
			Process mm2Process  = pb.redirectInput(ProcessBuilder.Redirect.from(inFile)).redirectError(ProcessBuilder.Redirect.to(new File("err.txt"))).start();
			//	Process mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("err.txt"))).start();
			//	OutputStream os = mm2Process.getOutputStream();
				reader =  SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));
				iterator = reader.iterator();
		 }
		@Override
		public boolean hasNext() {
			boolean res = iterator.hasNext();
			try{
			if(!res) reader.close();
			}catch(IOException exc){
				exc.printStackTrace();
			}
			return res;
		}

		@Override
		public SAMRecord next() {
			if(!iterator.hasNext()){
				try{
				reader.close();
				}catch(IOException exc){
					exc.printStackTrace();
				}
				return null;
			}
			return iterator.next();
		}
		 
	 }
	
	
	public static String src_tag = "SC";
	public static String pool_tag = "PT";
	
	 static class SamR implements Comparable{
		final int pos, ref, ind;
		public SamR(SAMRecord samRecord, int i) {
			ref = samRecord.getReferenceIndex();
			pos = samRecord.getAlignmentStart();
			ind = i;
			
		}
		public String toString(){
			return ref+","+pos+","+ind;
		}

		@Override
		public int compareTo(Object o) {
			SamR o1 = (SamR)o;
			int res = Integer.compare(ref, o1.ref);
			if(res==0) res = Integer.compare(pos, o1.pos);
			if(res==0) res = Integer.compare(ind, o1.ind);
			return res;
		}
		
	}
	
	public static class CombinedIterator implements Iterator<SAMRecord> {
		private final Iterator<SAMRecord>[] samIters;
		private final SAMRecord[] currentVals; 
		private final TreeSet<SamR> samR = new TreeSet<SamR>();
		
		final int len;
		int nullcount=0;
		public  CombinedIterator(Iterator<SAMRecord>[] samIters) {
			this.samIters = samIters;
			len = samIters.length; 
			currentVals = new SAMRecord[samIters.length];
			for(int i=0; i<len; i++){
				if(samIters[i].hasNext()){
					SAMRecord sr =   samIters[i].next();
					sr.setAttribute(SequenceUtils.src_tag,i);
					samR.add(new SamR(sr,i));
					currentVals[i] = sr;
					
				}
			}
			System.err.println("h");
		}
		
		@Override
		public boolean hasNext() {
			return nullcount<len;
//			return samR.size()>0;
		}
		

		@Override
		public SAMRecord next() {
			SamR  first = samR.first(); //next element
			samR.remove(first);
			int i = first.ind;
			SAMRecord sr = currentVals[i];
			sr.setAttribute(SequenceUtils.src_tag,i);
			if(samIters[i].hasNext()){
				//replace element removed if possible
				currentVals[i] = samIters[i].next();
				samR.add(new SamR(currentVals[i],i));
			}else{
				System.err.println("iterator"+i+" exhausted "+sr.getReferenceName()+" "+sr.getAlignmentStart());
				nullcount++;
				System.err.println(nullcount+" of "+len+" finished");
				currentVals[i] = null;
			}
			return sr;
		}
		
	}
	
	
	
	
	//chroms is string e.g 0:1:2  or 0,0,2400000:1,0,240000
public static Iterator<SAMRecord>  getCombined(Iterator<SAMRecord>[] samReaders){
		int len = samReaders.length;
		if(samReaders.length==1  && false){
			return samReaders[0];
		}else{
			return 		new CombinedIterator(samReaders);
		}
		
	}

public static String minimapIndex(File refFile, String mm2, String mem, boolean overwrite) throws IOException, InterruptedException {
	File indexFile = new File(refFile.getAbsolutePath()+".mmi");
	if(indexFile.exists() && !overwrite){
		System.out.println("minimap index exists");
		
	}else{
		ProcessBuilder pb = new ProcessBuilder(mm2, 
				"-I",
				mem,
				"-d",
				indexFile.toString(),
				refFile.toString()
				);
		//System.err.println(pb.toString());
			Process p =  pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
			p.waitFor();
	}
	return indexFile.getAbsolutePath();
}

public static Iterator<SAMRecord> getCombined(Iterator<SAMRecord>[] samIters, Collection[] reads, Integer max_reads,
		String chrToInclude) {
	final Map<Integer, int[]> chrom_indices_to_include = new HashMap<Integer, int[]>();
	final Set<String> reads_all = new HashSet<String>();
	for(int i=0; i<reads.length; i++){
		for(Iterator<String> it = reads[i].iterator(); it.hasNext();){
			reads_all.add(it.next());
		}
	}
//	Map<String, int[]> chromsToInclude = new HashMap<String, int[]>();
	if(chrToInclude!=null && !chrToInclude.equals("all")){
		String[] chri = chrToInclude.split(":");
		for(int i=0; i<chri.length; i++){
			String[] str = chri[i].split(",");
			int[] vals ;
			if(str.length==1) vals = new int[] {0,Integer.MAX_VALUE};
			else{
				vals = new int[] {Integer.parseInt(str[1]), Integer.parseInt(str[2])};
			}
			chrom_indices_to_include.put(Integer.parseInt(str[0]), vals);
		}
	}
//	if(chrom_indices_to_include.size()==0) chrom_indices_to_include=null;
	final Iterator<SAMRecord>  sr = getCombined(samIters);
	return new Iterator<SAMRecord>(){

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return sr.hasNext();
		}

		@Override
		public SAMRecord next() {
			SAMRecord nxt = sr.next();
			while(
					!(chrom_indices_to_include.containsKey(nxt.getReferenceIndex()) ||
							 reads_all.contains(nxt.getReadName())
							)
							){
				if(!sr.hasNext()) return null;
				nxt = sr.next();
			}
			return nxt;
		}
		
	};
}
}
