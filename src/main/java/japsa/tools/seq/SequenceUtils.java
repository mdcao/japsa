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
import java.util.Iterator;
import java.util.Map;

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
			
			return iterator.next();
		}
		 
	 }
	
	
	public static String src_tag = "SC";
	public static String pool_tag = "PT";
	
	public static  class CombinedIterator implements Iterator<SAMRecord> {
		private final Iterator<SAMRecord>[] samIters;
		int current_sam_index =0;
		int currentIndex =0; //relates to chromosome
		private final SAMRecord[] currentVals;
		private final boolean[] returned;
		private final int[] cnts ;
		int max;
		Collection<String>[] readList ;
		Map<Integer, int[]> chrs;
		public  CombinedIterator(Iterator<SAMRecord>[] samIters, int max, Collection<String>[]readList, Map<Integer, int[]> chrs) {
			this.samIters = samIters;
			this.readList = readList;
			currentVals = new SAMRecord[samIters.length];
			returned = new boolean [samIters.length];
			this.max = max;
			cnts = new int[samIters.length];
			this.chrs = chrs.size()>0 ? chrs : null;
		}
		/*public void close(){
			for(int i=0; i<samIters.length; i++){
				samIters[i].close();
			}
		}*/
		@Override
		public boolean hasNext() {
			for(int i=0; i<samIters.length; i++){
				if(samIters[i].hasNext() && cnts[i]<max) return true;
			}
		//	this.close();
			return false;
		}
		int pool_ind=-1;
		private SAMRecord next(int i){
			SAMRecord sr;
			if(currentVals[i]!=null && !returned[i]){
				sr = currentVals[i];
			}else{
				if(cnts[i]<max){
					
					//sr  = samIters[i].next();
				/*	while(!(
							sr==null || 
							(readList==null || readList.contains(sr.getReadName())) ||
							(chrs==null || chrs.contains(sr.getReferenceIndex()))
							)
							)*/
					inner: while(true)
					{
						sr  = samIters[i].next();
						if(sr==null) break inner;
						pool_ind =-1;
						if(readList!=null){
							inner1: for(int i2=0; i2<readList.length; i2++){
								if(readList[i2].contains(sr.getReadName())){
									pool_ind = i2;
									break inner1;
								}
							}
						}
						if(readList!=null && pool_ind>=0) break inner;
						if(readList==null &&  chrs!=null && contains(chrs,sr)) break inner;
						if(readList==null  && chrs==null) break inner;
					}
					if(sr!=null) {
						sr.setAttribute(pool_tag, pool_ind);
						sr.setAttribute(src_tag, i);
					}
					currentVals[i] = sr;
					returned[i] = false; 
				}else{
					sr = null;
				}
			}
			return sr;
		}

		private boolean contains(Map<Integer, int[]> chrs2, SAMRecord sr) {
			int[] obj = chrs2.get(sr.getReferenceIndex());
			if(obj==null) return false;
			else{
				if(sr.getAlignmentStart()>=obj[0] && sr.getAlignmentEnd() <=obj[1]) return true;
				else return false;
			}
		}
		@Override
		public SAMRecord next() {
			//int curr_index = current_sam_index;
			SAMRecord sr = next(current_sam_index);
			
			if(sr==null || sr.getReferenceIndex()>currentIndex ){
				int[] ref_inds = new int[samIters.length];
				int min_ind =-1;
				int minv = Integer.MAX_VALUE;
				for(int i=0; i<samIters.length; i++){
					SAMRecord sr_i = next(i);
					if(sr_i!=null) {
						ref_inds[i] = sr_i.getReferenceIndex(); //note this only advances if not returned;
						if(ref_inds[i]<minv){
							min_ind = i;
							minv = ref_inds[i];
						}
					}
				}
				current_sam_index = min_ind;
			}
			
			if(current_sam_index<0) return null;
			sr = this.currentVals[current_sam_index];
			if(sr!=null) this.currentIndex = sr.getReferenceIndex();
			cnts[current_sam_index]++;
			returned[current_sam_index] = true; 
			return sr;
			
		}
	}
	public static Iterator<SAMRecord>  getCombined(Iterator<SAMRecord>[] samReaders, Collection[] reads, int max_reads,String chroms){
		Map<Integer, int[]>chromIndices = new HashMap<Integer, int[]>();
		return getCombined(samReaders,reads, max_reads, chroms, chromIndices);
	}
	
	//chroms is string e.g 0:1:2  or 0,0,2400000:1,0,240000
public static Iterator<SAMRecord>  getCombined(Iterator<SAMRecord>[] samReaders, Collection[] reads, int max_reads,String chroms,
		Map<Integer, int[]>chromIndices ){
//	Map<Integer, int[] > chromIndices = null;
	if(chroms!=null && !chroms.equals("all")){
		
		String[] chri = chroms.split(":");
		for(int i=0; i<chri.length; i++){
			String[] str = chri[i].split(",");
			int[] vals ;
			if(str.length==1) vals = new int[] {0,Integer.MAX_VALUE};
			else{
				vals = new int[] {Integer.parseInt(str[1]), Integer.parseInt(str[2])};
			}
			chromIndices.put(Integer.parseInt(str[0]), vals);
		}
	}
		int len = samReaders.length;
	//	SAMRecordIterator[] samIters = new SAMRecordIterator[len];
		//for (int ii = 0; ii < len; ii++) {
		//	samIters[ii] = samReaders[ii].iterator();
		//}
		//if(reads==null && chrom_indices_to_include==null && samReaders.length==1){
		//	samIter = samIters[0];
	//	}else{
			return 		new CombinedIterator(samReaders, max_reads,reads, chromIndices);
		//}
		
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
}
