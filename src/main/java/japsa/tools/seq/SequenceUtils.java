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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.ProcessBuilder.Redirect;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

public class SequenceUtils {
	public static boolean deleteErr = true;
	public static SAMFileWriterFactory sfw = new SAMFileWriterFactory();
	{
		sfw.setUseAsyncIo(true);
	}
	public static void flip(SAMRecord sam, boolean switchFlag) {
		 if(true) throw new RuntimeException(" not clear effect on coordinates in read space");
		String sa = sam.getReadString();
		byte[]phredQs = sam.getBaseQualities();
		byte[] bases = sam.getReadBases();
		reverseArray(phredQs);
		reverseArray(bases);
		sam.setReadBases(bases);
		sam.setBaseQualities(phredQs);
		sa = SequenceUtil.reverseComplement(sa);
		sam.setReadString(sa);
		if(switchFlag){
			sam.setReadNegativeStrandFlag(!sam.getReadNegativeStrandFlag());
		}
	}
	public static void reverseArray(byte[] intArray) {
		 int size = intArray.length;
	        int i, k;
	        byte temp; 
	        for (i = 0; i < size / 2; i++) { 
	            temp = intArray[i]; 
	            intArray[i] = intArray[size - i - 1]; 
	            intArray[size - i - 1] = temp; 
	        } 
	 } 
	
	/** This processes input fastq */
	public static class PipeConnector implements Runnable {
		 final Collection<String> readsToInclude;
		final double q_thresh;
	//  private BufferedReader  _process1Output = null;
		  private Iterator<FastqRecord> _process1Output  = null;//
		  private BasicFastqWriter _process2Input  = null;
		  private SequenceOutputStream sos = null;
		  private int max;
		  final boolean fasta;
		  /**
		   * Initialize the PipeConnector
		   *
		   * @param  process1Output  The output stream from the first process (read as an InputStream)
		   * @param  process2Input   The input stream to the second process (written as an OutputStream)
		 * @param fasta 
		   */
		  public PipeConnector (Iterator<FastqRecord> process1Output, OutputStream process2Input, int max, Collection<String> readsToInclude, double q_thresh, boolean fasta) {
		    _process1Output = process1Output;
		    this.fasta = fasta;
		    this.readsToInclude = readsToInclude;
		    this.q_thresh = q_thresh;
		    if(fasta) sos = new SequenceOutputStream(process2Input);
		    else _process2Input  = new BasicFastqWriter(new PrintStream(process2Input));
		    this.max = max;
		  }
		// String st0,st1,st2,st3;
		  
		  /**
		   * Perform the copy operation in a separate thread
		   */
		  public void run () {
				try{
		 FastqRecord nxt = null;
		   boolean process = true;
		   for(int i =0; i<max; i++) {
		        nxt= _process1Output.next();
		        if (nxt == null){
		        	break;  // end of input stream
		        }
		     
		        if(readsToInclude==null){
		        	process=true;
		        }else{
		        	String readname = nxt.getReadName();
			        process = readsToInclude.remove(readname);
		        }
		        if(!fasta){
			        double q = SequenceUtils.getQual(nxt.getBaseQualities());
			        if(q<q_thresh) process = false;
		        }
		        if(process){
		        	if(fasta){
		        		String nme=">"+nxt.getReadName();
		        		this.sos.print(nme);sos.println();
		        		sos.print(nxt.getReadString()); sos.println();
		        		sos.flush();
		        	}else{
		        	_process2Input.write(nxt);
		        	_process2Input.flush();

		        	}
		        	//_process2Input.println();
		        }
		        if(readsToInclude!=null && readsToInclude.size()==0) break; // no more reads to include
		    }
		   
		System.err.println("finished piping input data "+fasta);
		   
		      //_process1Output.close();
		   if(fasta)sos.close();
		   else    _process2Input.close();
		  
		  }catch(IOException exc){
  			exc.printStackTrace();
  			
  		}
		  }
		}
	
	
	
public static void main(String[] args){
	try{
		SequenceUtils.apboa_path = "/home/lachlan/abPOA-v1.0.1/bin/abpoa";
		File f = new File("/home/lachlan/WORK/AQIP/japsa_species_typing/plasmids/aqip003.fastq/fastqs");
		SequenceUtils.makeConsensus(f, 4, true);
//		String refFile = "/home/lachlan/github/npTranscript/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz";
//		mm2_path="/home/lachlan/github/minimap2/minimap2";
//		String mm2_index = SequenceUtils.minimapIndex(new File(refFile),  false, true);
//		Iterator<SAMRecord>  sm = getSAMIteratorFromFastq(new String[]{"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4082025/ERR4082025_1.fastq.gz"},
//				mm2_index, 100, null,0, null);
///		while(sm.hasNext()){
	//		System.err.println(sm.next().getAlignmentStart());
		//}
	}catch(Exception exc){
		exc.printStackTrace();
	}
	}
	
	public static Iterator<SAMRecord> getSAMIteratorFromFastq(String[] url, String mm2Index, int maxReads, Collection<String>readsToInclude, 
			double q_thresh, Stack<File> bamOut) throws IOException{
		FastqToSAMRecord it =  new FastqToSAMRecord(url, mm2Index,maxReads, readsToInclude, q_thresh , bamOut!=null);
		if(bamOut!=null) bamOut.push(it.outputBAM);
		return it;
	}
	
	public static int mm2_threads=4;
	public static int max_secondary=10;
	public static String mm2_path="minimap2";
	public static String mm2_mem = "1000000000";
	public static String mm2Preset="splice";
	public static String mm2_splicing="-un";
	public static boolean secondary =true;
	
	public static String apboa_path="abpoa";
	
	public static void waitOnThreads(ExecutorService executor, int sleep) {
		if(executor==null) return ;
		if(executor instanceof ThreadPoolExecutor){
	    	while(((ThreadPoolExecutor) executor).getActiveCount()>0){
	    		try{
		    	System.err.println("OUTPUTS: awaiting completion "+((ThreadPoolExecutor)executor).getActiveCount());
		    	//Thread.currentThread();
				Thread.sleep(sleep);
	    		}catch(InterruptedException exc){
	    			exc.printStackTrace();
	    		}
	    	}
	    	}
		
	}
	
	static void printCommand(ProcessBuilder pb) {
		List<String> cmd = pb.command();
		StringBuffer sb  = new StringBuffer();
		for(int i=0; i<cmd.size(); i++){
			sb.append(cmd.get(i)+" ");
		}
		System.err.println(sb.toString());
		
	}
public static File makeConsensus(File file, int threads, boolean deleteFa) {
	File output = new File(file.getParentFile(), "consensus_output.fa");
	try{
	MultiAbpoa mab = new MultiAbpoa(file, threads,Alphabet.DNA16(), output, deleteFa);
	mab.run();
	
	}catch(Exception exc){
		exc.printStackTrace();
	}
	return output;
	}
	
	

	public static class MultiAbpoa{
		final ExecutorService executor;
		final SequenceOutputStream pw;
		final List<File> files = new ArrayList<File>();
		final Alphabet alph;
		MultiAbpoa(File file, int threads, Alphabet alph, File out, boolean del) throws FileNotFoundException{
			pw = new SequenceOutputStream(new FileOutputStream(out));
			executor = Executors.newFixedThreadPool(threads);
			this.alph = alph;
			File[] f = file.listFiles(new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.listFiles().length>0;
				}
				
			});
			for(int i=0; i<f.length; i++){
				files.addAll(Arrays.asList(f[i].listFiles()));
			}
			if(del){
				file.deleteOnExit();
				for(int i=0; i<f.length; i++) f[i].deleteOnExit();
				for(int i=0; i<files.size(); i++){
					files.get(i).deleteOnExit();
				}
				
				
			}
		}
		public void run() throws IOException{
			for(int i=0; i<files.size(); i++){
				try{
			
			executor.execute(new Abpoa(files.get(i), alph));
				}catch(Exception exc){
					System.err.println("prob with "+files.get(i));
					exc.printStackTrace();
				}
			}
			waitOnThreads(executor,100);
			executor.shutdown();
			pw.close();
		}
		
		synchronized void print(Sequence seq) throws IOException{
			seq.print(pw);
		}
		
		class Abpoa implements Runnable{
		//FastaReader br;
			Process proc;
			Alphabet alph;
			String name;
			String desc;
			File in;
			public Abpoa(File in, Alphabet alph) throws IOException{
				this.alph = alph;
				this.in = in;
				this.name = in.getName();
				this.desc = in.getParentFile().getName();//.replace(" ", "_");
				if(true){
					InputStream is = new FileInputStream(in);
					Sequence seq1 = FastaReader.read(is, alph);
					String desc1 = seq1.getDesc();
					desc= desc+" ; "+desc1;
					is.close();
				}
				//pw = new PrintWriter(new OutputStreamWriter(proc.getOutputStream()));
			}
			@Override
			public void run() {
				try{
					ProcessBuilder 	pb = new ProcessBuilder(apboa_path, 
							in.getAbsolutePath()
							);
					printCommand(pb);
					String extra = "";
					
					//System.err.println("HERE");
					//System.err.println(desc1);
					
					proc =  pb.redirectError(Redirect.INHERIT).start();//redirectError(ProcessBuilder.Redirect.to(new File("err_minimap2.txt"))).start();
				//	br  = new FastaReader(proc.getInputStream());
					Sequence seq = FastaReader.read(proc.getInputStream(),alph);
					
				
					
					seq.setName(name);seq.setDesc(desc+extra);
					print(seq);
				}catch(Exception exc){
					System.err.println("WARNING:  problem with building consensus for"+in.getAbsolutePath());
					//exc.printStackTrace();
				}
				
				// TODO Auto-generated method stub
				
			}
			
		}
		
	}
	
	 
	 static Iterator<FastqRecord> getFastaIterator(InputStream ins) throws IOException {
		 final BufferedReader br = new BufferedReader(new InputStreamReader(ins));//FastaReader fr = new FastaReader(ins);
		 String nxtLine1 = br.readLine();
		 if(nxtLine1==null){
			 System.err.println("warning, no entries");
			 return null;
		 }
		/* if(!br.hasNext()) {
			 throw new RuntimeException ("!!");
			 //
		//	 System.err.println(br.readLine());
			 return null;
		 }*/
		return new Iterator<FastqRecord>(){
			String nxtLine = nxtLine1;
			
			@Override
			public boolean hasNext() {
					boolean hasNext =  nxtLine!=null;
					
					return hasNext;
			}

			@Override
			public FastqRecord next() {
				FastqRecord fq = null;
				try{
				if(nxtLine==null) {
					br.close();
					return null;
				}
				
					int i1 = Math.max(nxtLine.indexOf(' '), nxtLine.length());
					String nme = nxtLine.substring(1, i1);
					StringBuffer sb = new StringBuffer();
					while((nxtLine=br.readLine())!=null){
						if(nxtLine.startsWith(">")) break;
						sb.append(nxtLine);
					}
					if(nxtLine==null) br.close();
				
				fq  =   new FastqRecord(nme,	sb.toString(),	null,null);
				}catch(IOException exc){
					exc.printStackTrace();
				}
				
				return fq;
			}};
		}
	
 static Iterator<FastqRecord> getFastqIterator(BufferedReader br) throws IOException {
	 String st0_ = br.readLine();
	 if(st0_==null) return null; // if file is empty
	Iterator<FastqRecord> it1 = new Iterator<FastqRecord>(){
		String st1,st2,st3;
		String st0 = st0_;
		@Override
		public boolean hasNext() {
			boolean hasNext =  st0!=null;
			if(!hasNext) try{
				br.close();
			}catch(IOException exc){
				exc.printStackTrace();
			}
			return hasNext;
		}

		@Override
		public FastqRecord next() {
			FastqRecord fq = null;
			try{
				if(st0==null) {
					br.close();
					return null;
				}
				
			st1 = br.readLine();
			
				st2 = br.readLine();
				st3 = br.readLine();
			
			fq  =   new FastqRecord(st0.split(" ")[0].substring(1),	st1,	st2,st3);
			 st0=br.readLine();
			
			}catch(IOException exc){
				exc.printStackTrace();
			}
			return fq;
		}};
		return it1;
	}
	
	
 static Iterator<FastqRecord> getFastqIterator(	final SamReader samR) throws IOException {
		
		return new Iterator<FastqRecord>(){
			Iterator<SAMRecord>sams = samR.iterator();
			@Override
			public boolean hasNext() {
				boolean nxt =  sams.hasNext();
				if(!nxt)this.close();
				return nxt;
			}

			public void close(){
				try{
					samR.close();
					}catch(IOException exc){
						exc.printStackTrace();
					}
			}
			@Override
			public FastqRecord next() {
				SAMRecord sam = sams.next();
				while(sam.isSecondaryOrSupplementary() && sams.hasNext()){
					sam = sams.next();
					if(sam==null) {
						close();
						return null;
					}
				}
				if(sam!=null && sam.isSecondaryOrSupplementary()) sam = null;
				if(sam==null){
					close();
					return null;
				}
				return  new FastqRecord(sam.getReadName(),	sam.getReadString(),"+",sam.getBaseQualityString());
			}
			
		};
	}
	
   
	private static class FastqToSAMRecord implements Iterator<SAMRecord> {
		// ProcessBuilder pb;
	//	 ExecutorService ex1 = Executors.newFixedThreadPool(1);
		InputStream is = null;
		 SamReader reader;
		SAMRecordIterator iterator;
		final String mm2Index;
		final double q_thresh;
		final String[] input;
		int max_per_file;
		final Collection<String> readsToInclude;
		 SAMTextWriter bfw = null;
		Log log = null;
		 ProgressLogger pl = null; 
		 int cntsAdded = 0;
		 
		 public boolean finished(){
			 return this.bfw!=null || pl.getCount()==cntsAdded;
		 }
		 public void refresh(){
			 try{
				 System.err.println("refreshing: analysed "+count+" records");
			 is.close();
			if(reader!=null) reader.close();
			if(iterator!=null) iterator.close();
			 iterator = null;
			 }catch(Exception exc){
				 exc.printStackTrace();
			 }
		 }
		
		public void close(){
			try{
			System.err.println("closing: analysed "+count+" records "+Arrays.asList(input));
			System.err.println("writing alignments "+bfw==null);
			this.refresh();
			if(bfw!=null) { // this is bfw is async
				while(!finished()){
					System.err.println("sleep "+100);
					 Thread.currentThread().sleep(100);
				}
				this.bfw.close();
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		
			//ex1.shutdown();
		}
	//	final boolean deleteFile;
		private void init(int k) {
			reader=null; iterator = null;
			try{
			ProcessBuilder pb;
			System.err.println("making builder");
			if(mm2_splicing==null) {
				if(mm2Preset==null){
					if(!secondary){
					pb = new ProcessBuilder(mm2_path, 	"-t","" + mm2_threads,"-a","--secondary=no",	"-I",	mm2_mem,	mm2Index,	"-");
					}else{
						pb = new ProcessBuilder(mm2_path, "-t","" + mm2_threads,	"-a","-I",mm2_mem,"-N",""+max_secondary,	mm2Index,	"-");
					}
				}else{
				pb = new ProcessBuilder(mm2_path, 	"-t","" + mm2_threads,"-ax",	mm2Preset,	"-I",	mm2_mem,	"-N",""+max_secondary,	mm2Index,"-");
				}
			}
			else pb = new ProcessBuilder(mm2_path, 
					"-t","" + mm2_threads,"-ax",mm2Preset,mm2_splicing,	"-I",	mm2_mem,	mm2Index,	"-");
			printCommand(pb);
			boolean bam = false;
				boolean fasta=false;
				SamReader samReader= null;
				SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
				if(inputFile==null){
					URL  url = new URL(input[k]);
					URLConnection urlc = url.openConnection();
					is= input[k].endsWith(".gz")  ? new GZIPInputStream(urlc.getInputStream()) : urlc.getInputStream();
				}else if(input[k].endsWith(".bam") || input[k].endsWith(".sam")){
					is=	new FileInputStream(inputFile[k]);
					bam = true;
					samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(is));
				}else{
					if(!inputFile[k].exists()) throw new RuntimeException("!!! does not exists");
					boolean gz = input[k].endsWith(".gz");
					fasta = input[k].endsWith(".fa") || input[k].endsWith(".fasta") || input[k].endsWith(".fa.gz") || input[k].endsWith(".fasta.gz");

					is= gz   ? new GZIPInputStream(new FileInputStream(inputFile[k])) : new FileInputStream(inputFile[k]);
					System.err.println("input file "+input[k]+" "+fasta);

				}
				Iterator<FastqRecord> fastqIt1 = bam ? 
						getFastqIterator(samReader):
					(fasta ? getFastaIterator(is)	 : getFastqIterator(new BufferedReader(new InputStreamReader(is))));
			//	BufferedIterator<FastqRecord> fastqIt = fastqIt1==null ? null :   new BufferedIterator<FastqRecord>(fastqIt1,5);
			if(fastqIt1!=null){
				
				if(false){
					File  tmp = new File("tmp."+System.currentTimeMillis()+"."+Math.random());
					this.print(tmp, fasta, fastqIt1);
					System.err.println("tmp : "+tmp.getAbsolutePath());
					tmp.deleteOnExit();
				}
				File minimaperr = new File("err_minimap."+System.currentTimeMillis()+"."+Math.random());
				if(deleteErr) minimaperr.deleteOnExit();
				Redirect error = ProcessBuilder.Redirect.to(minimaperr);//Redirect.INHERIT
				Process mm2Process =  pb.redirectInput(ProcessBuilder.Redirect.PIPE).redirectError(error).start();
			//	Process mm2Process =  pb.redirectInput(ProcessBuilder.Redirect.from(tmp)).redirectError(error).start();
				PipeConnector pc = new PipeConnector(fastqIt1, mm2Process.getOutputStream(), max_per_file, readsToInclude,q_thresh, fasta);
				Thread th = new Thread(pc);
				th.start();
				//.redirectError(Redirect.INHERIT).start();
				reader =  SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));
				iterator = reader.iterator();
				if(!iterator.hasNext()){
					System.err.println("WARNING: nothing to return");
				}
			//pc.run();
//			ex1.execute(run);
//			Thread th = new Thread(pc);
//			th.start();
			//	Process mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("err.txt"))).start();
			//	OutputStream os = mm2Process.getOutputStream();
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
				//pc.run();
				
		}
		
	 
		//** prints fastqIterator directly to file
		private void print(File tmp, boolean fasta, Iterator<FastqRecord> fastqIt1) throws IOException {
			
				if(fasta){
					SequenceOutputStream sos = new SequenceOutputStream(new FileOutputStream(tmp));
					while(fastqIt1.hasNext()){
						FastqRecord nxt = fastqIt1.next();
						String nme=">"+nxt.getReadName();
	        			sos.print(nme);sos.println();
	        		sos.print(nxt.getReadString()); sos.println();
					}
	        		sos.close();
				}else{
				FastqWriter pw = new BasicFastqWriter(new PrintStream(new FileOutputStream(tmp)));

				while(fastqIt1.hasNext()){
					pw.write(fastqIt1.next());
				}
				pw.close();
				}
		
	}



		//	static int id = 
		 File[] inputFile = null;
		 
		
		 public int count=0;
		 File outputBAM;
		public FastqToSAMRecord(String[] input, String mm2Index, int maxReads,Collection<String>readsToInclude, double q_thresh, boolean keepBam) throws IOException{
			this.mm2Index = mm2Index;
			this.q_thresh = q_thresh;
			this.readsToInclude = readsToInclude;
			 max_per_file = maxReads;
			if(max_per_file <0) max_per_file = Integer.MAX_VALUE;
		
			if(input[0].startsWith("ftp://")  || input[0].startsWith("file:/")){
				 this.input = input;
			}else{
				inputFile = new File[input.length];
				this.input = new String[input.length];
				for(int k=0; k<input.length; k++) {
					inputFile[k] = new File(input[k]);;
					this.input[k] = "file:/"+(inputFile[k]).getAbsolutePath();	
				}
			     
			}
			if(keepBam){
				this.outputBAM = new File(inputFile[0]+"."+System.currentTimeMillis()+".sam");

				outputBAM.deleteOnExit();
			
				//SAMFileHeader header  = new SAMFileHeader();
				//header.setSortOrder(SortOrder.unsorted);
				
				this.bfw = // (SAMTextWriter) sfw.makeSAMWriter(header, false, outputBAM);
						new SAMTextWriter(outputBAM);//sfw.makeSAMWriter(header, false, outputBAM);
				this.bfw.setSortOrder(SortOrder.unsorted, false);
				this.log = Log.getInstance(this.getClass());
				this.pl = new ProgressLogger(log);
			}
		 }
		
		
		
		@Override
		public boolean hasNext() {
			// if its null it has not been initialised
			boolean res = this.curr_index< this.input.length || iterator.hasNext();
			if(!res) {
				this.close();
			}
			return res;
		}

		int curr_index=0;
		
		
		
		@Override
		public SAMRecord next() {
			if(iterator !=null && !iterator.hasNext()){
				refresh();
			}
			if(iterator==null && curr_index < this.input.length) {
				init(curr_index);
				curr_index++;
			}
			//System.err.println(iterator.hasNext());
			 SAMRecord nxt =  iterator==null ? null :  iterator.next() ;
		//		System.err.println(nxt.getReadName()+" "+bfw==null);

			 if(this.bfw!=null && nxt!=null){
				 this.cntsAdded++;
				 bfw.addAlignment(nxt);
			 }
			 count++;
			 if(nxt==null) this.close();
			 return nxt;
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
	 
	 public static double getQual(byte[] b){
		// byte[] b = nxt.getBaseQualities();
			double sump = 0;
			for(int i=0; i<b.length; i++){
				sump+= Math.pow(10, -b[i]/10.0);
			}
			sump = sump/(double) b.length;
			double q = -10.0*Math.log10(sump);
			if(Double.isNaN(q)){
				q = -1;
			}
		//	System.err.println(q);
			return q;
	 }
	 public static Iterator<SAMRecord> getFilteredIterator(Iterator<SAMRecord> samIter, Collection<String> reads, int max_reads, double q_thresh, 
			 boolean flipName){
		 return new FilteredIterator(samIter, reads,max_reads, q_thresh, flipName);
	 }
	 
	 public static class FilteredIterator implements Iterator<SAMRecord>{
		 private final Iterator<SAMRecord> samIter;
		 final int max_reads;
		 final Collection<String> reads;
		 SAMRecord nxt;
		 final double qual_thresh;
		 int cnt=0;
		 boolean flipname;
		 public FilteredIterator(Iterator<SAMRecord>sam , Collection<String> reads, int max_reads, double qual_thresh, boolean flipname){
			 this.samIter= sam;
			 this.flipname = flipname;
			 this.reads = reads;
			 this.qual_thresh= qual_thresh;
			 this.max_reads = max_reads;
			 nxt = getNext();
			 if(nxt==null){
				 throw new RuntimeException(" iterator is empty");
			 }
		 }
		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return  nxt!=null && cnt < max_reads;
		}
public SAMRecord next(){
	SAMRecord nxt1 = nxt;
	cnt++;
	nxt = getNext();
	if(flipname){
		String refname = nxt1.getReferenceName();
		String readname = nxt1.getReadName();
		nxt1.setReferenceName(readname);
		nxt1.setReadName(refname);
		
	}
	return nxt1;
}
		
		public SAMRecord getNext() {
			
			while(samIter.hasNext()){
				SAMRecord nxt = samIter.next();
				String nme = nxt.getReadName();
				if(reads==null || reads.contains(nme)){
					byte[] b = nxt.getBaseQualities();
					if(b.length==0 || getQual(b)>=qual_thresh) {
						return nxt;
					}
				}
				//nxt=null;
			}
			return null;
		}
		 
	 }
	 
	 public static class SequentialIterator implements Iterator<SAMRecord> {
		 private final Iterator<SAMRecord>[] samIters;
		 final int len;
		 int curr=0;
		 public  SequentialIterator(Iterator<SAMRecord>[] samIters) {
				this.samIters = samIters;
				len = samIters.length; 
			}
		@Override
		public boolean hasNext() {
			return curr<samIters.length-1 || samIters[curr].hasNext();
		}

		@Override
		public SAMRecord next() {
			if(!samIters[curr].hasNext()){
				curr = curr+1;
				if(curr>=samIters.length) return null;
			}
			SAMRecord sr = samIters[curr].next();
			if(sr!=null) sr.setAttribute(SequenceUtils.src_tag,curr);
			return sr;
		}
		 
	 }
	
	public static class CombinedIterator implements Iterator<SAMRecord> {
		private final Iterator<SAMRecord>[] samIters;
		private final SAMRecord[] currentVals; 
		private final TreeSet<SamR> samR = new TreeSet<SamR>();
		
		final int len;
		int nullcount=0;
		final boolean sorted;
		public  CombinedIterator(Iterator<SAMRecord>[] samIters, boolean sorted) {
			this.samIters = samIters;
			len = samIters.length; 
			this.sorted = sorted;
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
	//		return nullcount<len;
			return samR.size()>0;
		}
		
		SamR previous = null;
		@Override
		public SAMRecord next() {
			SamR  first = samR.first(); //next element
			if(sorted && previous!=null && previous.compareTo(first)>0){
				System.err.println(previous);
				System.err.println(first);
				throw new  RuntimeException("not sorted");
			}
			samR.remove(first);
			int i = first.ind;
			SAMRecord sr = currentVals[i];
			sr.setAttribute(SequenceUtils.src_tag,i);
			if(samIters[i].hasNext()){
				//replace element removed if possible
				currentVals[i] = samIters[i].next();
				while(currentVals[i]!=null && currentVals[i].getReadUnmappedFlag()){
					currentVals[i] = samIters[i].hasNext() ? samIters[i].next() : null;
				}
				if(currentVals[i]!=null) samR.add(new SamR(currentVals[i],i));
				else nullcount++;
			}else{
				System.err.println("iterator"+i+" exhausted "+sr.getReferenceName()+" "+sr.getAlignmentStart());
				nullcount++;
				System.err.println(nullcount+" of "+len+" finished");
				currentVals[i] = null;
			}
			previous = first; 
			return sr;
		}
		
	}
	
	
	
	
	//chroms is string e.g 0:1:2  or 0,0,2400000:1,0,240000
public static Iterator<SAMRecord>  getCombined(Iterator<SAMRecord>[] samReaders, boolean sorted, boolean sequential){
		int len = samReaders.length;
		if(sequential) return new SequentialIterator(samReaders);
		if(samReaders.length==1  && false){
			return samReaders[0];
		}else{
			return 		new CombinedIterator(samReaders, sorted);
		}
		
	}

	
	
	
	
	
public static String minimapIndex(File refFile,   boolean overwrite, boolean saveSeqs) throws IOException, InterruptedException {
	final String mm2 = mm2_path;
	final String mem = mm2_mem;
	if(refFile.isDirectory()) throw new RuntimeException("is directory");
	String infile = refFile.getAbsolutePath();
	File indexFile = new File(infile+".mmi");
	File indexFile1 =  new File(infile.replaceAll(".fasta","").replaceAll(".fa", "")+".mmi");
	if(indexFile.exists() && !overwrite){
		System.out.println("minimap index exists");
		return indexFile.getAbsolutePath();
	}else if(indexFile1.exists()){
		System.out.println("minimap index exists");
		return indexFile1.getAbsolutePath();
	}else{
		System.err.println("making index "+saveSeqs+" "+mem);
		ProcessBuilder pb ;
		if(saveSeqs){
			pb= new ProcessBuilder(mm2, 
					"-I",
					mem,
					"-d",
					indexFile.toString(),
					refFile.toString()
					);
		}else{
		pb = new ProcessBuilder(mm2, 
				"--idx-no-seq",
				"-I",
				mem,
				"-d",
				indexFile.toString(),
				refFile.toString()
				);
		}
		//System.err.println(pb.toString());
			Process p =  pb.redirectError(ProcessBuilder.Redirect.to(new File("err_index.txt"))).start();
			p.waitFor();
	}
	return indexFile.getAbsolutePath();
}

public static Iterator<SAMRecord> getCombined(Iterator<SAMRecord>[] samIters, Collection[] reads, Integer max_reads,
		String chrToInclude, boolean sorted , boolean sequential) {
	final Map<Integer, int[]> chrom_indices_to_include = new HashMap<Integer, int[]>();
	final Set<String> reads_all = new HashSet<String>();
	if(reads==null){
		for(int i=0; i<reads.length; i++){
			for(Iterator<String> it = reads[i].iterator(); it.hasNext();){
				reads_all.add(it.next());
			}
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
	final Iterator<SAMRecord>  sr = getCombined(samIters, sorted, sequential);
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

/*
public static void annotateWithGenomeLength(File refFile, 
		HashMap<String, String> seq2Species, HashMap<String, Integer> seqToLen)  throws NumberFormatException, IOException{
//	File lenF = new File(refFile.getParentFile(),refFile.getName().replaceAll(".gz", "")+".len.txt.gz");
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(refFile))));
		String st = "";
		while((st = br.readLine())!=null){
			String[] str = st.split("\t");
			seqToLen.put(str[2],Integer.parseInt(str[4]));
		//	Node node = this.getNode(str[0]);
		//	if(node!=null){
		//	node.getIdentifier().setAttribute("length",Integer.parseInt(str[1]));
		//	}
		}
		br.close();

}*/

public static Collection<String> getReadList(String readList, boolean split) {
	if(readList==null || readList=="null") return null;
	Set<String> reads = new HashSet<String>();
	try{
	
	InputStream is = new FileInputStream(new File(readList));
	 if(readList.endsWith(".gz")) is = new GZIPInputStream(is);
	 BufferedReader br = new BufferedReader(new InputStreamReader(is));
	 String st = "";
	 while((st = br.readLine())!=null){
		 
		 reads.add(split ? st.split("\\s+")[0] : st);
	 }
	}catch(IOException exc){
		exc.printStackTrace();
	}
	return reads;
}


}
