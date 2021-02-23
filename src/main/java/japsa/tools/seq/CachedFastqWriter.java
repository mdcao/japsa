package japsa.tools.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Stack;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqEncoder;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

/** this enables splitting of output sequences into species specific bams */
public class CachedFastqWriter {
	
	//special thread for this class
	public static final ExecutorService writeCompressDirsExecutor  = Executors.newSingleThreadExecutor();

	static class FQWriter implements FastqWriter{
		private  PrintStream writer;
		public FQWriter(File outdir, String target ,boolean append) {
			try{
				File out = new File(outdir,target.replace('/', '_'));
				//System.err.println(out);
				OutputStream os = new FileOutputStream(out, append);
				if(target.endsWith(".gz")) os = new GZIPOutputStream(os);
				writer = new PrintStream(os);
			 
			}catch(Exception exc){
				exc.printStackTrace();
				
				System.exit(0);
			}
			
		}

		public FQWriter(FileOutputStream os) {
			writer = new PrintStream(os);
		}

		@Override
		public void close() {
			writer.println();
			writer.close();
			
		}

		@Override
		public void write(FastqRecord rec) {
			// TODO Auto-generated method stub
			 FastqEncoder.write(writer, rec);
		}
		
	}
	
	
  public static int max_seqs_per_cluster = Integer.MAX_VALUE;
  private File outdir;
  private  String species ;
  //boolean append;
  int thresh = 10;
  int seqs_printed=0;
  boolean append = false; // starts out as false, but then set to true after first batch is written
  public CachedFastqWriter(File outdir, String target) {
	this.species=target;
	this.outdir = outdir;
	}


   boolean lock = false;  // is stream locked for closing
public void printAll(){
	lock = true; 
	
	Runnable run = new Runnable(){
		@Override
		public void run() {
			//File out = new File()
			
			 FQWriter	so = new FQWriter(outdir, species+".fastq.gz", append);
 			while(stack.size()>0 && seqs_printed < max_seqs_per_cluster){
				 so.write(stack.pop());
					seqs_printed++;
			}
 			stack.clear();
			 so.close();
			lock = false;
			
		}
		
	};
	try{
		run.run();
//	writeCompressDirsExecutor.execute(run);
	}catch(Exception exc){
		exc.printStackTrace();
	}
	append=true;
}

  
  public void write(SAMRecord sam) throws IOException {
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
		FastqRecord repeat =  new FastqRecord(sam.getReadName(),	readSeq,	"",baseQ);
	  if(seqs_printed < max_seqs_per_cluster){
		stack.push(repeat);
		if(stack.size()==thresh){
			this.printAll();
		}
	  }
	}

Stack<FastqRecord>  stack= new Stack<FastqRecord>();


public void close(){
	if(stack.size()>0){
	this.printAll();
	}
	
}

public String entryname() {
	// TODO Auto-generated method stub
	return species;
}






  
 
  
  
}
