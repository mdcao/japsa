package japsa.tools.seq;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

/** this enables splitting of output sequences into species specific bams */
public class CachedFastqWriter implements CachedOutput{
	
	//special thread for this class
	//public static final ExecutorService writeCompressDirsExecutor  = Executors.newSingleThreadExecutor();
public static int MIN_READ_COUNT=20;
	//public static 	FastqWriterFactory fqFact = new FastqWriterFactory();

	
	
  public static int max_seqs_per_cluster = Integer.MAX_VALUE;
 private File outdir;
 private  String species ;
  //boolean append;
//  private File outF;
 // int thresh = 10;
  int seqs_printed=0;

 // FastqWriter base = null;
  //Stack<FastqRecord> baseStack = null;

   List<FastqWriter>fqw = null ;
   List<Stack<FastqRecord>>  stack = null;
   List<String> nmes = null;
  boolean separate = false;
  int total_count=0;
  boolean fasta;
  public CachedFastqWriter(File outdir, String species, boolean separateIntoContigs) {
	this.outdir = outdir;
	this.species = species.replace('/', '_');
	this.separate = separateIntoContigs;
			this.fqw = new ArrayList<FastqWriter>();
			this.stack = new ArrayList<Stack<FastqRecord>>();
			this.nmes = new ArrayList<String>();
	}

boolean print = false;
 
  public void write(SAMRecord sam)  {
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
	  String nme = sam.getReadName();
	  String ref = separate  ? sam.getReferenceName() : species;
		FastqRecord repeat =  new FastqRecord(nme,	readSeq,	"",baseQ);

	  FastqWriter fqw_i  ;
	  total_count++;
	  int  index =  this.nmes.indexOf(ref) ;
	 if(index<0){
			  index=nmes.size();
			  nmes.add(ref);
			  if(print){
				  File outdir1 = outdir;
	    		  if(separate){
	    			  outdir1 = new File(outdir, species);
	    		  }
	    		  File outf = new File(outdir1,ref.replace('|', '_')+".fq");
	    		  FastqWriter fqw_j = new BasicFastqWriter(outf);  
	    				  //fqFact.newWriter();
	    		  fqw.add(fqw_j); 
	    		  
			  }else{
				  fqw.add(null);
			  }
			  stack.add(new Stack<FastqRecord>());
		  }
		 fqw_i= this.fqw.get(index);
      if(fqw_i!=null){
    	  fqw_i.write(repeat);
      }
      else{
    	  
	    	  Stack<FastqRecord> st1 = stack.get(index);
	    	  st1.push(repeat);
	    	 
	    	  if(total_count>=MIN_READ_COUNT){ // opening up all the writers once one reaches minimum count
	    		  print = true;
	    		  this.outdir.mkdirs();
	    		  File outdir1 = outdir;
	    		  if(separate){
	    			  outdir1 = new File(outdir, species);
	    			  outdir1.mkdir();
	    		  }
	    		  for(int j=0; j<this.fqw.size(); j++){
	    			 FastqWriter fqw_j =  new BasicFastqWriter(new File(outdir1,ref.replace('|', '_')+".fq"));
		    		  this.fqw.set(j,  fqw_j);
		    		  Stack<FastqRecord> st2 = stack.get(j);
		    		  while(st2.size()>0){
		    			  fqw_j.write(st2.pop());
		    		  }
	    		  }
			}
    	  }
	  }
  

public void close(){
	for(int i=0; i<fqw.size(); i++){
	if(this.fqw.get(i)!=null) {
		fqw.get(i).close();
	}
	}
	
}







  
 
  
  
}
