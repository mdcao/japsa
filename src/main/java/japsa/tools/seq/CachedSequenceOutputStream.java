package japsa.tools.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/** this enables splitting of output sequences into species specific bams */
public class CachedSequenceOutputStream {
	
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
   List<SequenceOutputStream>fqw = null ;
   List<Stack<Sequence>>  stack = null;
   List<String> nmes = null;
  boolean separate = false;
  int total_count=0;
  boolean fasta;
  public CachedSequenceOutputStream(File outdir, String species, boolean separateIntoContigs) {
	this.outdir = outdir;
	this.species = species.replace('/', '_');
	this.separate = separateIntoContigs;
			this.fqw = new ArrayList<SequenceOutputStream>();
			this.stack = new ArrayList<Stack<Sequence>>();
			this.nmes = new ArrayList<String>();
	}

boolean print = false;
 static Alphabet alpha = Alphabet.getAlphabet("DNA");
  public void write(SAMRecord sam)  {
	  try{
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
	  String nme = sam.getReadName();
	  String ref = separate  ? sam.getReferenceName() : species;
		Sequence repeat =  new Sequence(alpha, readSeq,	nme);

	  SequenceOutputStream fqw_i  ;
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
	    		  SequenceOutputStream fqw_j =  new SequenceOutputStream(new FileOutputStream(new File(outdir1,ref.replace('|', '_')+".fa")));
	    		  fqw.add(fqw_j); 
	    		  
			  }else{
				  fqw.add(null);
			  }
			  stack.add(new Stack<Sequence>());
		  }
		 fqw_i= this.fqw.get(index);
      if(fqw_i!=null){
    	  repeat.writeFasta(fqw_i);
    	 // fqw_i.write(repeat);
      }
      else{
    	  
	    	  Stack<Sequence> st1 = stack.get(index);
	    	  st1.push(repeat);
	    	 
	    	  if(total_count>=MIN_READ_COUNT){ // opening up all the writers once one reaches minimum count
	    		  print = true;
	    		  File outdir1 = outdir;
	    		  if(separate){
	    			  outdir1 = new File(outdir, species);
	    			  outdir1.mkdir();
	    		  }
	    		  for(int j=0; j<this.fqw.size(); j++){
	    			 SequenceOutputStream fqw_j =  new SequenceOutputStream(new FileOutputStream(new File(outdir1,ref.replace('|', '_')+".fa")));
		    		  this.fqw.set(j,  fqw_j);
		    		  Stack<Sequence> st2 = stack.get(j);
		    		  while(st2.size()>0){
		    			  st2.pop().writeFasta(fqw_j);
		    		  }
	    		  }
			}
    	  }
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	  }
  

public void close(){
	for(int i=0; i<fqw.size(); i++){
		try{
	if(this.fqw.get(i)!=null) {
		fqw.get(i).close();
	}
		}catch(IOException exc){
			exc.printStackTrace();
		}
	}
	
}







  
 
  
  
}
