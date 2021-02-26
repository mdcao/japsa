package japsa.tools.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import japsa.seq.SequenceOutputStream;

/** this enables splitting of output sequences into species specific bams */
public class CachedFastqWriter extends CachedOutput{
  class Inner{
	  FastqWriter fqw;
	  Stack stack = new Stack();
	  final String nme1;
	 
	  Inner(String nme1){
		  this.nme1 = nme1;
	  }
	  public void push(Object fqw){
		  stack.push(fqw);
		  if(print ){
			  clear();
		  }
	  }
	  public void clear(){
		  if(fqw==null && stack.size()>0 && print) fqw =  new BasicFastqWriter(new File(outdir, nme1));
		  if(print){
		  while(stack.size()>0){
			  fqw.write((FastqRecord)stack.pop());
		  }
		  }
	  }
	  public void close(){
		  if(print){
			  this.clear();
			  if(fqw!=null) this.fqw.close();
		  }else{
			  stack.removeAllElements();
		  }
	  }
	  
  }
  List<Inner> l = new ArrayList<Inner>();
  final  Inner remainder;
 
  public CachedFastqWriter(File outdir, String species, boolean separateIntoContigs, boolean alignedOnly) {
	  super(outdir, species, separateIntoContigs, alignedOnly);
	  this.l = new ArrayList<Inner>();
	  this.remainder = new Inner("remainder.fq");
	}

    protected  String modify(String ref){
		 return ref.replace('|', '_')+".fq";
	 }

  public void write(SAMRecord sam)  {
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
	  String nme = sam.getReadName();
	  if(writeAlignedPortionOnly) {
		  int st = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			int end = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			if(st > 100){
		  		 this.remainder.push(new FastqRecord(nme+".L."+st, readSeq.substring(0,st-1),"", baseQ.substring(0,st-1) ));
		  }
		  if(end < sam.getReadLength()-100){
		  		 this.remainder.push(new FastqRecord(nme+".R."+end, readSeq.substring(end,sam.getReadLength()),"", baseQ.substring(end,sam.getReadLength()) ));
		  }
			//char strand = sam.getReadNegativeStrandFlag() ? '-' : '+';
	  	readSeq =  readSeq.substring(st-1,end); // because sam is 1-based
	  	nme = nme+" "+st+"-"+end+" "+sam.isSecondaryOrSupplementary();
  	  }
	  String ref = separate  ? sam.getReferenceName() : species;
	  FastqRecord repeat =  new FastqRecord(nme,	readSeq,	"",baseQ);
	  total_count++;
	  if(! print && total_count> MIN_READ_COUNT) {
		  print = true;
		  this.outdir.mkdirs();
	  }
	  int  index =  this.nmes.indexOf(ref) ;
	  if(index<0){
			  index=nmes.size();
			  nmes.add(ref);
			  this.l.add(new Inner(this.modify(ref)));
		  }
	 	 Inner inner = l.get(index);
	 	  inner.push(repeat);
		 
	  }
  

public void close(Map<String, Integer> species2Len){
	
	for(int i=0; i<l.size(); i++){
		l.get(i).close();
	}
	this.remainder.close();
	super.writeAssemblyCommand(species2Len);
	
}







  
 
  
  
}
