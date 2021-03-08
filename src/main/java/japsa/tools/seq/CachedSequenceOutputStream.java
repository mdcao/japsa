package japsa.tools.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/** this enables splitting of output sequences into species specific bams */
public class CachedSequenceOutputStream extends CachedOutput {
	 public  void  getOutFile(List<String> fi){
		  for(int i=0; i<l.size(); i++){
			  if(l.get(i).printed>0) fi.add(outdir+"/"+ l.get(i).nme1);
		  }
	  }
	 public int length(){
		 return l.size();
	 }
	class Inner{
		int printed=0;
		  SequenceOutputStream fqw_os;
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
			  try{
				  if(print){
					  if(fqw_os==null && stack.size()>0 ) {
						  fqw_os =  new SequenceOutputStream(new FileOutputStream((new File(outdir, nme1))));  
					  }
					
						  while(stack.size()>0){
							  printed++;
							 ((Sequence)stack.pop()).writeFasta(fqw_os);
						  }
					
				  }
			  }catch(IOException exc){
				  exc.printStackTrace();
			  }
		  }
		  public void close(){
			  try{
			  if(print){
				  this.clear();
				  if(fqw_os!=null) this.fqw_os.close();
			  }else{
				  stack.removeAllElements();
			  }
			  }catch(IOException exc){
				  exc.printStackTrace();
			  }
		  }
		  
	  }
	  List<Inner> l = new ArrayList<Inner>();
	  final  Inner remainder;
	
	  public CachedSequenceOutputStream(File outdir, String species, boolean separateIntoContigs, boolean alignedOnly) {
		super(outdir, species, separateIntoContigs, alignedOnly);
		this.l = new ArrayList<Inner>();
		this.remainder = new Inner("remainder.fa");
	}

 static Alphabet alpha = Alphabet.getAlphabet("DNA");
 protected   String modify(String ref){
	 return ref.replace('|', '_')+".fa";
 }

 public void write(SAMRecord sam, String annotation)  {
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
	  String nme = sam.getReadName();
	  if(writeAlignedPortionOnly) {
		  int st = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart());
			int end = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd());
			 if(st > 100){
		  		 this.remainder.push(new Sequence(alpha, readSeq.substring(0,st-1), nme+".L."+st));
		  }
		  if(end < sam.getReadLength()-100){
			  	this.remainder.push(new Sequence(alpha, readSeq.substring(end,sam.getReadLength()), nme+".R."+end));
		  }
		  readSeq =  readSeq.substring(st-1,end); // because sam is 1-based
	  		nme = nme+" "+st+"-"+end+" "+sam.isSecondaryOrSupplementary();
 	  }
	  String ref = separate  ? sam.getReferenceName() : species;
	  Sequence repeat =  new Sequence(alpha, readSeq,	nme+"_"+annotation);
	  total_count++;
	//  System.err.println(total_count);
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
