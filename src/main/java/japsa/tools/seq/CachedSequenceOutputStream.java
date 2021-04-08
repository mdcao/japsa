package japsa.tools.seq;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import htsjdk.samtools.SAMRecord;
import japsa.bio.np.RealtimeSpeciesTyping;
import japsa.bio.np.RealtimeSpeciesTyping.Interval;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;

/** this enables splitting of output sequences into species specific bams */
public class CachedSequenceOutputStream extends CachedOutput {
	
	public static int remainderThresh = Integer.MAX_VALUE;// limit to write remainder sequence
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
			  if(print  && stack.size()>=buffer){
				  clear();
			  }
		  }
		  public void clear(){
			  try{
				  if(print && (stack.size()+printed) >= MIN_READ_COUNT){
					  if(fqw_os==null && stack.size()>0 ) {
						  fqw_os =  new SequenceOutputStream(new FileOutputStream((new File(outdir, nme1+"_"+stack.size()))));  
					  }
					
						  while(stack.size()>0){
							  printed++;
							  Sequence seq = (Sequence)stack.pop();
							  try{
							 (seq).writeFasta(fqw_os);
							  }catch(Exception exc){
								  System.err.println("problem "+seq.getName());
								  System.err.println(seq.toString());
								  exc.printStackTrace();
							  }
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
	
	  public CachedSequenceOutputStream(File outdir, String species, boolean separateIntoContig, boolean alignedOnly) {
		super(outdir, species, separateIntoContig, alignedOnly);
		this.l = new ArrayList<Inner>();
		this.remainder = new Inner("remainder.fa");
	}

 static Alphabet alpha = Alphabet.getAlphabet("DNA");
 protected   String modify(String ref){
	 return ref.replace('|', '_')+".fa";
 }


 // interval is target interval in reference
 public void write(SAMRecord sam, String annotation, RealtimeSpeciesTyping.Interval interval)  {
	 boolean primary  = !sam.isSecondaryOrSupplementary();
	 
	 
	  String baseQ = sam.getBaseQualityString();
	  String readSeq = sam.getReadString();
	  String nme = sam.getReadName();
	  if(nme.length()>100){
		  System.err.println("h");
	  }
	  int stA = interval==null ? sam.getAlignmentStart() : Math.max(interval.start, sam.getAlignmentStart());
	  int endA = interval==null ? sam.getAlignmentEnd() : Math.min(interval.end, sam.getAlignmentEnd());
	  int lenA = endA - stA;
	  int st = sam.getReadPositionAtReferencePosition(stA, true);
		int end = sam.getReadPositionAtReferencePosition(endA, true);
		int len = end - st;
		  String desc = sam.getReferenceName()+":"+stA+","+endA;//

	  if(alignedOnly) {
			 if(primary && st > remainderThresh){
		  		 this.remainder.push(new Sequence(alpha, readSeq.substring(0,st-1), nme+".L."+st));
			 }
			  if(primary && end < sam.getReadLength()-remainderThresh){
				  	this.remainder.push(new Sequence(alpha, readSeq.substring(end,sam.getReadLength()), nme+".R."+end));
			  }
			  if(st <0 || end-1 >=readSeq.length()){
				  throw new RuntimeException("out of bounds");
			  }
		  readSeq =  readSeq.substring(st,end); // because sam is 1-based
	  		nme = nme+"."+st+"_"+end;
	 
 	  }else{
 		  desc = desc +" read:"+st+","+end;
 	  }
	  String ref = separate  ? sam.getReferenceName() : species;
	  if(interval!=null) ref = ref+"."+interval.start+"-"+interval.end+"."+interval.coverage;
	  Sequence repeat =  new Sequence(alpha, readSeq,	nme+"_"+annotation);
	  repeat.setDesc(desc);
	  total_count++;
	//  System.err.println(total_count);
	  if(! print && total_count>= MIN_READ_COUNT) {
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
//	super.writeAssemblyCommand(species2Len);
	
}








  
 
  
  
}
