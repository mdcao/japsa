package japsa.tools.seq;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMRecord;
import japsa.bio.np.RealtimeSpeciesTyping;

public abstract class CachedOutput {
	  public static int max_seqs_per_cluster = Integer.MAX_VALUE;
	  public static boolean writeAssemblyCommand = false;
	  final File outdir;
	  String species ;
	  int seqs_printed=0;
	  boolean separate = false;
	  int total_count=0;
	  boolean print = false;
	  public static int MIN_READ_COUNT=20;
	  public static int buffer = 50;
	//  boolean writeAlignedPortionOnly =false;
	  List<String> nmes = null;
	  List<Integer> lens = null;
	  final boolean alignedOnly;
	 public CachedOutput(File outdir, String species, boolean separateIntoContigs, boolean alignedOnly) {
		 this.alignedOnly = alignedOnly;
	//	 this.writeAlignedPortionOnly = alignedOnly;
		this.species = species.replace('/', '_').replace(" ", "_");
			this.separate = separateIntoContigs;
			this.nmes = new ArrayList<String>();
			 if(separateIntoContigs){
				 this.outdir = new File(outdir, this.species);
			  }else{
				  this.outdir = outdir;
			  }
			//  this.outdir.mkdirs();
			
	}
	public abstract void write(SAMRecord sam, String annotation, RealtimeSpeciesTyping.Interval interval);
	protected abstract String modify(String in);
	 public abstract void close(Map<String, Integer> species2Len);
	 
	 public abstract void  getOutFile(List<File> fi);
	 
	 public abstract int length();

	public void writeAssemblyCommand(Map<String, Integer> species2Len) {
		if(!writeAssemblyCommand) return;
		if(species2Len==null) return;
		if(print){
			try{
				File spec =outdir;
				File out;
				boolean isDir = spec.isDirectory();
				if(isDir){
					out = new File(spec,"flye.sh");
				}else{
					out = new File(spec+"_flye.sh");
				}
			PrintWriter pw  = new PrintWriter(new FileWriter(out));
		for(int i=0; i<this.nmes.size(); i++){
			String nmesi = nmes.get(i);
			String st;
			if(isDir){
				st= "flye -o flye_out --genome-size "+species2Len.get(nmes.get(i))+" --nano-raw "+species+"/"+this.modify(nmesi);
			}
			else {
				st = "flye -o flye_out --genome-size "+species2Len.get(nmes.get(i))+" --nano-raw "+this.modify(species);
			}
			pw.println(st);
		}
		pw.close();
			}catch(IOException exc){
				exc.printStackTrace();
			}
		}
		
	}
	public void close() {
	this.close(null);	// TODO Auto-generated method stub
		if(outdir.listFiles().length==0) outdir.delete();
	}
	public void writeAll(List<SAMRecord> records, List<String> resclasses) {
		if(records.size()>0){
		if(RealtimeSpeciesTyping.alignedOnly){
			for(int i=0; i<records.size(); i++){
				this.write(records.get(i), resclasses.get(i));
			}
		}else{
			StringBuffer sb = new StringBuffer();
			Iterator<String>res = resclasses.iterator();
			while(res.hasNext()) {
				sb.append(res.next());
				if(res.hasNext())sb.append("_");
			}
			this.write(records.get(0), sb.toString());
		}
		}
	}
	public  void write(SAMRecord samRecord, String string) {
		this.write(samRecord, string, null);
	}
	
}
