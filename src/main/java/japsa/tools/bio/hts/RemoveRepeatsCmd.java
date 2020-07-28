package japsa.tools.bio.hts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
/**
 * @author lachlancoin
 *
 */
@Deployable(
	scriptName = "jsa.hts.removeRepeats",
	scriptDesc = "Detecting repeats in long read sequencing data")
public class RemoveRepeatsCmd extends CommandLine {
	public RemoveRepeatsCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("resDir", "./results",  "Results directory", false);
		addString("excl", null,  "File with lines to exclude", false);
		addString("repeatsFile", null,  "Name of repeats file", false);
		addString("repeatsHeader", "chrom:chromStart:chromEnd:sequence:period",  "repeats header", false);
		addString("reference", null, "Name of reference genome",true);
		addString("chroms", "chr1:chr2:chr3:chr4:chr5:chr6:chr7:chr8:chr9:chr10:chr11:chr12:chr13:chr14:chr15:chr16:chr17:chr18:chr19:chr20:chr21:chr22:chrX:chrY:chrM", "chroms to retain",false);
		addInt("nrep",50,"Length of Ns to remove", true);
		addStdHelp();		
	} 
	static class Indel{
		int start;
		//String seq;
		int len;
		public String seq;
		Indel(int start, String seq){
			this.start = start;
			//this.seq = seq;
			this.len =  seq.length();
			this.seq = seq.substring(0,Math.min(seq.length(), 10));
		}
		public String toString(){
			//if(len==1) return start+"";
			//else
			return start+"-"+(start+len);
		}
		public int len() {
			// TODO Auto-generated method stub
			return len;
		}
		public int end() {
			// TODO Auto-generated method stub
			return start + len;
		}
	}
	static int nrep =50;
	static Pattern[] p = null;
	static void findPatt(Sequence chr, SortedMap<Integer, Indel> allMatches){
		if(nrep==0) return;
		String seqr = chr.toString();
		for(int i=0; i<p.length; i++){
			Matcher m = p[i].matcher(seqr);
			int tot =0;
			int prev =0;
			while (m.find()) {
			
				prev=m.start();
				int len =  m.end()-m.start();
				tot+=len;
				String subst = seqr.substring(m.start(), m.end());
				allMatches.put(m.start(),new Indel(m.start(),subst));
			}
			System.err.println(chr.getName()+" : "+i+" "+allMatches.size());//+" "+allMatches.firstKey()+" "+allMatches.lastKey());

		}
	}
	
	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new RemoveRepeatsCmd();		
		args = cmdLine.stdParseLine(args);		
		String reference = cmdLine.getStringVal("reference");		
		String repFile = cmdLine.getStringVal("repeatsFile");
		String resDir = cmdLine.getStringVal("resDir");
		String excl = cmdLine.getStringVal("excl");
		File exclFile = excl==null ? null: new File (excl);
		File resD =(new File(resDir)); 
		(resD).mkdir();
		RemoveRepeatsCmd.nrep = cmdLine.getIntVal("nrep");
		RemoveRepeatsCmd.mem = (Runtime.getRuntime().maxMemory()-1000000000)+"";
		System.err.println(RemoveRepeatsCmd.mem);
		Pattern p = Pattern.compile("N{"+nrep+",}");
		//Pattern p1 = Pattern.compile("\\p{javaLowerCase}{"+nrep+",}");
		List<String> chroms = Arrays.asList(cmdLine.getStringVal("chroms").split(":"));
		System.err.println(chroms);
		 RemoveRepeatsCmd.p = new Pattern[] {p};//,p1};
		 ArrayList<Sequence> genomes = SequenceReader.readAll(reference, Alphabet.DNA());
			String out = reference.replace(".fasta", "").replace(".fa", "").replace(".gz", "");
			//out = resDir+"/"+out;
		 for(int i=genomes.size()-1; i>=0; i--){
			 String nme = genomes.get(i).getName();
			 if(chroms!=null && chroms.indexOf(nme)<0 ){
				 genomes.remove(i);
				 System.err.println("removing "+nme);
			 }
		 }
		 RemoveRepeatsCmd.headers = cmdLine.getStringVal("repeatsHeader").split(":");
				 //"chrom:chromStart:chromEnd:sequence:period".split(":");
		// RemoveRepeatsCmd.headers="chrom:start:end:repeatUnit:period".split(":");//  unitNo  size    target  repeatUnit      #H:ID
		 
		System.err.println(genomes.size());
		RemoveRepeatsCmd.Inner inner = new RemoveRepeatsCmd.Inner(genomes, repFile==null ? null : new File(repFile), out, resD, exclFile);
		 if(repFile==null){
			 inner.removeN();
		 }else{
			 inner.removeRepeats();		
		 }
		
		
		//paramEst(bamFile, reference, qual);
	}
	static String mm2="/home/lachlan/github/minimap2/minimap2";
	//static int mm2Threads = 2;
	static String mem = "5g";
	static class SeqWriter{
		OutputStreamWriter os;
		PrintWriter coord_pw = null;
		File outf;
		final boolean mmi;
		public SeqWriter(String string, boolean mmi) throws FileNotFoundException, IOException {
			outf = new File( string);
			this.mmi = mmi;
			os = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf)));
			// coord_pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(string+"coord.bed"))));

		}
		public void writeln(String st) throws IOException{
			os.write(st);os.write("\n");
			if(coord_pw!=null) coord_pw.println(st);
		}
		public void close() throws IOException, InterruptedException{
			os.close();
			if(coord_pw!=null) coord_pw.close();
			//"$mm2 -I 8g -d $out $fa" 
			if(mmi){
			ProcessBuilder pb = new ProcessBuilder(mm2, 
					"-I",
					mem,
					"-d",
					outf.toString()+".mmi",
					outf.toString()
					);
			//System.err.println(pb.toString());
				Process p =  pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
				p.waitFor();
			}
			//return SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));
		}
		public void write(Sequence chr, int start, int end) throws IOException{
			this.write(chr, start, end, false, null); //not checking for Ns in repeat
		}
		public void write(Sequence chr, int start, int end, boolean check, String repeat_unit) throws IOException {
			if(coord_pw!=null) coord_pw.println(chr.getName()+"\t"+start+"\t"+end);
			String st = chr.subSequence(start, end).toString();
				
				if(repeat_unit!=null) {
					String nme = chr.getName()+"_"+start+"_"+end;
					os.write(">"+nme);
					os.write(" ");
					boolean prefix = st.startsWith(repeat_unit);
					os.write(prefix+" "+repeat_unit);
					os.write("\n");
				}
			if(check){
				Matcher m = p[0].matcher(st);
				if(m.matches()){
					System.err.println(start+"-"+end);
					throw new RuntimeException("should not match "+st);
				}
			}
		
			os.write(st);
			os.write("\n");

		}
	}
	static int threshold = 20;
	static String[] headers = "chrom:chromStart:chromEnd:sequence:period".split(":");
//	#bin    chrom   chromStart      chromEnd        name    period  copyNum consensusSize   perMatch        perIndel        score   A       C       G       T       entropy s
  static class Inner{
	
	   //just remove the Ns
	void removeN() throws IOException, InterruptedException{
		for(int i=0; i<genomes.size(); i++){
			Sequence seq = genomes.get(i);
			String chrom = seq.getName();
			String outf = resDir.getAbsolutePath()+"/"+chrom+"."+out+".no_repeats.fa.gz";
			String outf1 = resDir.getAbsolutePath()+"/"+chrom+"."+out+".repeats.fa.gz";
			indices = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(resDir.getAbsolutePath()+"/"+chrom+"."+out+".bed.gz")))));
			output =new SeqWriter(outf, true);
			output_repeats = new SeqWriter(outf1, false);
		//	System.err.println(i);
			System.err.println(seq.getName());
			allMatches.clear();
			findPatt(seq, allMatches);
			output.writeln(">"+seq.getName());
			 processN(allMatches,seq, output, output_repeats ,  indices,header_out,  0, 0);
			 indices.close();
				output.close();
				output_repeats.close();
		}
		
	}
	
	public Inner(ArrayList<Sequence> genomes, File repFile, String out, File resDir, File to_excl) throws FileNotFoundException, IOException {
		// TODO Auto-generated constructor stub
		for(int i=0; i<genomes.size(); i++){
			String chr = genomes.get(i).getName();
			m.put(chr,genomes.get(i));
			this.chrom_index.put(chr, i);
		}
		this.resDir = resDir;
		if(to_excl!=null && to_excl.exists()){
			BufferedReader br1= new BufferedReader(new InputStreamReader((new FileInputStream(to_excl))));
			String st = "";
			while((st = br1.readLine())!=null){
				this.excl.add(st);
			}
			br1.close();
		}
		this.out = out;
		this.genomes = genomes;
		chroms = new HashSet<String>(m.keySet());
		System.err.println(chroms);
		if(repFile!=null){
			br= new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(repFile))));
			String[] header_st = br.readLine().split("\t");
			List<String> header = Arrays.asList(header_st);
			 chr_ind = header.indexOf(headers[0]);
			 st_ind = header.indexOf(headers[1]);
			end_ind = header.indexOf(headers[2]);
			unit_ind = header.indexOf(headers[3]);
			period_ind = header.indexOf(headers[4]);
		}else{
			br=null;
			chr_ind = -1;
			st_ind = -1;
			end_ind= -1;
			unit_ind = -1;
			period_ind = -1;
		}
		out_st_r= new String[header_out.length];
		 out_st = new String[header_out.length];
		 Arrays.fill(out_st_r, "");
	}
	final Map<String,Sequence> m = new HashMap<String, Sequence>();
	final Map<String, Integer> chrom_index = new HashMap<String, Integer>();
	final ArrayList<Sequence> genomes;
	final Set<String> chroms;
	final BufferedReader br ;
	final String out;
	final File resDir; 
	String[] header_out = "chrom:start:end:start1:end1:offset:repeatlen:ID:sequence".split(":");
	final int chr_ind, st_ind, end_ind, unit_ind, period_ind;
	 PrintWriter indices; 
	 SeqWriter output, output_repeats;
	String[] out_st_r , out_st;
	TreeMap<Integer, Indel> allMatches = new TreeMap<Integer, Indel>();
	
	private void finishChrom(Sequence chr, int start, int end, int seqlen) throws IOException, InterruptedException{
		System.err.print("finishing.."+chr.getName());
		end = chr.length()+1;
		start = processN(allMatches.headMap(end),chr, output,output_repeats, indices, out_st_r, start, seqlen);
		if(end-start > threshold){
				output.write(chr, start, end) ;
			seqlen += (end-start);
		}
		indices.close();
		output.close();
		output_repeats.close();
		
		
	}
	
	
	
Set<String> excl = new HashSet<String>();
	
	 void removeRepeats() throws IOException, InterruptedException{	
		 Sequence chr = null;
		String st = "";
		int seqlen=0;
		int start =0;
		int end =-1;
		int period = 0;
		Set<String> done = new HashSet<String>();
		 // allMatches relates to NNNN sequence
		outer: while((st = br.readLine())!=null){
			if(excl.contains(st)) {
				System.err.println("excluding "+st);
				continue outer;
			}
			String[] str = st.split("\t");
			String chrom = str[chr_ind];
		
			if(!chroms.contains(chrom)) {
				continue outer;
			}
			if(chr==null || !chr.getName().equals(chrom)){
				if(chr!=null){
					finishChrom(chr, start, end, seqlen);
				}
				if(done.contains(chrom))throw new RuntimeException("already done "+chrom);
				int chromi = chrom_index.get(chrom);
				String outf = resDir.getAbsolutePath()+"/"+chrom+"."+out+".no_repeats.fa.gz";
				String outf1 = resDir.getAbsolutePath()+"/"+chrom+"."+out+".repeats.fa.gz";
				indices = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(resDir.getAbsolutePath()+"/"+chrom+"."+out+".bed.gz")))));
				output =new SeqWriter(outf, true);
				output_repeats = new SeqWriter(outf1, false);
				System.err.println("output : ");
				System.err.println(outf+"\n"+outf1);
				print(header_out, indices); indices.println();
				System.err.println(chrom);
				chr = m.remove(chrom);
				output.writeln(">"+chrom);
				allMatches.clear();
				findPatt(chr, allMatches);
				
				System.err.println(allMatches);
				done.add(chrom);
				seqlen=0;
				start=0;
			}
			if(chr!=null) output_repeats.write(chr, Integer.parseInt(str[st_ind]), Integer.parseInt(str[end_ind]), false, unit_ind<0 ? "NA" : str[unit_ind]);
			period = Integer.parseInt(str[period_ind]); //include one period of the repeat in the flanking

			end = Integer.parseInt(str[st_ind]) + period; // end of the flanking region  (start of repeat + 1 repeat unit)
			SortedMap<Integer, Indel> preceding = allMatches.headMap(end+1);
			if(preceding.size()>0){
				///if there is intervening NNN sequence we have to deal with that first
				start = processN(preceding,chr, output,output_repeats, indices,out_st_r,  start, seqlen);
			}
			
			if(end-start > threshold){ //only write intervening flanking (ie between repeats) if it at least 100bp
				if(chr!=null) output.write(chr, start, end);//.subSequence(start, end).toString()); output.write("\n");
				seqlen += (end-start);
			}
			out_st[0] = chrom;
			out_st[1]  = seqlen+"";
			out_st[2] = seqlen+"";
		//	out_st[3] = str[st_ind];
			out_st[3] = end+"";
			out_st[4] = str[end_ind];
			out_st[6] = period+"";
			out_st[7] = chrom+"_"+str[st_ind]+"_"+str[end_ind];
			out_st[8] = unit_ind<0 ? "NA":  str[unit_ind];
			out_st[5] = ""+(end - seqlen);
	
			
				print(out_st, indices); indices.println();
				start = Math.max(Integer.parseInt(str[end_ind]), start) ; ///start of next flanking
		}
		if(chr!=null){
			finishChrom(chr, start, end, seqlen);
		}
	}
	
  }
	
	static int processN(SortedMap<Integer, Indel> preceding, Sequence chr, SeqWriter os, SeqWriter out_repeats, PrintWriter indices,String[] out_st_r, int start, int seqlen  
			) throws IOException {
		
		if(preceding.size()==0) return start;
		String chrom = chr.getName();
		List<Integer> l = new ArrayList<Integer> (preceding.keySet());
		Iterator<Integer> it = l.iterator();
	//	
		while(it.hasNext()){
			Integer key = it.next();
			Indel indel = preceding.remove(key);
			int end1 = indel.start;
			int offset = start - seqlen;
			if(end1 - start> threshold){
				
				os.write(chr, start, end1);//.subSequence(start, end1).toString());
				//os.write("\n");
				seqlen += (end1-start);
			}
			
			out_st_r[0] = chrom;
			out_st_r[1]  = seqlen+"";
			out_st_r[2] = seqlen+"";
			out_st_r[3] = indel.start+"";
			out_st_r[4] = indel.end()+"";
			out_st_r[5] = ""+offset;
			out_st_r[6] = indel.len()+"";
			out_st_r[7] = chrom+"_"+indel.start+"_"+indel.end();
			out_st_r[8] = "N";
			out_repeats.write(chr, indel.start, indel.end(), false, indel.seq);
			print(out_st_r, indices); indices.println();
			start = indel.end();
		}
		return start;
		
	}

	
	static void print(String[] str, PrintWriter pw){
		for(int i=0; i<str.length; i++){
			pw.print(str[i]);
			if(i<str.length-1)pw.print("\t");
		}
	}

}
