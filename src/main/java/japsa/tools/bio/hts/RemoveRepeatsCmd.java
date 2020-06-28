package japsa.tools.bio.hts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
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

		addString("repeatsFile", null,  "Name of repeats file", true);
		addString("reference", null, "Name of reference genome",true);
		
		addStdHelp();		
	} 
	public static void main(String [] args) throws IOException, InterruptedException{		 		
		CommandLine cmdLine = new RemoveRepeatsCmd();		
		args = cmdLine.stdParseLine(args);		

		String reference = cmdLine.getStringVal("reference");		
		
		String repFile = cmdLine.getStringVal("repeatsFile");
	
		removeRepeats(repFile, reference);		


		//paramEst(bamFile, reference, qual);
	}
	
	static int threshold = 100;
	static void removeRepeats(String repFile, String refFile) throws IOException{	
		
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());
		Map<String,Sequence> m = new HashMap<String, Sequence>();
		for(int i=0; i<genomes.size(); i++){
			m.put(genomes.get(i).getName(),genomes.get(i));
		}
		//get the first chrom
	//	int currentIndex = 0;
		Sequence chr = null;
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(repFile))));
		List<String> header = Arrays.asList(br.readLine().split("\t"));
		int chr_ind = header.indexOf("chrom");
		int st_ind = header.indexOf("start");
		int end_ind = header.indexOf("end");
		String st = "";
		StringBuffer sb = new StringBuffer();
		PrintWriter indices = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream("conversion.txt.gz")))));
		SequenceOutputStream output = new SequenceOutputStream((new FileOutputStream(new File("compressed.fa.gz" ))));
		int seqlen=0;
		int start =0;
		int end =-1;
		Set<String> done = new HashSet<String>();
		while((st = br.readLine())!=null){
			String[] str = st.split("\t");
			//System.err.println(st);
			String chrom = str[chr_ind];
			if(chr==null || !chr.getName().equals(chrom)){
				if(sb.length()>0){
					Sequence seq = new Sequence(Alphabet.DNA(),sb.toString(), chr.getName());
					seq.writeFasta(output);
					sb.delete(0, sb.length());
					
				}
				if(done.contains(chrom))throw new RuntimeException("already done "+chrom);
				chr = m.remove(chrom);
				System.err.println(chrom);
				done.add(chrom);
				seqlen=0;
				start=0;
			
			}
			end = Integer.parseInt(str[st_ind]); // end of the flanking
			if(end-start > threshold){
				sb.append(chr.subSequence(start, end));
				seqlen += (end-start);
			}
				start = Integer.parseInt(str[end_ind]);
				indices.println(chrom+"\t"+seqlen+"\t"+start+"\t"+(start-seqlen));
		}
		indices.close();
		output.close();

	}
	//chrom   start   end     period  unitNo  size    target  repeatUnit      #H:ID
	//chr1    10000   10468   6       77.2    468     chr1:9000-11468 TAACCC  chr1_10000_10468

}
