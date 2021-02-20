package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class RemoveReadsFromFastq {

	public static void  main(String[] args){
		try{
			File in = new File("combined_nohuman.fastq.gz");
			File torem = new File("to_remove.txt");
			RemoveReadsFromFastq  rr = new RemoveReadsFromFastq(torem);
			File out = new File("fixed_nohuman.fastq");
			rr.filter(in, out);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

static BufferedReader getBR(File file)throws IOException{
	
	
	 BufferedReader br;
		if(file.getName().endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}
		else{
			br = new BufferedReader(new FileReader(file));
		}
		return br;
}
	Set<String> toremove = new HashSet<String>();
	
	RemoveReadsFromFastq(File toremovef) throws IOException {
		BufferedReader br =getBR(toremovef);
		String st = "";
		while((st = br.readLine())!=null){
			toremove.add("@"+st);
		}
		br.close();
	}
	
	//String torepl = "2018-12";
	//String torepl1 = "2018-03";
	
	//String runid="ff2ff8cf910b0531d154a48f1d0ea8a9667e5462";
	//String start_time="2018-03-26T17:23:19Z";
	
	//2018-12-06
	//2018-03-26  
	
	public void filter(File in, File out) throws IOException{
		BufferedReader br =getBR(in);
		PrintWriter pw = new PrintWriter(new FileWriter(out));
		PrintWriter removed = new PrintWriter(new FileWriter("removed"));
		PrintWriter kept = new PrintWriter(new FileWriter("kept"));
		String st1, st2, st3, st4;
		st1 = "";
		boolean cts = this.toremove.contains("@f0612e2f-25c4-47dc-b79c-d5ab26c9c2a1");
		while((st1 = br.readLine())!=null){
			String[] str = st1.split("\\s+");
			st2 = br.readLine();
			st3 = br.readLine(); 
			st4 = br.readLine();
			System.err.println(str[0]);
			if(toremove.contains(str[0])){
				System.err.println("removing");
				removed.println(str[0]);
			}else{
				kept.println(str[0]);
				pw.println(st1);
				pw.println(st2);
				pw.println(st3);
				pw.println(st4);
			}
			
		}
		br.close();
		pw.close();
		kept.close();
		removed.close();
	}
}
