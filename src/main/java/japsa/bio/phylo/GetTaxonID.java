package japsa.bio.phylo;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.io.input.ReversedLinesFileReader;
/** extract taxon ids matching speciesIndex from a list of assembly summary files */
public class GetTaxonID {
  Set<Integer> taxon_set = new HashSet<Integer>();
	Map<String, Integer> name2Taxa = new HashMap<String, Integer>();
	Map<String, Integer> name2Taxa4 = new HashMap<String, Integer>();
	Map<String, Integer> name2Taxa3 = new HashMap<String, Integer>();
	Map<String, Integer> name2Taxa2 = new HashMap<String, Integer>();



	Map<Integer, String> taxa2Sci = new HashMap<Integer, String>();
  public GetTaxonID(){
  }
  
  public String getSciName(final String specName){
	  Integer taxa = getTaxa(specName);
	  if(taxa!=null){
			 return taxa2Sci.get(taxa);
		 }
	  else return null;
  }
  
  public Integer getTaxa(final String specName){
	  String slug =  Slug.toSlug(specName, "");
		String slug4 =  Slug.toSlug(specName, 4,"");
		String slug3 =  Slug.toSlug(specName, 3,"");
		String slug2 =  Slug.toSlug(specName, 2,"");
	 Integer taxa = this.name2Taxa.get(slug);
	 if(taxa==null) taxa = name2Taxa4.get(slug4);
	 if(taxa==null) taxa = name2Taxa3.get(slug3);
	 if(taxa==null) taxa = name2Taxa2.get(slug2);
	 //if(specName.indexOf("229E-related")>=0){
	//	  System.err.println('h');
	 // }
	// err.println(specName+"->"+slug1+"->"+slug2+"->"+slug3+"->"+taxa);
	 return taxa;
  }
  void putTaxa(String nme, Integer taxa){
	 
	  String slug = Slug.toSlug(nme, "");
	  String slug4 = Slug.toSlug(nme, 4,"");
	  String slug3 = Slug.toSlug(nme, 3,"");
	  String slug2 = Slug.toSlug(nme, 2,"");
	  name2Taxa.put(slug, taxa);
	  name2Taxa4.put(slug4, taxa);
	  name2Taxa3.put(slug3, taxa);
	  name2Taxa2.put(slug2, taxa);
	 
  }
  PrintWriter err;
  
  public Map<Integer, Integer> nodeToParent = new HashMap<Integer,Integer >();
  
  public void addNodeDmp(File file) throws IOException{
	  BufferedReader br = getBR(file);
	  String st = "";
	  while((st = br.readLine())!=null){
		  String[] str = st.split("\\|");
			 nodeToParent.put(Integer.parseInt(str[0].trim()), Integer.parseInt(str[1].trim()));
		 
	  }
	  br.close();
  }
  
 
  private File  expand(Set<Integer> taxon_set2, File node_dmp, Set<Integer> set) throws IOException {
	  Set<Integer> todo = new HashSet<Integer>();
	todo.addAll(set);
	 // new in this run
	Set<Integer> done = new HashSet<Integer>();
	//newS.addAll(todo);
	Set<Integer> notP = new HashSet<Integer>();
	String st = "";
	String tme = System.currentTimeMillis()+"";
	File node_dmp1 = new File(node_dmp.getAbsolutePath()+"."+tme+".1");
	PrintWriter pw1 = new PrintWriter(new FileWriter(node_dmp1));
	//String prev = "";
	int ps = todo.size();
	for(int i=0; todo.size()>0; i++){
		//todo.clear();
		//todo.addAll(newS);
		Set<Integer> newS = new HashSet<Integer>();
		System.err.println("round "+i+" "+todo.size());
		System.err.println(todo);
		boolean reverse = i==0;
		Closeable br = reverse ? new ReversedLinesFileReader(node_dmp) : new BufferedReader(new FileReader(node_dmp));
		while((st = reverse ?  ((ReversedLinesFileReader)br).readLine(): ((BufferedReader)br).readLine())!=null){
			String[] str = st.split("\\|");
			Integer child = Integer.parseInt(str[0].trim());
			Integer parent =Integer.parseInt( str[1].trim());
		//	if(parent.intValue()==5439574){
			//	System.err.println(st);
			//}
			if(todo.contains(child) && ! done.contains(child)){
				todo.add(parent);
				set.add(parent);
				todo.remove(child);
				if(!done.contains(parent)) newS.add(parent);
				//todo.remove(child);
				done.add(child);
				pw1.println(st);
			}
			//prev = st;
		}
		br.close();
		notP.addAll(todo);
		todo.clear();
		todo.addAll(newS);
		if(todo.size()==ps) break;
		ps = todo.size();
	}
	//pw1.println(prev);
	
		pw1.close();
		ReversedLinesFileReader br1 = new ReversedLinesFileReader(node_dmp1);
		File node_dmp2 = new File(node_dmp.getAbsolutePath()+"."+tme+".2");
		PrintWriter pw2 = new PrintWriter(new FileWriter(node_dmp2));
		while((st = br1.readLine())!=null){
			pw2.println(st);
		}
		br1.close();
		pw2.close();
		node_dmp1.delete();
		node_dmp2.deleteOnExit();
		return node_dmp2;
	}
 

public  File name_dmp2;
  
  public GetTaxonID( File names_dmp, File node_dmp_, Set<Integer>taxon_set, File name_dmp2, boolean expand)  throws IOException{
		this.taxon_set = taxon_set;
		
	File node_dmp =   expand ? expand(taxon_set, node_dmp_, this.taxon_set) : node_dmp_;
	
	PrintWriter pw_2 = expand ?   new PrintWriter(new GZIPOutputStream(new FileOutputStream(name_dmp2))) : null;
	
		  BufferedReader br = getBR(names_dmp);
		  String st = "";
		  while((st = br.readLine())!=null){
			  String[] str = st.split("\t");
			  Integer taxa = Integer.parseInt(str[0]);
			  if(taxon_set.contains(taxa)){
				  if(pw_2!=null) pw_2.println(st);
				  	String nme = str[2];;
				  	String type = str[6];
				 putTaxa(nme, taxa);
				
				  if(type.startsWith("scientific")){
					  taxa2Sci.put(taxa, nme);
				  }
			  }
			
		  }
		  br.close();
	if(pw_2!=null)  pw_2.close();
	  this.addNodeDmp(node_dmp);
		// TODO Auto-generated constructor stub
	}
  
 

public void print(File out) throws IOException{
	  PrintWriter pw = new PrintWriter(new FileWriter(out));
	  for(Iterator<Integer> it = taxon_set.iterator(); it.hasNext();){
		  pw.println(it.next());
	  }
	  pw.close();
  }


public static BufferedReader getBR(File file)throws IOException{
	 BufferedReader br;
		if(file.getName().endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}
		else{
			br = new BufferedReader(new FileReader(file));
		}
		return br;
}
 public void process(File file)throws IOException {
	  BufferedReader br = getBR(file);
		String st;
		while((st = br.readLine())!=null){
		 String[] str = st.split("\t");
		 Integer taxa = this.processAlias(str, st);
		if(taxa!=null) {
			
			this.taxon_set.add(taxa);
		}
		}
}
 
 /* a few special cases */
 static String[] African_Cassava_Mosaic =  
		 (">AF|>AJ|>AM|>AY|>DQ|>EF|>EU|>FJ|>FM|>FN|>FR|>GQ|>GU|>HE|>HG|>HM|>HQ|>J0|>JF|>JN|>JX|>KC|>KF|>KM|>KP|>KR"
		 + "|>KT|>KU|>KX|>X1|>X6|>Z2|>Z8|>AF|>AJ|>AM|>AY|>DQ|>EF|>EU|>FJ|>FM|>FN|>FR|>GQ|"
		 + ">GU|>HE|>HG|>HM|>HQ|>J0|>JF|>JN|>JX|>KC|>KF|>KM|>KP|>KR|>KT|>KU|>KX|>X1|>X6|>Z2|>Z8|>JQ|>KJ").split("\\|");
 static String[][] spec = new String[][]{
		 ">chr	      :>HLA        :>Kqp".split(":"),
		 "Homo_sapiens:Homo_sapiens:Klebsiella_quasipneumoniae".split(":")
 	};
 	static{
 		for(int i=0; i<spec[0].length; i++){
 			spec[0][i] = spec[0][i].trim();
 		}
 		
 	}
 	Integer processAlias(String[] str, String st){
 	 //the fourth column should be taxa id
 		if(str[2].length()>0) return Integer.parseInt(str[2]);
 		String alias1 = collapse(str[1].split("\\s+"), 1);
 	 
 	/* int compg = alias1.indexOf(", complete genome");
 	 if(compg>=0){
 		 alias1 = alias1.substring(0, compg);
 	 }*/
	 if(  st.indexOf("GRCh38")>=0){
		 alias1= "Homo_sapiens";
		// str[0] = specName;
	 }else{
		 for(int i=0; i<African_Cassava_Mosaic.length; i++){
			 if(st.startsWith(African_Cassava_Mosaic[i])){
				 alias1="African_Cassava_Mosaic";
			 }
		 }
		 for(int i=0; i<spec[0].length; i++){
			 if(st.startsWith(spec[0][i])){
				 alias1=spec[1][i];
			 }
		 }
		
	 }
	 if(alias1.startsWith(">")) alias1 = alias1.substring(1);
	 
	return  this.getTaxa(alias1);
	
 	}
 	
 	 public static String collapse(String[] line, int start) {
 		int end = line.length; String string = " ";
 		StringBuffer sb = new StringBuffer(line[start]);

 		for(int i=start+1; i<end; i++){
 			sb.append(string+line[i]);
 		}
 		return sb.toString();
 	}
 
 
/*Map<String, String> slugToTaxon = new HashMap<String, String>();
//Map<String, String> slugToTaxonShort = new HashMap<String, String>();
  void processGenBank(String[] files) throws IOException{
	  for(int i=0; i<files.length; i++){
		  File file = new File(files[i]);
	  BufferedReader br;
		if(file.getName().endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}
		else{
			br = new BufferedReader(new FileReader(file));
		}
		String st = "";
		while((st = br.readLine())!=null){
			if(st.startsWith("#")) continue;
			String[] str = st.split("\t");
			String tax = str[6];
			String species = str[7];
			String slug = Slug.toSlug(species, "");
			String slugs = Slug.toSlug(species,2, "");
					
			slugToTaxon.put(slug, tax);
			slugToTaxonShort.put(slugs, tax);
			
			//
		}
	  }
  }*/


  
}
