package japsadev.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
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
/** extract taxon ids matching speciesIndex from a list of assembly summary files */
public class GetTaxonID {
  Set<String> taxon_set = new HashSet<String>();
	Map<String, String> name2Taxa = new HashMap<String, String>();
	Map<String, String> name2Taxa3 = new HashMap<String, String>();
	Map<String, String> name2Taxa2 = new HashMap<String, String>();


	Map<String, String> taxa2Sci = new HashMap<String, String>();
  public GetTaxonID(){
  }
  
  public String getName(String specName){
	  String slug1 =  Slug.toSlug(specName, "");
		String slug3 =  Slug.toSlug(specName, 4,"");
		String slug2 =  Slug.toSlug(specName, 3,"");
	 String taxa = this.name2Taxa.get(slug1);
	 if(taxa==null) taxa = name2Taxa3.get(slug3);
	 if(taxa==null) taxa = name2Taxa2.get(slug2);
	 if(taxa!=null){
		 specName = taxa2Sci.get(taxa);
	 }
	 return null;
  }
  void putTaxa(String nme, String taxa){
	  String slug = Slug.toSlug(nme, "");
	  String slug3 = Slug.toSlug(nme, 4,"");
	  String slug2 = Slug.toSlug(nme, 3,"");
	  name2Taxa.put(slug, taxa);
	   name2Taxa3.put(slug3, null);
	 // else  name2Taxa3.put(slug3, taxa);
	  //if(name2Taxa2.containsKey(slug2)) name2Taxa2.put(slug2, null);
	    name2Taxa2.put(slug2, taxa);
  }
  
  public GetTaxonID(File file, File names_dmp)  throws IOException{
	  if(file.exists()){
		  BufferedReader br = getBR(file);
		  String st = "";
		  while((st = br.readLine())!=null){
			  taxon_set.add(st.split("\\s+")[0]);
		  }
		  br.close();
	  }
	  if(names_dmp.exists()){
		  BufferedReader br = getBR(names_dmp);
		  String st = "";
		  while((st = br.readLine())!=null){
			  String[] str = st.split("\t");
			  String taxa = str[0];
			  //if(taxon_set.contains(taxa)){
				  	String nme = str[2];;
				  	String type = str[6];
				 putTaxa(nme, taxa);
				  if(type.startsWith("scientific")){
					  taxa2Sci.put(taxa, nme);
				  }
			  //}
			
		  }
		  br.close();
	  }
		// TODO Auto-generated constructor stub
	}
  
  void print(File out) throws IOException{
	  PrintWriter pw = new PrintWriter(new FileWriter(out));
	  for(Iterator<String> it = taxon_set.iterator(); it.hasNext();){
		  pw.println(it.next());
	  }
	  pw.close();
  }
public static void main(String[] args){
	  try{
		  String[] assembly_summary_inputs = "assembly_summary_genbank.txt.gz:assembly_summary_bacteria.txt.gz:assembly_summary_virus.txt.gz".split(":");
		  //these can be downloaded using DownloadTaxaFromNCBI
		  GetTaxonID gid = new GetTaxonID(new File("taxonid"), new File("taxdump/names.dmp"));
		  
		  if(true){
		 // gid.processGenBank(assembly_summary_inputs);
		  gid.process(new File("speciesIndex"), "add_taxons.txt", "missing_taxons.txt");
		  gid.print(new File("taxonid.new"));
		  }
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
  private void process(File file, String added_taxon, String missing_file)throws IOException {
	  BufferedReader br = getBR(file);
	
		PrintWriter pw = new PrintWriter(new FileWriter(added_taxon));
		PrintWriter missing = new PrintWriter(new FileWriter(missing_file));
		String st;
		while((st = br.readLine())!=null){
		 String[] str = st.split("\\s+");
		 String specName = str[0];
		 if(str[0].indexOf("GRCh38")>=0){
			 specName= "Homo_sapiens";
			 str[0] = specName;
		 }
		 String slug = Slug.toSlug(specName, "");
		 String taxa =  this.name2Taxa.get(slug);
		if(taxa!=null) this.taxon_set.add(taxa);
		else{
			  slug = Slug.toSlug(specName, 2,"");
			  taxa =  this.name2Taxa.get(slug);
			  if(taxa!=null) this.taxon_set.add(taxa);
		}
		
		}
		missing.close();
		pw.close();
	
}
Map<String, String> slugToTaxon = new HashMap<String, String>();
Map<String, String> slugToTaxonShort = new HashMap<String, String>();
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
			/*if(slug.startsWith("borrelia") && slug.indexOf("afz")>0){
				System.err.println(slug);
				System.err.println(slugs);
			}*/
			//
		}
	  }
  }
  
}
