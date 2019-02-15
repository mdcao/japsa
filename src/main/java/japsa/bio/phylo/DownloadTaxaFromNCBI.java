package japsadev.bio.phylo;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DownloadTaxaFromNCBI {
	 private static final Logger LOG = LoggerFactory.getLogger(DownloadTaxaFromNCBI.class);

public static void main(String[] args){
	try{
		runProcess(true);
}catch(Exception exc){
	exc.printStackTrace();
	
}
}


//"wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt";

public static void runProcess(boolean bacterial) throws Exception{

	String str = bacterial ? "bacteria" : "viral";
	File f = new File("assembly_summary.txt");
	String taxonid_file = "taxonid.txt";
	if(!f.exists()){
		String url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+str+"/assembly_summary.txt";
		{
			String cmd = "wget --tries=inf  --waitretry=60 --retry-connrefused -o download.log -O "+url;
			LOG.info("Running " + cmd);
			Process process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			LOG.info("Done " + cmd);
		}
		{
			
			String cmd = "awk -F \"\t\" '$12==\"Complete Genome\" && $11==\"latest\"{print $6}' assembly_summary.txt  > "+taxonid_file;
			LOG.info("Running " + cmd);
			Process process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			LOG.info("Done " + cmd);
		}
		

//upload taxonid_file to https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
//then download commontree.txt
//##mv ~/Downloads/commontree.txt  .
	}
}
}
