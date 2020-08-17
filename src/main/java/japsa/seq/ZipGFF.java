package japsa.seq;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;


/* *
 * 
 * author @lachlancoin
 * This class takes a gff file and makes a zip archive with each contig as separate entry
 */
public class ZipGFF {

	BufferedReader br;
	final boolean writeDirectToZip;
	File in;
	String next_st;
	String currChrom;
	int cnt =0;
//	Map<String, File> posm = new HashMap<String, File>();
	boolean haschr = false;
	//OutputStreamWriter pw;
	//final File outdir;
	
	 FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter osw;
	
	String getFile(String chrom1){
		String chrom = chrom1;
		if(haschr && !chrom.startsWith("chr")) chrom = "chr"+chrom1;
		if(!haschr && chrom.startsWith("chr")) chrom = chrom1.substring(3);
		//File f = this.posm.get(chrom);
		return chrom;
	}
	final File outzip;
	//final File output;
	
	final File inDir;
	int len;
	
	
	public ZipGFF(File in, File outzip,  boolean writeDirectToZip) throws IOException{
		String zippath = outzip.getAbsolutePath();
		int ind = zippath.indexOf(".zip");
		if(ind<0) throw new RuntimeException("ouzip should end with .zip");
		this.inDir = new File(zippath.substring(0, ind));
		this.outzip = new File(inDir.getAbsolutePath()+".zip");
		this.writeDirectToZip = writeDirectToZip;
		this.in = in;
	//	if(!in.exists()) throw new RuntimeException("input file needs to be .gz"+in.getAbsolutePath());
		if(inDir.exists() && inDir.listFiles().length>0) throw new RuntimeException("this file should not exist "+inDir);
    	inDir.mkdir();
    	len = inDir.getAbsolutePath().length()+1;
    	dest = new FileOutputStream(outzip);
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new  ZipOutputStream(new BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
		br = new BufferedReader(new InputStreamReader(in.getName().endsWith(".gz") ? new GZIPInputStream(new FileInputStream(in)) : new FileInputStream(in)));
		while((next_st = br.readLine())!=null){
			cnt++;
			if(!next_st.startsWith("#")){
				break;
			}
		};
		String[] st = next_st.split("\t");
		cnt++;
		if(st[0].startsWith("chr")) haschr=true;
		currChrom=st[0];
		
			osw = getWriter(currChrom, writeDirectToZip, true);
			//pw  = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf, true))));
			osw.write(next_st);
			osw.write("\n");
	}
	
	public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip, boolean append)throws IOException{
    	if(writeDirectToZip){
    	ZipEntry headings = new ZipEntry(entry);
	    outS.putNextEntry(headings);
	    return osw;
    	}else{
    		File f = new File(inDir, entry);
    		if(append && ! f.exists()) append=false;
    		return new OutputStreamWriter((new FileOutputStream(f,append)));
    	}
	        
    }
	public void closeWriter(OutputStreamWriter osw) throws IOException{
    	if(osw==this.osw){
    	   osw.flush();
           outS.closeEntry();
    	}else{
    		osw.close();
    	}
    }
	
	Set<String> done = new HashSet<String>();
	public void run() throws IOException{
		String currChromStr =currChrom+"\t";
			inner: while((next_st = br.readLine())!=null){
				cnt++;
				if(next_st.startsWith("#")) continue inner;
				if(!next_st.startsWith(currChromStr)){
					done.add(currChrom);
					closeWriter(osw);
					String[] st = next_st.split("\t");
					currChrom=st[0];
					
					
						if(writeDirectToZip && done.contains(currChrom)){
							outzip.deleteOnExit();
							throw new RuntimeException("not sorted gff "+currChrom+ " "+cnt);
						}
						osw  = getWriter(currChrom, writeDirectToZip, true);
					
					//posm.put(currChrom, outf);
					currChromStr = currChrom+"\t";
					
				}
				osw.write(next_st+"\n");
			}
			br.close();
			if(osw!=null) this.closeWriter(osw);
			this.writeAllFiles(inDir.listFiles());
			this.outS.close();

	}
	
	 public void writeHeader(File f, String newname) throws IOException{
	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			writeHeader(f1[i], f1[i].getAbsolutePath().substring(len));
	    		}
	    	}
	    	else{
	    		OutputStreamWriter osw1 = this.getWriter(newname, true,false);
	    		BufferedReader br = new BufferedReader(new FileReader(f));
	            String str = "";
	            while((str = br.readLine())!=null){
	        	  osw1.write(str);osw1.write("\n");
	           }
	          this.closeWriter(osw1);
	           br.close();
	    	}
	    }
	 
	   
	    private void writeAllFiles(File[] f) {
	    	try{
		    	//
		    	for(int i=0; i<f.length; i++){
		    		if(!f[i].getName().startsWith(".")){
		    			if(f[i]!=null && f[i].exists()){
		    				this.writeHeader(f[i],f[i].getAbsolutePath().substring(len));
		    				f[i].delete();
		    			}
		    		}
		    	}
		    	inDir.delete();
		    	}catch(Exception exc){
		    		System.err.println("problem with "+inDir);
		    		exc.printStackTrace();
		    	}
			
		}
}
