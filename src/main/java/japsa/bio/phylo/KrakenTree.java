package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;

/** written by Lachlan Coin to parse txt files from commontree */
public  class KrakenTree extends NCBITree {

	public KrakenTree(File file, String str) throws IOException {
		super(new File(file, str), true, true);
		// TODO Auto-generated constructor stub
	}
	
	public KrakenTree(BufferedReader file, String str) throws IOException {
		super(file, true, true);
		// TODO Auto-generated constructor stub
	}
	
	
	

	static FileFilter dirFilter = new FileFilter(){

		@Override
		public boolean accept(File dir) {
			return dir.isDirectory();
		}
		
	};
	private static File[] find(File file, String str) {
		if(str==null) return new File[] {file};
		List<File> l = new ArrayList<File>();
		find(file, str, l);//, file.getAbsolutePath(), file.getName());
		return l.toArray(new File[0]);
	}
	private static void  find(File file, String str, List<File>l){ // String prefix, String new_prefix) {
		File f1 = new File(file, str);
		if(f1.exists()){
			l.add(f1);//.getAbsolutePath().replace(prefix, new_prefix));
			File[] subdirs = file.listFiles(dirFilter);
			for(int i=0; i<subdirs.length; i++){
				find(subdirs[i], str, l);//, prefix, new_prefix);
			}
		}
	}
	
	
	
	
	public KrakenTree(File file) throws IOException {
		super(file, true, true);
		// TODO Auto-generated constructor stub
	}
	
	public String get(Node n, int target_level, String attr){
		 Integer level = ((Integer)n.getIdentifier().getAttribute("level")).intValue();
			 while(level>target_level){
				 n  = n.getParent();
				 if(n==null) return "NA";
				 level = ((Integer)n.getIdentifier().getAttribute("level")).intValue();
			 }
			 if(level==target_level) return n.getIdentifier().getAttribute(attr).toString();
			 else return "NA";
	}
	
	

	




	static Pattern p = Pattern.compile("[a-zA-Z]");
	 int getLevel(String line){
		 for(int i=0; i<line.length(); i++){
			 if(line.charAt(i)!=' ') {
				 return i-1;
			 }
		 }
		 return -1;
		 /*
		boolean mat =  Pattern.matches("cell", line);
		 Matcher m = p.matcher(line);
		 if(m.matches()){
			 m.find();
			 int st = m.start();
			 return st-1;
		 }else{
			 return -1;
		 }*/
	 }
	
	
	
	
}
