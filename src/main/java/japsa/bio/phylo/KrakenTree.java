package japsa.bio.phylo;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;

/** written by Lachlan Coin to parse txt files from commontree */
public  class KrakenTree extends NCBITree {

	public KrakenTree(File file, String str) throws IOException {
		super(find(file,str), true, true);
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
	public KrakenTree(File[] file) throws IOException {
		super(file, true, true);
		// TODO Auto-generated constructor stub
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
	
	
	static double bl = 0.04;
	@Override
	public void print(Node node, PrintStream pw, String[] count_tag, String[] format){
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName();
		// Integer[] counts = (Integer[]) id.getAttribute(count_tag);
		 Integer level = (int) Math.round((((Integer)id.getAttribute("level")).doubleValue()+1.0)/2.0);
		 if(node.isRoot()) level =0;
		String height = String.format("%5.3g", NodeUtils.getMinimumPathLengthLengthToLeaf(node)/bl).trim();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		Integer taxon = ((Integer)id.getAttribute("taxon"));
		double[] cssvals = (double[])id.getAttribute("cssvals");
		String[] cssvals1 = new String[cssvals.length];
		for(int i=0; i<cssvals.length; i++) cssvals1[i] = String.format("%5.3g",cssvals[i]).trim();
		StringBuffer sb = new StringBuffer();
		StringBuffer sb1 = new StringBuffer();
		Stack<String> l = new Stack<String>();
		Stack<String> l1 = new Stack<String>();
		if(!node.isRoot()){
			Node parent = node.getParent();
			while(!parent.isRoot()){
				l.push(parent.getIdentifier().getName());
				l1.push(parent.getIdentifier().getAttribute("taxon")+"");
				//l.push(p)
				//sb.append("->"+parent.getIdentifier().getName());
				parent = parent.getParent();
			}
			while(l.size()>0){
				sb.append(l.pop());
				sb1.append(l1.pop());
				if(l.size()>0){
					sb.append("->");
					sb1.append(",");
				}
			}
		}
		
		 //double height = node.getNodeHeight();
	//	 if(hex!=null) pw.print("\tcss="+hex);
		// if(alias!=null) pw.print("\talias="+alias);
		// if(alias1!=null) pw.print("\talias1="+alias1);
		//root    css=#000000ff   taxon=1 height=1.24
		//			header.append("name\tcolor\ttaxon\theight\tparents\ttaxon1\ttaxon2\ttaxon3\ttaxon4\ttaxon5");
		//header.append("name\tcolor\ttaxon\theight\tlevel\tcssvals\tparents\ttaxon1\ttaxon2\ttaxon3\ttaxon4\ttaxon5");

		 pw.print(nme+"\t"+hex+"\t"+taxon+"\t"+height+"\t"+level+"\t"+Arrays.asList(cssvals1)+"\t"+sb.toString()+"\t"+sb1.toString());
		 if(count_tag!=null){
			 for(int i=0; i<count_tag.length; i++){
				 pw.print(String.format(format[i],  id.getAttribute(count_tag[i])).replaceAll("\\s+",""));
			 }
		 }
		 //if(true) pw.print("\theight="+String.format("%5.3g", height).trim());

		 pw.println();
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
