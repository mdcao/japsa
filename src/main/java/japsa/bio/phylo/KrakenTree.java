package japsa.bio.phylo;

import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;

/** written by Lachlan Coin to parse txt files from commontree */
public  class KrakenTree extends NCBITree {


	
	
	void modAll(int i, int len) {
		for(int k=0; k<roots.size(); k++){
			Node root = roots.get(k);
			//if(root.getChildCount()>0){
				Iterator<Node> it = NodeUtils.preOrderIterator(root);
				while(it.hasNext()){
					Identifier id = it.next().getIdentifier();
					Integer cnt = (Integer) id.getAttribute(count_tag);
					if(cnt==null) cnt=0;
					int[] v = new int[len];
					v[i] = cnt;
					id.setAttribute(count_tag, v);
				}
			//}
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
				 level = ((Integer)n.getIdentifier().getAttribute("level")).intValue();
			 }
			 if(level==target_level) return n.getIdentifier().getAttribute(attr).toString();
			 else return "NA";
	}
	
	
	static double bl = 0.04;
	@Override
	public void print(Node node, PrintStream pw){
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName();
		 int[] counts = (int[]) id.getAttribute(NCBITree.count_tag);
		 if(counts!=null){
			 for(int i=0; i<counts.length; i++){
				 pw.print(counts[i]+"\t");
			 }
		 }
		//if(nme.indexOf("unclassified ssRNA")>=0){
		//	System.err.println("h");
		//}
		 Integer level = (int) Math.round((((Integer)id.getAttribute("level")).doubleValue()+1.0)/2.0);
		 
		String height = String.format("%5.3g", NodeUtils.getMinimumPathLengthLengthToLeaf(node)/bl).trim();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		Integer taxon = ((Integer)id.getAttribute("taxon"));	
		 //double height = node.getNodeHeight();
	//	 if(hex!=null) pw.print("\tcss="+hex);
		// if(alias!=null) pw.print("\talias="+alias);
		// if(alias1!=null) pw.print("\talias1="+alias1);
		 pw.print(nme+"\t"+taxon+"\t"+level+"\t"+hex+"\t"+height);
		 for(int i=1; i<=5; i+=2){
			  pw.print("\t"+get(node, i,"taxon")); 
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
