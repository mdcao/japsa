package japsa.bio.phylo;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;

public abstract class CommonTree {

	/**Returns String[][] res where res[0] is an array of taxonomy and res[1] is an array of css colors 
	  * includes current node all the way up to the root
	  * */
	

	 /* (non-Javadoc)
	 * @see japsadev.bio.phylo.CommonTree#print(java.io.File)
	 */

	
	public void print(Node node, PrintStream pw, String count_tag){
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName();
		//if(nme.indexOf("unclassified ssRNA")>=0){
		//	System.err.println("h");
		//}
		 Integer level = ((Integer)id.getAttribute("level")).intValue();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		Integer taxon = ((Integer)id.getAttribute("taxon"));	
		 double height = node.getNodeHeight();
		 pw.print(prefix+nme);
		 if(hex!=null) pw.print("\tcss="+hex);
		 if(alias!=null) pw.print("\talias="+alias);
		 if(alias1!=null) pw.print("\talias1="+alias1);
		 if(taxon!=null) pw.print("\ttaxon="+taxon);
		 if(true) pw.print("\theight="+String.format("%5.3g", height).trim());

		 pw.println();
	}
	
	public void print(File out) throws IOException{
		print(out, "------------------------------------\n", null, NCBITree.count_tag);
	}
	
	
	
	 public  void print(File out, String sep, String header, String count_tag) throws IOException{
		   PrintStream pw ;
		   if(out.getName().endsWith(".gz")){
			  pw = new PrintStream(new GZIPOutputStream(new FileOutputStream(out)));
		   }else{
			   pw = new PrintStream((new FileOutputStream(out)));
		   }
		//   PrintWriter pw = new PrintWriter(new FileWriter(out));
		 
		   pw.println(header);
		   for(int i=0; i<this.roots.size(); i++){
			   System.err.println(roots.get(i).getIdentifier().getName());
			 Iterator<Node> n = NodeUtils.preOrderIterator(roots.get(i));  
			
			inner: while(n.hasNext()){
				 Node node = n.next();
				 print(node, pw, count_tag);
				
			 }
			 pw.print(sep);
		   }
		   
		   pw.close();
	   }   
	
	abstract public Node getNode(String specName);
	List<Node> roots = new ArrayList<Node>();
	
	public Tree[] getTrees(){
		return tree;
	}
	public  Tree[] tree;
	
	 /* (non-Javadoc)
		 * @see japsadev.bio.phylo.CommonTree#getTaxonomy(java.lang.String)
		 */
		 
		public String[][] getTaxonomy(String in ){
		
			Node node = this.getNode(in);
			
			 List<String> tax = new ArrayList<String>();
			 List<String> css = new ArrayList<String>();
			 while(node!=null){
				 tax.add(node.getIdentifier().getName());
				 css.add((String) node.getIdentifier().getAttribute("css"));
				 if(node.isRoot()) node = null;
				 else node = node.getParent();
			 }
			 return new String[][] {tax.toArray(new String[0]), css.toArray(new String[0])};
		 }

}