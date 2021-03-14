package japsa.bio.phylo;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import japsa.bio.np.RealtimeSpeciesTyping;
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

	public void printKraken(Node node, PrintStream pw, int total){
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName().replaceAll("\\+-","");
	//	 System.err.println(nme);
		 Integer level = ((Integer)id.getAttribute("level")).intValue();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		 Integer taxon = ((Integer)id.getAttribute("taxon"));	
		 double height = node.getNodeHeight();
		 Integer[] nt = (Integer[] ) id.getAttribute(NCBITree.count_tag);
		 Integer[] nt1 = (Integer[] ) id.getAttribute(NCBITree.count_tag1);
		 String perc = String.format("%5.3g", 100*(double)nt[0]/ (double) total).trim();
		/* if(taxon==null){
			 taxon=0; nme="root";
		 }*/
		 Double coverage = (Double) id.getAttribute(RealtimeSpeciesTyping.fraction_covered);
		 if(coverage==null) coverage =Double.NaN;
		 pw.print(perc+"\t"+nt[0]+"\t"+nt1[0]+"\t"+level+"\t"+taxon+"\t"+nme+"\t"+String.format("%5.3g",coverage).trim());
		 pw.println();
	}
	public void print(Node node, PrintStream pw, String[] attributes, String[] format){
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
        if(attributes!=null){
        	for(int i=0; i<attributes.length; i++){
        		//System.err.println(format[i]);
        		Integer[] nt = (Integer[] ) id.getAttribute(attributes[i]);
        		//System.err.println(Arrays);
        		if(nt!=null){
        		String fmt =String.format(format[i],  nt[0]);//id.getAttribute(attributes[i])); 
        		pw.print("\t"+fmt.replaceAll("\\s+",""));
        		}
        	}
        }
		 pw.println();
	}
	public void print(File out) throws IOException{
		this.print(out, null, null, false);
	}
	public void print(File out, String[] attributes, String[] format, boolean krkn) throws IOException{
		print(out, "------------------------------------\n", null, attributes, format, krkn);
	}
	
	
	
	 public  void print(File out, String sep, String header, String[] count_tag, String[] format, boolean kraken_style) throws IOException{
		   PrintStream pw ;
		   
		   if(out.getName().endsWith(".gz")){
			  pw = new PrintStream(new GZIPOutputStream(new FileOutputStream(out)));
		   }else{
			   pw = new PrintStream((new FileOutputStream(out)));
		   }
		//   PrintWriter pw = new PrintWriter(new FileWriter(out));
		 
		   if(!kraken_style) pw.println(header);
		   int total =0;
			for(int i=roots.size()-1;i>=0; i--){
				Integer[] cnts = (Integer[]) roots.get(i).getIdentifier().getAttribute(NCBITree.count_tag);
				if(cnts!=null)		for(int j=0; j<cnts.length; j++) total+=cnts[j];
			}
		   for(int i=0; i<this.roots.size(); i++){
			   System.err.println(roots.get(i).getIdentifier().getName());
			 Iterator<Node> n = NodeUtils.preOrderIterator(roots.get(i));  
			
			inner: while(n.hasNext()){
				 Node node = n.next();
				 if(kraken_style) printKraken(node, pw, total);
				 else print(node, pw, count_tag, format);
				
			 }
			// pw.print(sep);
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