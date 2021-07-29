package japsa.bio.phylo;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;

import japsa.bio.np.RealtimeSpeciesTyping;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;

public abstract class CommonTree {

	/**Returns String[][] res where res[0] is an array of taxonomy and res[1] is an array of css colors 
	  * includes current node all the way up to the root
	 * @throws IOException 
	  * */
	

	 /* (non-Javadoc)
	 * @see japsadev.bio.phylo.CommonTree#print(java.io.File)
	 */

	public void printKraken(Node node, OutputStreamWriter pw, double total, boolean recursive) throws IOException{
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName().replaceAll("\\+-","");
	//	 System.err.println(nme);
		 Integer level = ((Integer)id.getAttribute("level")).intValue();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		 Integer taxon = ((Integer)id.getAttribute("taxon"));	
		// if(taxon==null || taxon <0) return;
		 double height = node.getNodeHeight();
		
		 Number[] nt = (Number[] ) id.getAttribute(NCBITree.count_tag);
		 Number[] nt1 = (Number[] ) id.getAttribute(NCBITree.count_tag1);
		 String perc = nt==null ? "null": String.format("%5.3g", 100*nt[0].doubleValue()/ (double) total).trim();
		/* if(taxon==null){
			 taxon=0; nme="root";
		 }*/
		 String cntString = nt==null ? "NA\tNA": String.format("%5.3g", nt[0].doubleValue()).trim()+"\t"+String.format("%5.3g", nt1[0].doubleValue()).trim();
		 Double coverage = (Double) id.getAttribute(RealtimeSpeciesTyping.fraction_covered);
		 if(coverage==null) coverage =Double.NaN;
		 pw.write(perc+"\t"+cntString+"\t"+level+"\t"+taxon+"\t"+nme+"\t"+String.format("%5.3g",coverage).trim());
		 pw.write("\n");
		 if(recursive){
			 for(int i=0; i<node.getChildCount(); i++){
				 printKraken(node.getChild(i), pw, total, recursive);
			 }
		 }
	}
	
	int getDistToLeaf(Node node) {
		if(node.isLeaf()) return 0;
		else{
			int val = Integer.MAX_VALUE;
			for(int i=0; i<node.getChildCount(); i++){
				int dist = 1+getDistToLeaf(node.getChild(i));
				if(dist < val) val = dist;
			}
			return val;
		}
	}
	
	public void print(Node node, OutputStreamWriter pw, String[] count_tag, String[] format, boolean recursive, boolean combined) throws IOException{
		if(combined) printCombined(node, pw, count_tag, format, recursive);
		else print1(node, pw, count_tag, format, recursive);
		
	}
	static double bl = 0.04;
	
	private void printCombined(Node node, OutputStreamWriter pw, String[] count_tag, String[] format, boolean recursive) throws IOException{
		 Identifier id  = node.getIdentifier();
		 String nme =id.getName();

		 if(nme.trim().startsWith("+-")){
			 nme = id.getName().replace("+-", "");
			 id.setName(nme);
		 }
		// Integer[] counts = (Integer[]) id.getAttribute(count_tag);
		 Integer level = (int) Math.round((((Integer)id.getAttribute("level")).doubleValue()+1.0)/2.0);
		 int level1 =getDistToLeaf(node);
		 if(node.isRoot()) level =0;
		String height = String.format("%5.3g", NodeUtils.getMinimumPathLengthLengthToLeaf(node)/bl).trim();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		Integer taxon = ((Integer)id.getAttribute("taxon"));
		double[] cssvals = (double[])id.getAttribute("cssvals");
		String[] cssvals1;
		if(cssvals!=null){
		cssvals1 = new String[cssvals.length];
		for(int i=0; i<cssvals.length; i++) cssvals1[i] = String.format("%5.3g",cssvals[i]).trim();
		}else{
			 cssvals1 = "0:0:0".split(":");
		}
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
if(nme.length()>0){
		 pw.write(nme+"\t"+hex+"\t"+taxon+"\t"+height+"\t"+level+"\t"+level1+"\t"+Arrays.asList(cssvals1)+"\t"+sb.toString()+"\t"+sb1.toString());
		 if(count_tag!=null){
			 for(int i=0; i<count_tag.length; i++){
				 Number[] num = (Number[])id.getAttribute(count_tag[i]);
				 double[] d = new double[num.length];
				 for(int k=0; k<d.length; k++){
					 d[k] = num[k].doubleValue();
					 pw.write("\t"+String.format(format[i],  d[k]).replaceAll("\\s+",""));
					 
				 }
				
			 }
		 }
		 //if(true) pw.print("\theight="+String.format("%5.3g", height).trim());

		 pw.write("\n");
}
		 if(recursive){
			 for(int i=0; i<node.getChildCount(); i++){
				 printCombined(node.getChild(i), pw, count_tag,  format, recursive);
			 }
		 }
	}
	
	private void print1(Node node, OutputStreamWriter pw, String[] attributes, String[] format, boolean recursive) throws IOException{
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
		 pw.write(prefix+nme);
		 if(hex!=null) pw.write("\tcss="+hex);
		 if(alias!=null) pw.write("\talias="+alias);
		 if(alias1!=null) pw.write("\talias1="+alias1);
		 if(taxon!=null) pw.write("\ttaxon="+taxon);
		 if(true) pw.write("\theight="+String.format("%5.3g", height).trim());
        if(attributes!=null){
        	for(int i=0; i<attributes.length; i++){
        		//System.err.println(format[i]);
        		Number[] nt = (Number[] ) id.getAttribute(attributes[i]);
        		//System.err.println(Arrays);
        		//do we go here? 
        		if(nt!=null){
        		String fmt =String.format(format[i],  nt[0]);//id.getAttribute(attributes[i])); 
        		pw.write("\t"+fmt.replaceAll("\\s+",""));
        		}
        	}
        }
		 pw.write("\n");
		 if(recursive){
			 for(int i=0; i<node.getChildCount(); i++){
				 print1(node.getChild(i), pw, attributes, format, recursive);
			 }
		 }
	}
	public void print(File out) throws IOException{
		OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(out));
		this.print(osw, null, null, false, false);
		osw.close();
	}
	public void print(OutputStreamWriter out, String[] attributes, String[] format, boolean krkn, boolean combined) throws IOException{
		print(out, "------------------------------------\n", null, attributes, format, krkn, combined);
	}
	
	
	
	 public  void print(OutputStreamWriter pw, String sep, String header, String[] count_tag, String[] format, boolean kraken_style, boolean combined) throws IOException{
		  // PrintStream pw 
/*		   if(out.getName().endsWith(".gz")){
			  pw = new PrintStream(new GZIPOutputStream(new FileOutputStream(out)));
		   }else{
			   pw = new PrintStream((new FileOutputStream(out)));
		   }*/
		//   PrintWriter pw = new PrintWriter(new FileWriter(out));
		 
		   if(!kraken_style && header!=null) pw.write(header);pw.write("\n");
		   double  total =0;
			for(int i=roots.size()-1;i>=0; i--){
				Object cnts =  roots.get(i).getIdentifier().getAttribute(NCBITree.count_tag);
				if(cnts!=null){
					if(cnts instanceof Number[] )	for(int j=0; j<((Number[])cnts).length; j++) total+=((Number[]) cnts)[j].doubleValue();
					else if(cnts instanceof Number )	total+=( ((Number) cnts).doubleValue());
					
				}
			}
		   for(int i=0; i<this.roots.size(); i++){
			   //System.err.println(roots.get(i).getIdentifier().getName());
			// Iterator<Node> n = NodeUtils.preOrderIterator(roots.get(i)); 
	//		 Node root = roots.get(i);
	//		 if(kraken_style) printKraken(root, pw, total);
	//		 else print(root, pw, count_tag, format);
			   
		//	inner: while(n.hasNext()){
				 Node node = roots.get(i);
				 if(kraken_style) printKraken(node, pw, total, true);
				 else print(node, pw, count_tag, format, true, combined);
				
			// }
			// pw.print(sep);
		   }
		   
		  // pw.close();
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