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


	
	public static void main(String[] args){
		try {
			FileFilter filter = new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return pathname.getName().endsWith("outreport");
				}
				
			};
			File[] f = (new File(".")).listFiles(filter);
			StringBuffer header = new StringBuffer();
			for(int i=0; i<f.length; i++){
				header.append(f[i].getName()+"\t");
			}
			header.append("name\ttaxon");
			KrakenTree[] kt = new KrakenTree[f.length];
			for(int i=0; i<kt.length; i++){
				kt[i] = new KrakenTree(f[i]);
				kt[i].modAll(i, kt.length);
				kt[i].print(new File(i+".combined.txt"),"",
						header.toString()
						);

			}
			
			NCBITree combined = kt[0];
			for(int j=1; j<kt.length; j++){
				combined.merge(kt[j], j);
			}
			combined.print(new File("combined.txt"),"", header.toString());
		///	combined.makeTrees();
			System.err.println(combined);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
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
		 Integer level = ((Integer)id.getAttribute("level")).intValue();
		 String hex = ((String)id.getAttribute("css"));		
		 String alias = ((String)id.getAttribute("alias"));	
		 String alias1 = ((String)id.getAttribute("alias1"));	
		 String prefix = ((String)id.getAttribute("prefix"));	
		Integer taxon = ((Integer)id.getAttribute("taxon"));	
		 //double height = node.getNodeHeight();
		 pw.print(prefix+nme);
	//	 if(hex!=null) pw.print("\tcss="+hex);
		// if(alias!=null) pw.print("\talias="+alias);
		// if(alias1!=null) pw.print("\talias1="+alias1);
		 if(taxon!=null) pw.print("\t"+taxon);
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
