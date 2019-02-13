package japsadev.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

/** written by Lachlan Coin to parse txt files from commontree */
public  class NCBITree {
 //Identifier attributes are css for hex value 
 //private static String default_source="src/test/resources/commontree.txt.css.gz";
 
	 public static void main(String[] args){
		   try{
			   NCBITree t = new NCBITree(new File(args[0]));
			  String[][] taxa =  t.getTaxonomy(t.tree[0].getExternalNode(0).getIdentifier().getName());
			  System.err.println(Arrays.asList(taxa[0]));
			  System.err.println(Arrays.asList(taxa[1]));
			t.print(new File("commontree.txt.out1.gz"));
			//  NCBITree t1 = new NCBITree("commontree.txt.out1");
			 // t1.print(new File("commontree.txt.out2"));*/
			//   System.err.println(t.tree.getExternalNodeCount());
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
	   }
	 
	 public static NCBITree read(File f) throws IOException{
		 return new NCBITree(f);
	 }
	 
	 /**Returns String[][] res where res[0] is an array of taxonomy and res[1] is an array of css colors 
	  * includes current node all the way up to the root
	  * */
	 public String[][] getTaxonomy(String in ){
		String slug =  Slug.toSlug(in);
		TreePos p = this.slugToPos.get(slug);
		Node node = null;
		if(p!=null){
			if(p.external) node = 	this.tree[p.tree_index].getExternalNode(p.node_index);
			else node = 	this.tree[p.tree_index].getInternalNode(p.node_index);
		}
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
	 
	 public static Tree[] readTree(File f) throws IOException{
		 NCBITree t = new NCBITree(f);
		 return t.tree;
	}
	
  public  Tree[] tree;
   
  private  static class TreePos{
	   public TreePos(int i, int j, boolean external) {
		   this.tree_index = i;
		   this.node_index = j;
		   this.external= external;
		// TODO Auto-generated constructor stub
	}
	   boolean external;
	int tree_index;
	   int node_index;
   }
   
   //this maps the slug to its tree and position in it.
   private Map<String,TreePos > slugToPos = new HashMap<String, TreePos>();
   
  
   public  void print(File out) throws IOException{
	   PrintStream pw ;
	   if(out.getName().endsWith(".gz")){
		  pw = new PrintStream(new GZIPOutputStream(new FileOutputStream(out)));
	   }else{
		   pw = new PrintStream((new FileOutputStream(out)));
	   }
	//   PrintWriter pw = new PrintWriter(new FileWriter(out));
	 
	   for(int i=0; i<tree.length; i++){
		 Iterator<Node> n = NodeUtils.depthFirstIterator(tree[i].getRoot());  
		
		inner: for(int j=0; n.hasNext()  ;j++){
		//	 System.err.println(i+" "+j);
			 Node node = n.next();
			 Identifier id  = node.getIdentifier();
			 String nme =id.getName();
			
			 int level = ((Integer)id.getAttribute("level")).intValue();
			 String hex = ((String)id.getAttribute("css"));		
			 String prefix = ((String)id.getAttribute("prefix"));		
			
			//System.err
			 pw.print(prefix+nme);
			 if(hex!=null) pw.println("\t"+hex);
			 else pw.println();
		 }
		 pw.println("------------------------------------");
	   }
	   
	   pw.close();
   }   
  


 
	private static final Pattern plusminus = Pattern.compile("[+|-][a-zA-Z\\[\\']");
	
	 private int getLevel(String nextLine){
		/*int a = nextLine.indexOf('-')+1;
		int b = nextLine.indexOf('+')+1;
		return Math.max(a, b);*/
		 Matcher matcher = plusminus.matcher(nextLine);
		 if(matcher.find())
		 	return  matcher.start()+1;
		else return 0;
		
	}
	
	
	
  
   private static Node make(String line_, int  level){
	   String[] lines = line_.split("\t");
	   String line = lines[0];
	   String name = line;
	   String prefix = "";
	   if(level>=0) {
		   name = line.substring(level, line.length());
		   prefix = line.substring(0,level);
	   }
	   Node n = new SimpleNode(name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
	   n.getIdentifier().setAttribute("prefix",prefix);
	   if(lines.length>1){
		   n.getIdentifier().setAttribute("css",lines[1]);
	   }
	   return n;
   }
	
	
	
	
protected NCBITree(File file) throws IOException {
	BufferedReader br;
	if(file.getName().endsWith(".gz")){
		br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
	}
	else{
		br = new BufferedReader(new FileReader(file));
	}
	String nextLine = br.readLine();
	List<Node> roots = new ArrayList<Node>();
	outer: while(nextLine !=null){
		Node parent; 
		Node root;
		int parentlevel;
		
		root= make(nextLine, 0);
		roots.add(root);
		parent = root;
		parentlevel =0;
		inner: while((nextLine = br.readLine())!=null){
			if(nextLine.startsWith("--")){
				nextLine=br.readLine();
				continue outer;
			}
			int level = getLevel(nextLine);
			//System.err.println(level+"->"+nextLine);
			while(level<=parentlevel){
				if(parent==root) {
					System.err.println("excluding  " +nextLine);
					continue inner;
				}
				parent = parent.getParent();
				parentlevel = ((Integer) parent.getIdentifier().getAttribute("level")).intValue();
			}
			Node n = make(nextLine, level);
			parent.addChild(n);
			//System.err.println(n+ "--->  "+parent);
			parentlevel = level;
			parent = n;
		}
		
		}
		this.tree = new Tree[roots.size()];
		for(int i=0; i<tree.length; i++){
			
			tree[i] = new SimpleTree(roots.get(i));
			
			int cnt =tree[i].getExternalNodeCount();
			for(int j=0; j<cnt; j++){
				String name = tree[i].getExternalNode(j).getIdentifier().getName();
				this.slugToPos.put(Slug.toSlug(name), new TreePos(i,j, true));
			}
			cnt =tree[i].getInternalNodeCount();
			for(int j=0; j<cnt; j++){
				String name = tree[i].getInternalNode(j).getIdentifier().getName();
				this.slugToPos.put(Slug.toSlug(name), new TreePos(i,j, false));
			}
		}
		
		br.close();
	}
	
	
	
	
	
	
	
}
