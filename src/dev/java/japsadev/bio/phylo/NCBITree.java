package japsadev.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

/** written by Lachlan Coin to parse txt files from commontree */
public  class NCBITree {
 //Identifier attributes are css for hex value 
 
  public  Tree[] tree;
   
   static class TreePos{
	   public TreePos(int i, int j) {
		   this.tree_index = i;
		   this.node_index = j;
		// TODO Auto-generated constructor stub
	}
	int tree_index;
	   int node_index;
   }
   
   //this maps the slug to its tree and position in it.
   Map<String,TreePos > slugToPos = new HashMap<String, TreePos>();
   
  
   
  
Set<Node> done = new HashSet<Node>();

 Node depthFirstNext(Node curr){

	 if(!curr.isLeaf()){
		 for(int i=0; i<curr.getChildCount(); i++){
			 if(!done.contains(curr.getChild(i))) return curr.getChild(i);
		 }
	 }
	 if(curr.isRoot()) return null;
	 return depthFirstNext(curr.getParent());
 }
   
   public  void print(File out) throws IOException{
	   PrintWriter pw = new PrintWriter(new FileWriter(out));
	   done.clear();
	   for(int i=0; i<tree.length; i++){
		   
		 Node root= tree[i].getRoot();
		   Set<String> done = new HashSet<String>();

		 Node node = root;
		inner: for(int j=0; node !=null  ;j++){
			 System.err.println(i+" "+j);
			 Identifier id  = node.getIdentifier();
			 String nme =id.getName();
			
			 int level = ((Integer)id.getAttribute("level")).intValue();
			 String hex = ((String)id.getAttribute("css"));		
			 String prefix = ((String)id.getAttribute("prefix"));		
			
			//System.err
			 pw.print(prefix+nme+"\t");
			 if(hex!=null) pw.println(hex);
			 else pw.println();
			 this.done.add(node);
			 node = depthFirstNext(node);
			 if(node==root) break inner;
		 }
	   }
	   pw.close();
   }
   
 
	private static final Pattern plusminus = Pattern.compile("[+|-][a-zA-Z\\[\\']");
	
	 int getLevel(String nextLine){
		/*int a = nextLine.indexOf('-')+1;
		int b = nextLine.indexOf('+')+1;
		return Math.max(a, b);*/
		 Matcher matcher = plusminus.matcher(nextLine);
		 if(matcher.find())
		 	return  matcher.start()+1;
		else return 0;
		
	}
	
	
   public static void main(String[] args){
	   try{
		   NCBITree t = new NCBITree(args[0]);
		  t.print(new File("test"));
		//   System.err.println(t.tree.getExternalNodeCount());
	   }catch(Exception exc){
		   exc.printStackTrace();
	   }
   }
   public static Node make(String line, int  level){
	   String name = line;
	   String prefix = "";
	   if(level>=0) {
		   name = line.substring(level, line.length());
		   prefix = line.substring(0,level);
	   }
	   Node n = new SimpleNode(name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
	   n.getIdentifier().setAttribute("prefix",prefix);
	   return n;
   }
	
	
	
	
public NCBITree(String file) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(file));
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
				this.slugToPos.put(Slug.toSlug(name), new TreePos(i,j));
			}
		}
		
		br.close();
	}
	
	
	
	
	
	public static NCBITree readTree(String f) throws IOException{
		 NCBITree t = new NCBITree(f);
		 
		 return t;
	}
	
}
