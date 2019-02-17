package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

/** written by Lachlan Coin to parse txt files from commontree */
public  class NCBITree extends CommonTree {
	
	  private static final Logger LOG = LoggerFactory.getLogger(NCBITree.class);
	

 //Identifier attributes are css for hex value 
 //private static String default_source="src/test/resources/commontree.txt.css.gz";
 

	 GetTaxonID gid = null;
	 
	
	 
	
	 
	/* private Node getNode(String in) {
		return this.getNode(this.getSlug1(in));
	}*/

	

	
	 
	/* private void putSlug1(String str, TreePos tp){
		String slug = Slug1.toSlug1(str);
		if(slugToPos.containsKey(slug)){
			System.err.println("already had "+str+ "->"+slugToPos.get(slug)+"\n"+tp);
		}
		this.slugToPos.put(slug, tp);
	 }
	 private TreePos getSlug1(String str){
		return  this.slugToPos.get(Slug1.toSlug1(str));
	 }*/
	 
	 public Integer getTaxa(String sciname){
		 return gid.getTaxa(sciname);
	 }
	 
	 public Node getNode(String sciName) {
		 if(useTaxaAsSlug){
			return  slugToNode1.get(getTaxa(sciName));
		 }else{
			 return slugToNode.get(sciName);
		 }
	 }
public Node getNode(Integer taxa) {
		 
		 Node n =  slugToNode1.get(taxa);
		
		 return n;
	 }
	
	 
	 	 
	 
	
	 /*str is a line from species index */
	 private void updateTree(String st, int lineno, double bl, PrintWriter missing){
		 String[] str = st.split("\\s+");
		// String specName = str[0];
		Integer taxa = gid.processAlias(str,st);
		String alias1 = gid.collapse(str, 1);
		String sciname = taxa==null ? null : gid.taxa2Sci.get(taxa);
		Node n;
		if(sciname==null && taxa==null)   {
			n = this.unclassified;
			Node n1 =  this.createFromSpeciesFile(str,alias1,  n, lineno, bl);
		//	System.err.println(n1.getIdentifier().getAttribute("level"));
			//System.err.println(n1.getIdentifier().getAttribute("prefix"));
			//System.err.println(n1.getIdentifier());
		}
		else{
			n = this.slugToNode1.get(taxa);//getNode( taxa);
			
			this.createFromSpeciesFile(str,alias1,  n, lineno, bl);
			//System.err.println(n.getIdentifier().getName());
		}
		
		// this.putSlug1(newnode);
	 }
	 
	

	void thinTree( Node n, int i){
		if(n.isLeaf() ) return;
		for(int i1=n.getChildCount()-1; i1>=0;  i1--){
			Node child = n.getChild(i1);
			thinTree(child,i1);
		}
		if(!n.isRoot() && n.getChildCount()==1 && !n.isLeaf() ){
			Node child  = n.getChild(0);
			System.err.println("removing1 "+n.getIdentifier().getName());
			n.getParent().setChild(i, child);
		}
		
	}
	
	 
 
/*
private Node getNode(TreePos tp) {
		// TODO Auto-generated method stub
	if(tp==null) return null;
		return tp.external ? tree[tp.tree_index].getExternalNode(tp.node_index): tree[tp.tree_index].getInternalNode(tp.node_index);
	}
*/



   
  private  static class TreePos{
	   public TreePos(int i, int j, boolean external) {
		   this.tree_index = i;
		   this.node_index = j;
		   this.external= external;
		// TODO Auto-generated constructor stub
	}
	   public String toString(){
		   return tree_index+","+node_index+","+external;
	   }
	   boolean external;
	int tree_index;
	   int node_index;
   }
   
   //private Map<String,Node > slugToNode = new HashMap<String, Node>();
  
  
  

//following better for new trees
 //private static final Pattern plusminus = Pattern.compile("[+-][a-zA-Z\\[\\']");
	private static final Pattern plusminus = Pattern.compile("[\\+|-][a-zA-Z\\[\\']");
	private static final Pattern plusminus1 = Pattern.compile("\\+-[1-9]");
	
	 private int getLevel(String line){
		/*int a = nextLine.indexOf('-')+1;
		int b = nextLine.indexOf('+')+1;
		return Math.max(a, b);*/
		 if(true) return line.indexOf("+-");
		 int level;
		 Matcher matcher = plusminus.matcher(line);
		 if(matcher.find()){
			 level =  matcher.start()+1;
		 }
		else{
			 Matcher matcher1 = plusminus1.matcher(line);
			 if(matcher.find()){
				 level =  matcher.start()+1;
			 }
			 else 
				 level= 0;
		}
		 if(level==0) {
			   String name = line.substring(level, line.length());
			   int pm_ind = name.indexOf("+-");
			   if(pm_ind>=0){
				   level = level + pm_ind+2;
				   //prefix = prefix+name.substring(0,pm_ind+2);
				  // name = name.substring(pm_ind)+2;
				  
			   }
		   }
		return level;
		   
		
	}
	
	//Map<String, String> alias = new HashMap<String, String>();
	
	 private Node createFromSpeciesFile(String[] line, String alias1, Node parent, int cnt, double bl){
		 int index =0;
		 String name = line[index];
		 String[] names = name.split("\\|");
		 if(names.length>3 ) name = names[3];
	if(name.startsWith(">")){
		name = name.substring(1);
	}
		 Node n = new SimpleNode(name,bl);
		 Identifier id = n.getIdentifier();
		 Identifier pid = parent.getIdentifier();
		 int level ;String prefix ;
		 int cc = 0;
		 if(parent!=null) cc = parent.getChildCount();
		
		 if(cc>0){
			level  = (Integer) parent.getChild(0).getIdentifier().getAttribute("level");
			prefix = (String) parent.getChild(0).getIdentifier().getAttribute("prefix");
		 }else{
			 level = ((Integer) pid.getAttribute("level"))+2;
			 if(parent.isRoot()){
				 prefix = "+-";
			 }else{
				 prefix = "  "+((String)pid.getAttribute("prefix"));
			 }
				 
		 }
		 
		 id.setAttribute("level",level);
		 id.setAttribute("prefix", "  "+prefix);
		// id.setAttribute("alias", line[0]);
		 id.setAttribute("alias1", alias1);
		 id.setAttribute("speciesIndex", cnt);
		 String css =(String) pid.getAttribute("css");
		 if(css!=null){
			 id.setAttribute("css", css);
		 }
		 parent.addChild(n);
		 this.putSlug1(n);
		 return n;
	 }
  

   
   /** moved all slug operations to GetTaxonID */
   final boolean useslug =false;
	final static int short_slug_length = 3;
	static String slug_sep = "";//"_";
   private String slug(String name, boolean shortslug){
	   if(!useslug) return name;
	   if(shortslug){
		   throw new RuntimeException("!!");
		  // return Slug.toSlug(name, this.short_slug_length, slug_sep);
	   }
	   else return Slug.toSlug(name, slug_sep);
   }
   
	 final boolean useTaxaAsSlug;
	 Map<String, Node> slugToNode = new HashMap<String, Node>();
     Map<Integer, Node> slugToNode1 = new HashMap<Integer, Node>();
 
 public void putSlug1( Node n){
	
	 if(useTaxaAsSlug){
	  Integer taxon = (Integer) n.getIdentifier().getAttribute("taxon");
	  boolean   contains = slugToNode1.containsKey(taxon);
		if(!contains) {
			slugToNode1.put(taxon, n);
		}
	 }else{
		 String name = n.getIdentifier().getName();
		 boolean   contains = slugToNode.containsKey(name);
			if(!contains) {
				slugToNode.put(name, n);
			}
	 }
   }
   
 private Node make( Integer taxon, Node child){
	String sci =  gid.taxa2Sci.get(taxon);
	Node n =  getNode(taxon);

	 if(n==null) {
		 n = new SimpleNode(sci, 0.1);
		 n.getIdentifier().setAttribute("taxon", taxon);
		 putSlug1(n);
	 }
	 if(child!=null){
		// String name = child.getIdentifier().getName();
		 //if(name.indexOf("Sclerophthora macrospora virus A")>=0){
			//   System.err.println("h");
		   //}
		 n.addChild(child);
	 }
	 return n;
 }
PrintWriter err;

	
	

private Node make(String line_, int  level, Node parent){
	  
	   String[] lines = line_.split("\t");
	   String line = lines[0];
	   String name = line;
	   err.println(name);err.flush();
	   String prefix = "";
	  // if(line_.indexOf("Sclerophthora macrospora virus A")>=0){
	//	   System.err.println("h");
	 //  }
	   
	   if(level>=0) {
		   name = line.substring(level, line.length());
		   prefix = line.substring(0,level);
	   }
	   int pm_ind = name.indexOf("+-");
	   if(pm_ind>=0){
		   prefix = prefix+name.substring(0,pm_ind+2);
		   name = name.substring(pm_ind+2);
		  
	   }
	   Node n = new SimpleNode(name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
	   n.getIdentifier().setAttribute("prefix",prefix);
	   Integer taxonvalue = null;
	  for(int i=1; i<lines.length; i++){
		   String[] v = lines[i].split("=");
		  Object value = v[1];
		   if(v[0].equals("taxon")){
			   taxonvalue = Integer.parseInt((String) value);
			   value = taxonvalue;
			   this.name2Taxa.put(name, taxonvalue);
		   }
		   n.getIdentifier().setAttribute(v[0],value);
	   }
	  if(name.equals("unclassified")){
		  n.getIdentifier().setAttribute("taxon", -1);
		  
	  }
	  
	   putSlug1(n);
	   if(parent!=null) parent.addChild(n);
	   return n;
   }
	
	
Node unclassified	;
/**
public NCBITree(File file, File taxonid, File taxdump) throws IOException {
	this(file, 
	(taxonid!=null && taxonid.exists()) ?  new GetTaxonID(taxonid, taxdump) : null);
} */



public NCBITree(GetTaxonID gid) throws IOException {	
	this.gid = gid;
	this.useTaxaAsSlug = true;
	BufferedReader br;
		
		for(Iterator<Integer> it = gid.taxon_set.iterator();it.hasNext();){
			Integer nxt = it.next();
			Integer parent= gid.nodeToParent.get(nxt);
			Node   n = make( nxt, null);
			List<Node>l = new ArrayList<Node>();
			l.add(n);
			inner: while(nxt!=null){
				Integer nextparent = gid.nodeToParent.get(parent);
				if(parent.equals(nxt)) parent = null;
				
				if(parent!=null){
					Node p = make(parent, n);
					l.add(p);
					n = p;
				}
				if(nextparent.equals("1")){
					nxt = null;
					parent = null;
				}else{
					nxt = parent;
					parent = nextparent;
				}
			}
			for(int i=0;  i<l.size(); i++){
				Node no = l.get(i);
				Identifier id = no.getIdentifier();
				if(id.getAttribute("level")==null){
					int level = l.size()-1-i;
					String prefix ; 
					if(i<l.size()-1){
						char[] c = new char[level+2];
						Arrays.fill(c, ' ');
						c[level]='+';
						c[level+1] = '-';	
						prefix = new String(c);
					}else{
						prefix = "";
					}
					id.setAttribute("level", level);
					id.setAttribute("prefix", prefix);
				}
			}
			Integer i0= (Integer) l.get(0).getIdentifier().getAttribute("taxon");
			//System.err.println(i0);;
			
			//System.err.println(l.get(0).getIdentifier());
			if(!roots.contains(n) ){
				//System.err.println("new root "+n.getIdentifier());;
				roots.add(n);
			}
			if(i0!=null){
			/*	if(i0.equals(new Integer(191289))){
				//	Node n1 = this.slugToNode.get(l.get(0).getIdentifier().getName());
					Node n2 = this.slugToNode1.get(i0);
					for(int i=0; i<l.size(); i++){
						System.err.println(l.get(i).getIdentifier());
					}
					System.err.println("roots");
					for(int i=0; i<roots.size(); i++){
						System.err.println(roots.get(i).getIdentifier());
					}
					System.err.println('h');
				}*/}
		}
	//}
}
Map<String, Integer> name2Taxa= new HashMap<String, Integer>();

public NCBITree(File file, boolean useTaxaAsAsslug) throws IOException {
//	this(f, null);
	    this.useTaxaAsSlug = useTaxaAsAsslug;
		err= new PrintWriter(new FileWriter(new File("error.txt")));
		BufferedReader br;
		if(file.getName().endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}
		else{
			br = new BufferedReader(new FileReader(file));
		}
		String nextLine = br.readLine();
		String[] nextLines = nextLine.split("\t");
		if(!nextLine.startsWith("unclassified")){
			unclassified = make("unclassified", 0, null);
			unclassified.getIdentifier().setAttribute("css", "#595959aa");
			roots.add(unclassified);
		}
	
	//assume unclassified are first
	
	
	outer: while(nextLine !=null){
		Node parent; 
		Node root;
		int parentlevel;
		
		root= make(nextLine, 0, null);
		if(root.getIdentifier().getName().equals("unclassified")) unclassified = root;
		roots.add(root);
		parent = root;
		parentlevel =0;
		inner: while((nextLine = br.readLine())!=null){
			nextLines = nextLine.split("\t");
			// if(nextLine.indexOf("Sclerophthora macrospora virus A")>=0){
			//	   System.err.println("h");
			 //  }
			if(nextLine.startsWith("--")){
				nextLine=br.readLine();
				nextLines = nextLine==null ? null : nextLine.split("\t");
				continue outer;
			}
			
			int level = getLevel(nextLines[0]);
			//System.err.println(level+"->"+nextLine);
			while(level<=parentlevel){
				if(parent==root) {
					Node n = make(nextLine, 2, unclassified);
					System.err.println("excluding  " +nextLine);
					continue inner;
				}
				parent = parent.getParent();
				parentlevel = ((Integer) parent.getIdentifier().getAttribute("level")).intValue();
			}
			Node n = make(nextLine, level, parent);
		
			
			//System.err.println(n+ "--->  "+parent);
			parentlevel = level;
			parent = n;
		}
		
		}
		br.close();
	   if(err!=null){
		   err.close();
	   }

	
		
	
	for(int i=1; i<roots.size(); i++){
		Node root = roots.get(i);
		if(root.getIdentifier().getName().equals("unclassified")){
			throw new RuntimeException("unclassified should be first entry, if it exists");
		}
	}
		makeTrees();
	}
		
	
	
	
	public void makeTrees(){
		System.err.println("making trees");
		this.tree = new Tree[roots.size()];
		System.err.println(tree.length);
	//	System.err.println(this.slugToNode.size());
		for(int i=0; i<tree.length; i++){
			tree[i] = new SimpleTree(roots.get(i));
			int cnt =tree[i].getInternalNodeCount();
			int cnt1 =tree[i].getExternalNodeCount();
		}
	}



	public void addSpeciesIndex(File speciesIndex) throws  IOException{
		if(speciesIndex!=null && speciesIndex.exists()){
			PrintWriter missing = new PrintWriter(new FileWriter("missing.txt"));
				 BufferedReader br1 = new BufferedReader(new FileReader(speciesIndex));
				 String st = "";
				for(int i=0; (st = br1.readLine())!=null; i++){
					updateTree(st, i, 0.0, missing);
				 }
				br1.close();
				missing.close();
				 System.err.println("done");
				
		}
	}

	@Override
	public Tree[] getTrees() {
		return this.tree;
	}
	
	
	
	
	
	
}
