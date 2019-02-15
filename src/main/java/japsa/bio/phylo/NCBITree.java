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
	 public static void main(String[] args){
		   try{
			   NCBITree t = new NCBITree(new File(args[0]), null, null);//, new File(args[1]));
			   
			   t.gid = new GetTaxonID(new File("taxonid"), new File("taxdump/names.dmp"));
				 
				t.addSpeciesIndex(new File(args[1]));
				t.gid.err.close();
			   System.err.println("here");
			  String[][] taxa =  t.getTaxonomy(t.tree[0].getExternalNode(0).getIdentifier().getName());
			  System.err.println(Arrays.asList(taxa[0]));
			  System.err.println(Arrays.asList(taxa[1]));
			t.print(new File(args[0]+".mod"));
			//  NCBITree t1 = new NCBITree("commontree.txt.out1");
			 // t1.print(new File("commontree.txt.out2"));*/
			//   System.err.println(t.tree.getExternalNodeCount());
		   }catch(Exception exc){
			   exc.printStackTrace();
		   }
	   }
	 
	 public static CommonTree read(File f) throws IOException{
		 return new NCBITree(f);
	 }
	 
	
	 
	/* private Node getNode(String in) {
		return this.getNode(this.getSlug1(in));
	}*/

	

	public static Tree[] readTree(File f) throws IOException{
		 NCBITree t = new NCBITree(f);
		 return t.tree;
	}
	 
	 public static NCBITree readTree(File f, File speciesIndex) throws IOException{
		 NCBITree t = new NCBITree(f);
		if(speciesIndex!=null){
			t.addSpeciesIndex(speciesIndex);
		}
		return t;
	}
	 
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
	 
	 public Node getNode(String specName, String taxa) {
		 Node n =  slugToNode.get(slug(specName, false));
		 if(n==null) {
			 n = slugToNode1.get(Integer.parseInt(taxa));
		 }
		 return n;
	 }
	 @Override
	public Node getNode(String specName) {
		 Node n;
		// String specName1 = null;
		/* if(gid!=null){
			 specName1 = this.gid.getSciName(specName);
			 if(specName1!=null){
				specName = specName1;
			 }
		 }*/
		 n =  slugToNode.get(slug(specName, false));
		 return n;
	}
	 
	 public Node getSlug1( String alias1){
		 Node n2 = getNode(alias1);
			 if(n2==null){
				 return slugToNode.get("unclassified");
			 }
			return n2;
	 }
	 
	 
	
	 /*str is a line from species index */
	 private void updateTree(String st, int lineno, double bl, PrintWriter missing){
		 String[] str = st.split("\\s+");
		// String specName = str[0];
		String taxa = gid.processAlias(str,st);
		String alias1 = gid.collapse(str, 1);
		String sciname = taxa==null ? null : gid.taxa2Sci.get(taxa);
		Node n;
		if(sciname==null && taxa==null)   {
			n = this.unclassified;
			Node n1 =  this.createFromSpeciesFile(str,alias1,  n, lineno, bl);
			
			System.err.println(n1.getIdentifier().getAttribute("level"));
			System.err.println(n1.getIdentifier().getAttribute("prefix"));
			System.err.println(n1.getIdentifier());
		}
		else{
			n = getNode(sciname, taxa);
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
   
   private Map<String,Node > slugToNode = new HashMap<String, Node>();
  
  
  

//following better for new trees
 //private static final Pattern plusminus = Pattern.compile("[+-][a-zA-Z\\[\\']");
	private static final Pattern plusminus = Pattern.compile("[\\+|-][a-zA-Z\\[\\']");
	private static final Pattern plusminus1 = Pattern.compile("\\+-[1-9]");
	
	 private int getLevel(String line){
		/*int a = nextLine.indexOf('-')+1;
		int b = nextLine.indexOf('+')+1;
		return Math.max(a, b);*/
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
 Map<Integer, Node> slugToNode1 = new HashMap<Integer, Node>();
 
 public boolean putSlug1( Node n){
	  String name = n.getIdentifier().getName();
	  Integer taxon = (Integer) n.getIdentifier().getAttribute("taxon");
	  boolean contains;
		   String slug = slug(name, false);
		   contains = slugToNode.containsKey(slug);
		    if(!contains) slugToNode.put(slug, n);
	 if(taxon !=null) {
		 slugToNode1.put(taxon, n);
		 if(taxon.intValue()==1739614){
			 System.err.println("h");
		 }
	 }
		return contains;
   }
   
 private Node make(String line_, Node child, Integer taxon){
	Node n =  slugToNode.get(slug(line_, false));

	 if(n==null) {
		 n = new SimpleNode(line_, 0.1);
		 n.getIdentifier().setAttribute("taxon", taxon);
		 putSlug1(n);
	 }
	 if(child!=null){
		 n.addChild(child);
	 }
	 return n;
 }
 
private Node make(String line_, int  level, Node parent){
	   String[] lines = line_.split("\t");
	   String line = lines[0];
	   String name = line;
	   String prefix = "";
	   if(level>=0) {
		   name = line.substring(level, line.length());
		   prefix = line.substring(0,level);
	   }
	   int pm_ind = name.indexOf("+-");
	   if(pm_ind>=0){
		   prefix = prefix+name.substring(0,pm_ind+2);
		   name = name.substring(pm_ind+2);
		  
	   }
	   /*String trimmed = name.trim();
	   int ind1 = name.indexOf(trimmed);
	   if(ind1>=0){
		   			prefix = prefix+name.substring(0, ind1); 
				   name = trimmed;
	   }
	   if(prefix.indexOf("+-")<0) {
		   prefix = prefix+"+-";
	   }*/
	   /*if(gid!=null){
		   //makes sure we use scientific name
			  String name1 = this.gid.getName(name);
			  if(name1!=null) name = name1;
		 }*/
	   Node n = new SimpleNode(name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
	   n.getIdentifier().setAttribute("prefix",prefix);
	  for(int i=1; i<lines.length; i++){
		   String[] v = lines[i].split("=");
		  Object value = v[1];
		   if(v[0].equals("taxon")) value = Integer.parseInt((String) value);
		   n.getIdentifier().setAttribute(v[0],value);
	   }
	   putSlug1(n);
	   if(parent!=null) parent.addChild(n);
	   return n;
   }
	
	
Node unclassified	;

public NCBITree(File file, File taxonid, File taxdump) throws IOException {
	this(file, 
	(taxonid!=null && taxonid.exists()) ?  new GetTaxonID(taxonid, taxdump) : null);
}



public NCBITree(File file, GetTaxonID gid) throws IOException {	
	this.gid = gid;
	BufferedReader br;
	if(file.getName().indexOf("nodes.dmp")>=0 && gid!=null){
		gid.addNodeDmp(file);
		System.err.println(gid.taxon_set.size());
		for(Iterator<String> it = gid.taxon_set.iterator();it.hasNext();){
			String nxt = it.next();
			String parent= gid.nodeToParent.get(nxt);
			Node   n = make(gid.taxa2Sci.get(nxt), null, Integer.parseInt(nxt));
			List<Node>l = new ArrayList<Node>();
			l.add(n);
			inner: while(nxt!=null){
				String nextparent = gid.nodeToParent.get(parent);
				if(parent==nxt) parent = null;
				
				if(parent!=null){
					Node p = make(gid.taxa2Sci.get(parent), n, Integer.parseInt(parent));
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
			//System.err.println(l.get(0).getIdentifier());
			if(!roots.contains(n) ){
				System.err.println("new root "+n.getIdentifier());;
				roots.add(n);
			}
		}
	}
}
	public NCBITree(File file) throws IOException {
//	this(f, null);
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
		System.err.println(this.slugToNode.size());
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
