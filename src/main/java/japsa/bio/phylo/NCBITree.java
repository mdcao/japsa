package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

/** written by Lachlan Coin to parse txt files from commontree */
public  class NCBITree extends CommonTree {
	
	  private static final Logger LOG = LoggerFactory.getLogger(NCBITree.class);
	



	 GetTaxonID gid = null;
	
	 // this modifies the arrays on the count tags, but will set counts to zero if no count information 
	 public void modAll(int i, int len) {
			String[] tags = new String[] {count_tag, count_tag1};
			for(int k=0; k<roots.size(); k++){
				Node root = roots.get(k);
				//if(root.getChildCount()>0){
					Iterator<Node> it = NodeUtils.preOrderIterator(root);
					while(it.hasNext()){
						Identifier id = it.next().getIdentifier();
						for(int j=0; j<tags.length; j++)
						{
							String count_tag2 = tags[j];
							Number[] cnt = (Number[]) id.getAttribute(count_tag2);
							if(cnt==null) cnt=new Number[] {0.0};
							Number[] v = new Number[len];
							Arrays.fill(v, 0);
							v[i] = cnt[0];
							id.setAttribute(count_tag2, v);
						}
					}
				//}
			}
			
		}
	 public void zeroCounts(int i, int len) {
			String[] tags = new String[] {count_tag, count_tag1};
			for(int k=0; k<roots.size(); k++){
				Node root = roots.get(k);
				//if(root.getChildCount()>0){
					Iterator<Node> it = NodeUtils.preOrderIterator(root);
					while(it.hasNext()){
						Identifier id = it.next().getIdentifier();
						for(int j=0; j<tags.length; j++)
						{
							String count_tag2 = tags[j];
							//Integer cnt = (Integer) id.getAttribute(count_tag2);
							//if(cnt==null) cnt=0;
							Number[] v = new Number[len];
							v[i] = 0;
							id.setAttribute(count_tag2, v);
						}
					}
				//}
			}
			
		}
	
	 
	 public Integer getTaxa(String sciname){
		 if(gid==null) return Integer.parseInt(sciname);
		 return gid.getTaxa(sciname);
	 }
	 
	 public Node getNode(String sciName) {
		 if(useTaxaAsSlug){
			return  slugToNode1.get(getTaxa(sciName));
		 }else{
			 
			 String slugn = this.slug(sciName,false);
			 return slugToNode.get(slugn);
		 }
	 }
public Node getNode(Integer taxa) {
		 
		 Node n =  slugToNode1.get(taxa);
		
		 return n;
	 }
	
	 
	 	 
	 
	
	 /*str is a line from species index */
	 private void updateTree(String st, int lineno, double bl, PrintWriter missing, int col_ind, Trie trie){
		 String[] str = st.split("\t");
		// String specName = str[0];
		Integer taxa = Integer.parseInt(str[col_ind]); //gid.processAlias(str,st);
//	  String alias1 = GetTaxonID.collapse(str[1].split("\\s+"), 1);
		//String sciname = taxa==null ? null : gid.taxa2Sci.get(taxa);
		 Node	n = this.slugToNode1.get(taxa);//getNode( taxa);
		 if(n==null){
			 taxa =  trie.find(str[0]);
			// System.err.println("putting at root "+st);
			 n = this.slugToNode1.get(taxa);
		 }
			 this.createFromSpeciesFile(str, n, lineno, bl);
			//System.err.println(n.getIdentifier().getName());
	
		
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
		//	System.err.println("removing1 "+n.getIdentifier().getName());
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
	
	 int getLevel(String line){
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
	
	 private Node createFromSpeciesFile(String[] line, Node parent, int cnt, double bl){
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
		 n.getIdentifier().setName(prefix+n.getIdentifier().getName());
		 id.setAttribute("taxon", Integer.parseInt(line[2]));
		 id.setAttribute("level",level);
		 id.setAttribute("prefix", "  "+prefix);
		// id.setAttribute("alias", line[0]);
	//	 id.setAttribute("alias1", alias1);
		 id.setAttribute("speciesIndex", cnt);
		 String css =(String) pid.getAttribute("css");
		 if(css!=null){
			 id.setAttribute("css", css);
		 }
		 parent.addChild(n);
		 n.setParent(parent);
		// String slug = this.slug(name, false);
		 //System.err.println("new node "+);
		 //this.slugToNode.put(slug, n);
		 this.putSlug1(n);
		 return n;
	 }
  

   
   /** moved all slug operations to GetTaxonID */
   final boolean useslug =true;
	final static int short_slug_length = 3;
	static String slug_sep = "";//"_";
 String slug(String name, boolean shortslug){
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
 
     public Node getSlug(String st){
    	 return slugToNode.get(slug(st, false));
     }
 public void putSlug1( Node n){
	
	// if(useTaxaAsSlug){
	  Integer taxon = (Integer) n.getIdentifier().getAttribute("taxon");
	  boolean   contains = slugToNode1.containsKey(taxon);
		if(!contains) {
			slugToNode1.put(taxon, n);
		}
	//	else{
		//	System.err.println("already contains "+taxon);
	//	}
	// }else{
		 String name = slug(n.getIdentifier().getName(), false);
		    contains = slugToNode.containsKey(name);
			if(!contains) {
				slugToNode.put(name, n);
			}
			else{
				//System.err.println("already contains "+name);
			}
	 //}
   }
 
 
 
 private Node make( Integer taxon, String sci, Node child){
	if(sci==null) sci =  gid.taxa2Sci.get(taxon);
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
		 child.setParent(n);
	 }
	 return n;
 }
//PrintWriter err;

	final boolean kraken;
	
	static boolean trim = false;

	
	public Node make( Node parent, String suffix){
		Identifier pid = parent.getIdentifier();
		   String name = pid.getName();
		   String name1 = name+suffix;
		//   String prefix = "";
		   Node n = new SimpleNode(name1, 0.1);
		   Identifier id =  n.getIdentifier();
		   parent.addChild(n);
		   n.setParent(parent);
			 Integer taxon = ((Integer)pid.getAttribute("taxon"));	
		  if(taxon!=null){
			  if(name2Taxa.containsKey(name1)){
				  throw new RuntimeException("!! "+name);
			  }
			  name2Taxa.put(name1, taxon);
			  id.setAttribute("taxon", taxon);
		  }
		  // int taxon1 = taxon+1;
		  id.setAttribute("level",(Integer)pid.getAttribute("level")+1);
		 id.setAttribute("prefix",(String)pid.getAttribute("prefix")+1);
		  this.slugToNode.put(slug(name1, false), n);
		   //putSlug1(n);
		
			 
		  
		   return n;
	   }
	
	
private Node make(String line_, int  level, Node parent, int index){
	  
	   String[] lines = line_.split("\t");
	   String line = lines[index];
	   String name = line;
	   int j=0;
	 
	 //  err.println(name);err.flush();
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
	   while(name2Taxa.containsKey(name)){
		   name = name+"_"+j;
		   j++;
	   }
	   Node n = new SimpleNode(trim ? name.trim() : prefix+name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
	   n.getIdentifier().setAttribute("prefix",prefix);
	   Integer taxonvalue = null;
	   if(kraken && !line.equals("unclassified")){
		   int taxon =  lines[index-1].equals("null") ? 0 : Integer.parseInt(lines[index-1]); 
		   if(name2Taxa.containsKey(name)){
				  throw new RuntimeException("!! "+name );
			  }
		   this.name2Taxa.put(name, taxon);
		 //  this.taxa2Node.put(taxon, n);
		   n.getIdentifier().setAttribute("taxon",taxon);
		  
		  // Integer[] l1 = new Integer[] {Integer.parseInt(lines[1]), Integer.parseInt(lines[2])};
		   n.getIdentifier().setAttribute(NCBITree.count_tag,new Number[] {Double.parseDouble(lines[1])});
		   n.getIdentifier().setAttribute(NCBITree.count_tag1,new Number[] {Double.parseDouble(lines[2])});

	   }else{
	  for(int i=1; i<lines.length; i++){
		   String[] v = lines[i].split("=");
		  Object value = v[1];
		   if(v[0].equals("taxon")){
			   taxonvalue = Integer.parseInt((String) value);
			   value = taxonvalue;
			   this.name2Taxa.put(name, taxonvalue);
			//   this.taxa2Node.put(taxonvalue, n);
		   }
		   n.getIdentifier().setAttribute(v[0],value);
	   }
	   }
	  if(name.equals("unclassified")){
		  n.getIdentifier().setAttribute("taxon", -1);
		  
	  }
	  
	   putSlug1(n);
	   if(parent!=null){
		   parent.addChild(n);
		   MergeKrakenCmd.checkDupl(parent);

	   }
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
	this.index=0;
	kraken =false;
	this.useTaxaAsSlug = true;
	BufferedReader br;
		
		for(Iterator<Integer> it = gid.taxon_set.iterator();it.hasNext();){
			Integer nxt = it.next();
			Integer parent= gid.nodeToParent.get(nxt);
			if(parent==null) {
				System.err.println("warning did not contain taxa "+nxt);
				continue;
			}
			Node   n = make( nxt,null,  null);
			List<Node>l = new ArrayList<Node>();
			l.add(n);
			inner: while(nxt!=null){
				Integer nextparent = gid.nodeToParent.get(parent);
				if(parent.equals(nxt)) parent = null;
				
				if(parent!=null){
					Node p = make(parent, null, n);
					l.add(p);
					n = p;
				}
				if(nextparent==null ){//|| nextparent.equals("1")){
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
				System.err.println("new root "+n.getIdentifier());;
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

final int index;

public NCBITree(File file, boolean useTaxaAsAsslug, boolean kraken) throws IOException {
this(getBR(file), useTaxaAsAsslug, kraken);
}

static BufferedReader getBR(File file) throws IOException{
	BufferedReader br;
	if(file.getName().endsWith(".gz")){
		br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
	}
	else{
		br = new BufferedReader(new FileReader(file));
	}
return br;
}

public NCBITree(BufferedReader br, boolean useTaxaAsAsslug, boolean kraken) throws IOException {
//	this(f, null);
	index = kraken ? 5 : 0;
	this.kraken = kraken;
	    this.useTaxaAsSlug = useTaxaAsAsslug;
	//	err= new PrintWriter(new FileWriter(new File("error.txt")));
	//	BufferedReader br;
		//inner1: for(int i=0; i<file.length; i++){
		
		String nextLine = br.readLine();
		
		if(nextLine.indexOf("unclassified")>=0 || nextLine.equals("null")){
			nextLine = br.readLine();
		}
		while(nextLine.startsWith(" ")) { // should use regex here
			nextLine = nextLine.substring(1);
		}
		String[] nextLines = nextLine.split("\t");
	/*	if(!nextLine.startsWith("unclassified")){
			unclassified = make("unclassified", 0, null,0);
			unclassified.getIdentifier().setAttribute("css", "#595959aa");
			roots.add(unclassified);
		}*/
	
	//assume unclassified are first
	
	
	outer: while(nextLine !=null){
		Node parent; 
		Node root;
		int parentlevel;
		
		root= make(nextLine, 0, null, index);
		if(root.getIdentifier().getName().equals("unclassified")) unclassified = root;
		roots.add(root);
		parent = root;
		parentlevel =0;
		inner: while((nextLine = br.readLine())!=null){
		/*	if(nextLine.indexOf("Viruses")>=0){
				System.err.println(nextLine);
			}*/
			nextLines = nextLine.split("\t");
			if(nextLine.startsWith("--")){
				nextLine=br.readLine();
				nextLines = nextLine==null ? null : nextLine.split("\t");
				continue outer;
			}
			
			int level = getLevel(nextLines[index]);
			Node n = null;
			//System.err.println(level+"->"+nextLine);
		inner2:	while(level<=parentlevel){
				if(parent==root) {
					n = make(nextLine, 2, root, index);
				//	System.err.println("excluding  " +nextLine);
					
					break inner2;
//					continue inner;
				}
				parent = parent.getParent();
				parentlevel = ((Integer) parent.getIdentifier().getAttribute("level")).intValue();
			}
			if(n==null) n = make(nextLine, level, parent, index);
		
			
		//	System.err.println(n.getIdentifier()+ "--->  "+parent.getIdentifier());
			parentlevel = level;
			parent = n;
		}
		
		}
		
	//   if(err!=null){
	//	   err.close();
	 //  }
		//for(int i=1; i<roots.size(); i++){
		//	Node root = roots.get(i);
		//	if(root.getIdentifier().getName().equals("unclassified")){
		//		throw new RuntimeException("unclassified should be first entry, if it exists");
		//	}
		//}
		if(!kraken)makeTrees(false);
		br.close();
		
	}
		
public static String count_tag = "cumul.css"; // this is cumulative
public  static String count_tag1 = "sep.css"; // this is for sep




public static double thresh=0.0001;
//static String count_below_tag = "count_below";
public void trim(double thresh_perc){
	double total =0;
	for(int i=roots.size()-1;i>=0; i--){
		if(roots.get(i)!=null){
		Number[] cnts = (Number[]) roots.get(i).getIdentifier().getAttribute(NCBITree.count_tag);
		for(int j=0; j<cnts.length; j++) total+=cnts[j].doubleValue();
		}
	}
	double thresh = Math.round(thresh_perc*total);
	//int thresh = Math.max(1, (int) );
	for(int i=roots.size()-1;i>=0; i--){
		if(roots.get(i)!=null){
		if(trimNode(roots.get(i), thresh_perc)){
			roots.remove(i);
		}
		}
	}
}
public void removeDupl(){
	for(int i=roots.size()-1;i>=0; i--){
		if(roots.get(i)!=null)
		removeDupl(roots.get(i));
	}
}
public void removeDupl(Node node){
	String nme = node.getIdentifier().getName().trim();
	while(node.getChildCount()>0 && node.getChild(0).getIdentifier().getName().trim().equals(nme) && 
			node.getChild(0).getChildCount()>0){
		Node child = node.getChild(0);
		node.removeChild(0);
		for(int j=child.getChildCount()-1; j>=0; j--){
			node.addChild(child.removeChild(j));
		}
	}
	for(int j=node.getChildCount()-1;j>=0; j--){
		removeDupl(node.getChild(j));
	}
}
public boolean  trimNode(Node node, double thresh){
	String count_tag2 = NCBITree.count_tag;
	Number[] v = (Number[])node.getIdentifier().getAttribute(count_tag2);
	double sum =0; 
	for(int i=0; i<v.length; i++){
		sum +=v[i].doubleValue();
	}
	if(sum<thresh){
	//	System.err.println("removing "+node.getIdentifier().getName()+" "+sum);
		return true;
	}
	else{
		for(int j=node.getChildCount()-1;j>=0; j--){
			if(trimNode(node.getChild(j), thresh)){
				node.removeChild(j);
			}
		}
	}
	return false;
}

//n is from the new tree. new_parent is from existing tree and will become the new parent

Set<String > added = new HashSet<String>();

public void merge(Node n, int pos){
	System.err.println(n.getIdentifier());
	final Integer taxon = (Integer) n.getIdentifier().getAttribute("taxon");
	Node node = this.slugToNode1.get(taxon);
	if(taxon==1224){
		System.err.println("h");
	}
	if(taxon==687329){
		System.err.println("h");
	}
	if(taxon==10239){
		System.err.println("h");
	}
	/*if(node==null){
		String sl = slug(n.getIdentifier().getName(),false);
		node = this.slugToNode.get(sl);
		if(true) throw new RuntimeException("!!");
	}*/
	//System.err.println(n.getIdentifier());
	if(node!=null){
		checkInTree(node);
		//System.err.println("already has "+n.getIdentifier());
		//this is adding in new samples
		{
			String count_tag2 = NCBITree.count_tag;
			Number[] v = (Number[])node.getIdentifier().getAttribute(count_tag2);
			Number[] v1 = (Number[])n.getIdentifier().getAttribute(count_tag2);
			for(int i=0; i<v.length; i++){
				v[i] = (v[i].doubleValue() +v1[i].doubleValue());
			}
		}
		{
			String count_tag2 = NCBITree.count_tag1;
			Number[] v = (Number[])node.getIdentifier().getAttribute(count_tag2);
			Number[] v1 = (Number[])n.getIdentifier().getAttribute(count_tag2);
			for(int i=0; i<v.length; i++){
				v[i] =v[i].doubleValue()+v1[i].doubleValue();
			}
		}
	}else{
		Node old_parent = n.getParent();
		Integer tax = (Integer) old_parent.getIdentifier().getAttribute("taxon");
		Node new_parent = this.slugToNode1.get(tax);
		if(new_parent==null){
				throw new RuntimeException("!!");
		}
		boolean contains = false;
		for(int i=0; i<new_parent.getChildCount(); i++){
			if(new_parent.getChild(i).equals(n)){
				contains =true;
			}else if(new_parent.getChild(i).getIdentifier().toString().equals(n.getIdentifier().toString())){
				throw new RuntimeException("!!");
			}
		}
		if(!contains) {
			String id = n.getIdentifier().getName().trim();
			if(added.contains(id)){
				throw new RuntimeException("!!");
			}else{
				added.add(id);
			}
			System.err.println("added "+id+"->"+new_parent.getIdentifier().getName());
			if(new_parent.isRoot()){
				if(new_parent!= this.roots.get(0)){
					System.err.println("problem");
				}
				System.err.println("h");
			}
			if(old_parent.isRoot()){
				System.err.println("h");
			}
			checkInTree(new_parent);
			new_parent.addChild(n);
		//	this.putChildren(n);
		}
		checkInTree(n);
	//	n.setParent(new_parent);
		this.slugToNode1.put(taxon,n);
		String sl = this.slug(n.getIdentifier().getName(), false);
		if(slugToNode.containsKey(sl)){
			throw new RuntimeException ("!!!");
		}
		else{
			this.slugToNode.put(sl, n);
		}
	}
	for(int i=0; i<n.getChildCount(); i++){
		merge(n.getChild(i),pos);
	}
	
}

public List<Node> checkInTree(Node node) {
	Node node1 = node;
	List<Node> l = new ArrayList<Node>();
	l.add(node);
	while(!node.isRoot()){
		node = node.getParent();
		l.add(node);
	}
	if(node!=this.roots.get(0)){
		System.err.println("problem for "+node1.getIdentifier());
		throw new RuntimeException("different root!! "+node.getIdentifier()+" "+node.getIdentifier());
	}
	return l;
}
public void merge(NCBITree tree1, int pos){
	this.name2Taxa.putAll(tree1.name2Taxa);
	if(tree1.roots.size()>1) {
		
		throw new RuntimeException(roots.size()+" roots : "+roots.get(0).getIdentifier()+ " vs "+roots.get(1).getIdentifier());
	}
	/*first check to see if the roots are actually sub nodes in this new tree
	for(int i=this.roots.size()-1;  i>=0; i--){
		Node root = this.roots.get(i);
		String id1 = root.getIdentifier().getName();
		Number[] num = (Number[])root.getIdentifier().getAttribute(count_tag);
		Number[] num1 = (Number[]) root.getIdentifier().getAttribute(count_tag1);
		String nme = root.getIdentifier().getName();
		System.err.println("root "+nme);
		if(nme.equals("unclassified")) continue;
		Node mtch = tree1.getSlug(nme);
		if(mtch==null || mtch.isRoot()) continue;
		boolean remove=false;
		if(mtch!=null ){
			inner: while(!mtch.isRoot()){
				Node mtchp = mtch.getParent();
				Integer taxon = (Integer) mtchp.getIdentifier().getAttribute("taxon");
				String sci = mtchp.getIdentifier().getName();
				Node newn = this.getNode(taxon);
				if(newn!=null){
					remove = true;
					newn.addChild(root);
					break inner;
				}else{
				 newn  = this.make(taxon, sci, root);
				}
				Number[] num_ = Arrays.copyOf(num, num.length);
				Number[] num1_ = new Number[num.length];
				Arrays.fill(num1_, 0.0);
				newn.getIdentifier().setAttribute(count_tag, num_);
				newn.getIdentifier().setAttribute(count_tag1, num1_);
				mtch =  mtchp;
				root = newn;
			}
			String id2 = root.getIdentifier().getName();
			System.err.println("switched "+id1+" "+id2);
			if(remove) roots.remove(i);
			else this.roots.set(i, root);
		}
	}*/
	for(int i=0; i<tree1.roots.size(); i++){
		Node root  = tree1.roots.get(i);
		//String nme = root.getIdentifier().getName();
		//if(nme.equals("unclassified")) continue;
		//Node mtch = this.getSlug(nme);
		/*for(int ij=0; ij<roots.size(); ij++){
			Node root_ij = roots.get(ij);
			if(root_ij!=null && root_ij.getIdentifier().getName().equals(nme)){
				mtch = roots.get(ij);
			}
		}*/
		
		/*if(mtch==null){
			String id = root.getIdentifier().getName().trim();
			if(added.contains(id)) {
				throw new RuntimeException("!!" + id);
			}
			System.err.println("added "+id);
			added.add(id);
			Node n = root;
			
			putChildren(n, false);
			mtch = n;
//			this.slugToNode.put(slug(id,false), root);
	//		this.slugToNode1.put((Integer)root.getIdentifier().getAttribute("taxon"),r);
			roots.add(root); // add new root
		}
		else{
		*/
			merge(root, pos);
			MergeKrakenCmd.checkDupl(root);
		//}
		/*	for(int j=0; j<root.getChildCount(); j++){
				merge(root.getChild(j), pos);
			}*/
		}
}
	
	
	/*private void putChildren(Node n, boolean recursive) {
		Integer taxon = (Integer) n.getIdentifier().getAttribute("taxon");
		if(slugToNode1.containsKey(taxon)){
			throw new RuntimeException("!!");
		}
		this.slugToNode1.put(taxon,n);
		
		String sl = this.slug(n.getIdentifier().getName(), false);
		this.slugToNode.put(sl, n);
		
		if(recursive){
		for(int i=0; i<n.getChildCount(); i++){
			putChildren(n.getChild(i), recursive);
		}
		}
}*/
	public NCBITree(File treein, boolean b) throws IOException{
		this(treein, b, false);
}

	public void makeTrees(boolean restrict){
		System.err.println("making trees");
		this.tree = new Tree[roots.size()];
		System.err.println(tree.length);
	//	System.err.println(this.slugToNode.size());
		for(int i=0; i<tree.length; i++){
			Node root = roots.get(i);
			if(root.getParent()!=null) throw new RuntimeException("!!");
			if(restrict){
				while(root.getChildCount()==1){
		
				Node root1 = root.getChild(0);
				String id = root.getIdentifier().getName();
				String id1 = root1.getIdentifier().getName();
				if(root1.getIdentifier().getName().equals(root.getIdentifier().getName())){
					 throw new RuntimeException("!!");
				}
				System.err.println(id+" "+id1);
				root = root1;
			}
			}
			System.err.println("making "+root.getIdentifier()+" "+root.getChildCount());
		//	System.err.println(root.getChild(0).getIdentifier()+", "+ root.getChild(1).getIdentifier());
		//	System.err.println(root);
			//Node following = NodeUtils.postorderSuccessor(root);
			tree[i] = new SimpleTree(root);
			System.err.println("done");
			int cnt1 =tree[i].getExternalNodeCount();
			//System.err.println(cnt+" "+cnt1);
		}
		
	}



	public void addSpeciesIndex(File speciesIndex, int col_ind, Trie trie) throws  IOException{
		if(speciesIndex!=null && speciesIndex.exists()){
			PrintWriter missing = new PrintWriter(new FileWriter("missing.txt"));
				 BufferedReader br1 = GetTaxonID.getBR(speciesIndex);
				 String st = "";
				for(int i=0; (st = br1.readLine())!=null; i++){
					updateTree(st, i, 0.0, missing, col_ind, trie);
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

	public void split() {
		List<String> mtch = Arrays.asList(new String[] {slug("unclassified", false), slug("root",false),  slug("cellular organisms",false),slug("bacteria", false)});
		int len = roots.size();
		List<Node> roots1 = new ArrayList<Node>();
		for(int i=0;i<len; i++){
			
			Node n = roots.get(i);
			if(roots.get(i)==null) continue;
			String nme  = slug(n.getIdentifier().getName().trim(),false);
			System.err.println(nme+" "+mtch.contains(nme));
			if(mtch.contains(nme)){
			if(n.getChildCount()>0){
				for(int j=n.getChildCount()-1;j>=0; j-- ){
					roots1.add(n.removeChild(j));
				}
			}
			}else{
				roots1.add(n);
			}
		}
		this.roots = roots1;
		
	}

	public void removeSingleNodes() {
		int len = roots.size();
		for(int i=len-1; i>=0; i--){
			int cnt = roots.get(i).getChildCount();
			int cnt1 = NodeUtils.getExternalNodes(roots.get(i)).length;
			System.err.println(roots.get(i).getIdentifier()+" "+cnt);
			if(cnt==0 || cnt1==1) { //|| !target.contains(roots.get(i).getIdentifier().getAttribute("taxon"))){
		//		System.err.println("removing "+roots.get(i).getIdentifier());
				roots.remove(i);
			}
		}
		
	}
	public void removeLeafNodes(String count_tag_, boolean add) {
		Node root = roots.get(0);
		Vector<Node> store = new Vector<Node>();
		NodeUtils.getExternalNodes(root, store);
		for(Iterator<Node> it = store.iterator(); it.hasNext();){
			Node n = it.next();
//			n.getIdentifier().getAttribute(count_tag);
			Number[] n1 = (Number[]) n.getIdentifier().getAttribute(count_tag_); // should be the sep tag
			Node p = n.getParent();
			Number[] p1 = (Number[]) n.getIdentifier().getAttribute(count_tag_); 
		//	System.err.println(Arrays.asList(n1));
	//		System.err.println(Arrays.asList(p1));
		//	System.err.println("....");
			if(add){
			for(int i=0; i<p1.length; i++){
				p1[i]=new Double(p1[i].doubleValue()+n1[i].doubleValue());
			}
			}
			
			inner: for(int j=0; j<p.getChildCount(); j++){
				if(p.getChild(j)==n){
					p.removeChild(j);
					break inner;
				}
			}
		}
		
	}
	
	
	
	
	
	
	
	
}
