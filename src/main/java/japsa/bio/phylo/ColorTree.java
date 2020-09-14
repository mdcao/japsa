package japsa.bio.phylo;

import java.awt.Color;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPOutputStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import mdsj.MDSJ;
import pal.misc.Identifier;
import pal.tree.ClockTree;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;
import pal.tree.TreeUtils;
//
//

public class ColorTree {
	 private static final Logger LOG = LoggerFactory.getLogger(ColorTree.class);
	//boolean slug = true;
	boolean species = true;
	double[][] distances ;
	Node[]  identifiers;
	
	
	//double[] startend ; // this partitions the hue space
static double maxlight = 85;




 static double rand =0;
 static Node[]  getInternalNodes(Node[] nodes){
	 Node anc = NodeUtils.getFirstCommonAncestor(nodes);
	 Set<Node> n = new HashSet<Node>();
	 for(int i=0; i<nodes.length; i++){
		 findAncs(anc, nodes[i], n);
	 }
	 return n.toArray(new Node[0]);
 }
 static double[][] getDist(Node[] nodes){
	 int cnt  = nodes.length;
	 double[][] X = new double[cnt][cnt];
	 for(int i=1; i<cnt; i++){
		 for(int j=0; j<i; j++){
			 double d = getDist(nodes[i], nodes[j]);
			 X[i][j] = d;
			 X[j][i] = d;
		 }
	 }
	 return X;
 }
 private static void findAncs(Node anc, Node n1, Set<Node>s){
		 while(n1!=anc){
			 n1 = n1.getParent();
			 s.add(n1);
		 }
	 }
 private static double getDistToAnc(Node anc, Node n1){
	double  cnt=0;
	 while(n1!=anc){
		 cnt+=n1.getBranchLength();
		 n1 = n1.getParent();
	 }
	 return cnt;
 }
 
 private static double getDist(Node nodeOne, Node nodeTwo) {
	Node anc = NodeUtils.getFirstCommonAncestor(nodeOne, nodeTwo);
	return getDistToAnc(anc, nodeOne) + getDistToAnc(anc, nodeTwo);
}

static double[][]getMatrix(Tree tree){
	 int cnt = tree.getExternalNodeCount();
	 double[] 	idist = new double[tree.getInternalNodeCount()];
	double epsilon = 0;
			boolean countEdges = false;
			
			double[][] X  = new double[cnt][cnt];
			for(int i=0; i<cnt; i++){
			
				Arrays.fill(idist, 0);
				double[] dist  = new double[cnt];
				int minind =0;
				double minval = 1e9;
			    TreeUtils.computeAllDistances(tree, i, dist, idist, countEdges, epsilon);
				for(int j=0; j<cnt; j++) {
					if(j!=i && dist[j]<minval){
						minind = j;
						minval = dist[j];
					}
					X[i][j] = dist[j];
					//X.set(i, j, dist[j]);
				}
				//String nme = 	tree.getExternalNode(i).getIdentifier().getName();
				//String nme1 = 	tree.getExternalNode(minind).getIdentifier().getName();
				//System.err.println(nme+" -> "+nme1);
			}
			return X;
 }
 void color() throws Exception{
	 if(tree!=null){
		 int cnt = tree.getExternalNodeCount();
			System.err.println("read tree with "+cnt + tree.getRoot().getIdentifier().getName());
		this.identifiers = new Node[tree.getExternalNodeCount()];
		for(int i=0; i<identifiers.length; i++){
			identifiers[i] = tree.getExternalNode(i);
		}
		 if(cnt>1){
				distances =  getMatrix(tree);
				internal  = NodeUtils.getInternalNodes(tree.getRoot(), true);
			
			
			}
		// double[] startend_ = this.startend;
		 double[][] X = distances;
		 System.err.println(X.length);
		 System.err.println(identifiers.length);
		 String outf = "distances.txt";
	 	 color(X, identifiers, outf, null);
	 	 this.colorInternal();
	 }
 }
 
 

 
// static double[] lightnessByDepth = new double[] {80,70,60,50,40,30,20,10,5,4,3,2,1}; 
 
public static void print(double[][] X, Node[] identifiers, String outf) throws IOException{
	PrintWriter distpw= new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf))));
	distpw.print("\t");
	 for(int i=0; i<X.length; i++){
		   distpw.print(identifiers[i].getIdentifier().getName()+"\t");
	 }
	 distpw.println();
   for(int i=0; i<X.length; i++){
	   distpw.print(identifiers[i].getIdentifier().getName()+"\t");
	   for(int j=0; j<X.length; j++){
		   distpw.print(String.format("%5.3g", X[i][j]).trim()+"\t");
	   }
	   distpw.println();
   }
	distpw.close();
}
//

public static void color(Node[] identifiers, String outf, double lightness) throws IOException{
	double[][] X = getDist(identifiers);
	color(X, identifiers, outf, lightness);
}

Map<Integer, List<Node>> levelToNode = new TreeMap<Integer, List<Node>>();

public void setLevels(){
	levelToNode.clear();
	setLevel(this.tree.getRoot(), 0);
}
public void setLevel(Node node, Integer level){
	List<Node> l = levelToNode.get(level);
	if(l==null){
		levelToNode.put(level, l = new ArrayList<Node>());
	}
	l.add(node);
	node.getIdentifier().setAttribute("level", level);
	for(int i=0; i<node.getChildCount(); i++){
		setLevel(node.getChild(i), level+1);
	}
}

public void colorEachLevel() throws Exception{
	setLevels();
	for(Iterator<Integer> it = levelToNode.keySet().iterator(); it.hasNext(); ){
		int nxti = it.next();
		Node[]  nxt = levelToNode.get(nxti).toArray(new Node[0]);

		if(nxt.length==1){
			setHSL(nxt[0],1,1,0);
		}else{
			ColorTree.color(nxt, null, 56);
		}
	}
}
public void colorRecurisvely(boolean even){
	Node root = tree.getRoot();
	setHSL(root,1,1,0);
	int max_depth = NodeUtils.getMaxNodeDepth(tree.getRoot());
	
	Node[] child = new Node[root.getChildCount()];
	
	for(int i=0; i<child.length; i++){
		child[i]  = root.getChild(i);
		child[i].getIdentifier().setAttribute("leaves", NodeUtils.getExternalNodes(child[i]).length);
	}
	Arrays.sort(child, new Comparator<Node>(){

		@Override
		public int compare(Node o1, Node o2) {
			return -1*Integer.compare((Integer)o1.getIdentifier().getAttribute("leaves"),(Integer) o2.getIdentifier().getAttribute("leaves"));
		}
		
	});
	
	double[] light_range = new double[] {0.0,1.0};
	double light_min = 0.0;
	double light_max ;
	double step = (light_range[1] - light_range[0])/ (double) child.length;
	List<Integer> order = new ArrayList<Integer>();
	double[] light = new double[child.length];
	for(int i=0; i<child.length; i++){
		light_max = light_min + step;
		
		light[i]  = (light_min + light_max)/2.0;
		light_min= light_max;
	}
	int mid1 =(int) Math.floor((double) child.length/2.0);
	int mid2 = mid1+1;
	//int off = child.length % 2-1;
	for(int i=0;i<child.length ;i++){
		child[i].getIdentifier().setAttribute("range", new double[] {0.0,0.95});
		if(i%2==0 && mid1>=0 || mid2==child.length){
			color(child[i],0, max_depth, light[mid1], even);
			mid1--;
		}else{
			color(child[i],0, max_depth, light[mid2], even);
			mid2++;
		}
	}
}

public static void color(Node node, int depth,int max_depth, double  lightness,  boolean even){
	
	double[] ranget = (double[]) node.getIdentifier().getAttribute("range");
//	double mid = (min + max)/2.0;
	if(depth==0){
		setHSL(node,1,1,0);
	}else{
		setHSL(node, (ranget[0]+ranget[1])/2.0, 0.5  + ((double)depth/(double) max_depth) * 0.5, lightness);
	}
	int cc = node.getChildCount();
	double min = ranget[0];
	double step = (ranget[1] - ranget[0])/(double)cc;
	//double parlen = NodeUtils.getExternalNodes(node).length;
	int sum=0;
	int[] leaves = new int[cc];
	for(int i=0; i<cc; i++){
		leaves[i] = NodeUtils.getLeafCount(node.getChild(i));
		sum+=leaves[i];
	}
	
	for(int i=0; i<cc; i++){
		Node child = node.getChild(i);
		if(!even){
		 step = ((double)leaves[i])/(double) sum * (ranget[1] - ranget[0]);
		}
		double max = min+step;
		child.getIdentifier().setAttribute("range", new double[] {min,max});
		color(node.getChild(i), depth+1,max_depth,  lightness, even);
		min = max;
	}
}

public  static void color(double[][] X, Node[] identifiers, String outf, Double light) throws IOException {
	if(outf!=null){
		print(X, identifiers, outf);
	}
	Node root = NodeUtils.getFirstCommonAncestor(identifiers);
	double maxheight = root.getNodeHeight();
	int cnt=X.length;    // number of data objects
	double[][] U=MDSJ.classicalScaling(X); // apply MDS
	double[] norms = new double[cnt];
	double[] thetas = new double[cnt];
	SortedSet<Double> norms_set = new TreeSet<Double>();
	SortedSet<Double> thetas_set = new TreeSet<Double>();

	for(int i=0; i<cnt; i++){
		//double r = eig1[i];//U.get(i, 0);
		//double g = eig2[i];//U.get(i, 1);
		double x = U[0][i];
		double y = U[1][i];
		double norm = Math.sqrt(Math.pow(x, 2)+ Math.pow(y, 2));
		double theta = - 1.0* Math.atan2(y, x)/Math.PI;
		System.err.println(identifiers[i].getIdentifier().getName()+":"+x+","+y+"->"+norm+","+theta);
		//double b = U.get(i, 2);
	//	norms[i] = Math.sqrt(Math.pow(r, 2) + Math.pow(g, 2)) + Math.random()*1e-5;;
	//	thetas[i] = (Math.atan2(g, r) ) + Math.random()*1e-5;
		norms[i] = norm +  Math.random()*1e-5;;
		thetas[i] =theta +  Math.random()*1e-5;
		norms_set.add(norms[i]);
		thetas_set.add(thetas[i]);
	}
	System.err.println(thetas_set);
	Map<String, double[]> colors1 = new HashMap<String, double[]>();
	double normsize = norms_set.size();
	double thetasize = thetas_set.size();
	//PrintWriter range = new PrintWriter(new FileWriter(new File("range.txt")));
	System.err.println("max depth "+maxheight);
	double theta_min = 0;
	double theta_max = 0.9;
	for(int i=0; i<cnt; i++){
		
		double norm = ((double)norms_set.headSet(norms[i]).size())/normsize;
		double t =((double)thetas_set.headSet(thetas[i]).size())/(thetasize-1);
		t = t*(theta_max-theta_min) + theta_min;
		// norm = norm*0.6  + 0.4;  //to make sure >0.4 saturation
		 double height =identifiers[i].getNodeHeight();
		 identifiers[i].getIdentifier().setAttribute("height", height);

		 double lightness = light==null ? 85*(1- height/maxheight) : light;
		 norm = 0.98;
		 setHSL(identifiers[i], t, norm, lightness);
		 
	}
	System.err.println("max depth " +maxheight);
	int sze1 = colors1.size();
	System.err.println(sze1);
	
}



private static void setHSL(Node node, double h, double s, double lightness) {
	String hexvalue = getHex(h * 360,s * 100,lightness*100,1.0);
	 node.getIdentifier().setAttribute("css", hexvalue);
	 
	 node.getIdentifier().setAttribute("cssvals", new double[] {h,s, lightness});
	
}
public void colorInternal(){
	double maxheight = tree.getRoot().getNodeHeight();
	//NodeUtils.postorderSuccessor(tree.getExternalNode(0));
	for(int i=internal.length-1; i>=0;  i--){
		Node n = internal[i];
		int cc = n.getChildCount();
		double[] d = new double[]{0,0,0};
		List<String> l = new ArrayList<String>();
		Identifier id = n.getIdentifier();
		for(int j =0; j<cc; j++){
			Node child = n.getChild(j);
			Identifier c_ident = child.getIdentifier();
			double[] d1 = (double[]) c_ident.getAttribute("cssvals");
			//System.err.println(d1);
			d[0]+=d1[0];
			d[1]+=d1[1];
		}
		d[0] = d[0]/cc;
		d[1] = d[1]/cc;
		//d[2] = ((Number)n.getIdentifier().getAttribute("level")).intValue();
		 double height =n.getNodeHeight();
		 d[2] = height;
		 double lightness = 85*(1- height/maxheight);
		 setHSL(identifiers[i], d[0] * 360, d[1]*100, lightness);
		//String hexvalue = getHex(, * 100,lightness,1.0);
		//id.setAttribute("css", hexvalue);
		id.setAttribute("height", height);

			id.setAttribute("cssvals", d);
	}
	//range.close();
}




//Map<String, String> colors = new HashMap<String, String>();
//Map<String, String> groups1 = new HashMap<String, String>();
//Map<String, List<String>> groups= new HashMap<String, List<String>>(); //maps slug to the group it belongs
//Map<String, double[]> colors1 = new HashMap<String, double[]>();

 
Node[] internal = new Node[0];
Tree tree; 



ColorTree(Tree tree_in, boolean clock) throws Exception{
//	setBL(tree_in.getRoot(), 1000, 0.1);
	tree_in.getRoot().getIdentifier().setAttribute("css" ,	"#000000ff");   // for homo sapiens
	if(tree_in.getExternalNodeCount()==1){
				tree_in.getExternalNode(0).getIdentifier().setAttribute("css" ,	"#00000000");   // for homo sapiens
	}else{
		this.tree = clock ?  new ClockTree(tree_in) : tree_in ;
		
	//this.startend= (new double[] {0,1});
	
	}
	
}

private static Identifier[] getIdentifiers(Tree stree) {
	Identifier[] ids = new Identifier[stree.getExternalNodeCount()];
	for(int i=0; i<ids.length; i++){
	ids[i] = stree.getIdentifier(i);
	}
	return ids;
}


private static void setBL(Node root, double d, double frac) {
	root.setBranchLength(d);
	if(!root.isLeaf()){
	int cnt = root.getChildCount();
	for(int i=0; i<cnt; i++){
		double dist = d*frac;
	//	dist = dist*(1+ Math.random()*rand);
		setBL(root.getChild(i), dist, frac);
	}
	}
	
}

public static String getHex(double h, double s, double l, double alpha){
	//System.err.println(String.format("%02x", 15));
	HSLColor color = new HSLColor((float) h, (float) s, (float) l, (float) alpha);
	Color rgb = color.getRGB();
	int alph = Math.round(color.getAlpha() * (float)255);
	String hex = String.format("#%02x%02x%02x%02x", rgb.getRed(), rgb.getGreen(), rgb.getBlue(), alph);  
	return hex;
  // return "#"+ Integer.parseInt(r, 16)+ Integer.toHexString(g)+ Integer.toHexString(b);
	
}

}
