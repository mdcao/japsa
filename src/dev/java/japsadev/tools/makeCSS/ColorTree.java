package japsadev.tools.makeCSS;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.Normalizer;
import java.text.Normalizer.Form;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import mdsj.MDSJ;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;
import pal.tree.TreeUtils;
//
//

public class ColorTree {
	boolean slug = true;
	boolean species = true;
	double[][] distances ;
	Identifier[]  identifiers;
	double[] startend ; // this partitions the hue space
static double maxlight = 85;
//https://journals.sagepub.com/doi/full/10.4137/EBO.S7565#_i6
public static void main(String[] args){
	try{
	System.err.println(getHex(240,100,30,1));
		ColorTree ct = new ColorTree(args[1] , false, args[0].equals("species"));
		ct.color();
		String default_hex = "#8888887F";
		
		//ct.addColorToIndex(args[2], Integer.parseInt(args[3]), "\\s+", default_hex, getHex(0,0,100,0));
		ct.printSlug(args[2]);

	}catch(Exception exc){
		exc.printStackTrace();
		
	}
}


public void printSlug(String out) throws Exception{
	PrintWriter pw = new PrintWriter(new FileWriter(out));
	for(Iterator<String> it = this.colors.keySet().iterator(); it.hasNext();){
		String key = it.next();
		String value = colors.get(key);
		pw.println(key+"\t"+value);
		
	}
	pw.close();
}

 static double rand =0;
 
 
 static double[][]getMatrix(Tree tree){
	
	 int cnt = tree.getExternalNodeCount();
		// int cnt1 = cnt;
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
	 double[] startend_ = this.startend;
	 double[][] X = distances;
	 System.err.println(X.length);
	 System.err.println(identifiers.length);
 	 color(X);
 }
 

 
// static double[] lightnessByDepth = new double[] {80,70,60,50,40,30,20,10,5,4,3,2,1}; 
 
public  void color(double[][] X) throws Exception {

	int cnt=X.length;    // number of data objects
	double[][] U=MDSJ.classicalScaling(X); // apply MDS


	//double maxnorm = 0;
	//double minnorm =1e5;
	//double maxtheta = 0;
	//double mintheta = 1e5;
	double[] norms = new double[cnt];
	double[] thetas = new double[cnt];
	SortedSet<Double> norms_set = new TreeSet<Double>();
	SortedSet<Double> thetas_set = new TreeSet<Double>();

	for(int i=0; i<cnt; i++){
		//double r = eig1[i];//U.get(i, 0);
		//double g = eig2[i];//U.get(i, 1);
		double r = U[0][i];
		double g = U[1][i];
		//double b = U.get(i, 2);
	//	norms[i] = Math.sqrt(Math.pow(r, 2) + Math.pow(g, 2)) + Math.random()*1e-5;;
	//	thetas[i] = (Math.atan2(g, r) ) + Math.random()*1e-5;
		norms[i] = r + Math.random()*1e-5;;
		thetas[i] = g+ Math.random()*1e-5;
		norms_set.add(norms[i]);
		thetas_set.add(thetas[i]);
	}
	Map<String, double[]> colors1 = new HashMap<String, double[]>();
	double normsize = norms_set.size();
	double thetasize = thetas_set.size();
	PrintWriter range = new PrintWriter(new FileWriter(new File("range.txt")));
	double maxdepth = 0;
	for(int i=0; i<cnt; i++){
		
		double h = ((double)norms_set.headSet(norms[i]).size())/normsize;
		double t = ((double)thetas_set.headSet(thetas[i]).size())/thetasize;
		 t = t*0.6  + 0.4;  //to make sure >0.4 saturatoin
		// range.println(h+" "+t);
		//String value = getHex(h * 360.0,t * 100.0, 80.0,1.0);
		String key =slug ?  Slug.toSlug(identifiers[i].getName()): identifiers[i].getName();
		if(colors1.containsKey(key)) {
			
			System.err.println("warning already has "+key);
		}
		double depth = ((Number)this.identifiers[i].getAttribute("level")).intValue();
		if(depth>maxdepth)maxdepth = depth;
		colors1.put(key, new double[] {h,t, depth});
		
	}
	System.err.println("max depth " +maxdepth);
	int sze1 = colors1.size();
	System.err.println(sze1);
	for(int i=internal.length-1; i>=0;  i--){
		Node n = internal[i];
		int cc = n.getChildCount();
		double[] d = new double[]{0,0,0};
		d[2] = ((Number)n.getIdentifier().getAttribute("level")).intValue();
		List<String> l = new ArrayList<String>();
		String ident = n.getIdentifier().getName();
		String key1 =slug ?  Slug.toSlug(ident): ident;
	
		for(int j =0; j<cc; j++){
			Node child = n.getChild(j);
			String c_ident = child.getIdentifier().getName();
			String key =slug ?  Slug.toSlug(c_ident): c_ident;
			l.add(key);
			double[] d1 = colors1.get(key);
			d[0]+=d1[0];
			d[1]+=d1[1];
			
		}
		d[0] = d[0]/cc;
		d[1] = d[1]/cc;
		
		
		if(colors1.containsKey(key1)){
			System.err.println("warning1 already has "+key1);
		}
		colors1.put(key1, d);
	}
	int sze2 = colors1.size();
	System.err.println(sze2);
	for(Iterator<String> it = colors1.keySet().iterator(); it.hasNext();){
		String key1 = it.next();
		double[] d = colors1.get(key1);
		double h = d[0];
		double t = d[1]; 
		double avgdepth =  (d[2]);
	
		double lightness = maxlight*(avgdepth/maxdepth);
		range.println(key1+""+h+" "+t+" "+avgdepth+" "+lightness);
		try{
			if(lightness>=0){
				String value = getHex(h * 360,t * 100,lightness,1.0);
				colors.put(key1, value);
		}
		}catch(Exception exc){
			System.err.println("problem with "+key1);
			exc.printStackTrace();
		}
	}
	range.close();
}




Map<String, String> colors = new HashMap<String, String>();
//Map<String, String> groups1 = new HashMap<String, String>();
//Map<String, List<String>> groups= new HashMap<String, List<String>>(); //maps slug to the group it belongs
Map<String, double[]> colors1 = new HashMap<String, double[]>();

 
Node[] internal = new Node[0];
Tree tree; 

ColorTree(String f, boolean split, boolean species) throws Exception{
		colors.put("grch38",	"#ffffff00");  // transparent for human
		slug = true;
		
		 tree =species ? NCBITree.readTree(f) : AntibioticTree.readTree(f);;
		int cnt = tree.getExternalNodeCount();
		System.err.println("read tree with "+cnt);
	//Tree tree = new ReadTree(f);
	setBL(tree.getRoot(), 100, 0.5);
	
	{
		distances =  getMatrix(tree);
		Identifier[] identifier = getIdentifiers(tree);
		//if(!species) this.getGroups(identifier, this.slug, groups);
		//else{
			internal  = NodeUtils.getInternalNodes(tree.getRoot(), true);
		//}
		this.identifiers = (identifier);
		this.startend= (new double[] {0,1});
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
		dist = dist*(1+ Math.random()*rand);
		setBL(root.getChild(i), dist, frac);
	}
	}
	
}

static class Slug{
	private static final Pattern NONLATIN = Pattern.compile("[^\\w-]");
	private static final Pattern WHITESPACE = Pattern.compile("[\\s_]");
	public static String toSlug(String input) {
	  Matcher matcher = WHITESPACE.matcher(input);
	 String nowhitespace = WHITESPACE.matcher(input).replaceAll("_");
	  String normalized = Normalizer.normalize(nowhitespace, Form.NFD);
	  String slug = NONLATIN.matcher(normalized).replaceAll("");
	  return slug.toLowerCase(Locale.ENGLISH);
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
