package japsadev.tools.makeCSS;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

/** written by Lachlan Coin to parse txt files from commontree */
public  class NCBITree {
 
 
   Tree tree;
   
   
 
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
		 
		   System.err.println(t.tree.getExternalNodeCount());
	   }catch(Exception exc){
		   exc.printStackTrace();
	   }
   }
   public static Node make(String line, int  level){
	   String name = line;
	   if(level>=0) name = line.substring(level, line.length());
	   Node n = new SimpleNode(name, 0.1);
	   n.getIdentifier().setAttribute("level",level);
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
			System.err.println(level+"->"+nextLine);
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
		Node newroot = roots.get(0);
		if(roots.size()>1){
			newroot =   make("root", -1);
			for(int i=0; i<roots.size(); i++){
				newroot.addChild(roots.get(i));
			}
		}
		this.tree = new SimpleTree(newroot);
		br.close();
	}
	
	
	
	
	
	public static Tree readTree(String f) throws IOException{
		 NCBITree t = new NCBITree(f);
		 return t.tree;
	}
	
}
