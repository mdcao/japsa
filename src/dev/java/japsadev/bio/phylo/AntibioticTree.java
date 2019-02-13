package japsadev.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

public class AntibioticTree implements CommonTree{
	Tree tree;
	Map<String, Node> nodes = new HashMap<String, Node>();
	
	
	 public  Node make(String line, int  level, Node parent){
		   String name = line;
		   Node n = nodes.get(line);
		   if(n==null){
		 //  if(level>=0) name = line.substring(level, line.length());
			   n = new SimpleNode(name, 0.1);
			   n.getIdentifier().setAttribute("level",level);
			   nodes.put(line,  n);
			   if(parent !=null) parent.addChild(n);
		   }
		   return n;
	   }
		
	 
	public AntibioticTree(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String nextLine = br.readLine();
		Node root = make("root", -1, null);
		
		outer: while((nextLine=br.readLine()) !=null){
			String[] str = nextLine.split("\\s+");
			Node parent = make(str[2], 1, root);
			Node child = make(str[1], 2, parent);
		}
		br.close();
		this.tree = new SimpleTree(root);
	}
	
	public static AntibioticTree readTree(File f) throws IOException{
		 return new AntibioticTree(f);
	}


	@Override
	public String[][] getTaxonomy(String in) {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void print(File out) throws IOException {
		// TODO Auto-generated method stub
		
	}


	@Override
	public Tree[] getTrees() {
		return new Tree[] {tree};
	}
}
