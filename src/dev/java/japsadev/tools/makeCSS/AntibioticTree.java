package japsadev.tools.makeCSS;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

public class AntibioticTree {
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
		
	 
	public AntibioticTree(String file) throws IOException {
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
	
	public static Tree[] readTree(String f) throws IOException{
		 AntibioticTree t = new AntibioticTree(f);
		 return new Tree[] {t.tree};
	}
}
