package japsadev.bio.hts.scaffold;

import java.util.ArrayList;

import japsadev.bio.hts.scaffold.Edge;
import japsadev.bio.hts.scaffold.Graph;
import japsadev.bio.hts.scaffold.Node;
import japsadev.bio.hts.scaffold.Path;
import japsadev.bio.hts.scaffold.Vertex;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

public class Path implements Comparable<Path>{
	ArrayList<Node> nodes;
	Graph graph;
	int length, deviation; //how this path differ to long read data (todo: by multiple-alignment??)
	public Path(){
		this.nodes=new ArrayList<Node>();
		graph=new Graph();
		length=0;
		deviation=Integer.MAX_VALUE;
	}
	
	public Path(Graph graph){
		this();
		associate(graph);
	} 
	public Path(Path p){
		this(p.graph);
		for(Node node:p.nodes)
			this.nodes.add(node);
		this.length=p.length;
	}
	/*
	 * @param String: a path as in contigs.paths of SPAdes output
	 * For example: 1+,2-,3+
	 */
	public Path(Graph graph, String paths){
		this(graph);
		paths=paths.replace(";", ""); //optimized it!
		String[] comps = paths.split(",");
		for(int i=0; i<comps.length; i++){
			String 	cur = comps[i];
			boolean curDir = cur.contains("+")?true:false;
			Vertex curComp = graph.getVertex(cur.substring(0,cur.length()-1));
			if(curComp == null){
				System.out.println("Could not find Vertex "+ cur.substring(0,cur.length()-1));
				break;
			}		
			if(!nodes.isEmpty()){
				Node lastNode = nodes.get(nodes.size()-1);
				Edge curEdge=new Edge(lastNode.getVertex(), curComp, lastNode.getDirection(), curDir);
				if(!graph.containsEdge(curEdge)){
					System.out.println(curEdge + " doesn't exist in the graph!");
					break;
				}
					
			}
			addNode(curComp, curDir);
		}
	}
	
	public void associate(Graph g){
		graph=g;
	}
	
	public void addNode(Vertex v, boolean dir){
		if(nodes.isEmpty())
			length=v.getSequence().length();
		else
			length+=v.getSequence().length()-Graph.getKmerSize();
		
		nodes.add(new Node(v,dir));
	}
	public void addNode(Node aNode){
		addNode(aNode.getVertex(), aNode.getDirection());
	}
	
	public ArrayList<Node> getComp(){
		return nodes;
	}
	public void setComp(ArrayList<Node> nodes){
		this.length=0;
		for(Node node:nodes)
			this.addNode(node);
	}
	
	public Path rc(){
		Path retval=new Path(graph);
		for(Node node:nodes){
			retval.nodes.add(0, node.getRC());
		}
		return retval;
	}
	
	public String toString(){
		return "P"+getID();
	}
	public String getID(){
		String retval="";
		for(Node node:nodes){
			retval+=node.toString();
		}
		return retval.trim();
	}
	 public Node removeLast(){
		 Node retval=nodes.remove(nodes.size()-1);
		 length-=retval.getSeq().length()-Graph.getKmerSize();
		 return retval;
	 }
	 
	 public Sequence spelling(){
		 SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  this.toString());

		 for(int i=0;i<nodes.size();i++){
			 Node aNode=nodes.get(i);
			 seq.append(aNode.getSeq().subSequence(0, (i==nodes.size()-1?aNode.getSeq().length():aNode.getSeq().length()-Graph.getKmerSize())));
		 }
		 
		 return seq.toSequence();
	 }
	 /**
	  * 
	  * @return Node: start node
	  */
	 public Node getStart(){
		 return nodes.get(0);
	 }
	/*
	 * @return Node: end node
	 */
	public Node getEnd(){
		return nodes.get(nodes.size()-1);			
	}
	
	public boolean isEmpty(){
		return nodes.isEmpty();
	}
	
	public int getDeviation(){
		return this.deviation;
	}
	public void setDeviation(int deviation){
		this.deviation=deviation;
	}

	@Override
	public int compareTo(Path o) {
		// TODO Auto-generated method stub
		return this.deviation-o.deviation;
	}
	
}
