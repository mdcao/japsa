package japsadev.bio.hts.scaffold;

import java.util.ArrayList;
import java.util.Arrays;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.util.Logging;

public class Path implements Comparable<Path>{
	ArrayList<Node> nodes;
	//ArrayList<Edge> edges;
	Graph graph;
	int length, deviation; //how this path differ to long read data (todo: by multiple-alignment??)
	int gapLen;
	public Path(){
		this.nodes=new ArrayList<Node>();
		//this.edges=new ArrayList<Edge>();
		graph=new Graph();
		gapLen = length = 0;
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
//		for(Edge edge:p.edges)
//			this.edges.add(edge);
		
		this.length = p.length;
		this.gapLen = p.gapLen;
	}
	/*
	 * @param String: a path as in contigs.paths of SPAdes output
	 * For example: 1+,2-,3+
	 */
	public Path(Graph graph, String paths){
		this(graph);
//		paths=paths.replace(";", ""); //optimized it!
//		String[] comps = paths.split(",");
//		for(int i=0; i<comps.length; i++){
//			String 	cur = comps[i];
//			boolean curDir = cur.contains("+")?true:false;
//			Vertex curComp = graph.getVertex(cur.substring(0,cur.length()-1));
//			if(curComp == null){
//				System.out.println("Could not find Vertex "+ cur.substring(0,cur.length()-1));
//				break;
//			}		
//			if(!nodes.isEmpty()){
//				Node lastNode = nodes.get(nodes.size()-1);
//				Edge curEdge=new Edge(lastNode.getVertex(), curComp, lastNode.getDirection(), curDir);
//				if(!graph.containsEdge(curEdge)){
//					System.out.println(curEdge + " doesn't exist in the graph!");
//					break;
//				}
//					
//			}
//			addNode(curComp, curDir);
//		}
		StringHelper.addPathFromSPAdes(this, paths);
	}
	
	public void associate(Graph g){
		graph=g;
	}
	
	public void addNode(Vertex v, boolean dir){
		length+=v.getSequence().length();
		if(!nodes.isEmpty()){
			Node end = getEnd();
			Edge e = graph.getEdge(end.v, v, end.dir, dir);
			if(e!=null)
				length+=e.getWeight();
			
			if(e.getWeight()>0)
				gapLen+=e.getWeight();
		}
		
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
			 
			 if(i<nodes.size()-1){
				 Node nextNode=nodes.get(i+1);
				 Edge e = graph.getEdge(aNode.v, nextNode.v, aNode.dir, nextNode.dir);
				 if(e==null){
					 System.err.println("Graph doesn't contain edge from " + aNode + " to " + nextNode);
					 break;
				 }
				 else if(e.getWeight()<=0){
					 int lenLeft = aNode.getSeq().length()+e.getWeight();
					 if(lenLeft > 0)
						 seq.append(aNode.getSeq().subSequence(0, lenLeft));
					 else if(ScaffoldGraph.verbose)
						 Logging.info( "Ignore " + aNode.toString() + " and " + nextNode.toString() + " overlap: " + lenLeft 
								 		+ "/" + aNode.getSeq().length() );
				 }
				 
				 else{
					 char[] chars = new char[e.getWeight()];
					 Arrays.fill(chars, 'N');
					 Sequence gap = new Sequence(Alphabet.DNA16(), new String(chars), "gaps");
					 seq.append(aNode.getSeq());
					 seq.append(gap);
				 }
				 					 
			 }else
				 seq.append(aNode.getSeq());

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
