package japsadev.bio.hts.newscarf;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;


import java.util.List;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.Path;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BidirectedPath extends Path{
	int deviation; //how this path differ to long read data (todo: by multiple-alignment??)
	private double coverage=-1; //representative coverage of this path (basically the lowest cov from its nodes)
    private static final Logger LOG = LoggerFactory.getLogger(BidirectedPath.class);
    
    @Override
    public void add(Edge edge) {
    	super.add(edge);
    	
    	double newCoverage = peekNode().getNumber("cov");
    	if(coverage<=0)
    		coverage=newCoverage;
    	else
    		coverage=Math.min(coverage, newCoverage);
    	
    }
    
	public BidirectedPath(){
		super();
	}
	public BidirectedPath(BidirectedPath p){
		super();
		if(p!=null && !p.empty()){
			setRoot(p.getRoot());
			for(Edge e:p.getEdgePath())
				add(e);
		}
		deviation=p.deviation;
	}
	
	//This constructor is only used to load in contigs.path from SPAdes
	//So no recursive path here (path contains all primitive edges)
	public BidirectedPath(BidirectedGraph graph, String paths){
		super();
		paths=paths.replace(";", ","); //optimized it!
		String[] comps = paths.split(",");
		if(comps.length<1)
			return;
		String curID = comps[0], nextID;
		boolean curDir = curID.contains("+")?true:false,
				nextDir;
		BidirectedNode curNode = graph.getNode(curID.substring(0,curID.length()-1)),
						nextNode;
		setRoot(curNode);
		for(int i=1; i<comps.length; i++){
			nextID = comps[i];
			nextDir = nextID.contains("+")?true:false;
			nextNode = graph.getNode(nextID.substring(0,nextID.length()-1));		
			BidirectedEdge curEdge=new BidirectedEdge(curNode, nextNode, curDir, !nextDir);
			add(curEdge);
			curDir=nextDir;
			curNode=nextNode;
		}
	}
	public BidirectedPath getReversedComplemented(){
		BidirectedPath rcPath = new BidirectedPath();
		rcPath.setRoot(this.peekNode());
		List<Edge> edges = this.getEdgePath();
		for(int i = edges.size()-1; i>=0; i--)
			rcPath.add(edges.get(i));
		return rcPath;
	}
	//It is not really ID because Path doesn't need an ID
	public String getId(){
		//need to make the Id unique for both sense and antisense spelling???
		BidirectedNode curNode = (BidirectedNode) getRoot();
		if(getEdgeCount()<1)
			return curNode.getId();

		String 	retval=curNode.getId(),
				curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode)?"+":"-";
		retval+=curDir;
		for(Edge e:getEdgePath()){
			curNode=e.getOpposite(curNode);
			retval+=","+curNode.getId();
			curDir=((BidirectedEdge) e).getDir(curNode)?"-":"+"; //note that curNode is target node
			retval+=curDir;
		}

		return retval.trim();
	}
	
	public String toString(){
		return "(" + getId() + ")";
	}
	 
	public Sequence spelling(){
	
		BidirectedNode curNode = (BidirectedNode) getRoot();
		Sequence curSeq = curNode.getAttribute("seq");
	
		if(getEdgeCount()<1)
			return curSeq;
		
		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  this.toString());
		boolean curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode);
		curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
	
		seq.append(curSeq.subSequence(0, curSeq.length()-BidirectedGraph.getKmerSize()));
		for(Edge edge:getEdgePath()){
			for(Edge e:((BidirectedEdge) edge).getPath().getEdgePath()){			
				curNode=e.getOpposite(curNode);
				curSeq= curNode.getAttribute("seq");
				curDir=!((BidirectedEdge) e).getDir(curNode);
				curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
		
				seq.append(curSeq.subSequence(0, curSeq.length()-(curNode==peekNode()?
						0:BidirectedGraph.getKmerSize())));
			}
		}
	 return seq.toSequence();
	}
	 /*
	  * Add a path to the current path. The path to be added must start with the last node
	  * of the current path.
	  */
	public void join(BidirectedPath bridge) {
		if(bridge==null || bridge.size() <=1)
			return;
		if(bridge.getRoot() != peekNode()){
			LOG.error("Cannot join path with disagreed first node " + bridge.getRoot().getId());
			return;
		}
		if(((BidirectedEdge) bridge.getEdgePath().get(0)).getDir((AbstractNode) bridge.getRoot())
			== ((BidirectedEdge) peekEdge()).getDir((AbstractNode) peekNode())){
			LOG.error("Conflict direction from the first node " + bridge.getRoot().getId());
			return;
		}
		//TODO: need a way to check coverage consistent

			
		for(Edge e:bridge.getEdgePath()){
			add(e);
		}
		
		coverage=Math.min(coverage, bridge.coverage);
	}
	
	public int getDeviation(){
		return this.deviation;
	}
	public void setDeviation(int deviation){
		this.deviation=deviation;
	}

	public double getCoverage(){
		return coverage;
	}
	/**
	 * 
	 * @return average depth of this path
	 */
//	public double averageCov(){
//		int len=0;
//		double res=0;
//		for(Node n:getNodePath()){
//			if(BidirectedGraph.isUnique(n)){
//				Sequence seq = (Sequence) n.getAttribute("seq");
//				double cov = Double.parseDouble(seq.getName().split("_")[5]);
//				len+=(n==getRoot())?seq.length():seq.length()-BidirectedGraph.getKmerSize();
//				res+=seq.length()*cov;
//			}
//		}
//		return res/len;
//	}

	public int length() {
		int retval = 0;
		for(Node n:getNodePath()){
			Sequence seq = (Sequence) n.getAttribute("seq");
			retval+=(n==getRoot())?seq.length():seq.length()-BidirectedGraph.getKmerSize();		
		}
		return retval;
	}
}
