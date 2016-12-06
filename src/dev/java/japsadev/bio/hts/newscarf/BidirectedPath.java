package japsadev.bio.hts.newscarf;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.util.Logging;

import java.util.List;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Path;
import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedPath extends Path{

	public BidirectedPath(){
		super();
	}
	public BidirectedPath(BidirectedGraph graph, String paths){
		super();
		paths=paths.replace(";", ""); //optimized it!
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
	
	public String getId(){
		//need to make the Id unique for both sense and antisense spelling???
		BidirectedNode curNode = (BidirectedNode) getRoot();
		if(getEdgeCount()<1)
			return curNode.getId()+"+";

		String 	retval=curNode.getId(),
				curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode)?"+":"-";
		retval+=curDir;
		for(Edge e:getEdgePath()){
			curNode=e.getOpposite(curNode);
			retval+=curNode.getId();
			curDir=((BidirectedEdge) e).getDir(curNode)?"-":"+"; //note that curNode is target node
			retval+=curDir;
		}

		return retval.trim();
	}
	
	 
	 public Sequence spelling(){

			BidirectedNode curNode = (BidirectedNode) getRoot();
			if(getEdgeCount()<1)
				return curNode.getAttribute("seq");
			
			SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  this.toString());
			Sequence curSeq = curNode.getAttribute("seq");
			boolean curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode);
			curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);

			seq.append(curSeq.subSequence(0, curSeq.length()-BidirectedGraph.getKmerSize()));
			for(Edge e:getEdgePath()){
				curNode=e.getOpposite(curNode);
				curSeq= curNode.getAttribute("seq");
				curDir=!((BidirectedEdge) e).getDir(curNode);
				curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);

				seq.append(curSeq.subSequence(0, curSeq.length()-(curNode==peekNode()?
						0:BidirectedGraph.getKmerSize())));
				
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
			Logging.error("Cannot join path with disagreed first node " + bridge.getRoot().getId());
			return;
		}
		if(((BidirectedEdge) bridge.getEdgePath().get(0)).getDir((AbstractNode) bridge.getRoot())
			== ((BidirectedEdge) peekEdge()).getDir((AbstractNode) peekNode())){
			Logging.error("Conflict direction from the first node " + bridge.getRoot().getId());
			return;
		}
		
		for(Edge e:bridge.getEdgePath()){
			add(e);
		}
	}
	 
}
