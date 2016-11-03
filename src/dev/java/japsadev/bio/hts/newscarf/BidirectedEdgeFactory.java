package japsadev.bio.hts.newscarf;

import org.graphstream.graph.EdgeFactory;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedEdgeFactory implements EdgeFactory<BidirectedEdge> {

	@Override
	public BidirectedEdge newInstance(String id, Node src, Node dst, boolean directed) {
		// TODO Auto-generated method stub
		return null;
	}
	public BidirectedEdge newInstance(String id, AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1) {
		return new BidirectedEdge(src, dst, dir0, dir1);
	}
}
