package japsadev.bio.hts.newscarf;

import java.security.AccessControlException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.*;
/**
 * Similar to {@link AdjacencyListNode}
 * 
 */
public class BidirectedNode extends AbstractNode {
	protected static final int INITIAL_EDGE_CAPACITY;
	protected static final double GROWTH_FACTOR = 1.1;

	static {
		String p = "org.graphstream.graph.node.initialEdgeCapacity";
		int initialEdgeCapacity = 32;
		try {
			initialEdgeCapacity = Integer.valueOf(System.getProperty(p, "32"));
		} catch (AccessControlException e) {
		}
		INITIAL_EDGE_CAPACITY = initialEdgeCapacity;
	}
	//edges are bidirected, here are 4 sub-types (name based on direction of the arrow relative to the corresponding node):
	protected static final byte OO_EDGE = 0b00; // src-->--<--dst
	protected static final byte OI_EDGE = 0b01; // src-->-->--dst
	protected static final byte IO_EDGE = 0b10; // src--<--<--dst
	protected static final byte II_EDGE = 0b11; // src--<-->--dst
	
	protected BidirectedEdge[] edges;//fast access to edges knowing direction from src
	protected int oStart, degree;

	protected HashMap<AbstractNode, List<BidirectedEdge>> neighborMap; //fast access to edges knowing dst
	// *** Constructor ***

	protected BidirectedNode(AbstractGraph graph, String id) {
		super(graph, id);
		edges = new BidirectedEdge[INITIAL_EDGE_CAPACITY];
		oStart = degree = 0;
		neighborMap = new HashMap<AbstractNode, List<BidirectedEdge>>(
				4 * INITIAL_EDGE_CAPACITY / 3 + 1);
	}

	// *** Helpers ***

	protected byte edgeType(BidirectedEdge e) {
		//return (byte) (e.getDir1()?0:1 + (e.getDir0()?0:1)<<1); //cool but less efficient
		if(e.getDir0())
			if(e.getDir1())
				return OI_EDGE;
			else 
				return OO_EDGE;
		else
			if(e.getDir1())
				return II_EDGE;
			else 
				return IO_EDGE;	
	}

	@SuppressWarnings("unchecked")
	protected BidirectedEdge locateEdge(Node opposite, byte type) {
		List<BidirectedEdge> l = neighborMap.get(opposite);
		if (l == null)
			return null;

		for (BidirectedEdge e : l) {
			if(type==edgeType(e))
				return e;
		}
		return null;
	}

	protected void removeEdge(int i) {
//		System.out.print("From node " + this.getId() + " remove edge number " + i);
//		System.out.println(" from a list of edges: ");
//		for(int j=0;j<edges.length;j++)
//			if(edges[j]!=null)
//				System.out.println("\t"+edges[j].getId());

		//remove from the hashmap
		AbstractNode opposite = edges[i].getOpposite(this);
		List<BidirectedEdge> l = neighborMap.get(opposite);
		l.remove(edges[i]);
		if (l.isEmpty())
			neighborMap.remove(opposite);
		//remove from the array
		if (i >= oStart) {
			edges[i] = edges[--degree];
			edges[degree] = null;
			return;
		}

		edges[i] = edges[--oStart];
		edges[oStart] = edges[--degree];
		edges[degree] = null;

	}

	// *** Callbacks ***

	@Override
	protected boolean addEdgeCallback(AbstractEdge edge) {
		AbstractNode opposite = edge.getOpposite(this);
		List<BidirectedEdge> l = neighborMap.get(opposite);
		if (l == null) {
			l = new LinkedList<BidirectedEdge>();
			neighborMap.put(opposite, l);
		}
		l.add((BidirectedEdge) edge);
		
		// resize edges if necessary
		if (edges.length == degree) {
			BidirectedEdge[] tmp = new BidirectedEdge[(int) (GROWTH_FACTOR * edges.length) + 1];
			System.arraycopy(edges, 0, tmp, 0, edges.length);
			Arrays.fill(edges, null);
			edges = tmp;
		}

		byte type = edgeType((BidirectedEdge) edge);

		if (type <= OI_EDGE) {
			edges[degree++] = (BidirectedEdge) edge;
			return true;
		}

		edges[degree++] = edges[oStart];
		edges[oStart++] =  (BidirectedEdge) edge;
		return true;
	}

	@Override
	protected void removeEdgeCallback(AbstractEdge edge) {
		// locate the edge first
		byte type = edgeType((BidirectedEdge) edge);
		int i = 0;
		if (type <= OI_EDGE)
			i = oStart;
		while (edges[i] != edge)
			i++;
		if(i < degree){ //only remove iff edge is found
			removeEdge(i);
		}
	}

	@Override
	protected void clearCallback() {
		Arrays.fill(edges, 0, degree, null);
		oStart = degree = 0;
	}

	// *** Access methods ***

	@Override
	public int getDegree() {
		return degree;
	}

	@Override
	public int getInDegree() {
		return oStart;
	}

	@Override
	public int getOutDegree() {
		return degree - oStart;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdge(int i) {
		if (i < 0 || i >= degree)
			throw new IndexOutOfBoundsException("Node \"" + this + "\""
					+ " has no edge " + i);
		return (T) edges[i];
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEnteringEdge(int i) {
		if (i < 0 || i >= getInDegree())
			throw new IndexOutOfBoundsException("Node \"" + this + "\""
					+ " has no entering edge " + i);
		return (T) edges[i];
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getLeavingEdge(int i) {
		if (i < 0 || i >= getOutDegree())
			throw new IndexOutOfBoundsException("Node \"" + this + "\""
					+ " has no edge " + i);
		return (T) edges[oStart + i];
	}

	// FIXME: I must override these stupid functions, let's just return random edge among 4 types!!!
	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdgeBetween(Node node) {
		return (T) locateEdge(node, IO_EDGE);
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdgeFrom(Node node) {
		return (T) locateEdge(node, OO_EDGE);
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdgeToward(Node node) {
		return (T) locateEdge(node, II_EDGE);
	}

	// *** Iterators ***

	protected class EdgeIterator<T extends Edge> implements Iterator<T> {
		protected int iPrev, iNext, iEnd;
		//0:in, 1:out, other(2):all
		protected EdgeIterator(int ori) {
			iPrev = -1;
			iNext = 0;
			iEnd = degree;
			if (ori==0)
				iEnd = oStart;
			else if(ori==1)
				iNext = oStart;
		}

		public boolean hasNext() {
			return iNext < iEnd;
		}

		@SuppressWarnings("unchecked")
		public T next() {
			if (iNext >= iEnd)
				throw new NoSuchElementException();
			iPrev = iNext++;
			return (T) edges[iPrev];
		}

		public void remove() {
			if (iPrev == -1)
				throw new IllegalStateException();
			AbstractEdge e = edges[iPrev];
			// do not call the callback because we already know the index
			//graph.removeEdge(e);
			((BidirectedGraph)graph).removeEdgeDup(e, true, e.getSourceNode() != BidirectedNode.this,
					e.getTargetNode() != BidirectedNode.this);
			removeEdge(iPrev);
			iNext = iPrev;
			iPrev = -1;
			iEnd--;
			
		}
	}

	@Override
	public <T extends Edge> Iterator<T> getEdgeIterator() {
		return new EdgeIterator<T>(2);
	}

	@Override
	public <T extends Edge> Iterator<T> getEnteringEdgeIterator() {
		return new EdgeIterator<T>(0);
	}

	@Override
	public <T extends Edge> Iterator<T> getLeavingEdgeIterator() {
		return new EdgeIterator<T>(1);
	}
}
