package japsadev.bio.hts.newscarf;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;


import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class BidirectedGraph extends AbstractGraph{
    private int kmer;
    static final int TOLERATE=500;
    
	public static final double GROW_FACTOR = 1.1;
	public static final int DEFAULT_NODE_CAPACITY = 128;
	public static final int DEFAULT_EDGE_CAPACITY = 1024;

	protected HashMap<String, AbstractNode> nodeMap;
	protected HashMap<String, AbstractEdge> edgeMap;

	protected AbstractNode[] nodeArray;
	protected AbstractEdge[] edgeArray;

	protected int nodeCount;
	protected int edgeCount;

	private Method setIndex;
	// *** Constructors ***
	/**
	 * Creates an empty graph.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 * @param strictChecking
	 *            If true any non-fatal error throws an exception.
	 * @param autoCreate
	 *            If true (and strict checking is false), nodes are
	 *            automatically created when referenced when creating a edge,
	 *            even if not yet inserted in the graph.
	 * @param initialNodeCapacity
	 *            Initial capacity of the node storage data structures. Use this
	 *            if you know the approximate maximum number of nodes of the
	 *            graph. The graph can grow beyond this limit, but storage
	 *            reallocation is expensive operation.
	 * @param initialEdgeCapacity
	 *            Initial capacity of the edge storage data structures. Use this
	 *            if you know the approximate maximum number of edges of the
	 *            graph. The graph can grow beyond this limit, but storage
	 *            reallocation is expensive operation.
	 */
	public BidirectedGraph(String id, boolean strictChecking, boolean autoCreate,
			int initialNodeCapacity, int initialEdgeCapacity) {
		super(id, strictChecking, autoCreate);
		// All we need to do is to change the node & edge factory
		setNodeFactory(new NodeFactory<BidirectedNode>() {
			public BidirectedNode newInstance(String id, Graph graph) {
				return new BidirectedNode((AbstractGraph) graph, id);
			}
		});

		setEdgeFactory(new EdgeFactory<BidirectedEdge>() {
			public BidirectedEdge newInstance(String id, Node src, Node dst, boolean directed) { //stupid??
				return new BidirectedEdge(id, (AbstractNode)src, (AbstractNode)dst);
			}
		});
		
	
		//setEdgeFactory(new BidirectedEdgeFactory());
		
		if (initialNodeCapacity < DEFAULT_NODE_CAPACITY)
			initialNodeCapacity = DEFAULT_NODE_CAPACITY;
		if (initialEdgeCapacity < DEFAULT_EDGE_CAPACITY)
			initialEdgeCapacity = DEFAULT_EDGE_CAPACITY;

		nodeMap = new HashMap<String, AbstractNode>(
				4 * initialNodeCapacity / 3 + 1);
		edgeMap = new HashMap<String, AbstractEdge>(
				4 * initialEdgeCapacity / 3 + 1);
		nodeArray = new AbstractNode[initialNodeCapacity];
		edgeArray = new AbstractEdge[initialEdgeCapacity];
		nodeCount = edgeCount = 0;
		
		try {
			setIndex = AbstractElement.class.getDeclaredMethod("setIndex", int.class);
			setIndex.setAccessible(true);
		} catch (Exception  e) {
			e.printStackTrace();
		} 
	}

	/**
	 * Creates an empty graph with default edge and node capacity.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 * @param strictChecking
	 *            If true any non-fatal error throws an exception.
	 * @param autoCreate
	 *            If true (and strict checking is false), nodes are
	 *            automatically created when referenced when creating a edge,
	 *            even if not yet inserted in the graph.
	 */
	public BidirectedGraph(String id, boolean strictChecking, boolean autoCreate) {
		this(id, strictChecking, autoCreate, DEFAULT_NODE_CAPACITY,
				DEFAULT_EDGE_CAPACITY);
	}

	/**
	 * Creates an empty graph with strict checking and without auto-creation.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 */
	public BidirectedGraph(String id) {
		this(id, true, false);
	}
	
	
	/**********************************************************************************
	 * ****************************Algorithms go from here*****************************
	 */
    public BidirectedGraph(){
    	this("Assembly graph",true,false,1000,100000);
        setKmerSize(127);//default kmer size used by SPAdes to assembly MiSeq data
    }
    public void loadFromFile(String graphFile){
        setAutoCreate(true);
        setStrict(false);
		//1. next iterate over again to read the connections
		SequenceReader reader;
		try {
			reader = new FastaReader(graphFile);
			Sequence seq;
			int shortestLen = 10000;
			while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
				if(seq.length()<shortestLen)
					shortestLen=seq.length();
				
				String[] adjList = seq.getName().split(":");
				String name = adjList[0];
				boolean dir0=name.contains("'")?false:true;
				
				name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim(); //EDGE_X_length_Y_cov_Z
				
				//FIXME: constructor is invoked by name but hashmap is based on label!!!
				String nodeID = name.split("_")[1];
				AbstractNode node = addNode(nodeID);
				node.setAttribute("name", name);
				
				if(dir0){
					seq.setName(name);
					//current.setSequence(seq);
					node.setAttribute("seq", seq);
				}
				if (adjList.length > 1){
					String[] nbList = adjList[1].split(",");
					for(int i=0; i < nbList.length; i++){
						String neighbor = nbList[i];
						// note that the direction is read reversely in the dest node
						boolean dir1=neighbor.contains("'")?true:false;
						neighbor=neighbor.replaceAll("[^a-zA-Z0-9_.]", "").trim();
						
						String neighborID = neighbor.split("_")[1];
						AbstractNode nbr = addNode(neighborID);
						
						BidirectedEdge e = new BidirectedEdge(node, nbr, dir0, dir1);
						//e.addAttribute("ui.label", e.getId());
						
						addEdge(e.getId(), node, nbr);
					}
				}
				
			}

			//rough estimation of kmer used
			if((shortestLen-1) != getKmerSize())
				setKmerSize(shortestLen-1);
			
			reader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
    }
	
    public int getKmerSize(){
    	return this.kmer;
    }
    public void setKmerSize(int kmer){
    	this.kmer=kmer;
    }
    
    /*
     * *********************************************************************************************************************
     */

	// *** Callbacks ***

	@Override
	protected void addEdgeCallback(AbstractEdge edge) {
		edgeMap.put(edge.getId(), edge);
		if (edgeCount == edgeArray.length) {
			AbstractEdge[] tmp = new AbstractEdge[(int) (edgeArray.length * GROW_FACTOR) + 1];
			System.arraycopy(edgeArray, 0, tmp, 0, edgeArray.length);
			Arrays.fill(edgeArray, null);
			edgeArray = tmp;
		}
		edgeArray[edgeCount] = edge;
		try {
			setIndex.invoke(edge, edgeCount++);
		} catch (Exception  e) {
			e.printStackTrace();
		} 

	}

	@Override
	protected void addNodeCallback(AbstractNode node) {
		nodeMap.put(node.getId(), node);
		if (nodeCount == nodeArray.length) {
			AbstractNode[] tmp = new AbstractNode[(int) (nodeArray.length * GROW_FACTOR) + 1];
			System.arraycopy(nodeArray, 0, tmp, 0, nodeArray.length);
			Arrays.fill(nodeArray, null);
			nodeArray = tmp;
		}
		nodeArray[nodeCount] = node;
		try {
			setIndex.invoke(node, nodeCount++);
		} catch (Exception  e) {
			e.printStackTrace();
		} 
	}

	@Override
	protected void removeEdgeCallback(AbstractEdge edge) {
		edgeMap.remove(edge.getId());
		int i = edge.getIndex();
		edgeArray[i] = edgeArray[--edgeCount];
		try {
			setIndex.invoke(edgeArray[i], i);
		} catch (Exception  e) {
			e.printStackTrace();
		} 
		
		edgeArray[edgeCount] = null;
	}

	@Override
	protected void removeNodeCallback(AbstractNode node) {
		nodeMap.remove(node.getId());
		int i = node.getIndex();
		nodeArray[i] = nodeArray[--nodeCount];
		try {
			setIndex.invoke(nodeArray[i], i);
		} catch (Exception  e) {
			e.printStackTrace();
		} 
		nodeArray[nodeCount] = null;
	}

	@Override
	protected void clearCallback() {
		nodeMap.clear();
		edgeMap.clear();
		Arrays.fill(nodeArray, 0, nodeCount, null);
		Arrays.fill(edgeArray, 0, edgeCount, null);
		nodeCount = edgeCount = 0;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdge(String id) {
		return (T) edgeMap.get(id);
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Edge> T getEdge(int index) {
		if (index < 0 || index >= edgeCount)
			throw new IndexOutOfBoundsException("Edge " + index
					+ " does not exist");
		return (T) edgeArray[index];
	}

	@Override
	public int getEdgeCount() {
		return edgeCount;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Node> T getNode(String id) {
		return (T) nodeMap.get(id);
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T extends Node> T getNode(int index) {
		if (index < 0 || index > nodeCount)
			throw new IndexOutOfBoundsException("Node " + index
					+ " does not exist");
		return (T) nodeArray[index];
	}

	@Override
	public int getNodeCount() {
		return nodeCount;
	}

	// *** Iterators ***

	protected class EdgeIterator<T extends Edge> implements Iterator<T> {
		int iNext = 0;
		int iPrev = -1;

		public boolean hasNext() {
			return iNext < edgeCount;
		}

		@SuppressWarnings("unchecked")
		public T next() {
			if (iNext >= edgeCount)
				throw new NoSuchElementException();
			iPrev = iNext++;
			return (T) edgeArray[iPrev];
		}

		public void remove() {
			if (iPrev == -1)
				throw new IllegalStateException();
			removeEdge(edgeArray[iPrev], true, true, true);
			iNext = iPrev;
			iPrev = -1;
		}
	}

	protected class NodeIterator<T extends Node> implements Iterator<T> {
		int iNext = 0;
		int iPrev = -1;

		public boolean hasNext() {
			return iNext < nodeCount;
		}

		@SuppressWarnings("unchecked")
		public T next() {
			if (iNext >= nodeCount)
				throw new NoSuchElementException();
			iPrev = iNext++;
			return (T) nodeArray[iPrev];
		}

		public void remove() {
			if (iPrev == -1)
				throw new IllegalStateException();
			removeNode(nodeArray[iPrev], true);
			iNext = iPrev;
			iPrev = -1;
		}
	}

	@Override
	public <T extends Edge> Iterator<T> getEdgeIterator() {
		return new EdgeIterator<T>();
	}

	@Override
	public <T extends Node> Iterator<T> getNodeIterator() {
		return new NodeIterator<T>();
	}

	/*
	 * For performance tuning
	 * 
	 * @return the number of allocated but unused array elements public int
	 * getUnusedArrayElements() { int count = 0; count += edgeArray.length -
	 * edgeCount; count += nodeArray.length - nodeCount; for (ALNode n :
	 * this.<ALNode> getEachNode()) count += n.edges.length - n.degree; return
	 * count; }
	 */
}
