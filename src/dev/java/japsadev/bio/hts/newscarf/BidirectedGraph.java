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
public class BidirectedGraph extends AdjacencyListGraph{
    private int kmer;
    static final int TOLERATE=500;

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
}
