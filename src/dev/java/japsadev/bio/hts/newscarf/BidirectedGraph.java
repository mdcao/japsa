package japsadev.bio.hts.newscarf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;


public class BidirectedGraph extends AdjacencyListGraph{
    static int kmer=127;
    static final int TOLERATE=500;
    static final int D_LIMIT=200;
    static final double ER_LOWERBOUND=.5,
    					ER_UPPERBOUND=1.5;
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
	
	//just to make AbstractGraph.removeEdge(AbstractEdge, boolean, boolean, boolean) visible
	protected void removeEdgeDup(AbstractEdge edge, boolean graphCallback,
			boolean sourceCallback, boolean targetCallback) {
		this.removeEdge(edge, graphCallback, sourceCallback, targetCallback);
	
	}
	protected BidirectedEdge addEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1){
		BidirectedEdge tmp = addEdge(BidirectedEdge.createID(src, dst, dir0, dir1), src, dst);
//		String s1=tmp.toString();
//		//tmp.setDir0(dir0);
//		//tmp.setDir1(dir1);
//		if(!s1.equals(tmp.toString()))
//			System.out.println(s1 + " ---> " + tmp);
		return tmp;
	}
	/**********************************************************************************
	 * ****************************Algorithms go from here*****************************
	 */
    public BidirectedGraph(){
    	this("Assembly graph",true,false,200000,1000000);
        setKmerSize(127);//default kmer size used by SPAdes to assembly MiSeq data
    }
    //TODO: read from ABySS assembly graph (graph of final contigs, not like SPAdes)
    public void loadFromFile(String graphFile) throws IOException{
        setAutoCreate(true);
        setStrict(false);
		//1. next iterate over again to read the connections
		SequenceReader reader = new FastaReader(graphFile);
		Sequence seq;
		int shortestLen = 10000;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length()<shortestLen)
				shortestLen=seq.length();
			
			String[] adjList = seq.getName().split(":");
			String name = adjList[0];
			boolean dir0=name.contains("'")?false:true;
			
			name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim(); //EDGE_X_length_Y_cov_Z
			
			String nodeID = name.split("_")[1];
			AbstractNode node = addNode(nodeID);
			node.setAttribute("name", name);
			
			if(dir0){
				seq.setName(name);
				node.setAttribute("seq", seq);
				double cov = Double.parseDouble(name.split("_")[5]);
				node.setAttribute("cov", cov);
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
					
					addEdge(node, nbr, dir0, dir1);
					//e.addAttribute("ui.label", e.getId());
				}
			}
			
		}

		//rough estimation of kmer used
		if((shortestLen-1) != getKmerSize())
			setKmerSize(shortestLen-1);
		
		reader.close();
    }
    /*
     * Read paths from contigs.path and reduce the graph
     */
    public void readPathsFromSpades(String paths) throws IOException{

		BufferedReader pathReader = new BufferedReader(new FileReader(paths));
		
		String s;
		//Read contigs from contigs.paths and refer themselves to contigs.fasta
		boolean flag=false;
		while((s=pathReader.readLine()) != null){
			if(s.contains("NODE")){
				flag=s.contains("'")?false:true;
				continue;
			}else if(flag){
				BidirectedPath path=new BidirectedPath(this, s);
//				System.out.println("Using path to reduce: " + path.getId());
//				System.out.println("Before reduce => Node: " + getNodeCount() + " Edge: " + getEdgeCount());
				
//				AbstractNode comp=

		    	this.reduce(path);


//				if(comp!=null){
//					System.out.println("Reverting node: " + comp.getId());
//					revert(comp);
//			        System.out.println("After revert => Node: " + getNodeCount() + " Edge: " + getEdgeCount());
//
//				}
			}	
				

		}
		pathReader.close();
    }
	
    public static int getKmerSize(){
    	return BidirectedGraph.kmer;
    }
    public static void setKmerSize(int kmer){
    	BidirectedGraph.kmer=kmer;
    }
    
    /**
     * 
     * @param p Path to be grouped as a virtually vertex
     */
    public AbstractNode reduce(BidirectedPath p){
    	//do nothing if the path has only one node
    	if(p==null || p.getEdgeCount()<1)
    		return null;
    	
    	//now only work with path containing more than 2 unique nodes
    	int uniqueCount=0;
    	for(Node n:p.getEachNode()){
    		if(isUnique(n))
    			uniqueCount++;
    	}
    	if(uniqueCount < 2)
    	{
    		System.out.println("ignore path with less than 1 unique contig!");
    		return null;
    	}
    	//add the new composite Node to the graph
    	//compare id from sense & anti-sense to get the unique one
    	AbstractNode comp = addNode(p.getId().compareTo(p.getReversedComplemented().getId())>0?
    								p.getReversedComplemented().getId():p.getId());
    	
    	comp.addAttribute("path", p);
    	comp.addAttribute("seq", p.spelling());
        comp.addAttribute("ui.label", comp.getId());
        comp.setAttribute("ui.style", "text-offset: -10;"); 
        comp.setAttribute("ui.class", "marked");
        try { Thread.sleep(100); } catch (Exception e) {}

    	//store unique nodes on p for removing
    	ArrayList<String> tobeRemoved=new ArrayList<String>();
    	for(Node n:p.getEachNode()){
    		if(isUnique(n))
    			tobeRemoved.add(n.getId());
    	}
    	BidirectedNode 	start = (BidirectedNode) p.getRoot(),
    					end = (BidirectedNode) p.peekNode();
    	boolean startDir = ((BidirectedEdge) p.getEdgePath().get(0)).getDir(start), 
    			endDir = ((BidirectedEdge) p.peekEdge()).getDir(end);
    	//set neighbors of the composite Node
    	Iterator<Edge> startEdges = startDir?start.getEnteringEdgeIterator():start.getLeavingEdgeIterator(),
    					endEdges = endDir?end.getEnteringEdgeIterator():end.getLeavingEdgeIterator();
    	while(startEdges.hasNext()){
    		BidirectedEdge e = (BidirectedEdge) startEdges.next();
    		BidirectedNode opNode = e.getOpposite(start);
    		boolean opDir = e.getDir(opNode);
    		//Edge tmp=
    		addEdge(BidirectedEdge.createID(comp, opNode, false, opDir), comp, opNode);//always into start node
    		//System.out.println("From " + start.getId() + ": " + tmp.getId() + " added!");
    	}
    	
    	while(endEdges.hasNext()){
    		BidirectedEdge e = (BidirectedEdge) endEdges.next();
    		BidirectedNode opNode = e.getOpposite(end);
    		boolean opDir = e.getDir(opNode);
    		//Edge tmp=
    		addEdge(BidirectedEdge.createID(comp, opNode, true, opDir), comp, opNode);//always out of end node
    	
    		//System.out.println("From " + end.getId() + ": " + tmp.getId() + " added!");

    	}

    	for(String nLabel:tobeRemoved){
    		//System.out.println("About to remove " + nLabel);
    		removeNode(nLabel);
    	}
    		
    	//TODO: remove bubbles...
    	return comp;
    }
    
    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param p Path to simplify the graph
     */
    public void reduce(BidirectedPath p, boolean keepTrack){
    	//do nothing if the path has only one node
    	if(p==null || p.getEdgeCount()<1)
    		return;
    	
    	//loop over the edges of path (like spelling())
    	BidirectedNode markerNode = null;
    	//search for an unique node as the marker. 
    	for(Edge e:p.getEdgePath()){
 
			if(markerNode==null){
				
			}
	    	//as the next marker coming, connect it to 
	    	//the previous marker by an edge bind with corresponding path. TODO:check conflict!!! 
			else{
				
			}
		}
    	
    	//remove appropriate edges
    	
    }
    
    /**
     * 
     * @param v Node to be reverted (1-level reverting)
     */
    public void revert(AbstractNode v){
    	System.out.println("Reverting...");
    	Path p=v.getAttribute("path");
    	if(p==null) return;
    	
    	BidirectedNode 	start = (BidirectedNode) p.getRoot(),
    					end = (BidirectedNode) p.peekNode();
    	boolean startDir = ((BidirectedEdge) p.getEdgePath().get(0)).getDir(start), 
    			endDir = ((BidirectedEdge) p.peekEdge()).getDir(end);
    	
    	//add back all neighbor edges of this composite vertex
    	Iterator<Edge> 	startEdges = v.getEnteringEdgeIterator(),
						endEdges = v.getLeavingEdgeIterator();
    	//add back all nodes from the path
		for(Node n:p.getNodeSet()){
			if(getNode(n.getId())!=null)
				continue;
			Node tmp = addNode(n.getId());
			tmp.addAttribute("seq", (japsa.seq.Sequence)n.getAttribute("seq"));
			tmp.addAttribute("name", (String)n.getAttribute("name"));
			tmp.addAttribute("path", (BidirectedPath)n.getAttribute("path"));

			//System.out.println("Adding back edge "+tmp.getId());
		}
		while(startEdges.hasNext()){
			BidirectedEdge e = (BidirectedEdge) startEdges.next();
			BidirectedNode opNode = e.getOpposite(v);
			boolean opDir = e.getDir(opNode);
			//Edge tmp = 
			addEdge(BidirectedEdge.createID(start, opNode, !startDir, opDir), start, opNode);
			//System.out.println("Adding back edge "+tmp.getId());
		}
		
		while(endEdges.hasNext()){
			BidirectedEdge e = (BidirectedEdge) endEdges.next();
			BidirectedNode opNode = e.getOpposite(v);
			boolean opDir = e.getDir(opNode);
			//Edge tmp = 
			addEdge(BidirectedEdge.createID(end, opNode, !endDir, opDir), end, opNode);
			//System.out.println("Adding back edge "+tmp.getId());
		}

    	//add back all edges from the path
		for(Edge e:p.getEdgeSet()){
			//Edge tmp = 
			addEdge(e.getId(), e.getSourceNode().getId(), e.getTargetNode().getId());
			//System.out.println("Adding back edge "+tmp.getId());
		}
    	//finally remove the composite node
    	removeNode(v);
	}
    
    /*
     * This function deduces a full path in this graph between 2 nodes aligned with a long read
     */
    protected ArrayList<BidirectedPath> getClosestPath(Alignment from, Alignment to, int distance){
    	BidirectedNode srcNode = from.node,
    					dstNode = to.node;
    	System.out.println("Looking for path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
    	BidirectedPath 	tmp = new BidirectedPath();
    	ArrayList<BidirectedPath>	possiblePaths = new ArrayList<BidirectedPath>(),
    								retval = new ArrayList<BidirectedPath>();
    	tmp.setRoot(srcNode);  	
    	
    	//traverse(tmp, dest, retval, distance+source.getSeq().length()+dest.getSeq().length());
    	traverse(tmp, dstNode, possiblePaths, distance, from.strand, to.strand);
    	//only get the best ones
    	if(possiblePaths.isEmpty())
    		return null;
    	double bestScore=possiblePaths.get(0).getDeviation();
    	for(int i=0;i<possiblePaths.size();i++){
    		BidirectedPath p = possiblePaths.get(i);
    		if(p.getDeviation() != bestScore)
    			break;
    		retval.add(p);
    	}
    	
    	return retval;
    	
    }
    private void traverse(BidirectedPath path, BidirectedNode dst, ArrayList<BidirectedPath> curResult, int distance, boolean srcDir, boolean dstDir){
    	BidirectedNode currentNode=(BidirectedNode) path.peekNode();
    	BidirectedEdge currentEdge;
    	boolean curDir;//direction to the next node, = ! previous'
    	
    	Iterator<BidirectedEdge> ite;
    	if(path.size() <= 1) //only root
			curDir=srcDir;//re-check
		else{
			currentEdge = (BidirectedEdge) path.peekEdge();
			curDir = !((BidirectedEdge) currentEdge).getDir(currentNode);
		}
		ite=curDir?currentNode.getLeavingEdgeIterator():currentNode.getEnteringEdgeIterator();

    	while(ite.hasNext()){
    		BidirectedEdge e = ite.next();
			path.add(e);

			if(e.getOpposite(currentNode)==dst && e.getDir(dst)!=dstDir && Math.abs(distance+getKmerSize()) < TOLERATE){
		    	BidirectedPath 	curPath=curResult.isEmpty()?new BidirectedPath():curResult.get(0), //the best path saved among all possible paths from the list curResult
		    					tmpPath=new BidirectedPath(path);
		    	tmpPath.setDeviation(Math.abs(distance+getKmerSize()));
		    	if(	Math.abs(distance+getKmerSize()) < curPath.getDeviation() )
		    		curResult.add(0, tmpPath);
		    	else
		    		curResult.add(tmpPath);
				
				System.out.println("Hit added: "+path.getId()+"(candidate deviation: "+Math.abs(distance+getKmerSize())+")");
			}else{
				int newDistance=distance-((Sequence) e.getOpposite(currentNode).getAttribute("seq")).length()+getKmerSize();
				if (newDistance+getKmerSize()<-TOLERATE){
					System.out.println("Stop following path with distance "+newDistance+" already! : "+path.getId());
				}else
					traverse(path, dst, curResult, newDistance, srcDir, dstDir);
			}
			path.popNode();
    	
    	}
    }
    
    /*
     * Find a path based on list of Alignments
     */
	public BidirectedPath pathFinding(ArrayList<Alignment> alignments) {
		if(alignments.size()<=1)
			return null;
		
		System.out.println("=================================================");
		for(Alignment alg:alignments)
			System.out.println("\t"+alg.toString());
		System.out.println("=================================================");
		//First bin the alignments into different overlap regions			
		//only considering useful alignments
		HashMap<Range,Alignment> allAlignments = new HashMap<Range,Alignment>();
		ArrayList<BidirectedPath> joinPaths = new ArrayList<BidirectedPath>();
		
		for(Alignment alg:alignments){
			if(alg.useful){
				Range range = new Range(alg.readAlignmentStart(),alg.readAlignmentEnd());
				allAlignments.put(range, alg);
			}
		}
		//now get all the bin in order
		List<Range> baseRanges=new ArrayList<Range>(allAlignments.keySet());
		List<List<Range>> rangeGroups = MetaRange.getOverlappingGroups(baseRanges);
		
		if(rangeGroups.size() < 2)
			return null;
		
		System.out.println("Binning ranges: ");
	    for(List<Range> group : rangeGroups){
	        System.out.println(group);
	    }

		//iterate all alignments in adjacent bins to find correct path
		ArrayList<Alignment> curGroup = new ArrayList<Alignment>(), 
							nextGroup = new ArrayList<Alignment>();
		List<Range> curRanges = rangeGroups.get(0);
		for(Range r:curRanges)
			curGroup.add(allAlignments.get(r));
		
		for(int i=1; i<rangeGroups.size();i++){
			List<Range> nextRanges = rangeGroups.get(i);
			nextGroup = new ArrayList<Alignment>();
			
			for(Range r:nextRanges)
				nextGroup.add(allAlignments.get(r));
			
			ArrayList<BidirectedPath> allPaths = new ArrayList<BidirectedPath>();
			
			
			for(Alignment curAlg:curGroup){
				for(Alignment nextAlg:nextGroup){
					int distance = nextAlg.readAlignmentStart()-curAlg.readAlignmentEnd();
					if(distance<D_LIMIT){
						ArrayList<BidirectedPath> bridges = getClosestPath(curAlg, nextAlg, distance);
						if(bridges!=null)
							allPaths.addAll(bridges);
					}
				}
			}
			//join all paths from previous to the new ones
			//TODO:optimize it
			if(joinPaths.isEmpty())
				joinPaths=allPaths;
			else{
				System.out.println("=====Current list of paths: " + joinPaths);
				System.out.println("=====Join to list of paths: " + allPaths);

				for(BidirectedPath p:joinPaths)
					for(BidirectedPath e:allPaths)
						p.join(e);			
				
			}
			curGroup=nextGroup;
		}
		
		for(BidirectedPath path:joinPaths)
			System.out.println("A member Path: " + path.toString() + " deviation: " + path.getDeviation());
		if(joinPaths.isEmpty())
			return null;
		else
			return joinPaths.get(0);
	}
	
    /*
     * Important function: determine if a node is able to be removed or not
     * TODO: re-implement it based on statistics of coverage also
     * 1. pick the least coverage ones among a path as the base
     * 2. global base
     */
    public static boolean isUnique(Node node){
    	boolean res = false;
    	if(node.getDegree()<=2){
//    		if(((Sequence)node.getAttribute("seq")).length() > 5000 || node.getDegree()==0)
    			res=true;
    	}
    		
    	return res;
    }

}
