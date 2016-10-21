package japsadev.bio.hts.metagenome;

import java.io.IOException;
import java.util.*;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

/**
 * This class models a simple, bidirected graph using an 
 * incidence list representation. Vertices are identified 
 * uniquely by their labels, and only unique vertices are allowed.
 * At most one unique Edge per vertex pair is allowed in this Graph.
 * 
 * @author Son Nguyen
 * @date August 20, 2016
 */
public class Graph {
    
    private HashMap<String, Vertex> vertices;
    private HashMap<Integer, Edge> edges;
    private int kmer;

    static final int TOLERATE=500;

    
    public Graph(){
        this.vertices = new HashMap<String, Vertex>();
        this.edges = new HashMap<Integer, Edge>();
        setKmerSize(127);//default kmer size used by SPAdes to assembly MiSeq data
    }
    
    
    public Graph(String graphFile) throws IOException{
    	this();
		//1. next iterate over again to read the connections
		SequenceReader reader = new FastaReader(graphFile);
		Sequence seq;
		int shortestLen = 10000;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length()<shortestLen)
				shortestLen=seq.length();
			
			String[] adjList = seq.getName().split(":");
			String name = adjList[0];
			boolean dir1=name.contains("'")?false:true;
			
			name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim();
			
			//FIXME: constructor is invoked by name but hashmap is based on label!!!
			Vertex current=new Vertex(name);
			if(getVertex(current.getLabel())!=null)
				current=getVertex(current.getLabel());
				
			addVertex(current, false);
			
			if(dir1){
				seq.setName(name);
				current.setSequence(seq);
				//System.out.println(current);
			}
			if (adjList.length > 1){
				String[] nbList = adjList[1].split(",");
				for(int i=0; i < nbList.length; i++){
					// create list of bridges here (distance=-kmer overlapped)
					String neighbor = nbList[i];
					boolean dir2=neighbor.contains("'")?false:true;
					neighbor=neighbor.replaceAll("[^a-zA-Z0-9_.]", "").trim();

					Vertex nbVertex=new Vertex(neighbor);
					if(getVertex(nbVertex.getLabel())!=null)
						nbVertex=getVertex(nbVertex.getLabel());
	
					addVertex(nbVertex, false);
					
					addEdge(current, nbVertex, dir1, dir2);
				}
			}
			
		}
		//rough estimation of kmer used
		if((shortestLen-1) != getKmerSize())
			setKmerSize(shortestLen-1);
		
		reader.close();
    }
    /**
     * This constructor accepts an ArrayList<Vertex> and populates
     * this.vertices. If multiple Vertex objects have the same label,
     * then the last Vertex with the given label is used. 
     * 
     * @param vertices The initial Vertices to populate this Graph
     */
    public Graph(ArrayList<Vertex> vertices){
        this.vertices = new HashMap<String, Vertex>();
        this.edges = new HashMap<Integer, Edge>();
        
        for(Vertex v: vertices){
            this.vertices.put(v.getLabel(), v);
        }
        setKmerSize(127);//default kmer size used by SPAdes to assembly MiSeq data

    }
    
    public int getKmerSize(){
    	return this.kmer;
    }
    public void setKmerSize(int kmer){
    	this.kmer=kmer;
    	Path.setK(kmer);
    }
    /**
     * This method adds am edge between Vertices one and two
     * and their corresponding direction of weight kmer, 
     * if no Edge between these Vertices already exists in the Graph.
     * 
     * @param one The first vertex to add
     * @param two The second vertex to add
     * @param d1 The direction on the side of vertex one
     * @param d2 The direction on the side of vertex two
     * @return true iff no Edge relating one and two exists in the Graph
     */
    public boolean addEdge(Vertex one, Vertex two, boolean d1, boolean d2){
        return addEdge(one, two, d1, d2, -kmer);
    }
    
    
    /**
     * Accepts two vertices, their directions and a weight, and adds the edge 
     * ({one, two}, {d1, d2}, weight) iff no Edge relating one and two 
     * exists in the Graph.
     * 
     * @param one The first Vertex of the Edge
     * @param two The second Vertex of the Edge
     * @param d1 The direction on the side of vertex one
     * @param d2 The direction on the side of vertex two
     * @param weight The weight of the Edge
     * @return true iff no Edge already exists in the Graph
     */
    public boolean addEdge(Vertex one, Vertex two, boolean d1, boolean d2, int weight){
       
        //ensures the Edge is not in the Graph
        Edge e = new Edge(one, two, d1, d2, weight);
        if(edges.containsKey(e.hashCode()) || edges.containsKey(e.getReversedRead().hashCode())){
            return false;
        }
       
        //and that the Edge isn't already incident to one of the vertices
        else if(one.containsNeighbor(e) || two.containsNeighbor(e.getReversedRead())){
            return false;
        }
            
        edges.put(e.hashCode(), e);
        one.addNeighbor(e);
        two.addNeighbor(e.getReversedRead());
        return true;
    }
    
    /**
     * 
     * @param e The Edge to look up
     * @return true iff this Graph contains the Edge e
     */
    public boolean containsEdge(Edge e){
        if(e.getOne() == null || e.getTwo() == null){
            return false;
        }
        
        return 	this.edges.containsKey(e.hashCode())
        		|| this.edges.containsKey(e.getReversedRead().hashCode());
    }
    
    
    /**
     * This method removes the specified Edge from the Graph,
     * including as each vertex's incidence neighborhood.
     * 
     * @param e The Edge to remove from the Graph
     * @return Edge The Edge removed from the Graph
     */
    public Edge removeEdge(Edge e){
       e.getOne().removeNeighbor(e);
       e.getTwo().removeNeighbor(e.getReversedRead());
       Edge rmEdge = this.edges.remove(e.hashCode());
       if (rmEdge==null)
    	   rmEdge = this.edges.remove(e.getReversedRead().hashCode());
       return rmEdge;
    }
    
    /**
     * 
     * @param vertex The Vertex to look up
     * @return true iff this Graph contains vertex
     */
    public boolean containsVertex(Vertex vertex){
        return this.vertices.get(vertex.getLabel()) != null;
    }
    
    /**
     * 
     * @param label The specified Vertex label
     * @return Vertex The Vertex with the specified label
     */
    public Vertex getVertex(String label){
        return vertices.get(label);
    }
    
    /**
     * This method adds a Vertex to the graph. If a Vertex with the same label
     * as the parameter exists in the Graph, the existing Vertex is overwritten
     * only if overwriteExisting is true. If the existing Vertex is overwritten,
     * the Edges incident to it are all removed from the Graph.
     * 
     * @param vertex
     * @param overwriteExisting
     * @return true iff vertex was added to the Graph
     */
    public boolean addVertex(Vertex vertex, boolean overwriteExisting){
        Vertex current = this.vertices.get(vertex.getLabel());
        if(current != null){
            if(!overwriteExisting){
                return false;
            }
            
            while(current.getNeighborCount() > 0){
                this.removeEdge(current.getNeighbor(0));
            }
        }
        
        
        vertices.put(vertex.getLabel(), vertex);
        return true;
    }
    
    /**
     * 
     * @param label The label of the Vertex to remove
     * @return Vertex The removed Vertex object
     */
    public Vertex removeVertex(String label){
        Vertex v = vertices.remove(label);
        
        while(v.getNeighborCount() > 0){
            this.removeEdge(v.getNeighbor((0)));
        }
        
        return v;
    }
    
    /**
     * 
     * @return Set<Vertex> All Graph's Vertex objects
     */
    public Set<Vertex> getVertices(){
        return new HashSet<Vertex>(this.vertices.values());
    }
    
    /**
     * 
     * @return Set<Edge> The Edges of this graph
     */
    public Set<Edge> getEdges(){
        return new HashSet<Edge>(this.edges.values());
    }
    
    /**
     * Find a path between two nodes within a given distance
     */
    public ArrayList<Path> DFS(Node source, Node dest, int distance){
    	System.out.println("Looking for path between " + source.toString() + " to " + dest.toString() + " with distance " + distance);
    	Path 	tmp = new Path();
    	ArrayList<Path>	retval = new ArrayList<Path>();
    	tmp.addNode(source);  	
    	
    	//traverse(tmp, dest, retval, distance+source.getSeq().length()+dest.getSeq().length());
    	traverse(tmp, dest, retval, distance);

    	return retval;
    }
    
    public void traverse(Path path, Node dest, ArrayList<Path> curResult, int distance){
    	Node source=path.getEnd();
    	assert source!=null:"Path null fault!";
    	
    	ArrayList<Edge> nList = source.getVertex().getNeighbors();
    	for(Edge e:nList){
    		if(e.getDOne()==source.getDirection()){
    			path.addNode(e.getTwo(), e.getDTwo());

    			if(e.getTwo()==dest.getVertex() && e.getDTwo()==dest.getDirection() && Math.abs(distance+getKmerSize()) < TOLERATE){

    		    	Path 	curPath=curResult.isEmpty()?new Path():curResult.get(0), //the best path saved among all possible paths from the list curResult
    		    			tmpPath=new Path();
    		    	tmpPath.setComp(path.getNodes());
    		    	tmpPath.setDeviation(Math.abs(distance+getKmerSize()));
    		    	if(	Math.abs(distance+getKmerSize()) < curPath.getDeviation() )
    		    		curResult.add(0, tmpPath);
    		    	else
    		    		curResult.add(tmpPath);
    				
    				System.out.println("Hit added: "+path+"(candidate deviation: "+Math.abs(distance+getKmerSize())+")");
    			}else{
    				int newDistance=distance-e.getTwo().getSequence().length()+getKmerSize();
    				if (newDistance+getKmerSize()<-TOLERATE){
    					System.out.println("Stop following path with distance "+newDistance+" already! : "+path);
    				}else
    					traverse(path, dest, curResult, newDistance);
    			}
    			path.removeLast();
    		}
    	}
    }
    /**
     * 
     * @param p Path to be grouped as a virtually vertex
     */
    public void reduce(Path p){
    	Vertex comp=new Vertex(p);
    	Node 	start = p.getStart(),
    			end = p.getEnd();
    	//set neighbors of the grouped Vertex
    	for(Edge e:start.getVertex().getNeighbors()){
    		if(e.getDOne()!=start.getDirection())
    			comp.addNeighbor(e);
    	}
    	for(Edge e:end.getVertex().getNeighbors()){
    		if(e.getDOne()!=end.getDirection())
    			comp.addNeighbor(e);
    	}
    	//remove unique nodes on p
    	for(Node n:p.getNodes()){
    		if(n.getVertex().isUnique())
    			removeVertex(n.getVertex().getLabel());
    	}
    	//add the new composite Vertex to the graph
    	addVertex(comp, true);
    	
    	//TODO: remove bubbles...
    }
    /**
     * 
     * @param v Vertex to be reverted (1-level reverting)
     */
    public void revert(Vertex v){
    	//TODO: revert to initial status by extracting a complex vertex into its initial components
    	Path p=v.getSubComps();
    	if(!containsVertex(v)||p==null) return;
    	//add back all vertices first
    	for(Node n:p.getNodes())
    		addVertex(n.getVertex(), false);
    	//then add back all neighbor edges of this composite vertex
    	for(Edge e:v.getNeighbors())
    		addEdge(v,e.getTwo(),e.getDOne(),e.getDTwo());
    	//finally add back all edges from the path
    	Node prev=p.getStart();
    	for(Node cur:p.getNodes()){
    		if(cur==p.getStart())
    			continue;
    		else{
    			addEdge(prev.getVertex(),cur.getVertex(),prev.getDirection(),cur.getDirection());
    			prev=cur;
    		}
    	}
    	
    	//remove the original composite vertex
    	removeVertex(v.getLabel());
    }
    public void printStats(){
    	System.out.println(vertices.size() + " vertices:");
    	for(String label:vertices.keySet())
    		System.out.print(label+", ");
    	System.out.println(edges.size() + " edges:");
    	for(Edge e:edges.values()){
    		System.out.print(e.toString()+", ");
    	}
    }
}


