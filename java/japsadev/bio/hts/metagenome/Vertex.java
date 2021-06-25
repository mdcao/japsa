package japsadev.bio.hts.metagenome;

import java.util.ArrayList;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;

/**
 * This class models a vertex for my string graph which actually corresponds to an edge in SPAdes's assembly graph. 
 * Label for this vertex is extracted from its full name and used as its index. 
 * For example, vertex named 'EDGE_1_length_1000_cov_50' is labeled as vertex 1.
 * This vertex's neighborhood is described by the Edges incident to it. 
 * 
 * @author Son Nguyen
 * @date August 20, 2016
 */
public class Vertex {

    private ArrayList<Edge> neighborhood;
    private String fullName, label;
    private double coverage;
    private Sequence seq=null;
    //a sub-graph (path) is equivalent to a Vertex recursively
    private Path components=null;
    public Vertex(){
    	fullName="";
    	label="";
    	coverage=0;
        this.neighborhood = new ArrayList<Edge>();
        this.seq=new Sequence(Alphabet.DNA5(), 0);
    }
    /**
     * 
     * @param name The unique label associated with this Vertex
     */
    public Vertex(String name){
    	this.fullName=name;
        this.label=extractID(name);
        this.coverage=extractCoverage(name);
        this.neighborhood = new ArrayList<Edge>();
        this.seq=new Sequence(Alphabet.DNA5(), 0);
    }
    
    public Vertex(String name, Sequence seq){
        this(name);
        this.seq = seq;
    }
    
    public Vertex(Path p){
    	this();
    	fullName=label=p.getID();
    	components=p;
    	coverage=p.averageCov();
    }
    /**
     * 
     * @return If this vertex is unique (only based on its degree)
     */
    public boolean isUnique(){
    	return (getNeighborCount() <= 2);
    }
    public Path getSubComps(){
    	return components;
    }
    /**
     * Extract ID from name: EDGE_xx_length_yy_cov_zz;
     * @param name The name of Edge in assembly graph that correspond to this Vertex
     */
    private String extractID(String name){
    	return name.split("_")[1];
    }
    /**
     * Extract coverage from name: EDGE_xx_length_yy_cov_zz;
     * @param name The name of Edge in assembly graph that correspond to this Vertex
     * @return coverage value
     */
    private double extractCoverage(String name){
    	double res=1.0;
    	try{
    	 res=Double.parseDouble(name.split("_")[5]);
    	}catch(Exception e){
    		e.printStackTrace();
    	}
    	return res;
    }
    public double getCoverage(){
    	return coverage;
    }
    /**
     * This method adds an Edge to the incidence neighborhood of this graph iff
     * the edge is not already present. 
     * 
     * @param edge The edge to add
     */
    public void addNeighbor(Edge edge){
        if(this.neighborhood.contains(edge)){
            return;
        }
        this.neighborhood.add(edge);
    }
    
    
    /**
     * 
     * @param other The edge for which to search
     * @return true iff other is contained in this.neighborhood
     */
    public boolean containsNeighbor(Edge other){
        return this.neighborhood.contains(other);
    }
    
    /**
     * 
     * @param index The index of the Edge to retrieve
     * @return Edge The Edge at the specified index in this.neighborhood
     */
    public Edge getNeighbor(int index){
        return this.neighborhood.get(index);
    }
    
    
    /**
     * 
     * @param index The index of the edge to remove from this.neighborhood
     * @return Edge The removed Edge
     */
    Edge removeNeighbor(int index){
        return this.neighborhood.remove(index);
    }
    
    /**
     * 
     * @param e The Edge to remove from this.neighborhood
     */
    public void removeNeighbor(Edge e){
        this.neighborhood.remove(e);
    }
    
    
    /**
     * 
     * @return int The number of neighbors of this Vertex
     */
    public int getNeighborCount(){
        return this.neighborhood.size();
    }
    /**
     * 
     * @return String The label of this Vertex
     */
    public String getLabel(){
        return this.label;
    }
    /**
     * 
     * @return String The full name of this Vertex
     */
    public String getName(){
        return this.fullName;
    }
    /**
     * 
     * @param Sequence A sequence
     */
    public void setSequence(Sequence seq){
   		this.seq = seq;
    }
    /**
     * 
     * @return Sequence The sequence of this Vertex
     */
    public Sequence getSequence(){
        return this.seq;
    }
    /**
     * 
     * @return String A String representation of this Vertex
     */
    public String toString(){
        return "Vertex " + label;
    }
    
    /**
     * 
     * @return The hash code of this Vertex's label
     */
    public int hashCode(){
        return this.label.hashCode();
    }
    
    /**
     * 
     * @param other The object to compare
     * @return true iff other instanceof Vertex and the two Vertex objects have the same label
     */
    public boolean equals(Object other){
        if(!(other instanceof Vertex)){
            return false;
        }
        
        Vertex v = (Vertex)other;
        return this.label.equals(v.label);
    }
    
    /**
     * 
     * @return ArrayList<Edge> A copy of this.neighborhood. Modifying the returned
     * ArrayList will not affect the neighborhood of this Vertex
     */
    public ArrayList<Edge> getNeighbors(){
        return new ArrayList<Edge>(this.neighborhood);
    }
    
}


