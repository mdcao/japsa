package japsadev.bio.hts.scaffold;

import japsa.util.Logging;

/**
 * This class models an bidirected Edge in my Graph implementation.
 * An Edge contains two vertices and a weight (distance between them).
 * A certain edge (v1,v2) can take one among 4 types: ++, --, +- and -+. Each
 * type corresponds to the way we read the DNA sequence in each read when traversing
 * this edge. 
 * For example: v1->---<-v2 or (v1,v2)+- spells out (v1 v2') and/or (v2 v1') as in SPAdes output.
 * This class also deviates from the expectations of the Comparable interface
 * in that a return value of 0 does not indicate that this.equals(other). The
 * equals() method only compares the vertices, while the compareTo() method 
 * compares the edge weights. This provides more efficient implementation for
 * checking uniqueness of edges, as well as the fact that two edges of equal weight
 * should be considered equitably in a path finding or spanning tree algorithm.
 * 
 * @author Son Nguyen
 * @date August 20, 2016
 */
public class Edge implements Comparable<Edge> {

    private Vertex one, two;
    private boolean dOne, dTwo;
    private int weight;
    
    /**
     * 
     * @param one The first vertex in the Edge
     * @param two The second vertex in the Edge
     */
    public Edge(Vertex one, Vertex two, boolean d1, boolean d2){
        this(one, two, d1, d2, -Graph.getKmerSize());
    }
    
    /**
     * 
     * @param one The first vertex in the Edge
     * @param two The second vertex of the Edge
     * @param weight The weight of this Edge
     */
    public Edge(Vertex one, Vertex two, boolean dOne, boolean dTwo, int weight){
        //this.one = (one.getLabel().compareTo(two.getLabel()) <= 0) ? one : two;
        //this.two = (this.one == one) ? two : one;
        this.one=one;
        this.two=two;
    	this.weight = weight;
        this.dOne=dOne;
        this.dTwo=dTwo;
    }
    
    
    /**
     * 
     * @param current
     * @return The neighbor of current along this Edge
     */
    public Vertex getNeighbor(Vertex current){
        if(!(current.equals(one) || current.equals(two))){
            return null;
        }
        
        return (current.equals(one)) ? two : one;
    }
    /**
     * Return the same Edge but reading the other way around 
     * just swap the order of its vertices upside down
     * @param 
     * @return the identical Edge
     */
    public Edge getReversedRead(){
    	return new Edge(this.two, this.one, !this.dTwo, !this.dOne, this.weight);
    }
    /**
     * 
     * @param current
     * @return The direction to spell *current* along this Edge
     */
    public boolean getDirection(Vertex current){
        assert (current.equals(one) || current.equals(two)):"Vertex doesn't belong to this Edge!";

        return (current.equals(one)) ? dOne : !dTwo;
    }
    
    /**
     * 
     * @return Vertex this.one
     */
    public Vertex getOne(){
        return this.one;
    }
    
    /**
     * 
     * @return Vertex this.two
     */
    public Vertex getTwo(){
        return this.two;
    }
    
    /**
     * 
     * @return boolean this.dOne
     */
    public boolean getDOne(){
        return this.dOne;
    }
    
    /**
     * 
     * @return boolean this.dTwo
     */
    public boolean getDTwo(){
        return this.dTwo;
    }
    /**
     * 
     * @return int The weight of this Edge
     */
    public int getWeight(){
        return this.weight;
    }
    
    
    /**
     * 
     * @param weight The new weight of this Edge
     */
    public void setWeight(int weight){
        this.weight = weight;
    }
    
    
    /**
     * Note that the compareTo() method deviates from 
     * the specifications in the Comparable interface. A 
     * return value of 0 does not indicate that this.equals(other).
     * The equals() method checks the Vertex endpoints, while the 
     * compareTo() is used to compare Edge weights
     * 
     * @param other The Edge to compare against this
     * @return int this.weight - other.weight
     */
    public int compareTo(Edge other){
        return this.weight - other.weight;
    }
    
    /**
     * 
     * @return String A String representation of this Edge
     */
    public String toString(){
        return "({" + one + (dOne?"":"'") + ", " + two + (dTwo?"":"'") + "}, " + weight + ")";
    }
    
    /**
     * 
     * @return int The hash code for this Edge 
     */
    public int hashCode(){
        return (one.getLabel() + (dOne?"+":"-") + two.getLabel() + (dTwo?"+":"-")).hashCode(); 
    }
    
    /**
     * 
     * @param other The Object to compare against this
     * @return true iff other is an Edge with the same Vertices as this
     */
    public boolean equals(Object other){
        if(!(other instanceof Edge)){
            return false;
        }
        
        Edge e = (Edge)other;
        
        return (e.one.equals(this.one) && e.two.equals(this.two) && (e.getDOne()==this.dOne) && (e.getDTwo()==this.dTwo))
        		|| (e.one.equals(this.two) && e.two.equals(this.one) && (e.getDOne()!=this.dOne) && (e.getDTwo()!=this.dTwo));
    }   
}


