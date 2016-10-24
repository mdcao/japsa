package japsadev.bio.hts.metagenome;

import java.util.ArrayList;
import java.util.HashSet;

/**
*
* @author Son Nguyen
* @date August 21, 2016
*/
public class DemoGraph {
   
   public static void main(String[] args){
       	Graph graph = new Graph();
//		try {
//			graph = new Graph("/home/s.hoangnguyen/Projects/scaffolding/data/spades_3.7/EcK12S-careful/assembly_graph.fastg");
//			HashSet<Vertex> vertices = (HashSet<Vertex>) graph.getVertices();
//			
//			//display the initial setup- all vertices adjacent to each other
//			for(Vertex vertex:vertices){
//				System.out.println(vertex);
//	           
//				for(int j = 0; j < vertex.getNeighborCount(); j++){
//					System.out.println(vertex.getNeighbor(j));
//				}
//	           
//				System.out.println();
//			}
//			
//			Node 	source=new Node(graph.getVertex("108"),false),
//					dest=new Node(graph.getVertex("201"),true);
//			int distance=2962;
//			
//			ArrayList<Path> paths=graph.DFS(source, dest, distance);
//	    	if(!paths.isEmpty()){
//	    		System.out.println("Paths found ("+distance+"):");
//	    		for(Path p:paths)
//	    			System.out.println(p.toString() + " d=" + p.getDeviation() );
//	    	}
//	    	else
//	    		System.out.println("Path not found ("+distance+")!");
//				
//		}catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
       
      Vertex 	v1= new Vertex("EDGE_1_length_1000_cov_9"),
    		  	v2= new Vertex("EDGE_2_length_200_cov_22"),
    		  	v3= new Vertex("EDGE_3_length_300_cov_12"),
    		  	v4= new Vertex("EDGE_4_length_500_cov_33"),
    		  	v5= new Vertex("EDGE_5_length_600_cov_21");
      graph.addVertex(v1, false);
      graph.addVertex(v2, false);
      graph.addVertex(v3, false);
      graph.addVertex(v4, false);
      graph.addVertex(v5, false);
      
      graph.addEdge(v1, v2, true, true);
      graph.addEdge(v1, v2, false, false);
      
      graph.addEdge(v2, v1, true, true);
      graph.addEdge(v2, v3, true, true);
      graph.addEdge(v2, v1, false, false);
      graph.addEdge(v2, v4, false, false);
      
      graph.addEdge(v3, v4, true, true);
      graph.addEdge(v3, v2, false, false);
      
      graph.addEdge(v4, v5, true, true);
      graph.addEdge(v4, v2, true, true);
      graph.addEdge(v4, v3, false, false);
      graph.addEdge(v4, v5, false, false);
      
      graph.addEdge(v5, v4, false, false);
      graph.addEdge(v5, v4, true, true);
      
      graph.printStats();
      
      Path 	p1=new Path(graph, "1+,2+,3+"),
    		p2=new Path(graph, "3+,4+");
      graph.reduce(p1);
      graph.printStats();
      graph.reduce(p2);
      graph.printStats();
   }
}


