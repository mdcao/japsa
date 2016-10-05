package japsadev.bio.hts.scaffold;


import japsadev.bio.hts.scaffold.RealtimeScaffolding;
import japsadev.bio.hts.scaffold.ScaffoldGraph;

/**
*
* @author Son Nguyen
* @date August 21, 2016
*/
public class DemoGraph {
   
   public static void main(String[] args){
//       	Graph graph;
		try {
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
		
			
			String prefix = "/home/s.hoangnguyen/Projects/scaffolding/data/test/EcK12S";
			String bamFile = "/home/s.hoangnguyen/Projects/scaffolding/data/test/EcK12S.sam";
	
			
			
			ScaffoldGraph.minContigLength = 300;
			ScaffoldGraph.minSupportReads = 1;	
			ScaffoldGraph.maxRepeatLength = 7500;		
			//ScaffoldGraph.marginThres = marginThres;
			ScaffoldGraph.verbose = true;
			ScaffoldGraph graph;
			
			RealtimeScaffolding rtScaffolding = null;
			rtScaffolding = new RealtimeScaffolding("/home/s.hoangnguyen/Projects/scaffolding/data/test/contigs.fasta", null, null, null, null, "-");
			graph = rtScaffolding.graph;
			graph.prefix=prefix;

			synchronized(graph){
				graph.readMore("/home/s.hoangnguyen/Projects/scaffolding/data/test/assembly_graph.fastg","/home/s.hoangnguyen/Projects/scaffolding/data/test/contigs.paths");
			}
			
			rtScaffolding.scaffolding(bamFile, 500, 5, 0, 1);
				
		}catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
       
      
   }
}


