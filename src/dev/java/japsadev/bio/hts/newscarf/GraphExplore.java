package japsadev.bio.hts.newscarf;

import java.io.IOException;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.*;
public class GraphExplore {
	public static String spadesFolder="/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/";
	
	public static void main(String args[]) {
    	try {
			new GraphExplore();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


    }

    public GraphExplore() throws IOException{
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 

        HybridAssembler ass = new HybridAssembler(spadesFolder+"EcK12S-careful/assembly_graph.fastg");
    	BidirectedGraph graph= ass.simGraph;
        graph.addAttribute("ui.quality");
        graph.addAttribute("ui.antialias");
        graph.addAttribute("ui.stylesheet", styleSheet);
        graph.display();
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());

        
        for (Node node : graph) {
            node.addAttribute("ui.label", node.getId());
            node.setAttribute("ui.style", "text-offset: -10;"); 
            if(BidirectedGraph.isUnique(node))
            	node.setAttribute("ui.class", "marked");
        }
        
        //explore(graph.getNode("A"));
        
        /*
         * Testing reduce function
         */
        try {
			ass.reduceFromSPAdesPaths(spadesFolder+"EcK12S-careful/contigs.paths");
			ass.assembly(spadesFolder+"bwa/EcK12S-careful.sam", 50);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
//        //remove dead-end nodes
//        for (Node node : graph) {
//            if(node.getInDegree() == 0 || node.getOutDegree() == 0)
//            	graph.removeNode(node);
//        }
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());

        /*
         * Testing BidirectedEdge id pattern
         */
//    	String pattern = "^\\[([0-9\\+\\-]*)\\]([oi])\\[([0-9\\+\\-]*)\\]([oi])$";
//        // Create a Pattern object
//        Pattern r = Pattern.compile(pattern);
//        // Now create matcher object.
//        String id="[3]i[4+8+]o";
//        Matcher m = r.matcher(id);
//         	
//        if(m.find()){
//        	System.out.println(m.group(1)+"|"+m.group(2)+"|"+m.group(3)+"|"+m.group(4));
//        } else
//        	System.out.println("Fuck");
    }

    public void explore(Node source) {
        Iterator<? extends Node> k = source.getBreadthFirstIterator();

        while (k.hasNext()) {
            Node next = k.next();
            next.setAttribute("ui.class", "marked");
            sleep();
        }
    }

    protected void sleep() {
        try { Thread.sleep(1000); } catch (Exception e) {}
    }

    protected String styleSheet =
        "node {" +
        "	fill-color: black;" +
        "}" +
        "node.marked {" +
        "	fill-color: red;" +
        "}" +
        "edge.marked {" +
        "	fill-color: red;" +
        "}";
    	
}