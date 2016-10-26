package japsadev.bio.hts.newscarf;

import java.util.Iterator;
import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;

public class GraphExplore {
    public static void main(String args[]) {
        new GraphExplore();
    }

    public GraphExplore() {
        //Graph graph = new SingleGraph("tutorial 1");
    	System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
    	
    	Graph graph= new MultiGraph("tutorial 1");
        graph.addAttribute("ui.quality");
        graph.addAttribute("ui.antialias");
        graph.addAttribute("ui.stylesheet", styleSheet);
        graph.setAutoCreate(true);
        graph.setStrict(false);
        graph.display();

        
        Edge 	ab=graph.addEdge("AB", "A", "B"),
        		ab2=graph.addEdge("AB2", "A", "B");      
  
        //ab.setAttribute("ui.style", "text-offset: -10;"); 
        ab.addAttribute("ui.label", ab.getId());    
        
        ab2.addAttribute("ui.label", ab2.getId());
        
        graph.addEdge("BC", "B", "C");
        graph.addEdge("CA", "C", "A");
        graph.addEdge("AD", "A", "D");
        graph.addEdge("DE", "D", "E");
        graph.addEdge("DF", "D", "F");
        graph.addEdge("EF", "E", "F");
        
        System.out.println(graph.getEdgeCount());
        for (Node node : graph) {
            node.addAttribute("ui.label", node.getId());
            node.setAttribute("ui.style", "text-offset: -10;"); 
        }
        
        //explore(graph.getNode("A"));
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
        "}";
    	
}