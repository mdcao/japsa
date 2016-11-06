package japsadev.bio.hts.newscarf;

import java.util.Iterator;
import org.graphstream.graph.*;
public class GraphExplore {
    public static void main(String args[]) {
    	new GraphExplore();

    }

    public GraphExplore(){
        //Graph graph = new SingleGraph("tutorial 1");
    	System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
    	
    	BidirectedGraph graph= new BidirectedGraph();
        graph.addAttribute("ui.quality");
        graph.addAttribute("ui.antialias");
        graph.addAttribute("ui.stylesheet", styleSheet);
        graph.display();
        
        graph.loadFromFile("/home/hoangnguyen/workspace/data/spades/EcK12S-careful/assembly_graph.fastg");

        System.out.println(graph.getEdgeCount());
        BidirectedNode n128 = graph.getNode("128");
        Iterator<BidirectedEdge> ite = n128.getEdgeIterator();
        while(ite.hasNext()){
        	System.out.println(ite.next());
        }
        
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