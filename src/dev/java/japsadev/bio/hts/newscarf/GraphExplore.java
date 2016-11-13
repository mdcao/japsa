package japsadev.bio.hts.newscarf;

import java.io.IOException;
import java.util.Iterator;
import org.graphstream.graph.*;
public class GraphExplore {
    public static void main(String args[]) {
    	new GraphExplore();

    }

    public GraphExplore(){
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
        
        /*
         * Testing reduce function
         */
        try {
			graph.readPathsFromSpades("/home/hoangnguyen/workspace/data/spades/EcK12S-careful/contigs.paths");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        /*
         * Testing BidirectedEdge id pattern
         */
//    	String pattern = "^\\[([0-9oi]*)\\]([oi])\\[([0-9oi]*)\\]([oi])$";
//        // Create a Pattern object
//        Pattern r = Pattern.compile(pattern);
//        // Now create matcher object.
//        String id="[3i56i324o67i]i[4o8i]o";
//        Matcher m = r.matcher(id);
//         	
//        if(m.find()){
//        	System.out.println(m.group(1)+"-"+m.group(2)+"-"+m.group(3)+"-"+m.group(4));
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
        "}";
    	
}