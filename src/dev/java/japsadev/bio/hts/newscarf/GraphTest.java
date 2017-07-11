package japsadev.bio.hts.newscarf;

import org.graphstream.graph.Graph;
import org.graphstream.graph.implementations.AdjacencyListGraph;

public class GraphTest {
    public static void main(String[] args) throws Exception {
        System.setProperty("org.graphstream.ui.renderer",
                "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
        Graph g = new AdjacencyListGraph("g");
        g.addNode("A").addAttribute("xyz", new double[] { 0, 0 });
        g.addNode("B").addAttribute("xyz", new double[] { 10, 10 });

        g.addEdge("AB", "A", "B", false)
                .addAttribute(
                        "ui.points",
                        (Object) new double[] { 0, 0, 0, 0, 5, 0, 5, 10, 0, 10,
                                10, 0 });

        g.addAttribute("ui.stylesheet", "edge {shape: polyline; }"); // or shape: cubic-curve

        g.display(false);
    }
}
