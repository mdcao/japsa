package japsadev.tools;
import java.io.IOException;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.ui.view.Viewer;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsadev.bio.hts.newscarf.Alignment;
import japsadev.bio.hts.newscarf.BidirectedGraph;
import japsadev.bio.hts.newscarf.HybridAssembler;

@Deployable(
		scriptName = "jsa.dev.newScarf", 
		scriptDesc = "New npscarf"
		)
public class NewScarfCmd extends CommandLine{
	public NewScarfCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("fastg", null, "Assembly graph fastg file",true);		
		addString("sam", null, "Sam file alignment of assembly graph to long reads",true);
		addInt("qual", 30, "Minimum quality of alignment to considered");
		addString("path", null, "SPAdes contigs path file");
		addStdHelp();
	}
		    	
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new NewScarfCmd ();
		args = cmdLine.stdParseLine(args);

		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		String fastgFile = cmdLine.getStringVal("fastg");
		String samFile = cmdLine.getStringVal("sam");
		String pathFile = cmdLine.getStringVal("path");
		
		String styleSheet =
			        "node {" +
			        "	fill-color: black; z-index: 0;" +
			        "}" +
			        "edge {" +
			        "	text-alignment: along;" +
			        "}" +
			        "node.marked {" +
			        "	fill-color: red;" +
			        "}" +
			        "edge.marked {" +
			        "	fill-color: red;" +
			        "}";
		System.setProperty("java.awt.headless", "false");
		HybridAssembler hbAss = new HybridAssembler(fastgFile);
		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
		//TODO: need to make this easier
		BidirectedGraph graph= hbAss.simGraph;
		
        //graph.addAttribute("ui.quality");
        //graph.addAttribute("ui.antialias");
        graph.addAttribute("ui.stylesheet", styleSheet);
        Viewer viewer = graph.display();
        // Let the layout work ...
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());

        
        for (Node node : graph) {
            node.addAttribute("ui.label", node.getId());
            node.setAttribute("ui.style", "text-offset: -10;"); 
            node.addAttribute("layout.weight", 10); 

            if(BidirectedGraph.isUnique(node))
            	node.setAttribute("ui.class", "marked");
        }


        try {
        	if(pathFile!=null)
        		hbAss.reduceFromSPAdesPaths(pathFile);
        	hbAss.assembly(samFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        

//        for (Edge edge: graph.getEdgeSet()){
//        	if(edge.hasAttribute("isReducedEdge"))
//        		edge.addAttribute("layout.weight", 10); 
//        }
        
        
        HybridAssembler.promptEnterKey();
        viewer.disableAutoLayout();
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
	}
}
