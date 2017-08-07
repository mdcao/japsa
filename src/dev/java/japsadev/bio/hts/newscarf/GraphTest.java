package japsadev.bio.hts.newscarf;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.Graph;
import org.graphstream.graph.implementations.AdjacencyListGraph;

import htsjdk.samtools.*;

public class GraphTest {
    public static void main(String[] args) throws Exception {
//        System.setProperty("org.graphstream.ui.renderer",
//                "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
//        Graph g = new AdjacencyListGraph("g");
//        g.addNode("A").addAttribute("xyz", new double[] { 0, 0 });
//        g.addNode("B").addAttribute("xyz", new double[] { 10, 10 });
//
//        g.addEdge("AB", "A", "B", false)
//                .addAttribute(
//                        "ui.points",
//                        (Object) new double[] { 0, 0, 0, 0, 5, 0, 5, 10, 0, 10,
//                                10, 0 });
//
//        g.addAttribute("ui.stylesheet", "edge {shape: polyline; }"); // or shape: cubic-curve
//
//        g.display(false);
    	
//    	String data="/home/hoangnguyen/workspace/data/spades/EcK12S-careful/assembly_graph.fastg"; //sony
//
//    	final SamFileValidator validator=new SamFileValidator(new PrintWriter(System.out),8000);
//        validator.setIgnoreWarnings(true);
//        validator.setVerbose(true,1000);
//        validator.setErrorsToIgnore(Collections.singletonList(SAMValidationError.Type.MISSING_READ_GROUP));
//        SamReaderFactory factory=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
//        SamReader samReader=factory.open(new File(data));
//        boolean ischeck=false;
//		SAMRecordIterator iter = samReader.iterator();
//		System.out.println(iter.next());
//        
//        samReader.close();
    	
    	Pattern versionPattern = Pattern.compile("^Version:\\s(\\d+\\.\\d+\\.\\d+).*");
    	String line="Version: 0.7.5a-r405";
		Matcher matcher =versionPattern.matcher(line);
		if (matcher.find()){
			System.out.println(line);
			System.out.println(matcher.group(1));
		}
    }
}
