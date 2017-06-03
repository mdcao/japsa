package japsadev.bio.hts.newscarf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class HybridAssembler {
    private static final Logger LOG = LoggerFactory.getLogger(HybridAssembler.class);
	
//	final BidirectedGraph origGraph;
	public BidirectedGraph simGraph; //original and simplified graph should be separated, no???
	
	public HybridAssembler(){
//		origGraph=new BidirectedGraph("batch");
		simGraph=new BidirectedGraph("real");
	}
	
	
	public HybridAssembler(String graphFile) throws IOException{
		this();
//		origGraph.loadFromFile(graphFile);
		simGraph.loadFromFile(graphFile);
	}
	
	
	public void assembly(String bamFile, int qual) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		SamReader reader;
		if ("-".equals(bamFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		//ReadFilling readFilling = null;
		ArrayList<Alignment> samList =  new ArrayList<Alignment>();;// alignment record of the same read;	
		BidirectedPath p = new BidirectedPath();
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			if (rec.getMappingQuality() < qual)
				continue;
			
			String refID = rec.getReferenceName().split("_")[1];
			Alignment myRec = new Alignment(rec, simGraph.getNode(refID)); //FIXME: optimize

			//////////////////////////////////////////////////////////////////
			// make list of alignments of the same (Nanopore) read. 

			//not the first occurrance				
			if (!readID.equals("") && !readID.equals(myRec.readID)) {		
				//Collections.sort(samList);
				//p=origGraph.pathFinding(samList);
				p=simGraph.pathFinding(samList); // the graph MUST be the same as from new Alignment(...)

				if(p!=null)
					System.out.println("Final path found: " + p.getId());
				reduce(p);
				samList = new ArrayList<Alignment>();
				//readID = myRec.readID;	
			}	
			readID = myRec.readID;
			samList.add(myRec); // FIXME: (optimize) insert sort here

		}// while
		iter.close();

		//outOS.close();
		reader.close();		
	
	}
    /*
     * Read paths from contigs.path and reduce the graph
     */
    public void reduceFromSPAdesPaths(String paths) throws IOException{

		BufferedReader pathReader = new BufferedReader(new FileReader(paths));
		
		String s;
		//Read contigs from contigs.paths and refer themselves to contigs.fasta
		boolean flag=false;
		while((s=pathReader.readLine()) != null){
			if(s.contains("NODE")){
				flag=s.contains("'")?false:true;
				continue;
			}else if(flag){
				//BidirectedPath path=new BidirectedPath(origGraph, s);
				BidirectedPath path=new BidirectedPath(simGraph, s);

		    	reduce(path);

//				if(comp!=null){
//					System.out.println("Reverting node: " + comp.getId());
//					revert(comp);
//			        System.out.println("After revert => Node: " + getNodeCount() + " Edge: " + getEdgeCount());
//
//				}
			}	
				

		}
		pathReader.close();
    }
	
    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param p Path to simplify the graph (from origGraph)
     * @param target Subjected graph for the simplification
     */
    private void reduce(BidirectedPath p){
    	//do nothing if the path has only one node
    	if(p==null || p.getEdgeCount()<1)
    		return;
    	
    	//loop over the edges of path (like spelling())
    	BidirectedNode 	markerNode = null,
    			curNodeFromSimGraph = (BidirectedNode) p.getRoot();
	
    	BidirectedPath curPath= null;
    	boolean markerDir=true, curDir;
    	
    	if(BidirectedGraph.isUnique(curNodeFromSimGraph)){
    		markerNode=curNodeFromSimGraph;
    		markerDir=((BidirectedEdge) p.getEdgePath().get(0)).getDir(markerNode);
    		curPath = new BidirectedPath();
    		curPath.setRoot(curNodeFromSimGraph);
    	}
    	

    	//search for an unique node as the marker. 
    	ArrayList<BidirectedEdge> 	tobeRemoved = new ArrayList<BidirectedEdge>(),
    								tobeAdded = new ArrayList<BidirectedEdge>();
    	for(Edge e:p.getEdgePath()){
    			
    		curNodeFromSimGraph=e.getOpposite(curNodeFromSimGraph);
    		   		
//    		curNodeFromSimGraph = simGraph.getNode(curNodeFromOrigGraph.getId()); //change back to Node belong to simGraph (instead of origGraph)
    		curDir=((BidirectedEdge) e).getDir(curNodeFromSimGraph);
    		
    		if(BidirectedGraph.isUnique(curNodeFromSimGraph)){
        		
				if(markerNode!=null){
					//this is when we have 1 jumping path (both ends are markers)
					curPath.add(e);	
//					LOG.info("Processing path {} with marker {}:{}:{} and curNode {}:{}:{}", curPath.getId(), markerNode.getId(), markerDir?"out":"in", markerNode.getGraph().getId(), curNodeFromSimGraph.getId(), curDir?"out":"in", curNodeFromSimGraph.getGraph().getId());
					//create an edge connect markerNode to curNode with curPath
					//Edge reducedEdge = simGraph.addEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);
					BidirectedEdge reducedEdge = new BidirectedEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);

//					if(reducedEdge!=null)
//						reducedEdge.addAttribute("path", new BidirectedPath(curPath));
				
					tobeAdded.add(reducedEdge);
					
					//loop over curPath to find out edges needed to be removed
					Node  	n0 = curPath.getRoot(),
							n1 = null;
					for(Edge ep:curPath.getEdgePath()){
						n1 = ep.getOpposite(n0);
						if(!BidirectedGraph.isUnique(n0) == BidirectedGraph.isUnique(n1)){
			    			tobeRemoved.add((BidirectedEdge)ep);
						}

//			    		if(!BidirectedGraph.isUnique(n1)){			    			
//			    			n1.setAttribute("cov", n1.getNumber("cov")-markerNode.getNumber("cov"));   		
//			    			LOG.info("...coverage of " + n1.getAttribute("name") + " now is " + n1.getNumber("cov"));
//			    		}
						
			    		n0=n1;
					}

				}
				
				
				markerNode=curNodeFromSimGraph;
        		markerDir=!curDir; //in-out, out-in
				curPath= new BidirectedPath();
				curPath.setRoot(curNodeFromSimGraph);
    		}
    		else{
    			if(markerNode!=null){
    				curPath.add(e);
    			}
    		}
    		
		}
    	
    	//remove appropriate edges
    	for(BidirectedEdge e:tobeRemoved){
    		LOG.info("REMOVING EDGE " + e.getId() + " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
    		LOG.info("before: \n\t" + simGraph.printEdgesOfNode(e.getNode0()) + "\n\t" + simGraph.printEdgesOfNode(e.getNode1()));
    		simGraph.removeEdge(e.getId());
    		LOG.info("after: \n\t" + simGraph.printEdgesOfNode(e.getNode0()) + "\n\t" + simGraph.printEdgesOfNode(e.getNode1()));
    	}
    	
    	//add appropriate edges
    	for(BidirectedEdge e:tobeAdded){
    		LOG.info("ADDING EDGE " + e.getId()+ " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
    		LOG.info("before: \n\t" + simGraph.printEdgesOfNode(e.getNode0()) + "\n\t" + simGraph.printEdgesOfNode(e.getNode1()));
    		
    		Edge reducedEdge = simGraph.addEdge(e.getSourceNode(),e.getTargetNode(),e.getDir0(),e.getDir1());
			if(reducedEdge!=null){
//				reducedEdge.addAttribute("path", new BidirectedPath(e.getAttribute("path")));
				reducedEdge.setAttribute("ui.style", "text-offset: -10;"); 
				reducedEdge.setAttribute("ui.class", "marked");
			}
    		LOG.info("after: \n\t" + simGraph.printEdgesOfNode(e.getNode0()) + "\n\t" + simGraph.printEdgesOfNode(e.getNode1()));

    	}

//		promptEnterKey();
    }
    
    @SuppressWarnings("resource")
	public static void promptEnterKey(){
    	   System.out.println("Press \"ENTER\" to continue...");
    	   Scanner scanner = new Scanner(System.in);
    	   scanner.nextLine();
    	}
    
    protected void sleep() {
        try { Thread.sleep(1000); } catch (Exception e) {}
    }

	public static void main(String[] argv) throws IOException{
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
		//TODO: need to make this easier
		hbAss.assembly(GraphExplore.spadesFolder+"bwa/EcK12S-careful.sam", 30);
	}
	
}
