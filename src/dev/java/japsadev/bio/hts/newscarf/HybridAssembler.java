package japsadev.bio.hts.newscarf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import org.graphstream.graph.Edge;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.util.Logging;

public class HybridAssembler {
	final BidirectedGraph origGraph;
	BidirectedGraph simGraph; //original and simplified graph should be separated, no???
	
	public HybridAssembler(){
		origGraph=new BidirectedGraph();
		simGraph=new BidirectedGraph();
	}
	
	
	public HybridAssembler(String graphFile) throws IOException{
		this();
		origGraph.loadFromFile(graphFile);
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
			Alignment myRec = new Alignment(rec, origGraph.getNode(refID)); //FIXME: optimize

			//////////////////////////////////////////////////////////////////
			// make list of alignments of the same (Nanopore) read. 

			//not the first occurrance				
			if (!readID.equals("") && !readID.equals(myRec.readID)) {		
				//Collections.sort(samList);
				p=origGraph.pathFinding(samList);
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
				BidirectedPath path=new BidirectedPath(origGraph, s);

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
     * @param p Path to simplify the graph (from base graph)
     * @param target Subjected graph for the simplification
     */
    private void reduce(BidirectedPath p){
    	//do nothing if the path has only one node
    	if(p==null || p.getEdgeCount()<1)
    		return;
    	
    	//loop over the edges of path (like spelling())
    	BidirectedNode 	markerNode = null,
    					curNode = (BidirectedNode) p.getRoot();

    	BidirectedPath curPath= null;
    	if(BidirectedGraph.isUnique(curNode)){
    		markerNode=curNode;
    		curPath = new BidirectedPath();
    		curPath.setRoot(curNode);
    	}
    	
    	boolean markerDir=true, curDir=true;
    	//search for an unique node as the marker. 
    	ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
    	for(Edge e:p.getEdgePath()){
    			
    		curNode=e.getOpposite(curNode);
    		curDir=((BidirectedEdge) e).getDir(curNode);
    		if(BidirectedGraph.isUnique(curNode)){
        		
				if(markerNode!=null){
					curPath.add(e);	
					//create an edge connect markerNode to curNode with curPath
					Edge reducedEdge = simGraph.addEdge(markerNode, curNode, markerDir, curDir);
					if(reducedEdge!=null){
						reducedEdge.addAttribute("path", new BidirectedPath(curPath));
						reducedEdge.setAttribute("ui.style", "text-offset: -10;"); 
						reducedEdge.setAttribute("ui.class", "marked");
					}
					//and do necessary things
					
				}
				markerNode=curNode;
        		markerDir=curDir;
				curPath= new BidirectedPath();
				curPath.setRoot(curNode);
    		}
    		else{
    			if(markerNode!=null){
    				curPath.add(e);
    				//curNode.setAttribute("cov", (double)curNode.getAttribute("cov")-(double)markerNode.getAttribute("cov"));
    			}
    		}
    		
    		if(BidirectedGraph.isUnique(curNode) || BidirectedGraph.isUnique(e.getOpposite(curNode)))
    			tobeRemoved.add(e);
		}
    	
    	//remove appropriate edges
    	for(Edge e:tobeRemoved)
    		simGraph.removeEdge(e);
    }
    
	public static void main(String[] argv) throws IOException{
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
		//TODO: need to make this easier
		hbAss.assembly(GraphExplore.spadesFolder+"bwa/EcK12S-careful.sam", 0);
	}
	
}
