package japsadev.bio.hts.scaffold;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class StringHelper {

	/* 
	 * @param String, Graph: name of the FASTG file and the Graph to build
	 * @return int: shortest contig length in the FASTG file, used for k-mer guess
	 */
	public static int buildGraphFromFastg(String graphFile, Graph g) throws IOException{
		int shortestLen = Integer.MAX_VALUE;
		SequenceReader reader = new FastaReader(graphFile);
		Sequence seq;

		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length()<shortestLen)
				shortestLen=seq.length();
			
			String[] adjList = seq.getName().split(":");
			String name = adjList[0];
			boolean dir1=name.contains("'")?false:true;
			
			name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim();
			
			Vertex current=new Vertex(name);
			if(g.getVertex(current.getLabel())!=null)
				current=g.getVertex(current.getLabel());
				
			g.addVertex(current, false);
			
			if(dir1){
				seq.setName(name);
				current.setSequence(seq);
				//System.out.println(current);
			}
			if (adjList.length > 1){
				String[] nbList = adjList[1].split(",");
				for(int i=0; i < nbList.length; i++){
					// create list of bridges here (distance=-kmer overlapped)
					String neighbor = nbList[i];
					boolean dir2=neighbor.contains("'")?false:true;
					neighbor=neighbor.replaceAll("[^a-zA-Z0-9_.]", "").trim();

					Vertex nbVertex=new Vertex(neighbor);
					if(g.getVertex(nbVertex.getLabel())!=null)
						nbVertex=g.getVertex(nbVertex.getLabel());
	
					g.addVertex(nbVertex, false);
					
					g.addEdge(current, nbVertex, dir1, dir2);
				}
			}
			
		}
		reader.close();
		
		return shortestLen;
	}
	/*
	 * @param Path, String: path to be built from a string representing paths from File contigs.path 
	 * in SPAdes output directory.
	 * 
	 */
	public static void addPathFromSPAdes(Path p, String paths){
		Graph graph = p.graph;
		paths=paths.replace(";", ""); //optimized it!
		String[] comps = paths.split(",");
		for(int i=0; i<comps.length; i++){
			String 	cur = comps[i];
			boolean curDir = cur.contains("+")?true:false;
			Vertex curComp = graph.getVertex(cur.substring(0,cur.length()-1));
			if(curComp == null){
				System.out.println("Could not find Vertex "+ cur.substring(0,cur.length()-1));
				break;
			}		
			if(!p.nodes.isEmpty()){
				Node lastNode = p.nodes.get(p.nodes.size()-1);
				Edge curEdge=new Edge(lastNode.getVertex(), curComp, lastNode.getDirection(), curDir);
				if(!graph.containsEdge(curEdge)){
					System.out.println(curEdge + " doesn't exist in the graph!");
					break;
				}
					
			}
			p.addNode(curComp, curDir);
		}
	}
	
	public static void buildGraphFromDot(String graphFile, Graph graph) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(graphFile));
	    
		for(String line; (line = br.readLine()) != null; ) {
	        if(line.contains("edge")){
	        	String[] toks = line.split("\\s");
	        	int d = getDistanceFromDotFileBlock(toks[1]);
	        	if(d < 0)
	        		Graph.setKmerSize(-d); //actually (k-1)mer in ABySS
	        }
	        if(line.contains("->")){
	        	String[] toks = line.split("\\s");
	        	if(toks.length < 3)
	        		continue; //smt wrong actually!
	        	//get rid of quote characters
	        	String 	source = toks[0].replaceAll("\"", ""),
	        			dest = toks[2].replaceAll("\"", "");
	        	boolean sourceDir = source.contains("+")?true:false,
	        			destDir = dest.contains("+")?true:false;
	        	String 	sourceID = source.substring(0, source.length()-1),
	        			destID = dest.substring(0, dest.length()-1);
	        	Vertex 	sourceVertex = graph.getVertex(sourceID),
	        			destVertex = graph.getVertex(destID);
	        	
	        	if(sourceVertex==null){
					sourceVertex = new Vertex(sourceID);
					graph.addVertex(sourceVertex, false);
	        	}else if(destVertex==null){
					destVertex = new Vertex(destID);
					graph.addVertex(destVertex, false);
	        	}
	        	
	        	if(toks.length > 3) //distance available
	        		graph.addEdge(sourceVertex, destVertex, sourceDir, destDir, getDistanceFromDotFileBlock(toks[3]));
	        	else
	        		graph.addEdge(sourceVertex, destVertex, sourceDir, destDir);
	        	
	        }	
	    }

		br.close();
	}
	/*
	 * @param String [d=-%d]
	 * @return overlap distance, must be negative
	 */
	static int getDistanceFromDotFileBlock(String block){
		int d = 0;
    	String pattern = "^\\[d=([-]?[0-9]*)\\]$";
    	Pattern r = Pattern.compile(pattern);
    	Matcher m = r.matcher(block);
    	if(m.find()){
    		//System.out.println("Pattern matched: " + m.group(1));
    		d=Integer.parseInt(m.group(1));
    		//do smt
    	} else{
    		System.err.println("Not a legal block for distance!");
    	}
		return d;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String blk="[d=-1234]";
		System.out.println(getDistanceFromDotFileBlock(blk));
	}

}
