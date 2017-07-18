package japsadev.bio.hts.newscarf;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Sequence;

public class BidirectedEdge extends AbstractEdge{
	protected boolean dir0, dir1;//true: outward, false: inward
	private int length=-BidirectedGraph.getKmerSize();//length of the edge (distance between tips of seqs represented by 2 nodes)
	private BidirectedPath path; //path that represents this edge. Should be a list if metagenome/multiploidy
	
	
    private static final Logger LOG = LoggerFactory.getLogger(BidirectedEdge.class);

	protected BidirectedEdge(String id, AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1) {
		// id fuck off!!! we'll make one for ourselves
		this(src,dst,dir0,dir1);
	}
	protected BidirectedEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1) {
		//this(createID(src,dst,dir0,dir1), src, dst);

		super(createID(src,dst,dir0,dir1),src,dst,false);
		this.dir0=dir0;
		this.dir1=dir1;
		
		path = new BidirectedPath();
		path.setRoot(src);
		path.add(this);
	}
	
	/* param id must have the form %s[o/i]%s[o/i], e.g. [1+2-]o[3]i
	 * the constructor will translate the id to the direction property 
	 * of the bidirected edge
	 */
//	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dest){
//		super(id, source, dest, false);
//    	String pattern = "^\\[([0-9\\+\\-]*)\\]([oi])\\[([0-9\\+\\-]*)\\]([oi])$";
//        // Create a Pattern object
//        Pattern r = Pattern.compile(pattern);
//        // Now create matcher object.
//        Matcher m = r.matcher(id);
//        String 	leftID, rightID, 
//        		leftDir, rightDir;
//        if(m.find()){
//        	leftID=m.group(1);
//        	leftDir=m.group(2);
//        	rightID=m.group(3);
//        	rightDir=m.group(4);
//        	if(source.getId().equals(leftID)){
//        		dir0=(leftDir.equals("o")?true:false);
//        		dir1=(rightDir.equals("o")?true:false);
//        	}else if(source.getId().equals(rightID)){
//        		dir0=(rightDir.equals("o")?true:false);
//        		dir1=(leftDir.equals("o")?true:false);
//        	}else{
//            	System.err.println("ID does not match");
//            	System.exit(1);
//            }
//        } else{
//        	System.err.println("Illegal ID for a bidirected edge (id must have the form node_id[o/i]node_id[o/i])");
//        	System.exit(1);
//        }
//
//	}
	
	/*
	 * To adapt Graph class from GraphStream library (called by newInstance())
	 */
	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dest){
		super(id, source, dest, false);
		
		assert (source.getGraph() == dest.getGraph()):"Nodes come from different graph " + source.getGraph().getId() + " and " + dest.getGraph().getId();
		BidirectedGraph g = (BidirectedGraph) source.getGraph();	
		path = new BidirectedPath(g,id);
		BidirectedNode 	n0 = (BidirectedNode) path.getRoot(),
						n1 = (BidirectedNode) path.peekNode();
		
		BidirectedEdge 	firstEdge = (BidirectedEdge) path.getEdgePath().get(0),
						lastEdge = (BidirectedEdge) path.peekEdge();
		if(n0.getId().equals(source.getId())){
			dir0 = firstEdge.getDir(n0);
			dir1 = lastEdge.getDir(n1);
		}else if(n0.getId().equals(dest.getId())){
			dir1 = firstEdge.getDir(n0);
			dir0 = lastEdge.getDir(n1);
			path = path.getReversedComplemented();
		}else{
			LOG.error("Path {} conflicts with src={} dst={}!", id, source.getId(), dest.getId());
			System.exit(1);
		}
	}
		
//	public static String createID(AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1){
//		String 	srcDes = "["+source.getId()+"]"+(dir0 ? "o":"i"),
//				dstDes = "["+dst.getId()+"]"+(dir1 ? "o":"i");
//		if(srcDes.compareTo(dstDes)<0)
//			return String.format("%s%s", srcDes, dstDes);
//		else
//			return String.format("%s%s", dstDes, srcDes);
//	}
	
	
	//should have smt like this
	public static BidirectedEdge makeEdge (AbstractNode src, AbstractNode dst, BidirectedPath path) {
		return new BidirectedEdge(path.getId(), src, dst);
	}

//	public BidirectedEdge getReversedComplemented(){
//		BidirectedEdge retval = makeEdge(getNode1(),getNode0(),path.getReversedComplemented());
//		
//		return retval;
//	}
	
	public void setPath(BidirectedPath path){
		this.path = path;
		//here set the new length due to that path
		Node 	curNode = path.getRoot();
		for(Edge e:path.getEdgePath()){
			curNode = e.getOpposite(curNode);
			if(curNode == path.peekNode())
				break;
			else{
				length+=((Sequence)curNode.getAttribute("seq")).length()-BidirectedGraph.getKmerSize();
			}
		}
	}
	public BidirectedPath getPath(){
		return path;
	}
	
	@Override
	public String toString() {
		return String.format("%s:%s-%s-%s-%s", getId(), source, (dir0?">":"<"), (dir1?"<":">"), target);
	}
	
	public void setDir0(boolean dir){
		this.dir0=dir;
	}
	public void setDir1(boolean dir){
		this.dir1=dir;
	}
	public boolean getDir0(){
		return dir0;
	}
	public boolean getDir1(){
		return dir1;
	}
	//if kmer!=127 we need to set initial lengths of edges again (easy+tedious way)
	public void changeKmerSize(int kmer){
		length=-kmer;
	}
	public int getLength(){
		return length;
	}
	
	//TODO: include the case of tandem repeats
	public boolean getDir(AbstractNode node){
		assert node==getSourceNode()||node==getTargetNode():"Node " + node.getId() + " does not belong to this edge src=" + getSourceNode().getId() + " dst=" + getTargetNode().getId();
		return node==getSourceNode()?getDir0():getDir1();
	}
	

	public static String createID(AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1){
		String 	srcDes = source.getId(),
				dstDes = dst.getId();
		if(srcDes.compareTo(dstDes)<0)
			return String.format("%s%s,%s%s", srcDes, (dir0?"+":"-"), dstDes, (dir1?"-":"+"));
		else
			return String.format("%s%s,%s%s", dstDes, (dir1?"+":"-"), srcDes, (dir0?"-":"+"));
	}
}
