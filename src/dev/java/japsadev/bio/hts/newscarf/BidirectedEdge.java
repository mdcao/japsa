package japsadev.bio.hts.newscarf;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedEdge extends AbstractEdge{
	protected boolean dir0, dir1;//true: outward, false: inward
	
	protected BidirectedEdge(String id, AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1) {
		// id fuck off!!! we'll make one for ourselves
		this(src,dst,dir0,dir1);
	}
	protected BidirectedEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1) {
		//this(createID(src,dst,dir0,dir1), src, dst);

		super(createID(src,dst,dir0,dir1),src,dst,false);
		this.dir0=dir0;
		this.dir1=dir1;
	}
	
	/* param id must have the form %s[o/i]%s[o/i], e.g. [1+2-]o[3]i
	 * the constructor will translate the id to the direction property 
	 * of the bidirected edge
	 * TODO: read & do smt with the composite node id
	 */
	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dest){
		super(id, source, dest, false);
    	String pattern = "^\\[([0-9\\+\\-]*)\\]([oi])\\[([0-9\\+\\-]*)\\]([oi])$";
        // Create a Pattern object
        Pattern r = Pattern.compile(pattern);
        // Now create matcher object.
        Matcher m = r.matcher(id);
        String 	leftID, rightID, 
        		leftDir, rightDir;
        if(m.find()){
        	leftID=m.group(1);
        	leftDir=m.group(2);
        	rightID=m.group(3);
        	rightDir=m.group(4);
        	if(source.getId().equals(leftID)){
        		dir0=(leftDir.equals("o")?true:false);
        		dir1=(rightDir.equals("o")?true:false);
        	}else if(source.getId().equals(rightID)){
        		dir0=(rightDir.equals("o")?true:false);
        		dir1=(leftDir.equals("o")?true:false);
        	}else{
            	System.err.println("ID does not match");
            	System.exit(1);
            }
        } else{
        	System.err.println("Illegal ID for a bidirected edge (id must have the form node_id[o/i]node_id[o/i])");
        	System.exit(1);
        }

	}
		
	public static String createID(AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1){
		String 	srcDes = "["+source.getId()+"]"+(dir0 ? "o":"i"),
				dstDes = "["+dst.getId()+"]"+(dir1 ? "o":"i");
		if(srcDes.compareTo(dstDes)<0)
			return String.format("%s%s", srcDes, dstDes);
		else
			return String.format("%s%s", dstDes, srcDes);
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
	//TODO: include the case of tandem repeats
	public boolean getDir(AbstractNode node){
		assert node==getSourceNode()||node==getTargetNode():"Node does not belong to this edge!";
		return node==getSourceNode()?getDir0():getDir1();
	}
}
