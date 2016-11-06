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
	
	/* param id must have the form %d[o/i]%d[o/i]
	 * the constructor will translate the id to the direction property 
	 * of the bidirected edge
	 */
	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dst){
		super(id, source, dst, false);
        String pattern = "\\b(\\d+)([oi])(\\d+)([oi])\\b";
        // Create a Pattern object
        Pattern r = Pattern.compile(pattern);
        // Now create matcher object.
        Matcher m = r.matcher(id);
        String 	
        		//srcID, dstID, 
        		srcDir, dstDir;
        if(m.find()){
        	srcDir=m.group(2);
        	dstDir=m.group(4);
        	dir0=(srcDir=="o"?true:false);
        	dir1=(dstDir=="o"?true:false);
        } else{
        	System.err.println("Illegal ID for a bidirected edge (id must have the form %d[o/i]%d[o/i])");
        	System.exit(1);
        }

	}
		
	public static String createID(AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1){
		String 	srcDes = source.getId()+(dir0 ? "o":"i"),
				dstDes = dst.getId()+(dir1 ? "o":"i");
		if(srcDes.compareTo(dstDes)<0)
			return String.format("%s%s", srcDes, dstDes);
		else
			return String.format("%s%s", dstDes, srcDes);
	}
	
	
	
	@Override
	public String toString() {
		return String.format("%s[%s%s%s%s]", getId(), source, (dir0?">":"<"), (dir1?"<":">"), target);
	}
	
	public boolean getDir0(){
		return dir0;
	}
	public boolean getDir1(){
		return dir1;
	}
	
}
