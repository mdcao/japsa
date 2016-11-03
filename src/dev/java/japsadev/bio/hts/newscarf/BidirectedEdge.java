package japsadev.bio.hts.newscarf;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedEdge extends AbstractEdge{
	protected boolean dir0, dir1;//true: outward, false: inward
	
	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1) {
		// id fuck off!!!
		this(source,dst,dir0,dir1);
	}
	protected BidirectedEdge(AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1) {
		super(createID(source,dst,dir0,dir1), source, dst, false);
		// TODO Auto-generated constructor stub
		this.dir0=dir0;
		this.dir1=dir1;
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
