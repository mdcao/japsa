package japsadev.bio.hts.newscarf;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedEdge extends AbstractEdge{
	protected boolean dir0, dir1;//true: outward, false: inward
	
	protected BidirectedEdge(String id, AbstractNode source, AbstractNode dst, boolean dir0, boolean dir1) {
		super(id, source, dst, false);
		// TODO Auto-generated constructor stub
		this.dir0=dir0;
		this.dir1=dir1;
	}
	
	@Override
	public String toString() {
		return String.format("%s[%s%s%s%s]", getId(), source, dir0 ? "+"
				: "-", dir1 ? "+" : "-", target);
	}
	
	public boolean getDir0(){
		return dir0;
	}
	public boolean getDir1(){
		return dir1;
	}
	
}
