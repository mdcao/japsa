package japsadev.bio.hts.newscarf;

import com.sun.org.apache.bcel.internal.generic.INSTANCEOF;

public class Range implements Comparable<Range>{
	int start, end;
	Range(){
		start=end=0;
	}
	Range(int from, int to){
		this.start=from;
		this.end=to;
	}
	
	public String toString(){
		return new String(start+" -> "+end);
	}
	@Override
	public int compareTo(Range o) {
		// TODO Auto-generated method stub
		if(end-o.start < 1.2*BidirectedGraph.getKmerSize())
			return -1;
		else if(start-o.end>-1.2*BidirectedGraph.getKmerSize())
			return 1;
		else
			return 0;
	}
	@Override
	public boolean equals(Object obj){
	    if (obj == null) {
	        return false;
	    }
	    if (!Range.class.isAssignableFrom(obj.getClass())) {
	        return false;
	    }
	    
		return compareTo((Range)obj)==0;
		
	}
	
	@Override
	public int hashCode() {
	    int hash = 3;
//	    hash = 53 * hash + (this.name != null ? this.name.hashCode() : 0);
//	    hash = 53 * hash + this.age;
	    return hash;
	}
}
