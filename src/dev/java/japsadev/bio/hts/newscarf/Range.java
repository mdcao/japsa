package japsadev.bio.hts.newscarf;

public class Range implements Comparable<Range>{
	int left, right;
	Range(){
		left=right=0;
	}
	Range(int left, int right){
		this.left=left;
		this.right=right;
	}
	
	public int getLeft(){
		return left;
	}
	public int getRight(){
		return right;
	}
	public void setLeft(int left){
		this.left=left;
	}
	public void setRight(int right){
		this.right=right;
	}
	
	public boolean isHomo(Range other){
		int order=compareTo(other);
		if(order==0) return true;
		
		boolean retval=false;
		Range 	ref=(order<0?this:other), 
				qry=(order<0?other:this);
		
		if(ref.right-qry.left > BidirectedGraph.ER_UPPERBOUND*BidirectedGraph.getKmerSize())
			retval=true;
		else 
			retval=(qry.right<=ref.right);
		
		return retval;
	}
	
    public void verifyHomo(Range other){
        if(!isHomo(other)){
            throw new IllegalStateException("Other range shouldn't belongs to the same group!");
        }
    }
    
	@Override
	public int compareTo(Range o) {
		// TODO Auto-generated method stub
		return left-o.left;
	}
	
	public String toString(){
		return new String(left+" -> "+right);
	}
    @Override
    public int hashCode() {
    	int retval=3;
    	retval=37*retval+left;
    	retval=37*retval+right;
        return retval;
    }

    @Override
    public boolean equals(Object obj) {
       if (!(obj instanceof Range))
            return false;
        if (obj == this)
            return true;

        Range rhs = (Range) obj;
        return left==rhs.getLeft()&&right==rhs.getRight();
    }
}
