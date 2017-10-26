package japsa.util;

public class Range implements Comparable<Range>{
   
	int left, right;
	Range(){
		left=right=0;
	}
	public Range(int left, int right){
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
	public void check(){
		assert right>left:"Invalid range: "+left + " to "+right;
	}
	public int getDistance(){
		return right-left;
	}
	public boolean merge(Range r){
		if(	r.left > right	|| r.right < left)
			return false;
		else{
			left=Math.min(left, r.left);
			right=Math.max(right, r.right);
			return true;
		}
	}
	@Override
	public int compareTo(Range o) {
		// TODO Auto-generated method stub
		return left-o.left;
	}
	
	public String toString(){
		return new String(left+"->"+right);
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
