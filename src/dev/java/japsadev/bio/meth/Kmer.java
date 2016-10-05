package japsadev.bio.meth;

public class Kmer implements Comparable<Kmer>{
	boolean isMeth=false;
	int coordinate=-1;
	boolean strand = true; // +/- correspond to true/false
	String containedSeq = "";
	String desc = "", sequence = "";
	
	double mean=0.0, stdv=0.0; //might be a class
	public Kmer(){
	}
	public Kmer(boolean isMeth, boolean strand, String container, int coordinate){
		this.isMeth = isMeth;
		this.strand = strand;
		this.containedSeq = container;
		this.coordinate = coordinate;
	}

	public void setSignal(double mean, double stdv){
		this.mean = mean;
		this.stdv = stdv;
	}
	public double getMeanSignal(){
		return mean;
	}
	public double getStdvSignal(){
		return stdv;
	}
	
	public void setCoordinate(int coor){
		coordinate = coor;
	}
	public int getCoordinate(){
		return coordinate;
	}
	
	public void setDesc(String str){
		desc = new String(str);
	}
	public String getDesc(){
		return desc;
	}
	
	public void setSeq(String str){
		sequence = new String(str);
	}
	public String getSeq(){
		return sequence;
	}
	
	public String toString(){
		return new String(containedSeq + ":" + coordinate +":" + (strand?'+':'-'));
	}

	@Override
	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
		return false;
		// object must be Test at this point
		Kmer test = (Kmer)obj;
		return strand == test.strand && coordinate == test.coordinate &&
		(containedSeq == test.containedSeq || (containedSeq != null && containedSeq.equals(test.containedSeq)));
	}

	public int hashCode()
	{
			int hash = 7;
			hash = 31 * hash + (strand?1:0);
			hash = 31 * hash + coordinate;
			hash = 31 * hash + (null == containedSeq ? 0 : containedSeq.hashCode());
			return hash;
	}
	@Override
	public int compareTo(Kmer kmer) {
		// TODO Auto-generated method stub
		return coordinate - kmer.coordinate;
	}
}
