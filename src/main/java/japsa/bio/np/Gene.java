package japsa.bio.np;

class Gene implements Comparable{
	String name;
	int st,end;
	int length;
	int cnt;  //number of times in genome;
	Gene(String[] str){
		this(str[0], Integer.parseInt(str[1]), Integer.parseInt(str[2]));
	}
	Gene(String name, int st, int end){
		this.st  = st;
		this.end = end;
		this.name = name;
		this.length = st - end;
	}
	
	@Override
	public boolean equals(Object obj){
		if(obj instanceof Gene){
			return ((Gene)obj).st == st &&  ((Gene)obj).end == end && ((Gene)obj).name == name;
		}else  return false;
	}

	@Override
	public int compareTo(Object obj) {
		int st1 = ((Gene)obj).st;
		if(st1 <st) return -1;
		else if(st1==st) return 0;
		else return 1;
	}
	
	
}