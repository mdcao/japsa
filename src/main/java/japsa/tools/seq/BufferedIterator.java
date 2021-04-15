package japsa.tools.seq;

import java.util.Iterator;
import java.util.Stack;

public class BufferedIterator<T> implements Iterator<T> {
	public BufferedIterator(Iterator<T> it, int buff){
		this.it1 = it;
		for(int i=0; i<buff && it.hasNext(); i++){
			fastq.push(it.next());
		}
	}
	public int stackSize(){
		return fastq.size();
	}
	Iterator<T> it1;
	
   Stack<T> fastq = new Stack<T>();
	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
		return fastq.size()>0 || it1.hasNext();
	}

	@Override
	public T next() {
		if(fastq.size()>0) return fastq.pop();
		else return it1.next();
	}

}
