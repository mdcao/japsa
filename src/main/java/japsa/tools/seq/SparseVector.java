package japsa.tools.seq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author Lachlan Coin
 *
 */


public class SparseVector{
	static Double zero = 0.0;
	
	SortedMap<Integer, Number> m = new TreeMap<Integer,Number>();
	private double valsum=0;
	//double valsumByKey=0;
	
	public void addZero(int pos){
		this.m.put(pos, zero);
	}
	
	public void addToEntry(Integer position, double i) {
		Number val = m.get(position);
		valsum+=i;
		//valsumByKey+=position.doubleValue()* (double) i;
		m.put(position, val==null ? i : val.doubleValue() + i);
	}
	public String toString(){
		return m.toString();
	}

	public List<Integer> keys() {
		List<Integer> l= new ArrayList<Integer>(this.m.keySet());
		Collections.sort(l);
		return l;
	}
	public List<Integer> keys(double thresh) {
		List<Integer> l= new ArrayList<Integer>();
		for(Iterator<Map.Entry<Integer, Number>> it = m.entrySet().iterator(); it.hasNext();){
			Map.Entry<Integer, Number> nxt = it.next();
			if(nxt.getValue().doubleValue()>=thresh){
				l.add(nxt.getKey());
			}
		}
		Collections.sort(l);
		return l;
	}
	public Iterator<Integer> keyIt(){
		return m.keySet().iterator();
	}

	public Number get(Integer val) {
		Number res =  this.m.get(val);
		if(res==null) return zero;
		else return res;
	}
	
	
	

	public void clear() {
		m.clear();
		this.valsum=0;
		
	}
	
	public Number getDepth(Integer i) {
		Number val = m.get(i);
		return val==null ? zero: val;
	}

	public double valsum() {
	return valsum;
	}

	public void update(Integer spec, double q) {
		Number qual = this.m.get(spec);
		if(qual==null  || q > qual.doubleValue()){
			this.m.put(spec,q);
		}
	}
	public SparseVector clone(){
		SparseVector sv =  new SparseVector();
		for(Iterator<Integer> it = this.m.keySet().iterator(); it.hasNext();){
			Integer key = it.next();
			this.m.put(key, this.m.get(key));
		}
		return sv;
	}

	public Set<Integer> keySet() {
		return m.keySet();
	}

	public int size() {
		return this.keySet().size();
	}

	public Map<Integer, Number> getMax() {
		Iterator<Map.Entry<Integer, Number>> it = this.m.entrySet().iterator();
		Entry<Integer, Number>  v = it.next();
		while( it.hasNext()){
			Entry<Integer, Number> nxt = it.next();
			if(nxt.getValue().doubleValue()>v.getValue().doubleValue()) v = nxt;
		}
		Map<Integer, Number> m1 = new HashMap<Integer, Number>();
		it = this.m.entrySet().iterator();
		while( it.hasNext()){
			Entry<Integer, Number> nxt = it.next();
			if(nxt.getValue().doubleValue()>=v.getValue().doubleValue()){
				m1.put(nxt.getKey(), nxt.getValue());
			}
			
		}
		return m1;
		
	}

	

	

	

	
}