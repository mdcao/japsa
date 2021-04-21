package japsa.tools.seq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import japsa.bio.np.RealtimeSpeciesTyping;

public class SparseVectorCollection{
		// each entry in svs represents a multimapping read
		public final double[] abundance;
		int len;
		public SparseVectorCollection(int numSpecies){
			this.len = numSpecies;
			this.abundance = new double[numSpecies];
			double val = 1.0/(double) numSpecies;
			Arrays.fill(abundance, val);
		}
		List<SparseVector> svs = new ArrayList<SparseVector>();
	//	List<Integer> single = new ArrayList<Integer>();
		public void add(SparseVector all_species) {
			Set<Integer> keys = all_species.keySet();
			if(keys.size()==1 && false){
		//		single.add(keys.iterator().next());
			}else{
				svs.add(all_species);
			}
		}
		
		public double logLike(){
			double sc =0;
			for(int i=0; i<this.svs.size(); i++){
				SparseVector sv = svs.get(i); // this represents a read
				double pr = 0;
				for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
					Integer j = it.next();
					pr+= abundance[j] * sv.get(j).doubleValue();
				}
				sc+=Math.log(pr);
			}
			return sc;
		}
		
		
		//update abundance
		//one round of normalisation
		public void maximisation(double[] v, double pseudo){
			double[] abund = new double[len];
			Arrays.fill(abund, pseudo);
			if(svs.size()==0){// && single.size()==0){
			//	System.err.println("no matches");
				return;//throw new RuntimeException("nothing");
			}
			double tot = pseudo*len;//+svs.size()+single.size();
			tot+=svs.size();//+ single.size();//  add 
			double score =0;
			for(int i=0; i<this.svs.size(); i++){
				SparseVector sv = svs.get(i);
				
				double totj = 0;
				double scorej=0;
				//split the read count according to the quality scores
				for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
					Integer j = it.next();
					double sc = abundance[j] * sv.get(j).doubleValue();
					totj+=sc;
					scorej+=sc;
				}
				if(totj>0){
					for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
						Integer j = it.next();
						double sc = abundance[j] * sv.get(j).doubleValue()/totj;
						abund[j] += sc;
					}
					score+=Math.log(scorej);
				}else{
					tot = tot -1;
					System.err.println("warning totj is zero "+i+" "+sv.toString());
				}
				
			}
			double diff = 0;
			for(int j=0; j<abund.length; j++){
				double newv = abund[j]/tot;
				diff += Math.abs(abundance[j] - newv);
				this.abundance[j] =newv ;
			}
			
			v[0] = score; v[1] = diff;
			if(zerovs.size()>0){ // if we keeping anything zero
				Iterator<Integer> iter = zerovs.iterator();
				while( iter.hasNext()){
					Integer nxt = iter.next();
					this.abundance[nxt] = RealtimeSpeciesTyping.epsilon;
					
				}
				this.renormalise();
			}
			
		}
		
		public double setZero(int i, double v1){
			double v = abundance[i];
			abundance[i] =v1;
			this.renormalise();
			return v;
		}
		
		public void renormalise() {
			double tot = 0;
			for(int i =0; i<this.abundance.length; i++){
				tot+=abundance[i];
			}
			for(int i =0; i<this.abundance.length; i++){
				abundance[i] = abundance[i]/tot;
			}
			this.check();
			
		}

		public Integer[] nonZero(double d) {
			List<Integer> l = new ArrayList<Integer>();
			for(int i=0; i<this.len; i++){
				if(abundance[i]>d && ! this.zerovs.contains(i)) l.add(i);
			}
			return l.toArray(new Integer[0]);
		}
public Stack<Integer> zerovs = new Stack<Integer>();
		
		public void maximisation(double[] v, double pseudo, int numrep, Integer j) {
			if(j!=null) zerovs.push(j);
			for(int i=0; i<numrep; i++){
				this.maximisation(v, pseudo);
			}
			if(j!=null) zerovs.pop();
		}

		public void check() {
			double tot =0;
			for(int i=0; i<this.abundance.length; i++){
				tot+=abundance[i];
			}
			if(Math.abs(tot-1.0)>1e-5){
				throw new RuntimeException("!!! "+tot);
			}
			
		}

		
	}