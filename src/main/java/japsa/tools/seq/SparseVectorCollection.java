package japsa.tools.seq;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.Stack;
import java.util.TreeMap;

import japsa.bio.np.RealtimeSpeciesTyping;

public class SparseVectorCollection{
		// each entry in svs represents a multimapping read
		private  final double[] abundance;
		public final SparseVector[] abundances;
		int len;
		int total_size=0;
		List<String> species;
		
		SortedMap<Number, String> getScores(SparseVector sv){
		//	SparseVector sv = this.svs[k].get(i);
			Iterator<Integer> it = sv.keyIt();
			SortedMap< Number, String> m = new TreeMap<Number, String>();
			while(it.hasNext()){
				Integer nxt = it.next();
				System.err.println(species.get(nxt)+" "+sv.get(nxt));
				if(species.get(nxt).startsWith("P")){
					System.err.println("HHH");
				}
				Number n =sv.get(nxt); 
				m.put(n,species.get(nxt) );
			}
			return m;
		}
		
		public SparseVectorCollection(List<String> species, boolean keepNames,String[] nmes){
			int numSpecies = species.size();
			this.species = species;
			this.src_names = nmes;
			this.num_sources =src_names.length;
			this.len = numSpecies;
			this.abundance = new double[numSpecies];
			this.abundances = new SparseVector[num_sources];
			double val = 1.0/(double) numSpecies;
			Arrays.fill(abundance, val);
			for(int j=0; j<num_sources; j++){
				abundance[j] = val;
				abundances[j] = new SparseVector();
			}
			this.keepNames = keepNames;
			svs = new List[num_sources];
			this.read_nmes = new List[num_sources];

			for(int i=0; i<svs.length; i++){
				svs[i] = new ArrayList<SparseVector>();
				this.read_nmes[i] = new ArrayList<String>();
			}
		}
		final String[] src_names;
		final List<SparseVector>[] svs;
	//	List<Integer> single = new ArrayList<Integer>();
		final List<String>[] read_nmes ; // read nes
		final public boolean keepNames; // whether to record the read names
		
		public void add(SparseVector all_species, String nme, int src_index) {
			Set<Integer> keys = all_species.keySet();
			if(keepNames){
				this.read_nmes[src_index].add(nme);
			}
			total_size++;
			if(keys.size()==1 && false){
		//		single.add(keys.iterator().next());
			}else{
				svs[src_index].add(all_species);
			}
		}
		final public  int num_sources;
		public double logLike(){
			double sc =0;
			for(int k=0; k<this.num_sources; k++)
			for(int i=0; i<this.svs[k].size(); i++){
				SparseVector sv = svs[k].get(i); // this represents a read
				double pr = 0;
				for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
					Integer j = it.next();
					pr+= abundance[j] * sv.score(j);
				}
				sc+=Math.log(pr);
			}
			return sc;
		}
		
		public void printMostLikely( List<String> species, Map<String, String>[] out) throws IOException{
			if(!this.keepNames) return;
			SparseVector sv1 = new SparseVector();
			for(int k=0; k<this.num_sources; k++){
			//	PrintWriter pw = new PrintWriter(new FileWriter(out[k]));

				for(int i=0; i<this.svs[k].size(); i++){
					 Map.Entry<Integer, Number>  entr = this.mostLikely(k,i, sv1);
					double prob = entr.getValue().doubleValue()/sv1.valsum();
					Integer key = entr.getKey();
					out[k].put(this.read_nmes[k].get(i), species.get(key)+"\t"+String.format("%5.3g",prob).trim());
				//	pw.println(this.read_nmes[k].get(i)+"\t"+this.src_names[k]+"\t"+species.get(key)+"\t"+String.format("%5.3g",prob).trim());
				}
				//pw.close();

			}
		}
		
		public Map.Entry<Integer, Number> mostLikely(int src_index, int index, SparseVector sv1){
			SparseVector sv = this.svs[src_index].get(index);
			sv1.clear();
			for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
				Integer j = it.next();
				sv1.addToEntry(j,abundance[j] * sv.score(j));
			}
			return sv1.getMaxEntry();
			
		}
		
		public void getAbundances(){
			for(int k=0; k<this.num_sources; k++){
				SparseVector abund = this.abundances[k];
			for(int i=0; i<this.svs[k].size(); i++){
				SparseVector sv = svs[k].get(i);
				double totj = 0;
			//	double scorej=0;
				//split the read count according to the quality scores
				for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
					Integer j = it.next();
					double sc = abundance[j] * sv.score(j);
					totj+=sc;
				//	scorej+=sc;
				}
				if(totj>0){
					for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
						Integer j = it.next();
						if(abundance[j]>0){
							double sc = abundance[j] * sv.score(j)/totj;
							if(sc>0) abund.addToEntry(j, sc);
						}
					}
					//score+=Math.log(scorej);
				}else{
				//	tot = tot -1;
					System.err.println("warning totj is zero "+i+" "+sv.toString());
				}
				abund.normalise();
			}
			}
			
			//SortedMap<Number,String > m = this.getScores(this.abundances[1]);
			//SortedMap<Number,String > m1 = this.getScores(this.svs[1].get(0));
		//	m1.get(m.lastKey());
		//	System.err.println(m1);
		//	System.err.println("calculated abundances");
		}
		
		//update abundance
		//one round of normalisation
		public void maximisation(double[] v, double pseudo){
			double[] abund = new double[len];
			
			Arrays.fill(abund, pseudo);
			if(total_size==0){// && single.size()==0){
				return;//throw new RuntimeException("nothing");
			}
			//Map <Integer, Number> max = svs[0].get(0).getMax();
			//int sze = max.size();
			double score =0;
			for(int k=0; k<svs.length; k++){
			for(int i=0; i<this.svs[k].size(); i++){
				SparseVector sv = svs[k].get(i);
				
				double totj = 0;
				double scorej=0;
				//split the read count according to the quality scores
				for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
					Integer j = it.next();
					double sc = abundance[j] * sv.score(j);
					totj+=sc;
					scorej+=sc;
				}
				if(totj>0){
					for(Iterator<Integer> it = sv.keyIt();it.hasNext();){
						Integer j = it.next();
						double sc = abundance[j] * sv.score(j)/totj;
						abund[j] += sc;
					}
					score+=Math.log(scorej);
				}else{
				//	tot = tot -1;
					System.err.println("warning totj is zero "+i+" "+sv.toString());
				}
				
			}
			}
	
			System.arraycopy(abund, 0, abundance, 0, abundance.length);
			
			
			v[0] = score; 
			if(zerovs.size()>0){ // if we keeping anything zero
				Iterator<Integer> iter = zerovs.iterator();
				while( iter.hasNext()){
					Integer nxt = iter.next();
					this.abundance[nxt] = RealtimeSpeciesTyping.epsilon;
					
				}
				this.renormalise(this.abundance);
			}else{
				this.renormalise(this.abundance);
			}
			/*	double diff = 0;
				for(int j=0; j<abund.length; j++){
				double newv = abund[j]/tot;
				diff += Math.abs(abundance[j] - newv);
				this.abundance[j] =newv ;
				}	
				
				v[1] = diff;*/
			
		}
		
		public double setZero(int i, double v1){
			double v = abundance[i];
			abundance[i] =v1;
			this.renormalise(this.abundance);
			return v;
		}
		
		static void renormalise(double[] abundance) {
			double tot = 0;
			for(int i =0; i<abundance.length; i++){
				tot+=abundance[i];
			}
			for(int i =0; i<abundance.length; i++){
				abundance[i] = abundance[i]/tot;
			}
			check(abundance);
			
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
		//	double prev=0;
			double[] abund = new double[this.abundance.length];
			inner: for(int i=0; i<numrep; i++){
				this.maximisation(v, pseudo);
				v[1] = diffAndTransfer(abund, abundance);
				
				if(v[1]<1e-5 && i>1)  break inner;
			}
			if(j!=null) zerovs.pop();
		}

		private double diffAndTransfer(double[] abund, double[] abundance) {
			double diff = 0;
			for(int i=0; i<abund.length; i++){
				diff+= Math.abs(abundance[i] - abund[i]);
				abund[i] = abundance[i];
			}
			return diff;
		}

		private int max_ind(double[] a) {
			 int ind =0;
			 for(int i=1; i<a.length;i++){
				 if(a[i]>=a[ind]){
					
					 ind = i;
				 }
			 }
			 return ind;
		}

		static void check(double[] abundance) {
			double tot =0;
			for(int i=0; i<abundance.length; i++){
				tot+=abundance[i];
			}
			if(Math.abs(tot-1.0)>1e-5){
				throw new RuntimeException("!!! "+tot);
			}
			
		}

		
	}