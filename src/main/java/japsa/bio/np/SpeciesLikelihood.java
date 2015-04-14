package japsa.bio.np;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class SpeciesLikelihood {

	//double genomeLength, logGenomeLength;
	double threshLength=0;
	Map<String, Integer>allgenes = new HashMap<String, Integer>();
	int totalGeneLength =0;
	int noGenes =0;

	SpeciesLikelihood(){

	}

	//double currSum=0;

	//first line says genome length;
	//remaining lines list gene start end

	public String species;
	SpeciesLikelihood(String species ,  Iterator<String> it, SpeciesLikelihood bg){
		this.species = species;
		try{
			//		 BufferedReader br = new BufferedReader(new FileReader(genes));
			String st = "";//br.readLine();
			//this.genomeLength = Double.parseDouble(st);
			while(it.hasNext() && (st = it.next())!=null){
				String[] str = st.split("\t");
				// if(str.length<=1){
				// System.err.println("skipping "+st);
				// continue;
				//}
				int len = 100;
				if(str.length>11){
					len =  Integer.parseInt(str[2]) - Integer.parseInt(str[1]);
				}			
				{
					totalGeneLength +=len;
					noGenes++;
					Integer v = allgenes.get(str[0]);
					if(v==null) allgenes.put(str[0],len);
					else allgenes.put(str[0],v+len);
				}{
					Integer v1 = bg.allgenes.get(str[0]);
					if(v1==null) bg.allgenes.put(str[0],len);
					else bg.allgenes.put(str[0],v1+len);
					bg.totalGeneLength+=len;
					bg.noGenes++;
				}
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	//	public List<Double> readLikelihoods = new ArrayList<Double>();

	public SpeciesLikelihood(String species, int st, int end, SpeciesLikelihood bg, boolean sample) {
		this(species, mkList(st,end, sample), bg);
	}

	public static int sample(int st, int end){
		return st +(int) Math.floor(Math.random() *(end-st));
	}

	private static Iterator<String> mkList(int st, int end, boolean sample) {
		List<String> l = new ArrayList<String>();
		for(int i = st; i<end; i++){
			int k= sample ? sample(st, end): i;
			int len = SpeciesLikelihood.sample(0, 200);
			l.add("gene"+k+"\t0\t"+len);
		}
		return l.iterator();
	}

	public SpeciesLikelihood(File file, SpeciesLikelihood bg)  throws Exception{
		this(file.getName(), getIterator(file), bg);
	}


	public SpeciesLikelihood(String spname, List<String> l, SpeciesLikelihood bg)  throws Exception{
		this(spname, l.iterator(), bg);
	}
	private static Iterator<String> getIterator(File file)  throws Exception{
		final BufferedReader br = new BufferedReader(new FileReader(file));
		List<String> l = new ArrayList<String>();
		String st = "";
		while((st = br.readLine())!=null){
			l.add(st);
		}
		return l.iterator();

	}

	//use -ve number for no overlap.
	//calculates likelihood of read matching overlap and not matching non-overlaps
	double likelihood(double readLength, String gene){
		//		double p = (geneLength + (readLength - threshLength))/(sum_k geneLengthk +(readLength - threshLength); 
		//	double q = 1 - p;
		//	double logpsum = 0;
		//	double logqsum = 0;
		Integer geneLength = allgenes.get(gene);
		double likelihood = 0;
		if(geneLength!=null){
			likelihood = (geneLength +(readLength - threshLength))/(this.totalGeneLength + noGenes*(readLength - threshLength));
		}
		return likelihood;
		//		readLikelihoods.add(likelihood);
		//		currSum+=;
		//		return sum;
	}

	public List<String> sampleGenes(int n) {
		List<String> res = new ArrayList<String>(0);
		List<String> genes = new ArrayList<String>(this.allgenes.keySet());
		for(int k=0; k<n; k++){
			int samp = sample(0,genes.size());
			res.add(genes.get(samp));
		}
		return res;

	}



}
