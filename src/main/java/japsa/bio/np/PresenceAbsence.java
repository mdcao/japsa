/**
 * LICENCE?
 */
package japsa.bio.np;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class PresenceAbsence {


	static boolean RESCALE = true;
	public static void main(String[] args){
		try{
			PresenceAbsence pa = new PresenceAbsence(new File(args[0]));
			pa.likelihood(100,"ABC10");
			pa.likelihood(100,"BAF20");
			// pa = new PresenceAbsence(new String[] {"1-5000", "3-5003"}, false);
			// List<String> genes = pa.spl[0].sampleGenes(50000);
			// for(int k=0; k<genes.size(); k++){

			//}
			//String[] nme = pa.species;
			double[] posterior = pa.calcPosterior();
			double[][] samp =pa.calcPosterior(100);
			double[][] ranges = pa.getRanges(samp);
			for(int i=0; i<posterior.length; i++){
				System.out.println(pa.spl[i].species+" "+posterior[i]+" "+ranges[i][0]+" "+ranges[i][1]);
			}

		}catch(Exception exc){
			exc.printStackTrace();
		}

	}

	public double[] calcPosterior() {
		calcPosterior(this.sumLogL, this.posterior);
		return this.posterior;
	}

	public SpeciesLikelihood[] spl;
	SpeciesLikelihood bg; //background null model;
	double[] posterior;

	double mixp = 0.8;
	double mixq = 1-mixp;

	List<Double>[] likelihoods;
	List<Double> bg_likelihood = new ArrayList<Double>();

	double[] sumLogL;
	double bg_sumLogL=0;

	double[] sample_sumLogL;
	double sample_bg_sumLogL=0;
	double[][] sample_posterior;

	static FilenameFilter fnFilter = 
		new FilenameFilter(){
		@Override
		public boolean accept(File arg0, String arg1) {
			return arg1.endsWith(".txt");
		}

	};

	PresenceAbsence(File dir){
		this(dir.listFiles(fnFilter));
	}

	PresenceAbsence(File[] f){
		this(f.length);
		try{
			for(int i=0; i<spl.length; i++){
				spl[i] = new SpeciesLikelihood(f[i], bg);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	PresenceAbsence(String[] speciesnames, List<String>[] genelists){
		this(genelists.length);
		try{
			for(int i=0; i<spl.length; i++){
				spl[i] = new SpeciesLikelihood(speciesnames[i], genelists[i], bg);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	/**
	 * Minh added based on the above constructor
	 * @param profileList
	 */
	public PresenceAbsence(ArrayList<RealtimeStrainTyping.GeneProfile> profileList){
		this(profileList.size());
		try{
			for(int i=0; i<profileList.size(); i++){
				spl[i] = new SpeciesLikelihood(profileList.get(i).strainID(), profileList.get(i).getGeneList().iterator(), bg);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}


	PresenceAbsence(String[] ranges, boolean sample){
		this(ranges.length);
		try{
			for(int i=0; i<ranges.length; i++){
				String[] str = ranges[i].split("-");
				spl[i] = new SpeciesLikelihood("species"+i,Integer.parseInt(str[0]),Integer.parseInt(str[1]),  bg, sample);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}



	@SuppressWarnings("unchecked")
	PresenceAbsence(int len){
		spl = new SpeciesLikelihood[len];
		bg = new SpeciesLikelihood();
		likelihoods = new List[len];
		sumLogL = new double[len];
		this.sample_sumLogL = new double[len];
		Arrays.fill(sumLogL, 0);
		posterior = new double[len];
		for(int i=0; i<len; i++){
			likelihoods[i] = new ArrayList<Double>();
		}
	}


	//double likelihoodAdded=0;

	public void likelihood(double readLength, String gene){
		double bg_lhood = bg.likelihood(readLength, gene);
		bg_sumLogL+=Math.log(bg_lhood);
		this.bg_likelihood.add(bg_lhood);
		for(int i=0; i<spl.length; i++){
			double likelihood = spl[i].likelihood(readLength, gene)*mixp + bg_lhood * mixq;
			double logLikelihood = Math.log(likelihood);
			this.sumLogL[i]+=logLikelihood;
			this.likelihoods[i].add(logLikelihood);

		}
	}



	public static void calcPosterior(double[] sumLogL, double[] posterior){
		double s = 0;
		double maxL = 0;
		int spllen = posterior.length;
		if(RESCALE){
			maxL = Double.NEGATIVE_INFINITY;
			for(int i=0; i<spllen; i++){
				if(sumLogL[i] > maxL){
					maxL = sumLogL[i];
				}
			}
		}
		for(int i=0; i<spllen; i++){
			posterior[i] = Math.exp(sumLogL[i]-maxL); 
			s+=posterior[i];
		}
		for(int i=0; i<spllen; i++){
			posterior[i] = posterior[i]/s;
		}

	}

	public double[][] getRanges(double[][] post, double confidentInterval){
		int nspec = post[0].length;
		int nrep = post.length;
		double[][] res = new double[nspec][2];
		double[] samples = new double[nrep];
		for(int k=0; k<nspec; k++){
			for(int j=0; j<nrep; j++){
				samples[j] = post[j][k];
			}
			Arrays.sort(samples);
			int index0 = (int) (nrep * (1.0 - confidentInterval) / 2);
			if (index0 < 1)
				index0 = 1;

			int index1 = nrep - index0 -1;
			if (index1 <= index0)
				throw new RuntimeException("Confident interval of " + confidentInterval + " does not work with sample = " + nrep);
			//LOG.info("YYY " + index0 + " " + index1 + " from " + nrep + " and " +confidentInterval );
			res[k][0] = samples[index0];
			res[k][1] = samples[index1];
		}
		return res;
	}

	public double[][] getRanges(double[][] post){
		int nspec = post[0].length;
		int nrep = post.length;
		double[][] res = new double[nspec][2];
		double[] samples = new double[nrep];
		for(int k=0; k<nspec; k++){
			for(int j=0; j<nrep; j++){
				samples[j] = post[j][k];
			}
			Arrays.sort(samples);
			res[k][0] = samples[0];
			res[k][1] = samples[nrep-1];
		}
		return res;
	}

	public double[][] calcPosterior(int numreps){
		sample_posterior = new double[numreps][this.spl.length];
		for(int k=0; k<numreps; k++){
			samplePosterior(sample_posterior[k], this.bg_likelihood.size());
		}
		return sample_posterior;
	}

	public void samplePosterior(double[] sample_post, int n) {
		Arrays.fill(sample_sumLogL, 0);
		//	double sum=0;
		//	int n = this.bg_likelihood.size();
		//if(n<5) throw new RuntimeException("must have at least 5 reads");
		for(int i=0; i<n; i++){
			int s = SpeciesLikelihood.sample(0, n);;
			if(s!=n){
				for(int k=0; k<this.spl.length; k++){
					this.sample_sumLogL[k]+=this.likelihoods[k].get(s);
				}
			}else{
				i--;
			}

		}
		calcPosterior(sample_sumLogL, sample_post);

	}

}
