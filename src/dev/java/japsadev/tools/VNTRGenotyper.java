package japsadev.tools;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VNTRGenotyper {
	public VNTRGenotyper(double downsample2) {
		this.downsample = downsample2;
	}

	/* test */
	public static void main(String[] args){
		VNTRGenotyper vg = new VNTRGenotyper(1);
		double certainty = 0.5;
		double[] ref = new double[] {5,10,15};
		double[] sample = new double[] {5,10,15} ;
		double[] depth  = new double[] {500};
		for(int i=0; i<ref.length; i++){
			for(int j=0; j<sample.length; j++){
				for(int k=0; k<depth.length; k++){
					vg.simulate(ref[i], sample[j], depth[k], certainty);
				}
			}
		}
//		double ref = 10;
//		double sample = 14;
	//	double depth = 100;
		//double mult = 1;
		
		
	}
	
	void simulate(double ref, double sample, double depth, double targetCertainty){
		setRef(ref, depth/2 , depth/2);
		setSample((sample/ref)*(depth/2), depth/2);
		
		
		double[] genos = new double[20];
		int half = genos.length/2;
		for(int i=0; i<genos.length; i++){
			genos[i] = sample - half +i;
		}
		double[] prob = new double[genos.length];
		probability(prob, genos);
		System.err.println("Simulated  ref:"+ref+" sample:"+sample+" depth:"+depth);
		int[] range = new int[2];
		
		Double[] mass = getconf( targetCertainty, range);
		double[] perc  = new double[] { (sample-genos[range[0]])/sample, (genos[range[1]]-sample)/sample};
		System.err.println(genos[range[0]] +" to "+genos[range[1]]+ "(+/-"+String.format("%2.2g", perc[0]*100)  +"%) certainty:"+String.format("%5.3g", mass[2])+"        max:"+mass[0]+ " maxp:"+String.format("%5.3g", mass[1]));
		System.err.println("---------------------------------");
		//for(int i=range[0]; i<=range[1]; i++){
		//	System.err.println(genos[i]+" => "+prob[i]);
		//}
	}
	/**NOTE:  THIS MODEL IS SPECIFICALLY DESIGNED TO COUNT READS.  IF YOU ARE COUNTING BASES IT IS CRITICAL TO DIVIDE THE COUNTS BY THE AVERAGE LENGTH OF THE READS
			THE REASON FOR THIS IS THAT THE BETABINOMIAL REFLECTS THE UNCERTAINTY IN THE ESTIMATE OF THE PROPORTION OF READS IN THE REFERENCE IN FLANKING VS NORMAL
		    IF YOU ARE COUNTING BASES IT GIVES AN ARTIFICIALLY HIGH DEGREE OF CERTAINTY IN THIS FRACTION.	
	
	**/
	static  BetaBinomial bb = new BetaBinomial();
	
	 
	 
	 private double[] genos, prob;
	 
	 /** just gets the array of genos and probs ready */
	 private boolean setGenos(){
		 double refAllele = this.genotype_reference;
		 double ratio = ( this.count_repeat_sample/ this.count_flanking_sample)/(this.count_repeat_ref/ this.count_flanking_ref);
		 ratio = Math.max(2*ratio, 4);  // make genos at least twice the expected
		 if(Double.isNaN(ratio) || Double.isInfinite(ratio)) return true;
		 genos = new double[(int) Math.max(40, Math.ceil(ratio*refAllele))];  // list of possible genotypes
		/* if(genos.length==0){
			 System.err.println(ratio);
			 System.err.println(refAllele);
			 throw new RuntimeException("this should not happen");
		 }*/
			double rem = refAllele - Math.floor(refAllele);
		//	int half = genos.length/2;
			for(int x=0; x<genos.length; x++){
				genos[x] = x+rem;
			}
			prob = new double[genos.length];
			return false;
	 }
	 
	public String getConf(double conf){
			boolean NA = this.setGenos();
			if(NA){
				return "NA,NA,NA";
			}
			double rsd = probability(prob, genos);
			int[] range = new int[2];
			Double[] mass = getconf(conf, range);
			return String.format("%5.3g,%5.3g,%5.3g",mass).replaceAll("\\s+", "")+String.format(",%5.3g",  rsd).trim(); 
		//	return String.format("%5.3g", mass[0]).trim()+"-"+String.format("%5.3g",  mass[1]).trim();
	 }
	
	public String getConf1(double perc){
		this.setGenos();
		probability(prob, genos);
		Double[] mass = getconf(perc);
		return String.format("%5.3g,%5.3g,%5.3g",mass).replaceAll("\\s+", ""); 
	}
	
	private Double[] getconf( double perc){
		int maxi=0;
		for(int i=1; i<prob.length; i++){
			if(prob[i] >prob[maxi]) maxi =i;
		}
		double mlgeno = genos[maxi];
		int range = (int) Math.floor(mlgeno * perc); // this assumes genotypes in steps of 1
		
		double sum=0;
		for(int i=maxi-range; i<=maxi+range; i++){
			sum +=prob[i];
		}
	//	double[] ranges = new double[] {maxi-range,maxi+range};
		return new Double[] {((double) mlgeno-range)/mult, ((double) mlgeno+range)/mult, sum};
		
	}
	
	
	private Double[] getconf( double mass, int[] range){
		int maxi=0;
		for(int i=1; i<prob.length; i++){
			if(prob[i] >prob[maxi]) maxi =i;
		}
		int i=0;
		int len = prob.length;
		double sum=prob[maxi];
		range[0] = maxi;
		range[1] = maxi;
		if(sum<mass){
		for(i=1; maxi-i>0 && maxi +i <len ; i++ ){
			double add = prob[maxi-i]+prob[maxi+i];
			sum += add;
			range[0] =maxi-i;
			range[1] = maxi+i;
			if(sum>=mass){
				break;
			}
		}
		}
		
		return new Double[] {genos[range[0]]/mult, genos[range[1]]/mult, sum,  genos[maxi]/mult,  prob[maxi]};
		
	}
	/*bdw is a factor that can be used to make the genotype calls more uncertain.  It does this by artificially decreasing the counts. So a value of 10, for example
	 * would decrease the counts (both in repeat and flanking) by a factor of 10.  Note that its best to change this on a log scale if you want to see any affect, i.e. going from 1 to 2
	 * has little effect, need to change from 1 to 10
	 *  */
	static double bdw = 1;
	
	double genotype_reference;
	double count_repeat_ref;
	double count_flanking_ref;
	double count_repeat_sample;
	double count_flanking_sample;
	double mult = 1;

 private final double downsample; // this is just to see the effect of randomly downsampling
 
 static double referencelevel =  1000.0; // this is an arbitrary constant to set the reference allele count to.  It controls the level of resolution in the VNTR length calls.
	void setRef(double genotype_reference, double  count_repeat_ref, double  count_flanking_ref){
		 mult =referencelevel/genotype_reference;
		this.genotype_reference = referencelevel;
		if(Math.abs(downsample-1)>1e-5){
			double n = count_repeat_ref + count_flanking_ref;
			double p= count_repeat_ref/n;
			double n_new =  Math.round(n/downsample);
			BinomialDistribution bin = new BinomialDistribution((int) n_new,p);
			double n1 = bin.sample();
			this.count_repeat_ref = n1;
			this.count_flanking_ref = n_new-n1;
		}else{
			this.count_repeat_ref = count_repeat_ref;
			this.count_flanking_ref = count_flanking_ref;
		}
	}
	
	void setSample(double  count_repeat_sample, double  count_flanking_sample){
		if(Math.abs(downsample-1)>1e-5){
			double n = count_repeat_sample + count_flanking_sample;
			double p= count_repeat_sample/n;
			double n_new =  Math.round(n/downsample);
			BinomialDistribution bin = new BinomialDistribution((int) n_new,p);
			double n1 = bin.sample();
			this.count_repeat_sample = n1;
			this.count_flanking_sample= n_new-n1;
		}else{
			this.count_flanking_sample = count_flanking_sample;
			this.count_repeat_sample = count_repeat_sample;
		}
		
	}
	
	
	/** This calculates the log likelihood */
	double loglikelihood(double genotype){
		double relative_cn = genotype/genotype_reference;
		double n = count_repeat_sample + count_flanking_sample;
		bb.set((relative_cn*count_repeat_ref)/bdw, count_flanking_ref/bdw,n);
		return bb.logdensity(count_repeat_sample);
		
	}
	
	double[] probability(int max_cn){
		double[] res = new double[max_cn+1];
		double[] genos = new double[max_cn+1];
		for(int i=0; i<genos.length; i++) genos[i] = i;
		probability(res, genos);
		return res;
	}
	
	
	double probability(double[] probs, double[] genos){
		double sum=0;
		//double maxv =0;
		double mean =0;
		for(int i=0; i<probs.length; i++){
			double v = Math.exp(loglikelihood(genos[i]));
		//	if(v>maxv) maxv = v;
			probs[i] = v;
			mean = mean + genos[i] * probs[i];
			sum+=v;
		}
		mean = mean/sum;
		//NEED TO CHECK THIS IS RIGHT - IDEA IS TO TRANSFORM LOG PROBS TO PROBS
		double var =0;
		for(int i=0; i<probs.length; i++){
			//probs[i] = Math.exp(probs[i])/sum;
			probs[i] =probs[i]/sum;
			var = var + Math.pow(genos[i] - mean,2)*probs[i];
		}
		double rsd = Math.sqrt(var)/mean;
		return rsd;
		//for(int j =0; j<logprobs.length; j++){
		//	sum+=Math.exp(logprobs[j]-maxv);
		//}
		
	}

	public String getConfs(double[] cI) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<cI.length; i++){
			sb.append(getConf(cI[i]));
			if(i<cI.length-1)sb.append(";");
		}
		return sb.toString();
	}
	

	
}
