package japsadev.tools;

public class VNTRGenotyper {
	/* test */
	public static void main(String[] args){
		VNTRGenotyper vg = new VNTRGenotyper();
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
	 private void setGenos(){
		 double refAllele = this.genotype_reference;
		 double ratio = ( this.count_repeat_sample/ this.count_flanking_sample)/(this.count_repeat_ref/ this.count_flanking_ref);
		 ratio = Math.max(2*ratio, 4);
		 genos = new double[(int) Math.max(40, Math.ceil(ratio*refAllele))];
			double rem = refAllele - Math.floor(refAllele);
		//	int half = genos.length/2;
			for(int x=0; x<genos.length; x++){
				genos[x] = x+rem;
			}
			prob = new double[genos.length];
	 }
	 
	public String getConf(double conf){
			this.setGenos();
			probability(prob, genos);
			int[] range = new int[2];
			Double[] mass = getconf(conf, range);
			return String.format("%5.3g,%5.3g,%5.3g",mass).replaceAll("\\s+", ""); 
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
	
	
	Double[] getconf( double mass, int[] range){
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
	void setRef(double genotype_reference, double  count_repeat_ref, double  count_flanking_ref){
		 mult = 100.0/genotype_reference;
		this.genotype_reference = 100;
		this.count_repeat_ref = count_repeat_ref;
		this.count_flanking_ref = count_flanking_ref;
		
	}
	
	void setSample(double  count_repeat_sample, double  count_flanking_sample){
		this.count_flanking_sample = count_flanking_sample;
		this.count_repeat_sample = count_repeat_sample;
	}
	
	
	/** This calculates the log likelihood */
	double likelihood(double genotype){
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
	
	
	void probability(double[] logprobs, double[] genos){
		double sum=0;
		//double maxv =0;
		for(int i=0; i<logprobs.length; i++){
			double v = likelihood(genos[i]);
		//	if(v>maxv) maxv = v;
			logprobs[i] = v;
			sum+=Math.exp(v);
		}
		//NEED TO CHECK THIS IS RIGHT - IDEA IS TO TRANSFORM LOG PROBS TO PROBS
		for(int i=0; i<logprobs.length; i++){
			logprobs[i] = Math.exp(logprobs[i])/sum;
		}
		//for(int j =0; j<logprobs.length; j++){
		//	sum+=Math.exp(logprobs[j]-maxv);
		//}
		
	}
	

	
}
