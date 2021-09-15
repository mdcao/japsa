package japsa.bio.np;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;

public class MultinomialCI {
	// THIS IS A PORT OF THE MULTINOMIALCI LIBRARY IN R. 
	
	 static  BigDecimal zero = new BigDecimal(0);
	  static BigDecimal one = new BigDecimal(1);
	  static BigDecimal two = new BigDecimal(2);
	  static BigDecimal three = new BigDecimal(3);
	  static BigDecimal four = new BigDecimal(4);  
	  static BigDecimal six= new BigDecimal(6);
	  static BigDecimal seven = new BigDecimal(7);
	  static BigDecimal twelve = new BigDecimal(12);
	 static  boolean verbose = false;
	 static boolean debug = false;
	int[] x = null;
	
	RealMatrix salida ;
	public static String getString(Number[] n){
		String format ="%5g";
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<n.length; i++){
			sb.append(n[i]==null ? "null ":String.format(format,  n[i].doubleValue()).trim()+ " ");
		}
		return sb.toString();
	}
	public static void main(String[] v){
		System.err.println(Double.MAX_VALUE);
		//int[] x = new int[] {1000,1,10};
	int[] x = new int[] {232,788,2920};
	//	int[] x = new int []{232000,788000};
		//System.err.println(Arrays.asList(moments(1,30)));
		//System.err.println(getString(moments(1,6000000)));
		if(true){
		MultinomialCI ci = new MultinomialCI(0.05);
		ci.assignCount(x);
		ci.eval();
		RealMatrix res = ci.salida;
		double[][]res1 = ci.tab();
		System.err.println("results1");
		for(int i =0; i<res.getRowDimension(); i++){
			System.err.println(res1[i][0]+" "+res1[i][1]);
		}
		}
	}
	public MultinomialCI(double alpha){
		this.alpha = alpha;
	}
	public  double alpha = 0.05;
	
	
	public static double ppois(int a, double lambda){
		 PoissonDistribution ppois = new PoissonDistribution(lambda); 
		 double d =   ppois.cumulativeProbability(a);
	 return d;
	}
	
	static class PoissonDistribution1 extends PoissonDistribution{
		public PoissonDistribution1(double p) throws NotStrictlyPositiveException {
			super(p);
		}
		public double cp(int x){
			return cumulativeProbability(x);
		}
		public double cp(int a, int b){
			return this.cumulativeProbability(b, a);
			//return cp(a).subtract(cp(b));
		}
		
	}
	
	public static Number[] moments(int c, int lambda){
		PoissonDistribution1 ppois = new PoissonDistribution1(lambda); 
	 int a=lambda+c;
		 int  b=lambda-c;
		   if(b<0){ 
		      b=0;
		   }
		  
		  double den = b>0 ? ppois.cp(a, b-1):  ppois.cp(a);
		  //Number[] mu=new Number[5];
		 Double[] mu=new Double[4];

		   // mom es global y se usa fuera de esta función
		 // BigDecimal lamb = new BigDecimal(lambda);
		   for(int r =1; r<=4; r++){
		      double  poisA=0;
		      double poisB=0;
		      if((a-r) >=0){ poisA=ppois.cp(a,a-r); }
		      if((a-r) < 0){ poisA=ppois.cp(a); }
		      if((b-r-1) >=0){ poisB=ppois.cp(b-1,b-r-1); }
		      if((b-r-1) < 0 && (b-1)>=0){ poisB=ppois.cp(b-1); }
		      if((b-r-1) < 0 && (b-1) < 0){ poisB=0; }
		      Number mu2 = 1-(poisA- poisB)/den;
		    // BigDecimal mu1 = lamb.pow(r);
		    // mu[r] =   mu1.multiply(new BigDecimal(mu2.doubleValue()));
		     mu[r-1] = Math.pow(lambda, r)*mu2.doubleValue();
		   }
		  if(debug) System.err.println(getString(mu));
		
		 
		  return moments(mu, den);
		   //mom[0]=mu[0];
		   //mom[1]=mu[1]+mu[0]-Math.pow(mu[0],2);
		  // mom[2]=mu[2]+mu[1]*(3-3*mu[0])+(mu[0]-3*Math.pow(mu[0],2)+2*Math.pow(mu[0],3));
		  // mom[3]=mu[3]+mu[2]*(6-4*mu[0])+mu[1]*(7-12*mu[0]+6*Math.pow(mu[0], 2))+mu[0]-4*Math.pow(mu[0], 2)
		//		   +6*Math.pow(mu[0],3)-3*Math.pow(mu[0],4);
		 //  mom[4]=den;
		   //return mom;
	}
	
 private static Number[] moments(Double[] mu, double den){
	  Number[] mom=new Number[5];
	 mom[0]=mu[0];
	   mom[1]=mu[1]+mu[0]-Math.pow(mu[0],2);
	   mom[2]=mu[2]+mu[1]*(3-3*mu[0])+(mu[0]-3*Math.pow(mu[0],2)+2*Math.pow(mu[0],3));
	   mom[3]=mu[3]+mu[2]*(6-4*mu[0])+mu[1]*(7-12*mu[0]+6*Math.pow(mu[0], 2))+mu[0]-4*Math.pow(mu[0], 2)
			   +6*Math.pow(mu[0],3)-3*Math.pow(mu[0],4);
	   mom[4]=den;
	   return mom;
 }
	private static BigDecimal[] momemntsBD(BigDecimal[] mu, double den) {
		   Number[] mom=new Number[6];
		 mom[1]=mu[1];
		   mom[2]=add(new BigDecimal[] {mu[2],mu[1],mu[1].pow(2)}, new boolean[] {true, true, false});
		   mom[3]=add(new BigDecimal[] {
				   mu[3],
				   mu[2].multiply(mu[1].multiply(three).subtract(three)),
				   mu[1],
					mu[1].pow(2).multiply(three),
					mu[1].pow(3).multiply(two)},
				   new boolean[] {
						   true, 
						   false, 
						   true, 
						   false, 
						   true
				   });
		   mom[4] = add(new BigDecimal[] {
				   mu[4],
				   mu[1].multiply(four).subtract(six).multiply(mu[3]),
				   mu[2].multiply(mu[1].pow(2).multiply(six).subtract(mu[1].multiply(twelve)).add(seven)),
				   mu[1],
				   mu[1].pow(2).multiply(four),
				   mu[1].pow(3).multiply(six),
				   mu[1].pow(4).multiply(three),
				   
		   }, new boolean[] {
				   true,
				   false,
				   true, 
				   true,
				   false,
				   true,
				   false
				   
		   });
			  mom[5]=den;
		 return   Arrays.asList(mom).subList(1, mom.length).toArray(new BigDecimal[0]);
	}
	static Number add(Number[] d, boolean[] add){
		
		Number res_ = d[0];
		BigDecimal res = (res_ instanceof BigDecimal) ? (BigDecimal) res_ : new BigDecimal(res_.doubleValue());
		
		for(int i=1; i<d.length; i++){
		//	if(d[i].doubleValue()<0) throw new RuntimeException("!!");
			BigDecimal di = (d[i] instanceof BigDecimal) ? (BigDecimal) d[i] : new BigDecimal(d[i].doubleValue());
			if(add==null || add[i])
				res = res.add(di);
			else res = res.subtract(di);
		}
		return res;
	}
	static double sq(double d){
		return Math.pow(d, 2);
	}
	static double cube(double d){
		return Math.pow(d, 3);
	}
	static double fourth(double d){
		return Math.pow(d, 4);
	}
	
	public static double[] colSums(RealMatrix r){
		int d = r.getColumnDimension();
		double [] res = new double[d];
		for(int i=0; i<d; i++){
			double[] v = r.getColumn(i);
			double sum =0; 
			for(int k=0; k<v.length; k++){
				sum +=v[k];
			}
			res[i] = sum;
		}
		return res;
	}
	public static double gamma(double x){
		return Math.exp(Gamma.logGamma(x));
	}
	public static double truncpoi(int c, int[] x, int n, int k){
		PoissonDistribution1 ppois = new PoissonDistribution1(n); 

	RealMatrix m=new Array2DRowRealMatrix(k,5);
	
		   for(int i=0; i<k; i++ ){
		   int lambda=x[i];
		   Number[] mom = moments(c,lambda);
		    for(int j=0; j<5; j++ ){
		     m.setEntry(i,j,mom[j].doubleValue());
		    }
		   }
		   for(int i=0;  i<k; i++){
		    m.setEntry(i,3,m.getEntry(i,3)-3*Math.pow(m.getEntry(i,1),2));
		   }
		  
		 
		   double[] s=colSums(m);
		   double s1=s[0];
		   double  s2=s[1];
		    double s3=s[2];
		    double  s4=s[3];
		 
		 
		     double probn=1.0/(ppois.cp(n,n-1));
		 
		   double z=((double) n-s1)/Math.sqrt(s2);
		   double  g1=s3/Math.pow(s2,1.5);
		  double  g2=s4/Math.pow(s2,2.0);
		 double  poly=1+g1*(Math.pow(z, 3)-3*z)/6.0+g2*(Math.pow(z, 4)-6*Math.pow(z, 2)+3)/24.0
				 +Math.pow(g1,2)*(Math.pow(z, 6)-15.0*Math.pow(z, 4)+45.0*Math.pow(z, 2)-15)/72.0;
		   double f=poly*Math.exp(-Math.pow(z, 2)/2)/(Math.sqrt(2)*gamma(0.5));
		  
		   double probx=1;
		   for(int i=0; i<k; i++){
		    probx=probx*m.getEntry(i,4);
		   }
		   return(probn*probx*f/Math.sqrt(s2));
	}
	public  void  eval (){
		  int n =0;
		//  double[] p = new double[x.length];
		  for(int i=0; i<x.length; i++){
			  n = n+x[i];
		  }
		 // for(int i=0; i<x.length; i++){
		//	  p[i] = x[i]/n;
		 // }
		 int k = x.length;
		 double c = 0;
		 double pold=0;
		 double p = 0;
		 for(int cc=1; cc<=n; cc++){
		   p = truncpoi(cc,x,n,k);
		//System.err.println(cc+" "+p);
		   if( (p > 1-alpha && pold < 1-alpha)) { c = cc; break; };
		   pold=p;
		  }
		 
	
		 double  delta=(1-alpha-pold)/(p-pold);
		 Array2DRowRealMatrix out=new Array2DRowRealMatrix(k,5);
		 Array2DRowRealMatrix num=new Array2DRowRealMatrix(k,1);
		  c=c-1;  
		 double vol1=1;
		 double vol2=1;
		  for(int i=0; i<k; i++ ){
		   num.setEntry(i, 0, i);
		 //  num[i,1]=i;
		   double obsp=(double) x[i]/(double) n;
		   out.setEntry(i, 0,obsp);
		   out.setEntry(i,1,obsp-c/(double) n);
		   out.setEntry(i,2,obsp+c/(double)n+2*delta/(double)n);
		   if(out.getEntry(i,1)<0){ out.setEntry(i,1,0); }
		   if(out.getEntry(i,2)>1){ out.setEntry(i,2,1); }
		   out.setEntry(i,3,obsp-c/(double)n-1/(double)n);
		   out.setEntry(i,4,obsp+c/(double)n+1/(double)n);
		   if(out.getEntry(i,1)<0){ out.setEntry(i,1,0); }
		   if(out.getEntry(i,2)>1){ out.setEntry(i,2,1); }
		   vol1=vol1*(out.getEntry(i,2)-out.getEntry(i,1));
		   vol2=vol2*(out.getEntry(i,4)-out.getEntry(i,3));
		   salida.setEntry(i,0,out.getEntry(i,1));
		   salida.setEntry(i,1, out.getEntry(i,2));
		   this.order[i] = new DoubleInt(out.getEntry(i, 2),i);
		  }
		  String[] c1=new String[] {"PROPORTION", "LOWER(SG)", "UPPER(SG)","LOWER(C+1)","UPPER(C+1)"};
		  double  cov=100*(1-alpha);
		  Double[] sg=new Double[x.length];
		  for(int i=0; i<sg.length; i++){
			 sg[i] =  ((double)x[i]+delta)/(double)n;
		  }
		  String[] c2=new String[] {"SG-midpoint"} ;
		  if(verbose){
		    System.err.println("-------------------------------------------------------------");
		    System.err.println("    "+cov+"% SIMULTANEOUS CONFIDENCE INTERVALS");
		    System.err.println("       BASED ON THE METHODS OF SISON AND GLAZ");
		    System.err.println("-------------------------------------------------------------");
		    System.err.println("C = "+c);
		    System.err.println("P(c+1) = "+p);
		    System.err.println("P(c)   = "+pold);
		    System.err.println("delta =  "+delta);
		    System.err.println("Volume(SG) = "+vol1);
		    System.err.println("Volume(C+1)= "+vol2);
		    System.err.println(Arrays.asList(c1));
		    for(int i=0; i<out.getRowDimension(); i++){
		    System.err.println(out.getRowVector(i));
		    }
		    //System.err.println(paste(c(num+out)));
		   // System.err.println(paste(c2));
		    //System.err.println(paste(c(num,sg)));
		    System.err.println(Arrays.asList(sg));
		  }
		}
	public void assignCount(double[] array) {
		this.x = new int[array.length];
		for(int j=0; j<array.length; j++){
			x[j] = (int) Math.round(array[j]);
		}
		salida = new Array2DRowRealMatrix(x.length,2);
		order = new DoubleInt[array.length];
	}
	public void assignCount(int[] array) {
		this.x = new int[array.length];
		for(int j=0; j<array.length; j++){
			x[j] = (int) Math.round(array[j]);
		}
		salida = new Array2DRowRealMatrix(x.length,2);
		order = new DoubleInt[array.length];
		
	}
	public double[][] tab() {
		return salida.getData();
	}
	static class DoubleInt implements Comparable{
		double v; int ind;
		DoubleInt(double v, int i){
			this.v = v; this.ind = i;
		}
		@Override
		public int compareTo(Object o) {
			return Double.compare(((DoubleInt)o).v,v);
		}
	}
	
	DoubleInt[] order;
	public int[] rank(){
		/*if(l!=null){
			for(int i=0; i<order.length; i++){
				order[i].v = l.get(i);
			}
		}*/
		Arrays.sort(order);;
		int[] res = new int[order.length];
		for(int i=0; i<res.length; i++){
			res[i] = order[i].ind;
		}
		return res;
	}
}
