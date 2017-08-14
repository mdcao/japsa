package japsa.bio.np;
import java.util.Arrays;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;

public class MultinomialCI {
	// THIS IS A PORT OF THE MULTINOMIALCI LIBRARY IN R. 
	int[] x = null;
	
	RealMatrix salida ;

	
	public static void main(String[] v){
		int[] x = new int[] {1000,1,10};
		MultinomialCI ci = new MultinomialCI(0.05);
		ci.x = x;
		ci.eval();
		RealMatrix res = ci.salida;
		for(int i =0; i<res.getRowDimension(); i++){
			System.err.println(res.getRowVector(i));
		}
	}
	public MultinomialCI(double alpha){
		this.alpha = alpha;
	}
	public  double alpha = 0.05;
	public  boolean verbose = false;
	public static double ppois(int a, double lambda){
		 PoissonDistribution ppois = new PoissonDistribution(lambda); // this should be lambda, not mean
		 double d =   ppois.cumulativeProbability(a);
	 return d;
	}
	public static double[] moments(int c, int lambda){
	 int a=lambda+c;
		 int  b=lambda-c;
		   if(b<0){ 
		      b=0;
		   }
		   
		   double den = b>0 ? ppois(a,lambda)-ppois(b-1,lambda):  ppois(a,lambda);
		   //if(b>0){ 
		      //den=poisson(lambda,a)-poisson(lambda,b-1);
		     //  den=ppois(a,lambda)-ppois(b-1,lambda);
		   //}
		   //if(b==0){
		      // den=poisson(lambda,a);
		     //  den=ppois(a,lambda);
		   //}
		   double[] mu=new double[4];
		   // mom es global y se usa fuera de esta función
		 double[] mom=new double[5];
		   for(int r =1; r<=4; r++){
		      double poisA=0;
		      double poisB=0;
		      if((a-r) >=0){ poisA=ppois(a,lambda)-ppois(a-r,lambda); }
		      if((a-r) < 0){ poisA=ppois(a,lambda); }
		      if((b-r-1) >=0){ poisB=ppois(b-1,lambda)-ppois(b-r-1,lambda); }
		      if((b-r-1) < 0 && (b-1)>=0){ poisB=ppois(b-1,lambda); }
		      if((b-r-1) < 0 && (b-1) < 0){ poisB=0; }
		      mu[r-1]=(Math.pow(lambda,r))*(1-(poisA-poisB)/den);
		      
		    
		   }
		   mom[0]=mu[0];
		   mom[1]=mu[1]+mu[0]-Math.pow(mu[0],2);
		   mom[2]=mu[2]+mu[1]*(3-3*mu[0])+(mu[0]-3*Math.pow(mu[0],2)+2*Math.pow(mu[0],3));
		   mom[3]=mu[3]+mu[2]*(6-4*mu[0])+mu[1]*(7-12*mu[0]+6*Math.pow(mu[0], 2))+mu[0]-4*Math.pow(mu[0], 2)
				   +6*Math.pow(mu[0],3)-3*Math.pow(mu[0],4);
		   mom[4]=den;
		   return mom;
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
	RealMatrix m=new Array2DRowRealMatrix(k,5);
		   for(int i=0; i<k; i++ ){
		   int lambda=x[i];
		   double[] mom = moments(c,lambda);
		    for(int j=0; j<5; j++ ){
		     m.setEntry(i,j,mom[j]);
		    }
		   }
		   for(int i=0;  i<k; i++){
		    m.setEntry(i,3,m.getEntry(i,3)-3*Math.pow(m.getEntry(i,1),2));
		   }
		   
		   //s1=m[+,1];
		   //s2=m[+,2];
		   //s3=m[+,3];
		   //s4=m[+,4];
		   double[] s=colSums(m);
		   double s1=s[0];
		   double  s2=s[1];
		    double s3=s[2];
		    double  s4=s[3];
		 
		 
		     double probn=1.0/(ppois(n,n)-ppois(n-1,n));
		 
		   double z=((double) n-s1)/Math.sqrt(s2);
		   double  g1=s3/Math.pow(s2,1.5);
		  double  g2=s4/Math.pow(s2,2.0);
		 double  poly=1+g1*(Math.pow(z, 3)-3*z)/6.0+g2*(Math.pow(z, 4)-6*Math.pow(z, 2)+3)/24.0
				 +Math.pow(g1,2)*(Math.pow(z, 6)-15.0*Math.pow(z, 4)+45.0*Math.pow(z, 2)-15)/72.0;
//		   poly=1+g1*(z^3-3*z)/6+g2*(z^4-6*z^2+3)/24
		    //     +g1^2*(z^6-15*z^4+45*z^2-15)/72;
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
		
	}
	public void assignCount(int[] array) {
		this.x = new int[array.length];
		for(int j=0; j<array.length; j++){
			x[j] = (int) Math.round(array[j]);
		}
		salida = new Array2DRowRealMatrix(x.length,2);
		
	}
	public double[][] tab() {
		return salida.getData();
	}
}
