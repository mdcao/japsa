package japsadev.lib.jMEF;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Vector;

public class Test {
	/**
	 * Main function.
	 * @param args
	 */
	public static void main(String[] args) {

		// Display
		String title = "";
		title += "+----------------------------------------+\n";
		title += "| Testing soft clustering & classical EM |\n";
		title += "+----------------------------------------+\n";
		System.out.print(title);

		// Variables
		int n = 12;
	
	    File file = new File("/home/sonhoanghguyen/Projects/scaffolding/repeat/porecamp_metaSpades.hist");
		ArrayList<PVector> vectors = new ArrayList<PVector>();
	    try {

	        Scanner sc = new Scanner(file);
	        int count=0;
	        while (sc.hasNext()) {
	            int 	length = sc.nextInt(); 
	            double	cov = sc.nextDouble();
	            
	            for(int i=0; i < length; i++){
	            	PVector v = new PVector(1);
	            	v.array[0]=cov;
	            	vectors.add(v);
	            	
		            //System.out.println(count++ + ":" + length + " | " + cov);
	            }
	        }
	        sc.close();
	    } 
	    catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		
		PVector[] points = vectors.stream().toArray(PVector[]::new);
		System.out.println("Starting to fit by kmean...");
		long start = System.currentTimeMillis();

		NumberFormat formatter = new DecimalFormat("#0.00000");
		// Draw points from initial mixture model and compute the n clusters
		Vector<PVector>[] clusters = KMeans.run(points, n);

		long end1 = System.currentTimeMillis();
		System.out.print("running time: " + formatter.format((end1 - start) / 1000d) + " seconds");

		// Bregman soft clustering for Gauss
		MixtureModel mmef;
		mmef = BregmanSoftClustering.initialize(clusters, new UnivariateGaussian());
		mmef = BregmanSoftClustering.run(points, mmef);
		System.out.println("Mixure model of Gaussian estimated using Bregman soft clustering \n" + mmef + "\n");
		long end2 = System.currentTimeMillis();
		System.out.print("running time: " + formatter.format((end2 - end1) / 1000d) + " seconds");

		// Bregman soft clustering for Poisson
		MixtureModel mmp;
		mmp = BregmanSoftClustering.initialize(clusters, new Poisson());
		mmp = BregmanSoftClustering.run(points, mmp);
		System.out.println("Mixure model of Poisson estimated using Bregman soft clustering \n" + mmp + "\n");
		
		long end3 = System.currentTimeMillis();
		System.out.print("running time: " + formatter.format((end3 - end2) / 1000d) + " seconds");
	}
		
}
