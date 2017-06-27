package japsadev.lib.jMEF;

import java.awt.image.BufferedImage;

import japsadev.lib.jMEF.Clustering.CLUSTERING_TYPE;

public class Tutorial3 {
	
	
	/**
	 * Main function.
	 * @param args
	 */
	public static void main(String[] args) {

		// Display
		String title = "";
		title += "+-----------------------------------------------+\n";
		title += "| Mixture simplification and image segmentation |\n";
		title += "+-----------------------------------------------+\n";
		System.out.print(title);
		
		// Variables
		int n = 32;
		int m = 8;
		
		// Image/texture information (to be changed to fit your configuration)
		String input_folder  = "/home/sonhoanghguyen/workspace/MixtureModels/Input/";
		String output_folder = "/home/sonhoanghguyen/workspace/MixtureModels/Output/";
		String image_name    = "Baboon";
		String image_path    = input_folder + image_name + ".png";
		String mixture_path  = String.format("%s%s_3D_%03d.mix", input_folder, image_name, n);
		
		// Read the input image
		System.out.print("Read input image                         : ");
		BufferedImage image = Image.readImage(image_path);
		System.out.println("ok");
		
		// Read or generate the mixture model
		System.out.print("Read/generate mixture model              : ");
		MixtureModel mm1 = Image.loadMixtureModel(mixture_path, image, 3, n);
		System.out.println("ok");
		
		// Compute the image segmentation based on the mixture mm1
		System.out.print("Segment image (mixture model)            : ");
		BufferedImage seg1 = Image.segmentColorImageFromMOG(image, mm1);
		Image.writeImage(seg1, String.format("%sTutorial3_%s_%03d.png", output_folder, image_name, n));
		System.out.println("ok");

		// Simplify mm1 in a mixture mm2 of m components and compute the image segmentation based on mm2
		System.out.print("Segment image (simplified mixture model) : ");
		MixtureModel  mm2  = BregmanHardClustering.simplify(mm1, m, CLUSTERING_TYPE.LEFT_SIDED);
		BufferedImage seg2 = Image.segmentColorImageFromMOG(image, mm2);
		Image.writeImage(seg2, String.format("%sTutorial3_%s_%03d.png", output_folder, image_name, m));
		System.out.println("ok");
	}

}
