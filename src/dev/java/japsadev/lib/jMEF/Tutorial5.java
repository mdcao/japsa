package japsadev.lib.jMEF;

import java.awt.image.BufferedImage;

public class Tutorial5{
	
	
	/**
	 * Main function.
	 * @param args
	 */
	public static void main(String[] args) {
		
		// Display
		String title = "";
		title += "+----------------------------------------+\n";
		title += "| Statistical images from mixture models |\n";
		title += "+----------------------------------------+\n";
		System.out.print(title);

		// Variables
		int n = 32;
		
		// Image/texture information (to be changed to fit your configuration)
		String input_folder  = "/home/sonhoanghguyen/workspace/MixtureModels/Input/";
		String output_folder = "/home/sonhoanghguyen/workspace/MixtureModels/Output/";
		String image_name    = "Baboon";
		String image_path    = input_folder + image_name + ".png";
		String mixture_path  = String.format("%s%s_5D_%03d.mix", input_folder, image_name, n); 
		
		// Read the input image
		System.out.print("Read input image             : ");
		BufferedImage image = Image.readImage(image_path);
		System.out.println("ok");
		
		// Read or generate the mixture model
		System.out.print("Read/generate mixture model  : ");
		MixtureModel f = Image.loadMixtureModel(mixture_path, image, 5, n);
		System.out.println("ok");
		
		// Creates and save the statistical image
		System.out.print("Create statistical image     : ");
		BufferedImage stat = Image.createImageFromMixtureModel(image.getWidth(), image.getHeight(), f);
		Image.writeImage(stat, String.format("%sTutorial5_%s_statistical_%03d.png", output_folder, image_name, n));
		System.out.println("ok");

		// Creates and save the ellipse image
		System.out.print("Create ellipse image         : ");
		BufferedImage ell = Image.createEllipseImage(image.getWidth(), image.getHeight(), f, 2);
		Image.writeImage(ell, String.format("%sTutorial5_%s_ellipses_%03d.png", output_folder, image_name, n));
		System.out.println("ok");
		
	}
	
}
