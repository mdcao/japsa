package japsadev.lib.jMEF;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Vector;

import javax.imageio.ImageIO;

public class Image {
	
	
	/**
	 * Converts an image into a point set of dimension 3 (RGB)
	 * @param   image  color image
	 * @return         point set
	 */
	public static PVector[] convertColorImageToPointSet3D(BufferedImage image){
		
		// Variables
		int       width     = image.getWidth();
		int       height    = image.getHeight();
		int       nb_pixels = width * height;
		PVector[] points    = new PVector[nb_pixels];

		// Conversion
		for (int row=0; row<height; row++)
			for (int col=0; col<width; col++){
				Color   c   = new Color(image.getRGB(col, row));
				PVector px  = new PVector(3);
				px.array[0] = c.getRed();
				px.array[1] = c.getGreen();
				px.array[2] = c.getBlue();	
				points[ row * width + col ] = px;
			}
		
		// Return
		return points;
	}
	
	
	/**
	 * Converts an image into a point set of dimension 5 (RGB + X + Y)
	 * @param   image  color image
	 * @return         point set
	 */
	public static PVector[] convertColorImageToPointSet5D(BufferedImage image){
		
		// Variables
		int       width     = image.getWidth();
		int       height    = image.getHeight();
		int       nb_pixels = width * height;
		PVector[] points    = new PVector[nb_pixels];

		// Conversion
		for (int row=0; row<height; row++)
			for (int col=0; col<width; col++){
				Color   c   = new Color(image.getRGB(col, row));
				PVector px  = new PVector(5);
				px.array[0] = c.getRed();
				px.array[1] = c.getGreen();
				px.array[2] = c.getBlue();
				px.array[3] = row;
				px.array[4] = col;			
				points[ row * width + col ] = px;
			}
		
		// Return
		return points;
	}
	
	
	/**
	 * Segment an image using a mixture model.
	 * @param  imgIn  input image
	 * @param  mm     mixture model
	 * @return        image segmentation
	 */
	public static BufferedImage segmentColorImageFromMOG(BufferedImage imgIn, MixtureModel mm){

		// Image initialization
		BufferedImage imgOut = new BufferedImage(imgIn.getWidth(), imgIn.getHeight(), BufferedImage.TYPE_INT_RGB);

		// Loop on pixels
		for (int row=0; row<imgIn.getHeight(); row++)
			for (int col=0; col<imgIn.getWidth(); col++){
							
				// Get the pixel
				Color c     = new Color(imgIn.getRGB(row, col));
				PVector px  = new PVector(3);
				px.array[0] = c.getRed();
				px.array[1] = c.getGreen();
				px.array[2] = c.getBlue();
				
				// Find and set the most probable class for the current point
				int    idx   = 0;
				double d_max = 0;
				for (int j=0; j<mm.size; j++){
					double d_tmp = mm.weight[j] * mm.EF.density(px, mm.param[j]);
					if (d_tmp>d_max){
						d_max = d_tmp;
						idx   = j;
					}
				}
				px = ((PVectorMatrix)mm.param[idx]).v;
				c  = new Color( (int)px.array[0], (int)px.array[1], (int)px.array[2] );
				imgOut.setRGB(row, col, c.getRGB());
			}
		
		// Return
		return imgOut;		
	}
	
	
	/**
	 * Reads an image.
	 * @param  imagePath image file to read
	 * @return           image
	 */
	public static BufferedImage readImage(String imagePath){
		BufferedImage image_in = null;
		try{
			image_in = ImageIO.read(new File(imagePath));
		}
		catch (IOException e) {
			    e.printStackTrace();
			    System.err.println("*** Error: Image file does not exist ***");
		}
		return image_in;
	}
	
	
	/**
	 * Writes an image.
	 * @param image      image to be written
	 * @param imagePath  image path
	 */
	public static void writeImage(BufferedImage image, String imagePath){
		try{
			ImageIO.write(image, "png", new File(imagePath));
		}
		catch (IOException e) {
			    e.printStackTrace();
		}
	}
	

	/**
	 * Computes the PSNR between two images.
	 * @param   i1  first image
	 * @param   i2  second image
	 * @return      PSNR(i1,i2)
	 */
	public static double PSNR(BufferedImage i1, BufferedImage i2){
		double mse  = 0;
		for (int r=0; r<i1.getHeight(); r++){
			for (int c=0; c<i1.getWidth(); c++){
				Color  c1 = new Color(i1.getRGB(c, r));
				Color  c2 = new Color(i2.getRGB(c, r));
				int    dr = c1.getRed()   - c2.getRed();
				int    dg = c1.getGreen() - c2.getGreen();
				int    db = c1.getBlue()  - c2.getBlue();
				mse      += dr*dr + dg*dg + db*db;
			}
		}
		mse /= (i1.getHeight() * i1.getWidth() * 3);
		return 10 * Math.log10( (255*255) / mse );
	}

	
	/**
	 * Load a mixture model from a file. If the mixture doesn't exist, the function create
	 * a mixture of multivariate Gaussians from the pixels of the image, and save this mixture.
	 * @param   path  file-path of the mixture model
	 * @param   image input image
	 * @param   d     mixture dimension, must be 3 or 5
	 * @param   n     number of components in the mixture model
	 * @return        a mixture of Gaussian of dimension d and of n components computed from the input image
	 */
	public static MixtureModel loadMixtureModel(String path, BufferedImage image, int d, int n){
		
		if (d!=3 && d!=5)
			throw new RuntimeException("Only dimension 3 and 5 are supported.");
		
		MixtureModel mm = MixtureModel.load(path);
		if (mm==null){
			PVector[]         px       = d==3 ? convertColorImageToPointSet3D(image) : convertColorImageToPointSet5D(image);
			Vector<PVector>[] clusters = KMeans.run(px, n);
			mm = BregmanSoftClustering.initialize(clusters, new MultivariateGaussian());
			mm = BregmanSoftClustering.run(px, mm);
			MixtureModel.save(mm, path);
		}
		else if (mm.getDimension()!=d) {
			throw new RuntimeException("Incorrect dimension.");
		}
		return mm;
	}
	
	
	/**
	 * Counts the minimum number of points assigned to image pixel.
	 * @param  tab    array
	 * @param  height height of the array
	 * @param  width  width of the array    
	 * @return        minimum number of points
	 */
	private static int min(int[][] tab, int height, int width){
		int min = Integer.MAX_VALUE;
		for (int y=0; y<height; y++)
			for (int x=0; x<width; x++)
				if (tab[y][x]<min)
					min=tab[y][x];
		return min;
	}

	
	/**
	 * Creates an image from a mixture of Gaussians.
	 * @param width  image width
	 * @param height image height
	 * @param mm     mixture model
	 * @return       statistical image
	 */
	public static BufferedImage createImageFromMixtureModel(int width, int height, MixtureModel mm){

		BufferedImage imgOut = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		double[][][] imgSum = new double[height][width][3];
		int[][]      imgCpt = new int[height][width];
		int n = 1024;
		int x, y, r, g, b;

		while (min(imgCpt, height, width)<10){
			
			// Draw point
			PVector[] pixels = mm.drawRandomPoints(n);
			
			// Fill the imgSum
			for (int i=0; i<n; i++){
				y = (int)pixels[i].array[3];
				x = (int)pixels[i].array[4];
				if (x>=0 && y>=0 && x<width && y<height){
					r = (int)pixels[i].array[0];
					g = (int)pixels[i].array[1];
					b = (int)pixels[i].array[2];
					if (r>=0 && g>=0 && b>=0  && r<255 && g<255 && b<255){
						imgSum[y][x][0] += r;
						imgSum[y][x][1] += g;
						imgSum[y][x][2] += b;
						imgCpt[y][x]++;
					}
				}
			}
		}
		
		// Normalize the colors
		for (y=0; y<height; y++){
			for (x=0; x<width; x++){
				r = (int)( imgSum[y][x][0] / imgCpt[y][x] );
				g = (int)( imgSum[y][x][1] / imgCpt[y][x] );
				b = (int)( imgSum[y][x][2] / imgCpt[y][x] );
				Color c = new Color( r, g, b );
				imgOut.setRGB( x, y, c.getRGB());
			}
		}
		
		// Return
		return imgOut;
	}
	

	/**
	 * Segment an image using a mixture model.
	 * @param  imgIn  input image
	 * @param  mm     mixture model
	 * @return        image segmentation
	 */
	public static BufferedImage createEllipseImage(int width, int height, MixtureModel f, double t){

		// Image initialization
		BufferedImage imgOut = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		for (int i=0; i<f.size; i++){

			PVector mean = (PVector)((PVectorMatrix)f.param[i]).v.clone();
			Color c      = new Color((int)mean.array[0], (int)mean.array[1], (int)mean.array[2]);

			for (int row=0; row<height; row++)
				for (int col=0; col<width; col++){
					
					PVector x  = (PVector)mean.clone();
					x.array[3] = row;
					x.array[4] = col;

					double v = (x.Minus(((PVectorMatrix)f.param[i]).v)).InnerProduct(((PVectorMatrix)f.param[i]).M.Inverse().MultiplyVectorRight(x.Minus(((PVectorMatrix)f.param[i]).v)));
					
					if (v<t)
						imgOut.setRGB(col, row, c.getRGB());
				}
		}

		// Return
		return imgOut;		
	}
}
