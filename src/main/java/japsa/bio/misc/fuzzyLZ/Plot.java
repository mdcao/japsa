/*
 * Copyright (c) David Powell <david@drp.id.au>
 * 
 * 
 * This file is part of FuzzyLZ
 * 
 * FuzzyLZ is a program orginally intended for the compression of DNA sequeces.
 * It can be viewed as a compression model like Lempel-Ziv 77, but instead of
 * exact matches, allowing matches that contain inserts/deletes/mismatches.
 *  
 */

package japsa.bio.misc.fuzzyLZ;

import java.io.*;

import japsa.bio.misc.common.*;

class Plot implements Serializable {
	private static final long serialVersionUID = 1L;

	int img[][];

	double scale;

	int numRows, numCols;

	int startRow;

	Plot(int rows, int columns, int maxColumns, int maxRows) {
		this(rows, columns, maxColumns, maxRows, 0);
	}

	Plot(int rows, int columns, int maxColumns, int maxRows, int startRow) {
		this.startRow = startRow;

		scale = MyMath.min2(1.0 * MyMath.min2(columns, maxColumns) / columns,
				1.0 * MyMath.min2(rows - startRow, maxRows) / rows);

		numRows = (int) (scale * (rows - startRow));
		numCols = (int) (scale * columns);
		img = new int[numRows + 1][numCols + 1];

	}

	void put(int row, int col, double r, double g, double b) {
		row = (int) (scale * (row - startRow));
		col = (int) (scale * col);
		img[row][col] = ((byte) (r * 255)) << 16 | ((byte) (g * 255)) << 8
				| ((byte) (b * 255));
	}

	void putMax(int row, int col, double r, double g, double b) {
		row = (int) (scale * (row - startRow));
		col = (int) (scale * col);

		byte r1 = (byte) MyMath.max2((img[row][col] >> 16) & 255, r * 255);
		byte g1 = (byte) MyMath.max2((img[row][col] >> 8) & 255, g * 255);
		byte b1 = (byte) MyMath.max2((img[row][col]) & 255, b * 255);

		img[row][col] = (r1) << 16 | (g1) << 8 | (b1);
	}

	void save(String fname, String comments) {
		if (FuzzyLZ.DEBUG >= 1)
			System.out.println("Writing image '" + fname + "'");
		try {
			File f = new File(fname);
			DataOutputStream out = new DataOutputStream(
					new BufferedOutputStream(new FileOutputStream(f)));

			out.writeBytes("P6\n");
			out.writeBytes("#Created by FuzzyLZ\n");
			if (comments.length() > 0)
				out.writeBytes("#" + comments + "\n");
			out.writeBytes(numCols + " " + numRows + "\n");
			out.writeBytes("255\n");

			for (int r = 0; r < numRows; r++) {
				for (int c = 0; c < numCols; c++) {
					out.writeByte((img[r][c] >> 16) & 255);
					out.writeByte((img[r][c] >> 8) & 255);
					out.writeByte((img[r][c]) & 255);
				}
			}

			out.close();
		} catch (Exception e) {
			System.err.println("Error writing file: " + e);
		}
		if (FuzzyLZ.DEBUG >= 1)
			System.out.println("Done image");
	}

	public static void main(String args[]) {
		Plot p = new Plot(100, 100, 50, 50);
		for (int i = 0; i < 100; i++) {
			p.put(i, i, 0, 255, 0);
			p.put(i, 99 - i, 0, 0, 255);
			p.put(50, i, 255, 0, 0);
			p.put(i, 50, 255, 0, 0);
		}
		p.save("out.ppm", "");
	}
}