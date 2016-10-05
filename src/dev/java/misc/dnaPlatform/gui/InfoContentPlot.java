/******************************************************************************
 * Copyright (C) 2006-2010 Minh Duc Cao                                        *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version. This program is distributed in the hope that it will be  *
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General    *
 * Public License for more details.                                            *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along with*
 * this program; if not, write to the Free Software  Foundation, Inc.,         *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                   *
 ******************************************************************************/

//This class is written by Julie Bernal and subsequently modified and maintained
//by Minh Duc Cao

package misc.dnaPlatform.gui;

import japsa.seq.JapsaFeature;

import javax.swing.*;

import misc.dnaPlatform.sequence.*;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import java.util.*;
import java.text.DecimalFormat;

/**
 * <p>
 * Title: InfoContentPlot
 * </p>
 * 
 * <p>
 * Description: This is a JPanel with fixed length used to plot long arrays of
 * doubles, the index of arrays are x coordinates and the values stored in
 * arrays are y coordinates. This class uses paint to calculate the section of
 * the plot that is visible in screen and display it.
 * </p>
 * 
 * 
 * @author Julie Bernal This class is modified by Hoang Anh Nguyen to add
 *         proteins and annotations
 * 
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public class InfoContentPlot extends JPanel {
	public static final long serialVersionUID = MainFrame.serialVersionUID;

	Grid plot;
	Yaxis yAxis;
	Xaxis xAxis;
	JPanel south = new JPanel();
	JPanel nothing = new JPanel();

	Vector<Graph> graphs = new Vector<Graph>();
	Graph selectedGraph = null;
	Vector<AnnotationSequenceData> featureList;
	JPanel mainPanel;

	/*
	 * string that represents the protein
	 */

	// char[] aminoAcidSequence=null;

	/* Just some colours for graphs */
	Color[] plotColors = { new Color(7, 61, 149), new Color(195, 98, 3),
			new Color(9, 114, 13), new Color(11, 79, 16),
			new Color(90, 8, 132), new Color(204, 31, 2) };

	Color[] featureColors = { Color.RED, Color.BLUE, Color.YELLOW, Color.GREEN };

	public InfoContentPlot(JPanel panel) {
		try {
			mainPanel = panel;
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception {
		/* initializing InfoContent objects */
		plot = new Grid(0, 4, 0.5, 0, 100, 1, 0, 200);
		yAxis = new Yaxis(0, 4, 0.5);
		xAxis = new Xaxis(0, 10, 1);

		setXrange(0, 200);

		/* arranging all the panels */
		setLayout(new BorderLayout());
		south.setLayout(new BorderLayout());

		nothing.setMaximumSize(new Dimension((int) yAxis.getMaximumSize()
				.getWidth(), (int) xAxis.getMaximumSize().getHeight()));
		nothing.setMinimumSize(new Dimension((int) yAxis.getMaximumSize()
				.getWidth(), (int) xAxis.getMinimumSize().getHeight()));
		nothing.setPreferredSize(new Dimension((int) yAxis.getPreferredSize()
				.getWidth(), (int) xAxis.getPreferredSize().getHeight()));

		nothing.setBackground(yAxis.getBackground());

		add(yAxis, BorderLayout.WEST);
		add(plot, BorderLayout.CENTER);
		add(south, BorderLayout.SOUTH);

		south.add(nothing, BorderLayout.WEST);
		south.add(xAxis, BorderLayout.CENTER);

		featureList = new Vector<AnnotationSequenceData>();
	}

	/**
	 * Calculates the coordinate in pixels of a value given the width and the
	 * range of x values displayed.
	 * 
	 * @param width
	 *            int, the length where values are displayed
	 * @param rangeValues
	 *            double, the number of values displayed
	 * @param value
	 *            double, value to calculate coordinate
	 * @return int, the coordinate of value given as parameter
	 */
	public int calculateXcoord(int width, double minXvalue, double maxXvalue,
			double value) {
		int xCoor = (int) ((value - minXvalue) / (maxXvalue - minXvalue) * width);
		if (xCoor < 0)
			return 0;

		else if (xCoor > width)
			return width;

		return xCoor;
	}

	/**
	 * Calculates value of given x coordinate given the width and the range of x
	 * values displayed.
	 * 
	 * @param width
	 *            int, the length where values are displayed
	 * @param minXvalue
	 *            double, minimum x value being displayed in plot
	 * @param maxXvalue
	 *            double, maximum x value being displayed in plot
	 * @param coordinate
	 *            int, coordinate we want to finc value for
	 * @return int, the value of coordinate given as a parameter
	 */
	public double calculateXvalue(int width, double minXvalue,
			double maxXvalue, int coordinate) {
		return (coordinate * (maxXvalue - minXvalue) / width) + minXvalue;
	}

	/**
	 * Calculates the y coordinate in pixels of a value given the height and the
	 * number of y values displayed.
	 * 
	 * @param height
	 *            int, the length where values are displayed
	 * @param minYvalue
	 *            double, minimum y value being displayed in plot
	 * @param maxYvalue
	 *            double, maximum y value being displayed in plot
	 * @param value
	 *            double, value to calculate coordinate
	 * @return int, the coordinate of value given as parameter
	 */
	public int calculateYcoord(int height, double minYvalue, double maxYvalue,
			double value) {
		int yCoor = height
				- (int) ((value - minYvalue) / (maxYvalue - minYvalue) * height);
		if (yCoor < 0)
			return 0;
		else if (yCoor > height)
			return height;

		return yCoor;
	}

	/**
	 * Calculates value of given y coordinate given the height and the range of
	 * y values displayed.
	 * 
	 * @param height
	 *            int, the length where values are displayed
	 * @param minYvalue
	 *            double, the minimum y value being displayed in plot
	 * @param maxYvalue
	 *            double, the maximum y value being displayed in plot
	 * @param coordinate
	 *            int, coordinate we want to finc value for
	 * @return int, the value of coordinate given as a parameter
	 */
	public double calculateYvalue(int height, double minYvalue,
			double maxYvalue, int coordinate) {
		return ((height - coordinate) * (maxYvalue - minYvalue) / height)
				+ minYvalue;
	}

	/**
	 * This function takes argument indicating whether or not to display mouse
	 * coordinates in plot
	 * 
	 * @param bool
	 *            boolean
	 */
	public void displayMouseCoord(boolean bool) {
		plot.showMouseCoordinates(bool);
		repaint();
	}

	/**
	 * Returns minimum x value that is displayed in the InfoContentPlot
	 * 
	 * @return double
	 */
	public double getMinXval() {
		return plot.getMinX();
	}

	/**
	 * Returns maximum x value that is displayed in the InfoContentPlot
	 * 
	 * @return double
	 */
	public double getMaxXval() {
		return plot.getMaxX();
	}

	/**
	 * Returns the minimum y value that is displayed in the InfoContentPlot
	 * 
	 * @return double
	 */
	public double getMinYval() {
		return plot.getMinY();
	}

	/**
	 * Returns the maximum y value that is diplayed in the InfoContentPlot
	 * 
	 * @return double
	 */
	public double getMaxYval() {
		return plot.getMaxY();
	}

	/**
	 * Changes the range of y values to be displated in the InfoContentPlot. The
	 * range of values for the Grid and the y axis have to be changed.
	 * 
	 * @param yMin
	 *            double
	 * @param yMax
	 *            double
	 */
	public void setYrange(double yMin, double yMax) {
		/* return if invalid y range of values */
		if (yMin > yMax)
			return;

		plot.setYrange(yMin, yMax);
		yAxis.setRange(yMin, yMax);

		/* fix the y scale */
		double divisions = (yMax - yMin) / 15;

		int i = 1;
		double scale = 0.1;
		while (scale < divisions) {
			if (i == 2) {
				i = 0;
				scale = (scale / 2) + (scale * 2);
			} else {
				i++;
				scale *= 2;
			}
		}

		plot.setYscale((double) scale);
		yAxis.setScale((double) scale);

		repaint();
	}

	/**
	 * Changes the range of x values to be displayed in the InfoContentPlot. The
	 * range of values for the Grid and the x axis have to be changed.
	 * 
	 * @param xMin
	 *            double
	 * @param xMax
	 *            double
	 */
	public void setXrange(double xMin, double xMax) {

		/* Return if new range is invalid */
		if (xMin > xMax || xMin < plot.getMinLimitX()) // || xMax >
														// plot.getMaxLimitX())
			return;

		plot.setXrange(xMin, xMax);
		xAxis.setRange(xMin, xMax);

		/* fix the y scale */
		double divisions = (xMax - xMin) / 15;

		int i = 1;
		double scale = 1;
		while (scale < divisions) {
			if (i == 2) {
				i = 0;
				scale = (scale / 2) + (scale * 2);
			} else {
				i++;
				scale *= 2;
			}
		}

		plot.setXscale(scale);
		xAxis.setScale(scale);

		repaint();
	}

	/**
	 * This function zooms in the plot. This is achieved by decreasing the range
	 * of values displayed in the plot and in the x axis.
	 */
	public void zoomIn() {
		double decrease = (plot.getMaxX() - plot.getMinX()) / 2;
		if (decrease > 1)
			setXrange(plot.getMinX(), plot.getMaxX() - decrease);
	}

	/**
	 * This function zooms out of the plot. This is achieved by increasing the
	 * range of values displayed in the plot and x axis.
	 */
	public void zoomOut() {
		double increase = (plot.getMaxX() - plot.getMinX()) / 2;
		setXrange(plot.getMinX(), plot.getMaxX() + increase);
	}

	/**
	 * Returns minimum limit x value for the plot
	 * 
	 * @return double
	 */
	public double getMinLimitX() {
		return plot.getMinLimitX();
	}

	/**
	 * Returns maximum limit x value for the plot
	 * 
	 * @return double
	 */
	public double getMaxLimitX() {
		return plot.getMaxLimitX();
	}

	/**
	 * Moves the plot to the value indicated
	 * 
	 * @param xValue
	 *            double
	 */
	public void moveXaxis(double value) {
		double xRange = plot.getMaxX() - plot.getMinX();
		setXrange(value, value + xRange);
	}

	/**
	 * This function is used to set character sequence of the plot It also set
	 * the aminoAcid array
	 * 
	 * public void setSequence (char[] japsa.seq, char[] aminoSeq) {
	 * this.aminoAcidSequence = aminoSeq; int seqLen =
	 * xAxis.setSequence(japsa.seq);
	 * 
	 * 
	 * if(plot.getMaxLimitX() < seqLen) plot.setMaxLimitX(seqLen);
	 * 
	 * repaint(); }
	 */

	/**
	 * This function is used to set character sequence of the plot
	 */
	public void setSequence(char[] seq) {
		int seqLen = xAxis.setSequence(seq);

		/* Change the x limit in the plot */
		if (plot.getMaxLimitX() < seqLen)
			plot.setMaxLimitX(seqLen);

		repaint();
	}

	/**
	 * Calculates maximum x value in the plot to be displayed. This value is the
	 * maxixmum lenght of all graphs and char sequence of x axis
	 */
	private int calculateMaxXvalue() {
		Iterator it = graphs.iterator();
		int maxLen = 0;

		while (it.hasNext()) {
			Graph gr = (Graph) it.next();
			SequenceData d = gr.getData();
			if (d.getData().length > maxLen)
				maxLen = d.getData().length;
		}

		for (int lidx = 0; lidx < featureList.size(); lidx++) {
			AnnotationSequenceData featureData = featureList.get(lidx);
			if (featureData.size() > 0) {
				int x = featureData.getFeature(featureData.size() - 1).getEnd();
				if (x > maxLen)
					maxLen = x;
			}
		}

		if (xAxis.getSequenceLenght() > maxLen)
			maxLen = xAxis.getSequenceLenght();

		return maxLen;

	}

	/**
	 * creates a graph with given DoubleSequenceData object and returns the
	 * index of this graph in vector graphs.
	 * 
	 * @param data
	 *            SequenceData
	 * @return int
	 */
	public void addGraph(DoubleSequenceData data) {
		Graph gr = new Graph(plotColors[graphs.size() % plotColors.length],
				data);
		graphs.add(gr);

		selectedGraph = gr;

		if (plot.getMaxLimitX() < data.getData().length)
			plot.setMaxLimitX(data.getData().length);

		repaint();
	}

	/**
	 * Res the graph in graphs at given index and repaints the plot
	 * 
	 * @param index
	 *            int
	 * @return Graph
	 */

	public boolean removeFeatures(AnnotationSequenceData seqData) {
		if (featureList.remove(seqData)) {

			// calculate new maximum x value
			plot.setMaxLimitX(calculateMaxXvalue());

			repaint();
			return true;
		}
		return false;
	}

	public boolean addFeatures(AnnotationSequenceData seqData) {
		if (featureList.contains(seqData)) {
			System.out.println("The features already added");
			return false;
		}
		if (featureList.size() >= featureColors.length) {
			System.out.println("No more space");
			return false;
		}

		featureList.add(seqData);

		// calculate new maximum x value
		plot.setMaxLimitX(calculateMaxXvalue());
		repaint();

		return true;
	}

	/* TODO: Eliminate string comparisons with SequenceData names and types */

	public void removeGraph(String name) {
		Iterator it = graphs.iterator();
		while (it.hasNext()) {
			Graph gr = (Graph) it.next();
			String grName = (gr.getData()).toString();
			if (name.equals(grName)) {
				graphs.removeElement(gr);
				if (selectedGraph == gr)
					selectedGraph = null;
				break;
			}
		}

		// calculate new maximum x value
		plot.setMaxLimitX(calculateMaxXvalue());

		repaint();
	}

	/**
	 * Sets the selected graph to be the graph corresponding to the
	 * DoubleSequenceData object with the given name.
	 * 
	 * @param readID
	 *            String
	 */
	public void selectGraph(DoubleSequenceData data) {
		Iterator it = graphs.iterator();
		while (it.hasNext()) {
			Graph gr = (Graph) it.next();
			if (data == gr.getData())
				selectedGraph = gr;
		}

		repaint();
	}

	/**
	 * 
	 * <p>
	 * Title: Grid
	 * </p>
	 * 
	 * <p>
	 * Description: This is the grid where graphs get drawn in the
	 * InfoContentPlot
	 * </p>
	 * 
	 * 
	 * @author Julie Bernal
	 * @version 1.0
	 */
	private class Grid extends JPanel {
		// public static final long serialVersionUID =
		// MainFrame.serialVersionUID;

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		/* holds minimum and maximum values to be displayed */
		double minY, maxY;
		double minX, maxX;

		double scaleY, scaleX;

		/* width and height of Plot */
		int width, height;

		/* to store mouse coordinates */
		int mouseX, mouseY;
		boolean showMouseCoord = false;

		/* holds limits for x values */
		double minLimitX, maxLimitX;

		public Grid(double minYvalue, double maxYvalue, double yScale,
				double minXvalue, double maxXvalue, double xScale,
				double minXlimit, double maxXlimit) {
			super();
			if (minYvalue <= maxYvalue) {
				minY = minYvalue;
				maxY = maxYvalue;
			}

			if (minXvalue <= maxXvalue) {
				minX = minXvalue;
				maxX = maxXvalue;
			}

			scaleY = yScale;
			scaleX = xScale;

			minLimitX = minXlimit;
			maxLimitX = maxXlimit;

			setBackground(Color.white);
			setBorder(BorderFactory.createLineBorder(Color.black));
		}

		public void paint(Graphics g) {
			width = getWidth();
			height = getHeight();

			g.setColor(getBackground());
			g.fillRect(0, 0, width, height);

			paintGrid(g);

			// painting graphs in the grid
			Iterator<Graph> i = graphs.iterator();
			while (i.hasNext()) {
				Graph gr = i.next();
				g.setColor(gr.getColor());
				paintGraph(g, gr.getData(), false);
			}

			// painting selected graph
			if (selectedGraph != null) {
				g.setColor(selectedGraph.getColor());
				paintGraph(g, selectedGraph.getData(), false);// true);
			}

			if (showMouseCoord)
				paintMouseCoordinates(g);

		}

		/**
		 * Paints a grid
		 * 
		 * @param g
		 *            Graphics
		 */
		private void paintGrid(Graphics g) {
			int x0 = calculateXcoord(width, minX, maxX, 0);
			int y0 = calculateYcoord(height, minY, maxY, 0);

			if (x0 < 0 || x0 > width)
				x0 = 0;

			if (y0 < 0 || y0 > height)
				y0 = 0;

			g.setColor(Color.lightGray);

			/* drawing grid from axis */
			for (double y = minY; y < maxY; y += scaleY) {
				int yCoor = calculateYcoord(height, minY, maxY, y);
				g.drawLine(0, yCoor, width, yCoor);
			}

			/*
			 * to have gray lines at x values multiples of scaleX, calculate
			 * xValue at coordinate 0 and add to coordinate until corresponding
			 * value is a multiple of scaleX
			 */

			double initialXval = minX;
			while (initialXval % scaleX > 0)
				initialXval++;

			for (double x = initialXval; x < maxX; x += scaleX) {
				int xCoor = calculateXcoord(width, minX, maxX, x);
				g.drawLine(xCoor, 0, xCoor, height);
			}

			/* draw axis */
			g.setColor(Color.black);
			if (x0 > 0 && x0 < width)
				g.drawLine(x0, 0, x0, height);

			if (y0 > 0 && y0 < height)
				g.drawLine(0, y0, width, y0);

		}

		/**
		 * Paints mouse coordinates
		 * 
		 * @param g
		 *            Graphics
		 */
		private void paintMouseCoordinates(Graphics g) {

			g.setColor(Color.black);

			int mX = mouseX, mY = mouseY;
			if (mouseX + 50 > width)
				mX = mouseX - 50;
			if (mouseY - 50 < 0)
				mY = mouseY + 50;

			double x = calculateXvalue(width, minX, maxX, mouseX);
			double y = calculateYvalue(height, minY, maxY, mouseY);

			String xVal = String.valueOf(x);
			if (xVal.indexOf('.') > -1)
				xVal = xVal.substring(0, xVal.indexOf('.'));

			String yVal = String.valueOf(y);
			if (yVal.indexOf('.') > -1 && yVal.indexOf('.') < yVal.length() - 4)
				yVal = yVal.substring(0, yVal.indexOf('.') + 3);

			// g.drawString("(" + mouseX + "," + mouseY + ")", mX, mY);
			g.drawString("(" + xVal + "," + yVal + ")", mX, mY - 20);
			// g.drawString("coords (" + xCoor + "," + yCoor + ")", mX, mY -
			// 40);
		}

		/**
		 * Sets whether or not mouse coordinates should be displayed in the plot
		 * 
		 * @param bool
		 *            boolean
		 */
		public void showMouseCoordinates(boolean bool) {
			showMouseCoord = bool;
		}

		/**
		 * Paints a DoubleSequenceData double array. Indexes in array are
		 * treated as x coordinates and values stored in array are treated as y
		 * coordinates. Only points visible in the panel are drawn. Graphs are
		 * drawn in the plot from point 1 while in array points start [0] so
		 * when x coordinates are calculated, they are calculated as x+1 for any
		 * x value.
		 * 
		 * @param g
		 *            Graphics
		 * @param data
		 *            SequenceData
		 */

		public void paintGraph(Graphics g, DoubleSequenceData data,
				boolean highlight) {
			double infoContent[] = data.getDoubleData();

			// minYval and maxYval are used to draw only one line at each
			// pixel even if there are many points to be plotted at that pixel
			// this is done by drawing a vertical line from the lowest y value
			// to
			// the highest y value
			double minYval = Double.POSITIVE_INFINITY;
			double maxYval = Double.NEGATIVE_INFINITY;

			int xCoor, lastXcoor = -1;

			int x = (int) minX;

			/* draw first line if minX > 0 */
			if (x > 0 && x < infoContent.length) {
				int xC1 = calculateXcoord(width, minX, maxX, x + 1);
				int xC2 = calculateXcoord(width, minX, maxX, x);
				int yC1 = calculateYcoord(height, minY, maxY, infoContent[x]);
				int yC2 = calculateYcoord(height, minY, maxY,
						infoContent[x - 1]);

				if (highlight) {
					int xPoints[] = { xC1 - 1, xC1 + 1, xC2 + 1, xC2 - 1 };
					int yPoints[] = { yC1 - 1, yC1 + 1, yC2 + 1, yC2 - 1 };

					g.fillPolygon(xPoints, yPoints, 4);
				}

				else
					g.drawLine(xC1, yC1, xC2, yC2);
			}

			lastXcoor = calculateXcoord(width, minX, maxX, x + 1);
			x++;
			/*
			 * draw graph lines while x is still visible in the screen and a
			 * point in infoContent[]
			 */
			while (x < maxX && x < infoContent.length) {
				xCoor = calculateXcoord(width, minX, maxX, x + 1);

				// Only draw one line at each pixel
				if (xCoor == lastXcoor) {
					if (infoContent[x] < minYval)
						minYval = infoContent[x];
					if (infoContent[x] > maxYval)
						maxYval = infoContent[x];
				} else {
					// If minYval and maxYval have been set draw a line at
					// lastXcoor
					if (minYval < maxYval) {
						g.drawLine(lastXcoor,
								calculateYcoord(height, minY, maxY, minYval),
								lastXcoor,
								calculateYcoord(height, minY, maxY, maxYval));

						minYval = Double.POSITIVE_INFINITY;
						maxYval = Double.NEGATIVE_INFINITY;
					}

					// drawing lines
					int lastYcoor = calculateYcoord(height, minY, maxY,
							infoContent[x - 1]);
					int yCoor = calculateYcoord(height, minY, maxY,
							infoContent[x]);
					if (highlight) {
						int xPoints[] = { lastXcoor - 1, lastXcoor + 1,
								xCoor + 1, xCoor - 1 };
						int yPoints[] = { lastYcoor - 1, lastYcoor + 1,
								yCoor + 1, yCoor - 1 };

						g.fillPolygon(xPoints, yPoints, 4);
					} else
						g.drawLine(lastXcoor, lastYcoor, xCoor, yCoor);
				}

				lastXcoor = xCoor;
				x++;
			}
		}

		/**
		 * returns minimum y value displayed in the plot
		 * 
		 * @return double
		 */
		public double getMinY() {
			return minY;
		}

		/**
		 * returns maximum y value displayed in the plot
		 * 
		 * @return double
		 */
		public double getMaxY() {
			return maxY;
		}

		/**
		 * Returns minimum x value displayed in the plot
		 * 
		 * @return double
		 */
		public double getMinX() {
			return minX;
		}

		/**
		 * Returns maximum x value displayed in the plot
		 * 
		 * @return double
		 */
		public double getMaxX() {
			return maxX;
		}

		/**
		 * Changes the range of y values to display in the plot
		 * 
		 * @param yMin
		 *            double, minimum y value to display
		 * @param yMax
		 *            double, maximum y value to display
		 */
		public void setYrange(double yMin, double yMax) {
			if (minY <= maxY) {
				minY = yMin;
				maxY = yMax;
			}
		}

		/**
		 * Changes the range of x values to display in the plot
		 * 
		 * @param xMin
		 *            double, minimum x value to display
		 * @param xMax
		 *            double, maximum x value to display
		 */
		public void setXrange(double xMin, double xMax) {
			if (minX <= maxX) {
				minX = xMin;
				maxX = xMax;
			}
		}

		/**
		 * Returns the x scale of the plot. The x scale is the increase value to
		 * draw lines in the grid and the x axis
		 * 
		 * @return double
		 */
		@SuppressWarnings("unused")
		public double getXscale() {
			return scaleX;
		}

		/**
		 * Returns the y scale of the plot. The y scale is the increase value to
		 * draw lines in the grid and the y axis
		 * 
		 * @return double
		 */
		@SuppressWarnings("unused")
		public double getYscale() {
			return scaleY;
		}

		/**
		 * Sets the x scale of the plot. The x scale is the increase value to
		 * draw lines in the grid and the x axis
		 * 
		 * @param scale
		 *            double
		 */
		public void setXscale(double scale) {
			scaleX = scale;
		}

		/**
		 * Sets the y scale of the plot. The y scale is the increase value to
		 * draw lines in the grid and the y axis
		 * 
		 * @param scale
		 *            double
		 */
		public void setYscale(double scale) {
			scaleY = scale;
		}

		/**
		 * Returns the minimum allowed value in x axis
		 * 
		 * @return double
		 */
		public double getMinLimitX() {
			return minLimitX;
		}

		/**
		 * Returns the maximum value allowed in x axis
		 */
		public double getMaxLimitX() {
			return maxLimitX;
		}

		/**
		 * Sets the maximum value allowed in x axis
		 * 
		 * @param maxLimit
		 *            double
		 */
		public void setMaxLimitX(double maxLimit) {
			maxLimitX = maxLimit;
		}

	}

	public void mouseMoved(int x, int y) {
		if (plot.showMouseCoord) {
			plot.mouseX = x - (plot.getX()) + 2;
			plot.mouseY = y - (plot.getY()) - 3;

			repaint();
		}
	}

	/**
	 * 
	 * <p>
	 * Title: Yaxis
	 * </p>
	 * 
	 * <p>
	 * Description: A Yaxis is a JPanel of width 30 that displays Y values of an
	 * InfoContentPlot
	 * </p>
	 * 
	 * @author Julie Bernal
	 * @version 1.0
	 */
	private class Yaxis extends JPanel {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		double minY, maxY, scaleY;

		int width, height;

		public Yaxis(double minYvalue, double maxYvalue, double yScale) {
			super();
			if (minYvalue <= maxYvalue) {
				minY = minYvalue;
				maxY = maxYvalue;
			}
			scaleY = yScale;
			setBackground(Color.lightGray);

			setMaximumSize(new Dimension(30, 32767));
			setMinimumSize(new Dimension(30,
					(int) ((maxY - minY) / yScale) * 10));
			setPreferredSize(new Dimension(30,
					(int) ((maxY - minY) / yScale) * 10));

		}

		/**
		 * Writes y values in y axis. Values written depend in the minimum and
		 * maximum y value to display and the y scale
		 * 
		 * @param g
		 *            Graphics
		 */
		public void paint(Graphics g) {
			height = getHeight();
			width = getWidth();

			g.setColor(getBackground());
			g.fillRect(0, 0, width, height);

			g.setColor(Color.black);

			/* painting y values */
			Font oldFont = g.getFont();
			g.setFont(new Font("Dialog", Font.PLAIN, 9));

			DecimalFormat df = new DecimalFormat("#0.0");
			// decimal format for big numbers:
			// decimal format for big numbers:
			if (scaleY >= 1000) {
				df = new DecimalFormat("#0.#####E0");
				g.setFont(new Font("Dialog", Font.PLAIN, 9));
			} else if (scaleY >= 1)
				df = new DecimalFormat("#0.#");

			for (double y = minY; y < maxY; y += scaleY) {
				int i = calculateYcoord(height, minY, maxY, y);
				// only draw numbers if they fit
				int x = width - g.getFontMetrics().stringWidth(df.format(y));

				g.drawString(df.format(y), x, i);
			}

			g.setFont(oldFont);

		}

		/**
		 * Changes the range of values being displayed in y axis
		 * 
		 * @param yMinValue
		 *            double, minimum y value to display
		 * @param yMaxValue
		 *            double, maximum y value to display
		 */
		public void setRange(double yMinValue, double yMaxValue) {
			if (yMinValue <= yMaxValue) {
				minY = yMinValue;
				maxY = yMaxValue;
			}
		}

		/**
		 * Sets the scale for the y axis
		 * 
		 * @param scale
		 *            double
		 */
		public void setScale(double scale) {
			scaleY = scale;
		}

	}

	/**
	 * 
	 * <p>
	 * Title: Xaxis
	 * </p>
	 * 
	 * <p>
	 * Description: A Xaxis is a class of height 30 that displays annotation of
	 * an InfoContentPlot
	 * </p>
	 * 
	 * @author Julie Bernal
	 * @version 1.0 Modified by Minh Duc Cao
	 */
	private class Xaxis extends JPanel {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		static final int START = 25, HEIGHT = 8, DISTANCE = 12;

		double minX, maxX, scaleX;
		int width, height;

		private char[] sequence;

		public Xaxis(double minXvalue, double maxXvalue, double xScale) {
			super();
			if (minXvalue <= maxXvalue) {
				minX = minXvalue;
				maxX = maxXvalue;
			}
			scaleX = xScale;
			setBackground(Color.lightGray);
			// the height of the Maxis is increase from 30 to 55
			// in order to have enough place to put 3 more strings: the protein
			// sequence.
			setMaximumSize(new Dimension(32767, 85));
			setMinimumSize(new Dimension(10, 85));
			setPreferredSize(new Dimension(10, 85));

			// add mouse listener
			// If a features is clicked, display the feature property
			addMouseListener(new MouseAdapter() {
				public void mousePressed(MouseEvent e) {
					int index = -1;
					int x_mouse = e.getX();
					int y_mouse = e.getY();
					for (int lidx = 0; lidx < featureList.size(); lidx++) {
						if (START + DISTANCE * lidx <= y_mouse
								&& START + DISTANCE * lidx >= y_mouse - HEIGHT) {
							index = lidx;
							break;
						}
					}

					if (index >= 0) {
						AnnotationSequenceData featureData = featureList
								.get(index);

						for (int i = 0; i < featureData.size(); i++) {
							JapsaFeature feature = featureData.getFeature(i);
							if (feature.getEnd() > minX
									&& feature.getStart() < maxX) {// overlap
								int begin_drawing = calculateXcoord(width,
										minX, maxX, feature.getStart());
								// width of the rectangle
								int length_drawing = calculateXcoord(width,
										minX, maxX, feature.getEnd())
										- begin_drawing;
								if (begin_drawing <= x_mouse
										&& begin_drawing + length_drawing >= x_mouse) {
									// pop up
									JOptionPane.showMessageDialog(mainPanel,
											feature.getProperty(),
											"Feature info",
											JOptionPane.PLAIN_MESSAGE);
								}
							}
						}
					}
				}
			});

		}

		/**
		 * Writes the values in x axis, these values depend in the minimum and
		 * maximum x values to display and the scale of x axis
		 * 
		 * @param g
		 *            Graphics
		 */
		public void paint(Graphics g) {
			height = getHeight();
			width = getWidth();

			g.setColor(getBackground());
			g.fillRect(0, 0, width, height);

			g.setColor(Color.black);

			/* painting x values */
			int pixelSpace = (int) (width / ((maxX - minX) / scaleX));
			Font oldFont = g.getFont();
			g.setFont(new Font("Dialog", Font.PLAIN, 9));

			DecimalFormat df = new DecimalFormat("#0.0");
			// decimal format for big numbers:
			if (scaleX >= 10000)
				df = new DecimalFormat("#0.#####E0");
			else if (scaleX >= 1)
				df = new DecimalFormat("#0.#");

			/*
			 * to have values at x values multiples of scaleX, calculate xValue
			 * at coordinate 0 and add to coordinate until corresponding value
			 * is a multiple of scaleX
			 */

			double initialXval = minX;
			while (initialXval % scaleX > 0)
				initialXval++;

			for (double x = initialXval; x < maxX; x += scaleX) {
				int i = calculateXcoord(width, minX, maxX, x);
				// only draw numbers if they fit
				if (g.getFontMetrics().stringWidth(df.format(x)) < (pixelSpace - 5)) {
					if (x == 0)
						g.drawString("0", i, 10);

					else {
						i -= g.getFontMetrics().stringWidth(df.format(x)) / 2;
						g.drawString(df.format(x), i, 10);
					}
				}
			}

			g.setFont(oldFont);

			/*
			 * draw sequence if it exists and the range of x values displayed in
			 * the screen is small enough to display characters in string
			 * sequence[0] should be painted as sequence[1]
			 */
			// changes made by Hoang Nguyen

			int pixelsUnit = (int) (width / (maxX - minX));
			if (sequence != null && pixelsUnit > 5) {
				for (int x = (int) minX; x < maxX && x < sequence.length; x++)
					g.drawString("" + sequence[x],
							calculateXcoord(width, minX, maxX, x + 1), 20);
			}// end if pixel unit < 5

			for (int lidx = 0; lidx < featureList.size(); lidx++) {
				AnnotationSequenceData featureData = featureList.get(lidx);

				for (int i = 0; i < featureData.size(); i++) {
					JapsaFeature feature = featureData.getFeature(i);
					if (feature.getEnd() > minX && feature.getStart() < maxX) {// overlap
						g.setColor(featureColors[lidx]);
						int begin_drawing = calculateXcoord(width, minX, maxX,
								feature.getStart());
						// System.out.println("=====" + width + "  " + minX +
						// "  " + maxX + " " + feature.getStart() + "  " +
						// begin_drawing);
						// width of the rectangle
						int length_drawing = calculateXcoord(width, minX, maxX,
								feature.getEnd()) - begin_drawing;
						g.drawRect(begin_drawing, START + DISTANCE * lidx,
								length_drawing, HEIGHT);
					}
				}

			}

		}

		/**
		 * Changes the range of x values being displayed in x axis
		 * 
		 * @param xMinValue
		 *            double, minimum x value to display
		 * @param xMaxValue
		 *            double, maximum x value to display
		 */
		public void setRange(double xMinValue, double xMaxValue) {
			if (xMinValue <= xMaxValue) {
				minX = xMinValue;
				maxX = xMaxValue;
			}
		}

		/**
		 * Sets the scale for the x axis
		 * 
		 * @param scale
		 *            double
		 */
		public void setScale(double scale) {
			scaleX = scale;
		}

		/**
		 * This funciton is used to give a character sequence to be displayed as
		 * x coordinate
		 * 
		 * @param japsa
		 *            .seq char[]
		 * @return int the lenght of the sequence added to x axis
		 */
		public int setSequence(char[] seq) {
			sequence = seq;
			return seq.length;
		}

		/**
		 * Returns the lenght of character sequence in x aXis or 0 if the
		 * sequence has not been set
		 * 
		 * @return int
		 */
		public int getSequenceLenght() {
			if (sequence == null)
				return 0;
			return sequence.length;
		}
	}

	/**
	 * 
	 * <p>
	 * Title: Graph
	 * </p>
	 * 
	 * <p>
	 * Description: A graph contains SequenceData and other information about
	 * how it is displayed. For example the colour of the graph
	 * </p>
	 */
	private class Graph {
		private Color myColor;
		DoubleSequenceData myData;

		public Graph(Color color, DoubleSequenceData data) {
			myColor = color;
			myData = data;
		}

		public Color getColor() {
			return myColor;
		}

		public DoubleSequenceData getData() {
			return myData;
		}
	}

}
