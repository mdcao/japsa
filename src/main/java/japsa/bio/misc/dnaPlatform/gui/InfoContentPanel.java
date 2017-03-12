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

package japsa.bio.misc.dnaPlatform.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.border.TitledBorder;

import japsa.bio.misc.dnaPlatform.sequence.*;

import javax.swing.border.Border;

/**
 * <p>
 * Title: InfoContentPanel
 * </p>
 * 
 * <p>
 * Description: This pabel contains a InfoContentPlot with a scrollbar and
 * buttons to manipulate the plot.
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * <p>
 * Company: Monash
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class InfoContentPanel extends JPanel {
	public static final long serialVersionUID = MainFrame.serialVersionUID;

	InfoContentPlot plot;
	JScrollBar plotScrollBar;

	Border border1 = BorderFactory
			.createLineBorder(new Color(100, 150, 200), 4);
	Border highlightBorder = new TitledBorder(border1, "");
	Border border = BorderFactory.createTitledBorder("");

	PlotPopupMenu graphPopup = new PlotPopupMenu(this);

	// tobe implemented later
	// PlotPopupMenu featurePopup= new PlotPopupMenu(this);

	// storing this reference to call method
	MainPanel myContainer;

	public InfoContentPanel(String name, MainPanel container) {
		setPlotTitle(name);
		myContainer = container;
		plot = new InfoContentPlot(myContainer);// Panel to draw the plots
		plotScrollBar = new JScrollBar(Scrollbar.HORIZONTAL, -20, 200, -20, 200);

		try {
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void setPlotTitle(String title) {
		border = BorderFactory.createTitledBorder(title);
		highlightBorder = new TitledBorder(border1, title);
	}

	private void jbInit() throws Exception {
		this.setLayout(new BorderLayout());
		plot.setBorder(border);

		plotScrollBar
				.addAdjustmentListener(new InfoContentPanel_plotScrollBar_adjustmentAdapter(
						this));

		this.addMouseMotionListener(new InfoContentPanel_this_mouseMotionAdapter(
				this));

		this.addMouseListener(new InfoContentPanel_this_mouseAdapter(this));

		add(plot, BorderLayout.CENTER);
		add(plotScrollBar, BorderLayout.SOUTH);
	}

	public void highlight(boolean b) {
		if (b)
			plot.setBorder(highlightBorder);
		else
			plot.setBorder(border);
	}

	/**
	 * This function returns the graph being displayed in this panel. This can
	 * be used to print the graph from main.
	 * 
	 * @return JPanel
	 */
	public JPanel getInfoContentGraph() {
		return plot;
	}

	/**
	 * Sets sequence of characters to be displayed in the plot under x
	 * coordinates.
	 * 
	 * @param sequence
	 *            char[]
	 */
	public void setSequence(char[] sequence) {
		plot.setSequence(sequence);// ,proteinConvertUtilities.convertToProtein(sequence));
	}

	/**
	 * This function adjusts the plotScrollBar according to the plot current
	 * minimum, maximum and limit x values
	 */
	private void adjustPlotScroll() {
		int extent = (int) (plot.getMaxXval() - plot.getMinXval());

		plotScrollBar.setValues((int) plot.getMinXval(), extent,
				(int) plot.getMinLimitX(), (int) plot.getMaxLimitX());
	}

	/**
	 * Is given a DoubleSequenceData object to display as a graph. First it gets
	 * or sets the name for DoubleSequenceData object. Then it adds this
	 * object's name to the graphBox and the object to the plot. Finally, the
	 * scroll bar of the plot is adjusted.
	 * 
	 * @param data
	 *            DoubleSequenceData
	 * @return boolean - indicating whether data was drawn in plot
	 */
	public boolean paintInfoContent(DoubleSequenceData data) {

		/* add graph to popup menu */
		if (graphPopup.addGraphToMenu(data.toString())) {
			/* displaying information content in plot */
			plot.addGraph(data);
			plot.selectGraph(data);
			adjustPlotScroll();

			return true;
		}
		return false;
	}

	public boolean addFeatures(AnnotationSequenceData data) {

		/* add graph to popup menu */

		plot.addFeatures(data);
		adjustPlotScroll();

		return true;
	}

	public boolean removeFeatures(AnnotationSequenceData data) {

		/* add graph to popup menu */

		plot.removeFeatures(data);
		adjustPlotScroll();

		return true;
	}

	/**
	 * Calls function zoomIn() of InfoContentPlot and adjusts the scrollbar for
	 * the information content plot.
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void zoomIn_actionPerformed() {
		plot.zoomIn();
		adjustPlotScroll();
	}

	/**
	 * Calls the function zoomOut() of InfoContentPlot and adjusts the scrollbar
	 * for the information content plot.
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void zoomOut_actionPerformed() {
		plot.zoomOut();
		adjustPlotScroll();
	}

	/**
	 * Whenever the value of the scrollbar for the plot is changed it calls
	 * function in InfoContentPlot to move display window of plot.
	 * 
	 * @param e
	 *            AdjustmentEvent
	 */
	public void plotScrollBar_adjustmentValueChanged(AdjustmentEvent e) {
		plot.moveXaxis(e.getValue());
	}

	/**
	 * Function to remove a graph from the graph box and plot
	 * 
	 * @param graphName
	 *            String
	 */
	public void removeGraph(String graphName) {
		if (graphName.equals(""))
			return;

		// remove graph from the plot and popup menu
		plot.removeGraph(graphName);
		graphPopup.removeGraphFromMenu(graphName);

		adjustPlotScroll();
	}

	/**
	 * When the gridRamgeButton is pressed then the GridRangeDialog is displayed
	 * for the user to change the grid range
	 * 
	 */
	public void showGridRangeDialog() {
		GridRangeDialog dlg = new GridRangeDialog(this.toString(),
				plot.getMinXval(), plot.getMaxXval(), plot.getMinYval(),
				plot.getMaxYval());
		Dimension dlgSize = dlg.getPreferredSize();
		Dimension frmSize = myContainer.getSize();
		java.awt.Point loc = getLocation();
		dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x,
				(frmSize.height - dlgSize.height) / 2 + loc.y);
		dlg.setModal(true);
		dlg.setVisible(true);

		/* update title if it has been changed */
		if (!dlg.getTitle().equals("")
				&& !dlg.getTitle().equals(this.toString())) {

			setPlotTitle(dlg.getTitle());
			plot.setBorder(highlightBorder);

			// if name changes then we have to update menu in MainFrame
			myContainer.updateInfoContentPlotMenu();
		}

		double yMin = dlg.getYminValue();
		double yMax = dlg.getYmaxValue();

		double xMin = dlg.getXminValue();
		double xMax = dlg.getXmaxValue();

		if (yMin < yMax)
			plot.setYrange(yMin, yMax);
		if (xMin < xMax)
			plot.setXrange(xMin, xMax);

		adjustPlotScroll();

	}

	/**
	 * When the user checks box to display coordinates then tell the plot to
	 * display mouse coordinates
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void showCoorBox_actionPerformed(ActionEvent e) {
		// plot.displayMouseCoord(showCoorBox.isSelected());
	}

	public String toString() {
		return ((TitledBorder) border).getTitle();
	}

	public void this_mouseMoved(MouseEvent e) {
		plot.mouseMoved(e.getX() - plot.getX(), e.getY() - plot.getY());
	}

	public void this_mousePressed(MouseEvent e) {
		myContainer.plotToolBar.selectInfoContentPanel(this);
		// Container.selectInformationContentPanel(this);Select info
	}

}

class InfoContentPanel_this_mouseAdapter extends MouseAdapter {
	private InfoContentPanel adaptee;

	InfoContentPanel_this_mouseAdapter(InfoContentPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void mousePressed(MouseEvent e) {
		adaptee.this_mousePressed(e);
		maybeShowPopup(e);
	}

	public void maybeShowPopup(MouseEvent e) {
		/* show popup menu if data in sequenceData */
		if (e.isPopupTrigger()) {
			adaptee.graphPopup.show(e.getComponent(), e.getX(), e.getY());
		}
	}
}

class InfoContentPanel_this_mouseMotionAdapter extends MouseMotionAdapter {
	private InfoContentPanel adaptee;

	InfoContentPanel_this_mouseMotionAdapter(InfoContentPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void mouseMoved(MouseEvent e) {
		adaptee.this_mouseMoved(e);
	}
}

class InfoContentPanel_plotScrollBar_adjustmentAdapter implements
		AdjustmentListener {
	private InfoContentPanel adaptee;

	InfoContentPanel_plotScrollBar_adjustmentAdapter(InfoContentPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void adjustmentValueChanged(AdjustmentEvent e) {
		adaptee.plotScrollBar_adjustmentValueChanged(e);
	}

}
