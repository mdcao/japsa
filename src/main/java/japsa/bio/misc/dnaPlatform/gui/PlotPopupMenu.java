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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * <p>
 * Title: DnaPlatform
 * </p>
 * 
 * <p>
 * Description: This is a java platform for visualization of compression results
 * obtained with different models
 * </p>
 * 
 * <p>
 * Copyright: Copyright (c) 2005
 * </p>
 * 
 * <p>
 * Company: Monash University
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class PlotPopupMenu extends JPopupMenu {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	InfoContentPanel plot;

	/* GUI elements for InfoContentPlot menu */
	JMenuItem zoomIn = new JMenuItem("zoom in");
	JMenuItem zoomOut = new JMenuItem("zoom out");
	JMenuItem gridRange = new JMenuItem("plot view");
	JMenu graph = new JMenu("graphs");
	JMenu graph_remove = new JMenu("remove");

	public PlotPopupMenu(InfoContentPanel plot) {
		this.plot = plot;

		try {
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception {
		gridRange.addActionListener(new MyActionAdapter(this));
		zoomIn.addActionListener(new MyActionAdapter(this));
		zoomOut.addActionListener(new MyActionAdapter(this));

		/* setting viewMenu */
		graph.add(graph_remove);
		add(gridRange);
		add(zoomIn);
		add(zoomOut);
		add(graph);

	}

	public void gridRange_actionPerformed(ActionEvent e) {
		plot.showGridRangeDialog();
	}

	public void zoomIn_actionPerformed(ActionEvent e) {
		plot.zoomIn_actionPerformed();
	}

	public void zoomOut_actionPerformed(ActionEvent e) {
		plot.zoomOut_actionPerformed();
	}

	/**
	 * This method adds a graph name to the popup menu if this name is not
	 * already in menu. It returns a boolean indicating whether or not new graph
	 * was added.
	 * 
	 * @param data
	 *            SequenceData
	 * @return String
	 */
	public boolean addGraphToMenu(String graphName) {

		for (int i = 0; i < graph_remove.getItemCount(); i++) {
			if ((graph_remove.getItem(i).getText()).equals(graphName)) {
				return false;
			}
		}

		JMenuItem item = new JMenuItem(graphName);
		item.addActionListener(new MyRemoveActionAdapter(this, item));
		graph_remove.add(item);
		return true;
	}

	/**
	 * This method removes given graph name from the menu
	 * 
	 * @param graphName
	 *            String
	 */
	public void removeGraphFromMenu(String graphName) {
		for (int i = 0; i < graph_remove.getItemCount(); i++) {
			if ((graph_remove.getItem(i).getText()).equals(graphName)) {
				graph_remove.remove(i);
			}
		}
	}

}

class MyActionAdapter implements ActionListener {
	private PlotPopupMenu adaptee;

	MyActionAdapter(PlotPopupMenu adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource().equals(adaptee.gridRange))
			adaptee.gridRange_actionPerformed(e);

		else if (e.getSource().equals(adaptee.zoomIn))
			adaptee.zoomIn_actionPerformed(e);

		else if (e.getSource().equals(adaptee.zoomOut))
			adaptee.zoomOut_actionPerformed(e);

	}
}

class MyRemoveActionAdapter implements ActionListener {
	private PlotPopupMenu adaptee;
	private JMenuItem removeItem;

	MyRemoveActionAdapter(PlotPopupMenu adaptee, JMenuItem removeGraphItem) {
		this.adaptee = adaptee;
		removeItem = removeGraphItem;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.plot.removeGraph(removeItem.getText());
		adaptee.graph_remove.remove(removeItem);
	}
}
