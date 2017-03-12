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

import java.util.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
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

public class PlotToolBar extends JToolBar {
	public static final long serialVersionUID = MainFrame.serialVersionUID;
	/* GUI elements for PlotMenuBar */
	JLabel plotLabel = new JLabel("Working with plot:    ");
	@SuppressWarnings("rawtypes")
	JComboBox plotBox = new JComboBox();

	JButton removeButton = new JButton("Remove Plot");

	// JButton viewButton = new JButton("View");

	JPanel zoomPanel = new JPanel();

	MainPanel myContainer;

	ImageIcon zoomInIcon;
	ImageIcon zoomOutIcon;
	JButton zoomOut;
	JButton zoomIn;

	JLabel coordLabel = new JLabel("(0,0)");

	// private InfoContentPanel selectedPanel =null;
	public PlotToolBar(MainPanel panel) {

		super();
		myContainer = panel;

		try {
			jbInit();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception {
		this.setFloatable(false);

		zoomInIcon = myContainer.loadImage("images/ZoomIn16.gif", "zoom in");
		zoomOutIcon = myContainer.loadImage("images/ZoomOut16.gif", "zoom out");

		zoomOut = new JButton(zoomOutIcon);
		zoomIn = new JButton(zoomInIcon);

		zoomIn.setFont(new java.awt.Font("Monospaced", Font.BOLD, 12));

		zoomIn.setMargin(new Insets(0, 4, 0, 4));
		zoomIn.addActionListener(new PlotToolBar_zoomIn_actionAdapter(this));

		zoomOut.setFont(new java.awt.Font("Monospaced", Font.BOLD, 12));
		zoomOut.setMargin(new Insets(0, 4, 0, 4));
		zoomOut.addActionListener(new PlotToolBar_zoomOut_actionAdapter(this));

		coordLabel.setBorder(BorderFactory.createLoweredBevelBorder());
		coordLabel.setMaximumSize(new Dimension(200, 20));
		coordLabel.setMinimumSize(new Dimension(100, 20));
		coordLabel.setPreferredSize(new Dimension(100, 24));
		coordLabel.setHorizontalAlignment(SwingConstants.CENTER);

		plotBox.setMaximumSize(new Dimension(100, 24));
		plotBox.setMinimumSize(new Dimension(100, 24));
		plotBox.setPreferredSize(new Dimension(100, 24));
		plotBox.addItemListener(new PlotToolBar_plotBox_itemAdapter(this));

		plotLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		zoomPanel.setMinimumSize(new Dimension(53, 24));
		zoomPanel.setPreferredSize(new Dimension(53, 24));

		zoomPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));
		zoomPanel.add(zoomIn);
		zoomPanel.add(zoomOut);

		this.setLayout(new GridLayout(1, 0));

		// removeButton.setMargin(new Insets(0, 4, 0, 4));
		removeButton.setMinimumSize(new Dimension(30, 24));
		removeButton.setPreferredSize(new Dimension(30, 24));
		removeButton
				.addActionListener(new PlotToolBar_removeButton_actionAdapter(
						this));
		// removeButton.setFont(new java.awt.Font("Monospaced", Font.BOLD, 12));

		add(plotLabel);
		add(plotBox);
		add(removeButton);
		add(zoomPanel);
		add(coordLabel);

		// Load images

	}

	/**
	 * method sets the selected panel to a given panel
	 * 
	 * @param InfoContentPanel
	 * 
	 */
	// public void setSelectedPanel(InfoContentPanel p){
	// this.selectedPanel=p;
	// System.out.println(">>>>");
	// }

	@SuppressWarnings("unchecked")
	public void updateInfoContentPanels(Vector<InfoContentPanel> panels) {
		plotBox.removeAllItems();

		Iterator<InfoContentPanel> it = panels.iterator();
		while (it.hasNext()) {
			plotBox.addItem((InfoContentPanel) it.next());
		}
	}

	// select InfocontentPanel
	// cause alot of exception
	public void selectInfoContentPanel(InfoContentPanel panel) {
		plotBox.setSelectedItem(panel);
	}

	public void plotBox_itemStateChanged(ItemEvent e) {
		myContainer.selectInformationContentPanel((InfoContentPanel) plotBox
				.getSelectedItem());
	}

	public void zoomIn_actionPerformed(ActionEvent e) {
		((InfoContentPanel) plotBox.getSelectedItem()).zoomIn_actionPerformed();
		// this.selectedPanel.zoomIn_actionPerformed();
	}

	public void removeButton_actionPerformed(ActionEvent e) {
		System.out.println(myContainer.plotList.size());
		myContainer.plotList.remove(plotBox.getSelectedIndex());
		System.out.println(myContainer.plotList.size());

		myContainer.repaint();
		// Remove the item from the drop down list
		plotBox.removeItem((plotBox.getSelectedItem()));
		// myContainer.plotList.remove(plotBox.getSelectedIndex());

		//
		// ((InfoContentPanel)).zoomIn_actionPerformed();
		// this.selectedPanel.zoomIn_actionPerformed();
	}

	public void zoomOut_actionPerformed(ActionEvent e) {
		((InfoContentPanel) plotBox.getSelectedItem())
				.zoomOut_actionPerformed();
		// this.selectedPanel.zoomOut_actionPerformed();
	}

}

class PlotToolBar_plotBox_itemAdapter implements ItemListener {
	private PlotToolBar adaptee;

	PlotToolBar_plotBox_itemAdapter(PlotToolBar adaptee) {
		this.adaptee = adaptee;
	}

	public void itemStateChanged(ItemEvent e) {
		adaptee.plotBox_itemStateChanged(e);
	}
}

class PlotToolBar_zoomIn_actionAdapter implements ActionListener {
	private PlotToolBar adaptee;

	PlotToolBar_zoomIn_actionAdapter(PlotToolBar adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.zoomIn_actionPerformed(e);
	}
}

class PlotToolBar_zoomOut_actionAdapter implements ActionListener {
	private PlotToolBar adaptee;

	PlotToolBar_zoomOut_actionAdapter(PlotToolBar adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.zoomOut_actionPerformed(e);
	}

}

class PlotToolBar_removeButton_actionAdapter implements ActionListener {
	private PlotToolBar adaptee;

	PlotToolBar_removeButton_actionAdapter(PlotToolBar adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.removeButton_actionPerformed(e);
	}
}
