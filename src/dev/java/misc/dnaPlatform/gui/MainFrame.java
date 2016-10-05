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

import java.awt.*;
import javax.swing.*;

/**
 * <p>
 * Title: MainFrame
 * </p>
 * 
 * <p>
 * Description: This is the main Frame of the DnaGUI tool. This class contains
 * CompressionModels, SequenceData objects and Functions. It also contains
 * Panels to display information content, information about SequenceData and a
 * panel to display other graphs created by compression models.
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
public class MainFrame extends JFrame {
	public static final long serialVersionUID = 1234567890;

	MainPanel mainPanel;

	public MainFrame() {
		init();
		// Center the window
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension frameSize = this.getSize();
		if (frameSize.height > screenSize.height) {
			frameSize.height = screenSize.height;
		}
		if (frameSize.width > screenSize.width) {
			frameSize.width = screenSize.width;
		}

		this.setLocation((screenSize.width - frameSize.width) / 2,
				(screenSize.height - frameSize.height) / 2);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		boolean packFrame = false;// Dont know why it is here--Minh

		// Validate frames that have preset sizes
		// Pack frames that have useful preferred size info, e.g. from their
		// layout
		if (packFrame) {
			this.pack();
		} else {
			this.validate();
		}

		this.setVisible(true);
	}

	public void init() {
		mainPanel = new MainPanel();
		this.add(mainPanel);
		this.setJMenuBar(mainPanel.mainMenu);
		this.setTitle("DNA Platform");
		this.setSize(new Dimension(900, 500));

	}
}
