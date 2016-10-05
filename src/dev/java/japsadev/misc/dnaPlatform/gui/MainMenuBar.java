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

package japsadev.misc.dnaPlatform.gui;

import javax.swing.*;

import japsadev.misc.dnaPlatform.compModel.*;
import japsadev.misc.dnaPlatform.function.*;
import japsadev.misc.dnaPlatform.sequence.*;

import java.awt.event.*;
import java.io.File;

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
public class MainMenuBar extends JMenuBar {
	public static final long serialVersionUID = MainFrame.serialVersionUID;
	// The current directory to read/write file
	MainPanel mainPanel;

	// File
	JMenu menuFile = new JMenu("File");
	JMenuItem menuFile_exit = new JMenuItem("Exit");
	JMenuItem menuFile_print = new JMenuItem("Print graph");
	JMenu menuFile_import = new JMenu("Import Sequence");
	// JMenuItemImport: inner class from MainPanel class
	JMenuItemImport[] menuFile_import_sequences;

	// Run
	JMenu menuRun = new JMenu("Run");
	JMenu menuRun_model = new JMenu("Run Model");
	JMenu menuRun_function = new JMenu("Run Function");
	JMenuItemFunction[] menuRun_runFunction;

	// View
	JMenu menuView = new JMenu("Plots");
	JMenu menuView_plots = new JMenu("Number of Plots");
	JMenuItem menuView_plots_add = new JMenuItem("Add plot");
	JMenu menuView_plots_remove = new JMenu("Remove plot");
	JMenu menuView_selectPlot = new JMenu("Select plot to work with");

	// Help
	JMenu menuHelp = new JMenu("Help");
	JMenuItem menuHelp_about = new JMenuItem("About");

	public MainMenuBar(MainPanel mainPanel) {
		this.mainPanel = mainPanel;
		initMenu();
	}

	/**
	 * This method creates all the JMenuItemImport's needed to import sequences
	 * into the tool. This method contains an array that holds all the types of
	 * sequences that can be read from files.
	 * 
	 * 
	 * @return JMenuItemImport[]
	 */

	private JMenuItemImport[] createJMenuItemsImport() {
		// This array contains all the types of sequences that can be
		// imported into the tool from files
		SequenceData[] types = { new DNASequenceData(),
				// new CharSequenceData(),
				new DoubleSequenceData(), new AnnotationSequenceData() };

		// create array to hold all JMenuItemImport items
		JMenuItemImport[] seqImport = new JMenuItemImport[types.length];

		// create JMenuItemImport elements in array depending on
		// the seqTypes available
		for (int i = 0; i < seqImport.length; i++)
			seqImport[i] = new JMenuItemImport(types[i], mainPanel);

		return seqImport;
	}

	// initialise the menu
	private void initMenu() {
		/**** Creating all the Functions and CompressionModels ****/
		Function[] functions = mainPanel.functions;

		/**** Setting up the Menu ****/
		// import menu
		menuFile_import.add(new JMenuItemFormatedImport(mainPanel));
		menuFile_import_sequences = createJMenuItemsImport();
		for (int i = 0; i < menuFile_import_sequences.length; i++)
			menuFile_import.add(menuFile_import_sequences[i]);

		// menu item for paste DNA sequence
		JMenuItem pasteSequence = new JMenuItem("Paste DNA");
		pasteSequence
				.addActionListener(new MainMenu_pasteSequence_actionAdapter(
						mainPanel));
		menuFile_import.add(pasteSequence);

		menuFile_print
				.addActionListener(new MainMenu_menuFile_print_actionAdapter(
						mainPanel));// in MainPanel class
		menuFile_exit
				.addActionListener(new MainMenu_menuFile_exit_actionAdapter(
						mainPanel));// in MainPanel class
		menuView_plots_add
				.addActionListener(new MainMenu_menuView_plots_add_actionAdapter(
						mainPanel));// in MainPanel class

		// adding JMenuItem to menuFile
		menuFile.add(menuFile_import);
		menuFile.add(menuFile_print);
		menuFile.addSeparator();
		menuFile.add(menuFile_exit);

		// creating menu items for functions and models
		menuRun_runFunction = new JMenuItemFunction[functions.length];

		for (int i = 0; i < functions.length; i++) {
			menuRun_runFunction[i] = new JMenuItemFunction(functions[i],
					mainPanel);

			if (functions[i] instanceof CompressionModel) {
				menuRun_model.add(menuRun_runFunction[i]);
			} else
				menuRun_function.add(menuRun_runFunction[i]);
		}

		menuRun.add(menuRun_model);
		menuRun.add(menuRun_function);

		menuView_plots.add(menuView_plots_add);
		menuView_plots.add(menuView_plots_remove);
		menuView.add(menuView_selectPlot);
		menuView.add(menuView_plots);

		menuHelp_about.addActionListener(new MainMenu_helpAbout_actionAdapter(
				mainPanel));
		menuHelp.add(menuHelp_about);

		add(menuFile);
		add(menuRun);
		add(menuView);
		add(menuHelp);
	}

	/**
	 * Exit application
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void menuFile_exit_actionPerformed(ActionEvent e) {
		int value = JOptionPane.showConfirmDialog(this, "Do you want to exit?",
				"DNA Tool", JOptionPane.OK_CANCEL_OPTION);

		if (value == JOptionPane.OK_OPTION)
			System.exit(0);
	}

	/**
	 * Print InfoContentGraph
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void menuFile_print_actionPerformed(ActionEvent e) {
		PrintUtilities.printComponent(mainPanel.infoContent
				.getInfoContentGraph());
	}
}

/**
 * <p>
 * Title: JMenuItemImport
 * </p>
 * 
 * <p>
 * Description: This class extens JMenuItem and has an action listener that
 * reads sequence data from a file and creates a SequenceData object using the
 * ReadFileFunction, which runs on the RunFunction thread.
 * 
 * </p>
 */

class JMenuItemImport extends JMenuItem {
	public static final long serialVersionUID = MainFrame.serialVersionUID;
	private SequenceData sequenceType;
	MainPanel mainPanel;

	JMenuItemImport(SequenceData sequenceClass, MainPanel mainPanel) {
		super(sequenceClass.toString());
		sequenceType = sequenceClass;
		this.mainPanel = mainPanel;

		this.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				importSequence(e);
			}
		});
	}

	private void importSequence(ActionEvent e) {
		JFileChooser fileChooser = new JFileChooser(mainPanel.currentDir);

		String sequenceFile = "";
		fileChooser.setDialogTitle("Import " + sequenceType.getProperty());

		if (JFileChooser.APPROVE_OPTION == fileChooser.showOpenDialog(this)) {
			File file = fileChooser.getSelectedFile().getAbsoluteFile();
			mainPanel.currentDir = file.getAbsoluteFile();

			// set file of current model:
			sequenceFile = file.getAbsolutePath();
			mainPanel.statusBar.setText("sequence file: " + file.toString());
		}
		if (sequenceFile != "" && sequenceFile != null) {

			// create SequenceData of the same type as sequenceType
			SequenceData data = sequenceType.getNewSequenceData();
			mainPanel.currentData = data;

			// reading file in separate thread with readFun

			ReadFileFunction readFun = new ReadFileFunction();
			readFun.setFile(sequenceFile);
			RunFunction runFun = new RunFunction(readFun, data, mainPanel);
			runFun.start();

			System.out.println("Created a sequence data of type"
					+ data.getClass());

		} else
			mainPanel.statusBar.setText("couldn't open file: " + sequenceFile);
	}
}

/**
 * <p>
 * Title: JMenuItemImport
 * </p>
 * 
 * <p>
 * Description: This class extens JMenuItem and has an action listener that
 * reads sequence data from a file and creates a SequenceData object using the
 * ReadFileFunction, which runs on the RunFunction thread.
 * 
 * </p>
 */

class JMenuItemFormatedImport extends JMenuItem {
	public static final long serialVersionUID = MainFrame.serialVersionUID;

	MainPanel mainPanel;

	JMenuItemFormatedImport(MainPanel mainPanel) {
		super("Import Rich Format");
		this.mainPanel = mainPanel;

		this.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				importSequence(e);
			}
		});
	}

	private void importSequence(ActionEvent e) {
		mainPanel.importSequence(e);
	}
}

/**
 * <p>
 * Title: MainFrame_paste_actionAdapter
 * <p>
 * this class implements ActionListener when use chooses paste sequence, it will
 * show a textArea user can paste text into it
 * 
 */
class MainMenu_pasteSequence_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	// constructor
	public MainMenu_pasteSequence_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	// acitonPerformed
	public void actionPerformed(ActionEvent e) {
		adaptee.pasteSequence_actionPerformed(e);
	}
}

/**
 * <p>
 * Title: MainMenu_readFormated_actionAdapter
 * <p>
 * this class implements ActionListener
 * 
 */

class MainMenu_readFormated_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	// constructor
	public MainMenu_readFormated_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	// actionPerformed
	public void actionPerformed(ActionEvent e) {
		adaptee.pasteSequence_actionPerformed(e);
	}
}

class MainMenu_menuView_plots_add_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_menuView_plots_add_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.createPlot();
	}
}

class MainMenu_menuFile_exit_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_menuFile_exit_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.menuFile_exit_actionPerformed(e);
	}
}

class MainMenu_menuFile_print_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_menuFile_print_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.menuFile_print_actionPerformed(e);
	}
}

class MenuFile_exit_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MenuFile_exit_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.menuFile_exit_actionPerformed(e);
	}
}

class MenuFile_print_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MenuFile_print_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.menuFile_print_actionPerformed(e);
	}
}

/*
 * This inner class implements ActionListerner When the event of this type
 * occurs (user has chosen to set the axis axis) the class performs
 * actionPerformed method to call the method popup_setXaxis_actionPerformed
 */
class MainMenu_popup_setXaxis_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_setXaxis_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.popup_setXaxis_actionPerformed(e);
	}
}

/*
 * class popup_property shows the property of the sequence when the user chooses
 * to show the sequence information author: Hoang Nguyen
 */
class MainMenu_popup_property_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_property_actionAdapter(MainPanel adaptee) {// ,SequenceData
																// sequenceData)
																// {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.popup_property_actionPerformed(e);
	}
}

class MainMenu_popup_redraw_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_redraw_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.popup_redraw_actionPerformed(e);
	}
}

class MainMenu_popup_remove_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_remove_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.popup_remove_actionPerformed(e);
	}
}

class MainMenu_popup_save_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_save_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		adaptee.popup_save_actionPerformed(e);
	}

}

class MainMenu_popup_saveAnnotation_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	MainMenu_popup_saveAnnotation_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	public void actionPerformed(ActionEvent e) {
		// TODO: Ask for annotation

		adaptee.popup_saveAnno_actionPerformed(e);
	}

}

/**
 * <p>
 * Title: MainFrame_helpAbout_actionAdapter
 * <p>
 * this class implements ActionListener when use chooses help menu, it will show
 * a textField contains information about the product
 * 
 */
class MainMenu_helpAbout_actionAdapter implements ActionListener {
	private MainPanel adaptee;

	// constructor
	public MainMenu_helpAbout_actionAdapter(MainPanel adaptee) {
		this.adaptee = adaptee;
	}

	// actionPerformed: override the parent class
	public void actionPerformed(ActionEvent e) {
		adaptee.helpMenu_about_actionPerformed(e);
	}
}

class JMenuItemSelectPlot extends JMenuItem {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	InfoContentPanel me;
	MainPanel mainPanel;

	JMenuItemSelectPlot(InfoContentPanel plot, MainPanel mainPanel) {
		super(plot.toString());
		me = plot;
		this.mainPanel = mainPanel;

		this.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				selectInformationContentPanel();

			}
		});
	}

	private void selectInformationContentPanel() {
		mainPanel.plotToolBar.selectInfoContentPanel(me);
		// mainPanel.selectInformationContentPanel(me);Select info
	}

}

class JMenuItemRemovePlot extends JMenuItem {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private InfoContentPanel me;
	MainPanel mainPanel;

	JMenuItemRemovePlot(InfoContentPanel plot, MainPanel mainPanel) {
		super(plot.toString());
		me = plot;
		this.mainPanel = mainPanel;

		this.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				removePlot();
			}

		});
	}

	private void removePlot() {

		for (int i = 0; i < mainPanel.mainMenu.menuView_selectPlot
				.getItemCount(); i++) {
			JMenuItemSelectPlot item = (JMenuItemSelectPlot) mainPanel.mainMenu.menuView_selectPlot
					.getItem(i);
			if (item.me == me) {
				mainPanel.mainMenu.menuView_selectPlot.remove(item);
				System.out.println("Removing " + item.me);
			}
		}

		mainPanel.mainMenu.menuView_plots_remove.remove(this);

		mainPanel.plotList.remove(me);

		if (mainPanel.infoContent == me) {
			if (mainPanel.plotList.size() > 0) {
				mainPanel.infoContent = mainPanel.plotList.firstElement();
				mainPanel.infoContent.highlight(true);
			} else
				mainPanel.infoContent = null;
		}

		mainPanel.plotSplitPanel.remove(mainPanel.plotsPanel);
		mainPanel.plotsPanel = mainPanel
				.addPlots(mainPanel.plotList.size() - 1);
		mainPanel.plotSplitPanel.add(mainPanel.plotsPanel);

		// mainPanel.plotToolBar.updateInfoContentPanels(mainPanel.plots);---------------------------------------------------

	}
}

class JMenuItemFunction extends JMenuItem {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Function fun;
	private MainPanel mainPanel;

	JMenuItemFunction(Function function, MainPanel mainPanel) {
		super(function.toString());
		fun = function;
		this.mainPanel = mainPanel;

		this.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				runFunction(e);
			}
		});
	}

	private void runFunction(ActionEvent e) {
		RunFunction runFunction = new RunFunction(fun, mainPanel.currentData,
				mainPanel);
		runFunction.start();
	}

	public Function getFunction() {
		return fun;
	}

}
