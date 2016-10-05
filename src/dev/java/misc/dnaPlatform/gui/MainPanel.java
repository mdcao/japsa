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

import misc.dnaPlatform.OptionsHandle;
import misc.dnaPlatform.compModel.*;
import misc.dnaPlatform.function.*;
import misc.dnaPlatform.sequence.*;

import java.awt.event.*;
import java.io.*;
import java.net.URL;
import java.util.*;

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
public class MainPanel extends JPanel {
	public static final long serialVersionUID = 1234567890;

	File currentDir = new File(".");
	MainMenuBar mainMenu;

	/*-- array to hold all functions that can be executed on SequenceData
	including compression models --*/
	Function[] functions;

	/*-- vector to hold SequenceData objects --*/
	Vector<SequenceData> sequenceData = new Vector<SequenceData>();
	SequenceData currentData = null; // used to select data to be modified

	/*--  Panels  --*/
	JPanel statusPane = new JPanel();// The panel to hold the status bar
	JPanel rightPanel = new JPanel();

	JSplitPane plotSplitPanel = new JSplitPane();
	JSplitPane plotsPanel = new JSplitPane();// Panel to contain plot
	JSplitPane centrePanel = new JSplitPane();
	JScrollPane extrasPane = new JScrollPane();

	DataTreePanel dataTree = new DataTreePanel();
	InfoContentPanel infoContent = null; // new InfoContentPanel("Plot 1");

	Vector<InfoContentPanel> plotList = new Vector<InfoContentPanel>();

	/* plot toolbar */
	PlotToolBar plotToolBar = new PlotToolBar(this);

	/*-- popup menu for SequenceData in dataTree --*/
	JPopupMenu popup = new JPopupMenu();
	JMenu popup_model = new JMenu("Run Model");
	JMenu popup_function = new JMenu("Run Function");
	JMenuItemFunction[] popup_runFunction;
	JMenuItem popup_remove = new JMenuItem("Remove");
	JMenuItem popup_redraw = new JMenuItem("Redraw");
	JMenuItem popup_save = new JMenuItem("Save");
	JMenuItem popup_save_ann = new JMenuItem("Save with annotation");
	JMenuItem popup_setXaxis = new JMenuItem("Set as x-axis");
	JMenuItem popup_property = new JMenuItem("Properties");

	JPanel buttonsPanel = new JPanel();

	/*--  progress label and bar  --*/
	JLabel statusBar = new JLabel("");
	JProgressBar progressBar = new JProgressBar();

	public MainPanel() {
		initPanel();
	}

	/**
	 * This method is where all compression models and functions to be used in
	 * the GUI are created. New functions and compression models that implement
	 * the Function or CompressionModel interface can be added to the GUI by
	 * creating them and adding them to the vector of functions this method
	 * returns.
	 * 
	 * @return Function[] containing all Functions and CompressionModels to be
	 *         used in the GUI
	 */

	private Function[] createFunctionsModels() {

		// adding compression models
		Function[] funcs = {
				new FuzzyModel(),
				new MarkovModel(),
				new ExpertCompressionModel(),
				// adding functions
				new AppendFunction(), new ReverseFunction(),
				new SelectFunction(), new DNAComplementFunction(),
				new SmoothingFunction(), new DifferenceFunction(),
				new NegateFunction(), new ThresholdFunction(),
				new FilterFeatureFunction() };// maybe don't add converter
												// here/add in xaxis

		return funcs;
	}

	/**
	 * Set up the centre panel, most things appear here
	 */
	private void initCentrePane() {

		extrasPane.setMaximumSize(new Dimension(32767, 2));// Extrapanel to
															// cover upo unused
															// space?
		extrasPane.setMinimumSize(new Dimension(19, 2));

		/**** Adding panels ****/
		plotsPanel = addPlots(0);// The current plot panel

		plotSplitPanel.add(plotsPanel, JSplitPane.LEFT);// = top
		plotSplitPanel.add(extrasPane, JSplitPane.RIGHT);// = down

		plotSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);
		plotSplitPanel.setDividerLocation(800);

		rightPanel.setLayout(new BorderLayout());
		rightPanel.add(plotToolBar, BorderLayout.NORTH);
		rightPanel.add(plotSplitPanel, BorderLayout.CENTER);

		centrePanel.add(dataTree, JSplitPane.LEFT);
		centrePanel.add(rightPanel, JSplitPane.RIGHT);
		centrePanel.setDividerLocation(180);
	}

	private void initStatusPane() {
		progressBar.setPreferredSize(new Dimension(200, 14));
		statusPane.add(statusBar);
		statusPane.add(progressBar);
	}

	/**
	 * Setting up all the popups
	 */
	private void initPopups() {// throws Exception {
		/**** Setting up popup menu ****/
		popup_remove.addActionListener(new MainMenu_popup_remove_actionAdapter(
				this));
		popup_redraw.addActionListener(new MainMenu_popup_redraw_actionAdapter(
				this));
		popup_save
				.addActionListener(new MainMenu_popup_save_actionAdapter(this));
		popup_save_ann
				.addActionListener(new MainMenu_popup_saveAnnotation_actionAdapter(
						this));

		popup_setXaxis
				.addActionListener(new MainMenu_popup_setXaxis_actionAdapter(
						this));
		popup_property
				.addActionListener(new MainMenu_popup_property_actionAdapter(
						this));

		// creating popup menu items for functions and models
		popup_runFunction = new JMenuItemFunction[functions.length];
		for (int i = 0; i < functions.length; i++) {
			popup_runFunction[i] = new JMenuItemFunction(functions[i], this);
			if (functions[i] instanceof CompressionModel) {
				popup_model.add(popup_runFunction[i]);
			} else
				popup_function.add(popup_runFunction[i]);
		}

		popup.add(popup_model);
		popup.add(popup_function);
		popup.add(popup_remove);
		popup.add(popup_redraw);
		popup.add(popup_save);
		popup.add(popup_save_ann);
		popup.add(popup_setXaxis);
		popup.add(popup_property);

	}

	private void initButtonPanel() {// throws Exception {
		/**************** The panel to add buttons ***************************/
		JPanel smallerPanel = new JPanel();

		smallerPanel.setLayout(new GridLayout(1, 3));
		buttonsPanel.setLayout(new BorderLayout());

		ImageIcon icon = loadImage("images/Open16.gif", "open new DNA file");//

		JButton openFileButton = new JButton(icon);
		// variable acts like an holder, inner class will access it
		final MainPanel mainPanel = this;
		// add action Listener to the openFileButton
		openFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				importSequence(e);
			}
		}

		);

		icon = loadImage("images/Import16.gif", "open window to paste DNA file");

		JButton pasteADNSequenceButton = new JButton(icon);
		pasteADNSequenceButton
				.addActionListener(new MainMenu_pasteSequence_actionAdapter(
						mainPanel));

		// this button is used to create a new plot
		JButton addNewPlotButton = new JButton("+");
		addNewPlotButton
				.addActionListener(new MainMenu_menuView_plots_add_actionAdapter(
						mainPanel));

		smallerPanel.add(openFileButton);
		smallerPanel.add(pasteADNSequenceButton);
		smallerPanel.add(addNewPlotButton);

		buttonsPanel.add(smallerPanel, BorderLayout.WEST);
		/********************* End of button panel ******************************/

	}

	private void initPanel() {// throws Exception {
		this.setSize(new Dimension(900, 500));
		this.setName(this.getName() + ".contentPane");
		this.setLayout(new BorderLayout() {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			/**
			 * This BorderLayout subclass maps a null constraint to CENTER.
			 * Although the reference BorderLayout also does this, some VMs
			 * throw an IllegalArgumentException.
			 */
			public void addLayoutComponent(Component comp, Object constraints) {
				if (constraints == null) {
					constraints = BorderLayout.CENTER;
				}
				super.addLayoutComponent(comp, constraints);
			}
		});

		/**** Creating all the Functions and CompressionModels ****/
		functions = createFunctionsModels();

		/* setting plot toolbar */

		initButtonPanel();
		initCentrePane();
		initStatusPane();
		initPopups();

		/**** add elements to contentPane ****/
		this.add(buttonsPanel, BorderLayout.NORTH);
		this.add(centrePanel, BorderLayout.CENTER);
		this.add(statusPane, java.awt.BorderLayout.SOUTH);

		/**** adding mouse listener to tree in dataTree panel ****/
		JTree tree = dataTree.getTree();
		tree.addMouseListener(new tree_popupListener());

		mainMenu = new MainMenuBar(this);
	}

	protected ImageIcon loadImage(String imgPath, String desc) {
		URL url = this.getClass().getResource(imgPath);

		if (url == null)
			return new ImageIcon(imgPath, desc);

		return new ImageIcon(url, desc);

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
	 * open a new frame : about frame do nothing just about the information of
	 * the program: from Monash
	 * 
	 */
	public void helpMenu_about_actionPerformed(ActionEvent e) {
		JFrame aboutFrame = new JFrame("about");
		JTextField about = new JTextField();
		about.setText("From Monash University ");
		about.setEditable(false);
		aboutFrame.add(about);
		aboutFrame.setSize(200, 100);
		aboutFrame.setLocation(300, 400);
		aboutFrame.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		aboutFrame.setVisible(true);

	}

	/**
	 * create a new frame window contains TextArea so that the user can cut and
	 * paste the sequence to the program instead of selecting the file
	 * 
	 * @param ACtionEvent
	 * 
	 */
	public void pasteSequence_actionPerformed(ActionEvent e) {
		final JFrame pasteSequenceFrame = new JFrame("Paste DNA Sequence");
		JPanel textPanel = new JPanel();
		// just a variable holder for the anonymous class
		final MainPanel mainPanel = this;
		JPanel buttonPanel = new JPanel();
		pasteSequenceFrame.setLayout(new BorderLayout());
		final JTextArea paste = new JTextArea(40, 52);
		paste.setEditable(true);
		JScrollPane scrollText = new JScrollPane(paste,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		textPanel.add(scrollText);
		// okie button
		JButton okieButton = new JButton("OK");
		JButton cancelButton = new JButton("Cancel");
		okieButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String text = paste.getText();
				DNASequenceData data = new DNASequenceData();
				data.readDataFromString(text);
				currentData = data;
				(new RunFunction(null, null, mainPanel)).attachSequence(data);

				System.out.println("Created a sequence data of type"
						+ data.getClass());
				pasteSequenceFrame.dispose();
			}
		}

		);

		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				pasteSequenceFrame.dispose();
			}
		});
		buttonPanel.add(okieButton);
		buttonPanel.add(cancelButton);
		pasteSequenceFrame.add(textPanel, BorderLayout.NORTH);
		pasteSequenceFrame.add(buttonPanel, BorderLayout.SOUTH);
		pasteSequenceFrame.setSize(600, 900);
		pasteSequenceFrame.setLocation(0, 0);
		pasteSequenceFrame
				.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		pasteSequenceFrame.setVisible(true);

	}

	/**
	 * Print InfoContentGraph
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void menuFile_print_actionPerformed(ActionEvent e) {
		PrintUtilities.printComponent(infoContent.getInfoContentGraph());
	}

	/**
	 * Delete selected data from dataTree and information content.
	 * 
	 * @param e
	 *            ActionEvent
	 */

	public void popup_remove_actionPerformed(ActionEvent e) {
		// TODO: Convert this to make it generic: remove, draw
		if (currentData == null)
			return;

		if (currentData instanceof DoubleSequenceData) {
			// go through all the plotList to remove
			Iterator<InfoContentPanel> it = plotList.iterator();
			while (it.hasNext()) {
				InfoContentPanel p = it.next();

				p.removeGraph(currentData.toString());
			}

			// if (infoContent != null)
			// infoContent.removeGraph(currentData.toString());
		}

		else if (currentData instanceof AnnotationSequenceData) {
			// go through all the plotList to remove
			Iterator<InfoContentPanel> it = plotList.iterator();
			while (it.hasNext()) {
				InfoContentPanel p = it.next();

				p.removeFeatures((AnnotationSequenceData) currentData);
			}

			// infoContent.removeGraph(currentData.toString()); ??
		}

		sequenceData.removeElement(currentData);
		currentData = null;

		/* remove history of sequenceData object from dataTree */
		dataTree.removeSelectedPath();

	}

	/**
	 * Save selected data from dataTree
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void popup_save_actionPerformed(ActionEvent e) {

		if (currentData == null)
			return;

		JFileChooser fileChooser = new JFileChooser(currentDir);
		fileChooser.setDialogTitle("Save File");
		// choose which type of file to save

		if (JFileChooser.APPROVE_OPTION == fileChooser.showOpenDialog(this)) {
			File file = fileChooser.getSelectedFile().getAbsoluteFile();

			currentDir = file.getAbsoluteFile();

			SaveFileFunction saveFunction = new SaveFileFunction(file);
			RunFunction runFun = new RunFunction(saveFunction, currentData,
					this);

			System.out.println("ready to run");
			runFun.start();
		}
	}

	/**
	 * Save selected data from dataTree
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void popup_saveAnno_actionPerformed(ActionEvent e) {

		if (currentData == null)
			return;

		JFileChooser fileChooser = new JFileChooser(currentDir);
		fileChooser.setDialogTitle("Save File");
		// choose which type of file to save

		if (JFileChooser.APPROVE_OPTION == fileChooser.showOpenDialog(this)) {
			File file = fileChooser.getSelectedFile().getAbsoluteFile();
			currentDir = file.getAbsoluteFile();

			SaveFormatFileFunction saveFunction = new SaveFormatFileFunction(
					file);
			RunFunction runFun = new RunFunction(saveFunction, currentData,
					this);

			System.out.println("ready to run");
			runFun.start();
		}
	}

	/**
	 * Function to redraw data in plot
	 * 
	 * @param e
	 *            ActionEvent
	 */
	public void popup_redraw_actionPerformed(ActionEvent e) {
		if (currentData == null || infoContent == null)
			return;
		statusBar.setText("");
		if (currentData instanceof DoubleSequenceData) {
			if (!infoContent.paintInfoContent((DoubleSequenceData) currentData))
				statusBar.setText("Graph is already in plot");
		} else if (currentData instanceof AnnotationSequenceData) {
			if (!infoContent.addFeatures((AnnotationSequenceData) currentData))
				statusBar.setText("Cant draw feature list");
		}
	}

	// set the Xaxis of the panel
	// draw the protein
	// draw the features(anottation)
	public void popup_setXaxis_actionPerformed(ActionEvent e) {
		if (currentData == null || infoContent == null)
			return;

		if (currentData instanceof CharSequenceData) {
			System.out.println("Changing sequence of plot");
			infoContent.setSequence(((CharSequenceData) currentData)
					.getCharData());
		}

		if (currentData instanceof AnnotationSequenceData) {
			if (!infoContent.addFeatures((AnnotationSequenceData) currentData))
				statusBar.setText("Cant draw feature list");
		}
	}

	/*
	 * shows the property of the current sequence
	 * 
	 * @param ActionEvent Hoang Nguyen
	 */
	public void popup_property_actionPerformed(ActionEvent e) {

		if (currentData == null)
			return;

		String message = currentData.getProperty();
		JOptionPane.showMessageDialog(this, message, "sequence info",
				JOptionPane.PLAIN_MESSAGE);

		/*
		 * if(dnaSeq.sequenceType==DNASequenceData.FASTA) {
		 * 
		 * message +=
		 * "DESCRIPTION LINE:"+annotation.getProperty("DESCRIPTION_LINE") +"\n";
		 * message += "DESCRIPTION:"+annotation.getProperty("DESCRIPTION")
		 * +"\n"; JOptionPane.showMessageDialog(this, message, "sequence info",
		 * JOptionPane.PLAIN_MESSAGE); } else//read from genebank
		 * if(dnaSeq.sequenceType==DNASequenceData.GENBANK) {
		 * 
		 * message += "LOCUS      :"+annotation.getProperty("LOCUS")+"\n";
		 * message += "SIZE       :"+annotation.getProperty("SIZE")+"\n";
		 * message += "TYPE       :"+annotation.getProperty("TYPE")+"\n";
		 * message += "CIRCULAR   :"+annotation.getProperty("CIRCULAR")+"\n";
		 * message += "DIVISION   :"+annotation.getProperty("DIVISION")+"\n";
		 * message += "MDAT       :"+annotation.getProperty("MDAT")+"\n";
		 * message += "SOURCE     :"+annotation.getProperty("SOURCE")+"\n";
		 * message += "DEFINITION :"+annotation.getProperty("DEFINITION") +"\n";
		 * JOptionPane.showMessageDialog(this, message, "sequence info",
		 * JOptionPane.PLAIN_MESSAGE); }
		 * 
		 * /*********************************************************
		 */
	}

	/**
	 * Displays a dialog to update options to run a function model.
	 * 
	 * @param handle
	 *            OptionsHandle
	 * @return boolean Whether the options for function were set in dialog
	 */
	protected boolean setParameters(OptionsHandle handle) {
		OptionsHandleDialog dlg = new OptionsHandleDialog((Frame) null, handle,
				sequenceData);
		Dimension dlgSize = dlg.getPreferredSize();
		Dimension frmSize = getSize();
		java.awt.Point loc = getLocation();
		dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x,
				(frmSize.height - dlgSize.height) / 2 + loc.y);
		dlg.setModal(true);
		dlg.setVisible(true);

		return dlg.isOptionsHandleSet();
	}

	/**
	 * This is a recursive function to add all plots from index to JScrollPanes
	 * and returns the Panel where all plots were added. When the index is equal
	 * to zero this function returns the first InfoContentPanel
	 * 
	 * @param index
	 *            Index to add plots from
	 */
	JSplitPane addPlots(int index) {

		// add plots in split panels and return final JPanel
		JSplitPane pan = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
		if (index < 0 || index >= plotList.size()) {
			pan.add(new JPanel(), JSplitPane.BOTTOM);
		}

		else if (index == 0) {
			pan.add(plotList.elementAt(index), JSplitPane.BOTTOM);
		}

		else {
			pan.add(plotList.elementAt(index), JSplitPane.BOTTOM);
			pan.add(addPlots(index - 1), JSplitPane.TOP);
		}

		return pan;

	}

	public static final int MAX_NUMBER_PLOT = 4;

	/**
	 * This method creates a new InformationContentPanel to be displayed by this
	 * tool
	 */
	public void createPlot() {
		if (plotList.size() < MAX_NUMBER_PLOT) {

			// find right name for the plot
			int number = 1;
			Iterator<InfoContentPanel> it = plotList.iterator();
			while (it.hasNext()) {
				InfoContentPanel p = it.next();
				String name = p.toString();

				if (name.indexOf("Plot ") != -1) {
					String n = name.substring(name.lastIndexOf(' ') + 1);

					if (number == Integer.parseInt(n))
						number++;
				}
			}

			InfoContentPanel panel = new InfoContentPanel("Plot " + number,
					this);

			plotList.add(panel);

			plotSplitPanel.remove(plotsPanel);
			plotsPanel = addPlots(plotList.size() - 1);
			plotSplitPanel.add(plotsPanel);

			// select plot
			plotToolBar.selectInfoContentPanel(panel);

			// this.selectInformationContentPanel(panel);

			/* update menus */
			mainMenu.menuView_plots_remove.add(new JMenuItemRemovePlot(panel,
					this));
			mainMenu.menuView_selectPlot.add(new JMenuItemSelectPlot(panel,
					this));
			plotToolBar.updateInfoContentPanels(plotList);
		} else
			statusBar.setText("A maximum of four plots can be displayed");
	}

	public void updateInfoContentPlotMenu() {
		mainMenu.menuView_selectPlot.removeAll();
		mainMenu.menuView_plots_remove.removeAll();

		Iterator<InfoContentPanel> it = plotList.iterator();
		while (it.hasNext()) {
			InfoContentPanel panel = it.next();

			mainMenu.menuView_plots_remove.add(new JMenuItemRemovePlot(panel,
					this));
			mainMenu.menuView_selectPlot.add(new JMenuItemSelectPlot(panel,
					this));
		}

		// plotToolBar.updateInfoContentPanels(plots);-----------------------------------------------------------
	}

	/**
	 * This method selects the informationContentPanel to work with
	 * 
	 * @param panel
	 *            InfoContentPanel
	 */
	public void selectInformationContentPanel(InfoContentPanel panel) {
		if (panel == null)
			return;

		// set all other graphs to no highlight
		Iterator<InfoContentPanel> it = plotList.iterator();
		while (it.hasNext()) {
			InfoContentPanel p = it.next();
			p.highlight(p == panel);
		}

		infoContent = panel;
		// select InfoContentPanel in plotToolBar
		// plotToolBar.selectInfoContentPanel(panel);

		// after doing this, reset the menuView_selectPlot
		// this.plotToolBar.setSelectedPanel(panel);
		// plotToolBar.selectInfoContentPanel(panel);

	}

	protected void importSequence(ActionEvent e) {
		JFileChooser fileChooser = new JFileChooser(currentDir);

		String sequenceFile = "";
		fileChooser.setDialogTitle("Import Rich Format File");

		if (JFileChooser.APPROVE_OPTION == fileChooser.showOpenDialog(this)) {
			File file = fileChooser.getSelectedFile().getAbsoluteFile();
			currentDir = file.getAbsoluteFile();

			// set file of current model:
			sequenceFile = file.getAbsolutePath();
			statusBar.setText("Sequence file: " + file.toString());
		}
		if (sequenceFile != "" && sequenceFile != null) {

			// reading file in separate thread with readFun

			ReadFormatFileFunction readFun = new ReadFormatFileFunction();
			readFun.setFile(sequenceFile);
			RunFunction runFun = new RunFunction(readFun, null, this);
			runFun.start();

			// System.out.println("Created a sequence data of type" +
			// data.getClass());

		} else
			this.statusBar.setText("couldn't open file: " + sequenceFile);
	}

	/**
	 * <p>
	 * Title: JMenuItemFunction
	 * </p>
	 * 
	 * <p>
	 * Description: This class extens JMenuItem and has an action listener that
	 * calls a RunFunction thread to run functions.
	 * </p>
	 */

	/**
	 * 
	 * <p>
	 * Title: DataTreePanel_tree_mouseAdapter
	 * </p>
	 * 
	 * <p>
	 * Description: This is a MouseAdapter for the JTree contained in
	 * DataTreePanel. It has been added in this MainFrame instead of
	 * DataTreePanel to handle mouse clicks in this class.
	 * </p>
	 * 
	 * @author Julie Bernal
	 * @version 1.0
	 */
	@SuppressWarnings("rawtypes")
	class tree_popupListener extends MouseAdapter {

		@SuppressWarnings("unchecked")
		public void mousePressed(MouseEvent e) {
			maybeShowPopup(e);

			// enable and disable items in main menu according to current data
			for (int i = 0; i < mainMenu.menuRun_runFunction.length; i++) {
				Function fun = mainMenu.menuRun_runFunction[i].getFunction();
				Class[] classes = fun.getTypeSequenceData();

				mainMenu.menuRun_runFunction[i].setEnabled(false);
				for (int j = 0; j < classes.length; j++) {
					if (currentData != null
							&& classes[j].isAssignableFrom(currentData
									.getClass())) {
						mainMenu.menuRun_runFunction[i].setEnabled(true);
					}
				}
			}

		}

		@SuppressWarnings("unchecked")
		public void mouseReleased(MouseEvent e) {
			maybeShowPopup(e);

			// enable and disable items in main menu according to current data
			for (int i = 0; i < mainMenu.menuRun_runFunction.length; i++) {
				Function fun = mainMenu.menuRun_runFunction[i].getFunction();
				Class[] classes = fun.getTypeSequenceData();

				mainMenu.menuRun_runFunction[i].setEnabled(false);
				for (int j = 0; j < classes.length; j++) {
					if (currentData != null
							&& classes[j].isAssignableFrom(currentData
									.getClass())) {
						mainMenu.menuRun_runFunction[i].setEnabled(true);
					}
				}
			}
		}

		/**
		 * Shows popup for SequenceData object selected elements in popup menu
		 * are enabled according to the type of SequenceData selected.
		 * 
		 * @param e
		 *            MouseEvent
		 */
		@SuppressWarnings("unchecked")
		private void maybeShowPopup(MouseEvent e) {
			if (dataTree == null)
				return;
			String data = dataTree.getFirstSelectedPath();
			if (data == null)
				return;
			else {
				// get SequenceData with the same name as selected SequenceData
				Iterator<SequenceData> it = sequenceData.iterator();
				while (it.hasNext()) {
					SequenceData d = it.next();
					if (data.equals(d.toString())) {
						// select data from sequenceData using name in tree
						currentData = d;
						// enable functions and models that can be applied to
						// current data
						for (int i = 0; i < popup_runFunction.length; i++) {
							Function fun = popup_runFunction[i].getFunction();
							Class[] classes = fun.getTypeSequenceData();

							popup_runFunction[i].setEnabled(false);
							for (int j = 0; j < classes.length; j++) {
								if (classes[j].isAssignableFrom(currentData
										.getClass())) {
									popup_runFunction[i].setEnabled(true);
								}
							}
						}

						popup_property.setEnabled(true);
						popup_save_ann.setEnabled(d instanceof DNASequenceData);

						if (d instanceof CharSequenceData) {
							popup_redraw.setEnabled(false);
							popup_setXaxis.setEnabled(true);
						} else if (d instanceof DoubleSequenceData) {
							popup_redraw.setEnabled(true);
							popup_setXaxis.setEnabled(false);
						} else if (d instanceof AnnotationSequenceData) {
							popup_redraw.setEnabled(true);
							popup_setXaxis.setEnabled(false);
						}

						/* show popup menu if data in sequenceData */
						if (e.isPopupTrigger()) {
							popup.show(e.getComponent(), e.getX(), e.getY());
						}
						break;

					}
				} // while

			}

		}
	}
}

/**
 * <p>
 * Title: RunFunction
 * </p>
 * 
 * <p>
 * Description: This class is a Thread that runs functions.
 * </p>
 */
class RunFunction extends Thread {
	Function myFunction;
	SequenceData myInputData;
	MainPanel mainPanel;

	public RunFunction(Function function, SequenceData inputData,
			MainPanel mainPanel) {
		myFunction = function;
		myInputData = inputData;
		this.mainPanel = mainPanel;
	}

	public void attachMultiSequence(Iterator<SequenceData> itrData) {
		int i = 1;
		while (itrData.hasNext()) {
			System.out.println(i++);
			SequenceData data = itrData.next();
			attachSequence(data);
		}

	}

	public void attachSequence(SequenceData data) {
		Iterator<SequenceData> it = mainPanel.sequenceData.iterator();
		int i = 1;
		while (it.hasNext()) {
			SequenceData d = it.next();
			if (d.getClass().equals(data.getClass())) {
				i++;
				String name = d.getSequenceName();
				String subName[] = name.split("\\D");
				if (subName.length > 0) {
					int n = Integer.parseInt(subName[0]);
					if (i == n)
						i = n + 1;
				}

			}
		}
		data.setSequenceName(i + "");
		// Attached into sequence
		mainPanel.sequenceData.add(data);
		mainPanel.dataTree.addDataHistory(data.getHistory());

		// Draw if need to
		if (data instanceof DoubleSequenceData && mainPanel.infoContent != null) {
			mainPanel.statusBar.setText("Creating graph");
			mainPanel.infoContent.paintInfoContent((DoubleSequenceData) data);
		}
	}

	/**
	 * Runs a function. The steps required to run a function are: 1. Set
	 * parameters for function 2. Call method to execute function 3. If
	 * inputData is DoubleSequenceData repaint infoContent 4. Display
	 * information about outpout data in dataTree panel
	 * 
	 * The SequenceData object created by CompressionModel is stored in
	 * sequenceData vector.
	 * 
	 * The progressBar is turned on and off to indicate execution of thread
	 * 
	 */
	public void run() {
		// if there is another thread stop this one
		if (mainPanel.progressBar.isIndeterminate()) {
			mainPanel.statusBar
					.setText("Cannot run more than one funcion/model at the same time.");
			System.err
					.println("Cannot run more than one function and model at the same time.");
			return;
		}

		else if (mainPanel.sequenceData == null) {
			mainPanel.statusBar.setText("Select a sequence");
			System.err
					.println("Select sequence before executing function/model.");
			return;
		}

		// indicate task is being performed
		mainPanel.progressBar.setIndeterminate(true);

		try {
			// 1. Set parameters for function if function has options
			OptionsHandle options = myFunction.getOptionsHandle();
			if (myFunction instanceof ReadFormatFileFunction) {
				Iterator<SequenceData> iterSeq = ((ReadFormatFileFunction) myFunction)
						.guessFormat();
				attachMultiSequence(iterSeq);

			} else if (options == null || mainPanel.setParameters(options)) {

				// 2. Call method to execute function for inputData
				mainPanel.statusBar.setText("Running " + myFunction);
				SequenceData data = myFunction.execute(options, myInputData);
				if (data != null) {
					attachSequence(data);
				}
			}

		} catch (Exception e) {
			System.err.println("Error executing function " + myFunction);
			System.err.println(e);
			mainPanel.statusBar.setText("Error running function");
			mainPanel.progressBar.setIndeterminate(false);
			System.out.println(e.getMessage());
		}

		mainPanel.statusBar.setText("");
		mainPanel.progressBar.setIndeterminate(false);
	}

}
