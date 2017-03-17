/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 17 Apr 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsadev.seq.nanopore;

import japsa.util.DynamicHistogram;
import japsa.util.JapsaException;
import japsa.util.Logging;

import java.awt.EventQueue;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;

import javax.swing.JRadioButton;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.JCheckBox;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.time.Second;
import org.jfree.data.time.TimeTableXYDataset;

/**
 * @author minhduc
 *
 */
public class NanoporeReaderWindow implements Runnable{

	private JFrame frmNanoporeReader;
	private int height = 50;
	private int topR = 100, topC = 100;
	//String downloadFolder;	

	TimeTableXYDataset dataSet;
	NanoporeReaderStream reader;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					NanoporeReaderWindow window = new NanoporeReaderWindow(new NanoporeReaderStream(),null);
					window.frmNanoporeReader.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 * @throws IOException 
	 */
	public NanoporeReaderWindow(NanoporeReaderStream r, TimeTableXYDataset dataset) throws IOException {
		reader = r;
		this.dataSet = dataset;

		initialize();		
		//frmNanoporeReader.pack();
		frmNanoporeReader.setVisible(true);
	}			

	/**
	 * Initialize the contents of the frame.
	 * @throws IOException 
	 */
	private void initialize() throws IOException {		
		frmNanoporeReader = new JFrame();
		frmNanoporeReader.setTitle("Nanopore Reader");
		frmNanoporeReader.setBounds(topC, topR, 1360, 714);
		frmNanoporeReader.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmNanoporeReader.getContentPane().setLayout(new BorderLayout(0, 0));

		final JPanel controlPanel = new JPanel();
		frmNanoporeReader.getContentPane().add(controlPanel, BorderLayout.WEST);
		controlPanel.setPreferredSize(new Dimension(330, height));
		controlPanel.setLayout(null);

		final JPanel inputPanel = new JPanel();
		inputPanel.setBounds(0, 8, 320, 120);
		inputPanel.setBorder(BorderFactory.createTitledBorder("Input"));
		controlPanel.add(inputPanel);
		inputPanel.setLayout(null);

		final JTextField txtDir = new JTextField(reader.folder == null?"":reader.folder);		
		txtDir.setBounds(10, 51, 300, 20);
		inputPanel.add(txtDir);		

		//final ButtonGroup group = new ButtonGroup();

		final JButton btnChange = new JButton("Change");
		btnChange.setBounds(28, 83, 117, 25);
		inputPanel.add(btnChange);

		final JCheckBox chckbxInc = new JCheckBox("Include fail folder",reader.doFail);
		chckbxInc.setBounds(153, 84, 159, 23);
		inputPanel.add(chckbxInc);

		JLabel lblNewLabel = new JLabel("Folder containing base-called reads");
		lblNewLabel.setBounds(10, 24, 300, 15);
		inputPanel.add(lblNewLabel);

		chckbxInc.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				reader.doFail = (e.getStateChange() == ItemEvent.SELECTED);
			}		
		});

		btnChange.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setDialogTitle("Select download directory");
				fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				fileChooser.setCurrentDirectory(new File(txtDir.getText()));

				int returnValue = fileChooser.showOpenDialog(null);
				if (returnValue == JFileChooser.APPROVE_OPTION) {
					reader.folder = fileChooser.getSelectedFile().getPath();
					txtDir.setText(reader.folder);	            
				}
			}
		});


		final JPanel outputPanel = new JPanel();
		outputPanel.setBounds(0, 140, 320, 188);
		outputPanel.setBorder(BorderFactory.createTitledBorder("Output"));
		controlPanel.add(outputPanel);
		outputPanel.setLayout(null);

		final JRadioButton rdbtnOut2Str = new JRadioButton("Output to output stream");		
		rdbtnOut2Str.setBounds(10, 22, 302, 23);
		outputPanel.add(rdbtnOut2Str);

		final JRadioButton rdbtnOut2File = new JRadioButton("Output to file");	
		rdbtnOut2File.setBounds(10, 48, 302, 23);
		outputPanel.add(rdbtnOut2File);

		final ButtonGroup group2 = new ButtonGroup();
		group2.add(rdbtnOut2Str);		
		group2.add(rdbtnOut2File);

		final JTextField txtOFile = new JTextField(reader.output);		
		txtOFile.setBounds(10, 79, 300, 20);
		outputPanel.add(txtOFile);


		final JButton btnFileChange = new JButton("Change");		
		btnFileChange.setBounds(26, 105, 117, 25);	
		outputPanel.add(btnFileChange);


		final JCheckBox chckbxStreamServer = new JCheckBox("Stream output to some servers");
		final JTextField txtStreamServers = new JTextField();

		chckbxStreamServer.setSelected(reader.streamServers!=null);
		chckbxStreamServer.setBounds(10, 132, 283, 23);
		outputPanel.add(chckbxStreamServer);


		txtStreamServers.setText(reader.streamServers != null?reader.streamServers:"");
		txtStreamServers.setEnabled(chckbxStreamServer.isSelected());
		txtStreamServers.setBounds(10, 163, 300, 19);
		outputPanel.add(txtStreamServers);
		txtStreamServers.setColumns(10);

		chckbxStreamServer.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				txtStreamServers.setEnabled(chckbxStreamServer.isSelected());		
			}
		});


		btnFileChange.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setDialogTitle("Select output file");
				//fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				fileChooser.setCurrentDirectory(new File(txtOFile.getText()));
				int returnValue = fileChooser.showOpenDialog(null);
				if (returnValue == JFileChooser.APPROVE_OPTION) {
					reader.output = fileChooser.getSelectedFile().getPath();
					txtOFile.setText(reader.output);	            
				}
			}
		});


		if ("-".equals(reader.output)){
			rdbtnOut2Str.setSelected(true);
			rdbtnOut2File.setSelected(false);
			btnFileChange.setEnabled(false);
			txtOFile.setEnabled(false);
			txtOFile.setText("");			
		}else{
			rdbtnOut2Str.setSelected(false);
			rdbtnOut2File.setSelected(true);
			btnFileChange.setEnabled(true);
			txtOFile.setEnabled(true);
		}

		rdbtnOut2Str.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (rdbtnOut2Str.isSelected()){
					btnFileChange.setEnabled(false);
					txtOFile.setEnabled(false);
				}else{
					btnFileChange.setEnabled(true);
					txtOFile.setEnabled(true);
				}
			}
		});


		rdbtnOut2File.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (rdbtnOut2Str.isSelected()){
					btnFileChange.setEnabled(false);
					txtOFile.setEnabled(false);
				}else{
					btnFileChange.setEnabled(true);
					txtOFile.setEnabled(true);
				}
			}
		});

		JPanel formatPanel = new JPanel();
		formatPanel.setBorder(BorderFactory.createTitledBorder("Output format"));
		formatPanel.setBounds(0, 338, 320, 55);
		controlPanel.add(formatPanel);

		final JRadioButton fqRadioButton = new JRadioButton("fastq");
		fqRadioButton.setBounds(46, 22, 72, 23);

		final JRadioButton faRadioButton = new JRadioButton("fasta");
		faRadioButton.setBounds(186, 22, 72, 23);
		formatPanel.setLayout(null);

		final ButtonGroup formatBtGroup = new ButtonGroup();

		formatBtGroup.add(fqRadioButton);
		formatBtGroup.add(faRadioButton);

		formatPanel.add(fqRadioButton);
		formatPanel.add(faRadioButton);

		if ("fasta".equals(reader.format)){
			fqRadioButton.setSelected(false);
			faRadioButton.setSelected(true);
		}else{
			faRadioButton.setSelected(false);
			fqRadioButton.setSelected(true);
		}

		fqRadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (fqRadioButton.isSelected())
					reader.format = "fastq";
				else 
					reader.format = "fasta";
			}
		});

		faRadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (fqRadioButton.isSelected())
					reader.format = "fastq";
				else 
					reader.format = "fasta";
			}
		});		



		final JPanel optionPanel = new JPanel();
		optionPanel.setBounds(0, 405, 320, 115);
		optionPanel.setBorder(BorderFactory.createTitledBorder("Options"));
		controlPanel.add(optionPanel);
		optionPanel.setLayout(null);

		final JCheckBox chckReads = new JCheckBox("Include template and complement reads",true);
		chckReads.setBounds(7, 23, 310, 23);
		optionPanel.add(chckReads);

		chckReads.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				reader.doLow = (e.getStateChange() == ItemEvent.SELECTED);
			}		
		});


		final JCheckBox chckbxAddAUnicqu = new JCheckBox("Add a unique number to read name",reader.number);
		chckbxAddAUnicqu.setBounds(8, 52, 304, 23);
		optionPanel.add(chckbxAddAUnicqu);

		chckbxAddAUnicqu.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				reader.number = (e.getStateChange() == ItemEvent.SELECTED);
			}		
		});



		final JLabel lblMinReadLength = new JLabel("Min read length");
		lblMinReadLength.setBounds(7, 83, 154, 15);
		optionPanel.add(lblMinReadLength);

		final JTextField txtMinLenth = new JTextField();
		txtMinLenth.setText(reader.minLength+"");
		txtMinLenth.setBounds(137, 77, 80, 21);
		optionPanel.add(txtMinLenth);

		///// Barcode analysis control panel
		/*****************************************************************************************************
		 *****************************************************************************************************
		 *****************************************************************************************************
		 */
		JPanel bcPanel = new JPanel();
		bcPanel.setBorder(BorderFactory.createTitledBorder("Demultiplex"));
		bcPanel.setBounds(0, 521, 320, 120);
		controlPanel.add(bcPanel);

		final JRadioButton yRadioButton = new JRadioButton("Yes");
		yRadioButton.setBounds(46, 22, 62, 23);

		final JRadioButton nRadioButton = new JRadioButton("No");
		nRadioButton.setBounds(186, 22, 62, 23);
		bcPanel.setLayout(null);

		final ButtonGroup bcBtGroup = new ButtonGroup();

		bcBtGroup.add(yRadioButton);
		bcBtGroup.add(nRadioButton);

		bcPanel.add(yRadioButton);
		bcPanel.add(nRadioButton);

		final JTextField txtBCFile = new JTextField(reader.getBCFileName());		
		txtBCFile.setBounds(10, 50, 300, 20);
		bcPanel.add(txtBCFile);


		final JButton btnBCFileChange = new JButton("Change");		
		btnBCFileChange.setBounds(26, 80, 117, 25);	
		bcPanel.add(btnBCFileChange);
		
		if (reader.getBCFileName() == null){
			yRadioButton.setSelected(false);
			nRadioButton.setSelected(true);
			btnBCFileChange.setEnabled(false);
			txtBCFile.setEnabled(false);
		}else{
			yRadioButton.setSelected(true);
			nRadioButton.setSelected(false);
			btnBCFileChange.setEnabled(true);
			txtBCFile.setEnabled(true);
		}

		yRadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (nRadioButton.isSelected()){
					btnBCFileChange.setEnabled(false);
					txtBCFile.setEnabled(false);
				}else{
					btnBCFileChange.setEnabled(true);
					txtBCFile.setEnabled(true);
				}
			}
		});

		nRadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (nRadioButton.isSelected()){
					btnBCFileChange.setEnabled(false);
					txtBCFile.setEnabled(false);
				}else{
					btnBCFileChange.setEnabled(true);
					txtBCFile.setEnabled(true);
				}
			}
		});		
		
		btnBCFileChange.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				JFileChooser fileChooser = new JFileChooser();
				fileChooser.setDialogTitle("Select barcode file");
				fileChooser.setCurrentDirectory(new File(txtBCFile.getText()));

				int returnValue = fileChooser.showOpenDialog(null);
				if (returnValue == JFileChooser.APPROVE_OPTION) {
					reader.updateDemultiplexFile(fileChooser.getSelectedFile().getPath());
					txtBCFile.setText(reader.getBCFileName());	            
				}
			}
		});

		
		
		//final JLabel lblGroup = new JLabel("Group");
		//lblGroup.setBounds(7, 113, 154, 15);
		//optionPanel.add(lblGroup);

		//final JTextField txtGroup = new JTextField();
		//txtGroup.setText(reader.group);
		//txtGroup.setBounds(137, 107, 150, 21);
		//optionPanel.add(txtGroup);
			

//		final JPanel lPanel = new JPanel();
//		lPanel.setBounds(10, 641, 320, 55);
//		controlPanel.add(lPanel);
//		lPanel.setLayout(null);		
//
//		final JButton btnStart = new JButton("Start");
//		btnStart.setBounds(28, 18, 117, 25);
//		btnStart.setEnabled(true);
//		lPanel.add(btnStart);
//
//
//		final JButton btnStop = new JButton("Stop");		
//		btnStop.addActionListener(new ActionListener() {
//			public void actionPerformed(ActionEvent e) {
//				reader.wait = false;
//
//				while (!reader.done){
//					try {
//						Thread.sleep(100);
//					} catch (InterruptedException ee) {					
//						ee.printStackTrace();
//					}
//				}
//
//				stillRun = false;
//				JOptionPane.showMessageDialog(null, "Done", "Information", JOptionPane.PLAIN_MESSAGE);
//
//			}
//		});
//		btnStop.setBounds(191, 18, 117, 25);
//		btnStop.setEnabled(false);
//		lPanel.add(btnStop);



		final JPanel mainPanel = new JPanel();
		frmNanoporeReader.getContentPane().add(mainPanel, BorderLayout.CENTER);
		//mainPanel.setBorder(BorderFactory.createTitledBorder("Statistics"));
		mainPanel.setLayout(null);

		final JPanel panelCounts = new JPanel();
		panelCounts.setBounds(12, 304, 428, 280);
		mainPanel.add(panelCounts);
		panelCounts.setLayout(null);

		/**********************************************************************
		 * New position for start/stop buttons
		 */
		final JPanel lPanel = new JPanel();
		lPanel.setBounds(12, 580, 320, 55);
		mainPanel.add(lPanel);
		lPanel.setLayout(null);		

		final JButton btnStart = new JButton("Start");
		btnStart.setBounds(40, 5, 117, 25);
		btnStart.setEnabled(true);
		lPanel.add(btnStart);


		final JButton btnStop = new JButton("Stop");		
		btnStop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				reader.wait = false;

				while (!reader.done){
					try {
						Thread.sleep(100);
					} catch (InterruptedException ee) {					
						ee.printStackTrace();
					}
				}

				stillRun = false;
				JOptionPane.showMessageDialog(null, "Done", "Information", JOptionPane.PLAIN_MESSAGE);

			}
		});
		btnStop.setBounds(203, 5, 117, 25);
		btnStop.setEnabled(false);
		lPanel.add(btnStop);

		//////////////////////////////////////////////////
		final JLabel lblFiles = new JLabel("Total fast5 files");
		lblFiles.setBounds(10, 10, 140, 20);
		panelCounts.add(lblFiles);

		txtTFiles = new JTextField("0");
		txtTFiles.setEditable(false);
		txtTFiles.setBounds(150, 10, 110, 20);
		panelCounts.add(txtTFiles);
		txtTFiles.setColumns(10);

		////////////////////////////////////////////////
		final JLabel lblpFiles = new JLabel("Pass files");
		lblpFiles.setBounds(10, 35, 140, 20);
		//lblpFiles.setBounds(63, 61, 68, 15);
		panelCounts.add(lblpFiles);

		txtPFiles = new JTextField("0");
		txtPFiles.setEditable(false);
		txtPFiles.setBounds(150, 35, 110, 20);
		panelCounts.add(txtPFiles);
		txtPFiles.setColumns(10);

		final JLabel lblFFiles = new JLabel("Fail files");
		lblFFiles.setBounds(10, 60, 140, 20);
		panelCounts.add(lblFFiles);

		txtFFiles = new JTextField("0");
		txtFFiles.setEditable(false);
		txtFFiles.setBounds(150, 60, 110, 20);
		panelCounts.add(txtFFiles);
		txtFFiles.setColumns(10);


		final JLabel lbl2DReads = new JLabel("2D reads");
		lbl2DReads.setBounds(10, 90, 110, 20);
		panelCounts.add(lbl2DReads);		

		final JLabel lblTempReads = new JLabel("Template reads");
		lblTempReads.setBounds(10, 115, 140, 20);
		panelCounts.add(lblTempReads);	


		final JLabel lblCompReads = new JLabel("Complement reads");
		lblCompReads.setBounds(10, 140, 140, 20);
		panelCounts.add(lblCompReads);

		txtTempReads= new JTextField("0");
		txtTempReads.setEditable(false);
		txtTempReads.setBounds(150, 115, 110, 20);
		panelCounts.add(txtTempReads);
		txtTempReads.setColumns(10);		

		txt2DReads= new JTextField("0");
		txt2DReads.setEditable(false);
		txt2DReads.setBounds(150, 90, 110, 20);
		panelCounts.add(txt2DReads);
		txt2DReads.setColumns(10);		

		txtCompReads= new JTextField("0");
		txtCompReads.setEditable(false);
		txtCompReads.setBounds(150, 140, 110, 20);
		panelCounts.add(txtCompReads);
		txtCompReads.setColumns(10);


		final JFreeChart chart = ChartFactory.createStackedXYAreaChart(
			"Read count",      // chart title
			"Time",             // domain axis label
			"Read number",                   // range axis label
			this.dataSet   
			);			

		final StackedXYAreaRenderer render = new StackedXYAreaRenderer();

		DateAxis domainAxis = new DateAxis();
		domainAxis.setAutoRange(true);
		domainAxis.setDateFormatOverride(new SimpleDateFormat("HH:mm:ss"));

		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setRenderer(render);
		plot.setDomainAxis(domainAxis);
		plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
		plot.setForegroundAlpha(0.5f);

		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		rangeAxis.setNumberFormatOverride(new DecimalFormat("#,###.#"));
		rangeAxis.setAutoRange(true);

		ChartPanel chartPanel = new ChartPanel(chart,	            
			450,
			280,
			450,
			280,
			450,
			280,
			true,
			true,  // properties
			true,  // save
			true,  // print
			true,  // zoom
			true   // tooltips
			);		


		chartPanel.setBounds(0, 12, 450, 280);
		mainPanel.add(chartPanel);		

		/////////////////////////////////////////////////////////////
		//Histogram
		histoLengthDataSet=new DynamicHistogram();
		histoLengthDataSet.prepareSeries("Read Length", 500, 0, 40000);
		//histoDataset.prepareSeries("2D", 50, 0, 50000);
		//histoDataset.prepareSeries("template", 50, 0, 50000);
		//histoDataset.prepareSeries("complement", 50, 0, 50000);		

		JFreeChart hisLengths=ChartFactory.createHistogram("Read length histogram","length","count",histoLengthDataSet,PlotOrientation.VERTICAL,true,true,false);
		ChartPanel hisPanel = new ChartPanel(hisLengths,	            
			450,
			280,
			450,
			280,
			450,
			280,
			true,
			true,  // properties
			true,  // save
			true,  // print
			true,  // zoom
			true   // tooltips
			);


		XYPlot hisPlot = (XYPlot) hisLengths.getPlot();
		hisPlot.getDomainAxis().setAutoRange(true);		
		hisPlot.getRangeAxis().setAutoRange(true);


		hisPanel.setBounds(452, 12, 450, 280);
		mainPanel.add(hisPanel);


		histoQualDataSet=new DynamicHistogram();
		histoQualDataSet.setType(HistogramType.SCALE_AREA_TO_1);

		histoQualDataSet.prepareSeries("2D", 100, 0, 30);
		histoQualDataSet.prepareSeries("complement", 100, 0, 30);
		histoQualDataSet.prepareSeries("template", 100, 0, 30);

		JFreeChart hisQual=ChartFactory.createXYLineChart("Quality","quality","frequency",histoQualDataSet,PlotOrientation.VERTICAL,true,true,false);
		ChartPanel hisQualPanel = new ChartPanel(hisQual,	            
			450,
			280,
			450,
			280,
			450,
			280,
			true,
			true,  // properties
			true,  // save
			true,  // print
			true,  // zoom
			true   // tooltips
			);


		XYPlot hisQualPlot = (XYPlot) hisQual.getPlot();
		hisQualPlot.getDomainAxis().setAutoRange(true);		
		hisQualPlot.getRangeAxis().setAutoRange(true);

		hisQualPlot.setForegroundAlpha(0.8F);
		//XYBarRenderer xybarrenderer = (XYBarRenderer)hisQualPlot.getRenderer();
		//xybarrenderer.setDrawBarOutline(false);


		hisQualPanel.setBounds(452, 300, 450, 280);
		mainPanel.add(hisQualPanel);



		btnStart.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				//1. Validate before running
				
				//validate input

				String _path = txtDir.getText().trim();				
				if (_path.equals("")){
					JOptionPane.showMessageDialog(null, "Please specify download directory", "Warning", JOptionPane.PLAIN_MESSAGE);
					txtDir.grabFocus();
					return;
				}

				File _file = new File(_path);
				if (!_file.exists()){
					JOptionPane.showMessageDialog(null, "Directory \"" + _path + "\" does not exist!", "Warning", JOptionPane.PLAIN_MESSAGE);
					txtDir.grabFocus();
					return;
				}
				reader.folder = _path;
				//validate output
				if (rdbtnOut2File.isSelected()){
					String _foutput = txtOFile.getText().trim();
					if (_foutput.equals("")){						
						JOptionPane.showMessageDialog(null, "Please specify output file", "Warning", JOptionPane.PLAIN_MESSAGE);
						txtOFile.grabFocus();
						return;
					}
					reader.output = _foutput;					
				}else
					reader.output = "-";//stream
					
				
				//validate stream
				if (chckbxStreamServer.isSelected()){
					if (txtStreamServers.getText().trim().equals("")){
						JOptionPane.showMessageDialog(null, "Please specify output address of a server", "Warning", JOptionPane.PLAIN_MESSAGE);
						txtStreamServers.grabFocus();
						return;
					}			
					reader.streamServers = txtStreamServers.getText().trim();
				}
				
				//validate barcode analysis
				if(yRadioButton.isSelected()){
					if(txtBCFile.getText().trim().equals("")){
						JOptionPane.showMessageDialog(null, "Please specify barcode file for demultiplex", "Warning", JOptionPane.PLAIN_MESSAGE);
						txtBCFile.grabFocus();
						return;
					}
				}


				String msg = reader.prepareIO();
				if (msg !=null){
					JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.PLAIN_MESSAGE);
					return;
				}

				//Start running
				txtDir.setEnabled(false);
				btnChange.setEnabled(false);
				chckbxInc.setEnabled(false);
				rdbtnOut2Str.setEnabled(false);
				rdbtnOut2File.setEnabled(false);
				txtOFile.setEnabled(false);
				btnFileChange.setEnabled(false);
				chckReads.setEnabled(false);
				chckbxAddAUnicqu.setEnabled(false);
				txtMinLenth.setEnabled(false);
//				txtGroup.setEnabled(false);

				btnStart.setEnabled(false);
				btnStop.setEnabled(true);

				reader.ready = true;
			}
		});	}

	JTextField txtCompReads, txtTempReads, txt2DReads;
	JTextField txtPFiles, txtFFiles, txtTFiles;
	DynamicHistogram histoLengthDataSet, histoQualDataSet;

	boolean stillRun = true;

	public void interupt(JapsaException e){
		this.stillRun = false;
		reader.wait = false;
		JOptionPane.showMessageDialog(null, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
	}

	public void run() {		
		int lastIndexLengths = 0;//, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;

		while(stillRun) {				                
			//synchronized(reader) {//avoid concurrent update					
			Second period = new Second();
			dataSet.add(period, reader.twoDCount,"2D");
			dataSet.add(period, reader.compCount,"complement");
			dataSet.add(period, reader.tempCount,"template");

			txtTFiles.setText(reader.fileNumber+"");	                
			txtPFiles.setText(reader.passNumber+"");
			txtFFiles.setText(reader.failNumber+"");

			txt2DReads.setText(reader.twoDCount+"");
			txtCompReads.setText(reader.compCount+"");
			txtTempReads.setText(reader.tempCount+"");

			int currentIndex = reader.lengths.size();

			if (currentIndex > lastIndexLengths){			
				int index = histoLengthDataSet.getSeriesIndex("Read Length");
				for (int i = lastIndexLengths; i < currentIndex;i++)
					histoLengthDataSet.addSeries(index, reader.lengths.get(i));

				lastIndexLengths = currentIndex;

				histoLengthDataSet.notifyChanged();
			}

			currentIndex = reader.qual2D.size();
			if (currentIndex > lastIndexQual2D){
				int index = histoQualDataSet.getSeriesIndex("2D");
				for (int i = lastIndexQual2D; i < currentIndex;i++)
					histoQualDataSet.addSeries(index, reader.qual2D.get(i));

				lastIndexQual2D = currentIndex;
				histoQualDataSet.notifyChanged();
			}

			currentIndex = reader.qualComp.size();
			if (currentIndex > lastIndexQualComp){
				int index = histoQualDataSet.getSeriesIndex("complement");
				for (int i = lastIndexQualComp; i < currentIndex;i++)
					histoQualDataSet.addSeries(index, reader.qualComp.get(i));

				lastIndexQualComp = currentIndex;
				histoQualDataSet.notifyChanged();
			}

			currentIndex = reader.qualTemp.size();
			if (currentIndex > lastIndexQualTemp){
				int index = histoQualDataSet.getSeriesIndex("template");
				for (int i = lastIndexQualTemp; i < currentIndex;i++)
					histoQualDataSet.addSeries(index, reader.qualTemp.get(i));

				lastIndexQualTemp = currentIndex;
				histoQualDataSet.notifyChanged();
			}

			try {
				Thread.sleep(1000);
			} catch (InterruptedException ex) {
				Logging.error(ex.getMessage());
			}
		}
	}
}
