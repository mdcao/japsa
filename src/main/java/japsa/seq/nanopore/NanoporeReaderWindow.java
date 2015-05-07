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
package japsa.seq.nanopore;

import japsa.seq.SequenceOutputStream;
import japsa.util.Logging;

import java.awt.EventQueue;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

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
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.time.Second;
import org.jfree.data.time.TimeTableXYDataset;

import com.sun.corba.se.impl.encoding.CodeSetConversion.BTCConverter;

/**
 * @author minhduc
 *
 */
public class NanoporeReaderWindow implements Runnable
{
	private JFrame frmNanoporeReader;
	private int height = 600;
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
		frmNanoporeReader.setBounds(topC, topR, 1238, 714);
		frmNanoporeReader.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmNanoporeReader.getContentPane().setLayout(new BorderLayout(0, 0));

		final JPanel controlPanel = new JPanel();
		frmNanoporeReader.getContentPane().add(controlPanel, BorderLayout.WEST);
		controlPanel.setPreferredSize(new Dimension(320, height));
		controlPanel.setLayout(null);

		final JPanel inputPanel = new JPanel();
		inputPanel.setBounds(0, 23, 400, 157);
		inputPanel.setBorder(BorderFactory.createTitledBorder("Input"));
		controlPanel.add(inputPanel);
		inputPanel.setLayout(null);

		final JRadioButton rdbtnInputStream = new JRadioButton("Read files from input Stream");


		rdbtnInputStream.setBounds(8, 24, 384, 23);
		inputPanel.add(rdbtnInputStream);

		final JRadioButton rdbtnF = new JRadioButton("Read files from download directory");
		rdbtnF.setBounds(8, 51, 289, 23);
		inputPanel.add(rdbtnF);

		if (reader.folder == null){
			rdbtnInputStream.setSelected(true);			
		}else{
			rdbtnF.setSelected(true);
		}

		final ButtonGroup group = new ButtonGroup();
		group.add(rdbtnInputStream);
		group.add(rdbtnF);


		final JTextField txtDir = new JTextField(reader.folder);		
		txtDir.setBounds(18, 82, 289, 20);
		inputPanel.add(txtDir);
		txtDir.setEditable(false);
		txtDir.getDocument().addDocumentListener(new DocumentListener(){
			public void changedUpdate(DocumentEvent e){
			
			}
			public void removeUpdate(DocumentEvent e){
				reader.folder=txtDir.getText();
			}
			public void insertUpdate(DocumentEvent e){
				reader.folder=txtDir.getText();
			}
		}
		);
		
		
		final JButton btnChange = new JButton("Change");
		btnChange.setBounds(28, 109, 117, 25);
		inputPanel.add(btnChange);

		final JCheckBox chckbxInc = new JCheckBox("Include fail folder",reader.doFail);
		chckbxInc.setBounds(173, 110, 207, 23);
		inputPanel.add(chckbxInc);

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



		rdbtnF.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				if (e.getStateChange() == ItemEvent.SELECTED){
					txtDir.setEditable(true);
					btnChange.setEnabled(true);
					chckbxInc.setEnabled(true);
					reader.folder = txtDir.getText();

				}else{
					txtDir.setEditable(false);
					btnChange.setEnabled(true);
					chckbxInc.setEnabled(false);
					reader.folder = null;
				}
			}           
		});

		rdbtnInputStream.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (rdbtnInputStream.isSelected()){
					btnChange.setEnabled(false);
					chckbxInc.setEnabled(false);
					txtDir.setEditable(false);					
				}else{
					btnChange.setEnabled(true);
					chckbxInc.setEnabled(true);
					txtDir.setEditable(true);
				}
			}
		});


		final JPanel outputPanel = new JPanel();
		outputPanel.setBounds(0, 192, 400, 169);
		outputPanel.setBorder(BorderFactory.createTitledBorder("Output"));
		controlPanel.add(outputPanel);
		outputPanel.setLayout(null);

		final JRadioButton rdbtnOut2Str = new JRadioButton("Output to output stream");		
		rdbtnOut2Str.setBounds(8, 25, 302, 23);
		outputPanel.add(rdbtnOut2Str);

		final JRadioButton rdbtnOut2File = new JRadioButton("Output to file");	
		rdbtnOut2File.setBounds(8, 51, 302, 23);
		outputPanel.add(rdbtnOut2File);

		final ButtonGroup group2 = new ButtonGroup();
		group2.add(rdbtnOut2Str);		
		group2.add(rdbtnOut2File);

		final JTextField txtOFile = new JTextField(reader.output);		
		txtOFile.setBounds(18, 82, 295, 20);
		outputPanel.add(txtOFile);

		//txtOFile.setEditable(false);
		txtOFile.getDocument().addDocumentListener(new DocumentListener(){
			public void changedUpdate(DocumentEvent e){
			}
			public void removeUpdate(DocumentEvent e){
				reader.output=txtOFile.getText();
			}
			public void insertUpdate(DocumentEvent e){
				reader.output=txtOFile.getText();
			}
		}
		);
		

		final JButton btnFileChange = new JButton("Change");		
		btnFileChange.setBounds(26, 109, 117, 25);	
		outputPanel.add(btnFileChange);


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
			txtOFile.setEditable(false);
			txtOFile.setText("");			
		}else{
			rdbtnOut2Str.setSelected(false);
			rdbtnOut2File.setSelected(true);
			btnFileChange.setEnabled(true);
			txtOFile.setEditable(true);
		}

		rdbtnOut2Str.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (rdbtnOut2Str.isSelected()){
					btnFileChange.setEnabled(false);
					txtOFile.setEditable(false);
					reader.output="-";
					txtOFile.setText("-");
				}else{
					btnFileChange.setEnabled(true);
					txtOFile.setEditable(true);
				}
			}
		});


		rdbtnOut2File.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (rdbtnOut2Str.isSelected()){
					btnFileChange.setEnabled(false);
					txtOFile.setEditable(false);
					reader.output="-";
					txtOFile.setText("-");
				}else{
					btnFileChange.setEnabled(true);
					txtOFile.setEditable(true);
				}
			}
		});


		final JPanel optionPanel = new JPanel();
		optionPanel.setBounds(0, 373, 400, 132);
		optionPanel.setBorder(BorderFactory.createTitledBorder("Options"));
		controlPanel.add(optionPanel);
		optionPanel.setLayout(null);

		final JCheckBox chckReads = new JCheckBox("Include template and complement reads",true);
		chckReads.setBounds(8, 23, 310, 23);
		optionPanel.add(chckReads);

		chckReads.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				reader.doLow = (e.getStateChange() == ItemEvent.SELECTED);
			}		
		});


		final JCheckBox chckbxAddAUnicqu = new JCheckBox("Add a unique number to read name",reader.number);
		chckbxAddAUnicqu.setBounds(8, 52, 373, 23);
		optionPanel.add(chckbxAddAUnicqu);

		chckbxAddAUnicqu.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e){
				reader.number = (e.getStateChange() == ItemEvent.SELECTED);
			}		
		});



		final JLabel lblMinReadLength = new JLabel("Min read length");
		lblMinReadLength.setBounds(8, 83, 154, 15);
		optionPanel.add(lblMinReadLength);

		final JTextField txtMinLenth = new JTextField();
		txtMinLenth.setText(reader.minLength+"");
		txtMinLenth.setBounds(137, 77, 71, 21);
		optionPanel.add(txtMinLenth);

		final JPanel lPanel = new JPanel();
		lPanel.setBounds(0, 508, 400, 83);
		controlPanel.add(lPanel);
		lPanel.setLayout(null);		

		final JButton btnStart = new JButton("Start");
		btnStart.setBounds(27, 36, 117, 25);
		btnStart.setEnabled(true);
		lPanel.add(btnStart);


		final JButton btnStop = new JButton("Stop");		
		btnStop.setBounds(179, 36, 117, 25);
		btnStop.setEnabled(false);
		lPanel.add(btnStop);


		final JPanel mainPanel = new JPanel();
		frmNanoporeReader.getContentPane().add(mainPanel, BorderLayout.CENTER);
		//mainPanel.setBorder(BorderFactory.createTitledBorder("Statistics"));
		mainPanel.setLayout(null);

		final JPanel panelCounts = new JPanel();
		panelCounts.setBounds(12, 304, 428, 280);
		mainPanel.add(panelCounts);
		panelCounts.setLayout(null);


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
		histoLengthDataSet.prepareSeries("Read Length", 200, 0, 50000);
		//histoDataset.prepareSeries("2D", 50, 0, 50000);
		//histoDataset.prepareSeries("template", 50, 0, 50000);
		//histoDataset.prepareSeries("complement", 50, 0, 50000);		

		final JFreeChart hisLengths=ChartFactory.createHistogram("Read length histogram","length","count",histoLengthDataSet,PlotOrientation.VERTICAL,true,true,false);
		final ChartPanel hisLengthPanel = new ChartPanel(hisLengths,	            
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


		hisLengthPanel.setBounds(452, 12, 450, 280);
		mainPanel.add(hisLengthPanel);


		histoQualDataSet=new DynamicHistogram();
		histoQualDataSet.setType(HistogramType.SCALE_AREA_TO_1);

		histoQualDataSet.prepareSeries("2D", 100, 0, 30);
		histoQualDataSet.prepareSeries("complement", 100, 0, 30);
		histoQualDataSet.prepareSeries("template", 100, 0, 30);

		final JFreeChart hisQual=ChartFactory.createXYLineChart("Quality","quality","frequency",histoQualDataSet,PlotOrientation.VERTICAL,true,true,false);
		final ChartPanel hisQualPanel = new ChartPanel(hisQual,	            
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

				if(btnStart.getText().equals("Restart")){
					if(reader.sos != null && !reader.output.equals("-")){
						try{
							reader.sos.close();
						}catch(IOException e1){
							e1.printStackTrace();
						}
					}
					reader=reader.reset(); 
					
					// reset the read count graph
					chart.setNotify(false);
					dataSet.clear();
					chart.setNotify(true);
					// reset the read length histogram
					hisLengths.setNotify(false);
					histoLengthDataSet.reset();
					hisLengths.setNotify(true);
					hisLengthPanel.restoreAutoBounds();
					// reset the quality histogram
					hisQual.setNotify(false);
					histoQualDataSet.reset();
					hisQual.setNotify(true);
					hisQualPanel.restoreAutoBounds();
					// set the flag to run
				}
				// Validate before running
				if (rdbtnF.isSelected()){
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
				}


				if (rdbtnOut2File.isSelected()){
					if (txtOFile.getText().trim().equals("")){						
						JOptionPane.showMessageDialog(null, "Please specify output file", "Warning", JOptionPane.PLAIN_MESSAGE);
						txtOFile.grabFocus();
						return;
					}
				}

				String msg = reader.prepareIO();
				if (msg !=null){
					JOptionPane.showMessageDialog(null, msg, "Warning", JOptionPane.PLAIN_MESSAGE);
					return;
				}
				//Start running

				rdbtnInputStream.setEnabled(false);
				txtDir.setEnabled(false);
				btnChange.setEnabled(false);
				chckbxInc.setEnabled(false);
				rdbtnF.setEnabled(false);
				rdbtnOut2Str.setEnabled(false);
				rdbtnOut2File.setEnabled(false);
				txtOFile.setEnabled(false);
				btnFileChange.setEnabled(false);
				chckReads.setEnabled(false);
				chckbxAddAUnicqu.setEnabled(false);
				txtMinLenth.setEnabled(false);
				reader.ready=true;
				start();
				btnStart.setEnabled(false);
				btnStop.setEnabled(true);
				
			}
		});

		btnStop.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				//Stop running

				rdbtnInputStream.setEnabled(true);
				txtDir.setEnabled(true);
				btnChange.setEnabled(true);
				chckbxInc.setEnabled(true);
				rdbtnF.setEnabled(true);
				rdbtnOut2Str.setEnabled(true);
				rdbtnOut2File.setEnabled(true);
				txtOFile.setEnabled(true);
				btnFileChange.setEnabled(true);
				chckReads.setEnabled(true);
				chckbxAddAUnicqu.setEnabled(true);
				txtMinLenth.setEnabled(true);
				btnStop.setEnabled(false);
				btnStart.setText("Restart");
				btnStart.setEnabled(true);
			
				reader.ready=false;
				//frmNanoporeReader.repaint();
			}

		});
	}

	private void start(){
		SwingWorker<Void, Void> worker= new SwingWorker<Void, Void>(){
			@Override
			protected Void doInBackground() throws Exception{
				reader.readFastq(pFolderName);
				reader.sos.flush();
				return null;
			}
		
		};
		worker.execute();
	}

	JTextField txtCompReads, txtTempReads, txt2DReads;
	JTextField txtPFiles, txtFFiles, txtTFiles;
	DynamicHistogram histoLengthDataSet, histoQualDataSet;
	String pFolderName=null;
	//boolean stillRun = true;

	public void run() {		
		int lastIndexLengths = 0, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;

		while(true) {	
			if (reader.ready) {
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

			}
			else{
				lastIndexLengths = lastIndexLengths2D = lastIndexLengthsComp =  lastIndexLengthsTemp = 0;
				lastIndexQual2D = lastIndexQualComp = lastIndexQualTemp = 0;
			}

			try {
				Thread.sleep(1000);
			} catch (InterruptedException ex) {
				Logging.error(ex.getMessage());
			}
		}//while
	}//run
}

