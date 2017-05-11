/*****************************************************************************
 * Copyright (c) Son Hoang Nguyen, IMB - UQ, All rights reserved.         *
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
 * 01/03/2017 - Son Hoang Nguyen: Created                                        
 *  
 ****************************************************************************/
package japsa.seq.nanopore;

import japsa.util.DynamicHistogram;
import japsa.util.JapsaException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.StandardCategoryItemLabelGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.time.Second;
import org.jfree.data.time.TimeTableXYDataset;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.embed.swing.SwingNode;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.RowConstraints;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Stage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class NanoporeReaderWindowFX extends Application{
    private static final Logger LOG = LoggerFactory.getLogger(NanoporeReaderWindowFX.class);


    TimeTableXYDataset 	allReadsCount = new TimeTableXYDataset(),
						demultiplexedStackReadsCount = new TimeTableXYDataset();
	DefaultCategoryDataset demultiplexedBarReadsCount = new DefaultCategoryDataset();
	DynamicHistogram 	histoLengthDataSet = new DynamicHistogram(),
						histoQualDataSet = new DynamicHistogram();
	
	static NanoporeReaderStream reader = new NanoporeReaderStream();

	public static void setReader(NanoporeReaderStream r){
		reader = r;
	}
    
    public void start(Stage primaryStage){  	 
    	
        running(primaryStage);
    }
    
    private void running(Stage stage){
    	//Start with a BorderPane as root
	    BorderPane border = new BorderPane();
	    // Put start/stop button here  
	    HBox hbox = addHBox();
	    border.setTop(hbox);
	    
	    // All the parameters setting to the left 
	    leftBox=addVBox(stage);
	    border.setLeft(leftBox);
	        
	    // Add a stack to the HBox in the top region with Restart button.
	    // Uncomment this and line 332 (buttonRestart.setDisability(true)) 
	    // to have the function
        // addStackPane(hbox, stage);  
        
        // Here the main content    
        tabPane = new TabPane();
        Tab mainTab = new Tab("Main",addMainGridPane()),
        	bcTab = new Tab("Barcode",addBarcodePane());
        tabPane.getTabs().addAll(mainTab,bcTab);
        
        bcTab.disableProperty().bind(barcodeCB.selectedProperty().not());     
        border.setCenter(tabPane);
        
     
        Scene scene = new Scene(border);
        stage.setScene(scene);
        stage.setTitle("npreader");
        stage.setOnCloseRequest(e -> {
            Platform.exit();
            System.exit(0);
        });
	    stage.show();
			
		new Thread(new Runnable(){

			@Override
			public void run() {
				while (!reader.ready){
					//LOG.info("NOT READY");
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {					
						e.printStackTrace();
					}			
				}
				// TODO Auto-generated method stub
				LOG.info("GO");

				updateData();
				try{
					reader.readFast5();
				}catch (JapsaException e){
					System.err.println(e.getMessage());
					e.getStackTrace();
					interupt(e);
				}catch (Exception e){
					e.printStackTrace();
				}

			}
			
		}).start();
	     
    }
    
    private void reset(){
    	buttonStart.setDisable(false);
    	buttonStop.setDisable(true);
    	leftBox.setDisable(false);
    		
    	allReadsCount = new TimeTableXYDataset();
    	demultiplexedStackReadsCount = new TimeTableXYDataset();
		histoLengthDataSet = new DynamicHistogram(); 
		histoQualDataSet = new DynamicHistogram();
    	 
    	NanoporeReaderStream newReader = new NanoporeReaderStream();
		//reader.getTime = time;
    	newReader.number = reader.number;
    	newReader.minLength = reader.minLength;
    	newReader.interval = reader.interval;
    	newReader.age = reader.age;
    	newReader.folder = reader.folder;		
    	newReader.doFail = reader.doFail;
    	newReader.output = reader.output;
    	newReader.format = reader.format.toLowerCase();
    	newReader.realtime = reader.realtime;
    	newReader.streamServers = reader.streamServers;
    	newReader.exhaustive = reader.exhaustive;
    	newReader.realtime = true;
		newReader.stats = true;//GUI implies stats
		newReader.ready = false;//wait for the command from GUI
		newReader.updateDemultiplexFile(reader.getBCFileName());
    	
		setReader(newReader);
		
		txtCompReads.setText("0");
		txtTempReads.setText("0");
		txt2DReads.setText("0");
		txtPFiles.setText("0");
		txtFFiles.setText("0"); 
		txtTFiles.setText("0");
    	
    }
    
    private void restart(Stage stage){
    	reset();
    	running(stage);
    }
    /*
     * Components from left pane
     */
    private VBox leftBox;
    private TabPane tabPane; 
    private Button 	buttonStart, buttonStop, buttonRestart,
    				inputBrowseButton, barcodeBrowseButton, outputBrowseButton;
    private TextField inputTF, barcodeTF, bcThresholdTF, outputTF, streamTF, minLenTF;
    private CheckBox failCB, exhautiveCB, barcodeCB, serversCB, saveDemultiplexToFilesOptCB, addNumberOptCB;
    private ComboBox<String> outputToCombo, outputFormatCombo;
    /*
     * Creates an HBox with two buttons for the top region
     */
        
    private HBox addHBox() {
 
        HBox hbox = new HBox();
        hbox.setPadding(new Insets(15, 12, 15, 12));
        hbox.setSpacing(10);   // Gap between nodes
        hbox.setStyle("-fx-background-color: #336699;");
 

        Image imageStart = new Image(getClass().getResourceAsStream("icons/start.png"));
        ImageView viewStart = new ImageView(imageStart); 
        viewStart.setFitWidth(20);
        viewStart.setFitHeight(20);
        buttonStart = new Button("Start", viewStart);
        buttonStart.setPrefSize(100, 20);
        buttonStart.setOnAction((event) -> {
			//1. Validate before running	
			//validate input
			String _path = inputTF.getText().trim();				
			if (_path.equals("")){
				FxDialogs.showWarning("File not found!", "Please specify download directory");
				inputTF.requestFocus();
				return;
			}

			File _file = new File(_path);
			if (!_file.isDirectory()){
				FxDialogs.showWarning("File not found!", "Directory \"" + _path + "\" does not exist!");
				inputTF.requestFocus();
				return;
			}
			reader.folder = _path;
			//validate output
			if (outputToCombo.getSelectionModel().getSelectedItem().toString().equals("to file")){
				String _foutput = outputTF.getText().trim();
				if (_foutput.equals("")){		
					FxDialogs.showWarning("File not found!", "Please specify output file");
					outputTF.requestFocus();
					return;
				} else if(new File(_foutput).exists()){
					String confirm = FxDialogs.showConfirm( "Output file already exists!", "Are you sure to overwrite the old file?", "No", "Yes");
					if(confirm.equals("No")){
						outputTF.requestFocus();
						return;
					}
				}
				reader.output = new File(_foutput).getAbsolutePath();			
				try{
					System.setProperty("usr.dir", Paths.get(reader.output).getParent().toString());
				}
				catch(NullPointerException | IllegalArgumentException | SecurityException e ){
					e.printStackTrace();
					FxDialogs.showWarning("Illegal output folder!", "Please specify another output destination");
					outputTF.requestFocus();
					return;
				}
			}else
				reader.output = "-";//stream
				
			
			//validate stream
			if (serversCB.isSelected()){
				if (streamTF.getText().trim().equals("")){
					FxDialogs.showWarning("Server(s) not found!", "Please specify output address of a server");
					streamTF.requestFocus();
					return;
				}			
				reader.streamServers = streamTF.getText().trim();
			}
			
			//validate barcode analysis
			if(barcodeCB.isSelected()){
				if(barcodeTF.getText().trim().equals("")){
					FxDialogs.showWarning("File not found!", "Please specify barcode file for demultiplex");
					barcodeTF.requestFocus();
					return;
				}
				reader.updateDemultiplexFile(barcodeTF.getText().trim());
			}
			try{
				reader.minLength = Integer.parseInt(minLenTF.getText().trim());
			}catch(NumberFormatException e){
				reader.minLength = 0;
				minLenTF.setText("0");
			}
			
			if(reader.dmplx!=null)
				try{
					reader.dmplx.setThreshold(Integer.parseInt(bcThresholdTF.getText().trim()));
				}catch(NumberFormatException e){
					bcThresholdTF.setText(Integer.toString(reader.dmplx.SCORE_THRES));
				}
			
			String msg = reader.prepareIO();
			if (msg !=null){
				FxDialogs.showWarning("Warning", msg);
				return;
			}

			//Start running
			leftBox.setDisable(true);

			buttonStart.setDisable(true);;
			buttonStop.setDisable(false);

			reader.ready = true;
		});
        
        Image imageStop = new Image(getClass().getResourceAsStream("icons/stop.png"));
        ImageView viewStop = new ImageView(imageStop); 
        viewStop.setFitWidth(20);
        viewStop.setFitHeight(20);
        buttonStop = new Button("Stop", viewStop);
        buttonStop.setPrefSize(100, 20);
        buttonStop.setDisable(true);
        buttonStop.setOnAction((event) -> {
			String confirm = FxDialogs.showConfirm( "STOP button just being hit...", "Do you really want to stop the process?", "No", "Yes");
			if(confirm.equals("No")){
				return;
			}
			
        	reader.wait = false;
        	buttonStop.setDisable(true);

			try {
				reader.flush();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			if(!reader.output.equals("-")){
				try {
					reader.close();
				} catch (IOException e){
					e.printStackTrace();
				}
			}

        	//buttonRestart.setDisable(false);

        });
        
        
        hbox.getChildren().addAll(buttonStart, buttonStop);
        
        return hbox;
    }
        
    private final int LeftPaneWidth=360;
    /*
     * Creates a VBox with a list of parameter settings
     */
    private VBox addVBox(Stage stage) {
        
        VBox vbox = new VBox();
        vbox.setPadding(new Insets(10)); // Set all sides to 10
        vbox.setSpacing(8);              // Gap between nodes
 
        final Text title = new Text("Settings");
        title.setFont(Font.font("Arial", FontWeight.BOLD, 15));
        vbox.getChildren().add(title);
        final Separator sep1 = new Separator();
        sep1.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(1, sep1);
        
        vbox.getChildren().add(addInputPane(stage));
        vbox.setSpacing(5);
        final Separator sep2 = new Separator();
        sep2.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(3, sep2);

        
        vbox.getChildren().add(addOutputPane(stage));
        vbox.setSpacing(5);
        final Separator sep3 = new Separator();
        sep3.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(5, sep3);
        
        vbox.getChildren().add(addOptionPane());
               
        
        return vbox;
    }
    private GridPane addInputPane(Stage stage) {
    	GridPane inputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label inputLabel = new Label("Input:");
    	inputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	inputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(inputLabel, 0,0);
    	inputPane.getChildren().add(inputLabel);
    	
    	failCB = new CheckBox("Include fail folder");
    	failCB.setSelected(reader.doFail);
    	failCB.selectedProperty().addListener(
            (obs_val,old_val,new_val) -> {
            	reader.doFail = new_val;
            });	
    	
    	GridPane.setConstraints(failCB, 2,0,2,1);
    	inputPane.getChildren().add(failCB);
    	
	
    	
    	inputTF = new TextField(reader.folder == null?"":reader.folder);
    	inputTF.setPromptText("Enter folder of basecalled reads...");
    	inputTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	//textField.setPrefWidth(250);
    	GridPane.setConstraints(inputTF, 0,1,4,1);
    	inputPane.getChildren().add(inputTF);
    	

    	inputBrowseButton = new ImageButton("icons/folder.png");
    	inputBrowseButton.setPrefSize(10, 10);
    	inputBrowseButton.setOnAction((event) -> {
    		DirectoryChooser chooser = new DirectoryChooser();
    		chooser.setTitle("Select basecalled raw data (fast5) directory");
    		File defaultDirectory = new File(inputTF.getText());
    		if(defaultDirectory.isDirectory())
    			chooser.setInitialDirectory(defaultDirectory);
    		File selectedDirectory = chooser.showDialog(stage);
    		if(selectedDirectory != null){
				reader.folder = selectedDirectory.getPath();
				inputTF.setText(reader.folder);	
    		}
        });
    	GridPane.setConstraints(inputBrowseButton, 4,1);
    	GridPane.setHalignment(inputBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(inputBrowseButton);
    	//inputPane.setGridLinesVisible(true);

    	barcodeCB = new CheckBox("Demultiplexing for barcode analysis");
    	barcodeCB.setSelected(reader.dmplx!=null);
    	GridPane.setConstraints(barcodeCB, 0,3,5,1);
    	inputPane.getChildren().add(barcodeCB);
    	
    	barcodeTF = new TextField(reader.getBCFileName() == null?"":reader.getBCFileName());
    	barcodeTF.setPromptText("Enter name of barcode sequences file...");
    	barcodeTF.setDisable(!barcodeCB.isSelected());
    	barcodeTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(barcodeTF, 0,4,4,1);
    	inputPane.getChildren().add(barcodeTF);
    	
    	barcodeBrowseButton = new ImageButton("icons/folder.png");
    	barcodeBrowseButton.setPrefSize(10, 10);
    	barcodeBrowseButton.setDisable(!barcodeCB.isSelected());
    	barcodeBrowseButton.setOnAction((event) -> {
    		FileChooser chooser = new FileChooser();
    		chooser.setTitle("Select barcode file");
    		File defaultFile = new File(barcodeTF.getText());
    		if(defaultFile.isFile())
    			chooser.setInitialFileName(defaultFile.getName());
    		chooser.setInitialDirectory(defaultFile.getParentFile());
    		chooser.setSelectedExtensionFilter(new ExtensionFilter("FASTA files", "*.fasta", "*.fna", "*.fa"));
    		File selectedFile = chooser.showOpenDialog(stage);
    		if(selectedFile != null){
				reader.updateDemultiplexFile(selectedFile.getPath());
				barcodeTF.setText(reader.getBCFileName());	
    		}
        });
    	
    	final Label label = new Label("Barcode matching threshold");
    	label.setDisable(!barcodeCB.isSelected());;
    	GridPane.setConstraints(label, 0,5,3,1);
    	inputPane.getChildren().add(label);
    	
    	bcThresholdTF = new TextField(reader.dmplx!=null?Integer.toString(reader.dmplx.SCORE_THRES):"");
    	bcThresholdTF.setPromptText("Enter minimum score...");
    	barcodeBrowseButton.setDisable(!barcodeCB.isSelected());
    	bcThresholdTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(bcThresholdTF, 3,5);
    	inputPane.getChildren().add(bcThresholdTF);
    	
    	GridPane.setConstraints(barcodeBrowseButton, 4,4);
    	GridPane.setHalignment(barcodeBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(barcodeBrowseButton);
    	
    	barcodeCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	if(!new_val)
                		reader.dmplx = null;
                	barcodeTF.setDisable(!new_val);
                	barcodeBrowseButton.setDisable(!new_val);
                	bcThresholdTF.setDisable(!new_val);
                	saveDemultiplexToFilesOptCB.setDisable(!new_val);
                	label.setDisable(!new_val);
                	tabPane.getSelectionModel().select(0);
                });	
    	
		return inputPane;
	}
    private GridPane addOutputPane(Stage stage) {
    	GridPane outputPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label outputLabel = new Label("Output:");
    	outputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	outputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(outputLabel, 0,0);
    	outputPane.getChildren().add(outputLabel);
    	
        outputToCombo = new ComboBox<String>();
        outputToCombo.getItems().addAll("to file", "to stdout");   
        outputToCombo.setValue(reader.output.equals("-")?"to stdout":"to file");
        outputToCombo.valueProperty().addListener((obs_val, old_val, new_val) -> {
        	if(new_val.trim().equals("to file")){
        		//outputTF.setText("");
        		outputTF.setDisable(false);
        		outputBrowseButton.setDisable(false);
        	} else{
        		outputTF.setText("-");
        		outputTF.setDisable(true);
        		outputBrowseButton.setDisable(true);
        	}
        });
        GridPane.setConstraints(outputToCombo, 1, 0, 2, 1);
        outputPane.getChildren().add(outputToCombo);
        
        outputFormatCombo = new ComboBox<String>();
        outputFormatCombo.getItems().addAll("fastq", "fasta");   
        outputFormatCombo.setValue(reader.format);
        outputFormatCombo.valueProperty().addListener((obs_val, old_val, new_val) -> {
        	reader.format = new_val;
        });
        GridPane.setConstraints(outputFormatCombo, 3, 0, 2, 1);
        outputPane.getChildren().add(outputFormatCombo);
    	
    	outputTF = new TextField(reader.output);
    	outputTF.setDisable(reader.output.equals("-"));
    	outputTF.setPromptText("Enter name for output file...");
    	outputTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(outputTF, 0,1,4,1);
    	outputPane.getChildren().add(outputTF);
    	

    	outputBrowseButton = new ImageButton("icons/folder.png");
    	outputBrowseButton.setPrefSize(10, 10);
    	outputBrowseButton.setDisable(reader.output.equals("-"));
    	outputBrowseButton.setOnAction((event) -> {
    		FileChooser fileChooser = new FileChooser();
    		fileChooser.setTitle("Save output to file");
    		File initFolder = new File(inputTF.getText());
    		if(initFolder.isDirectory())
    			fileChooser.setInitialDirectory(initFolder);
    		fileChooser.setInitialFileName("output."+reader.format);
    		File savedFile = fileChooser.showSaveDialog(stage);
    		if(savedFile != null){
    			reader.output = savedFile.getAbsolutePath();
    			outputTF.setText(reader.output);
    		}
        });
    	GridPane.setConstraints(outputBrowseButton, 4,1);
    	GridPane.setHalignment(outputBrowseButton, HPos.LEFT);
    	outputPane.getChildren().add(outputBrowseButton);

    	//init
//        if(reader.output.equals("-")){
//        	outputTF.setText("-");
//        	outputToCombo.setValue("to stdout");
//    		outputTF.setDisable(true);
//    		outputBrowseButton.setDisable(true);
//        }else{
//        	outputTF.setText(reader.output);
//        	outputToCombo.setValue("to file");
//    		outputTF.setDisable(false);
//    		outputBrowseButton.setDisable(false);
//        }
        
    	
    	serversCB = new CheckBox("Streaming output to server(s)");
    	serversCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	streamTF.setDisable(!new_val);
                });	
    	GridPane.setConstraints(serversCB, 0,4,3,1);
    	outputPane.getChildren().add(serversCB);
    	
    	streamTF = new TextField();
    	streamTF.setPromptText("address1:port1, address2:port2,...");
    	streamTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(streamTF, 0,5,4,1);
    	outputPane.getChildren().add(streamTF);
    	
    	if(reader.streamServers != null){
    		serversCB.setSelected(true);
    		streamTF.setText(reader.streamServers);
    	}else{
    		streamTF.setDisable(true);
    	}
    	//outputPane.setGridLinesVisible(true);
		return outputPane;
	}
    private GridPane addOptionPane() {
    	GridPane optionPane = createFixGridPane(LeftPaneWidth, 5);
    	
    	final Label optLabel = new Label("Other options:");
    	optLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	optLabel.setStyle("-fx-underline:true");
    	
    	GridPane.setConstraints(optLabel, 0,0,4,1);
    	optionPane.getChildren().add(optLabel);
    	
    	saveDemultiplexToFilesOptCB = new CheckBox("Save demultiplexed reads to separated files");
    	saveDemultiplexToFilesOptCB.setSelected(Demultiplexer.toPrint);
    	saveDemultiplexToFilesOptCB.setDisable(!barcodeCB.isSelected());;
    	saveDemultiplexToFilesOptCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	Demultiplexer.toPrint=new_val;
                });	
    	GridPane.setConstraints(saveDemultiplexToFilesOptCB, 0,2,4,1);
    	optionPane.getChildren().add(saveDemultiplexToFilesOptCB);

    	addNumberOptCB = new CheckBox("Assign unique number to every read name");
    	addNumberOptCB.setSelected(reader.number);
    	addNumberOptCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	reader.number=new_val;
                });	
    	GridPane.setConstraints(addNumberOptCB, 0,4,4,1);
    	optionPane.getChildren().add(addNumberOptCB);
    	
    	exhautiveCB = new CheckBox("Exhaustively watch-mode (Albacore)");
    	exhautiveCB.setSelected(reader.exhaustive);
    	exhautiveCB.selectedProperty().addListener(
            (obs_val,old_val,new_val) -> {
            	reader.exhaustive = new_val;
            });	
    	
    	GridPane.setConstraints(exhautiveCB, 0,6,4,1);
    	optionPane.getChildren().add(exhautiveCB);
    	
    	final Label label2 = new Label("Filter out read shorter than ");
    	GridPane.setConstraints(label2, 0,8,3,1);
    	optionPane.getChildren().add(label2);
    	
    	minLenTF = new TextField(Integer.toString(reader.minLength));
    	minLenTF.setPromptText("min.");
    	minLenTF.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.ENTER)  {
                buttonStart.requestFocus();
            }
    	});
    	GridPane.setConstraints(minLenTF, 3,8);
    	optionPane.getChildren().add(minLenTF);
    	
    	final Label label3 = new Label("bp");
    	GridPane.setConstraints(label3, 4,8);
    	optionPane.getChildren().add(label3);
    	
		return optionPane;
	}
    
    private GridPane createFixGridPane(int width, int ncols){
        GridPane gridpane = new GridPane();
        for (int i = 0; i < ncols; i++) {
            ColumnConstraints column = new ColumnConstraints(1.0*width/ncols);
            gridpane.getColumnConstraints().add(column);
        }
        gridpane.setPadding(new Insets(10, 10, 10, 10));
        gridpane.setVgap(5);
        gridpane.setHgap(5);
        return gridpane;
    }
    private GridPane createAutoresizeGridPane(int ncols, int nrows){
        GridPane gridpane = new GridPane();
        for (int i = 0; i < ncols; i++) {
            ColumnConstraints column = new ColumnConstraints();
            column.setPercentWidth(100/ncols);
            gridpane.getColumnConstraints().add(column);
        }
        for (int i = 0; i < nrows; i++) {
            RowConstraints row = new RowConstraints();
            row.setPercentHeight(100/nrows);
            gridpane.getRowConstraints().add(row);
        }
        gridpane.setPadding(new Insets(5, 5, 5, 5));
        gridpane.setVgap(5);
        gridpane.setHgap(5);
        return gridpane;
    }
	/*
     * Restart button locates here
     * @param hb HBox to add the stack to
     */
    private void addStackPane(HBox hb, Stage stage) {
 
        StackPane stack = new StackPane();
        
        Image imageStart = new Image(getClass().getResourceAsStream("icons/restart.png"));
        ImageView viewStart = new ImageView(imageStart); 
        viewStart.setFitWidth(20);
        viewStart.setFitHeight(20);
        buttonRestart = new Button("Restart", viewStart);
        buttonRestart.setOnAction(e ->{
        	restart(stage);
        });
        buttonRestart.setDisable(true);
        
        stack.getChildren().add(buttonRestart);
        stack.setAlignment(Pos.CENTER_RIGHT);
        // Add offset to right for question mark to compensate for RIGHT 
        // alignment of all nodes
        //StackPane.setMargin(buttonRestart, new Insets(0, 10, 0, 0));
        
        hb.getChildren().add(stack);
        HBox.setHgrow(stack, Priority.ALWAYS);
                
    }
     
    /*
     * Creates a grid for the center region with 2 columns and 2 rows
     */
    private GridPane addMainGridPane() {
 
        GridPane mainGrid = createAutoresizeGridPane(2,2);
        mainGrid.setStyle("-fx-background-color: #C0C0C0;");
        /*
         * Read count chart
         */
		final JFreeChart chart = ChartFactory.createStackedXYAreaChart(
				"Read count",      // chart title
				"Time",             // domain axis label
				"Read number",                   // range axis label
				allReadsCount   
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

		SwingNode chartSwingNode = new SwingNode();
		chartSwingNode.setContent(chartPanel);
		GridPane.setConstraints(chartSwingNode, 0,0);
		
		mainGrid.getChildren().add(chartSwingNode);	
			
		/*
		 * Read length histogram	
		 */
		
		//histoLengthDataSet=new DynamicHistogram();
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

		SwingNode lengthSwingNode = new SwingNode();
		lengthSwingNode.setContent(hisPanel);
		GridPane.setConstraints(lengthSwingNode, 1,0);
//		GridPane.setHalignment(lengthSwingNode, HPos.CENTER);
//		GridPane.setValignment(lengthSwingNode, VPos.CENTER);
		mainGrid.getChildren().add(lengthSwingNode);
		
		/*
		 * Quality histogram
		 */
		//histoQualDataSet=new DynamicHistogram();
		histoQualDataSet.setType(HistogramType.SCALE_AREA_TO_1);
		histoQualDataSet.prepareSeries("2D", 100, 0, 30);
		histoQualDataSet.prepareSeries("template", 100, 0, 30);
		histoQualDataSet.prepareSeries("complement", 100, 0, 30);


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
		
		SwingNode qualitySwingNode = new SwingNode();
		qualitySwingNode.setContent(hisQualPanel);
		GridPane.setConstraints(qualitySwingNode, 1,1);
		mainGrid.getChildren().add(qualitySwingNode);
       
		/*
         * Statistics field
         */
		GridPane countPane = createAutoresizeGridPane(3, 9);
		countPane.setPadding(new Insets(30, 30, 30, 30));
        countPane.setStyle("-fx-background-color: #AABBCC;");

        
        final Label lblFiles = new Label("Total fast5 files");
		GridPane.setConstraints(lblFiles, 0, 0);
		countPane.getChildren().add(lblFiles);

		//txtTFiles = new TextField("0");
		txtTFiles.setPrefWidth(100);
		GridPane.setConstraints(txtTFiles, 1, 0);
		countPane.getChildren().add(txtTFiles);

		final Label lblpFiles = new Label("Good-read fast5");
		GridPane.setConstraints(lblpFiles, 0, 1);
		countPane.getChildren().add(lblpFiles);

		//txtPFiles = new TextField("0");
		txtPFiles.setEditable(false);
		txtPFiles.setPrefWidth(100);
		GridPane.setConstraints(txtPFiles, 1, 1);
		countPane.getChildren().add(txtPFiles);


		final Label lblFFiles = new Label("Invalid fast5");
		GridPane.setConstraints(lblFFiles, 0, 2);
		countPane.getChildren().add(lblFFiles);

		//txtFFiles = new TextField("0");
		txtFFiles.setEditable(false);
		txtFFiles.setPrefWidth(100);
		GridPane.setConstraints(txtFFiles, 1, 2);
		countPane.getChildren().add(txtFFiles);



		final Label lbl2DReads = new Label("2D reads");
		GridPane.setConstraints(lbl2DReads, 0, 3);
		countPane.getChildren().add(lbl2DReads);		

		//txt2DReads= new TextField("0");
		txt2DReads.setEditable(false);
		txt2DReads.setPrefWidth(100);
		GridPane.setConstraints(txt2DReads, 1, 3);
		countPane.getChildren().add(txt2DReads);

		final Label lblTempReads = new Label("Template reads");
		GridPane.setConstraints(lblTempReads, 0, 4);
		countPane.getChildren().add(lblTempReads);

		//txtTempReads= new TextField("0");
		txtTempReads.setEditable(false);
		txtTempReads.setPrefWidth(100);
		GridPane.setConstraints(txtTempReads, 1, 4);
		countPane.getChildren().add(txtTempReads);

		final Label lblCompReads = new Label("Complement reads");
		GridPane.setConstraints(lblCompReads, 0, 5);
		countPane.getChildren().add(lblCompReads);
		
		//txtCompReads= new TextField("0");
		txtCompReads.setEditable(false);
		txtCompReads.setPrefWidth(100);
		GridPane.setConstraints(txtCompReads, 1, 5);
        countPane.getChildren().add(txtCompReads);

        
        GridPane.setConstraints(countPane, 0, 1);
        mainGrid.getChildren().add(countPane);
        
//        mainGrid.setGridLinesVisible(true);
        return mainGrid;
    }
    
    /*
     * Create a Grid Pane for barcode analysis
     */
    private GridPane addBarcodePane(){
	   GridPane mainGrid = createAutoresizeGridPane(1,2);
       mainGrid.setStyle("-fx-background-color: #C0C0C0;");
       mainGrid.setPadding(new Insets(5, 100, 5, 100));
       mainGrid.setVgap(5);
       mainGrid.setHgap(5);
       /*
        * Read count stack chart
        */
   		final JFreeChart stackChart = ChartFactory.createXYLineChart(
   				"",      // chart title
   				"Time",             // domain axis label
   				"Over-time read count",                   // range axis label
   				demultiplexedStackReadsCount,
   				PlotOrientation.VERTICAL,
   				true,
   				true,
   				false
   				);			

   		final StackedXYAreaRenderer stackRender = new StackedXYAreaRenderer();		
   				
   		DateAxis domainAxis = new DateAxis();
   		domainAxis.setAutoRange(true);
   		domainAxis.setDateFormatOverride(new SimpleDateFormat("HH:mm:ss"));

   		XYPlot plot = (XYPlot) stackChart.getPlot();
   		plot.setRenderer(stackRender);
   		plot.setDomainAxis(domainAxis);
   		plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
   		plot.setForegroundAlpha(0.5f);

   		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
   		rangeAxis.setNumberFormatOverride(new DecimalFormat("#,###.#"));
   		rangeAxis.setAutoRange(true);

   		ChartPanel stackChartPanel = new ChartPanel(stackChart,	            
   			500,
   			300,
   			500,
   			300,
   			500,
   			300,
   			true,
   			true,  // properties
   			true,  // save
   			true,  // print
   			true,  // zoom
   			true   // tooltips
   			);		

   		SwingNode stackChartSwingNode = new SwingNode();
   		stackChartSwingNode.setContent(stackChartPanel);
   		
   		GridPane.setConstraints(stackChartSwingNode, 0,0);
   		mainGrid.getChildren().add(stackChartSwingNode);	
   		
   		/*
   		 * Reads count bar chart
   		 */
   		final JFreeChart barChart = ChartFactory.createBarChart(
   				"", 
   				"Barcode", 
   				"In-time read count", 
   				demultiplexedBarReadsCount, 
   				PlotOrientation.VERTICAL, 
   				false, 
   				false, 
   				false);
   		final BarRenderer barRender = new CustomRenderer(ChartColor.createDefaultPaintArray());
   		barRender.setSeriesItemLabelGenerator(0, new StandardCategoryItemLabelGenerator());
   		barRender.setSeriesItemLabelsVisible(0, true);
//   	barRender.setBarPainter(new StandardBarPainter());
   		
   		CategoryPlot bar = barChart.getCategoryPlot();
   		bar.getDomainAxis().setVisible(false);
   		bar.setRenderer(barRender);
   		bar.setForegroundAlpha(0.5f);


   		ChartPanel barChartPanel = new ChartPanel(barChart,	            
   			500,
   			300,
   			500,
   			300,
   			500,
   			300,
   			true,
   			true,  // properties
   			true,  // save
   			true,  // print
   			true,  // zoom
   			true   // tooltips
   			);		

   		SwingNode barChartSwingNode = new SwingNode();
   		barChartSwingNode.setContent(barChartPanel);
   		
   		GridPane.setConstraints(barChartSwingNode, 0,1);
   		mainGrid.getChildren().add(barChartSwingNode);
   		
   		return mainGrid;
    }
    
    /******************************************************************************************
     * ** Here are variables and controls for all plots ***************************************
     ******************************************************************************************/
	final TextField txtCompReads= new TextField("0"), 
					txtTempReads= new TextField("0"), 
					txt2DReads= new TextField("0");
	final TextField txtPFiles= new TextField("0"), 
					txtFFiles= new TextField("0"), 
					txtTFiles= new TextField("0");

	//private static boolean stillRun = true;

	public static void interupt(JapsaException e){
		//stillRun = false;
		reader.wait = false;
		FxDialogs.showError("Unexpected errors happened!", e.getMessage());
	}
	
	private void updateData(){
		 
        final ScheduledExecutorService scheduler 
            = Executors.newScheduledThreadPool(1);
 
        scheduler.scheduleAtFixedRate(
                new Runnable(){         
            		int lastIndexLengths = 0;//, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
            		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;
                    @Override
                    public void run() {
                    	if(!reader.wait)
                    		scheduler.shutdown();
                    	
            			Second period = new Second();
            			allReadsCount.add(period, reader.twoDCount,"2D");
            			allReadsCount.add(period, reader.tempCount,"template");
            			allReadsCount.add(period, reader.compCount,"complement");

            
            			Demultiplexer myDmplx = reader.dmplx;
            			if(myDmplx!=null){
            				for(int i=0;i<myDmplx.readCount.length;i++){
            					demultiplexedStackReadsCount.add(period, myDmplx.readCount[i], myDmplx.barCodes.get(i).getName());
            					
            					demultiplexedBarReadsCount.setValue(myDmplx.readCount[i], "Count", myDmplx.barCodes.get(i).getName());
            				}
            				
            			}
            			
            			txtTFiles.setText(reader.getTotalFilesNumber()+"");	                
            			txtPFiles.setText(reader.getOKFilesNumber()+"");
            			txtFFiles.setText(reader.getSkippedFilesNumber()+"");
            
            			txt2DReads.setText(reader.twoDCount+"");
            			txtTempReads.setText(reader.tempCount+"");
            			txtCompReads.setText(reader.compCount+"");

            
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
            			currentIndex = reader.qualTemp.size();
            			if (currentIndex > lastIndexQualTemp){
            				int index = histoQualDataSet.getSeriesIndex("template");
            				for (int i = lastIndexQualTemp; i < currentIndex;i++)
            					histoQualDataSet.addSeries(index, reader.qualTemp.get(i));
            
            				lastIndexQualTemp = currentIndex;
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
                        
                    }
                }, 
                1, 
                1, 
                TimeUnit.SECONDS);     
	}

}
