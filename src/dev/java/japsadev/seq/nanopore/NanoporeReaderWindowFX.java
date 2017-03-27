package japsadev.seq.nanopore;
import japsa.util.DynamicHistogram;
import japsa.util.JapsaException;
import japsa.util.Logging;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

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

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingNode;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.HPos;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.geometry.Side;
import javafx.geometry.VPos;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Hyperlink;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.RowConstraints;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.paint.CycleMethod;
import javafx.scene.paint.LinearGradient;
import javafx.scene.paint.Stop;
import javafx.scene.shape.Rectangle;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.FileChooser.ExtensionFilter;
import javafx.stage.Stage;
public class NanoporeReaderWindowFX extends Application{
	
	static TimeTableXYDataset dataSet = new TimeTableXYDataset();
	static NanoporeReaderStream reader = new NanoporeReaderStream();
	public static void setData(TimeTableXYDataset data){
		dataSet = data;
	}
	public static void setReader(NanoporeReaderStream r){
		reader = r;
	}
	
    public static void main(String[] args) {
    	//set reader manually for testing purpose. comment it out for the final version
    	reader.folder = "/home/hoangnguyen/workspace/data/rawNanopore/downloads";
    	//reader.updateDemultiplexFile("/home/hoangnguyen/workspace/data/barcode/barcode-npreader/nnp_barcode.fasta");
		reader.realtime = true;
		System.setProperty("java.awt.headless", "false");
		reader.stats = true;//GUI implies stats
		reader.ready = false;
		
		new Thread(new Runnable(){

			@Override
			public void run() {
				while (!reader.ready){
					Logging.info("NOT READY");
					try {
						Thread.sleep(1000);
					} catch (InterruptedException e) {					
						e.printStackTrace();
					}			
				}
				// TODO Auto-generated method stub
				try{
					reader.readFast5();
				}catch (JapsaException e){
					System.err.println(e.getMessage());
					e.getStackTrace();
					interupt(e);
				}catch (Exception e){
					e.printStackTrace();;
				}
			}
			
		}).start();
		
        Application.launch(args);
    }
    
    public void start(Stage primaryStage){  	 
    	
    	//Start with a BorderPane as root
	    BorderPane border = new BorderPane();
	    // Put start/stop button here  
	    HBox hbox = addHBox();
	    border.setTop(hbox);
	    
	    // All the parameters setting to the left 
	    border.setLeft(addVBox(primaryStage));
	        
	    // Add a stack to the HBox in the top region
        addStackPane(hbox);  
        
        // Here the main content    
        border.setCenter(addGridPane());
        
 
        Scene scene = new Scene(border);
        primaryStage.setScene(scene);
        primaryStage.setTitle("npreader");
	    primaryStage.show();
    	

	     Task<Void> task = new Task<Void>() {
	         @Override protected Void call() throws Exception {
                 if (isCancelled()) return null;

                 Platform.runLater(new Runnable() {
                     @Override public void run() {
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
                 });
	             
	             return null;
	         }
	     };
        Thread th = new Thread(task);
        th.setDaemon(true);
        th.start();
    }
    /*
     * Components from left pane
     */
    private Button 	buttonStart, buttonStop,
    				inputBrowseButton, barcodeBrowseButton, outputBrowseButton;
    private TextField inputTF, barcodeTF, outputTF, streamTF, minLenTF;
    private CheckBox failCB, barcodeCB, serversCB, allReadsOptCB, addNumberOptCB;
    private ComboBox<String> outputToCombo, outputFormatCombo;
    /*
     * Creates an HBox with two buttons for the top region
     */
        
    private HBox addHBox() {
 
        HBox hbox = new HBox();
        hbox.setPadding(new Insets(15, 12, 15, 12));
        hbox.setSpacing(10);   // Gap between nodes
        hbox.setStyle("-fx-background-color: #336699;");
 

        Image imageStart = new Image(getClass().getResourceAsStream("/start.png"));
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
				}
				reader.output = _foutput;					
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
			}
			//TODO: set barcode tab visible here...

			String msg = reader.prepareIO();
			if (msg !=null){
				FxDialogs.showWarning("Warning", msg);
				return;
			}

			//Start running
			
			inputBrowseButton.setDisable(true);
			barcodeBrowseButton.setDisable(true);
			outputBrowseButton.setDisable(true);
			
			inputTF.setDisable(true);
			barcodeTF.setDisable(true);
			outputTF.setDisable(true);
			streamTF.setDisable(true);
			minLenTF.setDisable(true);
			
			failCB.setDisable(true);
			barcodeCB.setDisable(true);
			serversCB.setDisable(true);
			allReadsOptCB.setDisable(true);
			addNumberOptCB.setDisable(true);
			
			outputToCombo.setDisable(true);
			outputFormatCombo.setDisable(true);

			buttonStart.setDisable(true);;
			buttonStop.setDisable(false);

			reader.ready = true;
		});
        
        Image imageStop = new Image(getClass().getResourceAsStream("/stop.png"));
        ImageView viewStop = new ImageView(imageStop); 
        viewStop.setFitWidth(20);
        viewStop.setFitHeight(20);
        buttonStop = new Button("Stop", viewStop);
        buttonStop.setPrefSize(100, 20);
        buttonStop.setDisable(true);
        buttonStop.setOnAction((event) -> {
        	reader.wait = false;

//			while (!reader.done){
//				try {
//					Thread.sleep(100);
//				} catch (InterruptedException e) {					
//					e.printStackTrace();
//				}
//			}

		stillRun = false;
		//FxDialogs.showInformation("Process being terminated!", "Done");
		String confirm = FxDialogs.showConfirm( "Remember to save all plots before quit", "Terminate the GUI", "Yes", "No");
		if(confirm.equals("Yes")){
			Platform.exit();
		    System.exit(0);
		}
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
    	
    	GridPane.setConstraints(failCB, 2,0,3,1);
    	inputPane.getChildren().add(failCB);
    	
    	inputTF = new TextField(reader.folder == null?"":reader.folder);
    	inputTF.setPromptText("Enter folder of basecalled reads...");
    	//textField.setPrefWidth(250);
    	GridPane.setConstraints(inputTF, 0,1,4,1);
    	inputPane.getChildren().add(inputTF);
    	

    	inputBrowseButton = new ImageButton("/folder.png");
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
    	barcodeCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	barcodeTF.setDisable(!new_val);
                	barcodeBrowseButton.setDisable(!new_val);
                });	
    	GridPane.setConstraints(barcodeCB, 0,3,5,1);
    	inputPane.getChildren().add(barcodeCB);
    	
    	barcodeTF = new TextField(reader.getBCFileName() == null?"":reader.getBCFileName());
    	barcodeTF.setPromptText("Enter name of barcode sequences file...");
    	barcodeTF.setDisable(!barcodeCB.isSelected());
    	GridPane.setConstraints(barcodeTF, 0,4,4,1);
    	inputPane.getChildren().add(barcodeTF);
    	
    	barcodeBrowseButton = new ImageButton("/folder.png");
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
    	GridPane.setConstraints(barcodeBrowseButton, 4,4);
    	GridPane.setHalignment(barcodeBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(barcodeBrowseButton);
    	
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
        outputToCombo.valueProperty().addListener((obs_val, old_val, new_val) -> {
        	if(new_val.trim().equals("to file")){
        		outputTF.setText("");
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
    	
    	outputTF = new TextField(reader.output == null?"":reader.output);
    	outputTF.setPromptText("Enter name for output file...");
    	GridPane.setConstraints(outputTF, 0,1,4,1);
    	outputPane.getChildren().add(outputTF);
    	

    	outputBrowseButton = new ImageButton("/folder.png");
    	outputBrowseButton.setPrefSize(10, 10);
    	outputBrowseButton.setOnAction((event) -> {
    		FileChooser fileChooser = new FileChooser();
    		fileChooser.setTitle("Save output to file");
    		fileChooser.setInitialFileName("output."+reader.format);
    		File savedFile = fileChooser.showSaveDialog(stage);
    		if(savedFile != null){
    			reader.output = savedFile.getName();
    			outputTF.setText(reader.output);
    		}
        });
    	GridPane.setConstraints(outputBrowseButton, 4,1);
    	GridPane.setHalignment(outputBrowseButton, HPos.LEFT);
    	outputPane.getChildren().add(outputBrowseButton);

    	//init
        if(reader.output.equals("-")){
        	outputTF.setText("-");
        	outputToCombo.setValue("to stdout");
    		outputTF.setDisable(true);
    		outputBrowseButton.setDisable(true);
        }else{
        	outputTF.setText(reader.output);
        	outputToCombo.setValue("to file");
    		outputTF.setDisable(false);
    		outputBrowseButton.setDisable(false);
        }
        
    	
    	serversCB = new CheckBox("Streaming output to server(s)");
    	serversCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	streamTF.setDisable(!new_val);
                });	
    	GridPane.setConstraints(serversCB, 0,4,3,1);
    	outputPane.getChildren().add(serversCB);
    	
    	streamTF = new TextField();
    	streamTF.setPromptText("address1:port1, address2:port2,...");
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
    	
    	allReadsOptCB = new CheckBox("Including template and complement reads");
    	allReadsOptCB.setSelected(true);
    	allReadsOptCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	reader.doLow=new_val;
                });	
    	GridPane.setConstraints(allReadsOptCB, 0,2,4,1);
    	optionPane.getChildren().add(allReadsOptCB);

    	addNumberOptCB = new CheckBox("Assign unique number to every read name");
    	addNumberOptCB.setSelected(reader.number);
    	addNumberOptCB.selectedProperty().addListener(
                (obs_val,old_val,new_val) -> {
                	reader.number=new_val;
                });	
    	GridPane.setConstraints(addNumberOptCB, 0,4,4,1);
    	optionPane.getChildren().add(addNumberOptCB);
    	
    	final Label label2 = new Label("Filter out read shorter than ");
    	GridPane.setConstraints(label2, 0,6,3,1);
    	optionPane.getChildren().add(label2);
    	
    	minLenTF = new TextField(Integer.toString(reader.minLength));
    	minLenTF.setPromptText("Enter minimum length to consider...");
    	GridPane.setConstraints(minLenTF, 3,6);
    	optionPane.getChildren().add(minLenTF);
    	
    	final Label label3 = new Label("bp");
    	GridPane.setConstraints(label3, 4,6);
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
     * Uses a stack pane to create a help icon and adds it to the right side of an HBox
     * TODO: handle it!!!
     * @param hb HBox to add the stack to
     */
    private void addStackPane(HBox hb) {
 
        StackPane stack = new StackPane();
        Rectangle helpIcon = new Rectangle(30.0, 25.0);
        helpIcon.setFill(new LinearGradient(0,0,0,1, true, CycleMethod.NO_CYCLE,
            new Stop[]{
            new Stop(0,Color.web("#4977A3")),
            new Stop(0.5, Color.web("#B0C6DA")),
            new Stop(1,Color.web("#9CB6CF")),}));
        helpIcon.setStroke(Color.web("#D0E6FA"));
        helpIcon.setArcHeight(3.5);
        helpIcon.setArcWidth(3.5);
        
        final Text helpText = new Text("?");
        helpText.setFont(Font.font("Verdana", FontWeight.BOLD, 18));
        helpText.setFill(Color.WHITE);
        helpText.setStroke(Color.web("#7080A0")); 
        
        stack.getChildren().addAll(helpIcon, helpText);
        stack.setAlignment(Pos.CENTER_RIGHT);
        // Add offset to right for question mark to compensate for RIGHT 
        // alignment of all nodes
        StackPane.setMargin(helpText, new Insets(0, 10, 0, 0));
        
        hb.getChildren().add(stack);
        HBox.setHgrow(stack, Priority.ALWAYS);
                
    }
     
    /*
     * Creates a grid for the center region with four columns and three rows
     */
    private GridPane addGridPane() {
 
        GridPane mainGrid = createAutoresizeGridPane(2,2);
        mainGrid.setStyle("-fx-background-color: #C0C0C0;");
        /*
         * Read count chart
         */
		final JFreeChart chart = ChartFactory.createStackedXYAreaChart(
				"Read count",      // chart title
				"Time",             // domain axis label
				"Read number",                   // range axis label
				dataSet   
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
		
		SwingNode qualitySwingNode = new SwingNode();
		qualitySwingNode.setContent(hisQualPanel);
		GridPane.setConstraints(qualitySwingNode, 1,1);
		mainGrid.getChildren().add(qualitySwingNode);
        /*
         * Statistics field
         */
//        GridPane countPane = new GridPane();
//        countPane.setPadding(new Insets(10, 10, 10, 10));
//        countPane.setVgap(5);
//        countPane.setHgap(5);
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

		final Label lblpFiles = new Label("Pass files");
		GridPane.setConstraints(lblpFiles, 0, 1);
		countPane.getChildren().add(lblpFiles);

		//txtPFiles = new TextField("0");
		txtPFiles.setEditable(false);
		txtPFiles.setPrefWidth(100);
		GridPane.setConstraints(txtPFiles, 1, 1);
		countPane.getChildren().add(txtPFiles);


		final Label lblFFiles = new Label("Fail files");
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
    
    
    /******************************************************************************************
     * ** Here are variables and controls for all plots ***************************************
     ******************************************************************************************/
	final TextField txtCompReads= new TextField("0"), 
					txtTempReads= new TextField("0"), 
					txt2DReads= new TextField("0");
	final TextField txtPFiles= new TextField("0"), 
					txtFFiles= new TextField("0"), 
					txtTFiles= new TextField("0");
	final DynamicHistogram 	histoLengthDataSet = new DynamicHistogram(), 
						histoQualDataSet = new DynamicHistogram();

	private static boolean stillRun = true;

	private static void interupt(JapsaException e){
		stillRun = false;
		reader.wait = false;
		FxDialogs.showError("Unexpected errors happened!", e.getMessage());
	}
	
//TODO: fix this by concurrent/binding
//	private void run() {		
//		int lastIndexLengths = 0;//, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
//		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;
//
//		while(stillRun) {				                
//			//synchronized(reader) {//avoid concurrent update					
//			Second period = new Second();
//			dataSet.add(period, reader.twoDCount,"2D");
//			dataSet.add(period, reader.compCount,"complement");
//			dataSet.add(period, reader.tempCount,"template");
//
//			txtTFiles.setText(reader.fileNumber+"");	                
//			txtPFiles.setText(reader.passNumber+"");
//			txtFFiles.setText(reader.failNumber+"");
//
//			txt2DReads.setText(reader.twoDCount+"");
//			txtCompReads.setText(reader.compCount+"");
//			txtTempReads.setText(reader.tempCount+"");
//
//			int currentIndex = reader.lengths.size();
//
//			if (currentIndex > lastIndexLengths){			
//				int index = histoLengthDataSet.getSeriesIndex("Read Length");
//				for (int i = lastIndexLengths; i < currentIndex;i++)
//					histoLengthDataSet.addSeries(index, reader.lengths.get(i));
//
//				lastIndexLengths = currentIndex;
//
//				histoLengthDataSet.notifyChanged();
//			}
//
//			currentIndex = reader.qual2D.size();
//			if (currentIndex > lastIndexQual2D){
//				int index = histoQualDataSet.getSeriesIndex("2D");
//				for (int i = lastIndexQual2D; i < currentIndex;i++)
//					histoQualDataSet.addSeries(index, reader.qual2D.get(i));
//
//				lastIndexQual2D = currentIndex;
//				histoQualDataSet.notifyChanged();
//			}
//
//			currentIndex = reader.qualComp.size();
//			if (currentIndex > lastIndexQualComp){
//				int index = histoQualDataSet.getSeriesIndex("complement");
//				for (int i = lastIndexQualComp; i < currentIndex;i++)
//					histoQualDataSet.addSeries(index, reader.qualComp.get(i));
//
//				lastIndexQualComp = currentIndex;
//				histoQualDataSet.notifyChanged();
//			}
//
//			currentIndex = reader.qualTemp.size();
//			if (currentIndex > lastIndexQualTemp){
//				int index = histoQualDataSet.getSeriesIndex("template");
//				for (int i = lastIndexQualTemp; i < currentIndex;i++)
//					histoQualDataSet.addSeries(index, reader.qualTemp.get(i));
//
//				lastIndexQualTemp = currentIndex;
//				histoQualDataSet.notifyChanged();
//			}
//
//			try {
//				Thread.sleep(1000);
//			} catch (InterruptedException ex) {
//				Logging.error(ex.getMessage());
//			}
//		}
//	}
//    @Override
//    public void start(Stage primaryStage) {
//        primaryStage.setTitle("Tabs");
//        Group root = new Group();
//        Scene scene = new Scene(root, 400, 250, Color.WHITE);
//
//        TabPane tabPane = new TabPane();
//
//        BorderPane borderPane = new BorderPane();
//        for (int i = 0; i < 5; i++) {
//            Tab tab = new Tab();
//            tab.setText("Tab" + i);
//            HBox hbox = new HBox();
//            hbox.getChildren().add(new Label("Tab" + i));
//            hbox.setAlignment(Pos.CENTER);
//            tab.setContent(hbox);
//            tabPane.getTabs().add(tab);
//        }
//        // bind to take available space
//        borderPane.prefHeightProperty().bind(scene.heightProperty());
//        borderPane.prefWidthProperty().bind(scene.widthProperty());
//        
//        borderPane.setCenter(tabPane);
//        root.getChildren().add(borderPane);
//        primaryStage.setScene(scene);
//        primaryStage.show();
//    }

}
