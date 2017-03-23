package japsadev.seq.nanopore;
import japsa.util.DynamicHistogram;
import japsa.util.JapsaException;
import japsa.util.Logging;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;

import javax.swing.JOptionPane;

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
import javafx.stage.Stage;
public class NanoporeReaderWindowFX extends Application{
	
	static TimeTableXYDataset dataSet;
	static NanoporeReaderStream reader;
	public static void setData(TimeTableXYDataset data){
		dataSet = data;
	}
	public static void setReader(NanoporeReaderStream r){
		reader = r;
	}
	
    public static void main(String[] args) {
        Application.launch(args);
    }
    
    public void start(Stage primaryStage){  	 
	// Use a border pane as the root for scene
	    BorderPane border = new BorderPane();
	// Put start/stop button here  
	    HBox hbox = addHBox();
	    border.setTop(hbox);
	    
	// All the parameters setting to the left 
	    border.setLeft(addVBox());
	        
	// Add a stack to the HBox in the top region
        addStackPane(hbox);  
        
    // Here the main content    
        border.setCenter(addGridPane());
        
 
        Scene scene = new Scene(border);
        primaryStage.setScene(scene);
        primaryStage.setTitle("npreader");
	    primaryStage.show();
    	
    }
    
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
        Button buttonStart = new Button("Start", viewStart);
        buttonStart.setPrefSize(100, 20);
        buttonStart.setOnAction((event) -> {
        	
        });
        
        Image imageStop = new Image(getClass().getResourceAsStream("/stop.png"));
        ImageView viewStop = new ImageView(imageStop); 
        viewStop.setFitWidth(20);
        viewStop.setFitHeight(20);
        Button buttonStop = new Button("Stop", viewStop);
        buttonStop.setPrefSize(100, 20);
        buttonStop.setOnAction((event) -> {
        	
        });
        
        
        hbox.getChildren().addAll(buttonStart, buttonStop);
        
        return hbox;
    }
        
    private final int LeftPaneWidth=400;
    /*
     * Creates a VBox with a list of parameter settings
     */
    private VBox addVBox() {
        
        VBox vbox = new VBox();
        vbox.setPadding(new Insets(10)); // Set all sides to 10
        vbox.setSpacing(8);              // Gap between nodes
 
        Text title = new Text("Settings");
        title.setFont(Font.font("Arial", FontWeight.BOLD, 15));
        vbox.getChildren().add(title);
        final Separator sep1 = new Separator();
        sep1.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(1, sep1);
        
        vbox.getChildren().add(addInputPane());
        vbox.setSpacing(5);
        final Separator sep2 = new Separator();
        sep2.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(3, sep2);

        
        vbox.getChildren().add(addOutputPane());
        vbox.setSpacing(5);
        final Separator sep3 = new Separator();
        sep3.setMaxWidth(LeftPaneWidth);
        vbox.getChildren().add(5, sep3);
        
        vbox.getChildren().add(addOptionPane());
        
        vbox.getChildren().add(addBarcodePane());
        
        
        return vbox;
    }
    private GridPane addInputPane() {
		// TODO Auto-generated method stub
    	GridPane inputPane = createGridPane(LeftPaneWidth, 5);
    	
    	final Label inputLabel = new Label("Input:");
    	inputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	inputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(inputLabel, 0,0);
    	inputPane.getChildren().add(inputLabel);
    	
    	final CheckBox includeFailBox = new CheckBox("Include fail folder");
    	includeFailBox.setOnAction((event) -> {
    	    //reader.doFail = includeFailBox.isSelected();
    	    
    	});	
    	GridPane.setConstraints(includeFailBox, 2,0,3,1);
    	inputPane.getChildren().add(includeFailBox);
    	
    	final TextField textField = new TextField();
    	textField.setPromptText("Enter folder of basecalled reads...");
    	//textField.setPrefWidth(250);
    	GridPane.setConstraints(textField, 0,1,4,1);
    	inputPane.getChildren().add(textField);
    	

    	ImageButton inputBrowseButton = new ImageButton("/folder.png");
    	inputBrowseButton.setPrefSize(10, 10);
    	inputBrowseButton.setOnAction((event) -> {
        	
        });
    	GridPane.setConstraints(inputBrowseButton, 4,1);
    	GridPane.setHalignment(inputBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(inputBrowseButton);
    	//inputPane.setGridLinesVisible(true);

    	final CheckBox barcodeBox = new CheckBox("Demultiplexing for barcode analysis");
    	barcodeBox.setOnAction((event) -> {
    	    //reader.doFail = includeFailBox.isSelected();
    	    
    	});	
    	GridPane.setConstraints(barcodeBox, 0,3,5,1);
    	inputPane.getChildren().add(barcodeBox);
    	
    	final TextField barcodeInputField = new TextField();
    	barcodeInputField.setPromptText("Enter name of barcode sequences file...");
    	GridPane.setConstraints(barcodeInputField, 0,4,4,1);
    	inputPane.getChildren().add(barcodeInputField);
    	
    	ImageButton barcodeBrowseButton = new ImageButton("/folder.png");
    	barcodeBrowseButton.setPrefSize(10, 10);
    	barcodeBrowseButton.setOnAction((event) -> {
        	
        });
    	GridPane.setConstraints(barcodeBrowseButton, 4,4);
    	GridPane.setHalignment(barcodeBrowseButton,HPos.LEFT);
    	inputPane.getChildren().add(barcodeBrowseButton);
    	
		return inputPane;
	}
    private GridPane addOutputPane() {
		// TODO Auto-generated method stub
    	GridPane outputPane = createGridPane(LeftPaneWidth, 5);
    	
    	final Label outputLabel = new Label("Output:");
    	outputLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	outputLabel.setStyle("-fx-underline:true");
    	GridPane.setConstraints(outputLabel, 0,0);
    	outputPane.getChildren().add(outputLabel);
    	
        final ComboBox<String> outputDestComboBox = new ComboBox<String>();
        outputDestComboBox.getItems().addAll("to file", "to stdout");   
        outputDestComboBox.setValue("to file");
        GridPane.setConstraints(outputDestComboBox, 1, 0, 2, 1);
        outputPane.getChildren().add(outputDestComboBox);
        
        final ComboBox<String> outputFormatComboBox = new ComboBox<String>();
        outputFormatComboBox.getItems().addAll("fastq", "fasta");   
        outputFormatComboBox.setValue("fastq");
        GridPane.setConstraints(outputFormatComboBox, 3, 0);
        outputPane.getChildren().add(outputFormatComboBox);
    	
    	final TextField outputTextField = new TextField();
    	outputTextField.setPromptText("Enter name of output file...");
    	GridPane.setConstraints(outputTextField, 0,1,4,1);
    	outputPane.getChildren().add(outputTextField);
    	

    	ImageButton outputBrowseButton = new ImageButton("/folder.png");
    	outputBrowseButton.setPrefSize(10, 10);
    	outputBrowseButton.setOnAction((event) -> {
        	
        });
    	GridPane.setConstraints(outputBrowseButton, 4,1);
    	GridPane.setHalignment(outputBrowseButton, HPos.LEFT);
    	outputPane.getChildren().add(outputBrowseButton);


    	final CheckBox serversStreaming = new CheckBox("Streaming output to server(s)");
    	serversStreaming.setOnAction((event) -> {
    	    //reader.doFail = includeFailBox.isSelected();
    	    
    	});	
    	GridPane.setConstraints(serversStreaming, 0,4,3,1);
    	outputPane.getChildren().add(serversStreaming);
    	
    	final TextField outputStreamField = new TextField();
    	outputStreamField.setPromptText("address1:port1, address2:port2,...");
    	GridPane.setConstraints(outputStreamField, 0,5,4,1);
    	outputPane.getChildren().add(outputStreamField);
    	
    	//outputPane.setGridLinesVisible(true);
		return outputPane;
	}
    private GridPane addOptionPane() {
		// TODO Auto-generated method stub
    	GridPane optionPane = createGridPane(LeftPaneWidth, 5);
    	
    	final Label optLabel = new Label("Other options:");
    	optLabel.setFont(Font.font("Roman", FontWeight.BOLD, 12));
    	optLabel.setStyle("-fx-underline:true");
    	
    	GridPane.setConstraints(optLabel, 0,0,4,1);
    	optionPane.getChildren().add(optLabel);
    	
    	final CheckBox includeAllReadBox = new CheckBox("Including template and complement reads");
    	includeAllReadBox.setOnAction((event) -> {
    	    //reader.doFail = includeFailBox.isSelected();
    	    
    	});	
    	GridPane.setConstraints(includeAllReadBox, 0,2,4,1);
    	optionPane.getChildren().add(includeAllReadBox);

    	final CheckBox numGenerateBox = new CheckBox("Assign unique number to every read name");
    	numGenerateBox.setOnAction((event) -> {
    	    //reader.doFail = includeFailBox.isSelected();
    	    
    	});	
    	GridPane.setConstraints(numGenerateBox, 0,4,4,1);
    	optionPane.getChildren().add(numGenerateBox);
    	
    	final Label label2 = new Label("Filter out read shorter than ");
    	GridPane.setConstraints(label2, 0,6,3,1);
    	optionPane.getChildren().add(label2);
    	
    	final TextField minLenField = new TextField();
    	minLenField.setPromptText("1");
    	GridPane.setConstraints(minLenField, 3,6);
    	optionPane.getChildren().add(minLenField);
    	
    	final Label label3 = new Label("bp");
    	GridPane.setConstraints(label3, 4,6);
    	optionPane.getChildren().add(label3);
    	
		return optionPane;
	}
    private StackPane addBarcodePane() {
		// TODO Auto-generated method stub
    	StackPane barcodePane = new StackPane();
    	
    	
		return barcodePane;
	}
    
    private GridPane createGridPane(int width, int ncols){
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
	/*
     * Uses a stack pane to create a help icon and adds it to the right side of an HBox
     * 
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
        
        Text helpText = new Text("?");
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
 
        GridPane grid = new GridPane();
        grid.setHgap(10);
        grid.setVgap(10);
        grid.setPadding(new Insets(0, 10, 0, 10));
 
        // Category in column 2, row 1
        Text category = new Text("Sales:");
        category.setFont(Font.font("Arial", FontWeight.BOLD, 20));
        grid.add(category, 1, 0); 
        
        // Title in column 3, row 1
        Text chartTitle = new Text("Current Year");
        chartTitle.setFont(Font.font("Arial", FontWeight.BOLD, 20));
        grid.add(chartTitle, 2, 0);
        
        // Subtitle in columns 2-3, row 2
        Text chartSubtitle = new Text("Goods and Services");
        grid.add(chartSubtitle, 1, 1, 2, 1);
        
        // House icon in column 1, rows 1-2
        ImageView imageHouse = new ImageView(
                    new Image(getClass().getResourceAsStream("/house.png")));
        grid.add(imageHouse, 0, 0, 1, 2);
 
        // Left label in column 1 (bottom), row 3
        Text goodsPercent = new Text("Goods\n80%");
        GridPane.setValignment(goodsPercent, VPos.BOTTOM);
        grid.add(goodsPercent, 0, 2);
        
        // Chart in columns 2-3, row 3
        ImageView imageChart = new ImageView(
                    new Image(getClass().getResourceAsStream("/piechart.png")));
        grid.add(imageChart, 1, 2, 2, 1);
        
        // Right label in column 4 (top), row 3
        Text servicesPercent = new Text("Services\n20%");
        GridPane.setValignment(servicesPercent, VPos.TOP);
        grid.add(servicesPercent, 3, 2);
        
//            grid.setGridLinesVisible(true);
        return grid;
    }
    
    
    /******************************************************************************************
     * ** Here are variables and controls for all plots ***************************************
     ******************************************************************************************/
	DynamicHistogram histoLengthDataSet, histoQualDataSet;

	private boolean stillRun = true;

	private void interupt(JapsaException e){
		this.stillRun = false;
		reader.wait = false;
		JOptionPane.showMessageDialog(null, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
	}

	private void run() {		
		int lastIndexLengths = 0;//, lastIndexLengths2D = 0, lastIndexLengthsComp = 0, lastIndexLengthsTemp = 0;
		int lastIndexQual2D = 0, lastIndexQualComp = 0, lastIndexQualTemp = 0;

		while(stillRun) {				                
			//synchronized(reader) {//avoid concurrent update					
			Second period = new Second();
			dataSet.add(period, reader.twoDCount,"2D");
			dataSet.add(period, reader.compCount,"complement");
			dataSet.add(period, reader.tempCount,"template");

//			txtTFiles.setText(reader.fileNumber+"");	                
//			txtPFiles.setText(reader.passNumber+"");
//			txtFFiles.setText(reader.failNumber+"");
//
//			txt2DReads.setText(reader.twoDCount+"");
//			txtCompReads.setText(reader.compCount+"");
//			txtTempReads.setText(reader.tempCount+"");

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
