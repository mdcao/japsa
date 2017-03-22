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
import javafx.geometry.Pos;
import javafx.geometry.Side;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.Label;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.paint.Color;
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
    	
    	
    }

    
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
