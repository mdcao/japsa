package japsadev.tools;

import java.io.IOException;

import jaligner.matrix.MatrixLoaderException;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
//import japsadev.bio.hts.barcode.BarCode;
import japsadev.bio.hts.barcode.BarCodeAnalysis;
import japsadev.bio.hts.barcode.BarCodeAnalysisVerbose;

@Deployable(
		scriptName = "jsa.dev.barcode", 
		scriptDesc = "Clustering nanopore sequences based on barcode"
		)
public class BarCodeAnalysisCmd extends CommandLine{
	public BarCodeAnalysisCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("bcFile", null, "Barcode file",true);		
		addString("seqFile", null, "Nanopore sequences file",true);
		addString("scriptRun", null, "Invoke command script to run npScarf",true);
		addBoolean("verbose", false, "Verbose mode to show alignments (blossom62)");
		addStdHelp();
	}
	public static void main(String[] args) throws IOException, InterruptedException, MatrixLoaderException{
		CommandLine cmdLine = new BarCodeAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		String bcFile = cmdLine.getStringVal("bcFile");
		String script = cmdLine.getStringVal("scriptRun");
		String seqFile = cmdLine.getStringVal("seqFile");
		Boolean v = cmdLine.getBooleanVal("verbose");

		if(v){
			BarCodeAnalysisVerbose bc = new BarCodeAnalysisVerbose(bcFile,script);
			bc.clustering(seqFile);
		}else{
			BarCodeAnalysis bc = new BarCodeAnalysis(bcFile,script);
			bc.clustering(seqFile);
		}
		
	}
}
