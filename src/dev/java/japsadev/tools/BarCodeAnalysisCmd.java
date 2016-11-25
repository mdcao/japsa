package japsadev.tools;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
//import japsadev.bio.hts.barcode.BarCode;
import japsadev.bio.hts.barcode.BarCodeAnalysis;

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

		addStdHelp();
	}
	public static void main(String[] args){
		CommandLine cmdLine = new BarCodeAnalysisCmd ();
		args = cmdLine.stdParseLine(args);

		String bcFile = cmdLine.getStringVal("bcFile");
		String script = cmdLine.getStringVal("scriptRun");
		String seqFile = cmdLine.getStringVal("seqFile");



		BarCodeAnalysis bc;
		try {
			bc = new BarCodeAnalysis(bcFile,script);
			bc.clustering(seqFile);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 

	}
}
