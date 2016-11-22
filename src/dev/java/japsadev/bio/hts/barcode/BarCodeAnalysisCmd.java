package japsadev.bio.hts.barcode;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

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

		addStdHelp();
	}
	public static void main(String[] args){
		BarCode bc;
		try {
			bc = new BarCode("/home/hoangnguyen/workspace/poreFUME/inputData/pb_39.fasta");
			bc.clustering("/home/hoangnguyen/workspace/poreFUME/inputData/n.fasta.protein.homolog.fasta");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
	}
}
