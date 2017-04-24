package japsadev.obsolete.np;

import japsa.util.CommandLine;
import japsa.util.JapsaException;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "jsa.np.meth",
scriptDesc = "Detecting methylated bases using Oxford Nanopore sequencing signal")
public class BaseMethylationCmd extends CommandLine{
	public BaseMethylationCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("samFile", null, "Name of the .sam file when nanopore reads align with the reference");
		addString("refFile", null, "Name of the reference fasta file");
		addString("gffFile", null, "Name of the GFF file specifying methylation position of ref");
		addString("nanopore", null, "Name of folder containing nanopore raw reads in hdf5 format");
		addStdHelp();	
	}
	/**
	 * @param args
	 * @throws Exception 
	 * @throws JapsaException 
	 * @throws OutOfMemoryError 
	 */
	public static void main(String[] args) throws OutOfMemoryError, JapsaException, Exception {
		// TODO Auto-generated method stub
		CommandLine cmdLine = new BaseMethylationCmd();		
		args = cmdLine.stdParseLine(args);
		
		String refFile = cmdLine.getStringVal("refFile");
		String samFile = cmdLine.getStringVal("samFile");
		String gffFile = cmdLine.getStringVal("gffFile");
		String nanopore = cmdLine.getStringVal("nanopore");
		
		KmerMap map = new KmerMap(gffFile, "m6A");
//		map.scanUnmethylate(refFile);
//		map.scanAlignment(samFile);
//		map.scanHDF5(nanopore);
//		map.print();
		KmerMap.statsHDF5(nanopore);
	}

}
