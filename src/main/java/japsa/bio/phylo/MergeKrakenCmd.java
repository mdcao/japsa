package japsa.bio.phylo;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;

import japsa.tools.bio.hts.RepeatDetectionCmd;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

/**
 * @author lachlancoin
 *
 */
@Deployable(
		scriptName = "jsa.hts.mergeKraken",
		scriptDesc = "Merge kraken results")
public class MergeKrakenCmd extends CommandLine{
	public MergeKrakenCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("output", "combined.txt", "output file", false);
		addString("input", ".outreport", "input pattern", false);

		addStdHelp();
	}
	
	public static void main(String[] args){
		CommandLine cmdLine = new MergeKrakenCmd();		
		args = cmdLine.stdParseLine(args);	
		String outfile = cmdLine.getStringVal("output");
		String outreport = cmdLine.getStringVal("input");
		try {
			FileFilter filter = new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return pathname.getName().endsWith(outreport);
				}
				
			};
			File[] f = (new File(".")).listFiles(filter);
			StringBuffer header = new StringBuffer();
			for(int i=0; i<f.length; i++){
				header.append(f[i].getName()+"\t");
			}
			header.append("name\ttaxon");
			KrakenTree[] kt = new KrakenTree[f.length];
			for(int i=0; i<kt.length; i++){
				kt[i] = new KrakenTree(f[i]);
				kt[i].modAll(i, kt.length);
			//	kt[i].print(new File(i+".combined.txt"),"",
				//		header.toString()
					//	);

			}
			
			NCBITree combined = kt[0];
			for(int j=1; j<kt.length; j++){
				combined.merge(kt[j], j);
			}
			combined.print(new File(outfile),"", header.toString());
		///	combined.makeTrees();
			//System.err.println(combined);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}