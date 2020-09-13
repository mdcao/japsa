package japsa.bio.phylo;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
		addString("output0", "cumul.css", "output file for sum of all depths below", false);
		addString("output1", "sep.css", "output file for node specific", false);
		addBoolean("trim",false,"whether to trim the species name");
		addString("pattern", ".outreport", "input pattern", false);
		addString("dirs", ".", "input directories", false);
		addDouble("thresh", 0, "threshold", false);

		addStdHelp();
	}
	
	public static void main(String[] args){
		CommandLine cmdLine = new MergeKrakenCmd();		
		args = cmdLine.stdParseLine(args);	
		String outfile = cmdLine.getStringVal("output");
		String outfile1 = cmdLine.getStringVal("output1");
		String regex= cmdLine.getStringVal("pattern");
		String[] dirs = cmdLine.getStringVal("dirs").split(":");
NCBITree.trim = cmdLine.getBooleanVal("trim");
NCBITree.thresh = cmdLine.getDoubleVal("thresh");

		try {
			FileFilter filter = new FileFilter(){

				@Override
				public boolean accept(File pathname) {
//					return Pattern.matches(regex, pathname.getName());
					return pathname.getName().endsWith(regex);
				}
				
			};
			List<File>f = new ArrayList<File>();
			for(int i=0; i<dirs.length; i++){
				f.addAll(Arrays.asList((new File(dirs[i])).listFiles(filter)));
			}
		
			KrakenTree[] kt = new KrakenTree[f.size()];
			for(int i=0; i<kt.length; i++){
				kt[i] = new KrakenTree(f.get(i));
				kt[i].modAll(i, kt.length);
			//	kt[i].print(new File(i+".combined.txt"),"",
				//		header.toString()
					//	);

			}
			
			NCBITree combined = kt[0];
			for(int j=1; j<kt.length; j++){
				combined.merge(kt[j], j);
			}
			combined.trim(NCBITree.thresh);
			combined.removeDupl();
			combined.split(); combined.split();
			
			combined.removeSingleNodes(Arrays.asList(new Integer[] {2}));
		//	combined.roots;
			System.err.println(combined.roots.size());

			StringBuffer header = new StringBuffer();
			
			header.append("name\tcolor\ttaxon\theight\tlevel\tcssvals\tparents\ttaxon_parents");
			for(int i=0; i<f.size(); i++){
				File fi = f.get(i);
				header.append("\t");
				header.append(fi.getParentFile().getName()+"_"+fi.getName());
			}
			combined.makeTrees();
			for(int i=0; i<combined.tree.length; i++){
				//CSSProcessCommand.color(combined.tree);
			//	CSSProcessCommand.colorEachLevel(combined.tree);
				CSSProcessCommand.colorRecursive(combined.tree, true);
			}
			combined.print(new File(outfile),"", header.toString(), NCBITree.count_tag);
			combined.print(new File(outfile1),"", header.toString(), NCBITree.count_tag1);

			//System.err.println(combined);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}