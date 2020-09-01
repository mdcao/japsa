package japsa.bio.phylo;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import pal.tree.Node;
import pal.tree.NodeUtils;

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
			//Node root = combined.roots.get(1);
		//	Node n = NodeUtils.postorderSuccessor(root);
			//Node a = combined.taxa2Node.get("1812935");
			//Node b = combined.taxa2Node.get("2027919");
		//	System.err.println(a.getParent().getIdentifier());
		//	System.err.println(b.getParent().getIdentifier());
/*
			while(n!=root){
				System.err.println(n.getIdentifier().getAttribute("taxon"));
				n = NodeUtils.postorderSuccessor(n);
			}
*/			
		//NodeUtils.preOrderIterator(combined.roots.get(1));
		
			//combined.makeTrees();
			StringBuffer header = new StringBuffer();
			for(int i=0; i<f.length; i++){
				header.append(f[i].getName()+"\t");
			}
			header.append("name\ttaxon\tlevel\tcolor\ttaxon1\ttaxon2\ttaxon3");
			combined.makeTrees();

			combined.print(new File(outfile),"", header.toString());
			
			//System.err.println(combined);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}