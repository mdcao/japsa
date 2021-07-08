package japsa.bio.phylo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
		addString("output0", "cumul.css", "output file for sum of all depths below", false);
		addString("output1", "sep.css", "output file for node specific", false);
		addBoolean("trim",false,"whether to trim the species name");
		addString("pattern", null, "input pattern", false);//".outreport"
		addString("extra", null, "extra_input", false);
		addString("dirs", ".", "input directories", false);
		addString("todo",null, "path to file with list of input files",false);
		addDouble("thresh", 0, "threshold", false);

		addStdHelp();
	}
	
	static boolean removeLeafNodes = false;
	
	public static void main(String[] args){
		CommandLine cmdLine = new MergeKrakenCmd();		
		args = cmdLine.stdParseLine(args);	
		String outfile = cmdLine.getStringVal("output");
		String outfile1 = cmdLine.getStringVal("output1");
		String regex= cmdLine.getStringVal("pattern");
		String todo = cmdLine.getStringVal("todo");
		NCBITree.trim = cmdLine.getBooleanVal("trim");
		NCBITree.thresh = cmdLine.getDoubleVal("thresh");
	
		File outfileF = new File(outfile);
		File outfileF1 = new File(outfile1);
		outfileF.delete();
		outfileF1.delete();
		
		String ext = cmdLine.getStringVal("extra");
		
		String[] extra = ext!=null ? ext.split(":") : new String[0];
		List<String>f = new ArrayList<String>();
		try {
		if(todo!=null){
			try{
			BufferedReader br = new BufferedReader(new FileReader(todo));
			String st = "";
			while((st = br.readLine())!=null){
				if(st.startsWith("#")) continue;
				File f1 = new File(st);
				if(!f1.exists()) throw new RuntimeException(st+ " does not exist");
				f.add(st);
			}
			}catch(IOException exc){
				exc.printStackTrace();
			}
		}else{
		String[] dirs = cmdLine.getStringVal("dirs").split(":");
		


		
			FilenameFilter filter = new FilenameFilter(){

				@Override
				public boolean accept(File fir, String pathname) {
//					return Pattern.matches(regex, pathname.getName());
					return pathname.endsWith(regex);
				}

				
			};
			
			for(int i=0; i<dirs.length; i++){
				f.addAll(Arrays.asList((new File(dirs[i])).list(filter)));
			}
		}
		//Collections.sort(f);
			NCBITree combined = new KrakenTree(new File(f.get(0)),regex);
			MergeKrakenCmd.checkDupl(combined.roots.get(0));
			combined.modAll(0, f.size()+ extra.length);
			String sl = combined.slug("Torque teno virus",false);
			int sze = combined.roots.get(0).getChildCount();
			Node virus = combined.slugToNode.get(sl);
			for(int i=1; i<f.size(); i++){
				KrakenTree kt = new KrakenTree(new File(f.get(i)),regex);
				MergeKrakenCmd.checkDupl(kt.roots.get(0));

				if(!kt.roots.get(0).getIdentifier().toString().equals("root")){
					System.err.println("excluding "+f.get(i));
					continue;
				}
					kt.modAll(i, f.size()+ extra.length);
					try{
					combined.merge(kt, i);
					checkDupl(combined.roots.get(0));
					int sze1 = combined.roots.get(0).getChildCount();
					System.err.println(sze1);
					}catch (Exception exc){
						
						System.err.println("problem with "+f.get(i));
						exc.printStackTrace();
						
					}
				//	Node contains = kt.slugToNode.get(sl);
				//	boolean contains1 = combined.slugToNode.containsKey(sl);
				//	System.err.println(sl+"  "+contains!=null+ " "+contains1);
			}
			Node virus1 = combined.slugToNode.get(sl);

			int sze1 = combined.roots.get(0).getChildCount();
		
			
			//System.err.println(sl+" "+combined.slugToNode.containsKey(sl));
		
			int len = f.size();
			
			for(int i=0; i<extra.length; i++){
				BufferedReader br = new BufferedReader(new FileReader(extra[i]));
				String str = "";
				while(( str = br.readLine())!=null){
					String[] st = str.split(",");
					int len1 = st.length-1;
					Double  val = Double.parseDouble(st[len1]);
					int val1 = (int) Math.round(val*1e6);
					Node n = combined.getSlug(st[len1-1]);
					if(n==null){
						System.err.println("did not find "+st[len1-1]);
					}else{
						System.err.println("found "+n.getIdentifier());

					}
					//System.err.println("found "+n.getIdentifier());
					
				
					int[] cnts1 = ((int[]) n.getIdentifier().getAttribute(NCBITree.count_tag1));
					cnts1[len+i] = val1;
					
					while(n!=null){//add counts up the tree for cumulative
						int[] cnts_ = ((int[]) n.getIdentifier().getAttribute(NCBITree.count_tag));
						cnts_[len+i] = cnts_[len+i] + val1;
						n = n.getParent();
						
					}
				}
				
			}
		
			//combined.trim(NCBITree.thresh);
			//combined.removeDupl();
			//combined.split(); combined.split();
			
			//combined.removeSingleNodes(); //Arrays.asList(new Integer[] {2}));
		//	combined.roots;
			System.err.println(combined.roots.size());

			StringBuffer header = new StringBuffer();
			File currDir = new File(".");
			header.append("name\tcolor\ttaxon\theight\tlevel\tlevel1\tcssvals\tparents\ttaxon_parents");
			PrintWriter designF = new PrintWriter(new FileWriter(new File("design.csv")));
			designF.println("Name,Grp1,Grp2");
		
			for(int i=0; i<f.size(); i++){
				String nmei = f.get(i).replaceAll("./", "");
				header.append("\t");
				
			//	String	nmei = fi.getName();
				String grp = "all";
				String[] split = nmei.split("/"); 
			//	if(fi.getParentFile()!=null){
					//String prefix = fi.getParentFile().equals(currDir) ? "": fi.getParentFile().getName()+"/";
				grp = split[0];
				for(int ik=1; ik<split.length-1; ik++){
					grp  =grp+"/"+split[ik];
				}
					//nmei = nmei+fi.getName();
				//}
				nmei = nmei.replace("results.krkn", "");
				designF.println(nmei+","+(i+1)+","+grp);
				header.append(nmei);
			}
			designF.close();
			for(int i=0; i<extra.length; i++){
				File extraf = new File(extra[i]);
				header.append("\t");header.append(extraf.getName());
			}
			
			if(removeLeafNodes) combined.removeLeafNodes(NCBITree.count_tag1, true); /// this should really be true I think, but reflects a problem with generating the original kraken files
			
			combined.makeTrees(false);
		//	Node virus1_ = combined.slugToNode.get(sl);
			System.err.println("h");
			for(int i=0; i<combined.tree.length; i++){
				//CSSProcessCommand.color(combined.tree);
			//	CSSProcessCommand.colorEachLevel(combined.tree);
				CSSProcessCommand.colorRecursive(combined.tree, true);
			}
			
			combined.print(outfileF,"", header.toString(), new String[] {NCBITree.count_tag}, new String[] {"%5.3g"}, false);
			combined.print(outfileF1,"", header.toString(),new String[] {NCBITree.count_tag1}, new String[] {"%5.3g"}, false);

			//System.err.println(combined);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	static void checkDupl(Node node) {
		Set<String> s = new HashSet<String>();
		
		
		for(int i=0; i<node.getChildCount(); i++){
			String str1 = node.getChild(i).getIdentifier().getName();
			if(s.contains(str1)) {
				throw new RuntimeException("!!"+str1);
			}else{
				s.add(str1);
			}
			
		}
		
	}
}