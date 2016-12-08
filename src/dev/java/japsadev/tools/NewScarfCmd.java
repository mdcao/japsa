package japsadev.tools;
import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsadev.bio.hts.newscarf.HybridAssembler;

@Deployable(
		scriptName = "jsa.dev.newScarf", 
		scriptDesc = "New npscarf"
		)
public class NewScarfCmd extends CommandLine{
	public NewScarfCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addString("fastg", null, "Assembly graph fastg file",true);		
		addString("sam", null, "Sam file alignment of assembly graph to long reads",true);

		addStdHelp();
	}
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new NewScarfCmd ();
		args = cmdLine.stdParseLine(args);

		String fastgFile = cmdLine.getStringVal("fastg");
		String samFile = cmdLine.getStringVal("sam");
		HybridAssembler hbAss = new HybridAssembler(fastgFile);
		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
		//TODO: need to make this easier
		hbAss.assembly(samFile, 0);
	}
}
