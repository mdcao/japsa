package japsadev.tools;

import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import japsadev.bio.np.phage.ConsensusGenerator;

@Deployable(
        scriptName = "jsa.dev.consensus",
        scriptDesc = "Generate consensus sequences for grouping..."
)
public class ConsensusGenerateCmd extends CommandLine{
    public ConsensusGenerateCmd(){
        super();
        Deployable annotation = getClass().getAnnotation(Deployable.class);
        setUsage(annotation.scriptName() + " [options]");
        setDesc(annotation.scriptDesc());

        addString("input", null, "Name of the input FASTA file for all sequences", true);
        addString("list", null, "Name of the group file, each line containing id of component sequences from input", true);
        addString("output", "out.fasta", "Name of the output file, - for standard input");
        addString("aligner", "kalign", "Name of the aligner used for multiple-alignment in grouping phase");
        addBoolean("trim", false, "To trimg the flanking sequences");
        addStdHelp();
    }

    public static void main(String[] args) throws IOException {

        /*********************** Setting up script ****************************/
        CommandLine cmdLine = new ConsensusGenerateCmd();
        args = cmdLine.stdParseLine(args);
        /**********************************************************************/

        String 	input = cmdLine.getStringVal("input"),
                list = cmdLine.getStringVal("list"),
                output = cmdLine.getStringVal("output"),
                aligner = cmdLine.getStringVal("aligner");
        boolean trim = cmdLine.getBooleanVal("trim");
        try {
            ConsensusGenerator gen = new ConsensusGenerator();

            gen.generate(input, list, output, aligner, trim);

        } catch (Exception e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
    }

}



	/*RST*



	 
	  
	  
	  
	  
	  
	  
	  
	  
	*RST*/
	  

