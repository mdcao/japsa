package japsa.tools.bio.np;


import java.io.File;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import japsa.bio.np.ErrorCorrection;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;


@Deployable(
        scriptName = "jsa.np.consensus",
        scriptDesc = "Generate consensus sequences from a multi-FASTA"
)
public class ConsensusCmd extends CommandLine{
    public ConsensusCmd(){
        super();
        Deployable annotation = getClass().getAnnotation(Deployable.class);
        setUsage(annotation.scriptName() + " [options]");
        setDesc(annotation.scriptDesc());

        addString("input", null, "Name of the input FASTA file for all sequences", true);
        addString("prefix", "out", "Prefix for the output files");
        addString("aligner", "kalign", "Name of the aligner used for multiple-alignment in grouping phase");
        addStdHelp();
    }

    public static void main(String[] args) {

        /*********************** Setting up script ****************************/
        CommandLine cmdLine = new ConsensusCmd();
        args = cmdLine.stdParseLine(args);
        /**********************************************************************/

        String 	input = cmdLine.getStringVal("input"),
                prefix = cmdLine.getStringVal("prefix"),
                aligner = ErrorCorrection.msa = cmdLine.getStringVal("aligner");
        try {
			SequenceReader faiReader =  SequenceReader.getReader(input);
			Sequence seq=null;
			List<Sequence> seqs = new ArrayList<Sequence>();
			while((seq=faiReader.nextSequence(Alphabet.DNA()))!=null)
				seqs.add(seq);			
			faiReader.close();
			String tmpDir=System.getProperty("java.io.tmpdir");
			if(tmpDir.isEmpty())
				tmpDir=Files.createTempDirectory("tmp").toString();
			seq=ErrorCorrection.consensusSequence(seqs, tmpDir+File.separator+prefix, aligner);
			seq.setName(prefix);
			seq.setDesc("length="+seq.length());
			SequenceOutputStream out = new SequenceOutputStream(new FileOutputStream(prefix+".fasta",false));
			seq.writeFasta(out);
			out.close();
        } catch (Exception e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
    }

}

/*RST*
 
  
  
  
  
  
  
  
  
*RST*/