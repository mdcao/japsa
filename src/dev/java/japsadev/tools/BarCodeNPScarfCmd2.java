package japsadev.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Random;

import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;
import japsa.util.deploy.Deployable;

@Deployable(
		scriptName = "jsa.dev.barcodeScarf2", 
		scriptDesc = "Clustering nanopore sequences based on barcode"
		)
public class BarCodeNPScarfCmd2 extends CommandLine{
	public BarCodeNPScarfCmd2(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		addStdHelp();
	}
	public static void main(String[] args) throws IOException, InterruptedException{
		CommandLine cmdLine = new BarCodeNPScarfCmd2 ();
		//args = cmdLine.stdParseLine(args);

		//String bcFile = cmdLine.getStringVal("bcFile");
		//String seqFile = cmdLine.getStringVal("seqFile");
		String [] IDs = {
				"GN_093",
				"GN_101",
				"GN_132",
				"GN_106",
				"GN_096",
				"GN_133",
				"GN_092",
//				"GP_023"
		};
		//Set up these processes
		int nProcess = IDs.length;
		SequenceOutputStream [] outStrs = new SequenceOutputStream[nProcess];
		Process [] process = new Process[nProcess];

		for (int i = 0; i < nProcess;i++){
			ProcessBuilder pb = new ProcessBuilder("./script.py", IDs[i]);
			process[i]  = pb.start();
			Logging.info("Job for " + IDs[i]  + " started");
			outStrs[i] = new SequenceOutputStream(process[i].getOutputStream());
		}
		
		//Read data from standard in
		
		BufferedReader reader = SequenceReader.openFile("-");
		String line;
		Random random = new Random();
		while ((line = reader.readLine())!=null){
			//pretend index as the right barcode
			int index = random.nextInt(10);
			System.out.print("You entere " + line);
			if (index < nProcess){
				//Send data to the right process
				outStrs[index].print(line);
				//outStrs[index].flush();
				System.out.print(", it was sent to job " + index);
			}else
				System.out.print(", it wasnt sent anywhere");
			
		}
		System.out.println("Done all input");
		for (int i = 0; i < nProcess;i++){
			outStrs[i].close();
			process[i].waitFor();
		}
		System.out.println("Done every thing");
	}
}
