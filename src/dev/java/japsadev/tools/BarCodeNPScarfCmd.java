package japsadev.tools;

import java.io.IOException;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(
		scriptName = "jsa.dev.barcodeScarf", 
		scriptDesc = "Clustering nanopore sequences based on barcode"
		)
public class BarCodeNPScarfCmd extends CommandLine{
	public BarCodeNPScarfCmd(){
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);		
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc()); 

		//addString("barcode", null, "Barcode file",true);		
		//addString("seqFile", null, "Nanopore sequences file",true);

		addStdHelp();
	}
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new BarCodeNPScarfCmd ();
		//args = cmdLine.stdParseLine(args);

		//String bcFile = cmdLine.getStringVal("bcFile");
		//String seqFile = cmdLine.getStringVal("seqFile");

		String bwaExe = "";


		NPScarfThread thread = new NPScarfThread(1);
		thread.start();		

	}

	static class NPScarfThread extends Thread{
		int index;

		NPScarfThread(int index) throws IOException{
			this.index = index;			
		}		/* (non-Javadoc)
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			ProcessBuilder pb = new ProcessBuilder("./script.sh","" + index);			
			try {
				Process p = pb.start();				
				p.waitFor();				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}//run		
	}
}
