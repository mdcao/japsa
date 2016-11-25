package japsadev.bio.hts.barcode;

import java.io.IOException;
import java.io.OutputStream;

import japsa.bio.hts.scaffold.RealtimeScaffolding;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.Logging;

public class SampleData {
	String id;
	SequenceOutputStream streamToScaffolder, streamToFile;
	//RealtimeScaffolding rtScaffold;
	Process process=null;
	Sequence fBarcode, rBarcode;
	
	public SampleData(String id) throws IOException{
		this.id = id;
//		out = null;
//		rtScaffold = null;
		fBarcode = rBarcode = new Sequence(Alphabet.DNA(), 100);
		streamToFile = SequenceOutputStream.makeOutputStream(id+".fasta");
//		ProcessBuilder pb = new ProcessBuilder("./script.sh", id);
//		process  = pb.start();
//		out = new SequenceOutputStream(process.getOutputStream());
	}
	
	public SampleData(String id, OutputStream out, String sequenceFile, Sequence fBarcode, Sequence rBarcode) throws IOException, InterruptedException{
		this.id = id;
		//this.out = out;
		//rtScaffold = new RealtimeScaffolding(sequenceFile, null, null, null, null, "-");
		this.fBarcode = fBarcode;
		this.rBarcode = rBarcode;
	}
	
	public String getId(){
		return id;
	}
//	public OutputStream getOutputStream(){
//		return out;
//	}
//	public RealtimeScaffolding getScaffoldingControl(){
//		return rtScaffold;
//	}
	public Sequence getFBarcode(){
		return fBarcode;
	}
	public Sequence getRBarcode(){
		return rBarcode;
	}
	
	public void setId(String id){
		this.id = id;
	}
	public void setFBarcode(Sequence fBarcode){
		this.fBarcode = fBarcode;
	}
	public void setRBarcode(Sequence rBarcode){
		this.rBarcode = rBarcode;
	}
	
	public void passRead(Sequence seq) throws IOException{
		if(process==null){
			ProcessBuilder pb = new ProcessBuilder("./script.sh", id);
			process  = pb.start();
			streamToScaffolder = new SequenceOutputStream(process.getOutputStream());
			Logging.info("Process for sample " + getId() + " is started!");
		}
		
		if(process.isAlive() && seq.length() > 300){
			seq.writeFasta(streamToScaffolder);
			seq.writeFasta(streamToFile);
		}
		else
			return;
	}
	public boolean terminate(){
		int stat=-1;
		if(process!=null && process.isAlive()){
			try {
				streamToScaffolder.close();
				streamToFile.close();
				stat = process.waitFor();
			} catch (InterruptedException | IOException e) {
				// TODO Auto-generated catch block
				stat=-1;
				e.printStackTrace();
			}
		}
		if(stat==0) return true;
		else return false;
	}
//	public void setContigData(String spadesFile) throws IOException, InterruptedException{
//		rtScaffold = new RealtimeScaffolding(spadesFile, null, null, null, null, "-");
//	}
//	public void setOutputStream(OutputStream out){
//		this.out = out;
//	}
}
