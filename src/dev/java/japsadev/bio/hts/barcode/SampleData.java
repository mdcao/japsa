package japsadev.bio.hts.barcode;

import java.io.IOException;
import java.io.OutputStream;

import japsa.bio.hts.scaffold.RealtimeScaffolding;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;

public class SampleData {
	String id;
	OutputStream out;
	RealtimeScaffolding rtScaffold;
	Sequence fBarcode, rBarcode;
	
	public SampleData(){
		id = new String();
		out = null;
		rtScaffold = null;
		fBarcode = rBarcode = new Sequence(Alphabet.DNA(), 100);
	}
	
	public SampleData(String id, OutputStream out, String sequenceFile, Sequence fBarcode, Sequence rBarcode) throws IOException, InterruptedException{
		this.id = id;
		this.out = out;
		rtScaffold = new RealtimeScaffolding(sequenceFile, null, null, null, null, "-");
		this.fBarcode = fBarcode;
		this.rBarcode = rBarcode;
	}
	
	public String getId(){
		return id;
	}
	public OutputStream getOutputStream(){
		return out;
	}
	public RealtimeScaffolding getScaffoldingControl(){
		return rtScaffold;
	}
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
	public void setContigData(String spadesFile) throws IOException, InterruptedException{
		rtScaffold = new RealtimeScaffolding(spadesFile, null, null, null, null, "-");
	}
	public void setOutputStream(OutputStream out){
		this.out = out;
	}
}
