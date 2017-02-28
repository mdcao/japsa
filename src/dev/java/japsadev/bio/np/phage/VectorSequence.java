package japsadev.bio.np.phage;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Date;
import java.util.stream.Stream;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.HTSUtilities;
import japsa.util.Logging;

public class VectorSequence {
	static final int FLANKING=100;
	//Sequence plasmid;
	String plasmidFile; //plasmid fasta/q file that are already indexed by bwa
	int s5, e5,
		s3, e3;
	Process bwaProcess = null;
	SequenceOutputStream toBWA = null;
	SamReader reader = null;
	SAMRecordIterator iter = null;
	String currentName = "";
	Thread extracting=null;
	
	public VectorSequence(String seqFile, String bwaExe, int e5, int s3) throws IOException{
//		SequenceReader reader = SequenceReader.getReader(seqFile);
//		plasmid = reader.nextSequence(Alphabet.DNA());
		plasmidFile=seqFile;
		this.e5=e5;
		this.s5=this.e5-FLANKING;
		this.s3=s3;
		this.e3=this.s3+FLANKING;
//		reader.close();
//		File out = new File("/home/s.hoangnguyen/Projects/Phage/log.out");
		ProcessBuilder pb  = new ProcessBuilder(bwaExe, 
				"mem",
				"-t4",
				"-k11",
				"-W20",
				"-r10",
				"-A1",
				"-B1",
				"-O1",
				"-E1",
				"-L0",
				"-a",
				"-Y",
				plasmidFile,
				"-"
				);
		bwaProcess = pb.redirectError(Redirect.to(new File("/dev/null"))).start();
		toBWA = new SequenceOutputStream(bwaProcess.getOutputStream());
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		extracting = new Thread(){
			public void run(){
				try {					
						reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream())); // take the output from BWA
						iter = reader.iterator();

						boolean direction = true;
						int readFront = 0, readBack = 0;

						while (iter.hasNext()){
							SAMRecord record = iter.next();
							if (record.getReadString().length() < 10){
								continue;//while
							}
							
							if (!record.getReadName().equals(currentName)){
								currentName = record.getReadName();
								direction = true;
								readFront = 0;
								readBack = 0;
							}
							System.out.println("Read " + currentName + " length " + record.getReadLength() + ": (" + record.getAlignmentStart() + ", " + record.getAlignmentEnd() + ")");
							
						}
					}catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}		
		};
		
	}
	
	public Sequence extractInsertSequence(Sequence read){
		Sequence insert = null;

			try {
				read.writeFasta(toBWA);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
		//TODO: how to return the insert sequence here out of above thread????	
		return insert;
	}
	
	public static void main(String[] args){
		try {
			VectorSequence vector = new VectorSequence("/home/s.hoangnguyen/Projects/Phage/plasmid.fasta","bwa",1658,2735);
			vector.extracting.start();
			
			SequenceReader.readAll("/home/s.hoangnguyen/Projects/Phage/2d_1.fasta",Alphabet.DNA())
			.stream()
			.map(e -> vector.extractInsertSequence(e))
			.forEach(e -> System.out.println(e==null));
				
			vector.extracting.join();
			
			if(vector.bwaProcess.isAlive()){
				vector.reader.close();
				vector.toBWA.close();
				vector.bwaProcess.waitFor();
			}
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
	}
}
