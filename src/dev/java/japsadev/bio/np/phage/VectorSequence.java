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
	public VectorSequence(String seqFile, String bwaExe, int e5, int s3) throws IOException{
//		SequenceReader reader = SequenceReader.getReader(seqFile);
//		plasmid = reader.nextSequence(Alphabet.DNA());
		plasmidFile=seqFile;
		this.e5=e5;
		this.s5=this.e5-FLANKING;
		this.s3=s3;
		this.e3=this.s3+FLANKING;
//		reader.close();
		File out = new File("/home/s.hoangnguyen/Projects/Phage/log.out");
		bwaProcess  = new ProcessBuilder(bwaExe, 
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
				)					
				.redirectError(new File("/home/s.hoangnguyen/Projects/Phage/log.err"))
				.redirectOutput(out) // this is the BWA output for *reader*
//				.redirectInput(Redirect.INHERIT) //correspond to *toBWA*
				.start();
		toBWA = new SequenceOutputStream(bwaProcess.getOutputStream());
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

//		reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream())); // take the output from BWA
		reader = SamReaderFactory.makeDefault().open(SamInputResource.of(out)); // take the output from BWA

		iter = reader.iterator();
	}
	
	public Sequence extractInsertSequence(Sequence read){
		Sequence insert = null;
		try {
			//SequenceOutputStream stdout = new SequenceOutputStream(System.out);
			synchronized(this){
				read.writeFasta(toBWA);
			}

			//Logging.info("bwa started x");	
//			SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream()));
			boolean direction = true;
			int readFront = 0, readBack = 0;

//			SAMRecordIterator iter = reader.iterator();
			while (iter.hasNext()){
				SAMRecord record = iter.next();
				if (record.getReadString().length() < 10){
					System.out.println("== " + record.getReadName());
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
			//stdout.close();
//			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return insert;
	}
	
	public static void main(String[] args){
		try {
			VectorSequence vector = new VectorSequence("/home/s.hoangnguyen/Projects/Phage/plasmid.fasta","bwa",1658,2735);
			SequenceReader.readAll("/home/s.hoangnguyen/Projects/Phage/2d_1.fasta",Alphabet.DNA())
			.stream()
			.map(e -> vector.extractInsertSequence(e))
			.forEach(e -> System.out.println(e==null));
			
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
